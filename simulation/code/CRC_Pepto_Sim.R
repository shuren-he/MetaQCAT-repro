rm(list=ls())
library(MetaQCAT)

# Simulate and evaluate the FE/RE MetaQCAT components and omnibus tests.

# Locate the project from this script (Rscript/source) or the working directory.
# This keeps paths portable if the project folder is moved or renamed.
project.root <- local({
  source.files <- vapply(sys.frames(), function(frame) {
    if (is.null(frame$ofile)) "" else as.character(frame$ofile)[1L]
  }, character(1))
  script.args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  script.files <- c(rev(source.files[nzchar(source.files)]),
                    sub("^--file=", "", script.args))
  start.dirs <- unique(c(dirname(script.files), getwd()))
  for (start.dir in start.dirs) {
    if (!dir.exists(start.dir)) next
    candidate <- normalizePath(start.dir, winslash = "/", mustWork = TRUE)
    repeat {
      if (all(dir.exists(file.path(candidate,
                                   c("simulation/data"))))) {
        return(candidate)
      }
      parent <- dirname(candidate)
      if (identical(parent, candidate)) break
      candidate <- parent
    }
  }
  stop("Cannot locate the MetaQCAT project. Run from the project directory or source this script by its file path.")
})

data.dir <- file.path(project.root, "simulation", "data")

seed <- 2333
lineage.name <- c("Peptostreptococcaceae")
nthreads <- 125

load(file.path(data.dir, paste0("mod_gdm_lst_", lineage.name, ".RData")))
load(file.path(data.dir, "meta.rda"))
load(file.path(data.dir, "count.genus.rda"))

study.names <- c("AUS", "CHI", "FRA", "GER", "USA")
lst.num <- as.integer(table(factor(meta$Country, levels = study.names)))
n.taxa <- length(mod_gdm_lst[[1]]$col.order)
taxa.names <- paste0("taxa", seq_len(n.taxa))

stopifnot(length(mod_gdm_lst) == length(lst.num),
          length(taxa.names) == n.taxa)

# Use the rounded mean composition depth of each study for every simulated
# sample in that study, using the same three genera used to fit the GDM.
composition.genera <- c(
  "1257 Peptostreptococcus", "1505657 Intestinibacter",
  "NA unknown Peptostreptococcaceae"
)
stopifnot(n.taxa == length(composition.genera),
          all(composition.genera %in% colnames(count.genus)),
          identical(rownames(count.genus), make.names(meta$Sample_ID)))
observed.depth <- rowSums(count.genus[, composition.genera, drop = FALSE])
depth.pool <- lapply(study.names, function(study) {
  as.integer(observed.depth[meta$Country == study])
})
stopifnot(all(lengths(depth.pool) == lst.num))

# "Increase by Delta" means multiplying the selected proportion by 1 + Delta.
# The small-count signal draws a new delta for each affected sample.
delta.proportion.vec <- 0.75
count.prob.vec <- 0.05

n.perm <- 300
simulation.target <- "power"     # use "type1" for the null simulation
if (!simulation.target %in% c("power", "type1")) {
  stop("simulation.target must be 'power' or 'type1'.")
}
sim.time <- if (simulation.target == "power") 1000L else 5000L


stu.num <- length(lst.num)
output.dir <- file.path(project.root, "simulation", "outcome", "pos_part_res")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

safe_test_call <- function(sim.cnt, test.fun, expected.length, test.label,
                           ...){
  tryCatch({
    value <- as.numeric(test.fun(sim.cnt = sim.cnt, ...))
    if (length(value) != expected.length) {
      stop(
        sprintf(
          "%s returned %d values; expected %d",
          test.label, length(value), expected.length
        )
      )
    }
    list(value = value, error = NA_character_)
  }, error = function(e) {
    list(
      value = rep(NA_real_, expected.length),
      error = conditionMessage(e)
    )
  })
}

unpack_safe_results <- function(output, expected.length, test.label){
  values <- matrix(
    NA_real_, nrow = length(output), ncol = expected.length
  )
  errors <- rep(NA_character_, length(output))

  for (i in seq_along(output)) {
    item <- output[[i]]
    if (is.list(item) &&
        is.numeric(item$value) &&
        length(item$value) == expected.length) {
      values[i, ] <- item$value
      errors[i] <- item$error
    } else {
      errors[i] <- paste0(
        "Parallel worker returned an invalid result: ",
        paste(as.character(item), collapse = " ")
      )
    }
  }

  failed <- which(!is.na(errors))
  if (length(failed) > 0L) {
    message(
      test.label, ": ", length(failed), "/", length(output),
      " replicates failed and were replaced with NA."
    )
  }

  list(
    values = values,
    errors = if (length(failed) == 0L) {
      data.frame(
        test = character(), replicate = integer(), error = character()
      )
    } else {
      data.frame(
        test = test.label,
        replicate = failed,
        error = errors[failed],
        stringsAsFactors = FALSE
      )
    }
  )
}

r.GDM <- function(X, mod, SeqDepth){
  alpha <- mod$alpha.est
  phi <- mod$phi.est
  if (is.null(alpha) || is.null(phi)) {
    stop("The fitted study model does not contain alpha.est and phi.est.")
  }
  alpha <- as.matrix(alpha)
  N <- nrow(X)
  K <- ncol(alpha)
  if (length(SeqDepth) == 1L) SeqDepth <- rep(SeqDepth, N)
  if (length(SeqDepth) != N) {
    stop("The sequencing-depth vector must have one value per sample.")
  }
  if (ncol(X) != nrow(alpha) || length(phi) != K) {
    stop("The GDM design matrix and fitted parameters have incompatible dimensions.")
  }
  mu.mat <- plogis(X %*% alpha)
  phi.mat <- matrix(phi, nrow = N, ncol = K, byrow = TRUE)
  Z <- matrix(rbeta(N * K, mu.mat * phi.mat,
                    (1 - mu.mat) * phi.mat), N, K)
  P <- matrix(NA_real_, N, K + 1L)
  remaining <- rep(1, N)
  for (j in seq_len(K)) {
    P[, j] <- Z[, j] * remaining
    remaining <- remaining * (1 - Z[, j])
  }
  P[, K + 1L] <- remaining
  P / rowSums(P)
}

scenario_directions <- function(stu.num, type, direct = NULL){
  if (type == "fixed") {
    return(rep(1L, stu.num))
  }
  if (type != "random") {
    stop("type must be 'fixed' or 'random'.")
  }
  if (identical(direct, "single")) {
    # Heterogeneity approach 1: two randomly selected studies have no effect.
    ans <- rep(0L, stu.num)
    ans[sample(seq_len(stu.num), stu.num - 2L)] <- 1L
    return(ans)
  }
  if (identical(direct, "both")) {
    # Heterogeneity approach 2: reverse the effect in one randomly selected
    # study; all remaining studies retain the replicate's baseline direction.
    ans <- rep(1L, stu.num)
    ans[sample(seq_len(stu.num), 1L)] <- -1L
    return(ans)
  }
  stop("For type='random', direct must be 'single' or 'both'.")
}

multinomial_counts <- function(P, depth){
  ans <- t(vapply(seq_len(nrow(P)), function(i) {
    as.integer(rmultinom(1L, size = depth[i], prob = P[i, ])[, 1L])
  }, integer(ncol(P))))
  colnames(ans) <- colnames(P)
  ans
}

renormalize_count_row <- function(x, depth){
  scaled <- x / sum(x) * depth
  ans <- floor(scaled)
  remainder <- as.integer(depth - sum(ans))
  if (remainder > 0L) {
    idx <- order(scaled - ans, decreasing = TRUE)[seq_len(remainder)]
    ans[idx] <- ans[idx] + 1L
  }
  as.integer(ans)
}

simulate_crc_data <- function(sim.cnt, mod_gdm_lst, depth.pool, taxa.names,
                              lst.num, signal.size, signal.method, type,
                              direct = NULL, seed, correlated = FALSE){
  set.seed(seed + sim.cnt)
  stu.num <- length(lst.num)
  n.taxa <- length(taxa.names)
  if (length(mod_gdm_lst) != stu.num || length(depth.pool) != stu.num ||
      any(lengths(depth.pool) == 0L)) {
    stop("GDM models, observed depth pools, and studies must align.")
  }
  if (!all(vapply(mod_gdm_lst, function(mod) {
    length(mod$col.order) == n.taxa
  }, logical(1)))) {
    stop("All fitted GDM study models must contain the same number of taxa.")
  }
  if (!signal.method %in% c("proportion", "count")) {
    stop("signal.method must be 'proportion' or 'count'.")
  }
  if (signal.method == "proportion" &&
      (!is.finite(signal.size) || signal.size < 0)) {
    stop("The proportional increase Delta must be non-negative.")
  }
  if (signal.method == "count" &&
      (!is.finite(signal.size) || signal.size < 0 || signal.size > 1)) {
    stop("The binomial count-signal probability must be between 0 and 1.")
  }

  target.idx <- sample(seq_len(n.taxa), 1L)
  baseline.direction <- sample(c(-1L, 1L), 1L)
  directions <- baseline.direction * scenario_directions(stu.num, type, direct)
  count.trials <- as.integer(round(mean(unlist(depth.pool, use.names = FALSE))))
  dat.lst <- vector("list", stu.num)
  X.lst <- vector("list", stu.num)
  keep.lst <- vector("list", stu.num)
  depth.lst <- vector("list", stu.num)
  cluster.lst <- vector("list", stu.num)

  for (l in seq_len(stu.num)) {
    n <- lst.num[l]
    # GDM baseline proportions retain the study-specific fitted batch effects.
    # All samples in a study have its rounded observed mean composition depth.
    depth.all <- rep.int(as.integer(round(mean(depth.pool[[l]]))), n)
    P.all <- r.GDM(
      X = cbind(1, rep(0, n)),
      mod = mod_gdm_lst[[l]]$est.para,
      SeqDepth = depth.all
    )
    P.all <- P.all[, order(mod_gdm_lst[[l]]$col.order), drop = FALSE]
    colnames(P.all) <- taxa.names
    keep <- which(depth.all >= 1.5 * n.taxa)
    if (length(keep) == 0L) {
      stop(sprintf("Study %d has no samples after the depth filter.", l))
    }
    depth <- depth.all[keep]
    P <- P.all[keep, , drop = FALSE]
    if (correlated) {
      # Random pairs share an independent baseline draw from the same
      # study-specific generalized Dirichlet distribution. An unpaired final
      # sample (when n is odd) remains independent in a singleton cluster.
      shuffled <- sample.int(length(keep))
      pairs <- split(shuffled, ceiling(seq_along(shuffled) / 2L))
      cluster <- integer(length(keep))
      for (i in seq_along(pairs)) {
        members <- pairs[[i]]
        cluster[members] <- i
        if (length(members) == 2L) {
          shared <- r.GDM(
            X = matrix(c(1, 0), nrow = 1L),
            mod = mod_gdm_lst[[l]]$est.para,
            SeqDepth = depth[members[1L]]
          )
          shared <- shared[1L, order(mod_gdm_lst[[l]]$col.order)]
          P[members, ] <- (P[members, , drop = FALSE] +
            matrix(shared, nrow = 2L, ncol = n.taxa, byrow = TRUE)) / 2
        }
      }
      # Cluster-level labels keep paired samples in the same group.
      cluster.case <- rep(0L, length(pairs))
      cluster.case[sample.int(length(pairs), ceiling(length(pairs) / 2L))] <- 1L
      X <- cluster.case[cluster]
    } else {
      cluster <- seq_along(keep)
      X <- rep(0L, length(keep))
      X[sample.int(length(keep), ceiling(length(keep) / 2L))] <- 1L
    }
    affected <- if (directions[l] == 1L) X == 1L else X == 0L

    if (directions[l] != 0L && signal.size > 0 &&
        signal.method == "proportion") {
      P[affected, target.idx] <- P[affected, target.idx] * (1 + signal.size)
      P[affected, ] <- P[affected, , drop = FALSE] /
        rowSums(P[affected, , drop = FALSE])
    }

    Y <- multinomial_counts(P, depth)
    if (directions[l] != 0L && signal.size > 0 &&
        signal.method == "count") {
      affected.idx <- which(affected)
      delta <- rbinom(length(affected.idx), size = count.trials,
                      prob = signal.size)
      for (j in seq_along(affected.idx)) {
        i <- affected.idx[j]
        modified <- Y[i, ]
        modified[target.idx] <- modified[target.idx] + delta[j]
        Y[i, ] <- renormalize_count_row(modified, depth[i])
      }
    }

    stopifnot(all(rowSums(Y) == depth))
    dat.lst[[l]] <- Y
    X.lst[[l]] <- matrix(X, ncol = 1L)
    keep.lst[[l]] <- keep
    depth.lst[[l]] <- depth
    cluster.lst[[l]] <- cluster
  }

  list(OTU = dat.lst, X = X.lst,
       keep = keep.lst, depth = depth.lst, cluster = cluster.lst,
       correlated = correlated,
       target.taxon = taxa.names[target.idx],
       baseline.direction = baseline.direction, directions = directions)
}

ACAT<-function(Pvals,Weights=NULL){
  #### check if there is NA
  if (sum(is.na(Pvals))>0){
    stop("Cannot have NAs in the p-values!")
  }
  #### check if Pvals are between 0 and 1
  if ((sum(Pvals<0)+sum(Pvals>1))>0){
    stop("P-values must be between 0 and 1!")
  }
  #### check if there are pvals that are either exactly 0 or 1.
  is.zero<-(sum(Pvals==0)>=1)
  is.one<-(sum(Pvals==1)>=1)
  if (is.zero && is.one){
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero){
    return(0)
  }
  if (is.one){
    warning("There are p-values that are exactly 1!")
    return(1)
  }

  #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
  if (is.null(Weights)){
    Weights<-rep(1/length(Pvals),length(Pvals))
  }else if (length(Weights)!=length(Pvals)){
    stop("The length of weights should be the same as that of the p-values")
  }else if (sum(Weights<0)>0){
    stop("All the weights must be positive!")
  }else{
    Weights<-Weights/sum(Weights)
  }


  #### check if there are very small non-zero p values
  is.small<-(Pvals<1e-16)
  if (sum(is.small)==0){
    cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
  }else{
    cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
    cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
  }
  #### check if the test statistic is very large.
  if (cct.stat>1e+15){
    pval<-(1/cct.stat)/pi
  }else{
    pval<-1-pcauchy(cct.stat)
  }
  return(pval)
}

meta.qcat.test <- function(sim.cnt, mod_gdm_lst, depth.pool, taxa.names,
                           lst.num, signal.size, signal.method, type,
                           direct = NULL, seed, simulated = NULL){
  if (is.null(simulated)) {
    simulated <- simulate_crc_data(
      sim.cnt, mod_gdm_lst, depth.pool, taxa.names, lst.num, signal.size,
      signal.method, type, direct, seed
    )
  }
  dat.simu.poured.tmp <- simulated$OTU
  X.lst.tmp <- simulated$X
  cluster.ids <- if (isTRUE(simulated$correlated)) {
    lapply(simulated$cluster, factor)
  } else {
    NULL
  }
  permutation.strata <- if (isTRUE(simulated$correlated)) "between" else NULL
  fe.meta = MetaQCAT(OTU = dat.simu.poured.tmp, X= X.lst.tmp, X.index = 1, Tax=NULL,
                           meta.method = 'FE-MetaQCAT', class.id = cluster.ids,
                           permute.strata = permutation.strata,
                           min.depth=0, n.perm=n.perm, n.threads = 1)
  fe.zero.meta = MetaQCAT(OTU = dat.simu.poured.tmp, X= X.lst.tmp, X.index = 1, Tax=NULL,
                                meta.method = 'FE-MetaQCAT*', class.id = cluster.ids,
                                permute.strata = permutation.strata,
                                min.depth=0, n.perm=n.perm, n.threads = 1)
  re.meta = MetaQCAT(OTU = dat.simu.poured.tmp, X= X.lst.tmp, X.index = 1, Tax=NULL,
                           meta.method = 'RE-MetaQCAT', class.id = cluster.ids,
                           permute.strata = permutation.strata,
                           min.depth=0, n.perm=n.perm, n.threads = 1)
  re.zero.meta = MetaQCAT(OTU = dat.simu.poured.tmp, X= X.lst.tmp, X.index = 1, Tax=NULL,
                                meta.method = 'RE-MetaQCAT*', class.id = cluster.ids,
                                permute.strata = permutation.strata,
                                min.depth=0, n.perm=n.perm, n.threads = 1)
  return(c(fe.meta$pval[1],  fe.zero.meta$pval[1], ACAT(na.omit(c(fe.meta$pval[1],  fe.zero.meta$pval[1]))),
           fe.meta$pval[2],  fe.zero.meta$pval[2], ACAT(na.omit(c(fe.meta$pval[2],  fe.zero.meta$pval[2]))),
           re.meta$pval[1],  re.zero.meta$pval[1], ACAT(na.omit(c(re.meta$pval[1],  re.zero.meta$pval[1])))))

}

run_replicate_tests <- function(sim.cnt, mod_gdm_lst, depth.pool, taxa.names,
                                lst.num, signal.size, signal.method, type,
                                direct = NULL, seed, correlated = FALSE) {
  simulated <- tryCatch(
    simulate_crc_data(sim.cnt, mod_gdm_lst, depth.pool, taxa.names,
                      lst.num, signal.size, signal.method, type,
                      direct, seed, correlated = correlated),
    error = function(e) e
  )
  if (inherits(simulated, "error")) {
    return(list(value = rep(NA_real_, 9L),
                error = paste0("Data generation: ", conditionMessage(simulated))))
  }
  safe_test_call(
    sim.cnt = sim.cnt, test.fun = meta.qcat.test,
    expected.length = 9L, test.label = "MetaQCAT",
    mod_gdm_lst = mod_gdm_lst, depth.pool = depth.pool, taxa.names = taxa.names,
    lst.num = lst.num, signal.size = signal.size, signal.method = signal.method,
    type = type, direct = direct, seed = seed, simulated = simulated
  )
}

test.lst <- list(
  homogeneous = c(type = "fixed", direct = NA_character_),
  two.null.studies = c(type = "random", direct = "single"),
  one.flipped.study = c(type = "random", direct = "both")
)

signal.grid <- list(
  proportion = delta.proportion.vec,
  count = count.prob.vec
)

if (simulation.target == "type1") {
  # Under the null there is no signal, so one homogeneous configuration is
  # sufficient; the manuscript uses 5000 replicates for this calculation.
  test.lst <- test.lst["homogeneous"]
  signal.grid <- list(proportion = 0)
}

if (!isTRUE(getOption("crc.simulation.source.only", FALSE))) {
for (scenario.name in names(test.lst)) {
  # scenario.name = names(test.lst)[1]
  type <- unname(test.lst[[scenario.name]]["type"])
  direct <- unname(test.lst[[scenario.name]]["direct"])
  if (is.na(direct)) direct <- NULL
  for (signal.method in names(signal.grid)) {
    for (signal.size in signal.grid[[signal.method]]) {
      # signal.method = names(signal.grid)[1]
      # signal.size = signal.grid[[signal.method]][1]
      message("Scenario: ", scenario.name, "; signal: ", signal.method,
              "; ", if (signal.method == "count") "p" else "Delta",
              ": ", signif(signal.size, 5))
      asy.res <- perm.res <- NULL
      file.tag <- paste(
        lineage.name, simulation.target, signal.method, scenario.name,
        format(signal.size, scientific = FALSE), sep = "_"
      )

      output <- parallel::mclapply(
        seq_len(sim.time), run_replicate_tests,
        mod_gdm_lst = mod_gdm_lst, depth.pool = depth.pool,
        taxa.names = taxa.names, lst.num = lst.num,
        signal.size = signal.size,
        signal.method = signal.method, type = type, direct = direct,
        seed = seed, mc.cores = nthreads, mc.silent = TRUE
      )
      unpacked <- unpack_safe_results(output, 9L, "MetaQCAT")
      asy.res <- unpacked$values[, 1:3, drop = FALSE]
      perm.res <- unpacked$values[, 4:9, drop = FALSE]
      colnames(asy.res) <- c(
        "FE-MetaQCAT (asy)", "FE-MetaQCAT* (asy)", "FE-MetaQCAT-O (asy)"
      )
      colnames(perm.res) <- c(
        "FE-MetaQCAT (perm)", "FE-MetaQCAT* (perm)", "FE-MetaQCAT-O (perm)",
        "RE-MetaQCAT", "RE-MetaQCAT*", "RE-MetaQCAT-O"
      )

      asy.valid <- colSums(!is.na(asy.res))
      perm.valid <- colSums(!is.na(perm.res))
      print(ifelse(
        asy.valid > 0,
        colSums(asy.res < 0.05, na.rm = TRUE) / asy.valid,
        NA_real_
      ))
      print(ifelse(
        perm.valid > 0,
        colSums(perm.res < 0.05, na.rm = TRUE) / perm.valid,
        NA_real_
      ))

      write.csv(asy.res, row.names = FALSE,
                file = file.path(output.dir, paste0(file.tag, "_asy_MetaQCAT_org.csv")))
      write.csv(perm.res, row.names = FALSE,
                file = file.path(output.dir, paste0(file.tag, "_perm_MetaQCAT_org.csv")))

      failure.log <- unpacked$errors
      if (!is.null(failure.log) && nrow(failure.log) > 0L) {
        failure.log$scenario <- scenario.name
        failure.log$signal.method <- signal.method
        failure.log$signal.size <- signal.size
        write.csv(
          failure.log, row.names = FALSE,
          file = file.path(output.dir, paste0(file.tag, "_failures.csv"))
        )
      }
    }
  }
}
}
