rm(list=ls())
library(doRNG)
library(doParallel)
library(foreach)
library(HMP)
library(MetaQCAT)
library(MMUPHin)
library(miLineage)

# MetaQCAT() prints routine filtering notices with cat() and has no verbose
# argument. Create a local equivalent with only cat() calls removed; filtering
# and all statistical calculations remain unchanged.

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
lineage.name <- c("Enterobacteriaceae")
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

# The four genera used to fit this GDM exclude Hafnia. Use the rounded mean
# composition depth of each study for every simulated sample in that study.
composition.genera <- c(
  "NA Enterobacteriaceae gen. [C Escherichia/Shigella]",
  "544 Citrobacter", "570 Klebsiella", "547 Enterobacter"
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
delta.proportion.vec <- 1.5
count.prob.vec <- 0.005

n.perm <- 300
simulation.target <- "power"     # use "type1" for the null simulation
if (!simulation.target %in% c("power", "type1")) {
  stop("simulation.target must be 'power' or 'type1'.")
}
sim.time <- if (simulation.target == "power") 1000L else 5000L


stu.num <- length(lst.num)
output.dir <- file.path(project.root, "simulation", "outcome", "pos_part_res")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

Rarefy <- function (otu.tab, depth = min(rowSums(otu.tab)), seed = 1234){
  # Rarefaction function: downsample to equal depth
  #
  # Args:
  #		otu.tab: OTU count table, row - n sample, column - q OTU
  #		depth: required sequencing depth
  #
  # Returns:
  # 	otu.tab.rff: Rarefied OTU table
  #		discard: labels of discarded samples
  #
  otu.tab <- as.matrix(otu.tab)
  ind <- (rowSums(otu.tab) < depth)
  sam.discard <- rownames(otu.tab)[ind]
  otu.tab <- otu.tab[!ind, ]
  set.seed(seed)
  rarefy <- function(x, depth){
    y <- sample(rep(1:length(x), x), depth)
    y.tab <- table(y)
    z <- numeric(length(x))
    z[as.numeric(names(y.tab))] <- y.tab
    z
  }
  otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
  rownames(otu.tab.rff) <- rownames(otu.tab)
  colnames(otu.tab.rff) <- colnames(otu.tab)
  return(list(otu.tab.rff=otu.tab.rff, discard=sam.discard))
}

dm1.resample <- function(Y, case, n.perm=1000){

  group.data = list(grp1 = Y[case==0,], grp2 = Y[case==1,] )
  dm.stat = try(Xmcupo.sevsample(group.data)$`Xmcupo statistics`,
                silent = TRUE)
  if((!"try-error" %in% class(dm.stat)) & !is.na(dm.stat) ){

    n.dm = 0
    dm.acc = 0
    #print("simulated stat:")
    #set.seed(16)


    for(k in 1:n.perm){

      case.perm = sample(case)
      group.data.perm = list(grp1 = Y[case.perm==0,], grp2 = Y[case.perm==1,] )
      dm.stat.perm = try(
        Xmcupo.sevsample(group.data.perm)$`Xmcupo statistics`,
        silent = TRUE
      )

      if((!"try-error" %in% class(dm.stat.perm)) & !is.na(dm.stat.perm) ){

        n.dm = n.dm + 1
        if(dm.stat.perm >= dm.stat){
          dm.acc = dm.acc + 1
        }
      }

    }

    if(n.dm<n.perm/2){
      print("#replicate too small for dm test")
    }

    dm.Rpvalue = dm.acc/n.dm


  }else{

    dm.Rpvalue = NA

  }

  return(dm.Rpvalue)
}

F.test <- function(x){
  x <- as.numeric(na.omit(x))
  x <- x[is.finite(x) & x >= 0 & x <= 1]
  if (length(x) == 0L) return(NA_real_)
  if (any(x == 0)) return(0)
  x.stat = -2 * sum(log(x))
  return( 1 - pchisq(x.stat, df = 2 * length(x)) )
}

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

meta.tab <- as.data.frame(meta)
meta.tab$Sample_ID <- gsub("-", ".", meta$Sample_ID)
meta.tab$External_ID <- gsub("-", ".", meta.tab$External_ID)
rownames(meta.tab) <- paste0(meta.tab$Sample_ID,'.', meta.tab$External_ID)
meta.tab.lst <- split.data.frame(meta.tab, meta.tab$Country)[study.names]

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

p.comb.test <- function(sim.cnt, mod_gdm_lst, depth.pool, taxa.names,
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
  stu.num <- length(mod_gdm_lst)
  score.pvalue.beta = NULL
  score.pvalue.both = NULL
  score.pvalue.beta.perm  = NULL
  score.pvalue.both.perm = NULL
  dm.pvalue = NULL
  dm.pvalue.perm = NULL
  ZIGDM.pvalue = NULL
  ZIGDM.pvalue.perm = NULL
  for(j in 1:stu.num){
    Y = dat.simu.poured.tmp[[j]]
    Y.rare = Rarefy(Y)$otu.tab.rff
    X = X.lst.tmp[[j]]
    tmp.q1 = try(
      miLineage::QCAT(Y, X, 1, Tax = NULL, min.depth = 0,
                     n.perm = n.perm),
      silent = TRUE
    )
    tmp.q2 = try(
      miLineage::QCAT_GEE(Y.rare, X, 1, X, 1, Tax = NULL,
                          min.depth = 0, n.perm = n.perm),
      silent = TRUE
    )
    tmp.dm.perm = try(
      dm1.resample(Y, X, n.perm = n.perm),
      silent = TRUE
    )
    tmp.ZIGDM = try(
      miLineage::ZIGDM(
        Y, as.matrix(X, ncol = 1), NULL, NULL, Tax = NULL,
        X.index = 1, n.perm = n.perm
      ),
      silent = TRUE
    )
    Y = split.data.frame(Y,X)
    tmp.dm = try(HMP::Xmcupo.sevsample(Y), silent = TRUE)
    if(! "try-error" %in% class(tmp.q1)){
      score.pvalue.beta = append(score.pvalue.beta,tmp.q1$pval[1])
      score.pvalue.beta.perm = append(score.pvalue.beta.perm,tmp.q1$pval[2])
    }
    if(! "try-error" %in% class(tmp.q2)){
      score.pvalue.both = append(score.pvalue.both, tmp.q2$pval[1])
      score.pvalue.both.perm = append(score.pvalue.both.perm, tmp.q2$pval[2])
    }
    if(! "try-error" %in% class(tmp.dm)){
      dm.pvalue = append(dm.pvalue, tmp.dm$`p value`)
    }
    if(! "try-error" %in% class(tmp.dm.perm)){
      dm.pvalue.perm  = append(dm.pvalue.perm,tmp.dm.perm)
    }
    if(! "try-error" %in% class(tmp.ZIGDM)){
      ZIGDM.pvalue = append(ZIGDM.pvalue, tmp.ZIGDM$pval[1])
      ZIGDM.pvalue.perm = append(ZIGDM.pvalue.perm, tmp.ZIGDM$pval[2])
    }
  }
  score.pvalue.q1 = F.test(na.omit(score.pvalue.beta))
  score.pvalue.q1.perm = F.test(na.omit(score.pvalue.beta.perm))
  score.pvalue.q2 = F.test(na.omit(score.pvalue.both))
  score.pvalue.q2.perm = F.test(na.omit(score.pvalue.both.perm))
  dm.pvalue = F.test(na.omit(dm.pvalue))
  dm.pvalue.perm = F.test(na.omit(dm.pvalue.perm))
  ZIGDM.pvalue = F.test(na.omit(ZIGDM.pvalue))
  ZIGDM.pvalue.perm = F.test(na.omit(ZIGDM.pvalue.perm))
  return(c(score.pvalue.q1,  score.pvalue.q2, dm.pvalue, ZIGDM.pvalue,
           score.pvalue.q1.perm, score.pvalue.q2.perm, dm.pvalue.perm, ZIGDM.pvalue.perm))
}

cd.test <- function(sim.cnt, mod_gdm_lst, depth.pool, taxa.names, lst.num,
                    meta.tab.lst, signal.size, signal.method, type,
                    direct = NULL, seed, simulated = NULL){
  if (is.null(simulated)) {
    simulated <- simulate_crc_data(
      sim.cnt, mod_gdm_lst, depth.pool, taxa.names, lst.num, signal.size,
      signal.method, type, direct, seed
    )
  }
  dat.simu.poured.tmp <- simulated$OTU
  X.lst.tmp <- simulated$X
  count.genus = do.call(rbind,dat.simu.poured.tmp)
  meta.tab.lst.tmp <- Map(
    function(dat, keep) dat[keep, , drop = FALSE],
    meta.tab.lst, simulated$keep
  )
  meta.tab.tmp <- do.call(rbind, meta.tab.lst.tmp)
  X = as.matrix(do.call(c,X.lst.tmp))
  # MMUPHin must adjust for the simulated case/control assignment, not the
  # original CRC labels stored in meta.rda.
  meta.tab.tmp$Group <- factor(ifelse(X[, 1] == 1, "case", "control"))
  otu.tab <-t(count.genus)
  colnames(otu.tab) <- rownames(meta.tab.tmp)
  fit_adjust_batch <- MMUPHin::adjust_batch(feature_abd = otu.tab,
                                            batch = "Country",
                                            covariates = "Group",
                                            data = meta.tab.tmp,
                                            control = list(verbose = FALSE))
  taxa.MMUPHin <- fit_adjust_batch$feature_abd_adj
  Y = t(taxa.MMUPHin)
  Y.rare = Rarefy(Y)$otu.tab.rff
  q1.MMUPHin = miLineage::QCAT(Y, X, 1, Tax = NULL,min.depth = 0, n.perm = n.perm)
  q2.MMUPHin  = miLineage::QCAT_GEE(Y.rare, X, 1, X, 1, Tax = NULL,
                                    min.depth = 0, n.perm = n.perm)
  ZIGDM.MMUPHin = miLineage::ZIGDM(Y, as.matrix(X, ncol = 1),  NULL, NULL, Tax = NULL, X.index = 1, n.perm =  n.perm)
  dm.MMUPHin.res = dm1.resample(Y, X, n.perm = n.perm)
  Y = split.data.frame(Y,X)
  dm.MMUPHin = HMP::Xmcupo.sevsample(Y)
  return(c(q1.MMUPHin$pval[1], q2.MMUPHin$pval[1],dm.MMUPHin$`p value`,ZIGDM.MMUPHin$pval[1],
           q1.MMUPHin$pval[2], q2.MMUPHin$pval[2], dm.MMUPHin.res, ZIGDM.MMUPHin$pval[2]))
}

rd.test <-function(sim.cnt, mod_gdm_lst, depth.pool, taxa.names, lst.num,
                   signal.size, signal.method, type, direct = NULL, seed,
                   simulated = NULL){
  if (is.null(simulated)) {
    simulated <- simulate_crc_data(
      sim.cnt, mod_gdm_lst, depth.pool, taxa.names, lst.num, signal.size,
      signal.method, type, direct, seed
    )
  }
  dat.simu.poured.tmp <- simulated$OTU
  X.lst.tmp <- simulated$X
  Y = do.call(rbind,dat.simu.poured.tmp)
  X = as.matrix(do.call(c, X.lst.tmp))
  Y.rare = Rarefy(Y)$otu.tab.rff
  q1.pour = miLineage::QCAT(Y, X, 1, Tax = NULL,min.depth = 0, n.perm = n.perm)
  q2.pour  =miLineage::QCAT_GEE(Y.rare, X, 1, X, 1, Tax = NULL,
                                min.depth = 0, n.perm = n.perm)
  ZIGDM.pour = miLineage::ZIGDM(Y, as.matrix(X, ncol = 1), NULL, NULL, Tax = NULL, X.index = 1, n.perm = n.perm)
  dm.pour.res = dm1.resample(Y, X, n.perm = 1000)
  Y = split.data.frame(Y,X)
  dm.pour = HMP::Xmcupo.sevsample(Y)
  return(c(q1.pour$pval[1], q2.pour$pval[1], dm.pour$`p value`, ZIGDM.pour$pval[1],
           q1.pour$pval[2], q2.pour$pval[2], dm.pour.res, ZIGDM.pour$pval[2]))
}

run_replicate_tests <- function(sim.cnt, mod_gdm_lst, depth.pool, taxa.names,
                                lst.num, meta.tab.lst, signal.size,
                                signal.method, type, direct = NULL, seed,
                                correlated = FALSE) {
  test.functions <- list(
    meta.qcat.test = meta.qcat.test,
    p.comb.test = p.comb.test,
    cd.test = cd.test,
    rd.test = rd.test
  )
  test.lengths <- c(meta.qcat.test = 9L, p.comb.test = 8L,
                    cd.test = 8L, rd.test = 8L)
  simulated <- tryCatch(
    simulate_crc_data(sim.cnt, mod_gdm_lst, depth.pool, taxa.names,
                      lst.num, signal.size, signal.method, type,
                      direct, seed, correlated = correlated),
    error = function(e) e
  )
  if (inherits(simulated, "error")) {
    error <- paste0("Data generation: ", conditionMessage(simulated))
    return(lapply(test.lengths, function(n) {
      list(value = rep(NA_real_, n), error = error)
    }))
  }

  # The previous implementation regenerated identical data for each method.
  # Restore the post-generation RNG state before each test to preserve its
  # original permutation stream while sharing one generated data set.
  post.generation.rng <- get(".Random.seed", envir = .GlobalEnv)
  common.args <- list(
    sim.cnt = sim.cnt, mod_gdm_lst = mod_gdm_lst, depth.pool = depth.pool,
    taxa.names = taxa.names, lst.num = lst.num, signal.size = signal.size,
    signal.method = signal.method, type = type, direct = direct, seed = seed,
    simulated = simulated
  )
  results <- lapply(names(test.functions), function(test.name) {
    assign(".Random.seed", post.generation.rng, envir = .GlobalEnv)
    args <- common.args
    if (test.name == "cd.test") args$meta.tab.lst <- meta.tab.lst
    do.call(safe_test_call, c(
      list(test.fun = test.functions[[test.name]],
           expected.length = test.lengths[[test.name]],
           test.label = test.name),
      args
    ))
  })
  names(results) <- names(test.functions)
  results
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
      failure.records <- list()
      file.tag <- paste(
        lineage.name, simulation.target, signal.method, scenario.name,
        format(signal.size, scientific = FALSE), sep = "_"
      )

      test.lengths <- c(meta.qcat.test = 9L, p.comb.test = 8L,
                        cd.test = 8L, rd.test = 8L)
      output <- parallel::mclapply(
        seq_len(sim.time), run_replicate_tests,
        mod_gdm_lst = mod_gdm_lst, depth.pool = depth.pool,
        taxa.names = taxa.names, lst.num = lst.num,
        meta.tab.lst = meta.tab.lst, signal.size = signal.size,
        signal.method = signal.method, type = type, direct = direct,
        seed = seed, mc.cores = nthreads, mc.silent = TRUE
      )
      unpacked <- lapply(names(test.lengths), function(test.name) {
        unpack_safe_results(
          lapply(output, `[[`, test.name),
          test.lengths[[test.name]], test.name
        )
      })
      names(unpacked) <- names(test.lengths)
      failure.records <- lapply(unpacked, `[[`, "errors")
      meta <- unpacked$meta.qcat.test$values
      pooled <- unpacked$p.comb.test$values
      corrected <- unpacked$cd.test$values
      raw <- unpacked$rd.test$values
      asy.res <- cbind(
        meta[, 1:3, drop = FALSE], pooled[, 1:4, drop = FALSE],
        corrected[, 1:4, drop = FALSE], raw[, 1:4, drop = FALSE]
      )
      perm.res <- cbind(
        meta[, 4:9, drop = FALSE], pooled[, 5:8, drop = FALSE],
        corrected[, 5:8, drop = FALSE], raw[, 5:8, drop = FALSE]
      )

      colnames(asy.res) <- c(
        "FE-MetaQCAT (asy)", "FE-MetaQCAT* (asy)", "FE-MetaQCAT-O (asy)",
        "P-QCAT1 (asy)", "P-QCAT2 (asy)", "P-DM (asy)", "P-ZIGDM (asy)",
        "CD-QCAT1 (asy)", "CD-QCAT2 (asy)", "CD-DM (asy)", "CD-ZIGDM (asy)",
        "RD-QCAT1 (asy)", "RD-QCAT2 (asy)", "RD-DM (asy)", "RD-ZIGDM (asy)"
      )
      colnames(perm.res) <- c(
        "FE-MetaQCAT (perm)", "FE-MetaQCAT* (perm)", "FE-MetaQCAT-O (perm)",
        "RE-MetaQCAT", "RE-MetaQCAT*", "RE-MetaQCAT-O",
        "P-QCAT1 (perm)", "P-QCAT2 (perm)", "P-DM (perm)", "P-ZIGDM (perm)",
        "CD-QCAT1 (perm)", "CD-QCAT2 (perm)", "CD-DM (perm)", "CD-ZIGDM (perm)",
        "RD-QCAT1 (perm)", "RD-QCAT2 (perm)", "RD-DM (perm)", "RD-ZIGDM (perm)"
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
                file = file.path(output.dir, paste0(file.tag, "_asy_MetaQCAT_org2.csv")))
      write.csv(perm.res, row.names = FALSE,
                file = file.path(output.dir, paste0(file.tag, "_perm_MetaQCAT_org2.csv")))

      failure.log <- do.call(rbind, failure.records)
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
