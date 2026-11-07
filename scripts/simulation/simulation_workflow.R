# Shared, directly executable workflow used by the four short family entrypoints.
# Run from the repository root. No input or result files are distributed here.
source("R/project_paths.R", local = TRUE)
package.dir <- file.path(project.dir, "packages", "MetaQCAT")
required <- c("pkgload", "parallelly", "HMP", "miLineage", "MMUPHin")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
pkgload::load_all(package.dir, quiet = TRUE)
if (as.character(utils::packageVersion("MetaQCAT")) != "0.2.1" ||
    !"n.boot" %in% names(formals(MetaQCAT::MetaQCAT)))
  stop("Use the bundled MetaQCAT 0.2.1 wild-score bootstrap package.")

seed <- 2333L
sim.time <- as.integer(Sys.getenv("SIM_REPS",
  if (simulation.target == "power") "1000" else "5000"))
n.boot <- as.integer(Sys.getenv("SIM_BOOTS", "300"))
n.perm <- as.integer(Sys.getenv("SIM_PERMS", "300"))
n.perm.rd.dm <- as.integer(Sys.getenv("SIM_RD_DM_PERMS", "1000"))
nthreads <- as.integer(Sys.getenv("SIM_THREADS", "7"))
stopifnot(!anyNA(c(sim.time, n.boot, n.perm, n.perm.rd.dm, nthreads)),
          min(sim.time, n.boot, n.perm, n.perm.rd.dm, nthreads) > 0L)
if (.Platform$OS.type == "windows") nthreads <- 1L
nthreads <- min(nthreads, parallelly::availableCores())
min.prev <- 0.1
alpha <- 0.05
meta.methods <- c("FE-MetaQCAT", "FE-MetaQCAT*", "FE-MetaQCAT-O",
                  "RE-MetaQCAT", "RE-MetaQCAT*", "RE-MetaQCAT-O")
baseline.methods <- unlist(lapply(c("P", "CD", "RD"), function(strategy)
  paste0(strategy, "-", c("QCAT1", "QCAT2", "DM", "ZIGDM"), " (perm)")))
methods <- c(meta.methods, baseline.methods)

input.files <- file.path(project.dir, "data", "simu", c(
  paste0("mod_gdm_lst_", lineage.name, ".RData"), "meta.rda", "count.genus.rda"))
if (any(!file.exists(input.files)))
  stop("Restore the external simulation inputs listed in docs/SIMULATION_WORKFLOWS.md.")
for (file in input.files) load(file)
study.names <- c("AUS", "CHI", "FRA", "GER", "USA")
lst.num <- as.integer(table(factor(meta$Country, levels = study.names)))
n.taxa <- length(mod_gdm_lst[[1L]]$col.order)
taxa.names <- paste0("taxa", seq_len(n.taxa))
stopifnot(identical(lst.num, c(109L, 127L, 114L, 120L, 104L)),
          length(mod_gdm_lst) == length(study.names),
          n.taxa == length(composition.genera),
          all(composition.genera %in% colnames(count.genus)),
          identical(rownames(count.genus), make.names(meta$Sample_ID)))
observed.depth <- rowSums(count.genus[, composition.genera, drop = FALSE])
depth.pool <- lapply(study.names, function(study)
  as.integer(observed.depth[meta$Country == study]))
stopifnot(all(lengths(depth.pool) == lst.num))

source(file.path(project.dir, "scripts", "simulation", "simulation_helpers.R"),
       local = TRUE)
if (simulation.target == "power") {
  settings <- expand.grid(scenario = c("homogeneous", "one.flipped.study"),
    signal.method = c("proportion", "count"), stringsAsFactors = FALSE)
  settings$signal.size <- ifelse(settings$signal.method == "proportion",
                                 delta.proportion, count.probability)
  settings$correlated <- FALSE
} else {
  settings <- data.frame(scenario = c("independent", "correlated"),
    signal.method = "none", signal.size = 0,
    correlated = c(FALSE, TRUE), stringsAsFactors = FALSE)
}
output.dir <- Sys.getenv("SIM_OUTPUT", file.path(
  project.dir, "outcome", "simulation", simulation.target, lineage.name))

# Set options(crc.simulation.source.only = TRUE) to inspect the loaded functions
# without launching the run. The normal scripts execute the loops below.
if (!isTRUE(getOption("crc.simulation.source.only", FALSE))) {
  dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  source.files <- c(list.files(file.path(project.dir, "scripts", "simulation"),
                               pattern = "[.]R$", full.names = TRUE),
    list.files(file.path(package.dir, "R"), full.names = TRUE),
    list.files(file.path(package.dir, "src"),
               pattern = "[.](cpp|h)$", full.names = TRUE),
    file.path(package.dir, "DESCRIPTION"))
  configuration <- list(
    lineage = lineage.name, target = simulation.target, settings = settings,
    methods = methods, seed = seed, replicates = sim.time, n.boot = n.boot,
    n.perm = n.perm, n.perm.rd.dm = n.perm.rd.dm, threads = nthreads,
    min.prev = min.prev, alpha = alpha, planned.sample.sizes = lst.num,
    composition.genera = composition.genera, study.names = study.names,
    mod_gdm_lst = mod_gdm_lst, depth.pool = depth.pool,
    package.version = as.character(utils::packageVersion("MetaQCAT")),
    dependency.versions = vapply(required, function(package)
      as.character(utils::packageVersion(package)), character(1)),
    input.md5 = tools::md5sum(input.files),
    source.md5 = tools::md5sum(source.files), rng = RNGkind(),
    session = sessionInfo())
  saveRDS(configuration, file.path(output.dir, "configuration.rds"))
  all.rows <- summary.rows <- list()

  for (setting in seq_len(nrow(settings))) {
    scenario <- settings$scenario[setting]
    signal.method <- settings$signal.method[setting]
    signal.size <- settings$signal.size[setting]
    correlated <- settings$correlated[setting]
    type <- if (scenario == "one.flipped.study") "random" else "fixed"
    direct <- if (type == "random") "both" else NULL
    # At the null, either signal mechanism with size zero produces no signal.
    generation.signal <- if (signal.method == "none") "proportion" else signal.method
    tag <- paste(lineage.name, simulation.target, scenario, signal.method, sep = "_")
    message(tag, ": ", sim.time, " replicates; six bootstrap and 12 comparator methods")
    started <- proc.time()[["elapsed"]]
    results <- parallel::mclapply(seq_len(sim.time), run_simulation_replicate,
      signal.size = signal.size, signal.method = generation.signal,
      type = type, direct = direct, correlated = correlated,
      mc.cores = nthreads, mc.set.seed = FALSE, mc.silent = TRUE)
    wall.seconds <- proc.time()[["elapsed"]] - started
    if (!all(vapply(results, function(x)
      is.list(x) && is.data.frame(x$rows) && nrow(x$rows) == length(methods),
      logical(1))))
      stop("A worker terminated unexpectedly; no incomplete setting is summarized: ", tag)
    rows <- do.call(rbind, lapply(results, "[[", "rows"))
    rows$lineage <- lineage.name
    rows$simulation_target <- simulation.target
    rows$scenario <- scenario
    rows$signal.method <- signal.method
    rows$signal.size <- signal.size
    rows$setting_wall_time_sec <- wall.seconds

    # Failures remain NA in replicate tables, and count as nonrejections here.
    # Wilson intervals quantify finite-replicate Monte Carlo uncertainty.
    table <- do.call(rbind, lapply(methods, function(method) {
      r <- rows[rows$method == method, ]
      valid <- !r$failed & is.finite(r$resampling_p)
      valid.asy <- is.finite(r$asymptotic_p) & r$asymptotic_p >= 0 & r$asymptotic_p <= 1
      attempted <- nrow(r)
      rejected <- sum(r$resampling_p[valid] < alpha)
      rate <- rejected / attempted
      z <- qnorm(.975)
      den <- 1 + z^2 / attempted
      center <- (rate + z^2 / (2 * attempted)) / den
      half <- z * sqrt(rate * (1 - rate) / attempted +
                         z^2 / (4 * attempted^2)) / den
      data.frame(
        lineage = lineage.name, simulation_target = simulation.target,
        scenario = scenario, signal.method = signal.method, signal.size = signal.size,
        method = method, calibration = r$calibration[1L],
        rejection_rate = rate, n_rejections = rejected,
        mc_ci_95_lower = max(0, center - half), mc_ci_95_upper = min(1, center + half),
        failures = sum(!valid), failure_rate = mean(!valid),
        n_replicates = attempted, n_valid = sum(valid),
        rejection_rate_valid_only = if (any(valid)) rejected / sum(valid) else NA_real_,
        asymptotic_rejection_rate = if (any(valid.asy))
          sum(r$asymptotic_p[valid.asy] < alpha) / attempted else NA_real_,
        asymptotic_n_valid = sum(valid.asy),
        partial_study_failures = sum(r$partial_study_failure),
        n_boot_requested = r$n_boot_requested[1L],
        n_boot_completed_min = min(r$n_boot_completed),
        n_perm_configured = r$n_perm_configured[1L],
        n_perm_completed_min = if (all(is.na(r$n_perm_completed))) NA_integer_ else
          min(r$n_perm_completed, na.rm = TRUE),
        mean_method_runtime_sec = mean(r$method_runtime_seconds),
        mean_replicate_runtime_sec = mean(r$replicate_runtime_seconds),
        setting_wall_time_sec = wall.seconds, stringsAsFactors = FALSE)
    }))
    write.csv(rows, file.path(output.dir, paste0(tag, "_pvalues.csv")), row.names = FALSE)
    write.csv(table, file.path(output.dir, paste0(tag, "_summary.csv")), row.names = FALSE)
    write.csv(rows[rows$failed | rows$partial_study_failure, ],
              file.path(output.dir, paste0(tag, "_failures.csv")), row.names = FALSE)
    saveRDS(list(configuration = configuration, results = results, rows = rows),
            file.path(output.dir, paste0(tag, "_results.rds")))
    all.rows[[tag]] <- rows
    summary.rows[[tag]] <- table
    prefix <- paste0(lineage.name, "_", simulation.target)
    write.csv(do.call(rbind, all.rows), file.path(output.dir,
      paste0(prefix, "_pvalues.csv")), row.names = FALSE)
    write.csv(do.call(rbind, summary.rows), file.path(output.dir,
      paste0(prefix, "_summary.csv")), row.names = FALSE)
    saveRDS(list(configuration = configuration, pvalues = do.call(rbind, all.rows),
      summary = do.call(rbind, summary.rows)), file.path(output.dir,
        paste0(prefix, "_results.rds")))
    print(table[, c("method", "rejection_rate", "failures")], row.names = FALSE)
  }
}
