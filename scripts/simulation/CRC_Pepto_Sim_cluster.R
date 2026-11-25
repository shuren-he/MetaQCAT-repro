# Peptostreptococcaceae: four correlated-sample power settings.
# Run from the repository root. All 12 methods share each generated dataset.
# MetaQCAT: six wild-score bootstrap tests; P/CD/RD: two cluster-QCAT tests.
lineage.name <- "Peptostreptococcaceae"
simulation.target <- "power_cluster"
composition.genera <- c(
  "1257 Peptostreptococcus", "1505657 Intestinibacter",
  "NA unknown Peptostreptococcaceae")
delta.proportion <- 0.75
count.probability <- 0.05
source("scripts/simulation/simulation_workflow.R", local = TRUE)
