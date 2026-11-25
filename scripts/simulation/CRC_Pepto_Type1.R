# Peptostreptococcaceae: independent and correlated type I error.
# Run from the repository root; all methods share each generated dataset.
# Independent: 18 methods. Correlated: six MetaQCAT + six native cluster QCAT.
# To run only correlated nulls, set SIM_TYPE1_DESIGNS=correlated.
lineage.name <- "Peptostreptococcaceae"
simulation.target <- "type1"
composition.genera <- c(
  "1257 Peptostreptococcus", "1505657 Intestinibacter",
  "NA unknown Peptostreptococcaceae")
delta.proportion <- 0.75
count.probability <- 0.05
source("scripts/simulation/simulation_workflow.R", local = TRUE)
