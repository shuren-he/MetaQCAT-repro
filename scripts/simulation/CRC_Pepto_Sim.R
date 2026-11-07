# Peptostreptococcaceae: four power settings.
# Run from the repository root. All 18 methods share each generated dataset.
lineage.name <- "Peptostreptococcaceae"
simulation.target <- "power"
composition.genera <- c(
  "1257 Peptostreptococcus", "1505657 Intestinibacter",
  "NA unknown Peptostreptococcaceae")
delta.proportion <- 0.75
count.probability <- 0.05
source("scripts/simulation/simulation_workflow.R", local = TRUE)
