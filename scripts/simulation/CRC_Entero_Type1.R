# Enterobacteriaceae: independent and correlated type I error.
# Run from the repository root. All 18 methods share each generated dataset.
lineage.name <- "Enterobacteriaceae"
simulation.target <- "type1"
composition.genera <- c(
  "NA Enterobacteriaceae gen. [C Escherichia/Shigella]",
  "544 Citrobacter", "570 Klebsiella", "547 Enterobacter")
delta.proportion <- 1.5
count.probability <- 0.005
source("scripts/simulation/simulation_workflow.R", local = TRUE)
