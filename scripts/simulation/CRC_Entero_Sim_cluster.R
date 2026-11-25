# Enterobacteriaceae: four correlated-sample power settings.
# Run from the repository root. All 12 methods share each generated dataset.
# MetaQCAT: six wild-score bootstrap tests; P/CD/RD: two cluster-QCAT tests.
lineage.name <- "Enterobacteriaceae"
simulation.target <- "power_cluster"
composition.genera <- c(
  "NA Enterobacteriaceae gen. [C Escherichia/Shigella]",
  "544 Citrobacter", "570 Klebsiella", "547 Enterobacter")
delta.proportion <- 1.5
count.probability <- 0.005
source("scripts/simulation/simulation_workflow.R", local = TRUE)
