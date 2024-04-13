# plot_trees.R

# Load diversitree package
library(diversitree)
library(ape)
# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

set.seed(2)


# Convert command-line arguments to numeric vectors
y_test <- as.numeric(unlist(strsplit(args[1], ",")))
predicted_test <- as.numeric(unlist(strsplit(args[2], ",")))
max_taxa <- as.numeric(args[3])

# Define function to plot trees
plot_trees <- function(y_test, predicted_test, max_taxa) {
  
  t <- tree.musse(pars=y_test, x0=1, max.taxa=max_taxa)
  ape::write.tree(t, file = "target.nwk", append = TRUE)
  
  t <- tree.musse(pars=predicted_test, x0=1, max.taxa=max_taxa)
  ape::write.tree(t, file = "predicted.nwk", append = TRUE)
  
}

plot_trees(y_test, predicted_test, max_taxa)
