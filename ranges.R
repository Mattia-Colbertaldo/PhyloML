library(diversitree)


load("data/FamilyAllTrees.Rdata")
trees <- FamilyAllTrees

ranges <- list()

# loop over the trees
for(i in 1:length(FamilyAllTrees)){
  # get the tree
  tree <- FamilyAllTrees[[i]]
  # get the number of species in the tree
  n_species <- length(tree[[1]]$tip.label)
  # get the number of species in the data
  n_data <- tree[[2]]
  # store the results
  ranges[[i]] <- c(n_species, n_data)
}
ranges

completeTreesIndices <- which(sapply(ranges, function(x) x[1] == x[2]))
completeTrees <- FamilyAllTrees[completeTreesIndices]
completeTreesIndices
plot(completeTrees[[2]][[1]])

# save the complete trees
save(completeTrees, file = "data/completeTrees.Rdata")


t <- completeTrees[[1]][[1]]
plot(t)
lik <- make.bd(t)
p <- starting.point.bd(t)
fit <- find.mle(lik, p)
fit$lnLik
fit$par
coef(fit)

load("data/completeTrees.Rdata")
pars <- list()
threshold <- 0
cat("Running BD model\n")
for (tree in completeTrees) {
    t <- tree[[1]]
    # if the tree is not ultrametric, we need to make it ultrametric
    if (!is.ultrametric(t)) {
        t <- chronos(t)
    }
    tryCatch(
        {
            lik <- make.bd(t)
            p <- starting.point.bd(t)
            fit <- find.mle(lik, p)
            fit$lnLik
            if (fit$lnLik > threshold) {
                pars[[length(pars) + 1]] <- fit$par
            }
        },
        error = function(e) {
            cat("Error in tree\n")
        }
    )

}
pars
lambda_bd_range <- c(min(sapply(pars, function(x) x[1])), max(sapply(pars, function(x) x[1])))
mu_bd_range <- c(min(sapply(pars, function(x) x[2])), max(sapply(pars, function(x) x[2])))

# plot the couples lambda, mu
plot(sapply(pars, function(x) x[1]), sapply(pars, function(x) x[2]), xlab = "lambda", ylab = "mu")
# remove the outlier (the one with the highest lambda)
pars_filtered <- pars[-which(sapply(pars, function(x) x[1]) == max(sapply(pars, function(x) x[1])))]
# plot the couples lambda, mu
plot(sapply(pars, function(x) x[1]), sapply(pars, function(x) x[2]), xlab = "lambda", ylab = "mu")
lambda_bd_range <- c(min(sapply(pars_filtered, function(x) x[1])), max(sapply(pars_filtered, function(x) x[1])))
mu_bd_range <- c(min(sapply(pars_filtered, function(x) x[2])), max(sapply(pars_filtered, function(x) x[2])))
lambda_bd_range
mu_bd_range

# do runif to select a random lambda and mu 50000 times and plot the results
lambda_bd <- runif(5000, lambda_bd_range[1], lambda_bd_range[2])
mu_bd <- runif(5000, mu_bd_range[1], mu_bd_range[2])
plot(lambda_bd, mu_bd, xlab = "lambda", ylab = "mu")

# do the same with a for loop
p <- list()
for (i in 1:5000) {
    p[[i]] <- c(runif(1, lambda_bd_range[1], lambda_bd_range[2]), runif(1, mu_bd_range[1], mu_bd_range[2]))
}
p <- unlist(p)
p <- matrix(p, ncol = 2, byrow = TRUE)
plot(p[, 1], p[, 2], xlab = "lambda", ylab = "mu")

# create a bd tree and save it
t <- completeTrees[[1]][[1]]
plot(t)
ape::write.tree(t, file = paste(file.path("trees", "musse" , "musse_tree_"), 0, ".txt", sep=""))
ape::write.tree(tree_musse_simulations[[i]], file = paste("/trees/musse/musse_tree_", 0, ".nwk", sep=""))
