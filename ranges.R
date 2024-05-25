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

# get the minimum and maximum number of species in the complete trees
min(sapply(completeTrees, function(x) length(x[[1]]$tip.label)))
max(sapply(completeTrees, function(x) length(x[[1]]$tip.label)))

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
cat("Running bisseness model\n")
for (tree in completeTrees) {
    t <- tree[[1]]
    # if the tree is not ultrametric, we need to make it ultrametric
    if (!is.ultrametric(t)) {
        t <- chronos(t)
    }
    tryCatch(
        {
            lik <- make.bisseness(t)
            p <- starting.point.bisseness(t)
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







# for each tree, i want to guess the model that generated it (bd, musse, bisse, geosse, classe)
# i will use the likelihood of the model to do that

# load the trees: for each model there is a file with a lot of trees (bd.nwk, musse.nwk, bisse.nwk, geosse.nwk, classe_2.nwk)
# store in a list the correct model for each tree (the name of the file is the model)
# store in a list the trees
# store in a list the targets (the correct model for each tree)

# load the trees
trees <- list()
models <- c("bd", "musse", "bisse", "geosse", "classe_2")
for(model in models){
    cat("Loading ", model, " trees\n")
    trees[[model]] <- read.tree(paste("trees", paste(model, "_min.nwk", sep=""), sep="/"))
    cat("Loaded ", length(trees[[model]]), " trees\n")
}

# store the correct model for each tree
targets <- c()
for(model in models){
    targets <- c(targets, rep(model, length(trees[[model]])))
}

# store the trees
all_trees <- c()
for(model in models){
    all_trees <- c(all_trees, trees[[model]])
}

# for each tree, i want to guess the model that generated it (bd, musse, bisse, geosse, classe)
# i will use the likelihood of the model to do that (make.bd, make.musse, make.bisse, make.geosse, make.classe)

diversitree_functions <- list(make.bd, make.musse, make.bisse, make.geosse, make.classe)
#diversitree_starting_points <- list(starting.point.bd, starting.point.musse, starting.point.bisse, starting.point.geosse, starting.point.classe)
predictions <- list()
for(i in 1:length(all_trees)){
    tree <- all_trees[[i]]
    likelihoods <- c()
    for(j in 1:length(diversitree_functions)){
        states <- tree$tip.state
        cat("States: ", states, "\n")
        if (model == "classe_2"){
            lik <- make.classe(tree, tree$tip.state, k = 2)
        }
        else if (model == "musse"){
            lik <- make.musse(tree, tree$tip.state, k = 3)
        }
        else{
            lik <- diversitree_functions[[j]](tree, tree$tip.state)
        }
    }
    predictions[[i]] <- likelihoods

}

# get the model with the highest likelihood
predictions <- sapply(predictions, function(x) which.max(x))
predictions <- sapply(predictions, function(x) models[x])
predictions

# get the accuracy
accuracy <- sum(predictions == targets) / length(targets)
accuracy

# save the results
save(predictions, accuracy, file = "data/predictions.Rdata")