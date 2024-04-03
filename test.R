# Set parameters for the constant rate model
speciation_rate <- birth <- 0.1
extinction_rate <- death <- 0.05
n <- 50
t0 <- 50
# Total rate for coalescent model
total_rate <- speciation_rate - extinction_rate


## APE
## Tree Simulation Under the Time-Dependent Birth–Death Models

if (!requireNamespace("ape", quietly = TRUE)) {
  install.packages("ape")
}
library(ape)

ape_times <- data.frame(Function = c("rlineage", "rbdtree", "rphylo"),
                        Time = numeric(3))

par(mfrow = c(1, 3))
for (i in 1:3) {
  times <- c()
  for (j in 1:10) { # Run each function 10 times
    time <- system.time({
      plot_result <- switch(i,
                            rlineage = rlineage(birth, death, Tmax = 50, eps = 1e-6),
                            rbdtree = rbdtree(birth, death, Tmax = 50, eps = 1e-6),
                            rphylo = rphylo(n, birth, death, T0 = 50, fossils = FALSE, eps = 1e-06)
      )
    })
    times <- c(times, time[3])
  }
  ape_times[i, "Time"] <- mean(times)
}
colnames(ape_times) <- c("Function", "Time")

## Generate Random Trees

a <- rtree(n, rooted = TRUE, tip.label = NULL, br = runif, equiprob = FALSE)
b <- rtopology(n, rooted = FALSE, tip.label = NULL, br = runif)
c <- rcoal(n, tip.label = NULL, br = "coalescent")

random_trees_times <- data.frame(Function = c("rtree", "rtopology", "rcoal"),
                                 Time = numeric(3))

par(mfrow = c(1, 3))
for (i in 1:3) {
  times <- c()
  for (j in 1:10) { # Run each function 10 times
    time <- system.time({
      plot_result <- switch(i,
                            rtree = plot(a),
                            rtopology = plot(b),
                            rcoal = plot(c)
      )
    })
    times <- c(times, time[3])
  }
  random_trees_times[i, "Time"] <- mean(times)
}
colnames(random_trees_times) <- c("Function", "Time")

## PHYTOOLS and GEIGER
## Simulate pure-birth or birth-death stochastic tree or trees

install.packages("phytools")
library(phytools)
install.packages("geiger")
library(geiger)
install.packages("diversitree")
library(diversitree)
install.packages("TreeSim")
library(TreeSim)

phytools_geiger_times <- data.frame(Function = c("pbtree", "sim.bdtree rejection sampling", "sim.bdtree direct sampling"),
                                    Time = numeric(3))

par(mfrow = c(1, 3))
for (i in 1:3) {
  times <- c()
  for (j in 1:10) { # Run each function 10 times
    time <- system.time({
      plot_result <- switch(i,
                            pbtree = pbtree(b = birth, d = death, n = n, t = t0, type = "continuous"),
                            sim_bdtree_rejection = pbtree(b = birth, d = death, n = n, t = t0, type = "continuous", method = "direct"),
                            sim_bdtree_direct = sim.bdtree(b = birth, d = death, n = n, t = t0, extinct = FALSE, stop = c("taxa", "time"))
      )
    })
    times <- c(times, time[3])
  }
  phytools_geiger_times[i, "Time"] <- mean(times)
}
colnames(phytools_geiger_times) <- c("Function", "Time")

# Print the tables
library(knitr)
print("APE Functions:")
knitr::kable(ape_times)

print("Random Trees Functions:")
knitr::kable(random_trees_times)

print("Phytools and Geiger Functions:")
knitr::kable(phytools_geiger_times)














