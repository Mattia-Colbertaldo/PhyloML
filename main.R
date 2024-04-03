install.packages("phytools")
library(phytools)
install.packages("geiger")
library(geiger)
install.packages("diversitree")
library(diversitree)
install.packages("TreeSim")
library(TreeSim)
install.packages("TreeSimGM")
library(TreeSimGM)
#install.packages("laser") # Removed from CRAN
#library(laser)  # Removed from CRAN
#install.packages("PhyloSim") # Removed from CRAN
#library(PhyloSim) # Removed from CRAN
#install.packages("Rphylip") # Removed from CRAN
#library(Rphylip)  # Removed from CRAN

install.packages('DDD')
library(DDD)


'''
DiversiTree: Comparative Phylogenetic Analyses of Diversification
Simulate: Evolve Birth-Death Trees

Evolves one or more trees under the BiSSE (Binary State Speciation and Extinction), MuSSE
(Multi-State Speciation and Extinction), BiSSE-ness (BiSSE-node enhanced state shift), ClaSSE
(Cladogenetic State change Speciation and Extinction), or GeoSSE (Geographic State Speciation
and Extinction) model, or a simple character independent birth-death model. For the SSE models,
it simultaneously evolves a character that affects speciation and/or extinction, and the tree itself.
Usage
trees(pars, type=c("bisse", "bisseness", "bd", "classe", "geosse",
"musse", "quasse", "yule"), n=1, max.taxa=Inf, max.t=Inf,
include.extinct=FALSE, ...)
tree.bisse(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE,
x0=NA)
tree.musse(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE,
80 simulate
x0=NA)
tree.musse.multitrait(pars, n.trait, depth, max.taxa=Inf, max.t=Inf,
include.extinct=FALSE, x0=NA)
tree.quasse(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE, x0=NA,
single.lineage=TRUE, verbose=FALSE)
tree.bisseness(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE,
x0=NA)
tree.classe(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE,
x0=NA)
tree.geosse(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE,
x0=NA)
tree.bd(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE)
tree.yule(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE)
prune(phy, to.drop=NULL)
'''

# List of packages and corresponding functions
package_functions <- list(
  ape = c("rlineage", "rbdtree", "rphylo"),
  TreeSim = c("sim.bd.taxa", "sim.rateshift.taxa"),
  phytools = c("pbtree_rejection", "pbtree_direct"),
  geiger = c("sim.bdtree"),
  TreeSimGM = c("sim.taxa", "sim.age"),
  diversitree = c("trees", "tree.bisse", "tree.musse", 
                  "tree.quasse", "tree.bisseness", "tree.classe", 
                  "tree.geosse", "tree.bd"), #, "tree.yule"),  <-- yule = bd
  #Laser = c("sim.laser", "lik.laser"),
  #PhyloSim = c("make.bdtree", "sim.bdtree", "make.bdtree"),
  #Rphylip = c("rmtree", "rgenetrator", "rgentree"),
  DDD = c("dd_KI_sim", "dd_MS_sim", "dd_sim", "dd_SR_sim", "td_sim")
)

# call help for each function
for (package in names(package_functions)) {
  functions <- package_functions[[package]]
  library(package, character.only = TRUE)
  for (func in functions) {
    help(func)
  }
}

# Data frame to store execution times
all_times <- data.frame(Package = character(),
                        Function = character(),
                        Time = numeric(),
                        Species_Distribution = numeric(),
                        stringsAsFactors = FALSE)

speciation_rate <- birth <- 0.1
extinction_rate <- death <- 0.05
n <- 50
numbsim <- 1

# Set the number of time steps
num_steps <- t0 <- 50

num_events <- 2

# Generate vectors of speciation rates (lambda) and extinction rates (mu) for each event
lambda <- runif(num_events+1, min = 1, max = 2)  # Example: random uniform speciation rates between 1 and 2 for each time step
# Sample mu from a uniform distribution between 0 and lambda, otherwise sim.rateshift.taxa stucks
mu <- runif(num_events+1, min = 0, max = lambda)

# Generate vector of proportions of species surviving mass extinction event for each event
frac <- c(1, runif(num_events, min = 0.1, max = 1))  # Example: random uniform proportions between 0.1 and 1 for each time step

# Generate vector of mass extinction and rate shift times, starting from 0 and randomly distributed
rateshift_times <- c(0, sort(runif(num_events, min = 0, max = t0)))


lambda
mu
frac
rateshift_times

# Create directories for each package
packages <- names(package_functions)
for (pkg in packages) {
  dir.create(pkg, showWarnings = FALSE)
}

species_distribution <- data.frame()

# Loop over package and functions
for (package in names(package_functions)) {
  functions <- package_functions[[package]]
  library(package, character.only = TRUE)
  for (func in functions) {
    direct <- FALSE
    if(func == "pbtree_direct"){
      direct <- TRUE
    }
    times <- numeric()
    species_numbers <- numeric()
    plots_dir <- file.path(".", package, func)
    dir.create(plots_dir, showWarnings = FALSE)
    help(func)
    for (j in 1:10) { # Run each function 10 times
      tryCatch({
        cat(paste("Running package:", package, ", function:", func), "\n")
        num_ddmodels <- 1
        time <- system.time({
          if (package == "ape") {
            if (func == "rphylo") {
              result <- do.call(func, args = list(n, birth, death, T0 = t0, fossils = FALSE))
            } else {
              result <- do.call(func, args = list(birth, death, t0))
            }
          } else if (package == "TreeSim") {
            if (func == "sim.bd.taxa") {
              result <- do.call(func, args = list(n, numbsim, birth, death))[[1]]
            } else if (func == "sim.rateshift.taxa") {
              result <- do.call(func, args = list(n, numbsim, lambda, mu, frac, rateshift_times))[[1]]
            }
          } else if (package == "geiger") {
            result <- do.call(func, args = list(b = birth, d = death, n = t0, t = t0))
          } else if (package == "phytools") {
            if (func == "pbtree_rejection") {
              func <- "pbtree"
              result <- do.call(func, args = list(b = birth, d = death, n = t0, t = t0, type = "continuous", method="rejection"))
            } else if (func == "pbtree_direct") {
              func <- "pbtree"
              result <- do.call(func, args = list(b = birth, d = death, n = t0, t = t0, type = "continuous", method="direct"))
            }
          } else if (package == "diversitree") {
            # "tree.bisse", "tree.musse", "tree.quasse", "tree.bisseness", "tree.classe", "tree.geosse", "tree.bd", "tree.yule"
            pars <- c(.1, .15, .2, # lambda 1, 2, 3
                      .03, .045, .06, # mu 1, 2, 3
                      .05, 0, # q12, q13
                      .05, .05, # q21, q23
                      0, .05) # q31, q32
            if (func == "trees") {
              ## Simulate a tree under a constant rates birth-death model and look at
              ## the maximum likelihood speciation/extinction parameters:
              result <- trees(c(.1, .03), "bd", max.taxa=25)[[1]]
            } else if (func == "tree.bisse") {
              set.seed(1)
              pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
              result <- tree.bisse(pars, max.t=t0, x0=0)
              h <- history.from.sim.discrete(result, 0:1)
              result_plot <- plot(h, result)[[1]]
              
            } else if (func == "tree.musse") {
              pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
              result <- tree.musse(pars, n, x0=1)
              # h <- history.from.sim.discrete(result, 1:2)
              # result_plot <- plot(h, result)[[1]]
            } else if (func == "tree.quasse") {
              ## Example showing simple integration with two different backends,
              ## plus the splits.
              lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
              mu <- function(x) constant.x(x, 0.03)
              char <- make.brownian.with.drift(0, 0.025)
              result <- tree.quasse(c(lambda, mu, char), max.taxa=n, x0=0,
                                    single.lineage=FALSE, verbose=TRUE)
            } else if (func == "tree.bisseness") {
              pars <- c(0.1, 0.2, 0.03, 0.03, 0, 0, 0.1, 0, 0.1, 0)
              result <- tree.bisseness(pars, max.taxa = n, x0 = 0)
            } else if (func == "tree.classe") {
              pars <- c(0.1, 0.2, 0.03, 0.03, 0, 0, 0.1, 0, 0.1, 0)
              result <- tree.classe(pars, max.taxa = n, max.t = t0, include.extinct = FALSE)
            } else if (func == "tree.geosse") {
              pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
              names(pars) <- diversitree:::default.argnames.geosse()
              result <- tree.geosse(pars, max.t=4, x0=0)
              plot(result)
            } else if (func == "tree.bd") {
              pars <- c(0.1, 0.2, 0.03, 0.03, 0, 0, 0.1, 0, 0.1, 0)
              result <- tree.bd(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE)
            } else if (func == "tree.yule") {
              pars <- c(0.1, 0.2, 0.03, 0.03, 0, 0, 0.1, 0, 0.1, 0)
              result <- tree.yule(pars, max.taxa=Inf, max.t=Inf, include.extinct=FALSE)
              plot(result)
            }
          } else if (package == "DDD") {
            if (func == "dd_KI_sim") {
              result <- do.call(func, args = list(c(0.2, 0.1, 20, 0.1, 0.05, 30, 10), 10))[[1]]
              species_numbers <- c(species_numbers, length(result$tip.label))
            } else if (func == "dd_MS_sim") {
              num_ddmodels <- 4
              for(ddmodel in c(1.3, 2.1, 2.2, 2.3)){
                cat(paste("Running package:", package, ", function:", func, ", ddmodel:", ddmodel), "\n")
                result <- do.call(func, args = list(c(0.2, 0.05, 20, 0.1, 0.05, 40), 100, ddmodel))[[1]]
                species_numbers <- c(species_numbers, length(result$tip.label))/num_ddmodels
              }
            } else if (func == "dd_sim") {
              num_ddmodels <- 10
              for(ddmodel in c(1, 1.3, 2, 2.1, 2.2, 2.3, 3, 4, 4.1, 4.2)){
                cat(paste("Running package:", package, ", function:", func, ", ddmodel:", ddmodel), "\n")
                result <- do.call(func, args = list(c(0.2, 0.1, 20), 20, ddmodel))[[1]]
                species_numbers <- c(species_numbers, length(result$tip.label))/num_ddmodels
              }
            } else if (func == "dd_SR_sim") {
              num_ddmodels <- 9
              for(ddmodel in c(1, 1.3, 2, 2.1, 2.2, 3, 4, 4.1, 4.2)){
                cat(paste("Running package:", package, ", function:", func, ", ddmodel:", ddmodel), "\n")
                result <- do.call(func, args = list(c(0.4, 0.2, 20, 0.2, 0.1, 40, 5), 10, ddmodel))[[1]]
                species_numbers <- c(species_numbers, length(result$tip.label))/num_ddmodels
              }
            } else if (func == "td_sim") {
              num_ddmodels <- 11
              for(ddmodel in c(1, 1.3, 2, 2.1, 2.2, 2.3, 3, 4, 4.1, 4.2, 5)){
                cat(paste("Running package:", package, ", function:", func, ", ddmodel:", ddmodel), "\n")
                result <- do.call(func, args = list(c(0.2, 0.1, 20, 0.1, 0.05, 30, 10), 20, ddmodel))[[1]]
                species_numbers <- c(species_numbers, length(result$tip.label))/num_ddmodels
              }
            }
          }
          
          # Add conditions for other packages
          if (!is.function(result) && package != "diversitree") {
            
            filename <- file.path(plots_dir, paste0("plot_", j, ".png"))
            png(filename)
            dev.off()
          }
        })
        t <- time["elapsed"]/num_ddmodels
        times <- c(times, t)
        
        if(package != "DDD"){
          species_numbers <- c(species_numbers, length(result$tip.label))
        }
        
      }, error = function(e) {
        # Display error message
        cat(paste("Error occurred for package:", package, ", function:", func, "\n"), "\n")
        message(e)
        times <- c(times, NA)
      })
    }
    avg_time <- mean(times, na.rm = TRUE)
    mean_species <- mean(species_numbers, na.rm = TRUE)
    
    # Round time to 3 decimal places
    avg_time <- round(avg_time, 3)
    
    # Round species distribution to 0 decimal places
    mean_species <- round(mean_species, 1) 
    
    if(func == "pbtree"){
      if(direct){
        func <- "pbtree [direct sampling]"
      } else {
        func <- "pbtree [rejection sampling]"
      }
    }
    
    all_times <- rbind(all_times, data.frame(Package = package, Function = func, Time = avg_time, Species_Distribution = mean_species))
    
  }
}

# print the table
cat("Execution time for different functions across packages and number of species distribution:", "\n")
knitr::kable(all_times)




result <- trees(pars, type=c("bisse", "bisseness", "bd", "classe", "geosse",
                             "musse", "quasse", "yule"), n=1, max.taxa=n, max.t=Inf,
                include.extinct=FALSE)
plot(result[[1]])


'''

## Simulate a tree under a constant rates birth-death model and look at
## the maximum likelihood speciation/extinction parameters:
set.seed(1)
phy <- trees(c(.1, .03), "bd", max.taxa=25)[[1]]
lik <- make.bd(phy)
## By default, optimisation gives a lambda close to 0.1 and extremely
## small mu:
fit <- find.mle(lik, c(.1, .03))
coef(fit)
## The above optimisation uses the algorithm \link{nlm} for
## compatibility with ape s \link{birthdeath}. This can be slightly
## improved by using \link{optim} for the optimisation, which allows
## bounds to be specified:
fit.o <- find.mle(lik, c(.1, .03), method="optim", lower=0)
coef(fit.o)
logLik(fit.o) - logLik(fit) # slight improvement
'''


## ggplot
remotes::install_github("YuLab-SMU/ggtree")
install.packages("ggimage")
install.packages("TDbook")
library(ggtree)



speciation_rate <- birth <- 0.1
extinction_rate <- death <- 0.05
n <- 50
numbsim <- 1

# Set the number of time steps
num_steps <- t0 <- 50

num_events <- 2

# Generate vectors of speciation rates (lambda) and extinction rates (mu) for each event
lambda <- runif(num_events+1, min = 1, max = 2)  # Example: random uniform speciation rates between 1 and 2 for each time step
# Sample mu from a uniform distribution between 0 and lambda, otherwise sim.rateshift.taxa stucks
mu <- runif(num_events+1, min = 0, max = lambda)

# Generate vector of proportions of species surviving mass extinction event for each event
frac <- c(1, runif(num_events, min = 0.1, max = 1))  # Example: random uniform proportions between 0.1 and 1 for each time step

# Generate vector of mass extinction and rate shift times, starting from 0 and randomly distributed
rateshift_times <- c(0, sort(runif(num_events, min = 0, max = t0)))

# plot one function for each package

par(mfrow=c(1,2))

# Ape
library("ape")
tree <- rlineage(birth=0.1, death=0.05, Tmax = 50)
plot(tree)
ggtree(tree) + geom_tiplab() + labs(title="Ape:rlineage - constant birth death rate (0.1, 0.05)")
ggtree(tree, layout="circular") + geom_tiplab()



# TreeSim
library("TreeSim") # constant

num_events <- 20

# Generate vectors of speciation rates (lambda) and extinction rates (mu) for each event
lambda <- runif(num_events+1, min = 1, max = 2)  # Example: random uniform speciation rates between 1 and 2 for each time step
# Sample mu from a uniform distribution between 0 and lambda, otherwise sim.rateshift.taxa stucks
mu <- runif(num_events+1, min = 0, max = lambda)

# Generate vector of proportions of species surviving mass extinction event for each event
frac <- c(1, runif(num_events, min = 0.1, max = 0.3))  # Example: random uniform proportions between 0.1 and 1 for each time step

# Generate vector of mass extinction and rate shift times, starting from 0 and randomly distributed
rateshift_times <- c(0, sort(runif(num_events, min = 0, max = t0)))
tree <- sim.rateshift.taxa(30, 1, lambda, mu, frac, rateshift_times)[[1]]
plot(tree)
ggtree(tree) + geom_tiplab() + labs(title="TreeSim::sim.rateshift.taxa - mass extinction events")

# phytools
library("phytools") 
tree <- pbtree(birth=0.1, death=0.05, n=20, method="direct")
ggtree(tree) + geom_tiplab() + labs(title="Phytools::pbtree - constant birth death rate (0.1, 0.05) and taxa-stop criterion (20)")


# geiger
library("geiger")
tree <- sim.bdtree(b=0.1, d=0.05, n=50, t=50)
ggtree(tree)

# TreeSimGM
library("TreeSimGM")
tree <- sim.taxa(numbsim=1, n=50, waitsp=1)
ggtree(tree)

# diversitree
library("diversitree")

pars <- c(.1, .15, .2, # lambda 1, 2, 3
          .03, .045, .06, # mu 1, 2, 3
          .05, 0, # q12, q13
          .05, .05, # q21, q23
          0, .05) # q31, q32
set.seed(2)
phy <- tree.musse(pars, 30, x0=1)
## Extract history from simulated tree and plot
## (colours are 1: black, 2: red, 3: blue)
col <- c("blue", "orange", "red")
h <- history.from.sim.discrete(phy, 1:3)
plot(h, phy, cex=.7, col=col)
## The states are numbered 1:3, rather than 0:1 in bisse.
states <- phy$tip.state
table(states)


# DDD
library("DDD")
result <- dd_KI_sim(c(0.2, 0.1, 20, 0.1, 0.05, 30, 10), 10)
ggtree(result)

# bisse tree
pars <- c(0.1, 0.2, 0.03, 0.06, 0.01, 0.02)
set.seed(2)
phy <- tree.bisse(pars, max.t=50, x0=0)
col <- c("blue", "orange")
h <- history.from.sim.discrete(phy, 0:1)
plot(h, phy, cex = 0.6, col = col)


# diversitree:::default.argnames.classe(2)
#[1] "lambda111" "lambda112" "lambda122" "lambda211" "lambda212" "lambda222" "mu1"      
#[8] "mu2"       "q12"       "q21"
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01, 0.05, 0.05, 0.1, 0.1)
set.seed(2)
phy <- tree.classe(pars, max.taxa=50, max.t=50, include.extinct=FALSE)
plot(phy)




pars <- c(1.5, 0.5, 1.0, 0.7, 0.7, 2.5, 0.5)
names(pars) <- diversitree:::default.argnames.geosse()
set.seed(1)
phy <- tree.geosse(pars, max.t=4, x0=0)
statecols <- c("AB"="violet", "A"="blue", "B"="red")
plot(phy, tip.color=statecols[phy$tip.state+1], cex=0.5)

lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)
set.seed(1)
phy <- tree.quasse(c(lambda, mu, char), max.taxa=15, x0=0,
                   single.lineage=FALSE)
plot(phy)



# bisseness
pars <- c(0.1, 0.2, 0.03, 0.06, 0.01, 0.02, 0.1, 0.2, 0.1, 0.2)
set.seed(2)
phy <- tree.bisseness(pars, max.t=50, x0=0)
col <- c("blue", "orange")
h <- history.from.sim.discrete(phy, 0:1)
plot(h, phy, cex = 0.6, col = col)



library("ggtree")
tree <- rtree(50)
ggtree(tree)
ggtree(tree, layout="roundrect")
ggtree(tree, layout="slanted")
ggtree(tree, layout="ellipse")
ggtree(tree, layout="circular")
ggtree(tree, layout="fan", open.angle=120)
ggtree(tree, layout="equal_angle")
ggtree(tree, layout="daylight")
ggtree(tree, branch.length='none')
ggtree(tree, layout="ellipse", branch.length="none")
ggtree(tree, branch.length='none', layout='circular')
ggtree(tree, layout="daylight", branch.length = 'none')


# Matrix Lambda
Lambda <- matrix(c(lambda11, lambda12, lambda21, lambda22), nrow=2, ncol=2, dimnames=list(c("1", "2"), c("1", "2")))

# Matrix Q
Q <- matrix(c(q11, q12, q21, q22), nrow=2, ncol=2, dimnames=list(c("1", "2"), c("1", "2")))

# Vector mu
mu <- c(mu1, mu2)

result <- tree.classe(object, Lambda=Lambda, Q=Q, mu=mu)

