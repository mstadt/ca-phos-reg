# Run modeleqns.r with baseline initial conditions
# and parameters

# Run the model
library(deSolve)
library(rootSolve)
library(tidyverse)

# get relevant functions
source("setparams.r")
source("initconds.r")
source("modeqns.r")

# Get initial condition
IC <- unlist(init_conds())

# Get parameters
pars <- set_pars()

# Simulation time
t0 = 0 # start time
tf = 100 # final time
times = seq(t0,tf,1) # hours

# Model
mod = caphos_mod


# Steady state
ST <- stode(IC, time = 0, func = mod, parms = pars)
print(ST$y)

# Model simulation
out <- as.data.frame(lsoda(
                        IC,
                        times,
                        mod, 
                        pars,
                        rtol = 1e-8,
                        atol = 1e-8
                        ))
