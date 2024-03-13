# This driver is to run the ca bone model without any drug dosing

library(deSolve)
library(rootSolve)
library(tidyverse)

source("ca.bone.lib.r")

camod <- ca.bone.load.model()
camod <- ca.bone.derive.init(camod)


# set simulation times
times <- seq(0,500,1)

# initial condition
IC = unlist(camod$init[camod$cmt])

# Steady state solution
print('get SS')
ST <- stode(IC, time = 0, func = camod$model, parms = camod$param)
print(ST$y)

print('Running simulation')
## RUN THE MODEL
out <- as.data.frame(
                    lsoda(
                        IC, 
                        times, 
                        camod$model,
                        camod$param,
                        rtol=1e-10,
                        atol=1E-10,
                        ynames=F
                        )
                    )
print('Simulation finished')


## POST PROCESSING
out <- ca.bone.responses(out,camod)
# this outputs the model output to a file so that I can later postprocess
# write.csv(out, file = filename) 

#' Plot
print('Plotting results')
ggplot(out) +
  geom_line(aes(x=time,y=ECCPhos)) 
  labs(x="Time (hours)", y="ECCPhos") 

ggsave("eccphos.png",width = 8, height = 4, dpi = 300)


# Plot Calcium
ggplot(out) +
  geom_line(aes(x=time,y=P)) +
  labs(x="Time (hours)", y="Calcium conc")

  ggsave("caconc.png",width = 8, height = 4, dpi = 300)

# Plot FGF23
# ggplot(out) +
#   geom_line(aes(x=time,y=FGF23)) +
#   labs(x="Time  (hours)", y="FGF23") 
#   ggsave("fgf23.png",width = 8, height = 4, dpi = 300)


# Plot PTH
ggplot(out) +
  geom_line(aes(x=time,y=PTH)) +
  labs(x="Time (hours)", y="PTHconc") 
  ggsave("pth.png",width = 8, height = 4, dpi = 300)

print('done!')


