# This driver is to run the ca bone model without any drug dosing

############# Begin User Input #################
# # where to save the model output
# today <- Sys.Date()
# temp <- "out4daysTeri.csv"
# filename <- paste(today, temp, sep = '_')
############# End of User Input ################

library(deSolve)
library(lattice)
library(reshape)

source("ca.bone.lib.r")

camod <- ca.bone.load.model()
camod <- ca.bone.derive.init(camod)


## EVALUATION TIMES IN HOURS
times <- seq(0,4*24,0.5)


## TERIPARATIDE DOSING EVENTS (TIMES IN HOURS)
teri.times <- seq(12,4*24,24)
teri.dose.mcg <- 20
teri.dose <- teri.dose.mcg*1E6/4117.8



events <- data.frame(
                     var="TERISC",
                     time=teri.times,
                     value=teri.dose,
                     method="add"
                     )


## ADD DOSING EVENTS TO EVALUATION TIMES
times <- sort(unique(c(times,events$time)))


## RUN THE MODEL

out <- as.data.frame(
                     lsoda(
                           unlist(camod$init[camod$cmt]), 
                           times, 
                           camod$model,
                           camod$param,
                           rtol=1e-10,
                           atol=1E-10,
                           ynames=F,
                           events=list(data=events)
                           )
                     )


## POST PROCESSING
out <- ca.bone.responses(out,camod)
# this outputs the model output to a file so that I can later postprocess
write.csv(out, file = filename) 

## Optional plotting routine
## REQUIRES loading of reshape and lattice libraries
print(xyplot(value~time/24|variable,
       data=melt(out,measure.vars=c("BSAP","sCTx","PTHpM","PTHconc",camod$cmt),id.vars="time"),
       type='l',par.strip.text=list(cex=0.8),
       scales=list(y=list(relation='free')),
       xlab="Time (days)",
       ylab="DV"
       )
)


