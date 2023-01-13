#!/usr/bin/R

## ==============================================================================
## author          :Ghislain Vieilledent
## email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
## web             :https://ghislainv.github.io
## license         :GPLv3
## ==============================================================================

#= Import of R ayer functions
source("./ayer.dbh.R")
source("./ayer.dbh.comp.R")

#= Library
library(MCMCpack) # For Bayesian estimation of mortality rate in each bin for function ayer.dbh

#= Import of Katalin data
dat <- read.table(file="dat.txt",header=TRUE,sep="\t")
nobs <- nrow(dat)

# Call of the ayer.dbh function on test data
model.dbh <- ayer.dbh(mort=as.numeric(dat$mort),dbh=dat$dbh,ti=rep(1,nobs),init.bin=5)
model.dbh
model.dbh$deviance

# Call of the ayer.dbh.comp function on test data
model.dbh.comp <- ayer.dbh.comp(mort=as.numeric(dat$mort),dbh=dat$dbh,comp=dat$bigBA2,
                                ti=rep(1,nobs),
                                dbh0=model.dbh$dbh0,
                                init.bin.dbh=5,
                                init.bin.comp=10000)

model.dbh.comp$deviance
pdf("mortality-dbh-comp.pdf")
image(z=as.matrix(model.dbh.comp$model),col=gray((20:1)/32),
      x=seq(from=5.5,to=117.5,by=5),
      y=seq(from=0.5e+04,to=1.05e+05,by=1e+04),
      xlab="DBH (cm)",ylab="Competition")
dev.off()

# End