##
## ayer.dbh.R
##
## Author: Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
##
## Date: 24-11-2011
##

## ayer.dbh() function fits a U-shape semi-parametric mortality model
## as a function of the tree dbh. First, the function estimates the
## dbh at which tree mortality rate is minimal (dbh0). Second, the
## function estimates the dbh bins and the corresponding annual
## mortality rates to obtain a decreasing mortality with dbh before
## dbh0 and an increasing mortality with dbh beyond dbh0.

## Parameters:
##
## mort: vector of 0 or 1 for tree bernoulli mortality process (1=dead, 0=alive) between census c and c+1
## dbh: vector for tree diameter at census c
## ti: time interval in year between census c and c+1
## init.bin: initial bin width in cm (5 cm is an appropriate value)

## Value:
##
## dbh0: dbh for minimal tree mortality rate. There might be
## equivalent models while looking at the deviance and thus
## potentially several dbh0. In such cases, the model computes the
## mean of the dbh0 (dbh0.bar) and select the closest dbh value from
## the initial bins such that dbh0.final < dbh0.bar
##
## deviance: model deviance
##
## model: matrix describing the model
##
## Moreover, the function plots in the current working directory (1)
## the relationship between the deviance and the value of dbh0 and (2)
## the semi-parametric relationship between mortality rate and tree
## dbh

## References:
##
## Ayer, M.; Brunk, H. D.; Ewing, G. M.; Reid, W. T. & Silverman, E. An
## empirical distribution function for sampling with incomplete
## information The Annals of Mathematical Statistics, 1955, 26, 641-647
##  
## Vieilledent, G.; Courbaud, B.; Kunstler, G.; Dhote, J. F. & Clark,
## J. S. Biases in the estimation of size-dependent mortality models:
## advantages of a semiparametric approach Canadian Journal Of Forest
## Research-Revue Canadienne De Recherche Forestiere, 2009, 39, 1430-1443

#= inv.logit function
inv.logit <- function(x, min=0, max=1) {
  p <- exp(x)/(1+exp(x))
  p <- ifelse(is.na(p) & !is.na(x), 1, p) # fix problems with +Inf
  return(p*(max-min)+min)
}

#= Ayer function
ayer.dbh <- function (mort,dbh,ti,init.bin) {
 
#= Recompose data
data <- data.frame(dbh=dbh,ti=ti,mort=mort)
  
#= Compute minimal and maximal dbh for dbh interval
dbh.min <- floor(min(data$dbh)/init.bin)*init.bin
dbh.max <- ceiling(max(data$dbh)/init.bin)*init.bin

#= dbh breaks
data$dbh.breaks <- cut(data$dbh,breaks=seq(dbh.min,dbh.max,by=init.bin),right=FALSE)
Levels.dbh.breaks <- levels(as.factor(data$dbh.breaks))
data$dbh.breaks <- as.character(data$dbh.breaks)
n.dbh.breaks <- length(Levels.dbh.breaks)
theta.MLE <- rep(NA,n.dbh.breaks)

#= Matrix which ranges dbh
Mat.Range.dbh <- as.data.frame(matrix(NA,nrow=length(Levels.dbh.breaks),ncol=6))
for (i in 1:n.dbh.breaks) {
  Mat.Range.dbh[i,2] <- strsplit(as.character(Levels.dbh.breaks),split=c(","))[[i]][1]
  Mat.Range.dbh[i,3] <- strsplit(as.character(Levels.dbh.breaks),split=c(","))[[i]][2]
}
Mat.Range.dbh[,1]<-paste(Mat.Range.dbh[,2],Mat.Range.dbh[,3],sep=",")
colnames(Mat.Range.dbh)<-c("dbh.class","BeginClass","EndClass","theta.MLE","n.Alive","n.Dead")

#======================
#= Likelihood function
Binom.Likelihood <- function (theta,Y.dead,Y.alive) {
  logvrais <- sum(log(1-(1-theta)^Y.dead))+sum(log((1-theta)^Y.alive))
  return(-logvrais)
}

#===========================================
#= Mortality rate estimate for each dbh class
for (i in 1:n.dbh.breaks) {
  # Variables for considered dbh class
  data.sub <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
  Y.dead <- data.sub$ti[data.sub$mort==1]
  Y.alive <- data.sub$ti[data.sub$mort==0]
  # Maximizing loglikelihood
  Min.Minus.LL <- optim(par=0.01,
                        fn=Binom.Likelihood,
                        method="L-BFGS-B",
                        lower=0.000001,upper=0.999999,
                        Y.dead=Y.dead,Y.alive=Y.alive)
  Mat.Range.dbh$theta.MLE[Mat.Range.dbh$dbh.class==Levels.dbh.breaks[i]] <- max(0,Min.Minus.LL$par)
  # Dead and Alive
  Mat.Range.dbh[i,5] <- length(Y.alive)
  Mat.Range.dbh[i,6] <- length(Y.dead)
}

cat("\nComputing the model deviance for each value of dbh0 in the dbh range...\n")

#==========================================================
#= Searching for the DBH with minimal mortality rate (DBH0)
test.dbh0 <- seq(from=dbh.min,to=dbh.max,by=init.bin)
Deviance <- rep(NA,length(test.dbh0))
Stop.Values <- seq(from=3,to=length(Mat.Range.dbh$dbh.class)-1)
    
# Loops on DBH0 values
for (jj in 1:length(Stop.Values)) {

  Stop <- Stop.Values[jj]
  # Decrease for dbh<Stop
  Mat.Range.dbh.2 <- Mat.Range.dbh[-seq(from=Stop,to=length(Mat.Range.dbh$dbh.class)),]
  # Order will check for monotonicity
  Order <- vector()
  for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
    Order[j] <- Mat.Range.dbh.2$theta.MLE[j]<=Mat.Range.dbh.2$theta.MLE[j+1]|
              Mat.Range.dbh.2$theta.MLE[j]==0|
              Mat.Range.dbh.2$theta.MLE[j+1]==0
  }
  # Loops
  while (sum(Order)!=0) {
    for (k in 1:(length(Mat.Range.dbh.2[,1])-1)) {
        # We check for unordered bins between k and k+1
        if (Mat.Range.dbh.2$theta.MLE[k]<=Mat.Range.dbh.2$theta.MLE[k+1]|
            Mat.Range.dbh.2$theta.MLE[k]==0|
            Mat.Range.dbh.2$theta.MLE[k+1]==0) {
          # We replace the values of the bin in data
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k+1]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          # We modify Mat.Range.dbh.2 in consequence
          Mat.Range.dbh.2$EndClass[k] <- Mat.Range.dbh.2$EndClass[k+1]
          Mat.Range.dbh.2$dbh.class[k] <- paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k],sep=",")
          Mat.Range.dbh.2$n.Alive[k] <- Mat.Range.dbh.2$n.Alive[k]+Mat.Range.dbh.2$n.Alive[k+1]
          Mat.Range.dbh.2$n.Dead[k] <- Mat.Range.dbh.2$n.Dead[k]+Mat.Range.dbh.2$n.Dead[k+1]  
          Mat.Range.dbh.2 <- Mat.Range.dbh.2[-(k+1),]
          ############################################
          ### New estimate of theta for each dbh.class of Mat.Range.dbh.2
          Levels.dbh.breaks <- sort(unique(Mat.Range.dbh.2$dbh.class))
          n.dbh.breaks <- length(Levels.dbh.breaks)
          # Mortality rate estimate for each dbh class
          for (i in 1:n.dbh.breaks) {
            # Variables for considered dbh class
            data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
            Y.dead <- data.temp$ti[data.temp$mort==1]
            Y.alive <- data.temp$ti[data.temp$mort==0]
            if (length(data.temp$ti[data.temp$mort==1])==0) {
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- 0
            }
            if (length(data.temp$ti[data.temp$mort==0])>0) { 
              # Minimizing -loglikelihood
              Min.Minus.LL <- optim(par=0.01,
                                    fn=Binom.Likelihood,
                                    method="L-BFGS-B",
                                    lower=0.000001,upper=0.999999,
                                    Y.dead=Y.dead,Y.alive=Y.alive)
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- max(0,Min.Minus.LL$par)
            }
          }
          ############################
          # We actualize Order
          Order <- vector()
          if (length(Mat.Range.dbh.2$theta.MLE)>1) {
            for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
            Order[j] <- Mat.Range.dbh.2$theta.MLE[j]<=Mat.Range.dbh.2$theta.MLE[j+1]|
                      Mat.Range.dbh.2$theta.MLE[j]==0|
                      Mat.Range.dbh.2$theta.MLE[j+1]==0
            }
          }
          else { Order <- FALSE }
          break
        }
    }
  }
  
  MatRange.decrease <- Mat.Range.dbh.2
  
  # Increase for dbh>Stop
  Mat.Range.dbh.2 <- Mat.Range.dbh[-seq(from=1,to=Stop-1),]
  # Order will check for monotonicity
  Order <- vector()
  for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
    Order[j] <- Mat.Range.dbh.2$theta.MLE[j]>=Mat.Range.dbh.2$theta.MLE[j+1]|
              Mat.Range.dbh.2$theta.MLE[j]==0|
              Mat.Range.dbh.2$theta.MLE[j+1]==0
  }
  # Loops
  while (sum(Order)!=0) {
    for (k in 1:(length(Mat.Range.dbh.2[,1])-1)) {
        # We check for unordered bins between k and k+1
        if (Mat.Range.dbh.2$theta.MLE[k]>=Mat.Range.dbh.2$theta.MLE[k+1]|
            Mat.Range.dbh.2$theta.MLE[k]==0|
            Mat.Range.dbh.2$theta.MLE[k+1]==0) {
          # We replace the values of the bin in data
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k+1]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          # We modify Mat.Range.dbh.2 in consequence
          Mat.Range.dbh.2$EndClass[k] <- Mat.Range.dbh.2$EndClass[k+1]
          Mat.Range.dbh.2$dbh.class[k] <- paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k],sep=",")
          Mat.Range.dbh.2$n.Alive[k] <- Mat.Range.dbh.2$n.Alive[k]+Mat.Range.dbh.2$n.Alive[k+1]
          Mat.Range.dbh.2$n.Dead[k] <- Mat.Range.dbh.2$n.Dead[k]+Mat.Range.dbh.2$n.Dead[k+1]  
          Mat.Range.dbh.2 <- Mat.Range.dbh.2[-(k+1),]
          ############################################
          ### New estimate of theta for each dbh class
          Levels.dbh.breaks <- sort(unique(Mat.Range.dbh.2$dbh.class))
          n.dbh.breaks <- length(Levels.dbh.breaks)
          # Mortality rate estimate for each dbh class
          for (i in 1:n.dbh.breaks) {
            # Variables for considered dbhclass
            data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
            Y.dead <- data.temp$ti[data.temp$mort==1]
            Y.alive <- data.temp$ti[data.temp$mort==0]
            if (length(data.temp$ti[data.temp$mort==1])==0) {
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- 0
            }
            if (length(data.temp$ti[data.temp$mort==0])>0) { 
              # Minimizing -loglikelihood
              Min.Minus.LL <- optim(par=0.01,
                                    fn=Binom.Likelihood,
                                    method="L-BFGS-B",
                                    lower=0.000001,upper=0.999999,
                                    Y.dead=Y.dead,Y.alive=Y.alive) 
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- max(0,Min.Minus.LL$par)
            }
          }
          ############################
          # We actualize Order
          Order <- vector()
          if (length(Mat.Range.dbh.2$theta.MLE)>1) {
            for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
              Order[j] <- Mat.Range.dbh.2$theta.MLE[j]>=Mat.Range.dbh.2$theta.MLE[j+1]|
                        Mat.Range.dbh.2$theta.MLE[j]==0|
                        Mat.Range.dbh.2$theta.MLE[j+1]==0
            }
          }
          else { Order <- FALSE }
          break
        }
      }
  }
  
  MatRange.increase <- Mat.Range.dbh.2
  MatRange.all <- rbind(MatRange.decrease,MatRange.increase)
  
  ######## Deviance calcul #########
  Levels.dbh.breaks <- sort(unique(MatRange.all$dbh.class))
  n.dbh.breaks <- length(Levels.dbh.breaks)
  MatRange.all$LLikelihood <- rep(NA,n.dbh.breaks)
  # Mortality rate estimate for each dbh class
  for (i in 1:n.dbh.breaks) {
    data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
    Y.dead <- data.temp$ti[data.temp$mort==1]
    Y.alive <- data.temp$ti[data.temp$mort==0]
    MatRange.all$LLikelihood[MatRange.all$dbh.class==Levels.dbh.breaks[i]] <- 
      sum(log(1-(1-MatRange.all$theta.MLE[MatRange.all$dbh.class==Levels.dbh.breaks[i]])^Y.dead))+
        sum(log((1-MatRange.all$theta.MLE[MatRange.all$dbh.class==Levels.dbh.breaks[i]])^Y.alive))
  }
  Deviance[jj+2] <- -2*sum(MatRange.all$LLikelihood)
  
  #### data is reinitiate to its initial value
  data$dbh.breaks <- cut(data$dbh,breaks=seq(dbh.min,dbh.max,by=init.bin),right=FALSE)
  data$dbh.breaks <- as.character(data$dbh.breaks)


  ## Message
  cat("*")
}    

#===========================================
#= Particular cases

cat("\nParticular cases...")
    
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Stop=1: Only increasing mortality rate
Stop <- 1

# Increase for dbh>Stop
  Mat.Range.dbh.2 <- Mat.Range.dbh
  # Order will check for monotonicity
  Order <- vector()
  for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
      Order[j] <- Mat.Range.dbh.2$theta.MLE[j]>=Mat.Range.dbh.2$theta.MLE[j+1]|
      Mat.Range.dbh.2$theta.MLE[j]==0|
      Mat.Range.dbh.2$theta.MLE[j+1]==0
  }
  # Loops
  while (sum(Order)!=0) {
    for (k in 1:(length(Mat.Range.dbh.2[,1])-1)) {
        # We check for unordered bins between k and k+1
        if (Mat.Range.dbh.2$theta.MLE[k]>=Mat.Range.dbh.2$theta.MLE[k+1]|
            Mat.Range.dbh.2$theta.MLE[k]==0|
            Mat.Range.dbh.2$theta.MLE[k+1]==0) {
          # We replace the values of the bin in data
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k+1]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          # We modify Mat.Range.dbh.2 in consequence
          Mat.Range.dbh.2$EndClass[k] <- Mat.Range.dbh.2$EndClass[k+1]
          Mat.Range.dbh.2$dbh.class[k] <- paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k],sep=",")
          Mat.Range.dbh.2$n.Alive[k] <- Mat.Range.dbh.2$n.Alive[k]+Mat.Range.dbh.2$n.Alive[k+1]
          Mat.Range.dbh.2$n.Dead[k] <- Mat.Range.dbh.2$n.Dead[k]+Mat.Range.dbh.2$n.Dead[k+1]
          Mat.Range.dbh.2 <- Mat.Range.dbh.2[-(k+1),]
          ############################################
          ### New estimate of theta for each dbh.class
          Levels.dbh.breaks <- sort(unique(Mat.Range.dbh.2$dbh.class))
          n.dbh.breaks <- length(Levels.dbh.breaks)
          # Mortality rate estimate for each dbhclass
          for (i in 1:n.dbh.breaks) {
            # Variables for considered dbhclass
            data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
            Y.dead <- data.temp$ti[data.temp$mort==1]
            Y.alive <- data.temp$ti[data.temp$mort==0]
            if (length(data.temp$ti[data.temp$mort==1])==0) {
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- 0
            }
            if (length(data.temp$ti[data.temp$mort==1])>0) { 
              # Minimizing -loglikelihood
              Min.Minus.LL <- optim(par=0.01,
                                    fn=Binom.Likelihood,
                                    method="L-BFGS-B",
                                    lower=0.000001,upper=0.999999,
                                    Y.dead=Y.dead,Y.alive=Y.alive) 
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- Min.Minus.LL$par
            }
          }
          ############################
          # We actualize Order
          Order <- vector()
          if (length(Mat.Range.dbh.2$theta.MLE)>1) {
            for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
              Order[j] <- Mat.Range.dbh.2$theta.MLE[j]>=Mat.Range.dbh.2$theta.MLE[j+1]|
                        Mat.Range.dbh.2$theta.MLE[j]==0|
                        Mat.Range.dbh.2$theta.MLE[j+1]==0
            }
          }
          else { Order <- FALSE }
          break
        }
    }
  }
  
  MatRange.increase <- Mat.Range.dbh.2
  MatRange.all <- MatRange.increase
  
  ######## Deviance calcul #########
  Levels.dbh.breaks <- sort(unique(MatRange.all$dbh.class))
  n.dbh.breaks <- length(Levels.dbh.breaks)
  MatRange.all$LLikelihood <- rep(NA,n.dbh.breaks)
  # Mortality rate estimate for each dbhclass
  for (i in 1:n.dbh.breaks) {
    data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
    Y.dead <- data.temp$ti[data.temp$mort==1]
    Y.alive <- data.temp$ti[data.temp$mort==0]
    MatRange.all$LLikelihood[MatRange.all$dbh.class==Levels.dbh.breaks[i]] <- 
              sum(log(1-(1-MatRange.all$theta.MLE[MatRange.all$dbh.class==Levels.dbh.breaks[i]])^Y.dead))+
              sum(log((1-MatRange.all$theta.MLE[MatRange.all$dbh.class==Levels.dbh.breaks[i]])^Y.alive))
  }
  Deviance[1] <- -2*sum(MatRange.all$LLikelihood)
  
  #### data is reinitiate to its initial value
  data$dbh.breaks <- cut(data$dbh,breaks=seq(dbh.min,dbh.max,by=init.bin),right=FALSE)
  data$dbh.breaks <- as.character(data$dbh.breaks)
  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Stop=2: One bin then increase
Stop <- 2

  # Decrease for dbh<Stop
  Mat.Range.dbh.2 <- Mat.Range.dbh[-seq(from=Stop,to=length(Mat.Range.dbh$dbh.class)),]
  ############################################
  ### New estimate of theta for each dbh.class
  Levels.dbh.breaks <- sort(unique(Mat.Range.dbh.2$dbh.class))
  n.dbh.breaks <- length(Levels.dbh.breaks)
  # Mortality rate estimate for each dbhclass
  for (i in 1:n.dbh.breaks) {
    # Variables for considered dbhclass
    data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
    Y.dead <- data.temp$ti[data.temp$mort==1]
    Y.alive <- data.temp$ti[data.temp$mort==0]
    if (length(data.temp$ti[data.temp$mort==1])==0) {
      Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- 0
    }
    if (length(data.temp$ti[data.temp$mort==1])>0) { 
      # Minimizing -loglikelihood
      Min.Minus.LL <- optim(par=0.01,
                            fn=Binom.Likelihood,
                            method="L-BFGS-B",
                            lower=0.000001,upper=0.999999,
                            Y.dead=Y.dead,Y.alive=Y.alive) 
      Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- Min.Minus.LL$par
    }
  } 
  MatRange.decrease <- Mat.Range.dbh.2
  
  # Increase for dbh>Stop
  Mat.Range.dbh.2 <- Mat.Range.dbh[-seq(from=1,to=Stop-1),]
  # Order will check for monotonicity
  Order <- vector()
  for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
    Order[j] <- Mat.Range.dbh.2$theta.MLE[j]>=Mat.Range.dbh.2$theta.MLE[j+1]|
              Mat.Range.dbh.2$theta.MLE[j]==0|
              Mat.Range.dbh.2$theta.MLE[j+1]==0
  }
  # Loops
  while (sum(Order)!=0) {
    for (k in 1:(length(Mat.Range.dbh.2[,1])-1)) {
        # We check for unordered bins between k and k+1
        if (Mat.Range.dbh.2$theta.MLE[k]>=Mat.Range.dbh.2$theta.MLE[k+1]|
            Mat.Range.dbh.2$theta.MLE[k]==0|
            Mat.Range.dbh.2$theta.MLE[k+1]==0) {
          # We replace the values of the bin in data
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k+1]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          # We modify Mat.Range.dbh.2 in consequence
          Mat.Range.dbh.2$EndClass[k] <- Mat.Range.dbh.2$EndClass[k+1]
          Mat.Range.dbh.2$dbh.class[k] <- paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k],sep=",")
          Mat.Range.dbh.2$n.Alive[k] <- Mat.Range.dbh.2$n.Alive[k]+Mat.Range.dbh.2$n.Alive[k+1]
          Mat.Range.dbh.2$n.Dead[k] <- Mat.Range.dbh.2$n.Dead[k]+Mat.Range.dbh.2$n.Dead[k+1]
          Mat.Range.dbh.2 <- Mat.Range.dbh.2[-(k+1),]
          ############################################
          ### New estimate of theta for each dbh.class
          Levels.dbh.breaks <- sort(unique(Mat.Range.dbh.2$dbh.class))
          n.dbh.breaks <- length(Levels.dbh.breaks)
          # Mortality rate estimate for each dbhclass
          for (i in 1:n.dbh.breaks) {
            # Variables for considered dbhclass
            data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
            Y.dead <- data.temp$ti[data.temp$mort==1]
            Y.alive <- data.temp$ti[data.temp$mort==0]
            if (length(data.temp$ti[data.temp$mort==1])==0) {
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- 0
            }
            if (length(data.temp$ti[data.temp$mort==1])>0) { 
              # Minimizing -loglikelihood
              Min.Minus.LL <- optim(par=0.01,
                                    fn=Binom.Likelihood,
                                    method="L-BFGS-B",
                                    lower=0.000001,upper=0.999999,
                                    Y.dead=Y.dead,Y.alive=Y.alive) 
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- Min.Minus.LL$par
            }
          }
          ############################
          # We actualize Order
          Order <- vector()
          if (length(Mat.Range.dbh.2$theta.MLE)>1) {
            for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
              Order[j] <- Mat.Range.dbh.2$theta.MLE[j]>=Mat.Range.dbh.2$theta.MLE[j+1]|
                        Mat.Range.dbh.2$theta.MLE[j]==0|
                        Mat.Range.dbh.2$theta.MLE[j+1]==0
            }
          }
          else { Order <- FALSE }
          break
        }
    }
  }
  
  MatRange.increase <- Mat.Range.dbh.2
  MatRange.all <- rbind(MatRange.decrease,MatRange.increase)
  
  ######## Deviance calcul #########
  Levels.dbh.breaks <- sort(unique(MatRange.all$dbh.class))
  n.dbh.breaks <- length(Levels.dbh.breaks)
  MatRange.all$LLikelihood <- rep(NA,n.dbh.breaks)
  # Mortality rate estimate for each dbhclass
  for (i in 1:n.dbh.breaks) {
    data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
    Y.dead <- data.temp$ti[data.temp$mort==1]
    Y.alive <- data.temp$ti[data.temp$mort==0]
    MatRange.all$LLikelihood[MatRange.all$dbh.class==Levels.dbh.breaks[i]] <- 
              sum(log(1-(1-MatRange.all$theta.MLE[MatRange.all$dbh.class==Levels.dbh.breaks[i]])^Y.dead))+
              sum(log((1-MatRange.all$theta.MLE[MatRange.all$dbh.class==Levels.dbh.breaks[i]])^Y.alive))
  }
  Deviance[2] <- -2*sum(MatRange.all$LLikelihood)
  
  #### data is reinitiate to its initial value
  data$dbh.breaks <- cut(data$dbh,breaks=seq(dbh.min,dbh.max,by=init.bin),right=FALSE)
  data$dbh.breaks <- as.character(data$dbh.breaks)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Stop=length(Mat.Range.dbh$dbh.class)-1: Decrease then one bin
Stop <- length(Mat.Range.dbh$dbh.class)-1

# Decrease for dbh<Stop
  Mat.Range.dbh.2 <- Mat.Range.dbh[-seq(from=Stop+1,to=length(Mat.Range.dbh$dbh.class)),]
  # Order will check for monotonicity
  Order <- vector()
  for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
    Order[j] <- Mat.Range.dbh.2$theta.MLE[j]<=Mat.Range.dbh.2$theta.MLE[j+1]|
    Mat.Range.dbh.2$theta.MLE[j]==0|
    Mat.Range.dbh.2$theta.MLE[j+1]==0
  }
  # Loops
  while (sum(Order)!=0) {
    for (k in 1:(length(Mat.Range.dbh.2[,1])-1)) {
        # We check for unordered bins between k and k+1
        if (Mat.Range.dbh.2$theta.MLE[k]<=Mat.Range.dbh.2$theta.MLE[k+1]|
            Mat.Range.dbh.2$theta.MLE[k]==0|
            Mat.Range.dbh.2$theta.MLE[k+1]==0) {
          # We replace the values of the bin in data
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k+1]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          # We modify Mat.Range.dbh.2 in consequence
          Mat.Range.dbh.2$EndClass[k] <- Mat.Range.dbh.2$EndClass[k+1]
          Mat.Range.dbh.2$dbh.class[k] <- paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k],sep=",")
          Mat.Range.dbh.2$n.Alive[k] <- Mat.Range.dbh.2$n.Alive[k]+Mat.Range.dbh.2$n.Alive[k+1]
          Mat.Range.dbh.2$n.Dead[k] <- Mat.Range.dbh.2$n.Dead[k]+Mat.Range.dbh.2$n.Dead[k+1]
          Mat.Range.dbh.2 <- Mat.Range.dbh.2[-(k+1),]
          ############################################
          ### New estimate of theta for each dbh.class
          Levels.dbh.breaks <- sort(unique(Mat.Range.dbh.2$dbh.class))
          n.dbh.breaks <- length(Levels.dbh.breaks)
          # Mortality rate estimate for each dbhclass
          for (i in 1:n.dbh.breaks) {
            # Variables for considered dbhclass
            data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
            Y.dead <- data.temp$ti[data.temp$mort==1]
            Y.alive <- data.temp$ti[data.temp$mort==0]
            if (length(data.temp$ti[data.temp$mort==1])==0) {
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- 0
            }
            if (length(data.temp$ti[data.temp$mort==1])>0) { 
              # Minimizing -loglikelihood
              Min.Minus.LL <- optim(par=0.01,
                                    fn=Binom.Likelihood,
                                    method="L-BFGS-B",
                                    lower=0.000001,upper=0.999999,
                                    Y.dead=Y.dead,Y.alive=Y.alive) 
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- Min.Minus.LL$par
            }
          }
          ############################
          # We actualize Order
          Order <- vector()
          if (length(Mat.Range.dbh.2$theta.MLE)>1) {
            for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
              Order[j] <- Mat.Range.dbh.2$theta.MLE[j]<=Mat.Range.dbh.2$theta.MLE[j+1]|
              Mat.Range.dbh.2$theta.MLE[j]==0|
              Mat.Range.dbh.2$theta.MLE[j+1]==0
            }
          }
          else { Order <- FALSE }
          break
        }
    }
  }
  
  MatRange.decrease <- Mat.Range.dbh.2
  
  # Increase for dbh>Stop
  Mat.Range.dbh.2 <- Mat.Range.dbh[-seq(from=1,to=Stop),]
  ############################################
  ### New estimate of theta for each dbh.class
  Levels.dbh.breaks <- sort(unique(Mat.Range.dbh.2$dbh.class))
  n.dbh.breaks <- length(Levels.dbh.breaks)
  # Mortality rate estimate for each dbhclass
  for (i in 1:n.dbh.breaks) {
    # Variables for considered dbhclass
    data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
    Y.dead <- data.temp$ti[data.temp$mort==1]
    Y.alive <- data.temp$ti[data.temp$mort==0]
    if (length(data.temp$ti[data.temp$mort==1])==0) {
      Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- 0
    }
    if (length(data.temp$ti[data.temp$mort==1])>0) { 
      # Minimizing -loglikelihood
      Min.Minus.LL <- optim(par=0.01,
                            fn=Binom.Likelihood,
                            method="L-BFGS-B",
                            lower=0.000001,upper=0.999999,
                            Y.dead=Y.dead,Y.alive=Y.alive) 
      Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- Min.Minus.LL$par
    }
  }
           
  MatRange.increase <- Mat.Range.dbh.2
  MatRange.all <- rbind(MatRange.decrease,MatRange.increase)
  
  ######## Deviance calcul #########
  Levels.dbh.breaks <- sort(unique(MatRange.all$dbh.class))
  n.dbh.breaks <- length(Levels.dbh.breaks)
  MatRange.all$LLikelihood <- rep(NA,n.dbh.breaks)
  # Mortality rate estimate for each dbhclass
  for (i in 1:n.dbh.breaks) {
    data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
    Y.dead <- data.temp$ti[data.temp$mort==1]
    Y.alive <- data.temp$ti[data.temp$mort==0]
    MatRange.all$LLikelihood[MatRange.all$dbh.class==Levels.dbh.breaks[i]] <- 
              sum(log(1-(1-MatRange.all$theta.MLE[MatRange.all$dbh.class==Levels.dbh.breaks[i]])^Y.dead))+
              sum(log((1-MatRange.all$theta.MLE[MatRange.all$dbh.class==Levels.dbh.breaks[i]])^Y.alive))
  }
  Deviance[length(Mat.Range.dbh$dbh.class)] <- -2*sum(MatRange.all$LLikelihood)
  
  #### data is reinitiate to its initial value
  data$dbh.breaks <- cut(data$dbh,breaks=seq(dbh.min,dbh.max,by=init.bin),right=FALSE)
  data$dbh.breaks <- as.character(data$dbh.breaks)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Stop=length(Mat.Range.dbh$dbh.class): Only decrease
Stop <- length(Mat.Range.dbh$dbh.class)

# Decrease for dbh<Stop
  Mat.Range.dbh.2 <- Mat.Range.dbh
  # Order will check for monotonicity
  Order <- vector()
  for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
    Order[j] <- Mat.Range.dbh.2$theta.MLE[j]<=Mat.Range.dbh.2$theta.MLE[j+1]|
    Mat.Range.dbh.2$theta.MLE[j]==0|
    Mat.Range.dbh.2$theta.MLE[j+1]==0
  }
  # Loops
  while (sum(Order)!=0) {
    for (k in 1:(length(Mat.Range.dbh.2[,1])-1)) {
        # We check for unordered bins between k and k+1
        if (Mat.Range.dbh.2$theta.MLE[k]<=Mat.Range.dbh.2$theta.MLE[k+1]|
            Mat.Range.dbh.2$theta.MLE[k]==0|
            Mat.Range.dbh.2$theta.MLE[k+1]==0) {
          # We replace the values of the bin in data
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k+1]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          # We modify Mat.Range.dbh.2 in consequence
          Mat.Range.dbh.2$EndClass[k] <- Mat.Range.dbh.2$EndClass[k+1]
          Mat.Range.dbh.2$dbh.class[k] <- paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k],sep=",")
          Mat.Range.dbh.2$n.Alive[k] <- Mat.Range.dbh.2$n.Alive[k]+Mat.Range.dbh.2$n.Alive[k+1]
          Mat.Range.dbh.2$n.Dead[k] <- Mat.Range.dbh.2$n.Dead[k]+Mat.Range.dbh.2$n.Dead[k+1]
          Mat.Range.dbh.2 <- Mat.Range.dbh.2[-(k+1),]
          ############################################
          ### New estimate of theta for each dbh.class of Mat.Range.dbh.2
          Levels.dbh.breaks <- sort(unique(Mat.Range.dbh.2$dbh.class))
          n.dbh.breaks <- length(Levels.dbh.breaks)
          # Mortality rate estimate for each dbhclass
          for (i in 1:n.dbh.breaks) {
            # Variables for considered dbhclass
            data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
            Y.dead <- data.temp$ti[data.temp$mort==1]
            Y.alive <- data.temp$ti[data.temp$mort==0]
            if (length(data.temp$ti[data.temp$mort==1])==0) {
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- 0
            }
            if (length(data.temp$ti[data.temp$mort==1])>0) { 

              # Minimizing -loglikelihood
              Min.Minus.LL <- optim(par=0.01,
                                    fn=Binom.Likelihood,
                                    method="L-BFGS-B",
                                    lower=0.000001,upper=0.999999,
                                    Y.dead=Y.dead,Y.alive=Y.alive) 
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- Min.Minus.LL$par
            }
          }
          ############################
          # We actualize Order
          Order <- vector()
          if (length(Mat.Range.dbh.2$theta.MLE)>1) {
            for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
              Order[j] <- Mat.Range.dbh.2$theta.MLE[j]<=Mat.Range.dbh.2$theta.MLE[j+1]|
                        Mat.Range.dbh.2$theta.MLE[j]==0|
                        Mat.Range.dbh.2$theta.MLE[j+1]==0
            }
          }
          else { Order <- FALSE }
          break
        }
    }
  }
  
  MatRange.decrease <- Mat.Range.dbh.2
  MatRange.all <- MatRange.decrease
  
  ######## Deviance calcul #########
  Levels.dbh.breaks <- sort(unique(MatRange.all$dbh.class))
  n.dbh.breaks <- length(Levels.dbh.breaks)
  MatRange.all$LLikelihood <- rep(NA,n.dbh.breaks)
  # Mortality rate estimate for each dbhclass
  for (i in 1:n.dbh.breaks) {
    data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
    Y.dead <- data.temp$ti[data.temp$mort==1]
    Y.alive <- data.temp$ti[data.temp$mort==0]
    MatRange.all$LLikelihood[MatRange.all$dbh.class==Levels.dbh.breaks[i]] <- 
              sum(log(1-(1-MatRange.all$theta.MLE[MatRange.all$dbh.class==Levels.dbh.breaks[i]])^Y.dead))+
              sum(log((1-MatRange.all$theta.MLE[MatRange.all$dbh.class==Levels.dbh.breaks[i]])^Y.alive))
  }
  Deviance[length(Mat.Range.dbh$dbh.class)+1] <- -2*sum(MatRange.all$LLikelihood)
  
  #### data is reinitiate to its initial value
  data$dbh.breaks <- cut(data$dbh,breaks=seq(dbh.min,dbh.max,by=init.bin),right=FALSE)
  data$dbh.breaks <- as.character(data$dbh.breaks)

#== End of particular case
#================================================================================
    
#================================================================================
#= Graphic: deviance evolution with dbh0

cat("\nPlot graphic...")  
    
test.dbh0 <- seq(from=dbh.min,to=dbh.max,by=init.bin)
pdf(file="deviance-dbh0.pdf")
par(cex=1.4)
plot(test.dbh0,Deviance,
    type="l",
    xlab="dbh0 (cm)",
    ylab="Deviance")
Deviance.Best.Index <- ceiling(mean(which(Deviance==min(Deviance))))
abline(v=test.dbh0[Deviance.Best.Index],lty=2)
dev.off()

dbh0 <- test.dbh0[Deviance.Best.Index]
deviance0 <- Deviance[Deviance.Best.Index]

cat(paste("\ndbh0 = ",dbh0,sep=""))
cat(paste("\ndeviance = ",deviance0,sep=""))
    
#================================================================================
#= Model with identified dbh0

  cat("\nEstimates for the best model...\n")   

  Stop <- Deviance.Best.Index
  # Decrease for dbh<Stop
  Mat.Range.dbh.2 <- Mat.Range.dbh[-seq(from=Stop,to=length(Mat.Range.dbh$dbh.class)),]
  # Order will check for monotonicity
  Order <- vector()
  for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
    Order[j] <- Mat.Range.dbh.2$theta.MLE[j]<=Mat.Range.dbh.2$theta.MLE[j+1]|
              Mat.Range.dbh.2$theta.MLE[j]==0|
              Mat.Range.dbh.2$theta.MLE[j+1]==0
  }
  # Loops
  while (sum(Order)!=0) {
    for (k in 1:(length(Mat.Range.dbh.2[,1])-1)) {
        # We check for unordered bins between k and k+1
        if (Mat.Range.dbh.2$theta.MLE[k]<=Mat.Range.dbh.2$theta.MLE[k+1]|
            Mat.Range.dbh.2$theta.MLE[k]==0|
            Mat.Range.dbh.2$theta.MLE[k+1]==0) {
          # We replace the values of the bin in data
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k+1]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          # We modify Mat.Range.dbh.2 in consequence
          Mat.Range.dbh.2$EndClass[k] <- Mat.Range.dbh.2$EndClass[k+1]
          Mat.Range.dbh.2$dbh.class[k] <- paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k],sep=",")
          Mat.Range.dbh.2$n.Alive[k] <- Mat.Range.dbh.2$n.Alive[k]+Mat.Range.dbh.2$n.Alive[k+1]
          Mat.Range.dbh.2$n.Dead[k] <- Mat.Range.dbh.2$n.Dead[k]+Mat.Range.dbh.2$n.Dead[k+1]
          Mat.Range.dbh.2 <- Mat.Range.dbh.2[-(k+1),]
          ############################################
          ### New estimate of theta for each dbh.class of Mat.Range.dbh.2
          Levels.dbh.breaks <- sort(unique(Mat.Range.dbh.2$dbh.class))
          n.dbh.breaks <- length(Levels.dbh.breaks)
          # Mortality rate estimate for each dbhclass
          for (i in 1:n.dbh.breaks) {
            # Variables for considered dbhclass
            data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
            Y.dead <- data.temp$ti[data.temp$mort==1]
            Y.alive <- data.temp$ti[data.temp$mort==0]
            if (length(data.temp$ti[data.temp$mort==1])==0) {
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- 0
            }
            if (length(data.temp$ti[data.temp$mort==1])>0) { 
              # Minimizing -loglikelihood
              Min.Minus.LL <- optim(par=0.01,
		                       fn=Binom.Likelihood,
		                       method="L-BFGS-B",
		                       lower=0.000001,upper=0.999999,
		                       Y.dead=Y.dead,Y.alive=Y.alive) 
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- Min.Minus.LL$par
            }
          }
          ############################
          # We actualize Order
          Order <- vector()
          if (length(Mat.Range.dbh.2$theta.MLE)>1) {
            for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
            Order[j] <- Mat.Range.dbh.2$theta.MLE[j]<=Mat.Range.dbh.2$theta.MLE[j+1]|
                      Mat.Range.dbh.2$theta.MLE[j]==0|
                      Mat.Range.dbh.2$theta.MLE[j+1]==0
            }
          }
          else { Order <- FALSE }
          break
        }
    }
  }
  
  MatRange.decrease <- Mat.Range.dbh.2
  
  # Increase for dbh>Stop
  Mat.Range.dbh.2 <- Mat.Range.dbh[-seq(from=1,to=Stop-1),]
  # Order will check for monotonicity
  Order <- vector()
  for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
    Order[j] <- Mat.Range.dbh.2$theta.MLE[j]>=Mat.Range.dbh.2$theta.MLE[j+1]|
              Mat.Range.dbh.2$theta.MLE[j]==0|
              Mat.Range.dbh.2$theta.MLE[j+1]==0
  }
  # Loops
  while (sum(Order)!=0) {
    for (k in 1:(length(Mat.Range.dbh.2[,1])-1)) {
        # We check for unordered bins between k and k+1
        if (Mat.Range.dbh.2$theta.MLE[k]>=Mat.Range.dbh.2$theta.MLE[k+1]|
            Mat.Range.dbh.2$theta.MLE[k]==0|
            Mat.Range.dbh.2$theta.MLE[k+1]==0) {
          # We replace the values of the bin in data
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          data$dbh.breaks[data$dbh.breaks==Mat.Range.dbh.2$dbh.class[k+1]] <- 
                                paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k+1],sep=",")
          # We modify Mat.Range.dbh.2 in consequence
          Mat.Range.dbh.2$EndClass[k] <- Mat.Range.dbh.2$EndClass[k+1]
          Mat.Range.dbh.2$dbh.class[k] <- paste(Mat.Range.dbh.2$BeginClass[k],Mat.Range.dbh.2$EndClass[k],sep=",")
          Mat.Range.dbh.2$n.Alive[k] <- Mat.Range.dbh.2$n.Alive[k]+Mat.Range.dbh.2$n.Alive[k+1]
          Mat.Range.dbh.2$n.Dead[k] <- Mat.Range.dbh.2$n.Dead[k]+Mat.Range.dbh.2$n.Dead[k+1]
          Mat.Range.dbh.2 <- Mat.Range.dbh.2[-(k+1),]
          ############################################
          ### New estimate of theta for each dbh.class
          Levels.dbh.breaks <- sort(unique(Mat.Range.dbh.2$dbh.class))
          n.dbh.breaks <- length(Levels.dbh.breaks)
          # Mortality rate estimate for each dbhclass
          for (i in 1:n.dbh.breaks) {
            # Variables for considered dbhclass
            data.temp <- subset(data,data$dbh.breaks==Levels.dbh.breaks[i])
            Y.dead <- data.temp$ti[data.temp$mort==1]
            Y.alive <- data.temp$ti[data.temp$mort==0]
            if (length(data.temp$ti[data.temp$mort==1])==0) {
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- 0
            }
            if (length(data.temp$ti[data.temp$mort==1])>0) { 
              # Minimizing -loglikelihood
              Min.Minus.LL <- optim(par=0.01,
                                    fn=Binom.Likelihood,
                                    method="L-BFGS-B",
                                    lower=0.000001,upper=0.999999,
                                    Y.dead=Y.dead,Y.alive=Y.alive)
              Mat.Range.dbh.2$theta.MLE[Mat.Range.dbh.2$dbh.class==Levels.dbh.breaks[i]] <- Min.Minus.LL$par
            }
          }
          ############################
          # We actualize Order
          Order <- vector()
          if (length(Mat.Range.dbh.2$theta.MLE)>1) {
            for (j in 1:(length(Mat.Range.dbh.2[,1])-1)) {
              Order[j] <- Mat.Range.dbh.2$theta.MLE[j]>=Mat.Range.dbh.2$theta.MLE[j+1]|
                        Mat.Range.dbh.2$theta.MLE[j]==0|
                        Mat.Range.dbh.2$theta.MLE[j+1]==0
            }
          }
          else { Order <- FALSE }
          break
        }
    }
  }
  
  MatRange.increase <- Mat.Range.dbh.2
  MatRange.all <- rbind(MatRange.decrease,MatRange.increase)

#====================================================
#= Bayesian estimation of mortality rate for each bin

mat <- MatRange.all
nbin <- nrow(mat)
for (i in 1:nbin) {
  y <- c(rep(0,mat$n.Alive[i]),rep(1,mat$n.Dead[i])) 
  mod <- MCMClogit(y~1, verbose=0, tune=2, burnin=1000, mcmc=1000, thin=1)
  mat$theta.Bayes[i] <- mean(inv.logit(as.matrix(mod)[,1]))
  mat$theta.q1[i] <- quantile(inv.logit(as.matrix(mod)[,1]),0.025)
  mat$theta.q3[i] <- quantile(inv.logit(as.matrix(mod)[,1]),0.975)
}

#===========================
#= Graphic

# Some vectors to compute the width of each bin
l <- as.numeric(gsub("\\[","",mat$BeginClass))
u <- as.numeric(gsub(")","",mat$EndClass))

pdf("mortality-dbh.pdf")
par(mar=c(5,4,3,2),cex=1.4)
# barplot for MLE
x.at <- barplot(height=mat$theta.MLE,width=c(u-l),
                ylim=c(0,1),
                xlab="DBH (cm)",
                ylab="Annual mortality rate",
                col="transparent",
                border="black",
                axes=FALSE,
                space=0)
axis(1,at=c(dbh.min-dbh.min,dbh0-dbh.min,dbh.max-dbh.min),labels=c(dbh.min,paste("dbh0=",dbh0,sep=""),dbh.max))
axis(2,at=c(0,1),labels=c(0,1))

# posterior from Bayesian estimation
points(x.at,mat$theta.Bayes,pch=16)
lines(x.at,mat$theta.Bayes)
lines(x.at,mat$theta.q1,lty=2)
lines(x.at,mat$theta.q3,lty=2)
dev.off()

return (list(dbh0=dbh0,deviance=deviance0,model=mat))
    
} # End ayer.dbh function

#===============================================================================================================
#===============================================================================================================

