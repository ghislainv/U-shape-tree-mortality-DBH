##
## ayer.dbh.comp.R
##
## Author: Ghislain Vieilledent <ghislain.vieilledent@cirad.fr>
##
## Date: 12-01-2012
##

## ayer.dbh.comp() function fits a semi-parametric mortality model as
## a function of the tree dbh and neighborhood competition. The dbh at
## which tree mortality rate is minimal (dbh0) must be provided. The
## function estimates the dbh and competition bins and the
## corresponding annual mortality rates to obtain a decreasing
## mortality with dbh before dbh0, an increasing mortality with dbh
## beyond dbh0 and an increasing mortality with competition.

## Parameters:
##
## mort: vector of 0 or 1 for tree bernoulli mortality process (1=dead, 0=alive) between census c and c+1
## dbh: vector for tree diameter in cm at census c
## comp: vector for competition in m2.ha-1 in the neighborhood of the target tree at census c
## ti: time interval in year between census c and c+1
## dbh0: dbh in cm at which tree mortality rate is minimal
## init.bin.dbh: initial bin width in cm for the dbh (5 cm might be an appropriate value)
## init.bin.comp: initial bin width in m2.ha-1 for the competition (10 m^2.ha-1 cm might be an appropriate value)

## Value:
##
## deviance: model deviance
##
## living: matrix with number of living trees for each dbh and competition bin
##
## dead: matrix with number of dead trees for each dbh and competition bin
##
## model: matrix describing the model

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
##
## Vieilledent, G.; Courbaud, B.; Kunstler, G. & Dhôte,
## J.-F. Mortality of Silver Fir and Norway Spruce in the Western Alps
## – A semi-parametric approach combining size-dependent and
## growth-dependent mortality Annals of Forest Science, 2010, 67(3),
## 305


#= Function f.order.1 for [dbh.min,dbh0]
f.order.1 <- function (D.Ayer,G.Ayer,Mean.Combined.Post) {
  Order <- matrix(FALSE,nrow=D.Ayer,ncol=G.Ayer)
  for (j in 1:(D.Ayer-1)) {
    for (k in 1:(G.Ayer-1)) {
      Order[j,k] <- (Mean.Combined.Post[j,k]>Mean.Combined.Post[j,k+1]|
                     Mean.Combined.Post[j,k]<Mean.Combined.Post[j+1,k])
    }
  }
  for (j in 1:(D.Ayer-1)) {
    Order[j,G.Ayer] <- Mean.Combined.Post[j,G.Ayer]<Mean.Combined.Post[j+1,G.Ayer]
  }
  for (k in 1:(G.Ayer-1)) {
    Order[D.Ayer,k] <- Mean.Combined.Post[D.Ayer,k]>Mean.Combined.Post[D.Ayer,k+1]
  }
  return (Order)
}

#= Function f.order.2 for [dbh0,dbh.max]
f.order.2 <- function (D.Ayer,G.Ayer,Mean.Combined.Post) {
  Order <- matrix(FALSE,nrow=D.Ayer,ncol=G.Ayer)
  for (j in 1:(D.Ayer-1)) {
    for (k in 1:(G.Ayer-1)) {
      Order[j,k] <- (Mean.Combined.Post[j,k]>Mean.Combined.Post[j,k+1]|
                     Mean.Combined.Post[j,k]>Mean.Combined.Post[j+1,k])
    }
  }
  for (j in 1:(D.Ayer-1)) {
    Order[j,G.Ayer] <- Mean.Combined.Post[j,G.Ayer]>Mean.Combined.Post[j+1,G.Ayer]
  }
  for (k in 1:(G.Ayer-1)) {
    Order[D.Ayer,k] <- Mean.Combined.Post[D.Ayer,k]>Mean.Combined.Post[D.Ayer,k+1]
  }
  return(Order)
}

#= Function assessing the conditions on j for recombining bins for [dbh.min,dbh0]
f.conditions.j.1 <- function(Mean.Combined.Post,Matrix.Group,j,k) {
  nas <- ones <- not.mono <- nasd <- FALSE
  # 1st condition: for NAs
  if (is.na(Mean.Combined.Post[j,k])|is.na(Mean.Combined.Post[j+1,k])) {nas <- TRUE}
  if (!nas) {
    # 2sd condition: for ones
    if (((Mean.Combined.Post[j,k]==1|Mean.Combined.Post[j+1,k]==1)&
         (Matrix.Group[j,k]!=Matrix.Group[j+1,k]))) {ones <- TRUE}
    # 3rd condition: for monotonicity
    if (Mean.Combined.Post[j,k]<Mean.Combined.Post[j+1,k]) {not.mono <- TRUE}
  }
  # double NA
  if (is.na(Mean.Combined.Post[j,k])&is.na(Mean.Combined.Post[j+1,k])) {nasd <- TRUE}
  return(list(all=nas|ones|not.mono,nasd=nasd))
}

#= Function assessing the conditions on j for recombining bins for [dbh0,dbh.max]
f.conditions.j.2 <- function(Mean.Combined.Post,Matrix.Group,j,k) {
  nas <- ones <- not.mono <- nasd <- FALSE
  # 1st condition: for NAs
  if (is.na(Mean.Combined.Post[j,k])|is.na(Mean.Combined.Post[j+1,k])) {nas <- TRUE}
  if (!nas) {
    # 2sd condition: for ones
    if (((Mean.Combined.Post[j,k]==1|Mean.Combined.Post[j+1,k]==1)&
         (Matrix.Group[j,k]!=Matrix.Group[j+1,k]))) {ones <- TRUE}
    # 3rd condition: for monotonicity
    if (Mean.Combined.Post[j,k]>Mean.Combined.Post[j+1,k]) {not.mono <- TRUE}
  }
  # double NA
  if (is.na(Mean.Combined.Post[j,k])&is.na(Mean.Combined.Post[j+1,k])) {nasd <- TRUE}
  return(list(all=nas|ones|not.mono,nasd=nasd))
}

#= Function assessing the conditions on k for recombining bins
f.conditions.k <- function(Mean.Combined.Post,Matrix.Group,j,k) {
  nas <- ones <- not.mono <- nasd <- FALSE
  # 1st condition: for NAs
  if (is.na(Mean.Combined.Post[j,k])|is.na(Mean.Combined.Post[j,k+1])) {nas <- TRUE}
  if (!nas) {
    # 2sd condition: for ones
    if (((Mean.Combined.Post[j,k]==1|Mean.Combined.Post[j,k+1]==1)&
         (Matrix.Group[j,k]!=Matrix.Group[j,k+1]))) {ones <- TRUE}
    # 3rd condition: for monotonicity
    if (Mean.Combined.Post[j,k]>Mean.Combined.Post[j,k+1]) {not.mono <- TRUE}
  }
  # double NA
  if (is.na(Mean.Combined.Post[j,k])&is.na(Mean.Combined.Post[j,k+1])) {nasd <- TRUE}
  return(list(all=nas|ones|not.mono,nasd=nasd))
}

#= Likelihood function
Binom.Likelihood <- function (theta,Y.dead,Y.alive) {
  logvrais <- sum(log(1-(1-theta)^Y.dead))+sum(log((1-theta)^Y.alive))
  return(-logvrais)
}

#= Function ML.Binom
ML.Binom <- function (Y.dead,Y.alive) {
  if (length(Y.dead)==0 & length(Y.alive)==0) { # No data
    return(NA)
  }
  if (length(Y.dead)==0 & length(Y.alive)>0) { # Only living trees
    return(0)
  }
  if (length(Y.dead)>0 & length(Y.alive)==0) { # Only dead trees
    return(1)
  }
  if (length(Y.dead)>0 & length(Y.alive)>0) {
    Min.Minus.LL <- optim(par=0.01,
                          fn=Binom.Likelihood,
                          method="L-BFGS-B",
                          lower=0.000001,upper=0.999999,
                          Y.dead=Y.dead,Y.alive=Y.alive)
    return(Min.Minus.LL$par)
  }
}

#= Ayer function
ayer.dbh.comp <- function (mort,dbh,comp,ti,dbh0,init.bin.dbh,init.bin.comp) {

#= Recompose data
data <- data.frame(dbh=dbh,comp=comp,ti=ti,mort=mort)

#= Compute minimal and maximal dbh for dbh interval
dbh.min <- floor(min(data$dbh)/init.bin.dbh)*init.bin.dbh
dbh.max <- ceiling(max(data$dbh)/init.bin.dbh)*init.bin.dbh

#= Compute minimal and maximal comp for comp interval
comp.min <- floor(min(data$comp)/init.bin.comp)*init.bin.comp
comp.max <- ceiling(max(data$comp)/init.bin.comp)*init.bin.comp

#= Constructing classes
data$dbh.class <- cut(data$dbh,breaks=seq(from=dbh.min,to=dbh.max,by=init.bin.dbh),right=F)
data$comp.class <- cut(data$comp,breaks=seq(from=comp.min,to=comp.max,by=init.bin.comp),right=F)

#= Matrix of dead trees by class of dbh and comp
Matrix.Dead.Tot <- table(data$dbh.class[data$mort==1],data$comp.class[data$mort==1])

#= Matrix of living trees by class of dbh and comp
Matrix.Living.Tot <- table(data$dbh.class[data$mort==0],data$comp.class[data$mort==0])

#= Identifying rows for [dbh.min,dbh0]
nr <- (dbh0-dbh.min)/init.bin.dbh

#===============================================================================================
# For dbh on [dbh.min, dbh0]
#===============================================================================================

#= Useful objects
Matrix.Dead <- Matrix.Dead.Tot[1:nr,]
Matrix.Living <- Matrix.Living.Tot[1:nr,]
D.Ayer <- dim(Matrix.Living)[1]
G.Ayer <- dim(Matrix.Living)[2]
dbh.class.names <- rownames(Matrix.Living)
comp.class.names <- colnames(Matrix.Living)

#= Initial groups
data$group <- 0
group.id <- 1
for (j in 1:D.Ayer) {
  for (k in 1:G.Ayer) {
    data$group[data$dbh.class==dbh.class.names[j] & data$comp.class==comp.class.names[k]] <- group.id
    group.id <- group.id+1
  }
}
# Matrix of initial groups
Matrix.Group <- matrix(c(1:(D.Ayer*G.Ayer)),nrow=D.Ayer,ncol=G.Ayer,byrow=T)

#########################################################################################################################
########## Mean.Combined.Post is the matrix with the estimated mean mortality rate
########## for each starting class of G and D

Mean.Combined.Post <- matrix(NA,nrow=D.Ayer,ncol=G.Ayer)

for (j in 1:D.Ayer) {
  for (k in 1:G.Ayer) {
    # Variables for considered dbh class
    data.sub <- subset(data,data$dbh.class==dbh.class.names[j] & data$comp.class==comp.class.names[k])
    Y.dead <- data.sub$ti[data.sub$mort==1]
    Y.alive <- data.sub$ti[data.sub$mort==0]
    # Maximizing loglikelihood
    Mean.Combined.Post[j,k] <- ML.Binom(Y.dead,Y.alive)
  }
}

Mean.Combined.Post <- as.data.frame(Mean.Combined.Post)
rownames(Mean.Combined.Post) <- dbh.class.names
names(Mean.Combined.Post) <- comp.class.names

################################################################################
################################################################################
########## Ayer algorithm for combined mortality rate ##########################

#===================================
#= We actualize objects for checking
Order <- f.order.1(D.Ayer,G.Ayer,Mean.Combined.Post)     
Check.One <- sum(Mean.Combined.Post==1,na.rm=TRUE)
Check.NA <- sum(is.na(Mean.Combined.Post))

# Loops
cat("First part on [dbh.min,dbh0]\n")
Number.Iterst.Half <- 0
while ((sum(Order,na.rm=TRUE)!=0)|(Check.One!=0)|(Check.NA!=0)) {
  Number.Iterst.Half <- Number.Iterst.Half+1
  for (j in 1:(D.Ayer-1)) {
    for (k in 1:(G.Ayer-1)) {
      #====================================================================
      # We check for unordered bins between k and k+1 (until line D.Ayer-1)
      
      # Conditions for recombination
      cond <- f.conditions.k(Mean.Combined.Post,Matrix.Group,j,k)

      # Total
      if (cond$all) {

        #=========================================================
        # Observations for group[j,k+1] are affected to group[j,k]
        data$group[data$group==Matrix.Group[j,k+1]] <- Matrix.Group[j,k]
        Matrix.Group[Matrix.Group==Matrix.Group[j,k+1]] <- Matrix.Group[j,k]
        #= Considered.Group
        Considered.Group <- Matrix.Group[j,k]
        Matrix.True.False <- (Matrix.Group==Considered.Group)
        
        #==================================================
        # Computing mortality rate for the considered group

        #= Variables for considered dbh class
        data.sub <- subset(data,data$group==Considered.Group)
        Y.dead <- data.sub$ti[data.sub$mort==1]
        Y.alive <- data.sub$ti[data.sub$mort==0]
        
        #= Maximizing loglikelihood
        Mean.Combined.Post[Matrix.True.False] <- ML.Binom(Y.dead,Y.alive)

        #===============================
        #= We actualize objects for checking
        Order <- f.order.1(D.Ayer,G.Ayer,Mean.Combined.Post)     
        Check.One <- sum(Mean.Combined.Post==1,na.rm=TRUE)
        Check.NA <- sum(is.na(Mean.Combined.Post))

        #= restart from j=1, k=1
        if (!cond$nasd) {
          break
          break
        }
      }

      #======================================================================
      # We check for unordered bins between j and j+1 (until column G.Ayer-1)

      # Conditions for recombination
      cond <- f.conditions.j.1(Mean.Combined.Post,Matrix.Group,j,k)

      # Total
      if (cond$all) {
      
        #=========================================================
        # Observations for group[j+1,k] are affected to group[j,k]
        data$group[data$group==Matrix.Group[j+1,k]] <- Matrix.Group[j,k]
        Matrix.Group[Matrix.Group==Matrix.Group[j+1,k]] <- Matrix.Group[j,k] 
        #= Considered.Group
        Considered.Group <- Matrix.Group[j,k]
        Matrix.True.False <- (Matrix.Group==Considered.Group)
        
        #============================================
        #= Computing mortality rate for the considered group

        #= Variables for considered dbh class
        data.sub <- subset(data,data$group==Considered.Group)
        Y.dead <- data.sub$ti[data.sub$mort==1]
        Y.alive <- data.sub$ti[data.sub$mort==0]
        
        #= Maximizing loglikelihood
        Mean.Combined.Post[Matrix.True.False] <- ML.Binom(Y.dead,Y.alive)

        #===============================
        #= We actualize objects for checking
        Order <- f.order.1(D.Ayer,G.Ayer,Mean.Combined.Post)     
        Check.One <- sum(Mean.Combined.Post==1,na.rm=TRUE)
        Check.NA <- sum(is.na(Mean.Combined.Post))

        #= restart from j=1, k=1
        if (!cond$nasd) {
          break
          break
        }
      }

    } # End of k

    #=======================================================================
    # We check for unordered bins between j and j+1 (for last column G.Ayer)

    # Conditions for recombination
    cond <- f.conditions.j.1(Mean.Combined.Post,Matrix.Group,j,G.Ayer)

    # Total
    if (cond$all) {

      #=========================================================
      # Observations for group[j+1,G.Ayer] are affected to group[j,G.Ayer]
      data$group[data$group==Matrix.Group[j+1,G.Ayer]] <- Matrix.Group[j,G.Ayer]
      Matrix.Group[Matrix.Group==Matrix.Group[j+1,G.Ayer]] <- Matrix.Group[j,G.Ayer]
      #= Considered.Group
      Considered.Group <- Matrix.Group[j,G.Ayer]
      Matrix.True.False <- (Matrix.Group==Considered.Group)
        
      #============================================
      #= Computing mortality rate for the considered group

      #= Variables for considered dbh class
      data.sub <- subset(data,data$group==Considered.Group)
      Y.dead <- data.sub$ti[data.sub$mort==1]
      Y.alive <- data.sub$ti[data.sub$mort==0]
        
      #= Maximizing loglikelihood
      Mean.Combined.Post[Matrix.True.False] <- ML.Binom(Y.dead,Y.alive)

      #==============================
      #= We actualize objects for checking
      Order <- f.order.1(D.Ayer,G.Ayer,Mean.Combined.Post)     
      Check.One <- sum(Mean.Combined.Post==1,na.rm=TRUE)
      Check.NA <- sum(is.na(Mean.Combined.Post))

      #= restart from j=1
      if (!cond$nasd) {
        break
      }
    }
    
  } # End of j

  #=====================================================================
  # We check for unordered bins between k and k+1 (for last line D.Ayer)
  for (k in 1:(G.Ayer-1)) {

    # Conditions for recombination
    cond <- f.conditions.k(Mean.Combined.Post,Matrix.Group,D.Ayer,k)
    
    # Total
    if (cond$all) {

      #=========================================================
      # Observations for group[D.Ayer,k+1] are affected to group[D.Ayer,k]
      data$group[data$group==Matrix.Group[D.Ayer,k+1]] <- Matrix.Group[D.Ayer,k]
      Matrix.Group[Matrix.Group==Matrix.Group[D.Ayer,k+1]] <- Matrix.Group[D.Ayer,k]
      #= Considered.Group
      Considered.Group <- Matrix.Group[D.Ayer,k]
      Matrix.True.False <- (Matrix.Group==Considered.Group)
        
      #============================================
      #= Computing mortality rate for the considered group

      #= Variables for considered dbh class
      data.sub <- subset(data,data$group==Considered.Group)
      Y.dead <- data.sub$ti[data.sub$mort==1]
      Y.alive <- data.sub$ti[data.sub$mort==0]
    
      #= Maximizing loglikelihood
      Mean.Combined.Post[Matrix.True.False] <- ML.Binom(Y.dead,Y.alive)
      
      #==============================
      #= We actualize objects for checking
      Order <- f.order.1(D.Ayer,G.Ayer,Mean.Combined.Post)     
      Check.One <- sum(Mean.Combined.Post==1,na.rm=TRUE)
      Check.NA <- sum(is.na(Mean.Combined.Post))

      #= restart from k=1
      if (!cond$nasd) {
        break
      }
    }
  }
  print(paste("Current number of iterations:",Number.Iterst.Half,sep=" "))
}

Mean.Combined.Post.First.Half <- Mean.Combined.Post


#===============================================================================================
# For dbh on [dbh0, dbh.max]
#===============================================================================================

#= Useful objects
Matrix.Dead <- Matrix.Dead.Tot[(nr+1):nrow(Matrix.Dead.Tot),]
Matrix.Living <- Matrix.Living.Tot[(nr+1):nrow(Matrix.Living.Tot),]
D.Ayer <- dim(Matrix.Living)[1]
G.Ayer <- dim(Matrix.Living)[2]
dbh.class.names <- rownames(Matrix.Living)
comp.class.names <- colnames(Matrix.Living)

#= Groups
data$group <- 0
group.id <- 1
for (j in 1:D.Ayer) {
  for (k in 1:G.Ayer) {
    data$group[data$dbh.class==dbh.class.names[j] & data$comp.class==comp.class.names[k]] <- group.id
    group.id <- group.id+1
  }
}
# Matrix of initial groups
Matrix.Group <- matrix(c(1:(D.Ayer*G.Ayer)),nrow=D.Ayer,ncol=G.Ayer,byrow=T)

#########################################################################################################################
########## Mean.Combined.Post is the matrix with the estimated mean mortality rate
########## for each starting class of G and D

Mean.Combined.Post <- matrix(NA,nrow=D.Ayer,ncol=G.Ayer)

for (j in 1:D.Ayer) {
  for (k in 1:G.Ayer) {
    # Variables for considered dbh class
    data.sub <- subset(data,data$dbh.class==dbh.class.names[j] & data$comp.class==comp.class.names[k])
    Y.dead <- data.sub$ti[data.sub$mort==1]
    Y.alive <- data.sub$ti[data.sub$mort==0]
    # Maximizing loglikelihood
    Mean.Combined.Post[j,k] <- ML.Binom(Y.dead,Y.alive)
  }
}

Mean.Combined.Post <- as.data.frame(Mean.Combined.Post)
rownames(Mean.Combined.Post) <- dbh.class.names
names(Mean.Combined.Post) <- comp.class.names

################################################################################
################################################################################
########## Ayer algorithm for combined mortality rate ##########################

#===================================
#= We actualize objects for checking
Order <- f.order.2(D.Ayer,G.Ayer,Mean.Combined.Post) # Must be f.order.2 for [dbh0, dbh.max] !!
Check.One <- sum(Mean.Combined.Post==1,na.rm=TRUE)
Check.NA <- sum(is.na(Mean.Combined.Post))

# Loops
cat("Second part on [dbh0,dbh.max]\n")
Number.Iterst.Half <- 0
while ((sum(Order,na.rm=TRUE)!=0)|(Check.One!=0)|(Check.NA!=0)) {
  Number.Iterst.Half <- Number.Iterst.Half+1
  for (j in 1:(D.Ayer-1)) {
    for (k in 1:(G.Ayer-1)) {
      #====================================================================
      # We check for unordered bins between k and k+1 (until line D.Ayer-1)
      
      # Conditions for recombination
      cond <- f.conditions.k(Mean.Combined.Post,Matrix.Group,j,k)

      # Total
      if (cond$all) {

        #=========================================================
        # Observations for group[j,k+1] are affected to group[j,k]
        data$group[data$group==Matrix.Group[j,k+1]] <- Matrix.Group[j,k]
        Matrix.Group[Matrix.Group==Matrix.Group[j,k+1]] <- Matrix.Group[j,k]
        #= Considered.Group
        Considered.Group <- Matrix.Group[j,k]
        Matrix.True.False <- (Matrix.Group==Considered.Group)
        
        #==================================================
        # Computing mortality rate for the considered group

        #= Variables for considered dbh class
        data.sub <- subset(data,data$group==Considered.Group)
        Y.dead <- data.sub$ti[data.sub$mort==1]
        Y.alive <- data.sub$ti[data.sub$mort==0]
        
        #= Maximizing loglikelihood
        Mean.Combined.Post[Matrix.True.False] <- ML.Binom(Y.dead,Y.alive)

        #===============================
        #= We actualize objects for checking
        Order <- f.order.2(D.Ayer,G.Ayer,Mean.Combined.Post)     
        Check.One <- sum(Mean.Combined.Post==1,na.rm=TRUE)
        Check.NA <- sum(is.na(Mean.Combined.Post))

        #= restart from j=1, k=1
        if (!cond$nasd) {
          break
          break
        }
      }

      #======================================================================
      # We check for unordered bins between j and j+1 (until column G.Ayer-1)

      # Conditions for recombination
      cond <- f.conditions.j.2(Mean.Combined.Post,Matrix.Group,j,k) # Must be f.conditions.j.2 for [db0, dbh.max] !!

      # Total
      if (cond$all) {
      
        #=========================================================
        # Observations for group[j+1,k] are affected to group[j,k]
        data$group[data$group==Matrix.Group[j+1,k]] <- Matrix.Group[j,k]
        Matrix.Group[Matrix.Group==Matrix.Group[j+1,k]] <- Matrix.Group[j,k] 
        #= Considered.Group
        Considered.Group <- Matrix.Group[j,k]
        Matrix.True.False <- (Matrix.Group==Considered.Group)
        
        #============================================
        #= Computing mortality rate for the considered group

        #= Variables for considered dbh class
        data.sub <- subset(data,data$group==Considered.Group)
        Y.dead <- data.sub$ti[data.sub$mort==1]
        Y.alive <- data.sub$ti[data.sub$mort==0]
        
        #= Maximizing loglikelihood
        Mean.Combined.Post[Matrix.True.False] <- ML.Binom(Y.dead,Y.alive)

        #===============================
        #= We actualize objects for checking
        Order <- f.order.2(D.Ayer,G.Ayer,Mean.Combined.Post)     
        Check.One <- sum(Mean.Combined.Post==1,na.rm=TRUE)
        Check.NA <- sum(is.na(Mean.Combined.Post))

        #= restart from j=1, k=1
        if (!cond$nasd) {
          break
          break
        }      
      }

    } # End of k

    #=======================================================================
    # We check for unordered bins between j and j+1 (for last column G.Ayer)

    # Conditions for recombination
    cond <- f.conditions.j.2(Mean.Combined.Post,Matrix.Group,j,G.Ayer)

    # Total
    if (cond$all) {

      #=========================================================
      # Observations for group[j+1,G.Ayer] are affected to group[j,G.Ayer]
      data$group[data$group==Matrix.Group[j+1,G.Ayer]] <- Matrix.Group[j,G.Ayer]
      Matrix.Group[Matrix.Group==Matrix.Group[j+1,G.Ayer]] <- Matrix.Group[j,G.Ayer]
      #= Considered.Group
      Considered.Group <- Matrix.Group[j,G.Ayer]
      Matrix.True.False <- (Matrix.Group==Considered.Group)
        
      #============================================
      #= Computing mortality rate for the considered group

      #= Variables for considered dbh class
      data.sub <- subset(data,data$group==Considered.Group)
      Y.dead <- data.sub$ti[data.sub$mort==1]
      Y.alive <- data.sub$ti[data.sub$mort==0]
        
      #= Maximizing loglikelihood
      Mean.Combined.Post[Matrix.True.False] <- ML.Binom(Y.dead,Y.alive)

      #==============================
      #= We actualize objects for checking
      Order <- f.order.2(D.Ayer,G.Ayer,Mean.Combined.Post)     
      Check.One <- sum(Mean.Combined.Post==1,na.rm=TRUE)
      Check.NA <- sum(is.na(Mean.Combined.Post))

      #= restart from j=1
      if (!cond$nasd) {
        break
      }
    }
  } # End of j

  #=====================================================================
  # We check for unordered bins between k and k+1 (for last line D.Ayer)
  for (k in 1:(G.Ayer-1)) {

    # Conditions for recombination
    cond <- f.conditions.k(Mean.Combined.Post,Matrix.Group,D.Ayer,k)
    
    # Total
    if (cond$all) {

      #=========================================================
      # Observations for group[D.Ayer,k+1] are affected to group[D.Ayer,k]
      data$group[data$group==Matrix.Group[D.Ayer,k+1]] <- Matrix.Group[D.Ayer,k]
      Matrix.Group[Matrix.Group==Matrix.Group[D.Ayer,k+1]] <- Matrix.Group[D.Ayer,k]
      #= Considered.Group
      Considered.Group <- Matrix.Group[D.Ayer,k]
      Matrix.True.False <- (Matrix.Group==Considered.Group)
        
      #============================================
      #= Computing mortality rate for the considered group

      #= Variables for considered dbh class
      data.sub <- subset(data,data$group==Considered.Group)
      Y.dead <- data.sub$ti[data.sub$mort==1]
      Y.alive <- data.sub$ti[data.sub$mort==0]
    
      #= Maximizing loglikelihood
      Mean.Combined.Post[Matrix.True.False] <- ML.Binom(Y.dead,Y.alive)
      
      #==============================
      #= We actualize objects for checking
      Order <- f.order.2(D.Ayer,G.Ayer,Mean.Combined.Post)     
      Check.One <- sum(Mean.Combined.Post==1,na.rm=TRUE)
      Check.NA <- sum(is.na(Mean.Combined.Post))

      #= restart from k=1
      if (!cond$nasd) {
        break
      }
    }
  }
  print(paste("Current number of iterations:",Number.Iterst.Half,sep=" "))
}

Mean.Combined.Post.Second.Half <- Mean.Combined.Post

#= Combining the two parts
Mean.Post <- rbind(Mean.Combined.Post.First.Half,Mean.Combined.Post.Second.Half)

#======================
#= Deviance computation

#= Useful objects
D.Ayer <- dim(Matrix.Living.Tot)[1]
G.Ayer <- dim(Matrix.Living.Tot)[2]
dbh.class.names <- rownames(Matrix.Living.Tot)
comp.class.names <- colnames(Matrix.Living.Tot)

#= Deviance
deviance <- 0
for (j in 1:D.Ayer) {
  for (k in 1:G.Ayer) {
    data.sub <- data[data$dbh.class==dbh.class.names[j] & data$comp.class==comp.class.names[k],]
    Y.dead <- data.sub$ti[data.sub$mort==1]
    Y.alive <- data.sub$ti[data.sub$mort==0]
    deviance <- deviance + 2*Binom.Likelihood(theta=Mean.Post[j,k],Y.dead,Y.alive) # Binom.Likeihood return -logL
  }
}

return (list(deviance=deviance,living=Matrix.Living.Tot,dead=Matrix.Dead.Tot,model=Mean.Post))
    
} # End ayer.dbh.comp function

#===============================================================================================================
#===============================================================================================================


