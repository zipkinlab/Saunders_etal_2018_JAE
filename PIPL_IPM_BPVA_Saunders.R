########################################################################################
# Integrated population model (IPM) - Bayesian population viability analysis (BPVA) 
# for Great Lakes piping plovers, 1993 - 2016
# Sarah Saunders, Francesca Cuthbert, Elise Zipkin

# Adapted from original scripts by Marc Kéry & Michael Schaub (2016)
# Modified by S. Saunders, 2016 - 2017

########################################################################################

# Load data and libraries
library(jagsUI)

nyears <- 24	  # Number of years in analysis

#Load function to create a m-array based on capture-recapture data (CH)
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

########################################################################
# Capture-recapture data: m-array of juveniles (HY) and adults (AHY)
########################################################################

#First read in capture histories for birds marked as HY during 1993-2016
CH.J <- read.table("CH_HYmark16.txt")

#convert to matrix
CH.J <- data.matrix(CH.J)

#read in capture histories for birds marked as AHY during 1993-2016
CH.A <- read.table("CH_AHYmark16.txt")

#convert to matrix
CH.A <- data.matrix(CH.A)

#create two m-arrays, one for juveniles and one for adults
cap <- apply(CH.J, 1, sum)
ind <- which(cap >= 2)
CH.J.R <- CH.J[ind,] # Juvenile CH recaptured at least once
CH.J.N <- CH.J[-ind,] # Juvenile CH never recaptured

# Remove first capture
first <- numeric()
for (i in 1:dim(CH.J.R)[1]){
  first[i] <- min(which(CH.J.R[i,]==1))
}
CH.J.R1 <- CH.J.R
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R1[i,first[i]] <- 0
}

# Add grown-up juveniles to adults and create m-array
CH.A.m <- rbind(CH.A, CH.J.R1)
CH.A.marray <- marray(CH.A.m)

# Create CH matrix for juveniles, ignoring subsequent recaptures
second <- numeric()
for (i in 1:dim(CH.J.R1)[1]){
  second[i] <- min(which(CH.J.R1[i,]==1))
}
CH.J.R2 <- matrix(0, nrow = dim(CH.J.R)[1], ncol = dim(CH.J.R)[2])
for (i in 1:dim(CH.J.R)[1]){
  CH.J.R2[i,first[i]] <- 1
  CH.J.R2[i,second[i]] <- 1
}

# Create m-array for these
CH.J.R.marray <- marray(CH.J.R2)

# The last column should show the number of juveniles not recaptured again and should all be zeros, since all of them are released as adults
CH.J.R.marray[,dim(CH.J)[2]] <- 0

# Create the m-array for juveniles never recaptured and add it to the previous m-array
CH.J.N.marray <- marray(CH.J.N)
CH.J.marray <- CH.J.R.marray + CH.J.N.marray

#outputs: CH.A.marray and CH.J.marray
#convert outputs to names of m-arrays used in models

marray.j <- CH.J.marray
marray.a <- CH.A.marray

# Population count data, nesting PIPL pairs (1993-2016)
y <-  c(18,19,21,24,23,23,32,30,32,51,50,55,58,53,63,63,71,60,55,58,66,70,75,75)

# Productivity data (1993-2016)
J <- c(13,28,42,26,39,39,49,40,71,61,88,92,93,94,124,113,126,93,75,121,124,109,128,133) # Number of offspring/fledglings
R <- c(18,19,21,23,23,23,32,30,31,50,49,52,56,53,61,60,69,59,54,57,66,70,74,74) # Number of surveyed broods/brdg pairs contributing data

#########################################
# Specify model in BUGS language
#########################################

sink("imm.merlin.ipm.pva.jags")
cat("
    model {
    #----------------------------------------------------------------------------
    #  Integrated population model BPVA (10 yr predictions)
    #  - Stage structured model with 2 stages: juvenile and adult
    #  - Age at first breeding = 1 year
    #  - Prebreeding census, female-based
    #  - All vital rates assumed to be time-dependent (random env. stochasticity)
    #  - Includes env. stochasticity thru random time effects for all params
    #  - Explicit estimation of immigration as expected number of individuals
    #  - Merlin effect (latent abundance) on adult survival estimated by state-
    # space model
    #----------------------------------------------------------------------------
    
    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes
    n1 ~ dnorm(100, 0.001)I(0,)           # HY individuals
    nadSurv ~ dnorm(100, 0.001)I(0,)      # Adults >= 2 years
    nadimm ~ dnorm(100, 0.001)I(0,)       # Immigrants
    
    N1[1] <- round(n1)
    NadSurv[1] <- round(nadSurv)
    Nadimm[1] <- round(nadimm)
    Ntot[1] <- N1[1] + NadSurv[1] + Nadimm[1]
    
    # Mean demographic parameters (on appropriate scale)
    l.mphij ~ dnorm(0, 0.001)              
    l.mphia ~ dnorm(0, 0.001)
    l.mfec ~ dnorm(0, 0.001)
    b0.omm ~ dunif(0, 20)                 #expected number of immigrants              
    l.p ~ dnorm(0, 0.001)
    beta.phia ~ dnorm(0, 0.1)               
    
    #back transformation
    log.b0.omm <- log(b0.omm)
    
    
    # Precision of standard deviations of temporal variability
    sig.phij ~ dunif(0, 10)
    tau.phij <- pow(sig.phij, -2)
    sig.phia ~ dunif(0, 10)                  
    tau.phia <- pow(sig.phia, -2)
    sig.fec ~ dunif(0, 10)
    tau.fec <- pow(sig.fec, -2)
    sig.im ~ dunif(0, 10)
    tau.im <- pow(sig.im, -2)
    
    sig.obs ~ dunif(0.5, 50)
    tau.obs <- pow(sig.obs, -2)
    
    # Distribution of error terms (Bounded to help with convergence)
    for (t in 1:(nyears-1+K)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)	
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.fec[t] ~ dnorm(0, tau.fec)T(-5,5)
    epsilon.im[t] ~ dnorm(0, tau.im)T(-5,5)   
    }
    
    #-----------------------------------------------
    # 2. Constrain parameters (for temp variability)
    #-----------------------------------------------
    
    # Juv. apparent survival
    for (t in 1:(nyears-1+K)){
    logit(phij[t]) <- l.mphij + epsilon.phij[t]          
    
    # Adult apparent survival with merlin effect
    logit(phia[t]) <- l.mphia + beta.phia*N.cor[t] + epsilon.phia[t]    
    
    log(f[t]) <- l.mfec + epsilon.fec[t]                           # Productivity
    log(omega[t]) <- log.b0.omm + epsilon.im[t]                     # Immigration
    logit(p[t]) <- l.p                                    # Recapture probability
    }
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    
    mphij <- exp(l.mphij)/(1+exp(l.mphij))   # Mean juvenile survival probability
    mphia <- exp(l.mphia)/(1+exp(l.mphia))   # Mean adult survival probability
    mfec <- exp(l.mfec)                      # Mean productivity
    
    # Population growth rate (1993 to 2016)
    for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / (Ntot[t] + 0.0001)                             
    logla[t] <- log(lambda[t])
    }
    mlam <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)]))   # Geometric mean
    
    
    # Population growth rate (merlins)
    for (t in 1:(nyears-1+K)){
    lambda.mer[t] <- N.est[t+1] / (N.est[t] + 0.0001)    
    logla.mer[t] <- log(lambda.mer[t])
    }
    
    mlam.mer <- exp((1/(nyears-1+K))*sum(logla.mer[1:(nyears-1+K)])) 
    # Geometric mean (all years)
    
    mlampast.mer <- exp((1/(nyears-1))*sum(logla.mer[1:(nyears-1)]))   	
    # Geometric mean for 1993-2015
    
    mlamfut.mer <- exp((1/(K-1))*sum(logla.mer[nyears:(nyears-1+K)]))   	  
    # Geometric mean for 2016-2026 
    
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    for (t in 2:(nyears+K)){
    Ntot[t] <- NadSurv[t] + N1[t] + Nadimm[t]
    mean1[t] <- 0.5 * f[t-1] * phij[t-1] * Ntot[t-1]        
    N1[t] ~ dpois(mean1[t])
    NadSurv[t] ~ dbin(phia[t-1], Ntot[t-1])
    Nadimm[t] ~ dpois(omega[t-1])
    }
    
    # 4.1.2 Observation process
    for (t in 1:nyears){
    y[t] ~ dnorm(Ntot[t], tau.obs) 
    }
    
    # 4.2 Likelihood for capture-recapture data: CJS model (2 age classes)
    # Multinomial likelihood
    for (t in 1:(nyears-1)){
    marray.j[t,1:nyears] ~ dmulti(pr.j[t,], r.j[t])
    marray.a[t,1:nyears] ~ dmulti(pr.a[t,], r.a[t])
    }
    
    # m-array cell probabilities for juveniles
    for (t in 1:(nyears-1)){
    q[t] <- 1-p[t]
    # Main diagonal
    pr.j[t,t] <- phij[t]*p[t]
    # Above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.j[t,j] <- phij[t]*prod(phia[(t+1):j])*prod(q[t:(j-1)])*p[j]
    } #j
    
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.j[t,j] <- 0
    } #j
    # Last column
    pr.j[t,nyears] <- 1-sum(pr.j[t,1:(nyears-1)])
    } #t
    
    # m-array cell probabilities for adults
    for (t in 1:(nyears-1)){
    # Main diagonal
    pr.a[t,t] <- phia[t]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    
    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:(nyears-1)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t] * f[t]
    }
    
    #-------------------------------------------------------------------
    #  5. State-space model for merlin index (effect on adult survival)
    #-------------------------------------------------------------------
    
    # Priors and contraints
    logN.est[1] ~ dnorm(4.4, 0.01)            # Prior for initial population size
    
    mean.r ~ dnorm(0, 0.01)                           # Prior for mean grown rate
    
    sigma.proc ~ dunif(0, 1)                      # Prior for SD of state process
    sigma2.proc <- pow(sigma.proc, 2)
    tau.proc <- pow(sigma.proc, -2)
    
    sigma.obs ~ dunif(0, 1)                        # Prior for SD of obs. process
    sigma2.obs <- pow(sigma.obs, 2)
    t.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){                  # T is 34 years (24 + 10 prediction yrs)
    r[t] ~ dnorm(mean.r, tau.proc)
    logN.est[t+1] <- logN.est[t] + r[t]
    }
    
    # Observation process                     
    for (t in 1:T) {
    for (s in 1:S){
    x[t,s] ~ dnorm(logN.est[t], t.obs)
    }
    }
    
    # Population sizes on real scale
    for (t in 1:T) {                           
    N.est[t] <- exp(logN.est[t])
    N.cor[t] <- (N.est[t]-N.mean)/N.sd           #standardize to use as covariate
    }
    }
    ",fill = TRUE)
    sink()
    
###################################################################
    
# Load data----------------------------------------------------------
    
M <- read.table("merlins.txt",header=TRUE)    #Hawk Mtn. and Whitefish Pt. counts
    
#First, alter data input
hawk <- c(M$HM)           #  Hawk Mountain counts
white <- c(M$WP)          #  Whitefish Point counts
mat <- matrix(c(hawk, white), nrow=length(hawk))
    
#---------------------------------------------------------------------
    
# Bundle data
K <- 10                           # Number of years with predictions
nyears <- ncol(marray.j)          # Number of study years
N.mean = 114.9                    # Mean estimate of merlin population size
N.sd = 17.3                       # SD of merlin population size
    
#adjust merlin matrix for prediction years with NAs
v1 <- c(rep(NA, K))
v2 <- c(rep(NA, K))
mat.add <- matrix(c(v1, v2), nrow=length(v1))
mat.proj <-rbind(mat, mat.add)
    
jags.data <- list(nyears = nyears, marray.j = marray.j, marray.a = marray.a, y = y, J = J, R = R, r.j = rowSums(marray.j), r.a = rowSums(marray.a), K = K, x = log(mat.proj), T = nrow(mat.proj), S = ncol(mat.proj), N.mean = N.mean, N.sd = N.sd)
    
# Initial values
inits <- function(){list(l.mphij = rnorm(1, 0.2, 0.5), l.mphia = rnorm(1, 0.2, 0.5), l.mfec = rnorm(1, 0.2, 0.5), l.p = rnorm(1, 0.2, 1), sig.phij = runif(1, 0.1, 10), sig.phia = runif(1, 0.1, 10), sig.fec = runif(1, 0.1, 10),  n1 = round(runif(1, 1, 50), 0), nadSurv = round(runif(1, 5, 50), 0), beta.phia = runif(1, -1, 1), b0.omm = runif(1, 0, 10), sig.im = runif(1, 0.1, 10), nadimm = round(runif(1, 1, 50), 0), sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), sigma.obs = runif(1, 0, 1),logN.est = c(rnorm(1, 4.4, 0.1), rep(NA, (nrow(mat.proj) - 1))))}
    
# Parameters monitored
parameters <- c("phij", "phia", "f",  "p", "lambda", "mphij", "mphia", "mfec", "mlam", "mlam.mer", "mlampast.mer", "mlamfut.mer", "beta.phia", "sig.phij", "sig.phia", "sig.fec", "sig.obs", "omega", "sig.im", "N1", "NadSurv", "Ntot", "Nadimm", "b0.omm", "r", "mean.r", "sigma2.obs", "sigma2.proc", "N.cor", "N.est")
    
# MCMC settings 
ni <- 400000    
nt <- 10
nb <- 200000
nc <- 3
    
# Call JAGS from R (jagsUI)
ipm.pva <- jags(jags.data, inits, parameters, "imm.merlin.ipm.pva.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE, store.data=TRUE)

####################################################
#Plots of predicted pop sizes and demographic rates
###################################################

m1 <- min(c(y, ipm.pva$q2.5$Ntot))
m2 <- max(c(y, ipm.pva$q97.5$Ntot))
n.years <- length(ipm.pva$mean$Ntot)

par(mfrow = c(2, 3), mar = c(5, 5, 1, 1))
plot(0, 0, ylim = c(0, 200), xlim = c(0.5, n.years), ylab = "Piping plover population size", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
axis(2, las = 1)
axis(1, at = seq(1, n.years, 3), labels = seq(1, n.years, 3))
axis(1, at = 1:n.years, labels = rep("", n.years), tcl = -0.25)
polygon(x = c(1:n.years, n.years:1), y = c(ipm.pva$q2.5$Ntot, ipm.pva$q97.5$Ntot[n.years:1]), col = "gray90", border = "gray90")
points(y, type = "l", col = "black", lwd = 2)
points(ipm.pva$mean$Ntot, type = "l", col = "darkgoldenrod1", lwd = 2)
legend(1,200, legend = c("Observed", "Estimated"), lty = c(1, 1), lwd = c(2, 2), col = c("black", "darkgoldenrod1"), bty = "n", cex = 1)

plot(ipm.pva$mean$phij, ylim = range(c(0, 0.45)), xlim = c(0.5, n.years), ylab = "Juvenile survival", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
axis(2, las = 1)
axis(1, at = seq(1, n.years-1, 3), labels = seq(1, n.years-1, 3))
axis(1, at = 1:(n.years-1), labels = rep("", n.years-1), tcl = -0.25)
polygon(x = c(1:(n.years-1), (n.years-1):1), y = c(ipm.pva$q2.5$phij, ipm.pva$q97.5$phij [(n.years-1):1]), col = "gray90", border = "gray90")
points(ipm.pva$mean$phij, type = "l", col = "darkgoldenrod1", lwd = 2)

plot(ipm.pva$mean$f, ylim = range(c(0.5, 3.0)), xlim = c(0.5, n.years), ylab = "Fecundity", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
axis(2, las = 1)
axis(1, at = seq(1, n.years-1, 3), labels = seq(1, n.years-1, 3))
axis(1, at = 1:(n.years-1), labels = rep("", n.years-1), tcl = -0.25)
polygon(x = c(1:(n.years-1), (n.years-1):1), y = c(ipm.pva$q2.5$f, ipm.pva$q97.5$f [(n.years-1):1]), col = "gray90", border = "gray90")
points(ipm.pva$mean$f, type = "l", col = "darkgoldenrod1", lwd = 2)

plot(ipm.pva$mean$omega, ylim = range(c(ipm.pva$q2.5$omega, ipm.pva$q97.5$omega)), xlim = c(0.5, n.years), ylab = "Immigration (number indivs.)", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
axis(2, las = 1)
axis(1, at = seq(1, n.years-1, 3), labels = seq(1, n.years-1, 3))
axis(1, at = 1:(n.years-1), labels = rep("", n.years-1), tcl = -0.25)
polygon(x = c(1:(n.years-1), (n.years-1):1), y = c(ipm.pva$q2.5$omega, ipm.pva$q97.5$omega [(n.years-1):1]), col = "gray90", border = "gray90")
points(ipm.pva$mean$omega, type = "l", col = "darkgoldenrod1", lwd = 2)

plot(ipm.pva$mean$phia, ylim = range(c(0, 1.0)), xlim = c(0.5, n.years), ylab = "Adult survival", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
axis(2, las = 1)
axis(1, at = seq(1, n.years-1, 3), labels = seq(1, n.years-1, 3))
axis(1, at = 1:(n.years-1), labels = rep("", n.years-1), tcl = -0.25)
polygon(x = c(1:(n.years-1), (n.years-1):1), y = c(ipm.pva$q2.5$phia, ipm.pva$q97.5$phia [(n.years-1):1]), col = "gray90", border = "gray90")
points(ipm.pva$mean$phia, type = "l", col = "darkgoldenrod1", lwd = 2)

#adding merlin pop growth panel
plot(ipm.pva$mean$N.est, ylim = range(c(0, 250)), xlim = c(0.5, n.years), ylab = "Merlin population size", xlab = "Year", las = 1, col = "darkgreen", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
axis(2, las = 1)
axis(1, at = seq(1, n.years, 3), labels = seq(1, n.years, 3))
axis(1, at = 1:n.years, labels = rep("", n.years), tcl = -0.25)
polygon(x = c(1:n.years, n.years:1), y = c(ipm.pva$q2.5$N.est, ipm.pva$q97.5$N.est[n.years:1]), col = "gray90", border = "gray90")
points(ipm.pva$mean$N.est, type = "l", col = "forestgreen", lwd = 2)

#how often estimates above 78 in 2026?
mean(ipm.pva$sims.list$Ntot[,n.years]>78)

#below 72?
mean(ipm.pva$sims.list$Ntot[,n.years]<72)

#within CI of 72 to 78?
mean(ipm.pva$sims.list$Ntot[,n.years]>=72 & (ipm.pva$sims.list$Ntot[,n.years]<=78))

