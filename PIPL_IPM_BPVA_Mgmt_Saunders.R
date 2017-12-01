########################################################################################
# Integrated population model (IPM) - Bayesian population viability analysis (BPVA) 
# for Great Lakes piping plovers, 1993 - 2016
# Evaluating 4 management scenarios: merlin control, chick predator control, simultaneous
# control, and null scenario
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

sink("imm.merlin.ipm.pvaII.jags")
cat("
    model {
    #-------------------------------------------------------------------------
    #  Integrated population model BPVA (10 yr predictions)
    #  - Stage structured model with 2 stages: juvenile and adult
    #  - Age at first breeding = 1 year
    #  - Prebreeding census, female-based
    #  - All vital rates assumed to be time-dependent (random env. stochasticity)
    #  - Includes env. stochasticity thru random time effects for all params
    #  - Explicit estimation of immigration as expected number of individuals
    #  - Merlin effect on adult survival estimated by state-space model
    #  - Four management scenarios evaluated for efficacy
    #-------------------------------------------------------------------------
    
    #----------------------------------------
    # 1. Define the priors for the parameters
    #----------------------------------------
    
    # Initial population sizes
    n1 ~ dnorm(100, 0.001)I(0,)           # HY individuals
    nadSurv ~ dnorm(100, 0.001)I(0,)      # Adults >= 2 years
    nadimm ~ dnorm(100, 0.001)I(0,)       # Immigrants
    
    N1[1,1] <- round(n1)
    NadSurv[1,1] <- round(nadSurv)
    Nadimm[1,1] <- round(nadimm)
    Ntot[1,1] <- N1[1,1] + NadSurv[1,1] + Nadimm[1,1]
    
    # Mean demographic parameters (on appropriate scale)
    l.mphij ~ dnorm(0, 0.001)                
    l.mphia ~ dnorm(0, 0.001)
    l.mfec ~ dnorm(0, 0.001)
    b0.omm ~ dunif(0, 20)                     #expected number of immigrants              
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
    epsilon.phia[t,1] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.im[t] ~ dnorm(0, tau.im)T(-5,5)   
    }
    
    for (t in 1:(nyears+K)){
    epsilon.fec[t,1] ~ dnorm(0, tau.fec)T(-5,5)
    }
    
    #------------------------------------------------
    # 2. Constrain parameters (for temp variability)
    #------------------------------------------------
    
    #Scenario 1: no change    
    
    # Juvenile apparent survival
    for (t in 1:(nyears-1+K)){
    logit(phij[t]) <- l.mphij + epsilon.phij[t]          
    
    # Adult apparent survival with merlin effect
    logit(phia[t,1]) <- l.mphia + beta.phia*N.cor[t] + epsilon.phia[t,1]       
    
    log(omega[t]) <- log.b0.omm + epsilon.im[t]                     # Immigration
    logit(p[t]) <- l.p                                    # Recapture probability
    }
    
    for (t in 1:(nyears+K)){
    log(f[t,1]) <- l.mfec + epsilon.fec[t,1]                       # Productivity
    }
    
    #Scenario 2: increase of productivity
    
    for (t in 1:nyears){                          # Past: identical to scenario 1
    log(f[t,2]) <- log(f[t,1])
    epsilon.fec[t,2] <- epsilon.fec[t,1]
    }
    
    
    for (t in (nyears+1):(nyears+K)){                   # Future: increase by 20%
    log(f[t,2]) <- l.mfec + log(1.2) + epsilon.fec[t,2]
    epsilon.fec[t,2] ~ dnorm(0, tau.fec)T(-5,5)
    }
    
    #Scenario 3: reduction of mean merlin abundance on adult survival
    
    for (t in 1:(nyears-1)){
    logit(phia[t,2]) <- logit(phia[t,1])           # Past identical to scenario 1
    epsilon.phia[t,2] <- epsilon.phia[t,1]
    }
    
    # Future: new merlin covariate (20% fewer/yr)
    for (t in nyears:(nyears-1+K)){  
    logit(phia[t,2]) <- l.mphia + beta.phia*N.cor.new[t] + epsilon.phia[t,2]
    epsilon.phia[t,2] ~ dnorm(0, tau.phia)
    }
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    
    for (t in 1:(nyears+K)){
    N.tot[t,1] <- NadSurv[t,1] + N1[t,1] + Nadimm[t,1]  #Total population sizes
    N.tot[t,2] <- NadSurv[t,2] + N1[t,2] + Nadimm[t,2]  #Total population sizes
    N.tot[t,3] <- NadSurv[t,3] + N1[t,3] + Nadimm[t,3]  #Total population sizes
    N.tot[t,4] <- NadSurv[t,4] + N1[t,4] + Nadimm[t,4]  #Total population sizes
    }    
    
    mphij <- exp(l.mphij)/(1+exp(l.mphij))   # Mean juvenile survival probability
    mphia <- exp(l.mphia)/(1+exp(l.mphia))   # Mean adult survival probability
    mfec <- exp(l.mfec)                      # Mean productivity
    
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    
    # 4.1. Likelihood for population count data (state-space model)
    # 4.1.1 System process
    #Scenario 1: no change
    
    for (t in 2:(nyears+K)){
    mean1[t,1] <- 0.5 * f[t-1,1] * phij[t-1] * (NadSurv[t-1,1] + N1[t-1,1] +
    Nadimm[t-1,1])        
    N1[t,1] ~ dpois(mean1[t,1])
    NadSurv[t,1] ~ dbin(phia[t-1,1],(NadSurv[t-1,1] + N1[t-1,1] + Nadimm[t-1,1]))
    Nadimm[t,1] ~ dpois(omega[t-1])
    }
    
    #Scenario 2: increase of productivity
    
    # Past (same as scenario 1)
    for (t in 1:nyears){
    N1[t,2] <- N1[t,1]
    NadSurv[t,2] <- NadSurv[t,1]
    Nadimm[t,2] <- Nadimm[t,1]
    }
    
    #Future
    #use different value for productivity f[t,2]
    for (t in (nyears+1):(nyears+K)){
    N1[t,2] ~ dpois(0.5 * f[t-1,2] * phij[t-1] * (NadSurv[t-1,2] + N1[t-1,2] +
    Nadimm[t-1,2]))
    NadSurv[t,2] ~ dbin(phia[t-1,1],(NadSurv[t-1,2] + N1[t-1,2] + Nadimm[t-1,2]))
    Nadimm[t,2] ~ dpois(omega[t-1])
    }
    
    #Scenario 3: decrease mean merlin abundance
    
    #Past
    for (t in 1:nyears){
    N1[t,3] <- N1[t,1]
    NadSurv[t,3] <- NadSurv[t,1]
    Nadimm[t,3] <- Nadimm[t,1]
    }
    
    #Future
    #still use productivity from scenario 1
    #for adult survival, use a different value phia[t,2]
    for (t in (nyears+1):(nyears+K)){
    N1[t,3] ~ dpois(0.5 * f[t-1,1] * phij[t-1] * (NadSurv[t-1,3] + N1[t-1,3] +
    Nadimm[t-1,3]))
    NadSurv[t,3] ~ dbin(phia[t-1,2],(NadSurv[t-1,3] + N1[t-1,3] + Nadimm[t-1,3]))
    Nadimm[t,3] ~ dpois(omega[t-1])
    }
    
    #Scenario 4: combine both scenarios
    
    #Past
    for (t in 1:nyears){
    N1[t,4] <- N1[t,1]
    NadSurv[t,4] <- NadSurv[t,1]
    Nadimm[t,4] <- Nadimm[t,1]
    }
    
    #Future
    # use productivity from scenario 2
    # use adult survival from scenario 3
    for (t in (nyears+1):(nyears+K)){
    N1[t,4] ~ dpois(0.5 * f[t-1,2] * phij[t-1] * (NadSurv[t-1,4] + N1[t-1,4] +
    Nadimm[t-1,4]))
    NadSurv[t,4] ~ dbin(phia[t-1,2],(NadSurv[t-1,4] + N1[t-1,4] + Nadimm[t-1,4]))
    Nadimm[t,4] ~ dpois(omega[t-1])
    }
    
    # 4.1.2 Observation process 
    for (t in 1:nyears){
    y[t] ~ dnorm(NadSurv[t,1] + N1[t,1] + Nadimm[t,1], tau.obs)
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
    pr.j[t,j] <- phij[t]*prod(phia[(t+1):j,1])*prod(q[t:(j-1)])*p[j]
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
    pr.a[t,t] <- phia[t,1]*p[t]
    # above main diagonal
    for (j in (t+1):(nyears-1)){
    pr.a[t,j] <- prod(phia[t:j,1])*prod(q[t:(j-1)])*p[j]
    } #j
    # Below main diagonal
    for (j in 1:(t-1)){
    pr.a[t,j] <- 0
    } #j
    # Last column
    pr.a[t,nyears] <- 1-sum(pr.a[t,1:(nyears-1)])
    } #t
    
    # 4.3. Likelihood for productivity data: Poisson regression
    for (t in 1:(nyears)){
    J[t] ~ dpois(rho[t])
    rho[t] <- R[t] * f[t,1]
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
    N.cor[t] <- (N.est[t]-N.mean)/N.sd
    }
    
    for (t in nyears:T){         # new merlin abundance for future (20% fewer/yr)
    N.est.new[t] <- (N.est[t] - (N.est[t]*0.2))
    N.cor.new[t] <- (N.est.new[t]-N.mean)/N.sd                  # standardize
    }    
    }
    ",fill = TRUE)
sink()

###################################################################

# Load data--------------------------------------------------------

M <- read.table("merlins.txt",header=TRUE)   # Hawk Mtn. and Whitefish Pt. counts
# First, alter data input
hawk <- c(M$HM)
white <- c(M$WP)
mat <- matrix(c(hawk, white), nrow=length(hawk))

#---------------------------------------------------------------------

# Bundle data
K <- 10                           # Number of years with predictions
nyears <- ncol(marray.j)          # Number of study years
N.mean = 114.9
N.sd = 17.3

#adjust merlin matrix for prediction years with NAs
v1 <- c(rep(NA, K))
v2 <- c(rep(NA, K))
mat.add <- matrix(c(v1, v2), nrow=length(v1))
mat.proj <-rbind(mat, mat.add)

jags.data <- list(nyears = nyears, marray.j = marray.j, marray.a = marray.a, y = y, J = J, R = R, r.j = rowSums(marray.j), r.a = rowSums(marray.a), K = K, x = log(mat.proj), T = nrow(mat.proj), S = ncol(mat.proj), N.mean = N.mean, N.sd = N.sd)

# Initial values
inits <- function(){list(l.mphij = rnorm(1, 0.2, 0.5), l.mphia = rnorm(1, 0.2, 0.5), l.mfec = rnorm(1, 0.2, 0.5), l.p = rnorm(1, 0.2, 1), sig.phij = runif(1, 0.1, 10), sig.phia = runif(1, 0.1, 10), sig.fec = runif(1, 0.1, 10),  n1 = round(runif(1, 1, 50), 0), nadSurv = round(runif(1, 5, 50), 0), beta.phia = runif(1, -1, 1), b0.omm = runif(1, 0, 10), sig.im = runif(1, 0.1, 10), nadimm = round(runif(1, 1, 50), 0), sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), sigma.obs = runif(1, 0, 1),logN.est = c(rnorm(1, 4.4, 0.1), rep(NA, (nrow(mat.proj) - 1))))}

# Parameters monitored
parameters <- c("phij", "phia", "f", "p", "mphij", "mphia", "mfec", "beta.phia", "sig.phij", "sig.phia", "sig.fec", "sig.obs", "omega", "sig.im", "N1", "NadSurv", "N.tot", "Nadimm", "b0.omm", "r", "mean.r", "sigma2.obs", "sigma2.proc", "N.cor", "N.est")

# MCMC settings 
ni <- 400000    
nt <- 10
nb <- 200000
nc <- 3

# Call JAGS from R (jagsUI)
ipm.pva_Mgmt <- jags(jags.data, inits, parameters, " imm.merlin.ipm.pvaII.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel=TRUE, store.data=TRUE)

########################################################################
# Cumulative extinction probabilities and management scenarios graphs
#######################################################################

m1 <- min(c(y, ipm.pva_Mgmt$q2.5$N.tot))
m2 <- max(c(y, ipm.pva_Mgmt$q97.5$N.tot))
start <- nyears
end <- nyears + K

ext.threshold <- 15                         
extinct <- array(NA, dim = dim(ipm.pva_Mgmt$sims.list$N.tot))
extinct[,,1] <- ipm.pva_Mgmt$sims.list$N.tot[,,1] <= ext.threshold
extinct[,,2] <- ipm.pva_Mgmt$sims.list$N.tot[,,2] <= ext.threshold
extinct[,,3] <- ipm.pva_Mgmt$sims.list$N.tot[,,3] <= ext.threshold
extinct[,,4] <- ipm.pva_Mgmt$sims.list$N.tot[,,4] <= ext.threshold

co <- colorRampPalette(c("blue", "green"))(4)
n.years <- dim(extinct)[2]
par(mfrow=c(1,1))
plot(apply(extinct[,,1], 2, mean)[(nyears+1):n.years], type = "l", ylab = "Quasi-extinction probability", lwd = 2, xlab = "Year", frame = FALSE, axes = FALSE, col = co[1])
axis(1, at = 1:K, tck = -0.0125, labels = FALSE)
axis(1, at = c(1, 3, 5, 7, 9, 11, 13, 15), labels = c(1, 3, 5, 7, 9, 11, 13, 15), tck = -0.025)
axis(2, at = c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12), labels = c(0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12))
lines(apply(extinct[,,2], 2, mean)[(nyears+1):n.years], type = "l", ylab = "Quasi-extinction probability", lwd = 2, xlab = "Year", col = co[2])
lines(apply(extinct[,,3], 2, mean)[(nyears+1):n.years], type = "l", ylab = "Quasi-extinction probability", lwd = 2, xlab = "Year", col = co[3])
lines(apply(extinct[,,4], 2, mean)[(nyears+1):n.years], type = "l", ylab = "Quasi-extinction probability", lwd = 2, xlab = "Year", col = co[4])
legend("topleft", lty = rep(1, 4), lwd = rep(2, 4), col = co, legend = c("Null scenario", "Increased productivity", "Increased adult survival", "Both scenarios"), bty = "n")

#Given estimates are for each scenario the predicted population size and the extinction probability in K = 10 years (year 34).
ipm.pva_Mgmt$mean$N.tot[K+nyears,]
ipm.pva_Mgmt$q2.5$N.tot[K+nyears,]
ipm.pva_Mgmt$q97.5$N.tot[K+nyears,]
mean(extinct[,nyears+K,1])
mean(extinct[,nyears+K,2])
mean(extinct[,nyears+K,3])
mean(extinct[,nyears+K,4])

#To get further insight into the options to take, we compute the probability that the scenarios with management result in a larger 
#population size in K = 10 years compared to when no management is taken

mean(ipm.pva_Mgmt$sims.list$N.tot[,nyears+K,2] > ipm.pva_Mgmt$sims.list$N.tot[,nyears+K,1])
mean(ipm.pva_Mgmt$sims.list$N.tot[,nyears+K,3] > ipm.pva_Mgmt$sims.list$N.tot[,nyears+K,1])
mean(ipm.pva_Mgmt$sims.list$N.tot[,nyears+K,4] > ipm.pva_Mgmt$sims.list$N.tot[,nyears+K,1])

#Scenario 2
#how often estimates above 78 in 2026?
mean(ipm.pva_Mgmt$sims.list$N.tot[,n.years,2] > 78)

#below 72?
mean(ipm.pva_Mgmt$sims.list$N.tot[,n.years,2] < 72) 

#within CI of 72 to 78?
mean(ipm.pva_Mgmt$sims.list$N.tot[,n.years,2]>=72 & (ipm.pva_Mgmt$sims.list$N.tot[,n.years,2]<=78)) 
##############################################################################################
#scenario 3
#how often estimates above 78 in 2026?
mean(ipm.pva_Mgmt$sims.list$N.tot[,n.years,3]>78) 

#below 72?
mean(ipm.pva_Mgmt$sims.list$N.tot[,n.years,3]<72) 

#within CI of 72 to 78?
mean(ipm.pva_Mgmt$sims.list$N.tot[,n.years,3]>=72 & (ipm.pva_Mgmt$sims.list$N.tot[,n.years,3]<=78))
#################################################################################################
#scenario 4
#how often estimates above 78 in 2026?
mean(ipm.pva_Mgmt$sims.list$N.tot[,n.years,4]>78) 

#below 72?
mean(ipm.pva_Mgmt$sims.list$N.tot[,n.years,4]<72) 

#within CI of 72 to 78?
mean(ipm.pva_Mgmt$sims.list$N.tot[,n.years,4]>=72 & (ipm.pva_Mgmt$sims.list$N.tot[,n.years,4]<=78))
