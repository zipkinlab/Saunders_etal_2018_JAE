########################################################################################
# Integrated population model (IPM) for Great Lakes piping plovers, 1993 - 2016
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

#############################
#Merge juv and adult m-arrays 
#to create a single m-array (m)
############################

m <- rbind(CH.J.marray, CH.A.marray)

# Population count data, nesting PIPL pairs (1993-2016)
y <-  c(18,19,21,24,23,23,32,30,32,51,50,55,58,53,63,63,71,60,55,58,66,70,75,75)

# Productivity data (1993-2015)
J <- c(13,28,42,26,39,39,49,40,71,61,88,92,93,94,124,113,126,93,75,121,124,109,128) # Number of offspring/fledglings
R <- c(18,19,21,23,23,23,32,30,31,50,49,52,56,53,61,60,69,59,54,57,66,70,74) # Number of surveyed broods/brdg pairs contributing data

#########################################
# Specify model in BUGS language
#######################################

sink("pipl.ipm.merlin.jags")
cat("
    model {
    #-----------------------------------------------------------------------------------
    #  Integrated population model
    #  - Age structured model with 2 age classes: 
    #		HY and AHY
    #  - Age at first breeding = 1 year
    #  - Prebreeding census, female-based
    #  - All vital rates are assumed to be time-dependent
    #  - Includes env. stochasticity thru random time effects for all params
    #  - Explicit estimation of immigration
    #  - Merlin effect on adult survival only via state space model
    #-----------------------------------------------------------------------------------
    
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
    
    # Mean demographic parameters (on appropriate scale)
    l.mphij ~ dnorm(0, 0.001)              
    l.mphia ~ dnorm(0, 0.001)
    l.mfec ~ dnorm(0, 0.001)
    b0.omm ~ dunif(0, 20)                  #expected number of immigrants               
    l.p ~ dnorm(0, 0.001)
    beta.phia ~ dnorm(0, 0.1)              #uninformative prior for merlin effect on adult survival
    
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
    for (t in 1:(nyears-1)){
    epsilon.phij[t] ~ dnorm(0, tau.phij)T(-5,5)	
    epsilon.phia[t] ~ dnorm(0, tau.phia)T(-5,5)
    epsilon.fec[t] ~ dnorm(0, tau.fec)T(-5,5)
    epsilon.im[t] ~ dnorm(0, tau.im)T(-5,5)
    }
    
    #-------------------------
    # 2. Constrain parameters
    #-------------------------
    for (t in 1:(nyears-1)){
    logit(phij[t]) <- l.mphij + epsilon.phij[t]                            # Juv. apparent survival
    logit(phia[t]) <- l.mphia + beta.phia*N.cor[t] + epsilon.phia[t]       # Adult apparent survival
    log(f[t]) <- l.mfec + epsilon.fec[t]                                   # Productivity
    log(omega[t]) <- log.b0.omm + epsilon.im[t]                            # Immigration
    logit(p[t]) <- l.p                                                     # Recapture probability
    }
    
    #-----------------------
    # 3. Derived parameters
    #-----------------------
    mphij <- exp(l.mphij)/(1+exp(l.mphij))   # Mean juvenile survival probability
    mphia <- exp(l.mphia)/(1+exp(l.mphia))   # Mean adult survival probability
    mfec <- exp(l.mfec)                      # Mean productivity
    
    # Population growth rate
    for (t in 1:(nyears-1)){
    lambda[t] <- Ntot[t+1] / (Ntot[t] + 0.0001)                             
    logla[t] <- log(lambda[t])
    imrate[t] <- Nadimm[t+1] / Ntot[t]                     # Derived immigration rate
    }
    mlam <- exp((1/(nyears-1))*sum(logla[1:(nyears-1)]))   # Geometric mean
    
    #--------------------------------------------
    # 4. The likelihoods of the single data sets
    #--------------------------------------------
    # 4.1. Likelihood for population population count data (state-space model)
    # 4.1.1 System process
    for (t in 2:nyears){
    mean1[t] <- 0.5 * f[t-1] * phij[t-1] * Ntot[t-1]
    N1[t] ~ dpois(mean1[t])
    NadSurv[t] ~ dbin(phia[t-1], Ntot[t-1])
    Nadimm[t] ~ dpois(omega[t-1])
    }
    
    # 4.1.2 Observation process
    for (t in 1:nyears){
    Ntot[t] <- NadSurv[t] + N1[t] + Nadimm[t]
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
    #  5. State space model for merlin counts (effect on adult survival)
    #-------------------------------------------------------------------
    # Priors and contraints
    logN.est[1] ~ dnorm(4.4, 0.01)              # Prior for inital population size
    
    mean.r ~ dnorm(0, 0.01)                     # Prior for mean grown rate
    
    sigma.proc ~ dunif(0,1)                     # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc,2)
    tau.proc <- pow(sigma.proc,-2)
    
    sigma.obs ~ dunif(0,1)                      # Prior for sd of obs.process
    sigma2.obs <- pow(sigma.obs,2)
    t.obs <- pow(sigma.obs,-2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
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
    N.cor[t] <- (N.est[t]-N.mean)/N.sd         # standardize count to be used as covariate
    }
    }
    ",fill = TRUE)
sink()

###################################################################
# Load data
#-------------------------------------------------------------------

M <- read.table("merlins.txt",header=TRUE)

#First, alter data input
hawk <- c(M$HM)
white <- c(M$WP)
mat <- matrix(c(hawk, white), nrow=length(hawk))
#---------------------------------------------------------------------

# Bundle data
N.mean = 114.9
N.sd = 17.3
jags.data <- list(nyears = nyears, marray.j = marray.j, marray.a = marray.a, y = y, J = J, R = R, r.j = rowSums(marray.j), r.a = rowSums(marray.a), x=log(mat),T=nrow(mat), S=ncol(mat), N.mean = N.mean, N.sd = N.sd)

# Initial values
inits <- function(){list(l.mphij = rnorm(1, 0.2, 0.5), l.mphia = rnorm(1, 0.2, 0.5), l.mfec = rnorm(1, 0.2, 0.5), l.p = rnorm(1, 0.2, 1), sig.phij = runif(1, 0.1, 10), sig.phia = runif(1, 0.1, 10), sig.fec = runif(1, 0.1, 10),  n1 = round(runif(1, 1, 50), 0), nadSurv = round(runif(1, 5, 50), 0), beta.phia = runif(1, -1, 1), b0.omm = runif(1, 0, 10), sig.im = runif(1, 0.1, 10), nadimm = round(runif(1, 1, 50), 0), sigma.proc = runif(1, 0, 1), mean.r = rnorm(1), sigma.obs = runif(1, 0, 1),logN.est = c(rnorm(1, 4.4, 0.1), rep(NA, (nrow(mat) - 1))))}

# Parameters monitored
parameters <- c("phij", "phia", "f",  "p", "lambda", "mphij", "mphia", "mfec", "mlam", "beta.phia","sig.phij", "sig.phia", "sig.fec", "sig.obs", "N1", "NadSurv", "Ntot", "omega", "sig.im", "Nadimm", "b0.omm", "imrate", "r", "mean.r", "sigma2.obs", "sigma2.proc", "N.cor", "N.est")

# MCMC settings
ni <- 400000    
nt <- 10
nb <- 200000
nc <- 3

# Call JAGS from R
pipl.ipm.merlin <- jags(jags.data, inits, parameters, "pipl.ipm.merlin.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, parallel = TRUE, store.data = TRUE)

#----------------------------------------------------------------------------------------------------------------------------------------
# Code for fig of counts vs. ests, annual adult and juv surv. prob, fecundity, immigration [given merlin effect on adult survival]
par(mfrow = c(2, 2), cex.axis = 1, cex.lab = 1, las = 1, mar = c(5, 5, 1, 1), mgp=c(3, 1, 0))
lower <- upper <- numeric()
year <- 1993:2016
for (i in 1:nyears){
  lower[i] <- quantile(pipl.ipm.merlin$sims.list$Ntot[,i], 0.025)
  upper[i] <- quantile(pipl.ipm.merlin$sims.list$Ntot[,i], 0.975)}
m1 <- min(c(pipl.ipm.merlin$mean$Ntot, y, lower), na.rm = T)
m2 <- max(c(pipl.ipm.merlin$mean$Ntot, y, upper), na.rm = T)
plot(0, 0, ylim = c(0, m2), xlim = c(1, nyears), ylab = "Population size (pairs)", xlab = " ", col = "black", type = "l",  axes = F, frame = F)
axis(2)
axis(1, at = 1:nyears, labels = year)
polygon(x = c(1:nyears, nyears:1), y = c(lower, upper[nyears:1]), col = "grey85", border = "grey85")
points(y, type = "l", col = "grey30", lwd = 2)
points(pipl.ipm.merlin$mean$Ntot, type = "l", col = "cornflowerblue", lwd = 2)
legend(x = 0, y = 10, legend = c("Counts", "Estimates"), lty = c(1, 1),lwd = c(2, 2), col = c("grey30", "cornflowerblue"), bty = "n", cex = 1)

lower <- upper <- numeric()
T <- nyears-1
for (t in 1:T){
  lower[t] <- quantile(pipl.ipm.merlin$sims.list$phij[,t], 0.025)
  upper[t] <- quantile(pipl.ipm.merlin$sims.list$phij[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = pipl.ipm.merlin$mean$phij, x = (1:T)+0.5, xlim= c(1, 24), type = "b", pch = 16, ylim = c(0, 1.0), ylab = "Annual survival probability", xlab = "", axes = F, cex = 1.1, frame = F, lwd = 1.3)
axis(2)
axis(1, at = 1:(T+1), labels = 1993:2016)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)
segments(1, pipl.ipm.merlin$mean$mphij, T+1, pipl.ipm.merlin$mean$mphij, lty = 1, lwd = 1.4, col = "violetred")
segments(1, quantile(pipl.ipm.merlin$sims.list$mphij, 0.025), T+1, quantile(pipl.ipm.merlin$sims.list$mphij, 0.025), lty = 2, col = "violetred") 
segments(1, quantile(pipl.ipm.merlin$sims.list$mphij, 0.975), T+1, quantile(pipl.ipm.merlin$sims.list$mphij, 0.975), lty = 2, col = "violetred")
for (t in 1:T){
  lower[t] <- quantile(pipl.ipm.merlin$sims.list$phia[,t], 0.025)
  upper[t] <- quantile(pipl.ipm.merlin$sims.list$phia[,t], 0.975)}
points(y=pipl.ipm.merlin$mean$phia, x = (1:T)+0.5, type = "b", pch = 1, cex = 1.1, lwd = 1.3)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)
segments(1, pipl.ipm.merlin$mean$mphia, T+1, pipl.ipm.merlin$mean$mphia, lty = 1, lwd = 1.4, col = "violetred")
segments(1, quantile(pipl.ipm.merlin$sims.list$mphia, 0.025), T+1, quantile(pipl.ipm.merlin$sims.list$mphia, 0.025), lty = 2, col = "violetred") 
segments(1, quantile(pipl.ipm.merlin$sims.list$mphia, 0.975), T+1, quantile(pipl.ipm.merlin$sims.list$mphia, 0.975), lty = 2, col = "violetred")
legend(x = 0, y = 0.15, legend = c("Adults", "Juveniles"), pch = c(1, 16), bty = "n")

lower <- upper <- numeric()
T <- nyears-1
for (t in 1:T){
  lower[t] <- quantile(pipl.ipm.merlin$sims.list$f[,t], 0.025)
  upper[t] <- quantile(pipl.ipm.merlin$sims.list$f[,t], 0.975)}
plot(y=pipl.ipm.merlin$mean$f, x = (1:T), type = "b", pch = 16, ylim = c(0, 4), xlim=c(1,24), ylab = "Fecundity (fledgling / female)", xlab = "", axes = F, cex = 1.1, frame = F, lwd = 1.3)
axis(2)
axis(1, at = 1:(T+1), labels = 1993:2016)
segments((1:T), lower, (1:T), upper)
segments(1, pipl.ipm.merlin$mean$mfec, T, pipl.ipm.merlin$mean$mfec, lty = 1, lwd = 1.4, col = "violetred")
segments(1, quantile(pipl.ipm.merlin$sims.list$mfec, 0.025), T, quantile(pipl.ipm.merlin$sims.list$mfec, 0.025), lty = 2, col = "violetred") 
segments(1, quantile(pipl.ipm.merlin$sims.list$mfec, 0.975), T, quantile(pipl.ipm.merlin$sims.list$mfec, 0.975), lty = 2, col = "violetred")

lower <- upper <- numeric()
T <- nyears-1
for (t in 1:T){
  lower[t] <- quantile(pipl.ipm.merlin$sims.list$omega[,t], 0.025)
  upper[t] <- quantile(pipl.ipm.merlin$sims.list$omega[,t], 0.975)}
plot(y = pipl.ipm.merlin$mean$omega, x = (1:T)+0.5, xlim = c(1, 24), type = "b", pch = 16, ylim = c(0, 12), ylab = "Immigration (no. indivs)", xlab = "", axes = F, cex = 1.1, frame = F, lwd = 1.3)
axis(2)
axis(1, at = 1:(T+1), labels = 1993:2016)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)
segments(1, pipl.ipm.merlin$mean$b0.omm, T+1, pipl.ipm.merlin$mean$b0.omm, lty = 1, lwd = 1.4, col = "violetred")
segments(1, quantile(pipl.ipm.merlin$sims.list$b0.omm, 0.025), T+1, quantile(pipl.ipm.merlin$sims.list$b0.omm, 0.025), lty = 2, col = "violetred") 
segments(1, quantile(pipl.ipm.merlin$sims.list$b0.omm, 0.975), T+1, quantile(pipl.ipm.merlin$sims.list$b0.omm, 0.975), lty = 2, col = "violetred")


# Code for demo rates vs. pop growth and some descriptive statistics including MERLIN effect on adult survival
nyears <- 24
lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = nyears-1, ncol = 5)

for (i in 1:(nyears-1)){
  lambda.h[i] <- mean(pipl.ipm.merlin$sims.list$lambda[,i])
  lam.lower.h[i] <- quantile(pipl.ipm.merlin$sims.list$lambda[,i], 0.025)
  lam.upper.h[i] <- quantile(pipl.ipm.merlin$sims.list$lambda[,i], 0.975)
}

for (i in 1:(nyears-1)){
  Fitted.h[i,1] <- mean(pipl.ipm.merlin$sims.list$phij[,i])
  lower.h[i,1] <- quantile(pipl.ipm.merlin$sims.list$phij[,i], 0.025)
  upper.h[i,1] <- quantile(pipl.ipm.merlin$sims.list$phij[,i], 0.975)
}

for (i in 1:(nyears-1)){
  Fitted.h[i,2] <- mean(pipl.ipm.merlin$sims.list$phia[,i])
  lower.h[i,2] <- quantile(pipl.ipm.merlin$sims.list$phia[,i], 0.025)
  upper.h[i,2] <- quantile(pipl.ipm.merlin$sims.list$phia[,i], 0.975)
}

for (i in 1:(nyears-1)){
  Fitted.h[i,3] <- mean(pipl.ipm.merlin$sims.list$f[,i])
  lower.h[i,3] <- quantile(pipl.ipm.merlin$sims.list$f[,i], 0.025)
  upper.h[i,3] <- quantile(pipl.ipm.merlin$sims.list$f[,i], 0.975)
}

for (i in 1:(nyears-1)){
  Fitted.h[i,4] <- mean(pipl.ipm.merlin$sims.list$omega[,i])
  lower.h[i,4] <- quantile(pipl.ipm.merlin$sims.list$omega[,i], 0.025)
  upper.h[i,4] <- quantile(pipl.ipm.merlin$sims.list$omega[,i], 0.975)
}

####how correlated is merlin abundance with pop growth of plovers?

for (i in 1:(nyears-1)){
  Fitted.h[i,5] <- mean(pipl.ipm.merlin$sims.list$N.est[,i])
  lower.h[i,5] <- quantile(pipl.ipm.merlin$sims.list$N.est[,i], 0.025)
  upper.h[i,5] <- quantile(pipl.ipm.merlin$sims.list$N.est[,i], 0.975)
}


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 5, nrow = 60000)
for (i in 1:60000){
  correl.h[i,1] <- cor(pipl.ipm.merlin$sims.list$lambda[i,], pipl.ipm.merlin$sims.list$phij[i,])
  correl.h[i,2] <- cor(pipl.ipm.merlin$sims.list$lambda[i,], pipl.ipm.merlin$sims.list$phia[i,])
  correl.h[i,3] <- cor(pipl.ipm.merlin$sims.list$lambda[i,], pipl.ipm.merlin$sims.list$f[i,])
  correl.h[i,4] <- cor(pipl.ipm.merlin$sims.list$lambda[i,], pipl.ipm.merlin$sims.list$omega[i,])
  correl.h[i,5] <- cor(pipl.ipm.merlin$sims.list$lambda[i,], pipl.ipm.merlin$sims.list$N.est[i,1:23])
}

# Credible intervals of correlation coefficients
quantile(correl.h[,1], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,2], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,3], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,4], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,5], c(0.05, 0.5, 0.95), na.rm = TRUE)

# Compute the posterior modes of correlation coefficients
m <- density(correl.h[,1], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,2], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,3], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,4], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,5], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

# Probability that correlation coefficients (r) > 0
sum(correl.h[!is.na(correl.h[,1]),1]>0)/60000
sum(correl.h[!is.na(correl.h[,2]),2]>0)/60000
sum(correl.h[!is.na(correl.h[,3]),3]>0)/60000
sum(correl.h[!is.na(correl.h[,4]),4]>0)/60000
sum(correl.h[!is.na(correl.h[,5]),5]<0)/60000

# Plot retrospective fig
par(mfrow = c(3, 2), mar = c(5, 4, 1.5, 1), mgp=c(3, 1, 0), las = 1, cex = 0.9)
linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0.1, 0.5), ylim = c(0.6, 1.8), ylab = "Population growth rate", xlab = "Juvenile survival", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "darkblue")
text(x = 0.1, y = 1.75, "r = 0.42 (0.09, 0.63)", pos = 4, font = 3, cex = 0.8)
text(x = 0.1, y = 1.65, "P(r>0) = 0.98", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,2], type = "n", xlim = c(0.5, 1.0), ylim = c(0.6, 1.8),  ylab = "", xlab = "Adult survival", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,2], lam.lower.h, Fitted.h[,2], lam.upper.h, col = linecol)
segments(lower.h[,2], lambda.h, upper.h[,2], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,2], pch = 19, col = "darkblue")
text(x = 0.5, y = 1.75, "r = 0.34 (0.14, 0.57)", pos = 4, font = 3, cex = 0.8)
text(x = 0.5, y = 1.65, "P(r>0) = 0.99", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,3], type = "n", xlim = c(0.5, 3.0), ylim = c(0.6, 1.8), ylab = "Population growth rate", xlab = "Fecundity", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,3], lam.lower.h, Fitted.h[,3], lam.upper.h, col = linecol)
segments(lower.h[,3], lambda.h, upper.h[,3], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,3], pch = 19, col = "darkblue")
text(x = 0.5, y = 1.75, "r = 0.37 (0.07, 0.59)", pos = 4,  font = 3, cex = 0.8)
text(x = 0.5, y = 1.65, "P(r>0) = 0.98", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,4], type = "n", xlim = c(0, 10), ylim = c(0.6, 1.8),  ylab = "", xlab = "Immigration", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,4], lam.lower.h, Fitted.h[,4], lam.upper.h, col = linecol)
segments(lower.h[,4], lambda.h, upper.h[,4], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,4], pch = 19, col = "darkblue")
text(x = 0.0, y = 1.75, "r = 0.31 (-0.15, 0.77)", pos = 4, font = 3, cex = 0.8)   #0.61
text(x = 0.0, y = 1.65, "P(r>0) = 0.87", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,5], type = "n", xlim = c(50, 160), ylim = c(0.6, 1.8),  ylab = "Population growth rate", xlab = "Merlin abundance", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,5], lam.lower.h, Fitted.h[,5], lam.upper.h, col = linecol)
segments(lower.h[,5], lambda.h, upper.h[,5], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,5], pch = 19, col = "darkblue")
text(x = 50.0, y = 1.75, "r = -0.25 (-0.43, -0.08)", pos = 4, font = 3, cex = 0.8)
text(x = 50.0, y = 1.65, "P(r<0) = 0.99", pos = 4, font = 3, cex = 0.8)
