###################################################################%
#### Script of MSOMs used in the paper published in Ecol. Ent. ####%
###################################################################%

####
# Last modification 05/01/2021
# By B. Mourguiart
####


library(jagsUI)
load(here::here("Data/Data_analyse/data_model_article.RData"))


##### Model used in the paper ####


cat('
    model{
    
    #Define prior distributions for community-level model parameters
    omega ~ dunif(0,1)                            #probability of species inclusion
    
    a0.mean ~ dunif(0,1)                          #species level effect on occupancy
    mu.a0 <- log(a0.mean) - log(1-a0.mean)        #species level effect on occupancy on the logit scale
    mu.a1 ~ dnorm(0, 0.001)                        #site covariable 1 effect on occupancy
    mu.a2 ~ dnorm(0, 0.001)
    
    tau.a0 ~ dgamma(0.1,0.1)                      #variability of species level effect on occupancy
    tau.a1 ~ dgamma(0.1,0.1)                      #variability of site covariable 1 effect on occupancy
    tau.a2 ~ dgamma(0.1,0.1) 
    
    mu.b1 ~ dnorm(0, 0.001)
    
    mu.b0 ~ dnorm(0, 0.001)
    mu.b2 ~ dnorm(0, 0.001)
    mu.b3 ~ dnorm(0, 0.001)
    
    tau.b0 ~ dgamma(0.1,0.1)
    tau.b2 ~ dgamma(0.1,0.1)
    tau.b3 ~ dgamma(0.1,0.1)
    
    tau.b1 ~ dgamma(0.1,0.1)


    for (i in 1:(n+nzeroes)) {
    
    #Create priors for species i from the community level prior distributions
    w[i] ~ dbern(omega)                       #binary indicator: species i belongs or not to the community
    
    a0[i] ~ dnorm(mu.a0, tau.a0)                #species effect on occupancy of the species i  
    a1[i] ~ dnorm(mu.a1, tau.a1)               #site covariable 1 effect on occupancy of the species i
    a2[i] ~ dnorm(mu.a2, tau.a2)              
    
    b0[i] ~ dnorm(mu.b0, tau.b0)          #parameters for detection probabilities depend on mobility group
    b1[i] ~ dnorm(mu.b1, tau.b1)
    b2[i] ~ dnorm(mu.b2, tau.b2)
    b3[i] ~ dnorm(mu.b3, tau.b3)
    
    
    #Create a loop to estimate the Z matrix (true occurrence for species i 
    #at point j.      
    for (j in 1:J) {
    logit(psi[j,i]) <- a0[i] + 
		       a1[i]*covSite1[j] + a2[i]*pow(covSite1[j],2)  #altitude effect
		       
    mu.psi[j,i] <- psi[j,i]*w[i]          #site j could be occupied by species i only if it belongs to the community
    #if species i belongs to the community: mu.psi=psi, else: mu.psi=0
    Z[j,i] ~ dbern(mu.psi[j,i])           #binary indicator: species i presents or not in site j
    
    #Create a loop to estimate detection for species i at point j during 
    #sampling period k.      
    for (k in 1:K[j]) {  
    logit(p[j,k,i]) <-  b0[i] + #view
                        b1[i]*covDetection1[j,k] + #hearing
                        b2[i]*covDetection2[j,k] + #sweep netting
                        b3[i]*covDetection0[j,k]*covDetection3[j,k]  #sighting*grass height
                        
    mu.p[j,k,i] <- p[j,k,i]*Z[j,i]    #species i could be detected in site j during survey k only if site j is occupied
    #if j is occupied: mu.p=p, else: mu.p=0
    X[j,k,i] ~ dbern(mu.p[j,k,i])     #binary indicator observed: species i detected or not in site j during survey k
    
    #Create simulated dataset to calculate the Bayesian p-value
    Xnew[j,k,i] ~ dbern(mu.p[j,k,i])
    
    #Pearson residuals
    d[j,k,i]<-  abs(X[j,k,i] - mu.p[j,k,i])
    dnew[j,k,i]<- abs(Xnew[j,k,i]- mu.p[j,k,i])
    d2[j,k,i]<- pow(d[j,k,i],2)
    dnew2[j,k,i]<- pow(dnew[j,k,i],2)
    
    }   
    
    dsum[j,i]<- sum(d2[j,1:K[j],i]) 
    dnewsum[j,i]<- sum(dnew2[j,1:K[j],i])
    
    }
    }
    
    
    #Calculate the discrepancy measure
    p.fit<-sum(dsum[1:J,1:(n)])
    p.fitnew<-sum(dnewsum[1:J,1:(n)])
    
    #Sum all species observed (n) and unobserved species (n0) to find the 
    #total estimated richness
    n0 <- sum(w[(n+1):(n+nzeroes)])
    N <- n + n0
    
    #Create a loop to determine point level richness estimates for the 
    #whole community.
    for(j in 1:J){
      Nsite[j]<- sum(Z[j,1:(n+nzeroes)])
    }
}', file = here::here("Analyse/model_covarHerbe_articlev2.txt"))


#Load all the data needed to run the model

sp.data2 = list(n=n, nzeroes=50, J=J, K=rep(K,81), X=X,
                covDetection0=methodeV, covDetection1=methodeE, covDetection2=methodeF, 
                covDetection3=cov.detection_hauteur,
                covSite1=cov.sites_standard$altitude
)


#Specify the parameters to be monitored
sp.params2 = c('mu.a0','mu.a1','mu.a2', 
               'mu.b0','mu.b1', 'mu.b2', 'mu.b3',
               'tau.a0','tau.a1', 'tau.a2', 
               'tau.b0','tau.b1', 'tau.b2', 'tau.b3',
               'omega','a0', 'a1','a2', 
               'b0', 'b1', 'b2', 'b3',
               'Nsite', 'N', 'p.fit','p.fitnew') 

#Specify the initial values
sp.inits2 = function() {
  omegaGuess = runif(1, 0, 1)
  psi.meanGuess = runif(1, 0,1)
  list(w=c(rep(1, n), rbinom(50, size=1, prob=1)),
       a0.mean=runif(1,0,1),
       Z = array(1,dim=c(J, n+50)), 
       mu.a1=rnorm(1), mu.a2=rnorm(1), 
       mu.b0=rnorm(1),mu.b1=rnorm(1), mu.b2=rnorm(1), mu.b3=rnorm(1))
}



#Run the model

out_articlev2 <- autojags(data = sp.data2, 
                                inits = sp.inits2, 
                                parameters.to.save = sp.params2, 
                                model.file=here::here("Analyse/model_covarHerbe_articlev2.txt"),  
                                n.chains=3, 
                                n.adapt = 500,
                                n.burnin=15000,  
                                n.thin=15,  
                                iter.increment = 15000,
                                max.iter = 60000,
                                DIC=TRUE,
                                store.data = T,
                                Rhat.limit = 1.1,
                                parallel = T)
#save(out_articlev2, file="out_articlev2.RData")


#### Model for ESM 1 ####

cat('
    model{
    
    #Define prior distributions for community-level model parameters
    omega ~ dunif(0,1)                            #probability of species inclusion
    
    a0.mean ~ dunif(0,1)                          #species level effect on occupancy
    mu.a0 <- log(a0.mean) - log(1-a0.mean)        #species level effect on occupancy on the logit scale
    mu.a1 ~ dnorm(0, 0.001)                        #site covariable 1 effect on occupancy
    mu.a2 ~ dnorm(0, 0.001)
    mu.a3 ~ dnorm(0, 0.001)    

    tau.a0 ~ dgamma(0.1,0.1)                      #variability of species level effect on occupancy
    tau.a1 ~ dgamma(0.1,0.1)                      #variability of site covariable 1 effect on occupancy
    tau.a2 ~ dgamma(0.1,0.1) 
    tau.a3 ~ dgamma(0.1,0.1)
    
    mu.b1 ~ dnorm(0, 0.001)
    
    mu.b0 ~ dnorm(0, 0.001)
    mu.b2 ~ dnorm(0, 0.001)
    mu.b3 ~ dnorm(0, 0.001)
    
    tau.b0 ~ dgamma(0.1,0.1)
    tau.b2 ~ dgamma(0.1,0.1)
    tau.b3 ~ dgamma(0.1,0.1)
    
    tau.b1 ~ dgamma(0.1,0.1)


    for (i in 1:(n+nzeroes)) {
    
    #Create priors for species i from the community level prior distributions
    w[i] ~ dbern(omega)                       #binary indicator: species i belongs or not to the community
    
    a0[i] ~ dnorm(mu.a0, tau.a0)                #species effect on occupancy of the species i  
    a1[i] ~ dnorm(mu.a1, tau.a1)               #site covariable 1 effect on occupancy of the species i
    a2[i] ~ dnorm(mu.a2, tau.a2)              
    a3[i] ~ dnorm(mu.a3, tau.a3)
    
    b0[i] ~ dnorm(mu.b0, tau.b0)          #parameters for detection probabilities depend on mobility group
    b1[i] ~ dnorm(mu.b1, tau.b1)
    b2[i] ~ dnorm(mu.b2, tau.b2)
    b3[i] ~ dnorm(mu.b3, tau.b3)
    
    
    #Create a loop to estimate the Z matrix (true occurrence for species i 
    #at point j.      
    for (j in 1:J) {
    logit(psi[j,i]) <- a0[i] + 
		       a1[i]*covSite1[j] + a2[i]*pow(covSite1[j],2) + #altitude effect
		       a3[i]*covSite2[j]
		       
    mu.psi[j,i] <- psi[j,i]*w[i]          #site j could be occupied by species i only if it belongs to the community
    #if species i belongs to the community: mu.psi=psi, else: mu.psi=0
    Z[j,i] ~ dbern(mu.psi[j,i])           #binary indicator: species i presents or not in site j
    
    #Create a loop to estimate detection for species i at point j during 
    #sampling period k.      
    for (k in 1:K[j]) {  
    logit(p[j,k,i]) <-  b0[i] + #*covDetection0[j,k] + #view
                        b1[i]*covDetection1[j,k] + #hearing
                        b2[i]*covDetection2[j,k] + #sweep netting
                        b3[i]*covDetection0[j,k]*covDetection3[j,k] #sighting*grass height
                        
    mu.p[j,k,i] <- p[j,k,i]*Z[j,i]    #species i could be detected in site j during survey k only if site j is occupied
    #if j is occupied: mu.p=p, else: mu.p=0
    X[j,k,i] ~ dbern(mu.p[j,k,i])     #binary indicator observed: species i detected or not in site j during survey k
    
    #Create simulated dataset to calculate the Bayesian p-value
    Xnew[j,k,i] ~ dbern(mu.p[j,k,i])
    
    #Pearson residuals
    d[j,k,i]<-  abs(X[j,k,i] - mu.p[j,k,i])
    dnew[j,k,i]<- abs(Xnew[j,k,i]- mu.p[j,k,i])
    d2[j,k,i]<- pow(d[j,k,i],2)
    dnew2[j,k,i]<- pow(dnew[j,k,i],2)
    
    }   
    
    dsum[j,i]<- sum(d2[j,1:K[j],i]) 
    dnewsum[j,i]<- sum(dnew2[j,1:K[j],i])
    
    }
    }
    
    
    #Calculate the discrepancy measure
    p.fit<-sum(dsum[1:J,1:(n)])
    p.fitnew<-sum(dnewsum[1:J,1:(n)])
    
    #Sum all species observed (n) and unobserved species (n0) to find the 
    #total estimated richness
    n0 <- sum(w[(n+1):(n+nzeroes)])
    N <- n + n0
    
    #Create a loop to determine point level richness estimates for the 
    #whole community.
    for(j in 1:J){
      Nsite[j]<- sum(Z[j,1:(n+nzeroes)])
    }
}', file=here::here("Analyse/model_covarHerbe_bothv2.txt"))


#Load all the data
sp.data3 = list(n=n, nzeroes=50, J=J, K=rep(K,81), X=X,
                covDetection0=methodeV, covDetection1=methodeE, covDetection2=methodeF, 
                covDetection3=cov.detection_hauteur,
                covSite1=cov.sites_standard$altitude,
                covSite2=cov.sites_standard$hauteur
)

#Specify the parameters to be monitored
sp.params3 = c('mu.a0','mu.a1','mu.a2', 'mu.a3',
               'mu.b0','mu.b1', 'mu.b2', 'mu.b3',
               'tau.a0','tau.a1', 'tau.a2', 'tau.a3',
               'tau.b0','tau.b1', 'tau.b2', 'tau.b3',
               'omega','a0', 'a1','a2', 'a3',
               'b0', 'b1', 'b2', 'b3',
               'Nsite', 'N', 'p.fit','p.fitnew') 

#Specify the initial values
sp.inits3 = function() {
  omegaGuess = runif(1, 0, 1)
  psi.meanGuess = runif(1, 0,1)
  list(w=c(rep(1, n), rbinom(50, size=1, prob=1)),
       a0.mean=runif(1,0,1),
       Z = array(1,dim=c(J, n+50)), 
       mu.a1=rnorm(1), mu.a2=rnorm(1), mu.a3=rnorm(1),
       mu.b0=rnorm(1),mu.b1=rnorm(1), mu.b2=rnorm(1), mu.b3=rnorm(1))
}



#run

out_covarHerbe_bothv2 <- autojags(data = sp.data3, 
                                inits = sp.inits3, 
                                parameters.to.save = sp.params3, 
                                model.file=here::here("Analyse/model_covarHerbe_bothv2.txt"),  
                                n.chains=nc, 
                                n.adapt = 1000,
                                n.burnin=10000,  
                                n.thin=nt,  
                                iter.increment = 5000,
                                max.iter = 60000,
                                DIC=TRUE,
                                store.data = T,
                                Rhat.limit = 1.1,
                                parallel = T)
save(out_covarHerbe_bothv2, file=here::here("Data", "Data_analyse", "out_covarHerbe_bothv2.RData"))
