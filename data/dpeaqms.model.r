dpeaqms.model<- "model {     
   # The prior distribution on Tau, the precision of the log-intensity measurements, is
   # a gamma distribution with parameters a.tau and b.tau.
   # By default these are set so that the distribution on Tau is flat
   Tau  ~ dgamma(a.tau,b.tau)   
   Sigma <- 1/sqrt(Tau) 
     
   # Alpha is the log-intensity for the control group for a given peptide instance k
   # Its prior is normally distributed. Ideally this should be set so that it is flat
   # over the dynamic range of your mass spectrometer
   for (k in 1:Npeptides) {
       Alpha[k] ~ dnorm(a.alpha,b.alpha)
   }
   
  
   for (e in 1:E) {
     kappa[e,1] <- 0
     for (s in 2:ExpSamples[e]) {
        kappa[e,s] ~ dnorm(a.kappa,b.kappa)
     }
   }

   # Beta[g,j] is an binary parameter which is 0 if there is no differential expression 
   # between the group i and the control group (1)for the protein.
   # Gamma[i,protein] is the corresponding difference in the mean log-intensities 
   # between the group i and the control group
   # The prior distribution on Beta[i,x] is a Bernoulli distribution with parameter p[i]
   # The prior distribution on Gamma[i,x] is a Gaussian distribution with parameters c and d
   # By defaut c and d are set to 0 and 1 i.e. a standard normal distribution
   for (j in 1:P) {   
     p[1,j] <- 0.0  
     Beta[1,j]  <- 0.0
     Gamma[1,j] <- 0.0
     for (g in 2:G) {
       p[g,j] ~ dbeta(a.p,b.p)      
       Beta[g,j] ~ dbern(p[g,j]) 
       Gamma[g,j] ~ dnorm(a.gamma,b.gamma) ; 
     }  

   # The model for the data likelihood (see documentation)
   for (n in 1:N[j]) {                 
        intensity[offset[j]+n] ~ dnorm(kappa[experiment[offset[j]+n] , sample[offset[j]+n]-sampleoffset[experiment[offset[j]+n]]] +
                                             Alpha[peptide[offset[j]+n]] +
                                             Beta[group[offset[j]+n],j]*Gamma[group[offset[j]+n],j],
                                             Tau);
     }      
     
   }    
  }
  "
