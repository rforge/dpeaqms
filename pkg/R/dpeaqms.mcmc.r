# R function that runs an MCMC sampler to sample from the
# posterior density of the model parameters of a differential
# protein expression model based on intensity data from
# isobaric labeled MSMS data
dpeaqms.mcmc<-function(msmsdata, burnin=100000, samples=1000, thin=100,
                       a.alpha=10.0, b.alpha=1.0/9.0,
                       a.gamma=0.0,  b.gamma=1.0,
                       a.tau=1e-3, b.tau=1e-3,
                       a.kappa=0.0,b.kappa=1.0/9.0, 
                       a.p=1.0, b.p=19.0,
                       controlGroup=NULL, referenceSampleID=NULL, 
	               transform=T, infofilename=NULL, verbose=F) {
  
  # Load the R JAGS library
  require(rjags)
  
  # End model specification

  # Record run parameters to the specified info file
  if (!is.null(infofilename)) {
     infofile = file(infofilename, "w")
     line = paste('a.alpha = ' , a.alpha, '\n')
     cat(line,file=infofile)
     line = paste('b.alpha = ' , b.alpha, '\n')
     cat(line,file=infofile)
     line = paste('a.gamma = ' , a.gamma, '\n')
     cat(line,file=infofile)
     line = paste('b.gamma = ' , b.gamma, '\n')
     cat(line,file=infofile)
     line = paste('a.tau = ' , a.tau, '\n')
     cat(line,file=infofile)
     line = paste('b.tau = ' , b.tau, '\n')
     cat(line,file=infofile)
     line = paste('a.kappa = ' , a.kappa, '\n')
     cat(line,file=infofile)
     line = paste('b.kappa = ' , b.kappa, '\n')
     cat(line,file=infofile)
     line = paste('a.p = ' , a.p, '\n')
     cat(line,file=infofile)
     line = paste('b.p = ' , b.p, '\n')
     cat(line,file=infofile)
     line = paste('burnin runs = ' , burnin, '\n')
     cat(line,file=infofile)
     line = paste('sample runs = ' , samples*thin, '\n')
     cat(line,file=infofile)
     line = paste('thin factor = ' , thin, '\n')
     cat(line,file=infofile)
     close(infofile)
  } 

  proteinID    = msmsdata$proteinID
  intensity    = msmsdata$intensity
  peptideID    = msmsdata$peptideID
  sampleID     = msmsdata$sampleID
  groupID      = msmsdata$groupID
  experimentID = msmsdata$experimentID

  E = length(unique(experimentID))
 
  # Log transform the intensity data if transform=T
  if (transform) {
     intensity = log(intensity)
  }
 
  # Order the data by protein, peptide, replicate and finally sample 
  myord = order(proteinID, peptideID, sampleID)
  proteinID = proteinID[myord]
  peptideID = peptideID[myord]
  intensity = intensity[myord]
  sampleID = sampleID[myord]
  groupID = groupID[myord]
  experimentID  = experimentID[myord]

  # Factor experiment identifiers
  experimentID = factor(experimentID, levels=sort(unique.default(experimentID)))
  # Make a note of the experiment identifiers
  experimentLevels = levels(experimentID)
  # Numerically re-encode the experiment identifier vector
  levels(experimentID) <-seq(1,E)
  experimentID <- as.numeric(as.vector(experimentID))
  numberOfSamples = rep(0,E)
  for (i in seq(1,E)) {
    numberOfSamples[i] = length(unique(sampleID[experimentID==i]))
  }
  
  # Record the experiment offset i.e. where in the intensity vector
  # observations from the experiment start
  sampleoffset = rep(0,E)
  if (E > 1) {
   for (i in seq(2, E)) {
     sampleoffset[i] = sampleoffset[i-1]+numberOfSamples[i-1]
    }
  }
 
  
  # Factor ProteinID
  proteinID <-factor(proteinID, levels=unique.default(proteinID))
  # Store the orignal protein identifications
  pid = proteinID
  # Get the proteins (in the order in which they appear in the data file)  
  proteins  <- levels(proteinID)
  P <- length(proteins)
  # Numerically re-encode the proteins
  levels(proteinID)<-seq(1,length(levels(proteinID)))
  
  # Number of measurements for each protein and their "offset" in the list
  N <- vector(mode="integer",length=P) 
  m <- vector(mode="integer",length=P)
  offset <- vector(mode="integer" , length=P)
  moffset<- vector(mode="integer", length=P)
  # Factor SampleID
  sampleID    = factor(sampleID, levels=sort(unique.default(sampleID))) 
  sampleLevels = levels(sampleID)
  Nsamples = length(sampleLevels)
  
  # Numerically re-encode the sample identifiers
  samplenumericlevels <- seq(1,Nsamples)

  # If a reference sample (used for the sample normalization) is specified
  if (!is.null(referenceSampleID)) {
    soffset = 0 
    for (i in seq(1,E)) {
      # Use the specified sample identifier in each experiment as the reference     
      refSample = paste("Exp", i ,  "." , referenceSampleID, sep ="")
      # Make sure the reference samples specified is a valid samples
      x = match(refSample, sampleLevels)      
      if (is.na(x)) {
        print(paste("Error: Unable to find specified reference sample \"" , refSample , "\" in Experiment " , datafile[i], sep=""))
        return(NULL)
      }
      # Reorder the sample levels so that the specified sample identifers are the first
      # in each experiment
      temp = sampleLevels[soffset+1]
      sampleLevels[soffset+1] = sampleLevels[x]
      sampleLevels[x] = temp
      temp = samplenumericlevels[soffset+1]
      samplenumericlevels[soffset+1] = samplenumericlevels[x]
      samplenumericlevels[x] = temp   
      soffset = soffset + numberOfSamples[i]
      
    }      
  }

  if (verbose) {
    print(sampleLevels)
  }
  levels(sampleID) <-samplenumericlevels
  sampleID <- as.numeric(as.vector(sampleID))

  # Determine the peptide/replicate blocks 
  peptideID = factor(peptideID, levels=unique.default(peptideID))
  # Remove reserved operator symbols from peptide sequences and modifiers
  peptideNames = gsub('[(|)|\"]', '' , unique.default(peptideID))
 
  # Numerically encode the peptides
  Npeptides <- length(levels(peptideID))
  levels(peptideID) <- seq(1,length(levels(peptideID)))
  peptideID <- as.numeric(as.vector(peptideID))
  
  if (verbose) {
     print(paste("Processing " , Npeptides , " MSMS spectra" , sep=""))
  }

  # Numerically encode the group identifiers
  groupID = factor(groupID,levels=sort(unique(groupID)))
  groupLevels = levels(groupID)
  G = length(groupLevels)
  numericlevels <- seq(1,G)
  
  # If a control group is specified make sure it appears first in the list of groups
  if (!is.null(controlGroup)) {
    x = match(controlGroup, groupLevels)
    if (is.na(x)) {
      # If the control group doesn't exist print an error message and return a NULL
      print(paste("Error: Unable to find specified control group \"" , controlGroup , "\"" , sep=""))
      return(NULL)
    }
    temp = groupLevels[1]
    groupLevels[1] = groupLevels[x]
    groupLevels[x] = temp ;   
    temp = numericlevels[1]
    numericlevels[1] = numericlevels[x]
    numericlevels[x] = temp    
  }
  
  if (verbose) {
     print(groupLevels)
  }
   
  levels(groupID) <- numericlevels
  groupID <- as.numeric(as.vector(groupID))
  # Set the offset (i.e. what data point the protein observations start at) and number of observations per protein
  label = (pid==proteins[1])
  N[1] <- sum(label) 
  m[1] <- length(unique(peptideID[label]))
  offset[1] = 0
  moffset[1] = 0
  if (P > 1) {
   for (i in seq(2,P)) {
     prot = proteins[i] 
     label = (pid==prot)
     N[i] <- sum(label)      
     m[i] <- length(unique(peptideID[label]))
     offset[i] <- offset[i-1] + N[i-1]   
     moffset[i] <- moffset[i-1] + m[i-1]
   }
  }
# print(moffset)
  # Set data and prior parameter values for the model
  modeldata = list("intensity"=intensity, "experiment"=experimentID, "offset"=offset, 
                   "numberOfSamples"=numberOfSamples, "sampleoffset"=sampleoffset, 
                   "group"=groupID, "peptide"=peptideID, "sample"=sampleID,
                   "N"=N,"E"=E,"G"=G,"moffset"=moffset,"m"=m,"P"=P,
                   "a.alpha"=a.alpha, "b.alpha"=b.alpha,
                   "a.gamma"=a.gamma,"b.gamma"=b.gamma,
                   "a.kappa"=a.kappa,"b.kappa"=b.kappa,
                   "a.tau"=a.tau, "b.tau"=b.tau,
                   "a.p"=a.p, "b.p"=b.p
                   )

  # Parameter Initialization
  
  # Initializate "p" parameter
  
 
  # Initialize Beta parameter by simulating from its prior
  BetaInit = matrix(nrow=G, ncol=P)  
  p        = matrix(nrow=G, ncol=P)
  BetaInit[1,] = rep(0.0,P)
  
  # By definition the control group is non-differentially expressed w.r.t. itself   
  p[1,] = rep(0.0,P)
  for (g in seq(2,G)) {
    p[g,] = rbeta(P,a.p,b.p)
    BetaInit[g,] = rbinom(P,1,prob=p[g,])    
  }  

   if (verbose) {
   print(G)
   print(Npeptides)
  print(P)
  }
 # print(groupID)
 # print(intensity)
  # Initialize the alpha parameters
  # Initializing from the prior can cause problems with rjags so
  # initialize to the observed log intensities for the MSMS spectra
  # provided the intensities are not for samples from groups which
  # are differentially expressed
  AlphaInit = vector(length=Npeptides,mode="numeric")  
  for (i in seq(1,Npeptides)) {
    labelmatch = (peptideID==i)  
    # print(sum(labelmatch))
    # Should be length one if peptides are unique to proteins
    peptideprotein = unique(proteinID[peptideID==i])
    
    if (length(peptideprotein) > 1) {
      print(paste("Error: Peptide " , peptideNames[i] , " assigned to proteins:" , peptideprotein, sep=""))
      return(NULL)
    }
    groupmatch = (groupID==1)
    for (g in 2:G) {      
       if (BetaInit[g,peptideprotein[1]]==0)
         groupmatch = groupmatch | (groupID==g)
    }
    #print(intensity[labelmatch & groupmatch])
    #print(sum(labelmatch & groupmatch))
    AlphaInit[i] = mean(intensity[labelmatch & groupmatch])  
  }
  
 # if (verbose) {
 #    print(AlphaInit)
 # }

  # Initialize the "Gamma"
  GammaInit = matrix(nrow=G, ncol=P)
  GammaInit[1,] = rep(0.0,P)
  for (g in seq(2,G)) {
    for (i in seq(1,P)) {
      ObsG = c() 
      if (BetaInit[g,i] == 1) {      
        proteinpeptideinstancesInd = unique(peptideID[proteinID==i])     
        for (k in proteinpeptideinstancesInd) {     
          obsG = c(ObsG,mean(intensity[(peptideID==k) & (groupID==g)]) -  AlphaInit[k])
        }
        GammaInit[g,i] = mean(obsG)
        # If there are no observed intensities for samples in group g for a particular protein
        # spectra initialize Gamma for that protein and from the prior
        if (is.na(GammaInit[g,i])) {       
#          GammaInit[g,i] = rnorm(1,mean=a.gamma , sd=sqrt(1/b.gamma))
            GammaInit[g,i] = 0.0
        }
      }
      else {
        GammaInit[g,i] = 0.0
      }
    }
  }
  
  #if (verbose) { 
  #   print(GammaInit)

  # Initialize Tau 
  sum = 0.0
  for (j in seq(1,length(intensity))) {
    Y = intensity[j]
    sum = sum + (Y - AlphaInit[peptideID[j]] - (GammaInit[groupID[j],proteinID[j]] * BetaInit[groupID[j],proteinID[j]]))^2
  }
  # tau is the inverse variance i.e. the precision
  #print(paste("sum = ", sum))
  myvar = sum/length(intensity)
  tauinit = 1.0/myvar
  if (verbose) {
    print(paste("Var Init = " , myvar))
  }

  kappainit = matrix(nrow=E, ncol=max(numberOfSamples), 0)
 # for (e in seq(1,E)) {
   # kappainit[e,] = rnorm(numberOfSamples[e], mean=a.kappa, sd=sqrt(1/b.kappa))
 # }
  kappainit = matrix(nrow=E, ncol=max(numberOfSamples), 0)
 # print(kappainit)
  modelinits = list("Tau"=tauinit , "Alpha"=AlphaInit, "Gamma"=GammaInit, "Beta"=BetaInit, "p"=p , "kappa"=kappainit)
  
  tempfilename = tempfile()
  myc<-file(tempfilename, open="w")
  cat(dpeaqms.model,file=myc)
  close(myc)

  dpeaqms.jags.model <- jags.model(file=tempfilename, data=modeldata, inits=modelinits)

  # Run the MCMC chain for the specified number of burnin iterations
  update(dpeaqms.jags.model,burnin)
  # Sample the required number of samples from the MCMC chain
  dpeaqms.samples <- coda.samples(dpeaqms.jags.model,
                                  c("Beta", "Gamma" , "Sigma", "Alpha", "p", "kappa"),
                                  n.iter=samples*thin, thin=thin)
  # Return the results
  return(dpeaqms.samples)
  
}
