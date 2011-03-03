# R function that extracts MCMC sampled parameters from the
# output of the dpeaqms.mcmc function
dpeaqms.extract.sample<-function(msmsdata=NULL, mcmc_sample, numberOfMCMCSamples,
                                 controlGroup=NULL,  referenceSampleID=NULL,
	                         summaryOnly=F, outputprefix=NULL) {

 experimentID = msmsdata$experiment
  proteinID    = msmsdata$protein
  intensity    = msmsdata$intensity
  peptideID    = paste("Exp" , experimentID, ".", msmsdata$peptide, sep = "")
  sampleID     = paste("Exp" , experimentID, ".", msmsdata$sample , sep = "")
  groupID      = msmsdata$group
  

  E = length(unique(experimentID))
 
 
 
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
        print(paste("Error: Unable to find specified reference sample"))
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

  # If the full MCMC output is required for all parameters
  if (!summaryOnly) {
    pSampleString = paste("mcmc_sample[[1]][1:" , numberOfMCMCSamples," ,'p']", sep ='')
    # Print out the MCMC output for all parameters associated with the proteins in the datafile
    for (i in seq(1,P)) {
      if (!is.null(outputprefix)) {
        proteinoutputname = paste(outputprefix, "." , proteins[i], ".samples",sep='')
      }
      else {
        proteinoutputname = paste(proteins[i], ".samples", sep="")
      }
      PSampleString = ""
      kappaSampleString = "" 
      kappaString = "" 
      BetaSampleString = ""
      GammaSampleString = ""
      
      for (g in seq(2,G)) {        
        if (P > 1) {
           pGroupSampleString = paste("mcmc_sample[[1]][1:" , numberOfMCMCSamples," ,'p[", g,",",i , "]']", sep ='')           
           betaGroupSampleString  = paste("mcmc_sample[[1]][1:" , numberOfMCMCSamples," ,'Beta[", g,",",i , "]']", sep ='')
           gammaGroupSampleString = paste("mcmc_sample[[1]][1:" , numberOfMCMCSamples," ,'Gamma[", g,",",i , "]']", sep ='')
        }
        else {
           pGroupSampleString = paste("mcmc_sample[[1]][1:" , numberOfMCMCSamples," ,'p[", g,"]']", sep ='')
           betaGroupSampleString  = paste("mcmc_sample[[1]][1:" , numberOfMCMCSamples," , 'Beta[", g, "]']", sep ='')
           gammaGroupSampleString = paste("mcmc_sample[[1]][1:" , numberOfMCMCSamples," ,'Gamma[", g, "]']", sep ='')
        }
        if (BetaSampleString == "") {
           PSampleString = paste("p.",groupLevels[g],".",proteins[i], "=eval(", pGroupSampleString, ")", sep="")
           BetaSampleString  = paste("Beta.",groupLevels[g],
                                     "=eval(",betaGroupSampleString, ")", sep="")
           GammaSampleString = paste("Gamma.",groupLevels[g],
                                     "=eval(",gammaGroupSampleString, ")", sep="")
        }
        else {
           PSampleString = paste(PSampleString, ",p.",groupLevels[g],"." , proteins[i],"=eval(", pGroupSampleString, ")", sep="")
           BetaSampleString  = paste(BetaSampleString,",Beta.",groupLevels[g],
                                     "=eval(",betaGroupSampleString, ")", sep="")
           GammaSampleString = paste(GammaSampleString,",Gamma.",groupLevels[g],
                                     "=eval(",gammaGroupSampleString, ")", sep="")
        }               
      }
           
      SigmaSampleString = paste("mcmc_sample[[1]][1:" , numberOfMCMCSamples," ,'Sigma']", sep ='')      
      if (i == P) {
       # print(paste("Last Protein" ,proteins[i]))
        peptideInstances = unique(peptideID[(offset[i]+1):length(peptideID)])
       #print(peptideInstances)
      }
      else {
        peptideInstances = unique(peptideID[(offset[i]+1):(offset[i+1])])
      }
      ppeptides = length(peptideInstances)   
    
      blockEffectString = paste("Alpha." , peptideNames[peptideInstances[1]] , "=" , "mcmc_sample[[1]][1:" , numberOfMCMCSamples,",'Alpha[", peptideInstances[1], "]']", sep='')
      if (ppeptides > 1) {
        for (p in seq(2,ppeptides)) {
          blockEffectString = paste(blockEffectString , ",Alpha." ,peptideNames[peptideInstances[p]]  , "=mcmc_sample[[1]][1:" , numberOfMCMCSamples," ,'Alpha[", peptideInstances[p], "]']", sep='')
        }
      }
          
      proteinoutputsamplesString = paste("data.frame(" , PSampleString, "," ,  BetaSampleString , "," ,GammaSampleString, ", Sigma=eval(parse(text=SigmaSampleString)),", blockEffectString , ")", sep='')    
      print(proteinoutputname)      
      proteinoutputsamples = eval(parse(text=proteinoutputsamplesString))  
      write.table(proteinoutputnumberOfMCMCSamples,quote=FALSE, sep="\t", row.names=FALSE, file=proteinoutputname)
    }

   
    for (e in seq(1,E)) {     
     for (s in seq(1,numberOfSamples[e])) {
         if (E == 1) {
           kappaSampleString = paste("mcmc_sample[[1]][1:" , numberOfMCMCSamples," ,'kappa[" , s, "]']", sep ='')
         }
         else {
           kappaSampleString = paste("mcmc_sample[[1]][1:" , numberOfMCMCSamples," ,'kappa[" , e, "," , s, "]']", sep ='')
         }
         if ((s == 1) && (e==1)){
           kappaString = paste("kappa.", sampleLevels[sampleoffset[e]+s] , "=eval(" , kappaSampleString , ")" , sep="") ;
         }
         else {
           kappaString = paste(kappaString,",kappa.", sampleLevels[sampleoffset[e]+s] , "=eval(" , kappaSampleString , ")" , sep="") ; ;
         }
        
     }
    }
    kappaOutputString = paste("data.frame(" , kappaString,  ")", sep='')       
    kappaOutput = eval(parse(text=kappaOutputString))     
    if (!is.null(outputprefix)) {
      kappaOutputname = paste(outputprefix, ".kappa.samples",sep='')     
    }
    else {
      kappaOutputname = paste("kappa.samples", sep="")     
    }
    write.table(kappaOutput, quote=FALSE, sep="\t", row.names=FALSE, file=kappaOutputname) ;
  }

  # Write a summary of the posterior output to a file
  meanBeta  = matrix(0.0, nrow=G-1 , ncol=P)
  meanBetaGamma = matrix(0.0, nrow=G-1 , ncol=P)
  stdBetaGamma  = matrix(0,0, nrow=G-1 , ncol=P)
  upregulatedWrtCtl = matrix(groupLevels[1] , nrow=G-1,ncol=P)
 
  for (i in seq(1,P)) {
    theline = proteins[i]
    for (g in seq(2,G)) {  
      if (P > 1) {
      betaGammaEvalString =  paste("mean(mcmc_sample[[1]][1:", numberOfMCMCSamples,",'Beta[" , g,",",i , "]'] * ", 
                                        "mcmc_sample[[1]][1:",numberOfMCMCSamples,",'Gamma[" , g, ",",i , "]']" , ")", sep ='')

      betaEvalString   = paste("mean(mcmc_sample[[1]][1:", numberOfMCMCSamples,",'Beta[" , g,",",i , "]'])", sep ='')      
      }
      else {
      betaGammaEvalString =  paste("mean(mcmc_sample[[1]][1:", numberOfMCMCSamples,",'Beta[" , g , "]'] * ", 
                                        "mcmc_sample[[1]][1:",numberOfMCMCSamples,",'Gamma[" , g,  "]']" , ")", sep ='')

      betaEvalString   = paste("mean(mcmc_sample[[1]][1:", numberOfMCMCSamples,",'Beta[" , g, "]'])", sep ='')   
      }
      meanBeta[g-1,i]  = eval(parse(text=betaEvalString))
      meanBetaGamma[g-1,i] = eval(parse(text=betaGammaEvalString))
      theline = paste(theline, "\t" , meanBeta[g-1,i])
      theline = paste(theline, "\t" , meanBetaGamma[g-1,i])
      if (P > 1) {
      betaGammaEvalString =  paste("sd(mcmc_sample[[1]][1:", numberOfMCMCSamples,",'Beta[" , g,",",i , "]'] * ", 
                                        "mcmc_sample[[1]][1:",numberOfMCMCSamples,",'Gamma[" , g, ",",i , "]']" , ")", sep ='')
      }
      else {
      betaGammaEvalString =  paste("sd(mcmc_sample[[1]][1:", numberOfMCMCSamples,",'Beta[" , g, "]'] * ", 
                                        "mcmc_sample[[1]][1:",numberOfMCMCSamples,",'Gamma[" , g, "]']" , ")", sep ='')
      } 
      stdBetaGamma[g-1,i]  = eval(parse(text=betaGammaEvalString))
      theline = paste(theline, "\t" , stdBetaGamma[g-1,i])
      if (meanBetaGamma[g-1,i] > 0.0) {
        upregulatedWrtCtl[g-1,i] = groupLevels[g]
      }      
      theline = paste(theline, "\t" , upregulatedWrtCtl[g-1,i])
    }
 
  }


   resultsEvalString = "data.frame(Protein=proteins"   
   for (g in seq(2,G)) {
      resultsEvalString = paste(resultsEvalString, ",",
                                "mean.Beta.",  groupLevels[g], "=meanBeta[", g-1 , ",]," ,
                                "mean.Beta.x.Gamma.", groupLevels[g], "=meanBetaGamma[",g-1 , ",]," ,
                                "std.Beta.x.Gamma.",  groupLevels[g], "=stdBetaGamma[",g-1 , ",]," ,
                                "upregulated.", groupLevels[g], "=upregulatedWrtCtl[",g-1 , ",]", sep='')
                                    
    }
    resultsEvalString =paste(resultsEvalString,")",sep='')
    results = eval(parse(text=resultsEvalString))

    myord = order(results[,2], decreasing=T)
    if (!is.null(outputprefix)) {
       resultsoutputfile = paste(outputprefix, ".results", sep = "")
    }
    else {
       resultsoutputfile = "Run.Results"
    }
    write.table(results[myord,],quote=FALSE, sep="\t", row.names=FALSE, file=resultsoutputfile)
}
