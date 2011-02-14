# R function that extracts MCMC sampled parameters from the
# output of the dpeaqms.mcmc function
dpeaqms.extract.sample<-function(datafile, dpeaqms_mcmc_sample, 
                                 controlGroup=NULL,  referenceSampleID=NULL,
	                         samples, summaryOnly=F, outputprefix=NULL) {
	
  # Read in data from the specified file(s)
  # This data file should have six columns
  # Column 1 = Protein Identifier
  # Column 2 = Peptide Identifier
  # Column 3 = Peptide Replicate Identifier
  # Column 4 = Sample Identifier
  # Column 5 = Group Identifier
  # Column 6 = Intensity Measurement

  proteinColumn   = 1
  intensityColumn = 2
  peptideColumn   = 3
  replicateColumn = 4
  sampleColumn    = 5
  groupColumn     = 6
  
  proteinID = c()
  intensity = c()
  peptideID = c() 
  replicateID = c() 
  sampleID = c()
  groupID = c() 
  experimentID = c()
  expSamples = c()

  # Each experiment should have its own data file
  # the argument "datafile" is a list of experiment files
  Nexperiments = length(datafile)
 
  # Read in the MSMS quantitative data
  for (i in seq(1,Nexperiments)) {
     # Read the data to a temporary file
     temp <- read.table(datafile[i], header=T, as.is=T)
     proteinID = c(proteinID,temp[[proteinColumn]])
     intensity = c(intensity, temp[[intensityColumn]])
     # Tag the peptideID with an experiment number prefix and a "replicate" prefix
     peptideID = c(peptideID , paste(temp[[peptideColumn]], ".Exp" , i , "." , temp[[replicateColumn]], sep=""))
     # Tag the sample ID with an experiment number prefix
     sampleID  = c(sampleID, paste("Exp", i , "." ,temp[[sampleColumn]], sep =""))
     groupID = c(groupID, temp[[groupColumn]])
     # Keep track of the number of samples in each experiment
     expSamples = c(expSamples, length(unique(temp[[sampleColumn]])))
     experimentID = c(experimentID, rep(i,length(temp[[proteinColumn]])))
  }
  
  # Log transform the intensity data if transform=T
 
 
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
  levels(experimentID) <-seq(1,Nexperiments)
  experimentID <- as.numeric(as.vector(experimentID))

  # Record the experiment offset i.e. where in the intensity vector
  # observations from the experiment start
  sampleoffset = rep(0,Nexperiments)
  if (Nexperiments > 1) {
   for (i in seq(2, Nexperiments)) {
     sampleoffset[i] = sampleoffset[i-1]+expSamples[i-1]
    }
  }
 
  
  # Factor ProteinID
  proteinID <-factor(proteinID, levels=unique.default(proteinID))
  # Store the orignal protein identifications
  pid = proteinID
  # Get the proteins (in the order in which they appear in the data file)  
  proteins  <- levels(proteinID)
  Nproteins <- length(proteins)
  # Numerically re-encode the proteins
  levels(proteinID)<-seq(1,length(levels(proteinID)))
  
  # Number of measurements for each protein and their "offset" in the list
  N <- vector(mode="integer",length=Nproteins) 
  offset <- vector(mode="integer" , length=Nproteins)
  
  # Factor SampleID
  sampleID    = factor(sampleID, levels=sort(unique.default(sampleID))) 
  sampleLevels = levels(sampleID)
  Nsamples = length(sampleLevels)
  
  # Numerically re-encode the sample identifiers
  samplenumericlevels <- seq(1,Nsamples)

  # If a reference sample (used for the sample normalization) is specified
  if (!is.null(referenceSampleID)) {
    soffset = 0 
    for (i in seq(1,Nexperiments)) {
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
      soffset = soffset + expSamples[i]
      
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
  Ngroups = length(groupLevels)
  numericlevels <- seq(1,Ngroups)
  
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
  offset[1] = 0
  if (Nproteins > 1) {
   for (i in seq(2,Nproteins)) {
     prot = proteins[i] 
     label = (pid==prot)
     N[i] <- sum(label)      
     offset[i] <- offset[i-1] + N[i-1]   
   }
  }


  # If the full MCMC output is required for all parameters
  if (!summaryOnly) {
    pSampleString = paste("dpeaqms_mcmc_sample[[1]][1:" , samples," ,'p']", sep ='')
    # Print out the MCMC output for all parameters associated with the proteins in the datafile
    for (i in seq(1,Nproteins)) {
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
      
      for (g in seq(2,Ngroups)) {        
        if (Nproteins > 1) {
           pGroupSampleString = paste("dpeaqms_mcmc_sample[[1]][1:" , samples," ,'p[", g,",",i , "]']", sep ='')           
           betaGroupSampleString  = paste("dpeaqms_mcmc_sample[[1]][1:" , samples," ,'Beta[", g,",",i , "]']", sep ='')
           gammaGroupSampleString = paste("dpeaqms_mcmc_sample[[1]][1:" , samples," ,'Gamma[", g,",",i , "]']", sep ='')
        }
        else {
           pGroupSampleString = paste("dpeaqms_mcmc_sample[[1]][1:" , samples," ,'p[", g,"]']", sep ='')
           betaGroupSampleString  = paste("dpeaqms_mcmc_sample[[1]][1:" , samples," , 'Beta[", g, "]']", sep ='')
           gammaGroupSampleString = paste("dpeaqms_mcmc_sample[[1]][1:" , samples," ,'Gamma[", g, "]']", sep ='')
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
           
      SigmaSampleString = paste("dpeaqms_mcmc_sample[[1]][1:" , samples," ,'Sigma']", sep ='')      
      if (i == Nproteins) {
       # print(paste("Last Protein" ,proteins[i]))
        peptideInstances = unique(peptideID[(offset[i]+1):length(peptideID)])
       #print(peptideInstances)
      }
      else {
        peptideInstances = unique(peptideID[(offset[i]+1):(offset[i+1])])
      }
      ppeptides = length(peptideInstances)   
    
      blockEffectString = paste("Alpha." , peptideNames[peptideInstances[1]] , "=" , "dpeaqms_mcmc_sample[[1]][1:" , samples,",'Alpha[", peptideInstances[1], "]']", sep='')
      if (ppeptides > 1) {
        for (p in seq(2,ppeptides)) {
          blockEffectString = paste(blockEffectString , ",Alpha." ,peptideNames[peptideInstances[p]]  , "=dpeaqms_mcmc_sample[[1]][1:" , samples," ,'Alpha[", peptideInstances[p], "]']", sep='')
        }
      }
          
      proteinoutputsamplesString = paste("data.frame(" , PSampleString, "," ,  BetaSampleString , "," ,GammaSampleString, ", Sigma=eval(parse(text=SigmaSampleString)),", blockEffectString , ")", sep='')    
      print(proteinoutputname)      
      proteinoutputsamples = eval(parse(text=proteinoutputsamplesString))  
      write.table(proteinoutputsamples,quote=FALSE, sep="\t", row.names=FALSE, file=proteinoutputname)
    }

   
    for (e in seq(1,Nexperiments)) {     
     for (s in seq(1,expSamples[e])) {
         if (Nexperiments == 1) {
           kappaSampleString = paste("dpeaqms_mcmc_sample[[1]][1:" , samples," ,'kappa[" , s, "]']", sep ='')
         }
         else {
           kappaSampleString = paste("dpeaqms_mcmc_sample[[1]][1:" , samples," ,'kappa[" , e, "," , s, "]']", sep ='')
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
  meanBeta  = matrix(0.0, nrow=Ngroups-1 , ncol=Nproteins)
  meanBetaGamma = matrix(0.0, nrow=Ngroups-1 , ncol=Nproteins)
  stdBetaGamma  = matrix(0,0, nrow=Ngroups-1 , ncol=Nproteins)
  upregulatedWrtCtl = matrix(groupLevels[1] , nrow=Ngroups-1,ncol=Nproteins)
 
  for (i in seq(1,Nproteins)) {
    theline = proteins[i]
    for (g in seq(2,Ngroups)) {  
      if (Nproteins > 1) {
      betaGammaEvalString =  paste("mean(dpeaqms_mcmc_sample[[1]][1:", samples,",'Beta[" , g,",",i , "]'] * ", 
                                        "dpeaqms_mcmc_sample[[1]][1:",samples,",'Gamma[" , g, ",",i , "]']" , ")", sep ='')

      betaEvalString   = paste("mean(dpeaqms_mcmc_sample[[1]][1:", samples,",'Beta[" , g,",",i , "]'])", sep ='')      
      }
      else {
      betaGammaEvalString =  paste("mean(dpeaqms_mcmc_sample[[1]][1:", samples,",'Beta[" , g , "]'] * ", 
                                        "dpeaqms_mcmc_sample[[1]][1:",samples,",'Gamma[" , g,  "]']" , ")", sep ='')

      betaEvalString   = paste("mean(dpeaqms_mcmc_sample[[1]][1:", samples,",'Beta[" , g, "]'])", sep ='')   
      }
      meanBeta[g-1,i]  = eval(parse(text=betaEvalString))
      meanBetaGamma[g-1,i] = eval(parse(text=betaGammaEvalString))
      theline = paste(theline, "\t" , meanBeta[g-1,i])
      theline = paste(theline, "\t" , meanBetaGamma[g-1,i])
      if (Nproteins > 1) {
      betaGammaEvalString =  paste("sd(dpeaqms_mcmc_sample[[1]][1:", samples,",'Beta[" , g,",",i , "]'] * ", 
                                        "dpeaqms_mcmc_sample[[1]][1:",samples,",'Gamma[" , g, ",",i , "]']" , ")", sep ='')
      }
      else {
      betaGammaEvalString =  paste("sd(dpeaqms_mcmc_sample[[1]][1:", samples,",'Beta[" , g, "]'] * ", 
                                        "dpeaqms_mcmc_sample[[1]][1:",samples,",'Gamma[" , g, "]']" , ")", sep ='')
      } 
      stdBetaGamma[g-1,i]  = eval(parse(text=betaGammaEvalString))
      theline = paste(theline, "\t" , stdBetaGamma[g-1,i])
      if (meanBetaGamma[g-1,i] > 0.0) {
        upregulatedWrtCtl[g-1,i] = groupLevels[g]
      }      
      theline = paste(theline, "\t" , upregulatedWrtCtl[g-1,i])
    }
  #  if (verbose) {
  #     print(theline)
  #  }
  }


   resultsEvalString = "data.frame(Protein=proteins"   
   for (g in seq(2,Ngroups)) {
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
