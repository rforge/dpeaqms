# Converts a CSV file into an msms data frame usable by other functions that require 
# msms data
dpeaqms.csv2msms<-function(csvfiles=NULL, proteinColumn=1 , msmsidColumn=NULL,
                           reporterIonNames=c(126,127,128,129,130,131) ,                           
                           sampleGroupAssignment=c("A","A","A","B","B","B"),
                           includeUnknownProteins=T) {
 
  if (!is.null(csvfiles)) {
    proteinID    = c()
    intensity    = c()
    sampleID     = c()
    groupID      = c()
    msmsID       = c()
    experimentID = c()
    for (d in seq(1,length(csvfiles))) {
      csvdata <- read.table(csvfiles[d], header=T, sep=",", as.is=T)
      offset = length(experimentID)
      expsize = dim(csvdata)
      numObservations = expsize[1]*length(reporterIonNames)
      experimentID = c(experimentID,rep(d, numObservations))
      proteinID    = c(proteinID, rep("", numObservations))
      intensity    = c(intensity, rep(0.0,  numObservations))
      sampleID     = c(sampleID, rep("", numObservations))
      groupID      = c(groupID,rep("", numObservations))
      msmsID       = c(msmsID,rep("", numObservations))
      
      sampleNamesR = paste("X", reporterIonNames, sep="") ;
      ionIntensityColumns = match(sampleNamesR, colnames(csvdata))
   
      for (i in seq(1,length(ionIntensityColumns))) {
        if (is.nan(ionIntensityColumns[i])) {
          print(paste("Error: No column \"" , reporterIonNames[i] , "\" in " , csvfiles[d]))
          return(NULL)
        }
      }
   
      if (is.null(msmsidColumn)) {
        msmsName = seq(1,expsize[1])
      }
      else {
        msmsName = csvdata[,msmsidColumn[1]]
        for (z in seq(2,length(msmsidColumn))) {
          temp = csvdata[,msmsidColumn[z]]
          temp = gsub('[)$]', '' , temp)
          temp = gsub('[(|)|\"]', '.' , temp)
          temp = gsub('[ ]+','',temp)          
          msmsName = paste(msmsName, temp, sep=".")
        }
        msmsName = paste(msmsName, seq(1,expsize[1]), sep = ".")
      }
      proteinAccession = csvdata[,proteinColumn]
      proteinAccession = gsub('[-*+/]+', '.', proteinAccession)
    
      proteinAccession[proteinAccession==""] = "Unknown"
      experientID = rep(1,expsize[1])
      numSamples = length(reporterIonNames)
      for (i in seq(1,expsize[1])) {
        L = offset+(i-1)*numSamples+1
        U = offset+i*numSamples
      
        proteinID[L:U] = proteinAccession[i]
        for (j in seq(0,numSamples-1)) {
          intensity[L+j] = csvdata[i,ionIntensityColumns[j+1]]
        }
        sampleID[L:U]  = reporterIonNames
        groupID[L:U]   = sampleGroupAssignment
        msmsID[L:U]    = msmsName[i]
      }
    }
  }
  else {
    print("Error: No csvfiles specified")
    return (NULL) 
  }

  missingIntensity = is.na(intensity) | is.nan(intensity)
  experimentID = experimentID[!missingIntensity]
  proteinID    = proteinID[!missingIntensity]
  intensity    = intensity[!missingIntensity]
  sampleID     = sampleID[!missingIntensity]
  groupID      = groupID[!missingIntensity]
  msmsID       = msmsID[!missingIntensity]
  msmsdata<-data.frame("experiment"=experimentID,"protein"=proteinID,  "intensity"=intensity,
                       "sample"=sampleID, "group"=groupID, "msmsID"=msmsID)
  return(msmsdata)
}
