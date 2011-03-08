dpeaqms.create.msmsdata<-function(datafiles) {
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
  

  # Each experiment should have its own data file
  # the argument "datafiles" is a list of experiment files
  #Nexperiments = length(datafiles)
  E = length(datafiles)
  # Read in the MSMS quantitative data
  for (i in seq(1,E)) {
     # Read the data to a temporary file
     temp <- read.table(datafiles[i], header=TRUE, as.is=TRUE)
     proteinID = c(proteinID,temp[[proteinColumn]])
     intensity = c(intensity, temp[[intensityColumn]])
     # Tag the peptideID with an experiment number prefix and a "replicate" prefix
     peptideID = c(peptideID , paste(temp[[peptideColumn]], "." , temp[[replicateColumn]], sep=""))
     # Tag the sample ID with an experiment number prefix
     sampleID  = c(sampleID, temp[[sampleColumn]])
     groupID = c(groupID, temp[[groupColumn]])
     # Keep track of the number of samples in each experiment     
     experimentID = c(experimentID, rep(i,length(temp[[proteinColumn]])))
  }

  msmsdata <- data.frame("experiment"=experimentID, "protein"=proteinID, "peptide"=peptideID ,
                         "intensity"=intensity , "sample"=sampleID , "group"=groupID)
  return(msmsdata)
}
