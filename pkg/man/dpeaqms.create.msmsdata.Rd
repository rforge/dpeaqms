\name{dpeaqms.create.msmsdata}
\alias{dpeaqms.create.msmsdata}

\title{Creates and MSMS data frame for use in the functions
  \link{dpeaqms.mcmc} and \link{dpeaqms.extract.sample}}

\description{This function converts reporter ion intensities obtained
  from MSMS experiments on isobaric labeled samples into a data frame
  useable by the other functions in the dpeaqms package. The function
  expects the text files to be in tabular format with the fields
  protein, intensity, peptide sequence, peptide replicate, sample
  identifier and group identifier.}


\usage{
dpeaqms.create.msmsdata(datafiles) 
}


\arguments{
\item{datafiles}{A list containing the names (including paths if
  needed) of the files containing the reporter ion intensities
  obtained from MSMS experiments on isobaric labeled samples. One data
  file per experiment is expected.}
}

\keyword{mcmc, proteomics}