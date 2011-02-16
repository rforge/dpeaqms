\name{dpeaqms.extract.sample}
\alias{dpeaqms.extract.sample}

\title{Extract samples from the output of dpeaqms}

\description{This function will convert the output of the function
  dpeaqms.mcmc, i.e. a sample from the posterior distribution of model
  parameters, into a text file or series of text files if complete
  output is specified}


\usage{
dpeaqms.extract.sample(msmsdata, dpeaqms_mcmc_sample, 
                       controlGroup,  referenceSampleID,
	               samples, summaryOnly, outputprefix)
}

\arguments{
  \item{msmsdata}{The MS/MS data in tabular format consisting of one
    reporter ion per line with the fields corresponding to Protein,
    Intensity, Peptide Sequence , Replicate Number, Sample ID and Group
    ID in that order. Any additional fields will be ignored. NULL by default}
  \item{dpeaqms_mcmc_sample}{Sample out from dpeaqms.mcmc}
  \item{controlGroup}{The label of the group to be used as the
    control. By default the first group label is used as the control}
  \item{referenceSampleID}{A list of labels of the samples to be
    used as the calibration channel (one for each experiment). If none
    is specified the first sample label in each experiment is used.}
  \item{samples}{The number of samples to taken of the MCMC chain.}
  \item{summaryOnly}{If TRUE outputs all the sampled parameters}
  \item{outputprefix}{Prefix of the outputfile}
}

\keyword{mcmc, proteomics}
