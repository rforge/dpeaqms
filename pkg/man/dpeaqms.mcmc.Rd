\name{dpeaqms.mcmc}
\alias{dpeaqms.mcmc}

\title{dpeaqms MCMC Sampling function}
 
\description{A function that will take an MS/MS data object and run an
  MCMC algorithm to sample from the posterior distribution of the model
  (see dpeaqms.model) parameters}

\usage{
  dpeaqms.mcmc(msmsdata, burnin, samples, thin,
               a.alpha, b.alpha, a.gamma,  b.gamma,
               a.tau, b.tau,a.kappa, b.kappa, 
               a.p, b.p, controlGroup, referenceSampleID, 
	       transform, infofilename, verbose, seed)
}

\arguments{
  \item{msmsdata}{A data frame containing the MS/MS data. The data frame
    should have the columns: "experiment", "protein", "msmsID", "intensity",
    "sample" and "group".}
  \item{burnin}{The number of initial MCMC iterations to be discarded.}
  \item{samples}{The number of samples to take from the MCMC chain.}
  \item{thin}{This is the thinning interval on the sampled MCMC chain
    i.e. the interval between consecutive samples}
  \item{a.alpha}{The mean of the Gaussian prior distribution on the
    model parameter alpha (See model description)}
  \item{b.alpha}{The precision of the Gaussian prior distribution on the
    model parameter alpha (See model description)}
  \item{a.gamma}{The mean of the Gaussian prior distribution on the
    model parameter gamma (See model description)}
  \item{b.gamma}{The precision of the Gaussian prior distribution on the
    model parameter alpha (See model description)}
  \item{a.tau}{The shape parameter of the Gamma prior distribution on the
    model parameter tau (See model description)} 
  \item{b.tau}{The rate parameter of the Gamma prior distribution on the
    model parameter tau (See model description)}
  \item{a.kappa}{The mean of the Gaussian prior distribution on the
    model parameter kappa (See model description)}
  \item{b.kappa}{The precison of the Gaussian prior distribution on the
    model parameter kappa (See model description)}
  \item{a.p}{First parameter of the Beta hyperprior distribution on the
    parameter p}
  \item{b.p}{Second parameter on the Beta hyperprior distribution on the
  the model parameter p}
  \item{controlGroup}{The label of the group to be used as the
    control. By default the first group label is used as the control}
  \item{referenceSampleID}{A list of labels of the samples to be
    used as the calibration channel (one for each experiment). If none
    is specified the first sample label in each experiment is used.}
  \item{transform}{If TRUE uses the log of the reporter ion intensities
    otherwise leaves the intensities unchanged}
  \item{infofilename}{Specifies the name of the file to which to output
    the parameters used in the MCMC run}
  \item{verbose}{Print out information on the run progress}
  \item{seed}{Always the user to set the seed of the RNG used by JAGS
    for reproducible output}
}

\value{A sample from the posterior distribution in coda compatible format}

\examples{
  # Load the dpeaqms library
  library(dpeaqms)

  # Load the simulated dataset that comes with the package
  simdata <- dataset(dpeaqms.simulatedDataset1)

  # Sample from the posterior distribution of the model
  simdataSamples <- dpeaqms.mcmc(msmsdata=simdata, transform=FALSE,
                                 samples=1000, thin=100, burnin=1000)
}
\keyword{mcmc, proteomics}
