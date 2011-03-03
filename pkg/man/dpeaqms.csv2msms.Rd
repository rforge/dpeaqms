\name{dpeaqms.csv2msms}
\alias{dpeaqms.csv2msms}

\title{Converts a csv style file into an MSMS data frame that can be used by other functions in the package}
 
\description{This function will take a csv file and extract reporter ion intensity, an MSMS and a protein accession
information in order to create an MSMS data frame useable by the other functions in the dpeaqms package
}

\usage{
  dpeaqms.csv2msms(csvfile, proteinColumn, msmsidColumn,
                   reporterIonNames, sampleGroupAssignment,
                    includeUnknownProteins)

}

\arguments{
\item{csvfile}{The csvfile containing the reporter ion intensities, protein accession and msms information. The default is NULL.}
\item{proteinColumn}{The column number corresponding to the protein accessio. The default is 1 i.e. the first column.}
\item{msmsidColumn}{The column number(s) to use in order to construct. The default is none i.e. msmsIDs will simple be integers.}
\item{reporterIonNames}{A list containing the reporter ion intensity column names. By default these are those for a TMT 6-plex i.e. c(126,127,128,129,130,131)}
\item{sampleGroupAssignment}{The groups to which the samples labeled with the above reporter ions are assigned. By default assumes a TMT 6-plex and assigns the first three to a CTL group and the second three to a TRT group i.e. c("CTL","CTL","CTL","TRT","TRT","TRT")}.
\item{includeUnknownProteins}{This is a boolean which controls whether (TRUE) or not (FALSE) MSMS spectra with no known protein accessions are assigned to "Unknown" protein accession. The default is TRUE.}
}
