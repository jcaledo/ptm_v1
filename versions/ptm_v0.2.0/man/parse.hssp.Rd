% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evolutionary.R
\name{parse.hssp}
\alias{parse.hssp}
\title{Parse a HSSP File to Return Dataframes}
\usage{
parse.hssp(file, keepfiles = TRUE)
}
\arguments{
\item{file}{input hssp file.}

\item{keepfiles}{logical, if TRUE the dataframes will be saved in the working directory and we will keep the hssp file.}
}
\value{
Returns a dataframe corresponding to the profile.Rda described above.
}
\description{
Parses a HSSP file to return dataframes.
}
\details{
If the argument 'keepfiles' is not set to TRUE, the hssp file used to get the parsed dataframe will be removed. Otherwise, 4 dataframes will be saved:
\itemize{
\item{id_seq_list.Rda:}  {This block of information holds the metadata per sequence, and some alignment statistic. See https://swift.cmbi.umcn.nl/gv/hssp for a detailed description of the information that can be find in this block.}
\item{id_aln.Rda:}  {This dataframe contains the alignment itself (each sequence is a column). Additional information such as secondary structure, SASA, etc., is also found in this block.}
\item{id_profile.Rda:}  {This dataframe holds per amino acid type its percentage in the list of residues observed at that position. In addition, this dataframe also informs about the entropy at each position, as well as the number of sequences spanning this position (NOOC).}
\item{id_insertions.Rda:} {A dataframe with information regarding those sequences that contain insertions. See https://swift.cmbi.umcn.nl/gv/hssp for further details.}
}
}
\examples{
\dontrun{parse.hssp(file = './1u8f.hssp')}
}
\references{
Touw et al (2015) Nucl. Ac. Res. 43:D364-368.
}
\seealso{
msa(), custom.aln(), list.hom(), get.hssp(), shannon()
}
\author{
Juan Carlos Aledo
}
