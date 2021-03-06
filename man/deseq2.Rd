% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deseq2.R
\name{deseq2}
\alias{deseq2}
\title{DESeq2 Differential Analysis}
\usage{
deseq2(n, i, samp = matrix(), datatype)
}
\arguments{
\item{n}{integer. The number of samples for subsequent analysis.
Note that n equals 2 times the number of columns of matrix samp for 
datatype "luad" & "brca", while it equals the number of columns
of matrix samp for datatype "sarc".}

\item{i}{integer. The row index of matrix samp which is obtained from
the first-step or second-step resampling. Any integer is acceptable
if an analysis for the entire sample group is performed.}

\item{samp}{matrix obtained from first-step or second-step resampling.
It can be set as an empty matrix if analysis for the entire
sample group is desired.}

\item{datatype}{character indicating which type of data is involved. 
"luad", "brca" or "sarc".}
}
\value{
a list consisted of the result of DESeq2 analysis and filtered
        gene list. A criterion of abs(log2-FoldChange) >  1 and
        adjusted P-value < 0.05 is utilized during the filtering process.
}
\description{
This function performs DESeq2 gene differential analysis to
htseq counts files.
}
\examples{
deseq2(10, 1, samp=samp, datatype="sarc")
}
\keyword{DESeq2}
