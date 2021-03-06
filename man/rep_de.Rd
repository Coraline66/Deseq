% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rep_de.R
\name{rep_de}
\alias{rep_de}
\title{DESeq2 Analysis of Resampled Subsets}
\usage{
rep_de(n, sam, samp, datatype)
}
\arguments{
\item{n}{integer. The number of samples for subsequent analysis.
Note that n equals 2 times the number of columns of matrix samp for 
datatype "luad" & "brca", while it equals the number of columns
of matrix samp(element of samp if it is a list) for datatype "sarc".}

\item{sam}{matrix obtained from first-step resampling which indicates the 
new truth sets.}

\item{samp}{matrix or a list of matrices obtained from second-step resampling. 
It consists of the indices of subsets generated from the new truth sets.}

\item{datatype}{character indicating which type of data is involved. 
"luad", "brca" or "sarc".}
}
\value{
a list that is composed of the results of subgroup DESeq2 analysis.
        A criterion of abs(log2-FoldChange) >  1 and adjusted P-value < 0.05 
        is utilized during the filtering process.
}
\description{
This function performs DESeq2 differential analysis to a group 
             of subsets of the samples to investigate the effect of replicates.
}
\details{
Note that the length of exported list should be equal to the number of rows 
of matrix sam times the number of rows of matrix samp(element of samp if it is
a list). The results for the subsets of the first new truth set are stored in
"rep_de(n, sam, samp, datatype)[[k]], k=1, 2, ..., dim(samp[[1]])[1]", i.e., 
the first "dim(samp[[1]])[1]" elements of the returned result of function 
"rep_de()". Similarly, "rep_de(n, sam, samp, datatype)[[k]], 
k=dim(samp[[1]])[1]+1, ..., 2*(dim(samp[[1]])[1])" contains the information
corresponding to the second new truth set, etc.
}
\examples{
rep_de(n=4, sam=sam, samp=samp, datatype="luad")
}
