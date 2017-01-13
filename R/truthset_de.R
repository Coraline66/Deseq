#' Truth Set Redefinition
#'
#' This function imports a matrix attained from the first-step resampling,
#' compares the DESeq results with that of the original truth set and exports
#' a group of new truth sets(gene lists).
#' Note that N should be equal to the number of rows of matrix sam,
#' n should be equal to the 2 times the number of columns of matrix sam("luad",
#' "brca") and 1 time the number columns of matrix sam("sarc").
#' @param N  integer. The number of iterations for resampling the original
#'           truth set.
#' @param n  integer. The number of samples in the new truth set. 
#' @param sam  matrix obtained from first-step resampling.
#' @param datatype  character. "luad", "brca" or "sarc".
#' @description The function is used to obtain a new truth set for different 
#'              data types via a combination of certain variables.
#' @keywords  DESeq2
#' @return a list which is composed of N new truth sets for N iterations, 
#'         respectively.
#' @export
#' @examples
#' truthset_de(100, 170, sam, datatype="brca")

truthset_de <- function(N, n, sam, datatype) {
    # Apply DESeq2 differential analysis to the entire data set to get original
    # truth set.
    dir <- path.package(package="Deseq")
    files <- grep(datatype, list.files(paste0(dir, "/data", sep="")), value=TRUE)
    en_de <- deseq2(n=length(files), i=0, samp=matrix(), datatype=datatype)
    
    # Pick n out of the entire sample group for N iterations to construct the
    # new truth set.
    in_de <- list()
    for (i in 1:N) {
        in_de[[i]] <- deseq2(n, i, samp=sam, datatype)
    }
    
    # Find out the new truth sets for each iteration
    ts_de <- vector(mode="list")
    for (i in 1:N) {
        ts_de[[i]] <- intersect(en_de$genes, in_de[[i]]$genes)
    }
    
    # Return result
    return(ts_de)
}