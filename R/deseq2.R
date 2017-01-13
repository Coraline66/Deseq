#' DESeq2 Differential Analysis
#'
#' This function performs DESeq2 gene differential analysis to
#' htseq counts files.
#' @import mydata
#' @import DESeq2
#' @param n integer. The number of samples for subsequent analysis.
#'          Note that n equals 2 times the number of columns of matrix samp for 
#'          datatype "luad" & "brca", while it equals the number of columns
#'          of matrix samp for datatype "sarc". 
#' @param i  integer. The row index of matrix samp which is obtained from
#'           the first-step or second-step resampling. Any integer is acceptable
#'           if an analysis for the entire sample group is performed.
#' @param samp  matrix obtained from first-step or second-step resampling.
#'              It can be set as an empty matrix if analysis for the entire
#'              sample group is desired.
#' @param datatype  character indicating which type of data is involved. 
#'                  "luad", "brca" or "sarc".
#' @keywords  DESeq2
#' @return a list consisted of the result of DESeq2 analysis and filtered
#'         gene list. A criterion of abs(log2-FoldChange) >  1 and
#'         adjusted P-value < 0.05 is utilized during the filtering process.
#' @export
#' @examples
#' deseq2(10, 1, samp=samp, datatype="sarc")

deseq2 <- function(n, i, samp=matrix(), datatype) {
    # Load data
    datapath <- path.package(package="Deseq")
    files <- grep(datatype, list.files(paste0(datapath, "/data", sep="")), 
                  value=TRUE)
    dir <- paste0(datapath, "/data", sep="")
    sampleFiles <- list()
    # Construct sample counts matrix
    if ((datatype == "luad")|(datatype == "brca")) {
        if (n == length(files)) {
            sampleFiles <- files
            sampleCondition <- factor(rep(c("Normal", "Tumor"),
                                          each = (n/2)))
            sampleTable <- data.frame(sampleName = sampleFiles,
                                      filename = sampleFiles,
                                      condition = sampleCondition)
        }
        if (n < length(files)) {
            for (j in 1:(dim(samp)[2])) {
                sampleFiles <- append(sampleFiles, files[as.numeric(samp[i, j])])
            }
            for (j in 1:(dim(samp)[2])) {
                sampleFiles <- append(sampleFiles, files[as.numeric(samp[i, j])
                                                         +length(files)/2])
            }
            sampleFiles <- as.character(sampleFiles)
            sampleFiles <- sort(sampleFiles)
            sampleCondition <- factor(rep(c("Normal", "Tumor"),
                                          each = n/2))
            sampleTable <- data.frame(sampleName = sampleFiles,
                                      filename = sampleFiles,
                                      condition = sampleCondition)
        }
    }
    
    if (datatype == "sarc") {
        if (n == length(files)) {
            sampleFiles <- files
            sampleCondition <- factor(c(rep("DLP", length(grep("dl", files))),
                                        rep("LMS", length(grep("lm", files)))))
            sampleTable <- data.frame(sampleName = sampleFiles,
                                      filename = sampleFiles,
                                      condition = sampleCondition)
        }
        if (n < length(files)) {
            for (j in 1:(dim(samp)[2])) {
                sampleFiles <- append(sampleFiles, files[as.numeric(samp[i, j])])
            }
            sampleFiles <- as.character(sampleFiles)
            sampleFiles <- sort(sampleFiles)
            sampleCondition <- factor(rep(c("DLP", "LMS"),
                                          each = (n/2)))
            sampleTable <- data.frame(sampleName = sampleFiles,
                                      filename = sampleFiles,
                                      condition = sampleCondition)
        }
    }
    
    # Construct DESeq Dataset from the Matrix
    ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                           directory = dir,
                                           design = ~ condition)
    
    # Pre-filtering
    ddsHTSeq_fil <- ddsHTSeq[rowSums(counts(ddsHTSeq)) > 1,]
    
    # Differential Gene Expression Test
    ddsHTSeq_fil <- DESeq(ddsHTSeq_fil)
    res <- results(ddsHTSeq_fil)
    
    # Find out the DEGs with an Absolute Log2-Fold Change > 1
    # and an Adjusted P-value < 0.05
    # First We Remove the Rows with "NA" Adjusted P-values
    res_rmna <- res[!is.na(res$padj), ]
    res_padj <- res_rmna$padj < 0.05
    res_PADJ <- res_rmna[res_padj, ]
    res_padj_ordered <- res_PADJ[order(res_PADJ$padj), ]
    gene_padj <- rownames(res_padj_ordered)
    res_fc <- abs(res_rmna$log2FoldChange) > 1
    res_FC <- res_rmna[res_fc,]
    res_fc_ordered <- res_FC[order(abs(res_FC$log2FoldChange),
                                   decreasing = TRUE), ]
    gene_fc <- rownames(res_fc_ordered)
    res_ts <- res_rmna[res_fc & res_padj, ]
    ts_des <- rownames(res_ts)
    Nts_des <- length(ts_des)
    
    # Return Results
    if (Nts_des == 0) {
        ls <- list("genes"="empty", "number"=0)
    } else {
        ls <- list("result"=res, "genes"=ts_des, "number"=Nts_des)
    }
    return(ls)
}