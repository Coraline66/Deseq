#' Computing TPR and PPV
#'
#' This function computes TPR(true positive rate) and PPV(positive predictive value)
#' according to the results we obtained from DESeq2 analysis.
#' @import parallel
#' @param ts_de a list of new truth sets attained from function "truthset_de".
#' @param in_de a list of results of subgroup DESeq2 analyses.
#' @return a list of vectors which contain corresponding TPR values and PPV values.
#' @export
#' @examples
#' compute(ts_de, in_de)

compute <- function(ts_de, in_de) {
    # Use parallel processing to calculate TPR and PPV
    sen_de <- vector()
    ppv_de <- vector()
    no_cores <- detectCores() - 1
    cl <- makeCluster(no_cores)
    clusterExport(cl=cl, c("ts_de", "in_de", "sen_de", "ppv_de"), envir=environment())
    k <- 1:length(in_de)
    parLapply(cl, k,
              function(k) {
                  # Calculate the Numbers of "FP", "TP" for Each Gene List in "in_de".
                  # "TP" represents the intersection of new truth set and each genelist.
                  # "FP" is a resulting genelist subtracted by "TP".
                  tp <- intersect(in_de[[k]]$genes, ts_de)
                  fp <- setdiff(in_de[[k]]$genes, ts_de) 
                  Ntp <- length(tp)
                  Nfp <- length(fp)
                  # Compute Sensitivity(true positive rate)
                  sen_de[k] <- round(Ntp / length(ts_de), digits=3)
                  # Compute PPV = TP / (TP + FP)
                  ppv_de[k] <- round(Ntp / in_de[[k]]$number, digits=3)
                  if (in_de[[k]]$number == 0) {
                      ppv_de[k] <- 0
                  }
                  })
    stopCluster(cl)
    
    # Return result
    ls <- list("tpr" = sen_de, "ppv" = ppv_de)
    return(ls)
}