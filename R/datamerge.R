#' Merge TPR and PPV Data from Different Subgroups of Samples
#'
#' This function computes TPR(true positive rate) and PPV(positive predictive value)
#' according to the results we obtained from DESeq2 analysis.
#' @import parallel
#' @param ts_de a list of new truth sets attained from function "truthset_de".
#' @param in_de a list of results of subgroup DESeq2 analyses.
#' @return a list of vectors which contain corresponding TPR values and PPV values.
#' @export
#' @examples
#' datamerge(ts_de, in_de)

datamerge <- function(ls) {
    
}
# Combine data

files1 <- vector(mode = "list")
files2 <- vector(mode = "list")
files3 <- vector(mode = "list")
files4 <- vector(mode = "list")
case <- c("3.", "4.", "5.", "6.", "7.", "8.", "10.", "15.", "20.", "30.", 
          "40.", "50.", "60.", "70.", "80.")
for (i in 1:100) {
    files1[[i]] <- vector(mode = "list")
    files2[[i]] <- vector(mode = "list")
    files3[[i]] <- vector(mode = "list")
    files4[[i]] <- vector(mode = "list")
    for (j in 1:length(case)) {
        files1[[i]] <- append(files1[[i]], paste0("fdr_de", case[j], i, ".RData", sep=""))
        files2[[i]] <- append(files2[[i]], paste0("fdr_ed", case[j], i, ".RData", sep=""))
        files3[[i]] <- append(files3[[i]], paste0("sen_de", case[j], i, ".RData", sep=""))
        files4[[i]] <- append(files4[[i]], paste0("sen_ed", case[j], i, ".RData", sep=""))
    }
}
fdr_des <- vector(mode = "list")
fdr_edg <- vector(mode = "list")
sen_des <- vector(mode = "list")
sen_edg <- vector(mode = "list")
for (i in 1:100) {
    fdr_des[[i]] <- vector(mode = "list")
    fdr_edg[[i]] <- vector(mode = "list")
    sen_des[[i]] <- vector(mode = "list")
    sen_edg[[i]] <- vector(mode = "list")
    for (j in 1:length(case)) {
        setwd(path1)
        fdr_des[[i]] <- append(fdr_des[[i]], readRDS(files1[[i]][[j]])) 
        setwd(path2)
        fdr_edg[[i]] <- append(fdr_edg[[i]], readRDS(files2[[i]][[j]]))
        setwd(path3)
        sen_des[[i]] <- append(sen_des[[i]], readRDS(files3[[i]][[j]]))
        setwd(path4)
        sen_edg[[i]] <- append(sen_edg[[i]], readRDS(files4[[i]][[j]]))
    }
}
setwd("/Users/qux1/Documents/R/Rcode/BRCA100/data85_100")
for (i in 1:100) {
    fdr_de <- vector(mode = "list")
    fdr_ed <- vector(mode = "list")
    sen_de <- vector(mode = "list")
    sen_ed <- vector(mode = "list")
    for (j in 1:1500) {
        fdr_de[[j]] <- fdr_des[[i]][[j]]
        fdr_ed[[j]] <- fdr_edg[[i]][[j]]
        sen_de[[j]] <- sen_des[[i]][[j]]
        sen_ed[[j]] <- sen_edg[[i]][[j]]
    }
    saveRDS(fdr_de, file=paste0("fdr85_de100_", i, ".RData", sep=""))
    saveRDS(fdr_ed, file=paste0("fdr85_ed100_", i, ".RData", sep=""))
    saveRDS(sen_de, file=paste0("sen85_de100_", i, ".RData", sep=""))
    saveRDS(sen_ed, file=paste0("sen85_ed100_", i, ".RData", sep=""))
}
