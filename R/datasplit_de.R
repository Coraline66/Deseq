#' Save splitted data 
#'
#' This function loads data sets from package "mydata", splits them into 
#' separated text files and saves in "/data" subdirectory for DESeq2 analysis. 
#' @import mydata
#' @export
#' @examples
#' datasplit_de()

datasplit_de <- function() {
    data(luad, brca, sarc, files_luad, files_brca, files_sarc)
    dir <- path.package(package="Deseq")
    dir.create(path=paste0(dir, "/data", sep=""))
    for (i in 1:length(files_luad)) {
        write.table(luad[[i]], file=paste0(dir, "/data/", files_luad[i], sep=""),
                    row.names=FALSE, col.names=FALSE)
    }
    for (i in 1:length(files_brca)) {
        write.table(brca[[i]], file=paste0(dir, "/data/", files_brca[i], sep=""),
                    row.names=FALSE, col.names=FALSE)
    }
    for (i in 1:length(files_sarc)) {
        write.table(sarc[[i]], file=paste0(dir, "/data/", files_sarc[i], sep=""),
                    row.names=FALSE, col.names=FALSE)
    }
}