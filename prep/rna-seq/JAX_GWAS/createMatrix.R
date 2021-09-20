dir <- "/ssd/cory/AD_data/mouse_reads/JAX_GWAS"
files <- list.files(dir, pattern="*.txt")
length(files)  # 36

full.paths <- file.path(dir, files)


extractID <- function(s){
   match <- gregexpr("-[1-9][0-9][0-9][0-9]-", s)[[1]]
   stopifnot(match > 0)
   start <- as.numeric(match) + 1
   size <- attr(match, 'match.length')
   substring(s, start, start + size - 3)
   }

readSingleFile <- function(filename){
    full.path <- file.path(dir, filename)
    tbl <- read.table(full.path, sep="\t", nrow=-1, as.is=TRUE)
    mtx <- as.matrix(tbl[,2])
    rownames(mtx) <- tbl[,1]
    sampleId <-extractID(filename)
    sampleName <- sprintf("JAX_GWAS.%s", sampleId)
    colnames(mtx)[1] <- sampleName
    mtx
    }

count.matrices <- lapply(files, readSingleFile)
mtx <- do.call(cbind, count.matrices)
dim(mtx)
length(which(is.na(mtx)))

#------------------------------------------------------------
# make sure we have metadata for every sample
#------------------------------------------------------------
metadata.file <- "gwas_individual_metadata.csv"
tbl.md <- read.table(metadata.file,  sep=",", as.is=TRUE, nrow=-1, header=TRUE)
dim(tbl.md)
rownames(tbl.md) <- sprintf("JAX_GWAS.%s", tbl.md$individualID)

stopifnot(all(rownames(tbl.md) %in% colnames(mtx)))
stopifnot(all(colnames(mtx) %in% rownames(tbl.md)))

save(mtx, tbl.md, file="JAX_GWAS_matrixAndMetadata_readyForNoa.RData")
write.table(as.data.frame(mtx), file="JAX_GWAS_counts.tsv", sep="\t", quote=FALSE)
write.table(tbl.md, file="JAX_GWAS_metadata.tsv", sep="\t", quote=FALSE)

