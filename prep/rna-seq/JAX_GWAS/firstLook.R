dir <- "/ssd/cory/AD_data/mouse_reads/JAX_GWAS"
files <- list.files(dir, pattern="*.txt")
length(files)  # 36

f <- file.path(dir, files[1])
file.exists(f)
print(files[1])
# "Apoe-KO-4699-GES15-02415_CAGATC_AC6GK3ANXX_LaneALL_001_STARAligned.sortedByCoord.out_htseq_count.txt"
read.table(f, sep="\t", as.is=TRUE, nrow=3)
#                   V1   V2
# 1 ENSMUSG00000000001 1296
# 2 ENSMUSG00000000003    0
# 3 ENSMUSG00000000028   77


files <- list.files(dir, pattern="*.results")
length(files)  # 36

f <- file.path(dir, files[1])
file.exists(f)
print(files[1])
# "Apoe-KO-4699-GES15-02415_CAGATC_AC6GK3ANXX_LaneALL_001_rsem.genes.results"
dim(read.table(f, sep="\t", as.is=TRUE, nrow=3, header=TRUE)) # 7 columns

read.table(f, sep="\t", as.is=TRUE, nrow=3, header=TRUE)[, c(1,3,4,5,6,7)]

#              gene_id  length effective_length expected_count   TPM  FPKM
# 1 ENSMUSG00000000001 3262.00          3068.30           1160 22.50 11.03
# 2 ENSMUSG00000000003  799.50           605.80              0  0.00  0.00
# 3 ENSMUSG00000000028 1825.93          1632.22             61  2.22  1.09


metadata.file <- file.path(dir, "SYNAPSE_METADATA_MANIFEST.tsv")
colnames(read.table(metadata.file,  sep="\t", as.is=TRUE, nrow=3, header=TRUE))

as.data.frame(t(read.table(metadata.file,  sep="\t", as.is=TRUE, nrow=3, header=TRUE)[, -c(1:3)]))

tbl.md <- read.table(metadata.file,  sep="\t", as.is=TRUE, nrow=-1, header=TRUE)

# do any of the count filenames appear in tbl.md?
# choose a constrained example to start
# wildcard Clu*count.txt
#  Clu-4643-GES15-02439_CAGATC-_AC6R3HANXX_LaneALL_001_STARAligned.sortedByCoord.out_htseq_count.txt
#  Clu-4644-GES15-02440_ACTTGA-_AC6R3HANXX_LaneALL_001_STARAligned.sortedByCoord.out_htseq_count.txt
#  Clu-4648-GES15-02441_GATCAG-_AC6R3HANXX_LaneALL_001_STARAligned.sortedByCoord.out_htseq_count.txt
#  Clu-4655-GES15-02442_TAGCTT-_AC6R3HANXX_LaneALL_001_STARAligned.sortedByCoord.out_htseq_count.txt
#  Clu-4658-GES15-02443_GGCTAC-_AC6R3HANXX_LaneALL_001_STARAligned.sortedByCoord.out_htseq_count.txt
#  Clu-4659-GES15-02444_CTTGTA-_AC6R3HANXX_LaneALL_001_STARAligned.sortedByCoord.out_htseq_count.txt

dim(subset(tbl.md, individualID=="4659"))

subset(tbl.md, individualID=="4659")# [ -c(1:3)]
 [1] "path"               
 [2] "parent"             
 [3] "name"               
 [4] "synapseStore"       
 [5] "contentType"        
 [6] "used"               
 [7] "executed"           
 [8] "activityName"       
 [9] "activityDescription"
[10] "assay"              
[11] "dataType"           
[12] "isModelSystem"      
[13] "individualID"       
[14] "sex"                
[15] "fileFormat"         
[16] "isStranded"         
[17] "readLength"         
[18] "tissue"             
[19] "species"            
[20] "specimenID"         
[21] "dataSubtype"        
[22] "organ"              
[23] "grant"              
[24] "resourceType"       
[25] "study"              
[26] "individualIdSource" 
[27] "consortium"         
[28] "libraryPrep"        
[29] "platform"           
[30] "runType"            

subset(tbl.md, individualID=="4659")# [ -c(1:3)]
soi <- grep("^465", tbl.md$individualID, v=TRUE)

  # 3 htseq_count files, 3 samples
tbl.small <- subset(tbl.md, individualID %in% soi & contentType=="text/plain")
dim(tbl.small)
wdth(3000)
id.column <- grep("individualID", colnames(tbl.small))
id.column
tbl.small[, c(4:12, 14:30)]
unique(tbl.small[, c(4:12, 14:30)])


tbl.md <- read.table("gwas_individual_metadata.csv", sep=",", as.is=TRUE, header=TRUE)


