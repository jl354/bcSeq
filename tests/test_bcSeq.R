library(bcSeq)
#devtools::load_all("../")
#### Set the seed
set.seed(4523)

#### Generate barcode
lFName    <- "./libFile.fasta"
bases     <- c(rep('A', 4), rep('C',4), rep('G',4), rep('T',4))
numOfBars <- 7
Barcodes  <- rep(NA, numOfBars*2)
for (i in 1:numOfBars){
    Barcodes[2*i-1] <- paste0(">barcode_ID: ", i)
    Barcodes[2*i]   <- paste(sample(bases, length(bases)), collapse = '')
}
write(Barcodes, lFName)

#### Generate reads and phred score
rFName     <- "./readFile.fastq"
numOfReads <- 8
Reads      <- rep(NA, numOfReads*4)
for (i in 1:numOfReads){
    Reads[4*i-3] <- paste0("@read_ID_",i)
    Reads[4*i-2] <- Barcodes[2*sample(1:numOfBars,1,
                  replace=TRUE, prob=seq(1:numOfBars))]
    Reads[4*i-1] <- "+"
    Reads[4*i]   <- paste(rawToChar(as.raw(
                  33+sample(20:30, length(bases),replace=TRUE))),
                  collapse='')
}
write(Reads, rFName)

#### perform alignment 
ReadFile <- "./readFile.fastq"
BarFile  <- "./libFile.fasta"
outFile  <- "./countH.csv"

#### with default output for bcSeq_hamming
res <- bcSeq_hamming(ReadFile, BarFile, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = TRUE )
res <- read.csv(outFile, header=FALSE)
res

#### with return of alignment probability matrix to R
outFile  <- "./countH2.csv"
res <- bcSeq_hamming(ReadFile, BarFile, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = FALSE )
res 

#### with default output for bcSeq_edit
outFile  <- "./countE.csv"
res <- bcSeq_edit(ReadFile, BarFile, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = TRUE,
    gap_left = 2, ext_left = 1, gap_right = 2, ext_right = 1,
    pen_max = 7)
res <- read.csv(outFile, header=FALSE)
res

#### with return of alignment probability matrix to R
outFile  <- "./countE2.csv"
res <- bcSeq_edit(ReadFile, BarFile, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = FALSE,
    gap_left = 2, ext_left = 1, gap_right = 2, ext_right = 1,
    pen_max = 7)
res

#### user-defined probability model
comstomizeP <- function(m, x, y)
{
    x * (1 - log(2) + log(1 + m / (m + y) ) )
}
outFile = "comstomizeP.csv"
bcSeq_edit(ReadFile, BarFile, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = TRUE,
    gap_left = 2, ext_left = 1, gap_right = 2, ext_right = 1,
    pen_max = 7, userProb = comstomizeP)
