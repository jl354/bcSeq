\name{trim}
\alias{trim}
\title{Tool function to trim adaptor for read sequences.}
\usage{
trim(inputFile,outputFile, start, end)
}
\arguments{
\item{inputFile}{(string) filename for the library sequences, needs to be a fasta or fastq file.}

\item{outFile}{(string) output filename.}

\item{start}{(integer) starting position.}

\item{end}{(integer) ending position.}
}

\value{
\item{default}{ No objects are returned to R}
}
\description{
This a function for triming reads with adaptor. The sequencing within [start, end] will be written to 
the output file.}

\examples{
#### Generate barcodes
lFName    <- "./libFile.fasta"
bases     <- c(rep('A', 4), rep('C',4), rep('G',4), rep('T',4))
numOfBars <- 40
Barcodes  <- rep(NA, numOfBars*2)
for (i in 1:numOfBars){
    Barcodes[2*i-1] <- paste0(">barcode_ID: ", i)
    Barcodes[2*i]   <- paste(sample(bases, length(bases)), collapse = '')
}
write(Barcodes, lFName)

#### Generate reads and phred score
rFName     <- "./readFile.fastq"
numOfReads <- 800
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
outFile  <- "./readFile_trimed.fastq"
trim(rFName, outFile, 5,15)
}





