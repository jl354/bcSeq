\name{bcSeq_hamming}
\alias{bcSeq_hamming}
\title{Function to perform reads to barcode alignment tolerating hamming 
type error distance}
\usage{
bcSeq_hamming(sampleFile, libFile, outFile, misMatch=2, tMat = NULL, 
numThread = 4, count_only = TRUE, detail_info = FALSE)
}
\arguments{
\item{sampleFile}{(string) sample filename, needs to be a fastq file}

\item{libFile}{(string) library filename, needs to be a fasta or fastq file.}

\item{outFile}{(string) output filename.}

\item{misMatch}{(integer) the number of maximum mismatches or indels allowed in
the alignment.}

\item{tMat}{(two column dataframe) prior probability of a mismatch given a
sequence. The first column is the prior sequence, the second column is the
error rate. The default value for all prior sequences is 1/3.}

\item{numThread}{(integer) the number of threads for parallel computing,
default 4.}

\item{count_only}{(bool) option for function returns, default to be
\code{TRUE}.  If set to \code{FALSE}, a list contains a alignment probability
matrix between all the reads and barcodes, a read IDs vector, and barcode IDs
vector will be returned.  The row of the matrix is corresponding to the read
IDs and the column of the matrix is associated with the barcode IDs. Examples
of the probability matrix are provided in the vignettes files.}

\item{detail_info}{(bool) option for controlling function returns, default to
be \code{FALSE}. If set to \code{TRUE}, a file contain read indexes and library
indexes reads aligned will be created with file name \code{\$(outFile).txt}.}
}

\note{
The user need to perform the removing of any adaptor sequence before and after
the barcode for the fastq file.
}

\value{
\item{default}{ No objects are returned to R, instead, a csv count table is
created and written to files.  The .csv file contains two columns, the first
column is the sequences of the barcodes, and the second columns is the number
of reads that aligned to the barcodes.}

\item{count_only = FALSE}{ If set to \code{FALSE}, bcSeq will return list
contains a sparse matrix for alignment probabilities, a read IDs vector, and a
barcode IDs vector.  The rows of the matrix are corresponding to the read IDs
vector, and the columns is associated with the barcode IDs vector.} }

\description{ This a function for aligning CRISPR barcode reads to library, or
similar problems using hamming distance to evaluate the distance between a read
and a barcode for the error toleration.  }
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
outFile  <- "./count_hamming.csv"
#res <- bcSeq_hamming(rFName, lFName, outFile, misMatch = 2,
#    tMat = NULL, numThread = 1, count_only = TRUE)
}
