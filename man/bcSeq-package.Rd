\name{bcSeq-package}
\alias{bcSeq-package}
\alias{bcSeq-package}
\docType{package}
\title{Fast Sequence Alignment for High-throughput shRNA and CRISPR Screens
}
\description{
This \code{Rcpp}-based package implements highly efficient data structure and
algorithms for performing the alignment of short 'CRISPR' or shRNA screens
reads to library barcodes based on user specified mismatch, insert and delete
tolerance. Matching qualities are then evaluated based on \code{Phred} score. A
Bayers' classifier is employed to determine the originating barcode of a read.
We support user-defined probability model for evaluting matching qualities as
well as flexible output.The alignment also support multiple-thread to reduce
the processing time in the \code{C++} implementation.
}
\details{
\tabular{ll}{
Package: \tab bcSeq\cr
Type: \tab Package\cr
<<<<<<< HEAD
Version: \tab 1.0.1\cr
Date: \tab 2017-12-15\cr
=======
Version: \tab 1.3.2\cr
Date: \tab 2018-08-02\cr
>>>>>>> d40261c... adjust version number
License: \tab GPL-3\cr
}
Please see the example function calls below, or refer to the individual
function documentation or the included vignette for more information. The
package vignette serves as a tutorial for using this package. 
}

\author{
Jiaxing Lin, Jeremy Gresham, Tongrong Wang, So Young Kim, James Alvarez,
Jeffrey S. Damrauer, Scott Floyd, Joshua Granek, Andrew Allen, Cliburn Chan,
Jichun Xie, Kouros Owzar

Maintainer: Jiaxing Lin <jiaxing.lin@duke.edu>
}

\keyword{ package }
\seealso{
    \code{\link{Rcpp}}
}

\examples{
#### Generate barcodes
lFName    <- "./libFile.fasta"
bases     <- c(rep('A', 4), rep('C',4), rep('G',4), rep('T',4))
numOfBars <- 15
Barcodes  <- rep(NA, numOfBars*2)
for (i in 1:numOfBars){
    Barcodes[2*i-1] <- paste0(">barcode_ID: ", i)
    Barcodes[2*i]   <- paste(sample(bases, length(bases)), collapse = '')
}
write(Barcodes, lFName)

#### Generate reads and phred score
rFName     <- "./readFile.fastq"
numOfReads <- 100
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
outFile  <- "./counthamming.csv"
res <- bcSeq_hamming(rFName, lFName, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = TRUE)
outFile  <- "./countedit.csv"
res <- bcSeq_edit(rFName, lFName, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = TRUE,
    gap_left = 2, ext_left = 1, gap_right = 2, ext_right = 1,
    pen_max = 7, userProb = NULL) 
}
