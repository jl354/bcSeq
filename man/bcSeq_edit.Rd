\name{bcSeq_edit}
\alias{bcSeq_edit}
\title{Function to perform reads to barcode alignment tolerating edit
type error distance.}
\usage{
bcSeq_edit(sampleFile, libFile, outFile, misMatch = 2, tMat =
NULL, numThread = 4, count_only = TRUE, gap_left = 3,
ext_left = 1, gap_right = 3, ext_right = 1, pen_max =
6, userProb = NULL, detail_info = FALSE)
}
\arguments{
\item{sampleFile}{(string) sample filename, needs to be a fastq file.}

\item{libFile}{(string) library filename, needs to be a fasta or fastq file.}

\item{outFile}{(string) output filename.}

\item{misMatch}{(integer) the number of maximum mismatches or indels allowed in
the alignment.}

\item{tMat}{(two column dataframe) prior probability of a mismatch given a
sequence. The first column is the prior sequence, the second column is the
error rate. The default value for all prior sequences is 1/3.}

\item{numThread}{(integer) the number of threads for parallel computing,
default 4.}

\item{count_only}{(bool) option for controlling function returns, default to be
\code{TRUE}.  If set to \code{FALSE}, a list contains a alignment probability
matrix between all the reads and barcodes, a read IDs vector, and barcode IDs
vector will be returned. The row of the matrix is corresponding to the read IDs
and the column of the matrix is associated with the barcode IDs. Examples of
the probability matrix are provided in the vignettes file.}

\item{gap_left}{(double) Penalty score for delete a base for the reads.}

\item{ext_left}{(double) Penalty score for extending deletion of base for the
reads.}

\item{gap_right}{(double) Penalty score for delete a base for the barcodes.}

\item{ext_right}{(double) Penalty score for extending deletion of base for the
barcodes.}

\item{pen_max}{(double) Max penalty allowed for a alignment.}

\item{userProb}{ (function) a function to compute the alignment probability
with 3 arguments \code{userProb(max_pen, prob, pen_val)}, \code{max_pen} is
the max penalty allowed, \code{prob} is a vector the probability for match
and mismatch part between a read and a barcode for all the for all the possible
alignment forms (since there are multiple way to align a read to a barcode for
same edit distance), \code{pen_val} is a vector for the value of penalty for
all the possible alignment forms between a read and a barcode.  The purpose of
\code{userProb(max_pen, prob, pen_val)} is to provide a way to determine the
alignment pattern and probability.}

\item{detail_info}{(bool) option for controlling function returns, default to
be \code{FALSE}. If set to \code{TRUE}, a file contain read indexes and library
indexes reads aligned will be created with file name \code{\$(outFile).txt}. Not
available for user-defined probability model case.}
}

\note{ The user need to perform the removing of any adapter sequence before and
after the barcode for the fastq file.  }

\value{
\item{default}{ No objects are returned to R, instead, a csv count table is
created and written to files.  The .csv file contains two columns, the first
column is the sequences of the barcodes, and the second columns is the number
of reads that aligned to the barcodes.}

\item{count_only = FALSE}{ If set to \code{FALSE}, bcSeq will return list
contains a sparse matrix for alignment probabilities, a read IDs vector, and a
barcode IDs vector.  The rows of the matrix are corresponding to the read IDs
vector, and the columns is associated with the barcode IDs vector.} }

\description{
This a function for aligning CRISPR barcode reads to library, or similar
problems using edit distance to evaluate the distance between a read and a
barcode for the error toleration.  }

\examples{
#### Generate barcodes
lFName    <- "./libFile.fasta"
bases     <- c(rep('A', 4), rep('C',4), rep('G',4), rep('T',4))
numOfBars <- 20
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
outFile  <- "./count_edit.csv"
res <- bcSeq_edit(rFName, lFName, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = TRUE, userProb = NULL,
    gap_left = 2, ext_left = 1, gap_right = 2, ext_right = 1,
    pen_max = 7)

#### The user defined probability model function is a modeling 
#### of the alignment probability and the penalty encounted 
#### due to the edit distance between a read and a candidata 
#### barcode. The user define probability model function serves 
#### as a combined evaluation of the alignment quality by 
#### considering both the alignment probability and the edit 
#### distance penalty.

#### the user defined function has three arguements
#### (1) val: a double vector indicates the alignment probabilities 
####          between a read and its candidate barcodes based on 
####          edit distance.
#### (2) pens: a double vector indicates the penalties of edit distance
####          between a read and its candidate barcodes for the alignment
####          based on edit distance.
#### (3) m: a double scalor indicating the weight for alignment 
####        probabilities and edit distance penalty to determine 
####        the final alignment quality that can be used as 
####        as clasifier. 

#### User can also constuct a more complexed model by only keeping the 
#### function signature.

#### Example function in R
useP <-function(m, val, pens) { val * (1 - log(2) + log(1 + m / (m + pens) ) ) }

#### Example function in C++(can be ported to R using Rcpp packages)
library(Rcpp)
cppFunction(
'NumericVector cpp_fun(int m, NumericVector val, NumericVector pens) {
    int n = val.size();
    NumericVector out(n);
    for(int i = 0; i < n; ++i) {
        out[i] = val[i] * (1 - log(2) +
        log(1 + m / (m + pens[i]) ) );
    }
    return out;
}')
outFile  <- "./count_edit_2.csv"
res <- bcSeq_edit(rFName, lFName, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = TRUE, userProb = useP,
    gap_left = 2, ext_left = 1, gap_right = 2, ext_right = 1,
    pen_max = 7)
}
