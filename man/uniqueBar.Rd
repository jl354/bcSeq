\name{uniqueBar}
\alias{uniqueBar}
\title{Tool function to obtain unique barcode sequences from a library .}
\usage{
uniqueBar(inputFile,outputFile)
}
\arguments{
\item{inputFile}{(string) filename for the library sequences, needs to be a fasta or fastq file.}

\item{outFile}{(string) output filename.}
}

\value{
\item{default}{ No objects are returned to R}
}
\description{
This a function for removing the duplicated barcodes in the library file
that will be used for the bcSeq alignment.}

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
Barcodes <- rbind(Barcodes, Barcodes)
write(Barcodes, lFName)

outFile  <- "./libFile_unique.fasta"
uniqueBar(lFName, outFile)

}
