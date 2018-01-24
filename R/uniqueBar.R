uniqueBar <- function(inputFile, outputFile) {
    invisible(.Call('_bcSeq_uniqueBar', PACKAGE = 'bcSeq', inputFile, outputFile))
}

