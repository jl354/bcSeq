trim <- function(inputFile, outputFile, start, end) {
    invisible(.Call('_bcSeq_trim', PACKAGE = 'bcSeq', inputFile, outputFile, start, end))
}

uniqueBar <- function(inputFile, outputFile) {
    invisible(.Call('_bcSeq_uniqueBar', PACKAGE = 'bcSeq', inputFile, outputFile))
}

