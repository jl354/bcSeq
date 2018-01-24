trimRead <- function(inputFile, outputFile, start, end) {
    invisible(.Call('_bcSeq_trimRead', PACKAGE = 'bcSeq', inputFile, outputFile, start, end))
}


