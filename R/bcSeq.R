.bcSeq_hamming <- function(sampleFile, libFile, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = TRUE) {
    # check file status
    if (!file.exists(sampleFile))
    stop(paste0("Error! ", sampleFile, " does not exist!"))
    if (!file.exists(libFile))
    stop(paste0("Error! ", libFile, " does not exist!"))
    if (file.exists(outFile))
    stop(paste0("Error! ", outFile, " exists plese specify
        another name."))

    tMatSeq <- c("default")
    tMatProb <- c(0.3333)
    if (!is.null(tMat))
    {
        for (i in 1:nrow(tMat))
        {
            tMatSeq[i] <- tMat[i, 1]
            tMatProb[i] <- tMat[i, 2]
        }
    }

    tMatSeq <- as.vector(tMatSeq)
    tMatProb <- as.vector(tMatProb)
    res <- invisible(.Call("_bcSeq_CRISPR_matching", PACKAGE = "bcSeq",
        sampleFile, libFile, outFile, misMatch, tMatSeq, tMatProb,
        numThread, TRUE, count_only, 0, 0, 0, 0, 0))

    if (!count_only)
    {
        out <- do.call(Matrix::sparseMatrix, res[[2]])
        c(res[1], alignProb = out)
    }

}

.bcSeq_edit <- function(sampleFile, libFile, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = TRUE, gap_left = 3,
    ext_left = 1, gap_right = 3, ext_right = 1,
    pen_max = 6, userProb = NULL) {
    # check file status
    if (!file.exists(sampleFile))
    stop(paste0("Error! ", sampleFile, " does not exist!"))
    if (!file.exists(libFile))
    stop(paste0("Error! ", libFile, " does not exist!"))
    if (file.exists(outFile))
    stop(paste0("Error! ", outFile,
        " exists plese specify another name."))

    tMatSeq <- c("default")
    tMatProb <- c(0.3333)
    if (!is.null(tMat))
    {
        for (i in 1:nrow(tMat))
        {
            tMatSeq[i] <- tMat[i, 1]
            tMatProb[i] <- tMat[i, 2]
        }
    }

    tMatSeq <- as.vector(tMatSeq)
    tMatProb <- as.vector(tMatProb)

    res = NULL
    if( !is.null(userProb) )
    {
        res <- invisible(.Call("_bcSeq_CRISPR_user_matching", 
            PACKAGE = "bcSeq", sampleFile, libFile, outFile,
            misMatch, tMatSeq, tMatProb, numThread, count_only, gap_left,
            ext_left, gap_right, ext_right, pen_max, userProb))
    } else
    {
        res <- invisible(.Call("_bcSeq_CRISPR_matching", PACKAGE = "bcSeq",
            sampleFile, libFile, outFile, misMatch, tMatSeq,
            tMatProb, numThread, FALSE, count_only, gap_left, ext_left,
            gap_right, ext_right, pen_max))
    }

    if (!count_only)
    {
        out <- do.call(Matrix::sparseMatrix, res[[2]])
        c(res[1], out)
    }
}

.bcSeq_hamming_DNAString <- function(sampleFile, libFile, outFile, 
    misMatch = 2, tMat = NULL, numThread = 4, count_only = TRUE) {
    readSeq <- as.vector(as.character(sampleFile))
    readSeq_ids <- as.vector(names(sampleFile))
    readPhred <- as.vector(as.character(unlist(sampleFile@metadata)))
    libSeq <- as.vector(as.character(libFile))
    libSeq_ids <- as.vector(names(libFile))

    tMatSeq <- c("default")
    tMatProb <- c(0.3333)
    if (!is.null(tMat))
    {
        for (i in 1:nrow(tMat))
        {
            tMatSeq[i] <- tMat[i, 1]
            tMatProb[i] <- tMat[i, 2]
        }
    }

    tMatSeq <- as.vector(tMatSeq)
    tMatProb <- as.vector(tMatProb)
    res <- invisible(.Call("_bcSeq_CRISPR_matching_DNAString", 
        PACKAGE = "bcSeq",readSeq, readSeq_ids, readPhred,libSeq,libSeq_ids,
        outFile, misMatch, tMatSeq, tMatProb,
        numThread, TRUE, count_only, 0, 0, 0, 0, 0))

    if (!count_only)
    {
        out <- do.call(Matrix::sparseMatrix, res[[2]])
        c(res[1], alignProb = out)
    }

}

.bcSeq_edit_DNAString <- function(sampleFile, libFile, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = TRUE, gap_left = 3,
    ext_left = 1, gap_right = 3, ext_right = 1,
    pen_max = 6, userProb = NULL) {

    readSeq <- as.vector(as.character(sampleFile))
    readSeq_ids <- as.vector(names(sampleFile))
    readPhred <- as.vector(as.character(unlist(sampleFile@metadata)))
    libSeq <- as.vector(as.character(libFile))
    libSeq_ids <- as.vector(names(libFile))

    tMatSeq <- c("default")
    tMatProb <- c(0.3333)
    if (!is.null(tMat))
    {
        for (i in 1:nrow(tMat))
        {
            tMatSeq[i] <- tMat[i, 1]
            tMatProb[i] <- tMat[i, 2]
        }
    }

    tMatSeq <- as.vector(tMatSeq)
    tMatProb <- as.vector(tMatProb)

    res = NULL
    if( !is.null(userProb) )
    {
        res <- invisible(.Call("_bcSeq_CRISPR_user_matching_DNAString", 
            PACKAGE = "bcSeq",
            readSeq, readSeq_ids, readPhred,libSeq,libSeq_ids, outFile,
            misMatch, tMatSeq, tMatProb, numThread, count_only, gap_left,
            ext_left, gap_right, ext_right, pen_max, userProb))
    } else
    {
        res <- invisible(.Call("_bcSeq_CRISPR_matching_DNAString", 
            PACKAGE = "bcSeq",readSeq, readSeq_ids, readPhred,libSeq,
            libSeq_ids, outFile, misMatch, tMatSeq,
            tMatProb, numThread, FALSE, count_only, gap_left, ext_left,
            gap_right, ext_right, pen_max))
    }

    if (!count_only)
    {
        out <- do.call(Matrix::sparseMatrix, res[[2]])
        c(res[1], out)
    }
}


bcSeq_hamming <- function(sampleFile, libFile, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = TRUE)
{
    if(class(sampleFile) =="character")
    {
        .bcSeq_hamming(sampleFile=sampleFile, libFile=libFile,
            outFile =outFile, misMatch = misMatch,
            tMat = tMat, numThread = numThread, count_only = count_only)
    } else if(class(sampleFile) == "DNAStringSet")
    {
        .bcSeq_hamming_DNAString(sampleFile=sampleFile, libFile=libFile,
            outFile =outFile, misMatch = misMatch,
            tMat = tMat, numThread = numThread, count_only = count_only)
    }
}

bcSeq_edit <- function(sampleFile, libFile, outFile, misMatch = 2,
    tMat = NULL, numThread = 4, count_only = TRUE, gap_left = 3,
    ext_left = 1, gap_right = 3, ext_right = 1,
    pen_max = 6, userProb = NULL)
{
    if(class(sampleFile) == "character")
    {
        .bcSeq_edit(sampleFile, libFile, outFile, misMatch,
            tMat, numThread, count_only, gap_left,
            ext_left, gap_right, ext_right,
            pen_max, userProb)

    } else if(class(sampleFile) == "DNAStringSet")
    {
        .bcSeq_edit_DNAString(sampleFile, libFile, outFile, misMatch,
            tMat, numThread, count_only, gap_left,
            ext_left, gap_right, ext_right,
            pen_max, userProb)
    }
}
