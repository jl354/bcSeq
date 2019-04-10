.bcSeq_hamming <- function(sampleFile, libFile, outFile, misMatch = 2, 
    tMat = NULL, numThread = 4, count_only = TRUE, detail_info = FALSE) {
    # check file status
    if (!file.exists(sampleFile)) 
        stop(paste0(sampleFile, " does not exist!"))
    if (!file.exists(libFile)) 
        stop(paste0(libFile, " does not exist!"))
    if (file.exists(outFile)) 
        stop(paste0(outFile, " exists plese specify
        another name."))
    
    tMatSeq <- c("default")
    tMatProb <- c(0.3333)
    if (is.data.frame(tMat)) {
        for (i in 1:nrow(tMat)) {
            tMatSeq[i] <- tMat[i, 1]
            tMatProb[i] <- tMat[i, 2]
        }
    }
    
    tMatSeq <- as.vector(tMatSeq)
    tMatProb <- as.vector(tMatProb)
    res <- invisible(.Call("_bcSeq_CRISPR_matching", PACKAGE = "bcSeq", 
        sampleFile, libFile, outFile, misMatch, tMatSeq, tMatProb, 
        numThread, TRUE, count_only, 0, 0, 0, 0, 0, detail_info))
    
    if (!count_only) {
        out <- Matrix::sparseMatrix(i=(res[[2]]$i+1), j =(res[[2]]$j+1), 
            x=res[[2]]$x, dims = c(as.integer(max(res[[2]]$i)+1),
            as.integer(max(res[[2]]$j)+1)))
        #out <- do.call(Matrix::sparseMatrix, res[[2]])
        c(res[1], alignProb = out)
    }
    
}

.bcSeq_edit <- function(sampleFile, libFile, outFile, misMatch = 2, 
    tMat = NULL, numThread = 4, count_only = TRUE, gap_left = 3, 
    ext_left = 1, gap_right = 3, ext_right = 1, pen_max = 6, 
    userProb = NULL, detail_info = FALSE) {
    # check file status
    if (!file.exists(sampleFile)) 
        stop(paste0(sampleFile, " does not exist!"))
    if (!file.exists(libFile)) 
        stop(paste0(libFile, " does not exist!"))
    if (file.exists(outFile)) 
        stop(paste0(outFile, " exists plese specify another name."))
    
    tMatSeq <- c("default")
    tMatProb <- c(0.3333)
    if (is.data.frame(tMat)) {
        for (i in 1:nrow(tMat)) {
            tMatSeq[i] <- tMat[i, 1]
            tMatProb[i] <- tMat[i, 2]
        }
    }
    
    tMatSeq <- as.vector(tMatSeq)
    tMatProb <- as.vector(tMatProb)
    
    res <- NULL
    if (is.character(userProb)) {
        res <- invisible(.Call("_bcSeq_CRISPR_user_matching", 
            PACKAGE = "bcSeq", sampleFile, libFile, outFile, 
            misMatch, tMatSeq, tMatProb, numThread, count_only, 
            gap_left, ext_left, gap_right, ext_right, pen_max, 
            userProb))
    } else {
        res <- invisible(.Call("_bcSeq_CRISPR_matching", PACKAGE = "bcSeq", 
            sampleFile, libFile, outFile, misMatch, tMatSeq, 
            tMatProb, numThread, FALSE, count_only, gap_left, 
            ext_left, gap_right, ext_right, pen_max, detail_info))
    }
    
    if (!count_only) {
        out <- Matrix::sparseMatrix(i=(res[[2]]$i+1), j =(res[[2]]$j+1), 
            x=res[[2]]$x, dims = c(as.integer(max(res[[2]]$i)+1),
            as.integer(max(res[[2]]$j)+1)))
        #out <- do.call(Matrix::sparseMatrix, res[[2]])
        c(res[1], out)
    }
}

.bcSeq_hamming_DNAString <- function(sampleFile, libFile, outFile, 
    misMatch = 2, tMat = NULL, numThread = 4, count_only = TRUE, 
    detail_info = FALSE) {
    
    readSeq <- as.vector(as.character(sampleFile))
    readSeq_ids <- as.vector(names(sampleFile))
    readPhred <- as.vector(as.character(unlist(sampleFile@metadata)))
    libSeq <- as.vector(as.character(libFile))
    libSeq_ids <- as.vector(names(libFile))
    
    tMatSeq <- c("default")
    tMatProb <- c(0.3333)
    if (is.data.frame(tMat)) {
        for (i in 1:nrow(tMat)) {
            tMatSeq[i] <- tMat[i, 1]
            tMatProb[i] <- tMat[i, 2]
        }
    }
    
    tMatSeq <- as.vector(tMatSeq)
    tMatProb <- as.vector(tMatProb)
    res <- invisible(.Call("_bcSeq_CRISPR_matching_DNAString", 
        PACKAGE = "bcSeq", readSeq, readSeq_ids, readPhred, 
        libSeq, libSeq_ids, outFile, misMatch, tMatSeq, tMatProb, 
        numThread, TRUE, count_only, 0, 0, 0, 0, 0, detail_info))
    
    if (!count_only) {
        out <- Matrix::sparseMatrix(i=(res[[2]]$i+1), j =(res[[2]]$j+1), 
            x=res[[2]]$x, dims = c(as.integer(max(res[[2]]$i)+1),
            as.integer(max(res[[2]]$j)+1)))
        #out <- do.call(Matrix::sparseMatrix, res[[2]])
        c(res[1], alignProb = out)
    }
    
}

.bcSeq_edit_DNAString <- function(sampleFile, libFile, outFile, 
    misMatch = 2, tMat = NULL, numThread = 4, count_only = TRUE, 
    gap_left = 3, ext_left = 1, gap_right = 3, ext_right = 1, 
    pen_max = 6, userProb = NULL, detail_info = FALSE) {
    
    readSeq <- as.vector(as.character(sampleFile))
    readSeq_ids <- as.vector(names(sampleFile))
    readPhred <- as.vector(as.character(unlist(sampleFile@metadata)))
    libSeq <- as.vector(as.character(libFile))
    libSeq_ids <- as.vector(names(libFile))
    
    tMatSeq <- c("default")
    tMatProb <- c(0.3333)
    if (is.data.frame(tMat)) {
        for (i in 1:nrow(tMat)) {
            tMatSeq[i] <- tMat[i, 1]
            tMatProb[i] <- tMat[i, 2]
        }
    }
    
    tMatSeq <- as.vector(tMatSeq)
    tMatProb <- as.vector(tMatProb)
    
    res <- NULL
    if (is.character(userProb)) {
        res <- invisible(.Call("_bcSeq_CRISPR_user_matching_DNAString", 
            PACKAGE = "bcSeq", readSeq, readSeq_ids, readPhred, 
            libSeq, libSeq_ids, outFile, misMatch, tMatSeq, 
            tMatProb, numThread, count_only, gap_left, ext_left, 
            gap_right, ext_right, pen_max, userProb))
    } else {
        res <- invisible(.Call("_bcSeq_CRISPR_matching_DNAString", 
            PACKAGE = "bcSeq", readSeq, readSeq_ids, readPhred, 
            libSeq, libSeq_ids, outFile, misMatch, tMatSeq, 
            tMatProb, numThread, FALSE, count_only, gap_left, 
            ext_left, gap_right, ext_right, pen_max, detail_info))
    }
    
    if (!count_only) {
        out <- Matrix::sparseMatrix(i=(res[[2]]$i+1), j =(res[[2]]$j+1), 
            x=res[[2]]$x, dims = c(as.integer(max(res[[2]]$i)+1),
            as.integer(max(res[[2]]$j)+1)))
        #out <- do.call(Matrix::sparseMatrix, res[[2]])
        c(res[1], out)
    }
}


bcSeq_hamming <- function(sampleFile, libFile, outFile, misMatch = 2, 
    tMat = NULL, numThread = 4, count_only = TRUE, detail_info = FALSE) {
    if (is.character(sampleFile)) {
        .bcSeq_hamming(sampleFile = sampleFile, libFile = libFile, 
            outFile = outFile, misMatch = misMatch, tMat = tMat, 
            numThread = numThread, count_only = count_only, 
            detail_info = detail_info)
    } else {
        .bcSeq_hamming_DNAString(sampleFile = sampleFile, libFile = libFile, 
            outFile = outFile, misMatch = misMatch, tMat = tMat, 
            numThread = numThread, count_only = count_only, 
            detail_info = detail_info)
    }
}

bcSeq_edit <- function(sampleFile, libFile, outFile, misMatch = 2, 
    tMat = NULL, numThread = 4, count_only = TRUE, gap_left = 3, 
    ext_left = 1, gap_right = 3, ext_right = 1, pen_max = 6, 
    userProb = NULL, detail_info = FALSE) {
    if (is.character(sampleFile)) {
        .bcSeq_edit(sampleFile, libFile, outFile, misMatch, 
            tMat, numThread, count_only, gap_left, ext_left, 
            gap_right, ext_right, pen_max, userProb)
    } else {
        .bcSeq_edit_DNAString(sampleFile, libFile, outFile, 
            misMatch, tMat, numThread, count_only, gap_left, 
            ext_left, gap_right, ext_right, pen_max, userProb, 
            detail_info)
    }
}
