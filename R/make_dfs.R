# Make dataframe to combine with blast outputs

# Get the first value that is not NA
shiftNAs <- function(inmat, outList = TRUE, fill = NA, origDim = FALSE) {
        myList <- lapply(sequence(nrow(inmat)), function(x) {
                y <- inmat[x, ]
                y[!is.na(y)]})
        Len <- vapply(myList, length, 1L)
        if (isTRUE(origDim)) Ncol <- ncol(inmat) else Ncol <- max(Len)
        Nrow <- nrow(inmat)
        M <- matrix(fill, ncol = Ncol, nrow = Nrow)
        M[cbind(rep(sequence(Nrow), Len), sequence(Len))] <- 
                unlist(myList, use.names=FALSE)
        M
}

# Use the shiftNAs function above to format the dataframe for combining with blast outputs
make_model_df <- function(x){
        maf <- data.frame(cbind(x$MAF.x, x$MAF.y, x$MAF.x.x, x$MAF.y.y)) %>%
                shiftNAs()
        depth <- data.frame(cbind(x$AlignDepth.x, x$AlignDepth.y, x$AlignDepth.x.x, x$AlignDepth.y.y)) %>%
                shiftNAs()
        df <- x[c(1:7)]
        df$MAF <- maf[,1] # use first value
        df$depth <- depth[,1] # use first value
        df <- mutate(df, Query = paste(Contig_Name, SNP_Position, sep = '_'))
        df <- df
}