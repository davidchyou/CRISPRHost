summary_csv <- function(file_in, file_out) {
	pdata <- read.csv(file_in)
	pdata <- pdata[order(-1 * pdata$N_HOST_SPACER, pdata$VIRAL_GENOME, pdata$SPACER_ASSEMBLY),]
	data_out <- data.frame(VIRAL_GENOME = as.character(pdata$VIRAL_GENOME),
                           N_HOST_SPACER = as.numeric(pdata$N_HOST_SPACER),
                           E_VALUE = as.numeric(pdata$E_VALUE),
                           BITSCORE = as.numeric(pdata$BITSCORE),
                           PROP_NT_MATCHED_IN_HYBRID = as.numeric(round(pdata$PROP_NT_MATCHED_IN_HYBRID, 4)),
                           SPACER_ID = as.character(pdata$SPACER_ID),
                           KNOWN_SPECIES = as.character(pdata$KNOWN_SPECIES),
                           TARGET_ALIGN_START = as.numeric(pdata$TARGET_ALIGN_START),
                           SEQ_TARGET_ALIGNED = as.character(pdata$SEQ_TARGET_ALIGNED),
                           SPACER_ORI_ALIGN = as.character(pdata$SPACER_ORI_ALIGN))
    write.csv(data_out, file_out, row.names=FALSE)
    return(1)
}

reorder_csv <- function(file_in, file_out) {
	pdata <- read.csv(file_in)
	pdata <- pdata[order(-1 * pdata$N_HOST_SPACER, pdata$VIRAL_GENOME, pdata$SPACER_ASSEMBLY),]
	data_out <- data.frame(VIRAL_GENOME = as.character(pdata$VIRAL_GENOME),
                           N_HOST_SPACER = as.numeric(pdata$N_HOST_SPACER),
                           E_VALUE = as.numeric(pdata$E_VALUE),
                           BITSCORE = as.numeric(pdata$BITSCORE),
                           PROP_NT_MATCHED_IN_HYBRID = as.numeric(round(pdata$PROP_NT_MATCHED_IN_HYBRID, 4)),
                           SPACER_ID = as.character(pdata$SPACER_ID),
                           KNOWN_SPECIES = as.character(pdata$KNOWN_SPECIES),
                           TARGET_ALIGN_START = as.numeric(pdata$TARGET_ALIGN_START),
                           SEQ_TARGET_ALIGNED = as.character(pdata$SEQ_TARGET_ALIGNED),
                           SPACER_ORI_ALIGN = as.character(pdata$SPACER_ORI_ALIGN),
                           SPACER_ID = as.character(pdata$SPACER_ID),
                           ARRAY_SCORE = as.numeric(pdata$ARRAY_SCORE),
                           SPACER_ALIGN_START = as.numeric(pdata$SPACER_ALIGN_START),
                           SPACER_ALIGN_END = as.numeric(pdata$SPACER_ALIGN_END),
                           TARGET_ALIGN_END = as.numeric(pdata$TARGET_ALIGN_END),
                           SPACER_LENGTH = as.numeric(pdata$SPACER_LENGTH),
                           N_SPACER_NT_MISSED_5P = as.numeric(pdata$N_SPACER_NT_MISSED_5P),
                           N_SPACER_NT_MISSED_3P = as.numeric(pdata$N_SPACER_NT_MISSED_3P),
                           BLAST_PIDENT = as.numeric(pdata$BLAST_PIDENT),
                           N_NT_IN_HYBRID = as.character(pdata$N_NT_IN_HYBRID),
                           SEQ_SPACER_ALIGNED = as.character(pdata$SEQ_SPACER_ALIGNED),
                           TARGET_FLANK_5P = as.character(pdata$TARGET_FLANK_5P),
                           TARGET_FLANK_3P = as.character(pdata$TARGET_FLANK_3P),
                           ORIGINAL_SPACER = as.character(pdata$ORIGINAL_SPACER), 
                           SPACER_ASSEMBLY = as.character(pdata$SPACER_ASSEMBLY))
    write.csv(data_out, file_out, row.names=FALSE)
    return(1)
}
          