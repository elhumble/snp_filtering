# Generate files for illumina and affymetrix to calculate probe design scores

# For some swap454 SNPs, the corresponding base pair in the genome does not match 
# either the reference or alternative allele identified by the SNP caller. These 
# functions also include a step to remove these.

getADTfile <- function(x){
  
  require(dplyr)
  #x <- all.snps_infinium[c(1,2,7,8,6,36)]
  x$seq <- toupper(x$seq) # make upper case for string matching
  infinium_snps <- data.frame(Ref = substr(x$seq, x$SNPenclose, x$SNPenclose)) # get ref allele according to allele in seq
  
  # get alt alleles according to allele in seq
  for (i in 1:length(x$seq)){
    if (infinium_snps$Ref[i] == x$SNP_RefAllele[i]){
      infinium_snps$Alt[i] <- as.character(x$SNP_AltAllele[i])}
    else {
      infinium_snps$Alt[i] <- as.character(x$SNP_RefAllele[i])
    }
  }
  
  # alphabetical
  infinium_snps <- t(apply(infinium_snps, 1, sort))
  x <- data.frame(x, infinium_snps)
  
  # check errors
  check_snp <- function(df) {
    good <- TRUE
    if ((df["X1"] != df["SNP_RefAllele"]) & (df["X1"] != df["SNP_AltAllele"])) good <- FALSE
    if ((df["X2"] != df["SNP_RefAllele"]) & (df["X2"] != df["SNP_AltAllele"])) good <- FALSE
    good
  }
  
  snp_goodness <- apply(x, 1, check_snp)
  x$goodness <- snp_goodness
  x <- filter(x, goodness == TRUE) # filter where seq does not agree with caller, 30149
  
  for (i in 1:length(x$seq)){
    substr(x$seq[i], x$SNPenclose[i], x$SNPenclose[i]) <- "]"
  }
  
  to <- paste("[", x$X1, "/", x$X2, "]", sep = "")
  
  for(i in 1:length(x$seq))
    x$seq[i]<-gsub("]", to[i], x$seq[i])
  
  x <- mutate(x, SNPid = paste(Contig_Name, SNP_Position, sep = '_'))
  
  adt_df <- data.frame(Locus_Name = x$SNPid,
                       Target_Type = "SNP",
                       Sequence = x$seq,
                       Chromosome = 0,
                       Coordinate = 0,
                       Genome_Build_Version = 0,
                       Source = "UNKNOWN",
                       Source_Version = 0,
                       Sequence_Orientation = "UNKNOWN",
                       Plus_Minus = "UNKNOWN",
                       Force_Infinium_I = "FALSE")
  
  # adt_df <- adt_df
  write.csv(adt_df, "data/processed/WGGT_Seq_FurSeal.csv", quote = F, row.names = F)
  
}

getPconvertfile <- function(x){
  
  require(dplyr)
  x$seq <- toupper(x$seq) # make upper case for string matching
  axiom_snps <- data.frame(Ref = substr(x$seq, x$SNPenclose, x$SNPenclose))
  
  # get alt alleles according to allele in seq
  for (i in 1:length(x$seq)){
    if (axiom_snps$Ref[i] == x$SNP_RefAllele[i]){
      axiom_snps$Alt[i] <- as.character(x$SNP_AltAllele[i])}
    else {
      axiom_snps$Alt[i] <- as.character(x$SNP_RefAllele[i])
    }
  }
  
  # alphabetical
  axiom_snps <- t(apply(axiom_snps, 1, sort))
  x <- data.frame(x, axiom_snps)
  
  # check errors
  check_snp <- function(df) {
    good <- TRUE
    if ((df["X1"] != df["SNP_RefAllele"]) & (df["X1"] != df["SNP_AltAllele"])) good <- FALSE
    if ((df["X2"] != df["SNP_RefAllele"]) & (df["X2"] != df["SNP_AltAllele"])) good <- FALSE
    good
  }
  
  snp_goodness <- apply(x, 1, check_snp)
  x$goodness <- snp_goodness
  x <- filter(x, goodness == TRUE) # filter where seq does not agree with caller, 31580
  
  for (i in 1:length(x$seq)){
  substr(x$seq[i], x$SNPenclose[i], x$SNPenclose[i]) <- "]"
  }
          
  to <- paste("[", x$X1, "/", x$X2, "]", sep = "")
          
  for(i in 1:length(x$seq))
  x$seq[i]<-gsub("]", to[i], x$seq[i])
  
  x <- mutate(x, SNPid = paste(Contig_Name, SNP_Position, sep = '_'))
  
  axiom_df <- data.frame(Organism = "Arctocephalus_gazella",
                         Locus_Name = x$SNPid,
                         SEQ = x$seq,
                         SNP_PRIORITY = 0,
                         CHR = "unknown",
                         CHR_TYPE = "autosomal",
                         SNP_VAL = 0)
  
  # adt_df <- adt_df
  write.table(axiom_df, "data/processed/axiomdesign_Seq_FurSeal.txt",
              quote = F, col.names = T, row.names = F, sep = "\t")
  
}