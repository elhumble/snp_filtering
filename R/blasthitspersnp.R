# Function that takes blast output and summarises on SNP by SNP 
# basis for model prediction

# x = blast output format 6, y = snp metadata - should include columns: Query, Seq

blasthitspersnp <- function(x, y){
  
  df <- data.frame()
  
  colnames(x) <- c("Query","Subject","PercentIdentity", "AlignmentLength", 
                           "Mismatches", "Gap_Opening", "QueryStart", "QueryEnd", 
                           "SubjectStart", "SubjectEnd", "E.value", "BitScore")
  
  snps <- merge(x, y, by = "Query", sort=F)
  snps$Query<-gsub("Seq", " ", snps$Query)
  
  snps$Count <- 1
  require(reshape)
  molten.snps <- melt.data.frame(snps, id.vars = c("Query"),
                                 measure.vars = c("Count"))
  melt.snps <- cast(molten.snps, formula = Query ~ variable, sum)
  snps$Count <- NULL 
  snps$SNP.Poly <- NULL 
  
  snps_tophit <- snps[!duplicated(snps$Query),]
  df <- merge(snps_tophit, melt.snps, by = "Query", sort = F)
  
}