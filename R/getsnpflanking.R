# get flanking seq positions for stuff like bedtools 
# x = snp data, y = contig lengths                   


get_illuminaflanks <- function(x, y){
  snp.pos <- x[1:2]
  snp.pos$one <- snp.pos$SNP_Position - 61
  snp.pos$two <- snp.pos$SNP_Position + 60
  snp.pos <- merge(snp.pos, y, by = "Contig_Name", sort = F)
  
  for (i in 1:length(snp.pos$two)){
    if (snp.pos$two[i] > snp.pos$Length[i]){
      snp.pos$two[i] <- snp.pos$Length[i]}
    if (snp.pos$one[i] < 1){
      snp.pos$one[i] <- 1
    }
  }
  
   snp.pos$SNPenclose <- (snp.pos$SNP_Position - snp.pos$one)
   snp.pos <- snp.pos
}


get_axiomflanks <- function(x, y){
  snp.pos <- x[1:2]
  snp.pos$one <- snp.pos$SNP_Position - 36   
  snp.pos$two <- snp.pos$SNP_Position + 35
  snp.pos <- merge(snp.pos, y, by = "Contig_Name", sort = F)
  
  for (i in 1:length(snp.pos$two)){
    if (snp.pos$two[i] > snp.pos$Length[i]){
      snp.pos$two[i] <- snp.pos$Length[i]}
    if (snp.pos$one[i] < 1){
      snp.pos$one[i] <- 1
    }
  }
  
  snp.pos$SNPenclose <- (snp.pos$SNP_Position - snp.pos$one)
  snp.pos <- snp.pos # [-5]
}