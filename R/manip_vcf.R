# Extract data from vcf, add some important variables and filter where MAF == 0

# see http://goo.gl/a6dpQ6 for vcf interpretation
# depth per allele by sample -- no. of ref bases, no. of alt bases. Unfiltered for quality
# depth of coverage -- no. of reads at given variant. Filtered for quality

manipvcf <- function(x){
        library(VariantAnnotation)
        library(reshape2)
        AD <- (geno(x)$AD[,]) #creates a list of AD values (depthperallelepersample) for each contig
        DP <- (geno(x)$DP) #creates a matrix of DP (coverage) values for each contig
      
        ADref <- sapply(AD, "[[", 1) #extract the reference allele
        ADalt <- sapply(AD, "[[", 2) #extract the alternate allele
        ADalt<-transform(ADalt) #df
        ADref<-transform(ADref)
      
        ADalt$contig <- rownames(ADalt)
        ADref$contig <- rownames(ADref)
      
        DP <- as.data.frame(DP)
        DP$contig <- rownames(DP)
      
        colnames(ADalt) <- c("ALT", "contig")
        colnames(ADref) <- c("REF", "contig")
        colnames(DP) <- c("Depth", "contig")
      
        AD <- merge(ADref, ADalt)
        data <- merge(AD, DP)
      
        data$contig <- gsub("_v1", "v1", data$contig)
        contig <- colsplit(data$contig, "[:_/]", c("Contig_Name", "SNP_Position", "SNP_RefAllele", "SNP_AltAllele"))
        data <- cbind(contig, data)
        data <-data[-5]
        colnames(data)[c(5:7)] <- c("SNP_Count_RefAllele", "SNP_Count_AltAllele", "AlignDepth")
      
        data$SNP_Freq_RefAllele <- data$SNP_Count_RefAllele/(data$SNP_Count_RefAllele + data$SNP_Count_AltAllele)
        data$SNP_Freq_AltAllele <- data$SNP_Count_AltAllele/(data$SNP_Count_RefAllele + data$SNP_Count_AltAllele)
      
        for(i in 1:length(data$SNP_Freq_AltAllele)) { 
          if(data$SNP_Freq_AltAllele[i] < data$SNP_Freq_RefAllele[i]){
            data$MAF[i] = data$SNP_Freq_AltAllele[i]
            } else data$MAF[i] = data$SNP_Freq_RefAllele[i]
        }
      
        data$Contig_Name <- gsub("v1", "_v1", data$Contig_Name)
        data$Contig_Name <- as.factor(data$Contig_Name)
        data$SNP_Position <- as.numeric(data$SNP_Position)
        
        #remove SNPs where no reads mapped against reference SNP exactly- altfreq = 1.00?
        data <- subset(data, data$MAF != 0)

}

