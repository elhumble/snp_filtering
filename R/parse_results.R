# Functions to extract main results included in Table 1 & Table 2

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Get overall outcomes, broken down by share
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Determining outcomes by share number

get_overall_outcomes <- function(infinium_df, axiom_df, 
                                 infinium_model_pass, infinium_filter_pass,
                                 axiom_model_pass, axiom_filter_pass){
        
        # Make empty matrix
        share_outcomes <- matrix(ncol = 5,
                         nrow = 16)

        # Fill matrix
        share_outcomes[,1] <- c(rep(c(1,2,3,4), 4))
        share_outcomes[,2] <- c(rep("inf", 8), rep("ax", 8))
        share_outcomes[,3] <- c(rep("mod", 4), rep("filt", 4))
        share_outcomes <- data.frame(share_outcomes)
        names(share_outcomes) <- c("share", "tech","method", 
                                   "initial_share", "success_share")
        head(share_outcomes)

        # Number of SNPs in final infinium and axiom datasets prior to predictive modelling

        share_outcomes$initial_share <- c(nrow(subset(infinium_df, infinium_df$share == 1)),
                                  nrow(subset(infinium_df, infinium_df$share == 2)),
                                  nrow(subset(infinium_df, infinium_df$share == 3)),
                                  nrow(subset(infinium_df, infinium_df$share == 4)),
                                  nrow(subset(infinium_df, infinium_df$share == 1)),
                                  nrow(subset(infinium_df, infinium_df$share == 2)),
                                  nrow(subset(infinium_df, infinium_df$share == 3)),
                                  nrow(subset(infinium_df, infinium_df$share == 4)),
                                  nrow(subset(axiom_df, axiom_df$share == 1)),
                                  nrow(subset(axiom_df, axiom_df$share == 2)),
                                  nrow(subset(axiom_df, axiom_df$share == 3)),
                                  nrow(subset(axiom_df, axiom_df$share == 4)),
                                  nrow(subset(axiom_df, axiom_df$share == 1)),
                                  nrow(subset(axiom_df, axiom_df$share == 2)),
                                  nrow(subset(axiom_df, axiom_df$share == 3)),
                                  nrow(subset(axiom_df, axiom_df$share == 4)))

        # Number of SNPs in each dataset passing predictive modelling

        share_outcomes$success_share <- c(nrow(subset(infinium_model_pass, infinium_model_pass$share == 1)),
                                  nrow(subset(infinium_model_pass, infinium_model_pass$share == 2)),
                                  nrow(subset(infinium_model_pass, infinium_model_pass$share == 3)),
                                  nrow(subset(infinium_model_pass, infinium_model_pass$share == 4)),
                                  nrow(subset(infinium_filter_pass, infinium_filter_pass$share == 1)),
                                  nrow(subset(infinium_filter_pass, infinium_filter_pass$share == 2)),
                                  nrow(subset(infinium_filter_pass, infinium_filter_pass$share == 3)),
                                  nrow(subset(infinium_filter_pass, infinium_filter_pass$share == 4)),
                                  nrow(subset(axiom_model_pass, axiom_model_pass$share == 1)),
                                  nrow(subset(axiom_model_pass, axiom_model_pass$share == 2)),
                                  nrow(subset(axiom_model_pass, axiom_model_pass$share == 3)),
                                  nrow(subset(axiom_model_pass, axiom_model_pass$share == 4)),
                                  nrow(subset(axiom_filter_pass, axiom_filter_pass$share == 1)),
                                  nrow(subset(axiom_filter_pass, axiom_filter_pass$share == 2)),
                                  nrow(subset(axiom_filter_pass, axiom_filter_pass$share == 3)),
                                  nrow(subset(axiom_filter_pass, axiom_filter_pass$share == 4)))

# Success rate across SNPs share sets

share_outcomes <- share_outcomes %>%
        mutate(success_rate = success_share / initial_share * 100)

share_outcomes <- share_outcomes

}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Determine how callers are represented in the final SNP datasets
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

get_caller_outcomes <- function(infinium_model_pass, infinium_filter_pass, 
                                axiom_model_pass, axiom_filter_pass){
        
# Make empty results dataframe
caller_outcomes <- matrix(ncol = 9,
                          nrow = 16)

# Fill matrix
caller_outcomes[,1] <- c(rep(c(1,2,3,4), 4))
caller_outcomes[,3] <- c(rep("inf", 8), rep("ax", 8))
caller_outcomes[,4] <- c(rep("mod", 4), rep("filt", 4))
caller_outcomes[,5] <- share_outcomes$success_rate
caller_outcomes <- data.frame(caller_outcomes)
names(caller_outcomes) <- c("share", "total", "tech","method", "success_rate", "bowtie", "bwa", "newbler", "swap454")

# Get total number of putatively succesfull SNPS in final datasets
# Broken down by share, technology, filtering approach (modelling vs simple filtering)

caller_outcomes$total <- c(nrow(subset(infinium_model_pass, infinium_model_pass$share == 1)),
                                  nrow(subset(infinium_model_pass, infinium_model_pass$share == 2)),
                                  nrow(subset(infinium_model_pass, infinium_model_pass$share == 3)),
                                  nrow(subset(infinium_model_pass, infinium_model_pass$share == 4)),
                                  nrow(subset(infinium_filter_pass, infinium_filter_pass$share == 1)),
                                  nrow(subset(infinium_filter_pass, infinium_filter_pass$share == 2)),
                                  nrow(subset(infinium_filter_pass, infinium_filter_pass$share == 3)),
                                  nrow(subset(infinium_filter_pass, infinium_filter_pass$share == 4)),
                                  nrow(subset(axiom_model_pass, axiom_model_pass$share == 1)),
                                  nrow(subset(axiom_model_pass, axiom_model_pass$share == 2)),
                                  nrow(subset(axiom_model_pass, axiom_model_pass$share == 3)),
                                  nrow(subset(axiom_model_pass, axiom_model_pass$share == 4)),
                                  nrow(subset(axiom_filter_pass, axiom_filter_pass$share == 1)),
                                  nrow(subset(axiom_filter_pass, axiom_filter_pass$share == 2)),
                                  nrow(subset(axiom_filter_pass, axiom_filter_pass$share == 3)),
                                  nrow(subset(axiom_filter_pass, axiom_filter_pass$share == 4)))

# Extract the subsets as variables

# subset final infinium datasets
i_one_m <- subset(infinium_model_pass, infinium_model_pass$share == 1)
i_two_m <- subset(infinium_model_pass, infinium_model_pass$share == 2)
i_three_m <- subset(infinium_model_pass, infinium_model_pass$share == 3)
i_four_m <- subset(infinium_model_pass, infinium_model_pass$share == 4)

i_one_f <- subset(infinium_filter_pass, infinium_filter_pass$share == 1)
i_two_f <- subset(infinium_filter_pass, infinium_filter_pass$share == 2)
i_three_f <- subset(infinium_filter_pass, infinium_filter_pass$share == 3)
i_four_f <- subset(infinium_filter_pass, infinium_filter_pass$share == 4)

# subset final axiom datasets

a_one_m <- subset(axiom_model_pass, axiom_model_pass$share == 1)
a_two_m <- subset(axiom_model_pass, axiom_model_pass$share == 2)
a_three_m <- subset(axiom_model_pass, axiom_model_pass$share == 3)
a_four_m <- subset(axiom_model_pass, axiom_model_pass$share == 4)

a_one_f <- subset(axiom_filter_pass, axiom_filter_pass$share == 1)
a_two_f <- subset(axiom_filter_pass, axiom_filter_pass$share == 2)
a_three_f <- subset(axiom_filter_pass, axiom_filter_pass$share == 3)
a_four_f <- subset(axiom_filter_pass, axiom_filter_pass$share == 4)

# Function to determine the SNPS in the final dataset that were called by
# BOWTIE2, BWA, NEWBLER & SWAP454 respectively. Broken down by share.

# Must specify exact subsets
get_callers <- function(callerID, one, two, three, four){
        out <- c(length(which(one$Query %in% callerID)),
                 length(which(two$Query %in% callerID)),
                 length(which(three$Query %in% callerID)),
                 length(which(four$Query %in% callerID)))
        out <- out
}

# Modify caller IDs
newblerID <- gsub(" ", "_", newblerID)
swap454ID <- gsub(" ", "_", swap454ID)
bwaID <- gsub(" ", "_", bwaID)
bowtieID <- gsub(" ", "_", bowtieID)

# Get values and enter into dataframe
caller_outcomes$bowtie <- c(get_callers(bowtieID, 
                                                i_one_m, i_two_m, i_three_m, i_four_m),
                            get_callers(bowtieID, 
                                                i_one_f, i_two_f, i_three_f, i_four_f),
                            get_callers(bowtieID, 
                                                a_one_m, a_two_m, a_three_m, a_four_m),
                            get_callers(bowtieID, 
                                                a_one_f, a_two_f, a_three_f, a_four_f))

caller_outcomes$bwa <- c(get_callers(bwaID, 
                                             i_one_m, i_two_m, i_three_m, i_four_m), 
                         get_callers(bwaID, 
                                             i_one_f, i_two_f, i_three_f, i_four_f),
                         get_callers(bwaID, 
                                             a_one_m, a_two_m, a_three_m, a_four_m),
                         get_callers(bwaID, 
                                             a_one_f, a_two_f, a_three_f, a_four_f))

caller_outcomes$newbler <- c(get_callers(newblerID, 
                                                 i_one_m, i_two_m, i_three_m, i_four_m), 
                             get_callers(newblerID, 
                                                 i_one_f, i_two_f, i_three_f, i_four_f),
                             get_callers(newblerID, 
                                                 a_one_m, a_two_m, a_three_m, a_four_m),
                             get_callers(newblerID, 
                                                 a_one_f, a_two_f, a_three_f, a_four_f))

caller_outcomes$swap454 <- c(get_callers(swap454ID, 
                                                 i_one_m, i_two_m, i_three_m, i_four_m), 
                             get_callers(swap454ID, 
                                                 i_one_f, i_two_f, i_three_f, i_four_f),
                             get_callers(swap454ID, 
                                                 a_one_m, a_two_m, a_three_m, a_four_m),
                             get_callers(swap454ID, 
                                                 a_one_f, a_two_f, a_three_f, a_four_f))

caller_outcomes <- caller_outcomes

}


