library(tidyverse)
library(dplyr)

#########################################################################################################################
# Define some functions
#########################################################################################################################

# A function to combine rows in df that represent the same sample (sample name = Passage) 
coalesce_by_column <- function(df) {
  return(dplyr::coalesce(!!! as.list(df)))
}

# A function to compute the symmetric difference (opposite of intersection) between two vectors 
sym_diff <- function(a,b) unique(c(setdiff(a,b), setdiff(b,a)))

#########################################################################################################################
# Define some parameters for this analysis 
#########################################################################################################################

# List of mutations of interest
lineage_specific_mutations <- c("polymerase_PA M579I",	
                                "polymerase_PA N222S",	
                                "polymerase_PA S395N",	
                                "polymerase_PB1 M290V",
                                "polymerase_PB1 M339I",	
                                "polymerase_PB1 T46A",	
                                "polymerase_PB1 V285I",	
                                "polymerase_PB2 E180K",	
                                "polymerase_PB2 E191K",	
                                "polymerase_PB2 K189R",	
                                "polymerase_PB2 T491M",	
                                "polymerase_PB2 Y488C")

# Minimum depth for considering allele frequency at a position 
min_depth <- 20

#########################################################################################################################
setwd('/Volumes/lizso_backup_drive/Lieber_Flu_final/Figure-5H_Ferret-Transmission/')
hyb_all <- read_csv('rava_ferret_final_tables_by_lineage/lineage4_and_lineage5_final.csv')

# remove spaces in column names
names(hyb_all)<-str_replace_all(names(hyb_all), c(" " = "." , "," = "" ))

# get list of all samples
all_samples <- unique(hyb_all$Passage)

hyb_lin_specific <- hyb_all %>% filter(Amino.Acid.Change %in% lineage_specific_mutations)

# list all samples that don't have lineage-specific mutations
lin_mut_samples <- unique(hyb_lin_specific$Passage)
no_lin_muts <- sym_diff(all_samples, lin_mut_samples)

# remove "polymerase_" and the spaces from Amino.Acid.Change values for legibility
hyb_lin_specific$Amino.Acid.Change <- str_replace_all(hyb_lin_specific$Amino.Acid.Change, 
                                                      c("polymerase_" = "", " " = ":"))

# remove unused columns from the dataframe
temp <- hyb_lin_specific %>% select(-one_of('Position', 
                                            'Change', 
                                            'Protein', 
                                            'NucleotideChange', 
                                            'LetterChange',
                                            'Syn',
                                            'Sample')) 

# a slightly hacky way to add samples without lineage mutationss to the final table for downstream analysis. Chose 
# one AA change at random to initialize the rows (PB2:Y488C). The rest will get filled in with after pivot_wider(). 
temp <- temp %>% add_row(Passage = no_lin_muts, Amino.Acid.Change = "PB2:Y488C")

# calculate the number of reads with the mutant allele
temp$mut_allele_reads <- temp$AF * temp$Depth / 100

# Sum the depth and number of mutant allele reads across replicate entries (sum the data from first and second sequencing runs)
temp <- temp %>% group_by(Passage, Amino.Acid.Change) %>% summarise(across(c(mut_allele_reads, Depth), sum))

# Re-calculate the % AF using the summed Depth and number of mutant allele reads
temp$AF <- temp$mut_allele_reads / temp$Depth * 100
temp <- temp %>% select(-one_of('mut_allele_reads'))

# Find samples with low depth at loci of interest
isd <- temp %>% filter(Depth <= min_depth)
isd <- unique(isd$Passage)

# add columns for each lineage-specific mutation (Rows are filled in with null if data not already there)
final_hyb_lin_specific <- temp %>% 
  pivot_wider(names_from = Amino.Acid.Change,
              values_from = c(AF,Depth))

# Reorder columns so that AF and Depth for each mutation are next to each other and mutations are grouped by segment
# TODO: would rather this not be hardcoded
final_hyb_lin_specific <- final_hyb_lin_specific[, c(1,2,7,4,9,3,8,5,10,6,11)]

# Filter out the samples with insufficient depth, put them in a separate df
insufficient_depth <- final_hyb_lin_specific %>% filter(Passage %in% isd)
final_hyb_lin_specific <- final_hyb_lin_specific %>% filter(!(Passage %in% isd))


write.csv(final_hyb_lin_specific, file = "ferret_RAVA_output_to_figure_input/ferret_summary_lineage_mutations.csv")

write.csv(no_lin_muts, file = "ferret_RAVA_output_to_figure_input/samples_without_lineage_mutations.csv")

write.csv(insufficient_depth, file = "ferret_RAVA_output_to_figure_input/ferret_insufficient_depth.csv")