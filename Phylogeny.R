#Remember to make into relative path
setwd("C:/statistics/R/repositories/BIOL432_Final/")

#OTUtable <- read.delim("otu_table.txt", header = TRUE, sep = "\t", 
                       #check.names = FALSE, row.names= "#OTU ID")
MetaData<- read.delim("mapping.txt", row.names = "X.SampleID")
OTUtable<- read.csv("otu_table.csv", row.names = "OTU_ID")

library(rentrez)
library(BiocManager)
library(Biostrings)

#How do relevant microbial communities of each dung beetle species 
#differ in terms of phylogenetic relationships?


#Replace accession numbers with those associated with OTUs found to 
#be relevant for species classification (top 10 bacteria from 
#random Forest $importance in question 1)
accNum <- paste0("KX459", seq(698, 718))

fetched <- lapply(accNum, function(accNum) {
  entrez_fetch(db = "nucleotide", id = accNum, rettype = "fasta", retmode = "text")})

formSeq <- unlist(fetched)

#Write initial sequences
writeLines(formSeq, "ncbi_sequences.fasta")

#Load
seq <- readDNAStringSet("ncbi_sequences.fasta")

align <- muscle::muscle(seq, quiet=T)

#Formatting
alignDF <- as.data.frame(as.matrix(align))

alignString <- apply(alignDF, 1, paste, collapse = "")

alignSeq <- DNAStringSet(alignString)

#Write aligned
writeXStringSet(alignSeq, "aligned_sequences.fasta")

library(ape)

alignSeq <- read.dna("aligned_sequences.fasta", format = "fasta")

dist <- dist.dna(alignSeq)

phylo <- nj(dist)

library(ggtree)
library(ggplot2)

#isolate as specific of an ID as possible for each bacteria
#example
#phylo$tip.label <- sapply(strsplit(phylo$tip.label, " "), function(x) paste(x[2], x[3]))

tip_labels <- phylo$tip.label

#Extract the OTU numbers using a regular expression
otu_numbers <- sub(".*OTU#(\\d+).*", "OTU#\\1", tip_labels)

#Update the tip labels to show only the OTU numbers
phylo$tip.label <- otu_numbers 

#which OTUs are important across both species.

# Assuming OTU importance scores are stored in a dataframe 'otu_importance'
otu_importance <- data.frame(OTU = c("OTU#10", "OTU#100",   "OTU#107",  "OTU#1007", "OTU#101",  "OTU#1015", "OTU#102", "OTU#1027",
                                     "OTU#1029", "OTU#103",  "OTU#1037", "OTU#104",  "OTU#1045", "OTU#105",  "OTU#106", 
                                     "OTU#107",  "OTU#1073", "OTU#1077", "OTU#1079", "OTU#108",  "OTU#1083", "OTU#1087"), 
                             importance = c(0.25, 0.75, 0.45, 0.85, 0.002, 0.75, 0.92, 0.87, 0.78, 0.34, 0.19, 0.89, 0.43, 0.78, 0.67, 0.11, 0.98, 0.86, 0.54, 0.58, 0.73, 0.45))

matched_otu_importance <- data.frame(label = phylo$tip.label) # Extract tree tips
matched_otu_importance <- merge(matched_otu_importance, otu_importance, 
                                by.x = "label", by.y = "OTU", all.x = TRUE)

#Branch thickness by importance
library(dplyr)
phylo_data <- full_join(fortify(phylo), matched_otu_importance, by = "label")

phylo_data <- phylo_data %>% mutate(branch_size = ifelse(is.na(importance), 0.5, importance * 3))

ggtree(phylo, aes(size = branch_size), data = phylo_data) + 
  scale_size_continuous(range = c(0.2, 2)) +  # Adjust thickness range
  theme_tree()

#Node colours by importance                       
ggtree(phylo) %<+% matched_otu_importance +  # Attach importance to the tree
  geom_tippoint(aes(color = importance), size = 3) +  # Color tree tips based on importance
  scale_color_gradient(low = "lightblue", high = "red") +
  theme_tree()

#og
ggtree(phylo) +
  geom_tiplab() +  # Add tip labels to show the sample names
  theme_tree() +   # Adjust the tree theme
  ggtitle("Phylogenetic Tree (Neighbor Joining)")

#Example phylogeny interpretation:
#Possible Functional Specialization: The clustering of Dysgonomonas 
#suggests a strong adaptation to the beetle gut, as this genus is often
#linked to carbohydrate breakdown and fermentation.

#The presence of Bacteroidales and Chryseobacterium in different branches 
#may indicate microbiome diversification for different metabolic functions.

#Distinct Clades Suggest Functional or Evolutionary Differences: Some bacterial 
#groups, such as Dysgonomonas, cluster together, possibly indicating 
#shared functional roles (e.g., gut symbiosis or digestion aid).

#Other clades (e.g., Bacteroidales, Butyrivibrio, Lachnospiraceae) might
#represent separate ecological roles in the beetle gut or dung environment.

#Cannot colour by species. Maybe edit branch thickness to visualize most important otus from random forest?