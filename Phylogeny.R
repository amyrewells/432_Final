#Remember to make into relative path
setwd("C:/statistics/R/repositories/BIOL432_Final/")

#OTUtable <- read.delim("otu_table.txt", header = TRUE, sep = "\t", 
                       #check.names = FALSE, row.names= "#OTU ID")
MetaData<- read.delim("mapping.txt", row.names = "X.SampleID")
OTUtable<- read.csv("otu_table.csv", row.names = "OTU_ID")

library(dplyr)
library(tidyr)

OTUtable_no_last_col <- OTUtable[,-ncol(OTUtable)]

# Transpose the modified OTU table
OTUtransClean <- as.data.frame(t(OTUtable_no_last_col))

nrow(MetaData)
nrow(OTUtransClean)

# Merge metadata and OTU table by columns
data_combined <- cbind(MetaData, OTUtransClean)

metaFinal <- data_combined %>%
  separate(Description, into = c("Species", "Temp", "Life_Stage"), 
           sep = " ", extra = "merge", fill = "right") %>%
  mutate(Species = ifelse(is.na(Life_Stage), Species, paste(Species, Temp)),
    Life_Stage = ifelse(is.na(Life_Stage), Temp, Life_Stage)  
  ) %>%
  select(-Temp, -BarcodeSequence, -LinkerPrimerSequence, -Life_Stage)  

metaFinal_colnames <- as.character(colnames(metaFinal))
selected_values <- c("10", "100", "1007", "101", "1015", "102", "1027", "1029", "103", "1037", 
                     "104", "1045", "105", "106", "107", "1073", "1077", "1079", "108", "1083", "1087")

species_column <- "Species"  # Replace with the actual name for the species column

# First, select the columns for the OTUs specified in selected_values
selected_otu_columns <- colnames(metaFinal)[colnames(metaFinal) %in% selected_values]

# Then, combine the species column with the selected OTU columns
selected_otus <- metaFinal[, c(species_column, selected_otu_columns)]

Trans<- as.data.frame(t(selected_otus))

nrow(Trans)
nrow(MetaData)

Trans <- head(Trans, -1)

Merged<- merge(Trans, MetaData)

transMerge<- as.data.frame(t(Merged))

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
library(tidyr)
library(dplyr)

BLAH<- as.data.frame(t(selected_otus))
# Set the first row as column names
colnames(BLAH) <- BLAH[1,]

# Remove the first row (which is now the column names row)
BLAH <- BLAH[-1,]

anyDuplicated(colnames(BLAH))
colnames(BLAH) <- make.names(colnames(BLAH), unique = TRUE)

library(tibble)

colnames(BLAH)[which(colnames(BLAH) == "OTU")] <- "OTU_existing"

BLAH <- rownames_to_column(BLAH, var = "OTU")

# Function to generate sequential accession numbers
generate_accession_numbers <- function(prefix, start, end) {
  seq_numbers <- seq(start, end)  # Generate the sequence of numbers
  paste0(prefix, seq_numbers)     # Append prefix to each number
}

# Define the accession range
prefix <- "KX"        # Common prefix
start_num <- 459698   # Starting number
end_num <- 460822     # Ending number

# Generate the accession numbers
accession_numbers <- generate_accession_numbers(prefix, start_num, end_num)

# Convert to data frame
ncbi_data <- data.frame(Accession = accession_numbers, OTU = NA) # OTU will be added later

# Save to CSV for later use
write.csv(ncbi_data, "ncbi_accessions.csv", row.names = FALSE)

# Read FASTA file
fasta_file <- "ncbi_sequences.fasta"
sequences <- readDNAStringSet(fasta_file)

# Extract sequence names (headers)
seq_names <- names(sequences)

# Convert to a data frame
fasta_metadata <- data.frame(Title = seq_names, stringsAsFactors = FALSE)

# Extract Accession Number (first word in the title)
fasta_metadata$Accession <- sub("^>([A-Z0-9.]+).*", "\\1", fasta_metadata$Title)

# Extract OTU ID (pattern "OTU#X")
fasta_metadata$OTU <- sub(".*(OTU#[0-9]+).*", "\\1", fasta_metadata$Title)

# Convert the OTUs in fasta_metadata to just numbers by removing the "otu" prefix
fasta_metadata$OTU <- sub("OTU#", "", fasta_metadata$OTU)


# Now you can merge with BLAH (assuming BLAH's OTU is in numeric format)
finalData <- merge(BLAH, fasta_metadata, by = "OTU")


# Remove 'title' and 'accession' columns before reshaping
finalData_cleaned <- finalData %>%
  select(-Title, -Accession)  # Exclude the 'title' and 'accession' columns

# Reshape the data (pivoting species columns into rows)
finalData_long <- finalData_cleaned %>%
  pivot_longer(cols = -OTU, names_to = "species", values_to = "OTU_ID")

df_long <- finalData %>%
  pivot_longer(
    cols = -c(OTU, Title, Accession),  # Keep OTU, Title, and Accession as identifiers
    names_to = "Species", 
    values_to = "Count")
# Make sure `phylo$tip.label` contains the correct identifiers
phylo$tip.label <- df_long$OTU  # Replace with the actual OTU IDs

# Create a data frame that links OTUs to species
# Assuming 'tree_data' has a column 'OTU' for OTU IDs and 'Species' for species names
OTU_to_species_mapping <- data.frame(OTU = df_long$OTU, Species = df_long$Species)

# Match the tree's tip labels to species using the OTU_to_species_mapping
species_for_tree <- OTU_to_species_mapping$Species[match(phylo$tip.label, OTU_to_species_mapping$OTU)]

# Now assign the species names to the tree's tip labels (or use them for coloring)
phylo$tip.label <- species_for_tree
# Make tip labels unique
phylo$tip.label <- make.unique(phylo$tip.label)

# Now you can set row names in tree_data
rownames(tree_data) <- phylo$tip.label

tree_data <- data.frame(Species = species_for_tree)
rownames(tree_data) <- phylo$tip.label  # Ensure that rownames match the tip labels of the tree

# Check tree structure
str(phylo)
# Look at the number of tips and internal nodes
length(phylo$tip.label)  # should match the number of species in your metadata
phylo$Nnode  # should be 2 * (number of tips) - 1
# Check the first few tip labels of the tree and your metadata
head(phylo$tip.label)  # Should match the OTUs/species
head(tree_data$Species)  # Should match the Species column in your metadata

# Check if all the species in metadata are in the tree
all(tree_data$Species %in% phylo$tip.label)  # This should return TRUE

# Check edges and tips
nrow(phylo$edge)  # should be 2 * (Nnode + number of tips) - 1

ggtree(phylo) +
  geom_tiplab(aes(color = Species), size = 3) +  # Color tip labels by species
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) +  # Customize colors
  theme_tree() +  # Tree-friendly theme
  ggtitle("Phylogenetic Tree Colored by Species")  # Add t

# Check that species are correctly assigned
head(phylo$tip.label)

# Visualize the tree with species coloring
library(ggtree)
ggtree(phylo, aes(color = tip.label)) + 
  geom_tiplab(aes(color = tip.label), size = 2) + 
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) + 
  ggtitle("Phylogenetic Tree Colored by Species")





    
# Assuming 'phylo$tip.label' contains OTU IDs, we'll create a data frame with tip labels and species
tree_tips <- data.frame(Species = phylo$tip.label, stringsAsFactors = FALSE)

tree_data <- cbind(tree_tips, df_long[, -1])  # Remove the OTU column from df_long to avoid duplication

#isolate as specific of an ID as possible for each bacteria
#example
#phylo$tip.label <- sapply(strsplit(phylo$tip.label, " "), function(x) paste(x[2], x[3]))
#fasta_metadata$OTU <- sapply(strsplit(fasta_metadata$OTU, " "), function(x) paste(x[2], x[3]))


# Ensure tree_data has the correct columns
tree_data <- tree_data %>% select(OTU, Species)

# Convert OTU column to match tip labels in the tree
tree_data$OTU <- as.character(tree_data$OTU)

# Ensure that the row names in tree_data are in the same order as phylo$tip.label
tree_data <- tree_data[match(phylo$tip.label, tree_data$OTU), ]
#################
Species <- df_long$Species
OTUs <- df_long$Count
# Now assign the correct tip labels to the phylo object
phylo$tip.label <- df_long$Species

head(phylo$tip.label)  # Check the first few tip labels of the tree
head(Species)  # Check the first few species in the metadata

SpecPhy <- Species[match(phylo$tip.label, Species)] 
tree_data <- data.frame(Species = SpecPhy)

any(duplicated(phylo$tip.label))  # Check if there are any duplicates
# Make tip labels unique by adding a suffix to duplicates
phylo$tip.label <- make.unique(phylo$tip.label)

# Now assign row names to tree_data
rownames(tree_data) <- phylo$tip.label

library(ggtree)
ggtree(phylo, aes(color = Species)) + 
  geom_tiplab(aes(color = Species), size = 2) + 
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) + 
  ggtitle("Phylogenetic Tree Colored by Species")

ggtree(phylo, layout = "circular", aes(color = Species)) + 
  geom_tiplab(size = 2, aes(angle = angle)) + 
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) + 
  ggtitle("Phylogenetic Tree Colored by Species")

# Number of tips (should match length of tip labels)
length(phylo$tip.label)

# Number of internal nodes
phylo$Nnode

# Check the number of edges
nrow(phylo$edge)

#############

length(phylo$tip.label)  # Number of tips in the phylogenetic tree
nrow(tree_data)  # Number of rows in tree_data


rownames(tree_data) <- phylo$tip.label 
head(phylo$tip.label)
Species<- tree_data$Species

Groups<- split(phylo$tip.label, Species)
ggtree(phylo,layout="circular",aes(colour=Group)) +geom_tiplab(size=2,aes(angle=angle))
# Check if the tree data is attached correctly to the phylo object
head(phylo$tip.data)  # Check the metadata attached to the tree
# Check the number of tips and nodes
length(phylo$tip.label)  # Number of tips
phylo$Nnode  # Number of internal nodes
nrow(phylo$edge)  # Number of edges in the tree (should be 2 * (Nnode + number of tips) - 1)

# Assuming tree_data$Species corresponds to the species information, and it matches the order of phylo$tip.label
Species <- tree_data$Species.1

# Ensure the order of Species matches the order of tip labels in the tree
# If necessary, reorder the Species to match the tip labels
Species <- Species[match(phylo$tip.label, tree_data$Count)]  # Assuming OTU is the same as phylo$tip.label

# Create a grouped list based on Species
Groups <- split(phylo$tip.label, Species)

# Plot the tree and color tips by species (Groups)
ggtree(phylo, layout = "circular", aes(color = Species)) +
  geom_tiplab(size = 2, aes(angle = angle)) +
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) +  # Customize colors
  ggtitle("Circular Phylogenetic Tree Colored by Species")  # Add a title
Tree <- ggtree(phylo) +
  geom_tiplab(aes(color = phylo$tip.data$Species), size = 3) +  # Color by species
  scale_color_manual(values = c("red", "blue", "green", "purple", "orange")) +  # Customize colors for species
  theme_tree() +  # Apply a tree-friendly theme
  ggtitle("Phylogenetic Tree Colored by Species")  # Add title


write.csv(finalData, "finalData.csv", row.names = FALSE)





# Example phylogeny, assuming OTUs are row names in phylo object
phylo$tip.label <- rownames(phylo$tip.label)


#figure out how to colour by species? How to link accession numbers to metadata...
