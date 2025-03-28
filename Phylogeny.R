library(rentrez)
library(BiocManager)
library(Biostrings)
library(msa)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(tidyverse)

# ---- Step 1: Load MetaData, OTUtable, and Importance Values ----
MetaData <- read.delim("./data/input/mapping.txt", row.names = "X.SampleID")
OTUtable <- read.csv("./data/input/otu_table.csv", row.names = "OTU_ID")
Importance<- read.csv("./data/top_10.csv")

# ---- Step 2: Fetch Sequences from NCBI ----
accNum <- c("KX459959", "KX459698", "KX460680", "KX459898", 
            "KX460175", "KX460089", "KX460557", "KX460234", "KX459830", "KX460616")

fetched <- lapply(accNum, function(accNum) {
  entrez_fetch(db = "nucleotide", id = accNum, rettype = "fasta", retmode = "text")})
formSeq <- unlist(fetched)
writeLines(formSeq, "ncbi_sequences.fasta")

# ---- Step 3: Align Sequences ----
seq <- readDNAStringSet("ncbi_sequences.fasta")
align <- msa(seq, method = "Muscle", verbose = FALSE)
writeXStringSet(DNAStringSet(as.character(align)), "aligned_sequences.fasta")

# ---- Step 4: Distance Matrix ----
alignSeq <- read.dna("aligned_sequences.fasta", format = "fasta")
dist <- dist.dna(alignSeq)
dist_matrix <- as.matrix(dist)
phylo <- nj(dist)

# ---- Step 5: Extract and Clean Label Info ----
label_info <- data.frame(
  original_label = phylo$tip.label,
  Accession = sub(" .*", "", phylo$tip.label),                        
  Bacteria = sub("^.*?\\d\\.\\d (.*?) clone.*", "\\1", phylo$tip.label),
  OTU = sub(".*(OTU#\\d+).*", "\\1", phylo$tip.label))

label_info$Bacteria <- gsub("^Uncultured ", "", label_info$Bacteria) 

label_info[9, 3] <- "Bacteria_OTU#10"
label_info[10, 3] <- "Bacteria_OTU#5678"

# ---- Step 6: Replace labels with bacteria names ----
rownames(dist_matrix) <- label_info$Bacteria
colnames(dist_matrix) <- label_info$Bacteria

# ---- Step 7: Reshape for ggplot ----
dist_long <- as.data.frame(dist_matrix, .name_repair = "minimal") %>%
  rownames_to_column("Bacteria1") %>%
  pivot_longer(-Bacteria1, names_to = "Bacteria2", values_to = "Distance")

# ---- Step 8: Visualize Distances ----
Heatmap<- ggplot(dist_long, aes(x = Bacteria1, y = Bacteria2, fill = Distance)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "black", high = "lightblue") +
  labs(title = "Phylogenetic Distance Heatmap of Microbial Taxa",
       x = "Bacteria", y = "Bacteria") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
        axis.text.y = element_text(face = "italic"),
        panel.grid = element_blank())

print(Heatmap)
#ggsave("./Figures/Heatmap.png", plot = Heatmap, width = 8, height = 6, dpi = 300)

# ---- Step 9: Match OTUs with Cleaned Importance Data ----
otu_clean<- Importance[,-c(2, 3, 5)]
colnames(otu_clean)[1] <- "OTU"
colnames(otu_clean)[2] <- "Importance"
otu_clean$OTU <- sub("X\\.OTU\\.ID\\.(\\d+)", "OTU#\\1", otu_clean$OTU)

# ---- Step 10: Merge Importance Data ----
matched_otu_importance <- merge(label_info, otu_clean, by = "OTU", all.x = TRUE)

# ---- Step 11: Prepare Phylogenetic Data ----
phylo_data <- fortify(phylo) %>%
  filter(isTip == TRUE) %>%  # Only keep tips for plotting
  left_join(matched_otu_importance, by = c("label" = "original_label"))

# ---- Step 12: Plotting The Tree ----
Tree <- ggtree(phylo) %<+% phylo_data +  
  geom_tippoint(aes(color = Importance), size = 3) +  
  geom_tiplab(aes(label = Bacteria), size = 2.8, align = TRUE, hjust = 1, vjust= -0.5, offset = 0.1) +
  scale_color_gradient(low = "red", high = "green") + 
  theme_tree() 

print(Tree)
#ggsave("./Figures/Tree.png", plot = Tree, width = 8, height = 6, dpi = 300)