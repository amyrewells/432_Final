library(rentrez)
library(BiocManager)
library(Biostrings)
library(msa)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)

# ---- Step 1: Load MetaData and OTUtable ----
MetaData <- read.delim("data/mapping.txt", row.names = "X.SampleID")
OTUtable <- read.csv("data/otu_table.csv", row.names = "OTU_ID")
write.csv(top_17, "Top17_OTU.csv")
Importance<- read.csv("Top17_OTU.csv")

# ---- Step 2: Fetch Sequences from NCBI ----
accNum <- c("KX459721", "KX460758", "KX459915", "KX460102", "KX459847", 
            "KX460465", "KX459737", "KX460199", "KX459813")
fetched <- lapply(accNum, function(accNum) {
  entrez_fetch(db = "nucleotide", id = accNum, rettype = "fasta", retmode = "text")
})
formSeq <- unlist(fetched)
writeLines(formSeq, "ncbi_sequences.fasta")

# ---- Step 3: Align Sequences ----
seq <- readDNAStringSet("ncbi_sequences.fasta")
seq <- seq[!duplicated(names(seq))]  # Remove duplicated sequences
align <- msa(seq, method = "Muscle", verbose = FALSE)
writeXStringSet(DNAStringSet(as.character(align)), "aligned_sequences.fasta")

# ---- Step 4: Build the Phylogenetic Tree ----
alignSeq <- read.dna("aligned_sequences.fasta", format = "fasta")
dist <- dist.dna(alignSeq)
phylo <- nj(dist)

# ---- Step 5: Extract Bacteria Names and OTU Numbers ----
label_info <- data.frame(
  original_label = phylo$tip.label,
  Accession = sub(" .*", "", phylo$tip.label),                        # Extract Accession numbers
  Bacteria = sub("^.*?\\d\\.\\d (.*?) clone.*", "\\1", phylo$tip.label),  # Extract Bacteria name
  OTU = sub(".*(OTU#\\d+).*", "\\1", phylo$tip.label)                 # Extract OTU numbers
)

# ---- Step 6: Clean Bacteria Names (Remove "Uncultured") ----
label_info$Bacteria <- gsub("^Uncultured ", "", label_info$Bacteria)  # Remove "Uncultured" from names

# ---- Step 7: Remove Duplicates ----
# Keeping the most important occurrence if there are duplicates
label_info <- label_info %>%
  group_by(Bacteria) %>%
  mutate(Bacteria = ifelse(duplicated(Bacteria), paste0(Bacteria, "_", row_number()), Bacteria)) %>%
  ungroup()

# ---- Step 8: Match OTUs with Cleaned Importance Data ----
otu_clean<- Importance[-c(1, 5, 7, 8, 10, 13, 14, 15), -c(2, 3, 5)]
colnames(otu_clean)[1] <- "OTU"
colnames(otu_clean)[2] <- "Importance"
otu_clean$OTU <- sub("X\\.OTU\\.ID\\.(\\d+)", "OTU#\\1", otu_clean$OTU)

# ---- Step 9: Merge Importance Data ----
matched_otu_importance <- merge(label_info, otu_clean, by = "OTU", all.x = TRUE)

# ---- Step 10: Prepare Phylogenetic Data ----
phylo_data <- fortify(phylo) %>%
  filter(isTip == TRUE) %>%  # Only keep tips for plotting
  left_join(matched_otu_importance, by = c("label" = "original_label"))

# ---- Step 11: Plotting The Tree ----
p <- ggtree(phylo) %<+% phylo_data +  
  geom_tippoint(aes(color = Importance), size = 3) +  
  geom_tiplab(aes(label = Bacteria), size = 2.5, align = TRUE, hjust = -0.1) +  
  scale_color_gradient(low = "lightblue", high = "red", na.value = "grey") +
  theme_tree() +
  ggtitle("Phylogenetic Tree")

# Display the plot
print(p)