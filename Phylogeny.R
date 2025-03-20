#Remember to make into relative path
setwd("C:/statistics/R/repositories/BIOL432_Final/")

OTUtable <- read.delim("otu_table.txt", header = TRUE, sep = "\t", 
                       check.names = FALSE, row.names= "#OTU ID")
MetaData<- read.delim("mapping.txt")

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
phylo$tip.label <- sapply(strsplit(phylo$tip.label, " "), function(x) paste(x[2], x[3]))

ggtree(phylo) +
  geom_tiplab() +  # Add tip labels to show the sample names
  theme_tree() +   # Adjust the tree theme
  ggtitle("Phylogenetic Tree (Neighbor Joining)")

#figure out how to colour by species? How to link accession numbers to metadata...