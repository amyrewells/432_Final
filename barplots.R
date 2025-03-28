# Version 2.2 code edited and moved to top

# Libraries used:
library(tidyverse)
library(stringr)
library(purrr)

# Step 1: Load your data
otu_table <- read_csv("otu_table.csv")
top10 <- read_csv("./Random Forest/top_10.csv")
mapping <- read_tsv("mapping.txt", comment = "#", 
                    col_names = c("SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Description"))

# Step 2: Extract numeric OTU IDs from the top 17 list
top10 <- top10 %>%
  mutate(OTU_ID = as.numeric(str_extract(`...1`, "\\d+$")))

top10_tax <- top10 %>%
  left_join(otu_table %>% select(1,63), by = "OTU_ID")

# Step 4: Extracting the most specific available taxonomic label
top10_tax <- top10_tax %>%
  mutate(Taxon_Label = map_chr(top10_tax[[7]], function(lineage) {
    if (is.na(lineage)) return("Unknown")
    tax_levels <- str_split(lineage, ";")[[1]]
    tax_levels <- rev(tax_levels)  # from specific to general
    cleaned <- str_replace_all(tax_levels, ".*__", "") %>% str_trim()
    specific <- cleaned[!cleaned %in% c("uncultured", "unclassified", "", NA)]
    if (length(specific) > 0) specific[1] else "Unknown"
  })) %>%
  mutate(Taxon_Label = str_replace_all(Taxon_Label, "_", " "))

# Step 3: Subset OTU table to top 17 OTUs and drop Consensus Lineage
top_otus <- otu_table %>%
  filter(OTU_ID %in% top10$OTU_ID) %>%
  select(-`Consensus Lineage`)

otu_rel_abund <- top_otus %>%
  column_to_rownames("OTU_ID") %>%  # Convert OTU_ID to row names
  mutate(across(everything(), as.numeric)) %>%  # Convert all columns to numeric
  rownames_to_column("OTU_ID") %>%  # Convert row names back to a column
  pivot_longer(-OTU_ID, names_to = "SampleID", values_to = "Abundance") %>%  # Reshape to long format
  mutate(OTU_ID = as.numeric(OTU_ID))  # Convert OTU_ID to numeric

# Step 5: Extract species info from mapping file
mapping <- mapping %>%
  mutate(Species = str_extract(Description, "E\\. intermedius|E\\. triangulatus")) %>% 
  filter(!str_detect(Species, "untreated dung"))

# Step 6: Merge abundance data with species metadata
otu_with_species <- otu_rel_abund %>%
  left_join(mapping %>% select(SampleID, Species), by = "SampleID")

# Step 7: Calculate mean abundance for each OTU by species
sum_abundance <- otu_with_species %>%
  group_by(OTU_ID, Species) %>%
  summarise(Abundances = sum(Abundance), .groups = "drop")

abundanceData <- sum_abundance %>%
  filter(!is.na(sum_abundance[[2]]))

top10_tax_long <- top10_tax %>%
  pivot_longer(cols = starts_with("Taxon_Label"),  # Assuming Taxon_Label columns start with 'Taxon_Label'
               names_to = "Taxon_Type", 
               values_to = "Taxon_Label") %>%
  select(OTU_ID, Taxon_Label)  # Keep only the relevant columns

final <- abundanceData %>%
  left_join(top10_tax_long %>% select(OTU_ID, Taxon_Label), by = "OTU_ID")

# Step 8: Plot the grouped bar chart
ggplot(final, aes(x = Taxon_Label, y = Abundances, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "OTU ID",
    y = "Total Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  scale_fill_manual(values = c("#66c2a5", "#fc8d62")) +
  theme_minimal(base_size = 14)
  
####################################
#old code
# Step 1: Load your data
otu_table <- read_csv("otu_table.csv")
top10 <- read_csv("./Random Forest/top_10.csv")
mapping <- read_tsv("mapping.txt", comment = "#", 
                    col_names = c("SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Description"))

# Step 2: Extract numeric OTU IDs from top 17
top10 <- top10 %>%
  mutate(OTU_ID = as.numeric(str_extract(`...1`, "\\d+$")))

# Step 3: Subset OTU table to top 17 and remove taxonomy column
top_otus <- otu_table %>%
  filter(OTU_ID %in% top10$OTU_ID) %>%
  select(-`Consensus Lineage`)

# Step 4: Convert OTU counts to relative abundance per sample
otu_rel_abund <- top_otus %>%
  column_to_rownames("OTU_ID") %>%
  mutate(across(everything(), as.numeric)) %>%
  sweep(2, colSums(.), FUN = "/") %>%
  rownames_to_column("OTU_ID") %>%
  pivot_longer(-OTU_ID, names_to = "SampleID", values_to = "RelAbundance") %>%
  mutate(OTU_ID = as.numeric(OTU_ID))

# Step 5: Extract species info from mapping
mapping <- mapping %>%
  mutate(Species = str_extract(Description, "E\\. intermedius|E\\. triangulatus"))

# Step 6: Merge abundance data with species metadata
otu_with_species <- otu_rel_abund %>%
  left_join(mapping %>% select(SampleID, Species), by = "SampleID")

# Step 7: Calculate mean relative abundance by OTU and species
mean_abundance <- otu_with_species %>%
  group_by(OTU_ID, Species) %>%
  summarise(MeanRelAbundance = mean(RelAbundance, na.rm = TRUE), .groups = "drop")

# Step 8: Order OTUs by total abundance for prettier sorting
mean_abundance <- mean_abundance %>%
  group_by(OTU_ID) %>%
  mutate(Total = sum(MeanRelAbundance)) %>%
  ungroup() %>%
  mutate(OTU_ID = fct_reorder(as.factor(OTU_ID), Total))

# Step 9: Final horizontal bar plot
ggplot(mean_abundance, aes(x = OTU_ID, y = MeanRelAbundance, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  coord_flip() +
  labs(
    x = "OTU ID",
    y = "Mean Relative Abundance",
    title = "Top 17 OTUs by Mean Abundance in Each Dung Beetle Species"
  ) +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62")) +
  theme_minimal(base_size = 14)

##############################################
#old code
# Step 1: Loading the data
top10 <- read_csv("./Random Forest/top_10.csv")
otu_table <- read_csv("otu_table.csv")

# Step 2: Extracting numeric OTU ID from '...1' column
top10 <- top10 %>%
  mutate(OTU_ID = as.numeric(str_extract(`...1`, "\\d+$")))

# Step 3: Joining Consensus Lineage from OTU table
top10_tax <- top10 %>%
  left_join(otu_table %>% select(OTU_ID, `Consensus Lineage`), by = "OTU_ID")

# Step 4: Extracting the most specific available taxonomic label
top10_tax <- top10_tax %>%
  mutate(Taxon_Label = map_chr(`Consensus Lineage`, function(lineage) {
    if (is.na(lineage)) return("Unknown")
    tax_levels <- str_split(lineage, ";")[[1]]
    tax_levels <- rev(tax_levels)  # from specific to general
    cleaned <- str_replace_all(tax_levels, ".*__", "") %>% str_trim()
    specific <- cleaned[!cleaned %in% c("uncultured", "unclassified", "", NA)]
    if (length(specific) > 0) specific[1] else "Unknown"
  })) %>%
  mutate(Taxon_Label = str_replace_all(Taxon_Label, "_", " "))

# Step 5: Creating the plot
ggplot(top10_tax, aes(x = reorder(Taxon_Label, MeanDecreaseAccuracy), y = MeanDecreaseAccuracy)) +
  geom_bar(stat = "identity", fill = "#2c7fb8") +
  coord_flip() +
  labs(
    x = "Bacterial Taxon",
    y = "Mean Decrease Accuracy",
    title = "Top OTUs by Importance in Species Classification"
  ) +
  theme_minimal(base_size = 14)



#Progress?
top10 <- read.csv("./Random Forest/top_10.csv")
otu_table <- read.csv("otu_table.csv", row.names = 1)
mapping <- read_tsv("mapping.txt", comment = "#", 
                    col_names = c("SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Description"))

mapping <- mapping %>%
  mutate(Species = str_extract(Description, "E\\. intermedius|E\\. triangulatus"))

# Step 2: Extracting numeric OTU ID from 'X' column
top10 <- top10 %>%
  mutate(OTU_ID = as.numeric(str_extract(`X`, "\\d+$")))

# Step 3: Joining Consensus Lineage from OTU table
top10_tax <- top10 %>%
  left_join(otu_table %>% select(1,63), by = "OTU_ID")

# Step 4: Extracting the most specific available taxonomic label
top10_tax <- top10_tax %>%
  mutate(Taxon_Label = map_chr(top10_tax[[7]], function(lineage) {
    if (is.na(lineage)) return("Unknown")
    tax_levels <- str_split(lineage, ";")[[1]]
    tax_levels <- rev(tax_levels)  # from specific to general
    cleaned <- str_replace_all(tax_levels, ".*__", "") %>% str_trim()
    specific <- cleaned[!cleaned %in% c("uncultured", "unclassified", "", NA)]
    if (length(specific) > 0) specific[1] else "Unknown"
  })) %>%
  mutate(Taxon_Label = str_replace_all(Taxon_Label, "_", " "))

top10_tax <- top10_tax %>%
  left_join(mapping %>% select(Species), by = "OTU_ID")

#changed otu ids to rows so now this wont work
top_otus <- otu_table %>%
  filter(OTU_ID %in% top10$OTU_ID) %>%
  select(-62)

top_otus$sum_abundance <- apply(top_otus[, -1], 1, sum, na.rm = TRUE)
print(top_otus[ ,63])

top_otus <- top_otus[, -c(2:62)]

New<- as.data.frame(t(otu_table))
col_names <- paste("X.OTU.ID.", colnames(New), sep="")
colnames(New) <- col_names
New$SampleID <- rownames(New)
merged<- cbind(, mapping)

top10_tax <- top10_tax %>%
  select(OTU_ID, everything())

top10_tax <- top10_tax %>%
  arrange(OTU_ID)

top_otus<- top_otus %>%
  left_join(top10_tax %>% select(Taxon_Label), by = "OTU_ID")

bound<- cbind(top_otus, top10_tax)

otu_table2<- otu_table[, -ncol(otu_table)]

otu_table2<- as.data.frame(t(otu_table))
col_names <- paste("OTU_ID", colnames(otu_table2), sep="")
colnames(otu_table2) <- col_names
otu_table2$SampleID <- rownames(otu_table2)


# Load the tibble package
library(tibble)
library(dplyr)

otu_with_species <- OTUnew %>%
  left_join(mapping %>% select(SampleID, Species), by = "SampleID") %>% 
  filter(!str_detect(Species, "untreated dung"))

ggplot(OTUnew, aes(x = OTU_ID, y = MeanRelAbundance, fill = Species)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  coord_flip() +
  labs(
    x = "OTU ID",
    y = "Mean Relative Abundance",
    title = "Top 17 OTUs by Mean Abundance in Each Dung Beetle Species"
  ) +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62")) +
  theme_minimal(base_size = 14)
