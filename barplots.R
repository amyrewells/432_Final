# 1st bar plot: Top OTUs by importance in Species Classification:
# Libraries used:
library(tidyverse)
library(stringr)
library(purrr)

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

#2nd barplot (there's 2 versions of it):

# Load libraries
library(tidyverse)

# Step 1: Load your data
otu_table <- read_csv("otu_table.csv")
top10 <- read_csv("./Random Forest/top_10.csv")
mapping <- read_tsv("mapping.txt", comment = "#", 
                    col_names = c("SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Description"))

# Step 2: Extract numeric OTU IDs from the top 17 list
top10 <- top10 %>%
  mutate(OTU_ID = as.numeric(str_extract(`...1`, "\\d+$")))

# Step 3: Subset OTU table to top 17 OTUs and drop Consensus Lineage
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

# Step 5: Extract species info from mapping file
mapping <- mapping %>%
  mutate(Species = str_extract(Description, "E\\. intermedius|E\\. triangulatus"))

# Step 6: Merge abundance data with species metadata
otu_with_species <- otu_rel_abund %>%
  left_join(mapping %>% select(SampleID, Species), by = "SampleID")

# Step 7: Calculate mean abundance for each OTU by species
mean_abundance <- otu_with_species %>%
  group_by(OTU_ID, Species) %>%
  summarise(MeanRelAbundance = mean(RelAbundance, na.rm = TRUE), .groups = "drop")

# Step 8: Plot the grouped bar chart
ggplot(mean_abundance, aes(x = factor(OTU_ID), y = MeanRelAbundance, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "OTU ID",
    y = "Mean Relative Abundance",
    title = "Top 17 OTUs by Mean Abundance in Each Dung Beetle Species"
  ) +
  scale_fill_manual(values = c("#66c2a5", "#fc8d62")) +
  theme_minimal(base_size = 14)


# 2nd version:

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

#Update:
top10 <- read_csv("./Random Forest/top_10.csv")
otu_table <- read_csv("otu_table.csv")
mapping <- read_tsv("mapping.txt", comment = "#", 
                    col_names = c("SampleID", "BarcodeSequence", "LinkerPrimerSequence", "Description"))

# Step 2: Extracting numeric OTU ID from '...1' column
top10 <- top10 %>%
  mutate(OTU_ID = as.numeric(str_extract(`...1`, "\\d+$")))

# Step 3: Joining Consensus Lineage from OTU table
top10_tax <- top10 %>%
  left_join(otu_table %>% select(OTU_ID, `Consensus Lineage`), by = "OTU_ID")

top10_tax <- top10_tax %>%
  left_join(mapping %>% select(Species), by = "OTU_ID")

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

mapping <- mapping %>%
  mutate(Species = str_extract(Description, "E\\. intermedius|E\\. triangulatus"))
  
top_otus <- otu_table %>%
  filter(OTU_ID %in% top10$OTU_ID) %>%
  select(-`Consensus Lineage`)

top_otus$mean_abundance <- apply(top_otus[, -1], 1, mean, na.rm = TRUE)

# Check the result
print(otu_data)


# Step 6: Merge abundance data with species metadata
OTUnew<- as.data.frame(t(top_otus))
# Set the first row as column names
colnames(OTUnew) <- OTUnew[1, ]

# Remove the first row (now that it's been set as column names)
OTUnew <- OTUnew[-1, ]

# Load the tibble package
library(tibble)
library(dplyr)

# Convert row names into a column
OTUnew <- rownames_to_column(OTUnew, var = "SampleID")

otu_with_species <- OTUnew %>%
  left_join(mapping %>% select(SampleID, Species), by = "SampleID") %>% 
  filter(!str_detect(Species, "untreated dung"))

# Reshape the data to make OTU numbers a single column
otu_with_species_long <- otu_with_species %>%
  pivot_longer(cols= 2:11,  # Select columns that start with "OTU"
               names_to = "OTU_ID",     # Name of the new column for OTU numbers
               values_to = "Abundance")     # Name of the new column for OTU values

otu_with_species_long_unique <- otu_with_species_long %>%
  distinct()

mean_abundance <- otu_with_species %>%
  group_by(OTU_ID, Species) %>%
  summarise(MeanRelAbundance = mean(RelAbundance, na.rm = TRUE), .groups = "drop")

# Step 8: Order OTUs by total abundance for prettier sorting
mean_abundance <- mean_abundance %>%
  group_by(OTU_ID) %>%
  mutate(Total = sum(MeanRelAbundance)) %>%
  ungroup() %>%
  mutate(OTU_ID = fct_reorder(as.factor(OTU_ID), Total))

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
