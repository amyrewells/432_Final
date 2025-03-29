# Libraries used:
library(tidyverse)
library(stringr)
library(purrr)

# Step 1: Load your data
otu_table <- read_csv("./data/input/otu_table.csv")
top10 <- read_csv("./data/top_10.csv")
mapping <- read_tsv("./data/input/mapping.txt", comment = "#", 
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

# Step 5: Subset OTU table to top 17 OTUs and drop Consensus Lineage
top_otus <- otu_table %>%
  filter(OTU_ID %in% top10$OTU_ID) %>%
  select(-`Consensus Lineage`)

otu_rel_abund <- top_otus %>%
  column_to_rownames("OTU_ID") %>%  # Convert OTU_ID to row names
  mutate(across(everything(), as.numeric)) %>%  # Convert all columns to numeric
  rownames_to_column("OTU_ID") %>%  # Convert row names back to a column
  pivot_longer(-OTU_ID, names_to = "SampleID", values_to = "Abundance") %>%  # Reshape to long format
  mutate(OTU_ID = as.numeric(OTU_ID))  # Convert OTU_ID to numeric

# Step 6: Extract species info from mapping file
mapping <- mapping %>%
  mutate(Species = str_extract(Description, "E\\. intermedius|E\\. triangulatus")) %>% 
  filter(!str_detect(Species, "untreated dung"))

# Step 7: Merge abundance data with species metadata
otu_with_species <- otu_rel_abund %>%
  left_join(mapping %>% select(SampleID, Species), by = "SampleID")

# Step 8: Calculate mean abundance for each OTU by species
sum_abundance <- otu_with_species %>%
  group_by(OTU_ID, Species) %>%
  summarise(Abundances = sum(Abundance), .groups = "drop")

abundanceData <- sum_abundance %>%
  filter(!is.na(sum_abundance[[2]]))

# Step 9: Adding taxon labels to data for plotting
top10_tax_long <- top10_tax %>%
  pivot_longer(cols = starts_with("Taxon_Label"),  # Assuming Taxon_Label columns start with 'Taxon_Label'
               names_to = "Taxon_Type", 
               values_to = "Taxon_Label") %>%
  select(OTU_ID, Taxon_Label)  # Keep only the relevant columns

final <- abundanceData %>%
  left_join(top10_tax_long %>% select(OTU_ID, Taxon_Label), by = "OTU_ID")

# Step 10: Plotting bar chart
barplot<- ggplot(final, aes(x = Taxon_Label, y = Abundances, fill = Species)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Bacteria",
    y = "Total Abundance") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  scale_fill_manual(values = c("#66c2a5", "#fc8d62")) +
  theme_minimal(base_size = 14)

# Save
save(barplot, file= "./Figures/barplot.RData")  