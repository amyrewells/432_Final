# This code transforms the absolute abundance OTU data from otu_table.csv 
# into relative abundance OTU data, and saves as otu_rel.csv

library(dplyr)

absOTU <- read.csv("./data/otu_table.csv", check.names = FALSE, row.names= "OTU_ID")

numeric_cols <- absOTU %>%
  select(where(is.numeric)) %>%
  colnames()

# Convert count data to relative abundance (for each OTU per sample)
relOTU <- absOTU %>%
  mutate(across(all_of(numeric_cols), ~ (. / sum(.)) * 1000))

write.csv(relOTU, "./data/otu_rel.csv", row.names = TRUE)
