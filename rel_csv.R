library(dplyr)
setwd("C:/statistics/R/repositories/BIOL432_Final/")

absOTU <- read.csv("otu_table.csv", check.names = FALSE, row.names= "OTU_ID")

numeric_cols <- absOTU %>%
  select(where(is.numeric)) %>%
  colnames()

# Convert count data to relative abundance (for each OTU per sample)
relOTU <- absOTU %>%
  mutate(across(all_of(numeric_cols), ~ (. / sum(.)) * 1000))

write.csv(relOTU, "otu_rel.csv", row.names = TRUE)
