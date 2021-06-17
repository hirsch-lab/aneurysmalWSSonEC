library(dplyr)

original_endo <- read.csv("data/original/20171127endothelial.csv",header=TRUE,sep=";",stringsAsFactors = F)
str(original_endo)

# remove 0 flow
publication_endo <- original_endo %>%
  filter(Flow != 0)
str(publication_endo)

# Save
write.csv(publication_endo, file = "data/raw_data_endothelial.csv", row.names = FALSE)
