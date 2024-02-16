#Manipulating data
#Working directory "C:/Users/MOSES A/Documents/Data_R"

#Load libraries

library(dplyr)
library(tidyverse)
library(Seurat)
library(readr)

library(GEOquery)
#read data

cancer_data <- read_csv("breast_cancer/data/GSE183947_fpkm.csv")
dim(cancer_data)

#Get meta data

gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)

metadata <- pData(phenoData(gse[[1]]))

#Sub-setting meta data  and renaming columns and rows

metadata_subset <- select(metadata, c(1, 10,11, 17))

metadata_modified <- metadata_subset %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

#reshapping data

head(cancer_data)


cancer_data_long <- cancer_data %>%
  rename(gene = ...1) %>%
  gather(key = "samples", value = "FPKM", -gene) #Converting from Y

#Join data frame metadata and gene expression data
cancer_data_long <- cancer_data_long %>%
  left_join(., metadata_modified, by = c("samples" = "description"))


#exploring the data

cancer_data_long %>%
  filter(gene == "BRCA1" | gene == "BRCA2" )%>%
  group_by(gene, tissue) %>%
  summarise(mean_FPKM = mean(FPKM), meadian_FPKM = median(FPKM)) %>%
  arrange(mean_FPKM)


#Visualization of the Data

library(ggplot2)

#barplot

cancer_data_long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = samples, y = FPKM, fill = tissue)) + #., means used the value of the pipe function
  geom_col()

#Density Distribution of expression between tumor and normal tissue in BRCA1 gene
cancer_data_long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density(alpha = 0.3)

#boxplot

cancer_data_long %>%
  filter(gene == "BRCA1") %>%
  ggplot(., aes(x = metastasis, y = FPKM, fill = tissue)) +
  geom_boxplot()
  

#Scatter plot. Comparing correlation between the expression of two genes

cancer_data_long %>%
  filter(gene == "BRCA1" | gene == "BRCA2") %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, color = tissue)) +
  geom_point() +
  geom_smooth(method ="lm", se = FALSE)

#Heatmap Multiple comparision and expressioon of gene

gene_of_interest <- c("BRCA1", "BRCA2", "TP53", "ALK", "MYCN")

cancer_data_long %>%
  filter(gene %in% gene_of_interest) %>%
  ggplot(., aes(x= samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red")


#if you want to save....  ggsave(save the plot, filename= intended filename, width, height)






