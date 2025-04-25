# Title: Variant Analysis
# Author: Ariana Silva

# Load necessary libraries
# readxl: for reading Excel files
# dplyr: for data manipulation
library(readxl)
library(dplyr)

# Load the Excel file
data <- read_excel("/Users/arianasilva/Documents/GenomicaHumana/FAM5.xlsx")

# Inspect the structure and first few rows of the data
str(data)
head(data)

# Rename ambiguous genotype columns to clearer names
data <- data %>% 
  rename(
    `UASib1:ZYG` = `UASib:ZYG...51`,
    `UASib2:ZYG` = `UASib:ZYG...52`
  )

# Quality Filtering
# Keep variants with sufficient read depth, allelic balance, and genotype quality
# Remove synonymous variants
filtering <- data %>%
  filter(
    TR >= 10,               # Total reads ≥ 10
    RATIO >= 0.20,          # Allelic balance ≥ 20%
    GQ >= 30,               # Genotype quality ≥ 30
    ANNOTATION != "synonymous SNV"  # Remove synonymous variants
  )
dim(filtering)  # Check how many variants remain

# Convert population frequencies to numeric format
freq_filtering <- filtering %>%
  mutate_at(vars(TGP_FREQ, 
                 ESP_FREQ, 
                 EVE_ALT_FREQ), as.numeric)

# Frequency Filtering
# Keep only rare variants MAF < 1% or those not found in databases
freq_filtered <- freq_filtering %>%
  filter(
    TGP_FREQ < 0.01 | is.na(TGP_FREQ),
    ESP_FREQ < 0.01 | is.na(ESP_FREQ),
    EVE_ALT_FREQ < 0.01 | is.na(EVE_ALT_FREQ)
  )
dim(freq_filtered)

# Prediction Filtering 
# Keep variants predicted as damaging by SIFT and not benign by PolyPhen2
pred_filtered <- freq_filtered %>%
  filter(
    (SIFT_SCORE < 0.05 | SIFT_SCORE %in% c(".", NA)),
    (PPH2_PRED != "benign" | PPH2_PRED %in% c(".", NA))
  )
dim(pred_filtered)

# Prepare Genotype Data
# Keep only the genotype state for each individual
data_filtered <- pred_filtered %>%
  mutate(
    `Proband:ZYG` = sub(".*:", "", `Proband:ZYG`),
    `Father:ZYG` = sub(".*:", "", `Father:ZYG`),
    `Mother:ZYG` = sub(".*:", "", `Mother:ZYG`),
    `UASib1:ZYG` = sub(".*:", "", `UASib1:ZYG`),
    `UASib2:ZYG` = sub(".*:", "", `UASib2:ZYG`)
  )

# Recessive model filtering
# Identify variants that are homozygous in the proband but not in family members
data_recessive <- data_filtered %>% 
  filter(
    `Proband:ZYG` == "hom", 
    `Father:ZYG` != "hom", 
    `Mother:ZYG` != "hom",
    `UASib1:ZYG` != "hom", 
    `UASib2:ZYG` != "hom"
  ) 
dim(data_recessive)
View(data_recessive)

# Recessive X-linked model filtering
# Homozygous in proband, mother is a carrier, and father/siblings are not affected
data_xlinked <- data_filtered %>%
  filter(
    `Proband:ZYG` == "hom",
    (`Mother:ZYG` == "het" | `Mother:ZYG` == "na"),
    is.na(`Father:ZYG`) | `Father:ZYG` == "na",
    `UASib1:ZYG` != "hom",
    `UASib2:ZYG` != "hom"
  )
dim(data_xlinked)

# De novo mutation filtering
# Variant present in proband but absent in both parents
data_denovo <- data_filtered %>%
  filter(
    (`Proband:ZYG` == "het" | `Proband:ZYG` == "hom"),
    `Father:ZYG` == "na",
    `Mother:ZYG` == "na"
    # Uncomment the sibling filters if needed
    # `UASib1:ZYG` == "na",
    # `UASib2:ZYG` == "na"
  )
dim(data_denovo)
View(data_denovo)

# Combine candidate variants from recessive and de novo models for easier search in OMIM
variants <- merge(data_recessive, data_denovo, all = TRUE)
dim(variants)
