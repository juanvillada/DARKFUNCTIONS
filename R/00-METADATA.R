library("dplyr")
source("00-DARKFUNCTIONS.R")


##### IMPORT AND PROCESS DATA #####

# Set working directory
setwd("~/Helsinki/analyses/")

# Read data
metadata <- read.table("00_METADATA/metadata_clean.tsv", header = T, sep = "\t")

# Create list of samples
SAMPLES <- metadata %>% 
  select(Sample) %>% 
  unlist %>% 
  as.vector

SAMPLES.FULL <- metadata %>% 
  na.omit %>% 
  select(Sample) %>% 
  unlist %>% 
  as.vector

SAMPLES.ORG <- metadataSub(SAMPLES.FULL) %>% 
  subset(Layer == "Organic") %>% 
  select(Sample) %>% 
  unlist %>% 
  as.vector

SAMPLES.MIN <- metadataSub(SAMPLES.FULL) %>% 
  subset(Layer == "Mineral") %>% 
  select(Sample) %>% 
  unlist %>% 
  as.vector

# Create list of numeric variables
NUM.VARS = metadata %>% 
  select(c(5:18)) %>% 
  names

# Transform numeric variables to log
metadata[NUM.VARS] <- log(metadata[NUM.VARS])


##### ORDINATION #####

library("vegan")

# Prepare data
all.ord <- metadataSub(SAMPLES.FULL) %>% 
  select(NUM.VARS) %>% 
  rda

organic.ord <- metadataSub(SAMPLES.ORG) %>% 
  select(NUM.VARS) %>% 
  rda

mineral.ord <- metadataSub(SAMPLES.MIN) %>% 
  select(NUM.VARS) %>% 
  rda

# Plot
png("00_METADATA/metadata-PCA.png", width = 1600, height = 1600, res = 150)
plotOrdinationLayer(all.ord, metadataSub(SAMPLES.FULL))
plot(envfit(all.ord, metadataSub(SAMPLES.FULL)[NUM.VARS]), col = "#e0a8d0")
dev.off()

png("00_METADATA/metadata-ORGANIC-PCA.png", width = 1600, height = 1600, res = 150)
plotOrdinationHabitat(organic.ord, metadataSub(SAMPLES.ORG))
plot(envfit(organic.ord, metadataSub(SAMPLES.ORG)[NUM.VARS]), col = "#e0a8d0")
dev.off()

png("00_METADATA/metadata-MINERAL-PCA.png", width = 1600, height = 1600, res = 150)
plotOrdinationHabitat(mineral.ord, metadataSub(SAMPLES.MIN))
plot(envfit(mineral.ord, metadataSub(SAMPLES.MIN)[NUM.VARS]), col = "#e0a8d0")
dev.off()

