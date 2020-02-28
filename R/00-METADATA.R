library("tidyverse")
source("~/Data/Scripts/DARKFUNCTIONS/R/00-DARKFUNCTIONS.R")


##### IMPORT AND PROCESS DATA #####

# Read data
metadata <- read_delim("00_METADATA/metadata.tsv", delim = "\t") %>% 
  select(Sample, Layer, Vegetation, SoilWater, pH, SOM, SnowDepth, Elevation, CH4, CO2, N2O) %>% 
  mutate(Sample = as.factor(Sample)) %>% 
  mutate(Layer = as.factor(Layer)) %>% 
  mutate(Vegetation = as.factor(Vegetation))

# Transform variables to log
metadata <- metadata %>% 
  mutate_at(vars(NUM.VARS), log) %>% 
  mutate_at(vars(PROXY.VARS), log)

# Keep only samples with full metadata
metadata <- metadata %>% 
  drop_na(NUM.VARS, PROXY.VARS)


##### ORDINATION #####

library("vegan")

# Ordinate

## All samples
ord <- metadata %>% 
  select(NUM.VARS) %>%
  rda

## Only organic samples
SAMPLES.ORG <- metadata %>% 
  filter(Layer == "organic") %>% 
  pull(Sample) %>% 
  as.vector

ord.organic <- metadata %>% 
  filter(Sample %in% SAMPLES.ORG) %>% 
  select(NUM.VARS, PROXY.VARS) %>%
  rda

# Plot

## All samples
png("00_METADATA/metadata-PCA.png", width = 1600, height = 1600, res = 150)
plotOrdinationLayer(ord, metadata)
metadata %>% 
  select(NUM.VARS) %>% 
  envfit(ord, .) %>% 
  plot(col = "#e0a8d0")
dev.off()

## Only organic samples
png("00_METADATA/metadata-PCA-ORG.png", width = 1600, height = 1600, res = 150)
plotOrdinationLayer(ord.organic, metadata %>% filter(Sample %in% SAMPLES.ORG), ORDIELLIPSE = F, LEGEND = F)
metadata %>% 
  filter(Sample %in% SAMPLES.ORG) %>% 
  select(NUM.VARS, PROXY.VARS) %>% 
  envfit(ord.organic, .) %>% 
  plot(col = "#e0a8d0")
metadata %>% 
  filter(Sample %in% SAMPLES.ORG) %>% 
  select(FLUX.VARS) %>% 
  envfit(ord.organic, .) %>% 
  plot(col = "#b3a5cb")
dev.off()


##### CORRELATIONS #####

library("ggplot2")
library("PerformanceAnalytics")

# Prepare data
correl.organic <- metadata %>% 
  filter(Sample %in% SAMPLES.ORG) %>% 
  select(NUM.VARS, PROXY.VARS, FLUX.VARS)

# All against all
png("00_METADATA/metadata-correlation.png", width = 1600, height = 1600, res = 150)
chart.Correlation(correl.organic, histogram = F)
dev.off()

# Only proxy variables
correl.organic.snow <- gather(correl.organic, key = "key", value = "value", -SnowDepth)
correl.organic.elevation <- gather(correl.organic, key = "key", value = "value", -Elevation)

png("00_METADATA/metadata-correlation-snow.png", width = 1600, height = 1600, res = 150)
ggplot(correl.organic.snow, aes(x = SnowDepth, y = value, color = key, fill = key)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(vars(key), scale = "free") +
  theme(axis.title.y = element_blank())
dev.off()

png("00_METADATA/metadata-correlation-elevation.png", width = 1600, height = 1600, res = 150)
ggplot(correl.organic.elevation, aes(x = Elevation, y = value, color = key, fill = key)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(vars(key), scale = "free") +
  theme(axis.title.y = element_blank())
dev.off()



##### TEST 19 CLUSTERS #####

clusters19 <- metadata %>% 
  filter(Sample %in% SAMPLES.ORG) %>% 
  column_to_rownames("Sample") %>% 
  select(SoilWater, SOM, pH, SnowDepth)

# UPGMA
upgma <- dist(clusters19, method = "euclidean")
upgma <- hclust(upgma, method = "average")
upgma.groups <- cutree(upgma, k = 19) 

# PCA
library(vegan)

pca <- rda(clusters19)

png("00_METADATA/19_clusters.png", width = 1600, height = 1600, res = 150)
par(mar = c(1, 1, 1, 1))
ordipointlabel(pca, choices = c(1,2),
               col = c("forestgreen", "navy"),
               pch = 1:2, font = c(1,2))
ordihull(pca, upgma.groups)
dev.off()

