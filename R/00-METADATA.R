library("tidyverse")
library("DARKFUNCTIONS.R")
setwd("~/Data/Helsinki/analyses/")


##### IMPORT AND PROCESS DATA #####

# Read data
metadata <- readMetadata("00_METADATA/metadata.tsv") %>% 
  drop_na(NUM.VARS, PROXY.VARS)


##### BOXPLOT #####

library("ggplot2")

# Prepare data
boxplot <- metadata %>% 
  select(all_of(NUM.VARS), all_of(CAT.VARS), all_of(FLUX.VARS)) %>% 
  mutate_at(vars(NUM.VARS), exp) %>% 
  gather(key = "key", value = "Value", -all_of(CAT.VARS))

# Plot
png("00_METADATA/metadata-boxplot.png", width = 1600, height = 1600, res = 150)
plotBoxplot(boxplot, Vegetation) +
  facet_grid(rows = vars(key), scales = "free")


##### CORRELATIONS #####

library("ggplot2")
library("PerformanceAnalytics")

# Prepare data
correl.min <- metadata %>% 
  filter(Layer == "mineral") %>% 
  select(all_of(NUM.VARS), all_of(PROXY.VARS))

correl.org <- metadata %>% 
  filter(Layer == "organic") %>% 
  select(all_of(NUM.VARS), all_of(PROXY.VARS))

# All against all
png("00_METADATA/metadata-correl-min.png", width = 1600, height = 1600, res = 150)
chart.Correlation(correl.min, histogram = F)

png("00_METADATA/metadata-correl-org.png", width = 1600, height = 1600, res = 150)
chart.Correlation(correl.org, histogram = F)

# Only proxy variables
correl.min.snow <- correl.min %>% 
  gather(key = "key", value = "value", -SnowDepth)

correl.min.elev <- correl.min %>% 
  gather(key = "key", value = "value", -Elevation)

correl.org.snow <- correl.org %>% 
  gather(key = "key", value = "value", -SnowDepth)

correl.org.elev <- correl.org %>% 
  gather(key = "key", value = "value", -Elevation)

png("00_METADATA/metadata-correl-min-snow.png", width = 1600, height = 1600, res = 150)
ggplot(correl.min.snow, aes(x = SnowDepth, y = value, color = key, fill = key)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(vars(key), scale = "free") +
  theme(axis.title.y = element_blank())

png("00_METADATA/metadata-correl-min-elevation.png", width = 1600, height = 1600, res = 150)
ggplot(correl.min.elev, aes(x = Elevation, y = value, color = key, fill = key)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(vars(key), scale = "free") +
  theme(axis.title.y = element_blank())

png("00_METADATA/metadata-correl-org-snow.png", width = 1600, height = 1600, res = 150)
ggplot(correl.org.snow, aes(x = SnowDepth, y = value, color = key, fill = key)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(vars(key), scale = "free") +
  theme(axis.title.y = element_blank())

png("00_METADATA/metadata-correl-org-elevation.png", width = 1600, height = 1600, res = 150)
ggplot(correl.org.elev, aes(x = Elevation, y = value, color = key, fill = key)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(vars(key), scale = "free") +
  theme(axis.title.y = element_blank())


##### ORDINATION #####

library("vegan")

# All samples
ord <- metadata %>% 
  select(NUM.VARS) %>%
  rda

png("00_METADATA/metadata-PCA.png", width = 1600, height = 1600, res = 150)
plotOrdination(ord, metadata, "Layer")
metadata %>% 
  select(all_of(NUM.VARS)) %>% 
  envfit(ord, .) %>% 
  plot(col = "#e0a8d0")

## Mineral samples
SAMPLES.MIN <- metadata %>% 
  filter(Layer == "mineral") %>% 
  pull(Sample) %>% 
  as.vector

ord.min <- metadata %>% 
  filter(Sample %in% SAMPLES.MIN) %>% 
  select(all_of(NUM.VARS)) %>%
  rda

png("00_METADATA/metadata-PCA-min.png", width = 1600, height = 1600, res = 150)
plotOrdination(ord.min, metadata %>% filter(Sample %in% SAMPLES.MIN), "Vegetation")
metadata %>% 
  filter(Sample %in% SAMPLES.MIN) %>% 
  select(NUM.VARS, PROXY.VARS) %>% 
  envfit(ord.min, .) %>% 
  plot(col = "#e0a8d0")

## Organic samples
SAMPLES.ORG <- metadata %>% 
  filter(Layer == "organic") %>% 
  pull(Sample) %>% 
  as.vector

ord.org <- metadata %>% 
  filter(Sample %in% SAMPLES.ORG) %>% 
  select(NUM.VARS) %>%
  rda

png("00_METADATA/metadata-PCA-org.png", width = 1600, height = 1600, res = 150)
plotOrdination(ord.org, metadata %>% filter(Sample %in% SAMPLES.ORG), "Vegetation")
metadata %>% 
  filter(Sample %in% SAMPLES.ORG) %>% 
  select(NUM.VARS, PROXY.VARS) %>% 
  envfit(ord.org, .) %>% 
  plot(col = "#e0a8d0")
metadata %>% 
  filter(Sample %in% SAMPLES.ORG) %>% 
  select(FLUX.VARS) %>% 
  envfit(ord.org, .) %>% 
  plot(col = "#b3a5cb")


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

