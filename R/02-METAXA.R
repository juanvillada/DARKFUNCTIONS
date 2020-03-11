# Initialize
library("tidyverse")
source("~/Data/Scripts/DARKFUNCTIONS/R/00-DARKFUNCTIONS.R")


##### IMPORT AND PROCESS DATA #####

# Create list of samples
SAMPLES <- read_lines("SAMPLES.txt")

# Read metadata
metadata <- readMetadata() %>% 
  filter(Sample %in% SAMPLES)

# Read data
phylum <- read_delim("02_METAXA/summary_level_2.txt", delim = "\t")
genus  <- read_delim("02_METAXA/summary_level_6.txt", delim = "\t")

# Split taxonomy
phylum <- phylum %>% 
  separate(Taxa, c("Domain", "Phylum"), ";")

genus <- genus %>% 
  separate(Taxa, c("Domain", "Phylum", "Class", "Order", "Family", "Genus"), ";")

# Keep only bacteria and archaea
phylum <- phylum %>% 
  filter(Domain %in% c("Bacteria", "Archaea"))

genus <- genus %>% 
  filter(Domain %in% c("Bacteria", "Archaea"))

# Split into counts and taxonomic hierarchy
phylum.hier <- phylum %>% 
  select(-SAMPLES)

genus.hier <- genus %>% 
  select(-SAMPLES)

phylum <- phylum %>% 
  select(SAMPLES)

genus <- genus %>% 
  select(SAMPLES)

# Reorder samples
SAMPLES <- metadata %>% 
  pull(Sample) %>% 
  as.vector

phylum <- phylum %>% 
  select(SAMPLES)

genus <- genus %>% 
  select(SAMPLES)

# Compute number of reads per sample
nreads <- tibble(Sample = SAMPLES,
                 Reads = colSums(phylum)) %>% 
  mutate(Proportion = (Reads / 2000000) * 100)

# Transform to relative abundance
phylum.rel <- phylum %>% 
  mutate_all(function (x) x / sum(x))

genus.rel <- genus %>% 
  mutate_all(function (x) x / sum(x))

# Summarise means and standard deviations
phylum.sum <- phylum.rel %>% 
  summariseMeans(phylum.hier) %>% 
  filterUnclassified("Phylum")

phylum.sum.barren <- phylum.rel %>% 
  select(metadata %>% filter(Ecosystem == "barren") %>% pull(Sample)) %>% 
  summariseMeans(phylum.hier) %>% 
  filterUnclassified("Phylum")

phylum.sum.heath <- phylum.rel %>% 
  select(metadata %>% filter(Ecosystem == "heathland") %>% pull(Sample)) %>% 
  summariseMeans(phylum.hier) %>% 
  filterUnclassified("Phylum")

phylum.sum.wet <- phylum.rel %>% 
  select(metadata %>% filter(Ecosystem == "wetland") %>% pull(Sample)) %>% 
  summariseMeans(phylum.hier) %>% 
  filterUnclassified("Phylum")

genus.sum <- genus.rel %>% 
  summariseMeans(genus.hier) %>% 
  filterUnclassified("Genus")

genus.sum.barren <- genus.rel %>% 
  select(metadata %>% filter(Ecosystem == "barren") %>% pull(Sample)) %>% 
  summariseMeans(genus.hier) %>% 
  filterUnclassified("Genus")

genus.sum.heath <-genus.rel %>% 
  select(metadata %>% filter(Ecosystem == "heathland") %>% pull(Sample)) %>% 
  summariseMeans(genus.hier) %>% 
  filterUnclassified("Genus")

genus.sum.wet <- genus.rel %>% 
  select(metadata %>% filter(Ecosystem == "wetland") %>% pull(Sample)) %>% 
  summariseMeans(genus.hier) %>% 
  filterUnclassified("Genus")


##### BARPLOT #####

library("ggplot2")
library("forcats")

# Prepare data
phylum.bar <- bind_cols(Phylum = phylum.hier %>% select(Phylum), 
                        phylum.rel)

# Reorder taxa by abundance
phylum.bar <- phylum.bar %>% 
  arrange(match(Phylum, phylum.sum %>% pull(Phylum)))

# Melt and merge with metadata
phylum.bar <- phylum.bar %>% 
  gather(key = "Sample", value = "Count", -Phylum) %>% 
  full_join(metadata %>% select(CAT.VARS), by = "Sample") %>% 
  dplyr::rename(Taxa = Phylum) %>% 
  mutate(Sample = as.factor(Sample))

# Plot
png("02_METAXA/METAXA-PHYLUM-barplot.png", width = 1600, height = 1200, res = 150)
plotBarplot(phylum.bar) +
  scale_fill_manual(values = c("#ffbb9b", "#47e0ff", "#fbc387", "#64cbff", "#f0c683", "#6fb1ec", "#cfdd90", "#ce9cde", "#b9f5b5", "#faa0ca",
                               "#7ff9e1", "#ffaba1", "#59f0fc", "#d1a377", "#3ab9dc", "#ffddae", "#46d0c3", "#ffbdce", "#74dcb7", "#e9ccff",
                               "#91ba76", "#c3a1bf", "#e2ffb9", "#97acd4", "#73bc87", "#f0e3ff", "#68bb9e", "#d89d9a", "#c1ffdd", "#bcaa77",
                               "#bae8ff", "#fff4c2", "#75b8ae", "#eeffe2", "#91b3a0", "#c5fff0", "#a7af8c", "#b4e1e2"))
devClose()


##### BOXPLOT RICHNESS #####

library("ggplot2")

# Remove unclassified
genus.rich <- bind_cols(genus.hier, genus.rel) %>% 
  filterUnclassified("Genus") %>% 
  select(SAMPLES)

# Compute richness
genus.rich <- genus.rich %>% 
  computeRich()

# Plot
png("02_METAXA/METAXA-GENUS-richness.png", width = 1600, height = 1200, res = 150)
plotBoxplot(genus.rich, Ecosystem)
devClose()


##### HEATMAP #####

library("pheatmap")

# Prepare data
phylum.map <- bind_cols(phylum.hier, phylum.rel)
genus.map  <- bind_cols(genus.hier , genus.rel )

# Remove unclassified
phylum.map <- phylum.map %>% 
  filterUnclassified("Phylum")

genus.map <- genus.map %>% 
  filterUnclassified("Genus") 

# For the genus levelkeep only top 50 genera
genus.map <- genus.map %>% 
  arrange(desc(apply(genus.map %>% select(SAMPLES), 1, mean))) %>% 
  head(50)

# Reorder taxa
phylum.map <- phylum.map %>% 
  arrange(Domain)

genus.map  <- genus.map %>% 
  arrange(Phylum, Genus)

# Transform to data.frame
phylum.map <- data.frame(phylum.map, row.names = rownames(phylum.map), check.names = F)
genus.map  <- data.frame(genus.map , row.names = rownames(genus.map ), check.names = F)

# Plot
plotHeatmap(phylum.map,
            filename = "02_METAXA/METAXA-PHYLUM-heatmap.png",
            gaps_row = table(phylum.map$Domain) %>% as.vector %>% cumsum,
            gaps_col = table(metadata$Ecosystem) %>% as.vector %>% cumsum,
            annotation_row = phylum.map %>% select(Domain),
            annotation_colors = list(Ecosystem = c(barren = "#dfc3f8", heathland = "#eca6c1", wetland = "#f9b99f"),
                                     Layer = c(mineral = "#b7d8ff", organic = "#98c699"),
                                     Domain = c(Archaea = "#8befff", Bacteria = "#acaf79")),
            labels_row = phylum.map %>% pull(Phylum) %>% as.character)
devClose()

plotHeatmap(genus.map,
            filename = "02_METAXA/METAXA-GENUS-heatmap.png",
            gaps_row = table(genus.map$Phylum) %>% as.vector %>% cumsum,
            gaps_col = table(metadata$Ecosystem) %>% as.vector %>% cumsum,
            annotation_row = genus.map %>% select(Phylum),
            annotation_colors = list(Ecosystem = c(barren = "#dfc3f8", heathland = "#eca6c1", wetland = "#f9b99f"),
                                     Layer = c(mineral = "#b7d8ff", organic = "#98c699"),
                                     Phylum = c(Acidobacteria = "#b3ad91", Actinobacteria = "#cdb4fc",
                                                Bacteroidetes = "#e6e294", Chloroflexi = "#49bbc0",
                                                Euryarchaeota = "#ffddb6", Firmicutes = "#88f0cb",
                                                Gemmatimonadetes = "#f0ffc6", Planctomycetes = "#b3f4ff",
                                                Proteobacteria = "#7ac18a", Spirochaetes = "#d8ffee",
                                                Verrucomicrobia = "#93b48d")),
            labels_row = genus.map %>% pull(Genus) %>% as.character)
devClose()


##### ORDINATION #####

library("vegan")

# Prepare data
genus.ord <- bind_cols(genus.hier, genus.rel ) %>% 
  filterUnclassified("Genus") %>% 
  select(SAMPLES) %>% 
  sqrt

# NP-MANOVA
adonis(genus.ord %>% t ~ Layer*Ecosystem, metadata, permutations = 9999, method = "bray")

# NMDS
genus.mds <- genus.ord %>% 
  t %>% 
  metaMDS(distance = "bray")

genus.mds.min <- genus.ord %>% 
  select(metadata %>% filter(Layer == "mineral") %>% pull(Sample)) %>% 
  t %>% 
  metaMDS(distance = "bray")

genus.mds.org <- genus.ord %>% 
  select(metadata %>% filter(Layer == "organic") %>% pull(Sample)) %>% 
  t %>% 
  metaMDS(distance = "bray")

# Plot
png("02_METAXA/METAXA-NMDS.png", width = 1600, height = 1600, res = 150)
plotOrdination(genus.mds, metadata, "Ecosystem")
devClose()

png("02_METAXA/METAXA-NMDS-min.png", width = 1600, height = 1600, res = 150)
plotOrdination(genus.mds.min, metadata %>% filter(Layer == "mineral"), "Ecosystem")
devClose()

png("02_METAXA/METAXA-NMDS-org.png", width = 1600, height = 1600, res = 150)
plotOrdination(genus.mds.org, metadata %>% filter(Layer == "organic"), "Ecosystem")
metadata %>% 
  filter(Layer == "organic") %>% 
  select(FLUX.VARS) %>% 
  envfit(genus.mds.org, .) %>% 
  plot(col = "#b3a5cb")
devClose()


##### NEGATIVE BINOMIAL ANALYSES (VEGETATION) #####

library("DESeq2")
library("pheatmap")

# Run models
genus.bin.veg <- runDeseq(genus, metadata, "Ecosystem")

# Get normalised counts
genus.bin.veg.norm <- counts(genus.bin.veg, normalized = T) %>% 
  as_tibble

# Get results
genus.bin.veg <- results(genus.bin.veg) %>% 
  as_tibble

# Plot

## Prepare data
genus.bin.veg.map <- bind_cols(genus.bin.veg.norm, genus.bin.veg, genus.hier)

## Keep only p < 0.05, reorder by log2FoldChange and remove unclassified
genus.bin.veg.map <- genus.bin.veg.map %>% 
  arrange(Phylum) %>% 
  filterDeseq(0.05) %>% 
  filterUnclassified("Genus")

## Transform to data.frame
genus.bin.veg.map <- data.frame(genus.bin.veg.map, row.names = rownames(genus.bin.veg.map), check.names = F)

## Plot
plotHeatmapBin(genus.bin.veg.map, 
               filename = "02_METAXA/METAXA-heatmap-bin-ecosystem.png",
               annotation_row = genus.bin.veg.map %>% select(Phylum),
               annotation_colors = list(Ecosystem = c(barren = "#dfc3f8", heathland = "#eca6c1", wetland = "#f9b99f"),
                                        Layer = c(mineral = "#b7d8ff", organic = "#98c699"),
                                        Phylum = c(Acidobacteria = "#b3ad91", Actinobacteria = "#cdb4fc",
                                                   Bacteroidetes = "#e6e294", Chloroflexi = "#49bbc0",
                                                   Euryarchaeota = "#ffddb6", Firmicutes = "#88f0cb",
                                                   Gemmatimonadetes = "#f0ffc6", Planctomycetes = "#b3f4ff",
                                                   Proteobacteria = "#7ac18a", Spirochaetes = "#d8ffee",
                                                   Verrucomicrobia = "#93b48d")),
               labels_row = genus.bin.veg.map %>% pull(Genus) %>% as.character)
devClose()


##### NEGATIVE BINOMIAL ANALYSES (FLUX) #####

library("DESeq2")
library("pheatmap")

# Run models
genus.bin.ch <- runDeseq(genus %>% select(metadata %>% filter(Layer == "organic") %>% pull(Sample)), metadata %>% filter(Layer == "organic"), "CH4")
genus.bin.co <- runDeseq(genus %>% select(metadata %>% filter(Layer == "organic") %>% pull(Sample)), metadata %>% filter(Layer == "organic"), "CO2")
genus.bin.no <- runDeseq(genus %>% select(metadata %>% filter(Layer == "organic") %>% pull(Sample)), metadata %>% filter(Layer == "organic"), "N2O")

# Get normalised counts
genus.bin.ch.norm <- counts(genus.bin.ch, normalized = T) %>% 
  as_tibble

genus.bin.co.norm <- counts(genus.bin.co, normalized = T) %>% 
  as_tibble

genus.bin.no.norm <- counts(genus.bin.no, normalized = T) %>% 
  as_tibble

# Get results
genus.bin.ch <- results(genus.bin.ch) %>% 
  as_tibble

genus.bin.co <- results(genus.bin.co) %>% 
  as_tibble

genus.bin.no <- results(genus.bin.no) %>% 
  as_tibble

# Plot

## Prepare data
genus.bin.ch.map <- bind_cols(genus.bin.ch.norm, genus.bin.ch, genus.hier)
genus.bin.co.map <- bind_cols(genus.bin.co.norm, genus.bin.co, genus.hier)
genus.bin.no.map <- bind_cols(genus.bin.no.norm, genus.bin.no, genus.hier)

## Keep only p < 0.05, reorder by log2FoldChange and remove unclassified
genus.bin.ch.map <- genus.bin.ch.map %>% 
  arrange(Phylum) %>% 
  filterDeseq(0.05) %>% 
  filterUnclassified("Genus")

genus.bin.co.map <- genus.bin.co.map %>% 
  arrange(Phylum) %>% 
  filterDeseq(0.05) %>% 
  filterUnclassified("Genus")

genus.bin.no.map <- genus.bin.no.map %>% 
  arrange(Phylum) %>% 
  filterDeseq(0.05) %>% 
  filterUnclassified("Genus")

## Transform to data.frame
genus.bin.ch.map <- data.frame(genus.bin.ch.map, row.names = rownames(genus.bin.ch.map), check.names = F)
genus.bin.co.map <- data.frame(genus.bin.co.map, row.names = rownames(genus.bin.co.map), check.names = F)
genus.bin.no.map <- data.frame(genus.bin.no.map, row.names = rownames(genus.bin.no.map), check.names = F)

# Plot
plotHeatmapBin2(genus.bin.ch.map, "CH4",
                filename = "02_METAXA/METAXA-heatmap-bin-CH4.png",
                annotation_row = genus.bin.ch.map %>% select(Phylum),
                annotation_colors = list(Phylum = c(Acidobacteria = "#b3ad91", Actinobacteria = "#cdb4fc",
                                                    Bacteroidetes = "#e6e294", Chloroflexi = "#49bbc0",
                                                    Euryarchaeota = "#ffddb6", Firmicutes = "#88f0cb",
                                                    Gemmatimonadetes = "#f0ffc6", Planctomycetes = "#b3f4ff",
                                                    Proteobacteria = "#7ac18a", Spirochaetes = "#d8ffee",
                                                    Verrucomicrobia = "#93b48d", Nitrospirae = "#a3ada1")),
                labels_row = genus.bin.ch.map %>% pull(Genus) %>% as.character)
devClose()

plotHeatmapBin2(genus.bin.co.map, "CO2",
                filename = "02_METAXA/METAXA-heatmap-bin-CO2.png",
                annotation_row = genus.bin.co.map %>% select(Phylum),
                annotation_colors = list(Phylum = c(Acidobacteria = "#b3ad91", Actinobacteria = "#cdb4fc",
                                                    Bacteroidetes = "#e6e294", Chloroflexi = "#49bbc0",
                                                    Euryarchaeota = "#ffddb6", Firmicutes = "#88f0cb",
                                                    Gemmatimonadetes = "#f0ffc6", Planctomycetes = "#b3f4ff",
                                                    Proteobacteria = "#7ac18a", Spirochaetes = "#d8ffee",
                                                    Verrucomicrobia = "#93b48d", Nitrospirae = "#a3ada1")),
                labels_row = genus.bin.co.map %>% pull(Genus) %>% as.character)
devClose()


##### OLD (WITH NUMERICAL VARIABLES) #####

##### ORDINATION

# Create list of samples with full metadata
SAMPLES.FULL <- metadata %>% 
  select(Sample, NUM.VARS) %>% 
  drop_na %>% 
  pull(Sample) %>% 
  as.vector

# NP-MANOVA
adonis(genus.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ ., 
       metadata %>% filter(Sample %in% SAMPLES.FULL) %>% select(NUM.VARS), permutations = 9999, method = "bray")

# dbRDA
genus.rda <- ordiR2step(dbrda(genus.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ 1, metadata %>% filter(Sample %in% SAMPLES.FULL) %>% select(NUM.VARS), distance = "bray"), 
                        dbrda(genus.rel %>% select(SAMPLES.FULL) %>% t %>% sqrt ~ ., metadata %>% filter(Sample %in% SAMPLES.FULL) %>% select(NUM.VARS), distance = "bray"), "forward")

png("02_METAXA/METAXA-dbRDA-layer.png", width = 1600, height = 1600, res = 150)
plotOrdinationLayer(genus.rda, metadata %>% filter(Sample %in% SAMPLES.FULL))
text(genus.rda, display = "bp", col = "#e0a8d0")
devClose()

png("02_METAXA/METAXA-dbRDA-vegetation.png", width = 1600, height = 1600, res = 150)
plotOrdinationVeg(genus.rda, metadata %>% filter(Sample %in% SAMPLES.FULL))
text(genus.rda, display = "bp", col = "#e0a8d0")
devClose()


##### NEGATIVE BINOMIAL MODELLING

library("DESeq2")

# Run models
genus.bin.ph <- runDeseq(genus %>% select(SAMPLES.FULL), metadata %>% filter(Sample %in% SAMPLES.FULL), "pH")
genus.bin.wt <- runDeseq(genus %>% select(SAMPLES.FULL), metadata %>% filter(Sample %in% SAMPLES.FULL), "SoilWater")
genus.bin.om <- runDeseq(genus %>% select(SAMPLES.FULL), metadata %>% filter(Sample %in% SAMPLES.FULL), "SOM")

# Get normalised counts
genus.norm <- counts(genus.bin.ph, normalized = T) %>% 
    as_tibble

# Get results
genus.bin.ph <- results(genus.bin.ph) %>% 
  as_tibble

genus.bin.wt <- results(genus.bin.wt) %>% 
  as_tibble

genus.bin.om <- results(genus.bin.om) %>% 
  as_tibble


##### HEATMAP NEGATIVE BINOMIAL MODELLING

library("pheatmap")

# Prepare data
genus.bin.ph.map <- bind_cols(genus.norm, genus.bin.ph, genus.hier)
genus.bin.wt.map <- bind_cols(genus.norm, genus.bin.wt, genus.hier)
genus.bin.om.map <- bind_cols(genus.norm, genus.bin.om, genus.hier)

# Keep only p < 0.05, reorder by log2FoldChange and remove unclassified
genus.bin.ph.map <- filterDeseq(genus.bin.ph.map)
genus.bin.wt.map <- filterDeseq(genus.bin.wt.map)
genus.bin.om.map <- filterDeseq(genus.bin.om.map)

# Transform to data.frame
genus.bin.ph.map <- data.frame(genus.bin.ph.map, row.names = rownames(genus.bin.ph.map), check.names = F)
genus.bin.wt.map <- data.frame(genus.bin.wt.map, row.names = rownames(genus.bin.wt.map), check.names = F)
genus.bin.om.map <- data.frame(genus.bin.om.map, row.names = rownames(genus.bin.om.map), check.names = F)

# Plot
plotHeatmapBin(genus.bin.ph.map, "pH",
               filename = "02_METAXA/METAXA-pH-heatmap.png",
               annotation_row = genus.bin.ph.map["Change"],
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
               labels_row = as.character(paste(genus.bin.ph.map$Phylum,
                                               genus.bin.ph.map$Genus, sep = " / ")))

plotHeatmapBin(genus.bin.wt.map, "SoilWater",
               filename = "02_METAXA/METAXA-SoilWater-heatmap.png",
               annotation_row = genus.bin.wt.map["Change"],
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
               labels_row = as.character(paste(genus.bin.wt.map$Phylum,
                                               genus.bin.wt.map$Genus, sep = " / ")))

plotHeatmapBin(genus.bin.om.map, "SOM",
               filename = "02_METAXA/METAXA-SOM-heatmap.png",
               annotation_row = genus.bin.om.map["Change"],
               annotation_colors = list(Change = c(Negative = "#dbbec1", Positive = "#b0d4cd")),
               labels_row = as.character(paste(genus.bin.om.map$Phylum,
                                               genus.bin.om.map$Genus, sep = " / ")))

