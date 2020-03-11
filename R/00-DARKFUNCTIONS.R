###  ### ### ### ### ### ### # ###
### DARKFUNCTIONS METAGENOMICS ###
### Auxiliary functions        ###
###  ### ### ### ### ### ### # ###


##### INITIALIZE #####

# Set working directory
setwd("~/Data/Helsinki/analyses/")

# Create list of variables
CAT.VARS <- c("Sample", "Layer", "Ecosystem", "Vegetation")
NUM.VARS <- c("SoilWater", "pH", "SOM")
PROXY.VARS <- c("SnowDepth", "Elevation")
FLUX.VARS <- c("CH4", "CO2", "N2O")


##### GENERAL #####

readMetadata <- function () {
  read_delim("00_METADATA/metadata.tsv", delim = "\t") %>%
    select(all_of(CAT.VARS), all_of(NUM.VARS), all_of(PROXY.VARS), all_of(FLUX.VARS)) %>%
    mutate(Sample = as.factor(Sample)) %>%
    mutate(Layer = as.factor(Layer)) %>%
    mutate(Ecosystem = as.factor(Ecosystem)) %>%
    mutate(Vegetation = as.factor(Vegetation)) %>%
    mutate_at(vars(NUM.VARS), log) %>%
    mutate_at(vars(PROXY.VARS), log)
}

summariseMeans <- function (DATA, METADATA) {
  bind_cols(METADATA,
            Mean = apply(DATA, 1, mean),
            SD = apply(DATA, 1, sd)) %>% 
    arrange(desc(Mean)) %>% 
    as_tibble
}

filterUnclassified <- function (x, y) {
  x %>% 
    filter(!!rlang::sym(y) != "") %>% 
    filter(!grepl("Unclassified", !!rlang::sym(y))) %>% 
    filter(!grepl("Incertae Sedis", !!rlang::sym(y)))
}

computeRich <- function (x) {
  tibble(Sample = SAMPLES,
         value = apply(x, 2, function (y) sum(y > 0))) %>% 
    full_join(metadata %>% select(CAT.VARS), by = "Sample")
}

runDeseq <- function (x, y, z) {
  x %>%
    DESeqDataSetFromMatrix(., y, as.formula(paste("~", z))) %>%
    DESeq(fitType = "local", quiet = T)
}

filterDeseq <- function (x, y) {
  x %>%
    filter(pvalue < y) %>%
    mutate(Change = ifelse(log2FoldChange < 0, "Negative", "Positive")) %>%
    arrange(Change)
}

readTax <- function () {
  bind_rows(read_delim("05_MAGs/GTDB/gtdbtk.bac120.summary.tsv", delim = "\t") %>% 
              select(user_genome, classification),
            read_delim("05_MAGs/GTDB/gtdbtk.ar122.summary.tsv", delim = "\t") %>% 
              select(user_genome, classification)) %>% 
    filter(user_genome %in% MAGS) %>% 
    arrange(match(user_genome, MAGS)) %>% 
    mutate(classification = gsub("[a-z]__", "", classification)) %>% 
    separate(classification, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
}

selectLevel3 <- function (LEVEL3) {
  MODULES <- modules %>% 
    filter(`Level 3` == LEVEL3) %>% 
    pull(`Level 4`) %>% 
    as.vector
  
  MAGS <- modules %>%
    filter(`Level 4` %in% MODULES) %>% 
    select(all_of(MAGS)) %>% 
    mutate_all(sum) %>% 
    slice(1) %>% 
    gather %>% 
    filter(value > 0) %>% 
    pull(key) %>% 
    as.vector
  
  TAXONOMY <- taxonomy %>% 
    filter(user_genome %in% MAGS) %>% 
    arrange(Domain, Phylum, Class, Order, Family, Genus, Species) %>% 
    column_to_rownames("user_genome")
  
  MAGS <- rownames(TAXONOMY) %>% 
    as.vector
  
  PRESENCE <- modules %>% 
    filter(`Level 4` %in% MODULES) %>% 
    select(all_of(MAGS))
  
  MODULES.HIER <- modules %>% 
    filter(`Level 4` %in% MODULES) %>% 
    select(starts_with("Level"))
  
  COVERAGE <- coverage %>% 
    select(Ecosystem, all_of(MAGS))
  
  PRESENCE <- data.frame(PRESENCE, row.names = rownames(PRESENCE), check.names = F)
  
  COVERAGE <- data.frame(COVERAGE, row.names = rownames(COVERAGE), check.names = F)
  
  return(list(LEVEL3 = LEVEL3, 
              MODULES = MODULES, 
              MAGS = MAGS, 
              TAXONOMY = TAXONOMY,
              PRESENCE = PRESENCE, 
              MODULES.HIER = MODULES.HIER,
              COVERAGE = COVERAGE))
}

selectLevel4 <- function (LEVEL4) {
  LEVEL4 <- KEGG.modules %>% 
    filter(str_detect(`Level 4`, LEVEL4)) %>% 
    pull(`Level 4`) %>% 
    as.vector
  
  MAGS <- genes %>%
    filter(`Level 4` %in% LEVEL4) %>% 
    select(all_of(MAGS)) %>% 
    mutate_all(sum) %>% 
    slice(1) %>% 
    gather %>% 
    filter(value > 0) %>% 
    pull(key) %>% 
    as.vector
  
  TAXONOMY <- taxonomy %>% 
    filter(user_genome %in% MAGS) %>% 
    arrange(Domain, Phylum, Class, Order, Family, Genus, Species) %>% 
    column_to_rownames("user_genome")
  
  MAGS <- rownames(TAXONOMY) %>% 
    as.vector
  
  PRESENCE <- genes %>% 
    filter(`Level 4` == LEVEL4) %>%
    select(all_of(MAGS))
  
  MODULES.HIER <- genes %>% 
    filter(`Level 4` == LEVEL4) %>%
    select(starts_with("Level"), KO, gene)
  
  COVERAGE <- coverage %>% 
    select(Ecosystem, all_of(MAGS))
  
  PRESENCE <- data.frame(PRESENCE, row.names = rownames(PRESENCE), check.names = F)
  
  COVERAGE <- data.frame(COVERAGE, row.names = rownames(COVERAGE), check.names = F)
  
  return(list(LEVEL4 = LEVEL4, 
              MAGS = MAGS, 
              TAXONOMY = TAXONOMY,
              PRESENCE = PRESENCE, 
              MODULES.HIER = MODULES.HIER,
              COVERAGE = COVERAGE))
}


##### PLOTTING #####

devClose <- function () {
  while (!is.null(dev.list())) {
    dev.off()
  }
}

plotBoxplot <- function (DATA, MAPPING) {
  ggplot(DATA, aes(x = !! enquo(MAPPING), y = value)) +
    geom_boxplot(aes(fill = forcats::fct_rev(Layer)), outlier.shape = NA, position = position_dodge(preserve = "single")) +
    geom_point(aes(fill = forcats::fct_rev(Layer)), shape = 16, position = position_jitterdodge()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.title = element_blank(),
          axis.line = element_blank(), panel.border = element_rect(fill = NA), strip.background = element_rect(fill = "grey90")) +
    scale_fill_manual(name = "Layer", values = c("#98c699", "#b7d8ff"))
}

plotOrdination <- function (DATA, METADATA, MAPPING) {
  attach(METADATA, warn.conflicts = F)
  
  if (MAPPING == "Layer") {
    COLORS <- c("#b7d8ff", "#98c699")
  }
  if (MAPPING == "Vegetation") {
    COLORS <- c("#dfc3f8", "#beefc1", "#eca6c1", "#61b7d9", "#f9b99f")
  }
  if (MAPPING == "Ecosystem") {
    COLORS <- c("#dfc3f8", "#eca6c1", "#f9b99f")
  }
  
  plot(DATA, display = "sites", type = "n")
  
  if (MAPPING == "Layer") {
    points(DATA, display = "sites", pch = c(16, 17)[Layer], col = COLORS[Layer])
    ordiellipse(DATA, groups = Layer, kind = "sd", draw = "polygon", col = COLORS)
    legend("bottomleft", legend = levels(metadata$Layer), bty = "n", col = COLORS, pch = c(16, 17))
  }
  if (MAPPING == "Vegetation") {
    points(DATA, display = "sites", pch = c(16, 17)[Layer], col = COLORS[Vegetation])
    ordiellipse(DATA, groups = Vegetation, kind = "sd", draw = "polygon", col = COLORS)
    legend("bottomleft", legend = levels(metadata$Layer), bty = "n", pch = c(16, 17))
    legend("bottomright", legend = levels(metadata$Vegetation), bty = "n", col = COLORS, pch =15)
  }
  if (MAPPING == "Ecosystem") {
    points(DATA, display = "sites", pch = c(16, 17)[Layer], col = COLORS[Ecosystem])
    ordiellipse(DATA, groups = Ecosystem, kind = "sd", draw = "polygon", col = COLORS)
    legend("bottomleft", legend = levels(metadata$Layer), bty = "n", pch = c(16, 17))
    legend("bottomright", legend = levels(metadata$Ecosystem), bty = "n", col = COLORS, pch =15)
  }
  
  detach(METADATA)
}

plotBarplot <- function (DATA) {
  ggplot(DATA, aes(x = Samplx, y = Count, fill = fct_reorder(Taxa, Count, mean))) +
    geom_bar(color = "black", size = 0.1, stat = "identity", position = "stack") +
    facet_grid(cols = vars(Vegetation), scales = "free_x", space = "free_x") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank(),
          axis.line = element_blank(), legend.key.size = unit(8, "pt"),
          legend.spacing.x = unit(4, "pt"), legend.title = element_blank(), panel.border = element_rect(fill = NA),
          strip.background = element_rect(fill = "grey90")) +
    ylab("Percentage of mapped reads") +
    guides(fill = guide_legend(ncol = 1)) +
    scale_y_continuous(expand = c(0,0), labels = scales::percent)
}

plotHeatmap <- function (DATA, ...) {
  ORDER <- metadata %>%
    arrange(Ecosystem, Layer) %>%
    pull(Sample) %>%
    as.vector

  ANNOTATION_COL <- data.frame(Layer = metadata %>% select(Layer), Ecosystem = metadata %>% select(Ecosystem), row.names = metadata %>% pull(Sample))

  pheatmap(DATA %>% select(ORDER) %>% sqrt, color = colorRampPalette(colors = c("#DEEBF7", "#4292C6", "#08306B"))(100),
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F,
           annotation_col = ANNOTATION_COL, ...)
}

plotHeatmapBin <- function (DATA, ...) {
  ORDER <- metadata %>%
    arrange(Ecosystem) %>%
    pull(Sample) %>%
    as.vector
  
  ANNOTATION_COL <- data.frame(metadata %>% select(Ecosystem), row.names = metadata %>% pull(Sample))
  
  pheatmap(DATA %>% select(ORDER), color = colorRampPalette(colors = c("red", "black", "green"))(100),
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, scale = "row",
           gaps_col = table(metadata$Ecosystem) %>% as.vector %>% cumsum,
           annotation_col = ANNOTATION_COL,
           gaps_row = subset(DATA, Change == "Negative") %>% nrow, ...)
}

plotHeatmapBin2 <- function (DATA, MAPPING, ...) {
  METADATA <- metadata %>% 
    filter(Layer == "organic")
  
  ORDER <- METADATA %>%
    arrange(!! rlang::sym(MAPPING)) %>%
    pull(Sample) %>%
    as.vector
  
  ANNOTATION_COL <- data.frame(METADATA %>% select(MAPPING), row.names = METADATA %>% pull(Sample))
  
  pheatmap(DATA %>% select(ORDER), color = colorRampPalette(colors = c("red", "black", "green"))(100),
           border_color = NA, cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, scale = "row",
           annotation_col = ANNOTATION_COL,
           gaps_row = subset(DATA, Change == "Negative") %>% nrow, ...)
}


plotHeatmapKEGG <- function (DATA, RATIO) {
  p1 <- pheatmap(DATA$PRESENCE, silent = T,
                 color = c("white", "orange"),
                 cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, legend = F,
                 labels_row = DATA$MODULES.HIER %>% pull(`Level 4`) %>% as.character,
                 labels_col = paste(DATA$TAXONOMY %>% pull(Phylum), rownames(DATA$TAXONOMY), sep = " | ") %>% as.character)
  
  p2 <- pheatmap(DATA$COVERAGE %>% select(-Ecosystem) %>% sqrt, main = DATA$LEVEL3, silent = T, 
                 color = colorRampPalette(colors = c("white", "yellow", "red"))(100), 
                 cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, legend = F, 
                 annotation_col = DATA$TAXONOMY %>% select(Phylum),
                 labels_row = DATA$COVERAGE %>% pull(Ecosystem) %>% as.character,
                 show_colnames = F)
  
  plot_grid(as.grob(p2), as.grob(p1), nrow = 2, align = "v", rel_heights = RATIO)
}

plotHeatmapKEGG2 <- function (DATA, RATIO) {
  p1 <- pheatmap(DATA$PRESENCE, silent = T,
                 color = c("white", "orange"), 
                 cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, legend = F,
                 labels_row = paste(DATA$MODULES.HIER %>% pull(KO), DATA$MODULES.HIER %>% pull(gene), sep = " | ") %>% as.character,
                 labels_col = paste(DATA$TAXONOMY %>% pull(Phylum), rownames(DATA$TAXONOMY), sep = " | ") %>% as.character)
  
  p2 <- pheatmap(DATA$COVERAGE %>% select(-Ecosystem) %>% sqrt, main = DATA$LEVEL4, silent = T, 
                 color = colorRampPalette(colors = c("white", "yellow", "red"))(100), 
                 cellheight = 10, cellwidth = 10, cluster_cols = F, cluster_rows = F, legend = F, 
                 annotation_col = DATA$TAXONOMY %>% select(Phylum),
                 labels_row = DATA$COVERAGE %>% pull(Ecosystem) %>% as.character,
                 show_colnames = F)
  
  plot_grid(as.grob(p2), as.grob(p1), nrow = 2, align = "v", rel_heights = RATIO)
}
