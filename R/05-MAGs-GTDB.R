library("tidyverse")
library("ape")
library("phytools")
library("treeio")
library("ggtree")
library("DARKFUNCTIONS.R")
setwd("~/Data/Helsinki/analyses/")


##### IMPORT AND PROCESS DATA #####

# Create list of MAGs
MAGs <- read_delim("05_MAGs/renamed_MAGs.txt", delim = "\t", col_names = F) %>% 
  pull(X2)

# Read data
bac_tree <- read.tree("05_MAGs/GTDB/gtdbtk.bac120.classify.tree")

arc_tree <- read.tree("05_MAGs/GTDB/gtdbtk.ar122.classify.tree" )

bac_MAGs_meta <- readTax() %>% 
  filter(Domain == "Bacteria")

arc_MAGs_meta <- readTax() %>% 
  filter(Domain == "Archaea") %>% 
  filter(! MAG %in% c("heath_MAG_0016", "heath_MAG_0018"))

arc_GTDB_meta <- read_delim("05_MAGs/GTDB/ar122_metadata_r89.tsv", delim = "\t") %>% 
  select(taxon = accession, taxonomy = gtdb_taxonomy, name = ncbi_organism_name) %>% 
  mutate(taxonomy = gsub("[a-z]__", "", taxonomy)) %>% 
  separate(taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

bac_GTDB_meta <- read_delim("05_MAGs/GTDB/bac120_metadata_r89.tsv", delim = "\t") %>% 
  select(taxon = accession, taxonomy = gtdb_taxonomy, name = ncbi_organism_name) %>% 
  mutate(taxonomy = gsub("[a-z]__", "", taxonomy)) %>% 
  separate(taxonomy, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")


##### ARCHAEA TREE #####

# Prepare view data
arc_view_data <- bind_cols(contig = arc_MAGs_meta %>% pull(MAG),
                           name = arc_MAGs_meta %>% pull(MAG),
                           arc_MAGs_meta) %>% 
  mutate_all(list(~na_if(., "")))

# Subset tree
arc_tree <- keep.tip(arc_tree, arc_view_data %>% pull(contig))

# Export
write.tree(arc_tree, "05_MAGs/GTDB/gtdb.archaea.tre")
write_delim(arc_view_data, "05_MAGs/GTDB/archaea_view_data.txt", delim = "\t")


##### BACTERIA TREE #####

# Prepare view data
bac_view_data <- bind_cols(contig = bac_MAGs_meta %>% pull(MAG),
                           name = bac_MAGs_meta %>% pull(MAG),
                           bac_MAGs_meta) %>% 
  mutate_all(list(~na_if(., "")))

# Subset tree
bac_tree <- keep.tip(bac_tree, bac_view_data %>% pull(contig))

# Export
write.tree(bac_tree, "05_MAGs/GTDB/gtdb.bacteria.tre")
write_delim(bac_view_data, "05_MAGs/GTDB/bacteria_view_data.txt", delim = "\t")


##### SELECTED TREES #####

# N-fixers: archaea
NFIXERS_ARC <- c("bog_MAG_0006", "bog_MAG_0007", "bog_MAG_0017", "bog_MAG_0043")

NFIXERS_ARC_TIPS <- vector()
for (MAG in NFIXERS_ARC) {
  res <- getSisters(arc_tree, node = MAG, mode = "number")
  res <- tree_subset(arc_tree, res, levels_back = 2)
  res <- res$tip.label
  
  NFIXERS_ARC_TIPS <- c(NFIXERS_ARC_TIPS, res)
}

nfixers_arc_tree <- keep.tip(arc_tree, NFIXERS_ARC_TIPS)

nfixers_arc_annot <- bind_rows(arc_GTDB_meta %>% subset(taxon %in% nfixers_arc_tree$tip.label),
                               arc_MAGs_meta %>% subset(MAG %in% nfixers_arc_tree$tip.label) %>% rename(taxon = MAG)) %>% 
  mutate(name = ifelse(! name %in% NA, paste(name, taxon, sep = "; "), taxon))

nfixers_arc_tree %>%
  ggtree(layout = "rectangular") %<+%
  nfixers_arc_annot +
  geom_tiplab(aes(label = name %>% as.character, color = Genus))+
  ggplot2::xlim(0, 2) +
  theme(legend.position="right") +
  ggtitle("Nitrogen fixers - archaea")

ggsave("05_MAGs/PLOTS/GTDB-nfixers-archaea.pdf")

# N-fixers: bacteria
NFIXERS_BAC <- c("bog_MAG_0046", "bog_MAG_0056", "bog_MAG_0061", "bog_MAG_0063", "bog_MAG_0075", "bog_MAG_0151", "bog_MAG_0153", 
                 "bog_MAG_0155", "bog_MAG_0156", "bog_MAG_0162", "bog_MAG_0168", "bog_MAG_0169", "bog_MAG_0173", "bog_MAG_0195", "bog_MAG_0205", 
                 "bog_MAG_0212", "bog_MAG_0235", "bog_MAG_0264", "bog_MAG_0282", "bog_MAG_0284", "bog_MAG_0285", "bog_MAG_0287", "bog_MAG_0288", 
                 "bog_MAG_0312")

NFIXERS_BAC_TIPS <- vector()
for (MAG in NFIXERS_BAC) {
  res <- getSisters(bac_tree, node = MAG, mode = "number")
  res <- tree_subset(bac_tree, res, levels_back = 1)
  res <- res$tip.label
  
  NFIXERS_BAC_TIPS <- c(NFIXERS_BAC_TIPS, res)
}

nfixers_bac_tree <- keep.tip(bac_tree, NFIXERS_BAC_TIPS)

nfixers_bac_annot <- bind_rows(bac_GTDB_meta %>% subset(taxon %in% nfixers_bac_tree$tip.label),
                               bac_MAGs_meta %>% subset(MAG %in% nfixers_bac_tree$tip.label) %>% rename(taxon = MAG)) %>% 
  mutate(name = ifelse(! name %in% NA, paste(name, taxon, sep = "; "), taxon)) %>% 
  mutate(color = paste(Phylum, Class, sep = ";"))

nfixers_bac_tree %>%
  ggtree(layout = "rectangular") %<+%
  nfixers_bac_annot +
  geom_tiplab(aes(label = name %>% as.character, color = color))+
  ggplot2::xlim(0, 4) +
  theme(legend.position="right") +
  ggtitle("Nitrogen fixers - bacteria")
ggsave("05_MAGs/PLOTS/GTDB-nfixers-bacteria.pdf")

# Methanotrophs: archaea
CH4TROPHS_ARC <- c("bog_MAG_0179")

CH4TROPHS_ARC_TIPS <- vector()
for (MAG in CH4TROPHS_ARC) {
  res <- getSisters(arc_tree, node = MAG, mode = "number")
  res <- tree_subset(arc_tree, res, levels_back = 3)
  res <- res$tip.label
  
  CH4TROPHS_ARC_TIPS <- c(CH4TROPHS_ARC_TIPS, res)
}

ch4trophs_arc_tree <- keep.tip(arc_tree, CH4TROPHS_ARC_TIPS)

ch4trophs_arc_annot <- bind_rows(arc_GTDB_meta %>% subset(taxon %in% ch4trophs_arc_tree$tip.label),
                                 arc_MAGs_meta %>% subset(MAG %in% ch4trophs_arc_tree$tip.label) %>% rename(taxon = MAG)) %>% 
  mutate(name = ifelse(! name %in% NA, paste(name, taxon, sep = "; "), taxon))

ch4trophs_arc_tree %>%
  ggtree(layout = "rectangular") %<+%
  ch4trophs_arc_annot +
  geom_tiplab(aes(label = name %>% as.character, color = Genus))+
  ggplot2::xlim(0, 1) +
  theme(legend.position="right") +
  ggtitle("Methanotrophs - archaea")

ggsave("05_MAGs/PLOTS/GTDB-ch4trophs-archaea.pdf")

# Methanotrophs: bacteria
CH4TROPHS_BAC <- c("bog_MAG_0057", "bog_MAG_0100", "bog_MAG_0095", "heath_MAG_0055", "bog_MAG_0195", "bog_MAG_0056")

CH4TROPHS_BAC_TIPS <- vector()
for (MAG in CH4TROPHS_BAC) {
  res <- getSisters(bac_tree, node = MAG, mode = "number")
  res <- tree_subset(bac_tree, res, levels_back = 2)
  res <- res$tip.label
  
  CH4TROPHS_BAC_TIPS <- c(CH4TROPHS_BAC_TIPS, res)
}

ch4trophs_bac_tree <- keep.tip(bac_tree, CH4TROPHS_BAC_TIPS)

ch4trophs_bac_annot <- bind_rows(bac_GTDB_meta %>% subset(taxon %in% ch4trophs_bac_tree$tip.label),
                                 bac_MAGs_meta %>% subset(MAG %in% ch4trophs_bac_tree$tip.label) %>% rename(taxon = MAG)) %>% 
  mutate(name = ifelse(! name %in% NA, paste(name, taxon, sep = "; "), taxon)) %>% 
  mutate(color = paste(Phylum, Class, sep = ";"))

ch4trophs_bac_tree %>%
  ggtree(layout = "rectangular") %<+%
  ch4trophs_bac_annot +
  geom_tiplab(aes(label = name %>% as.character, color = color)) +
  ggplot2::xlim(0, 3.5) +
  theme(legend.position="right") +
  ggtitle("Methanotrophs - bacteria")

ggsave("05_MAGs/PLOTS/GTDB-ch4trophs-bacteria.pdf")
