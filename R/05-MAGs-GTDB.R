library("tidyverse")
library("ape")
library("DARKFUNCTIONS.R")
setwd("~/Data/Helsinki/analyses/")


##### IMPORT AND PROCESS DATA #####

# Create list of MAGs
MAGs <- read_delim("05_MAGs/renamed_MAGs.txt", delim = "\t", col_names = F) %>% 
  pull(X2) %>% 
  as.vector

# Read data
bac_tree <- read.tree("05_MAGs/GTDB/gtdbtk.bac120.classify.tree")

arc_tree <- read.tree("05_MAGs/GTDB/gtdbtk.ar122.classify.tree" )

bac_MAGs_meta <- readTax() %>% 
  filter(Domain == "Bacteria")

arc_MAGs_meta <- readTax() %>% 
  filter(Domain == "Archaea") %>% 
  filter(! MAG %in% c("heath_MAG_0016", "heath_MAG_0018"))


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
