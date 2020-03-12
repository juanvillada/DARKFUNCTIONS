library("tidyverse")
# library("lazyeval")
library("DARKFUNCTIONS.R")
setwd("~/Data/Helsinki/analyses/")

# Define group
GROUP <- "ACIDOBACTERIA" #OR
GROUP <- "ARCHAEA"

# Read data
ANI <- read_delim(paste("05_MAGs", GROUP, "PHYLOGENOMICS/PYANI/ANIb_percentage_identity.txt", sep = "/"), delim = "\t") %>% 
  filter(key != "Escherichia_coli_str__K_12_substr__MG1655_GCF_002843685_1") %>% 
  select(-"Escherichia_coli_str__K_12_substr__MG1655_GCF_002843685_1")

length <- read_delim(paste("05_MAGs", GROUP, "PHYLOGENOMICS/PYANI/ANIb_alignment_lengths.txt", sep = "/"), delim = "\t") %>% 
  filter(key != "Escherichia_coli_str__K_12_substr__MG1655_GCF_002843685_1") %>% 
  select(-"Escherichia_coli_str__K_12_substr__MG1655_GCF_002843685_1")

coverage <- read_delim(paste("05_MAGs", GROUP, "PHYLOGENOMICS/PYANI/ANIb_alignment_coverage.txt", sep = "/"), delim = "\t") %>% 
  filter(key != "Escherichia_coli_str__K_12_substr__MG1655_GCF_002843685_1") %>% 
  select(-"Escherichia_coli_str__K_12_substr__MG1655_GCF_002843685_1")

annotation <- read_delim(paste("05_MAGs", GROUP, "PHYLOGENOMICS/view_data.txt", sep = "/"), delim = "\t") %>% 
  filter(contig != "Escherichia_coli_str__K_12_substr__MG1655_GCF_002843685_1")

# Put everything in same order
TAXA <- annotation %>% 
  pull(contig)

ANI <- ANI %>% 
  filter(key %in% TAXA) %>% 
  arrange(match(key, TAXA)) %>% 
  select(all_of(TAXA))

length <- length %>% 
  filter(key %in% TAXA) %>% 
  arrange(match(key, TAXA)) %>% 
  select(all_of(TAXA))

coverage <- coverage %>% 
  filter(key %in% TAXA) %>% 
  arrange(match(key, TAXA)) %>% 
  select(all_of(TAXA))

# Create list of MAGs
MAGs <- TAXA[grepl("[bog|heath]_MAG", TAXA)]

# Summarise
summary <- list(ANI = bind_cols(Query = annotation %>% select(contig), ANI) %>% 
                  gather(key = "Reference", value = "ANI", -contig) %>% 
                  rename(Query = contig),
                length = bind_cols(Query = annotation %>% select(contig), length) %>% 
                  gather(key = "Reference", value = "Length", -contig) %>% 
                  rename(Query = contig),
                coverage = bind_cols(Query = annotation %>% select(contig), coverage) %>% 
                  gather(key = "Reference", value = "Coverage", -contig) %>% 
                  rename(Query = contig))

summary <- bind_cols(Query = summary[["ANI"]] %>% select(Query),
                     Reference = summary[["ANI"]] %>% select(Reference),
                     ANI = summary[["ANI"]] %>% select(ANI),
                     Length = summary[["length"]] %>% select(Length),
                     Coverage = summary[["coverage"]] %>% select(Coverage))

summary <- summary %>% 
  filter(Query %in% MAGs) %>% 
  filter(! Reference %in% MAGs) %>% 
  filter(Coverage >= 0.2) %>% 
  arrange(Query, desc(ANI))

write_delim(summary, paste("05_MAGs", GROUP,"PHYLOGENOMICS/ANI_parsed.txt", sep = "/"), delim = "\t")


##### HEATMAP #####

library("pheatmap")

ANI.map <- bind_cols(ANI,
                     contig = annotation %>% select(contig),
                     name = annotation %>% select(name),
                     division = annotation %>% select(division))

ANI.map <- data.frame(ANI.map, row.names = rownames(ANI.map))

pheatmap(ANI.map %>% select(all_of(TAXA)),
         # filename = paste("05_MAGs", GROUP, "PHYLOGENOMICS/ANI.pdf", sep = "/"),
         color = colorRampPalette(c("yellow", "yellow", "red"))(100),
         cellwidth = 5, cellheight = 5,
         cluster_rows = F, cluster_cols = F,
         annotation_row = ANI.map %>% select(division),
         annotation_col = data.frame(ANI.map %>% select(division), row.names = ANI.map %>% pull(contig)),
         show_colnames = F, show_rownames = F,
         annotation_names_row = F, annotation_names_col = F)
devClose()


##### TEST #####

test <- data.frame(ANI, name = annotation$contig) %>% 
  melt %>% 
  rename(QUERY = name, REFERENCE = variable, ANI = value)

test <- merge(annotation %>% select(contig, division), test, by.y = "QUERY", by.x =  "contig", sort = F) %>% 
  rename(QUERY = contig, QUERY.DIV = division)

test <- merge(annotation %>% select(contig, division), test, by.y = "REFERENCE", by.x =  "contig", sort = F) %>% 
  rename(REFERENCE = contig, REFERENCE.DIV = division)

test <- test %>% 
  mutate(COMP = as.factor(ifelse(QUERY.DIV == REFERENCE.DIV, "SAME", "DIFFERENT")))

t.test(ANI ~ COMP, test)

boxplot(test$ANI~test$COMP)

test.GP3 <- test %>% 
  filter(QUERY.DIV == "GP3" & REFERENCE.DIV == "GP3")

boxplot(test.GP3$ANI)
