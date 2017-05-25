# This file generates .csv files that contain data on genomes that were detected to contain prophages
# Alison requested this data on 2017-05-25 during our weekly meeting

library(tidyverse)

compressed1 <- select(refseq_by_prophage, Genome.ID) %>%
  unique()

compressed2 <- select(viromes_by_prophage, Genome.ID) %>%
  unique()

metadata <- read.csv("../data/genome_list.csv", header = TRUE) %>%
  mutate(Genome.ID = IMG.Genome.ID)

Ref_data <- inner_join(compressed1, metadata, by = "Genome.ID") %>%
  select(Genome.Name = Genome.Name...Sample.Name,
         IMG.Genome.ID,
         Phylum,
         Class,
         Order,
         Family,
         Genus,
         Species,
         Genome.Size = Genome.Size.....assembled,
         GOLD.Sequencing.Strategy)

Vir_data <- inner_join(compressed2, metadata, by = "Genome.ID") %>%
  select(Genome.Name = Genome.Name...Sample.Name,
         IMG.Genome.ID,
         Phylum,
         Class,
         Order,
         Family,
         Genus,
         Species,
         Genome.Size = Genome.Size.....assembled,
         GOLD.Sequencing.Strategy)

write.csv(Ref_data, file = "../tables/RefSeq_genomes_with_prophages.csv", row.names = FALSE)

write.csv(Vir_data, file = "../tables/Viromes_genomes_with_prophages.csv", row.names = FALSE)