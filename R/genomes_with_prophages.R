# This file generates .csv files that contain data on genomes that were detected to contain prophages
# Alison requested this data on 2017-05-25 during our weekly meeting

library(tidyverse)

compressed1 <- select(refseq_by_prophage, Genome.ID) %>%
  unique()

compressed2 <- select(viromes_by_prophage, Genome.ID) %>%
  unique()

metadata <- read.csv("./data/genome_list.csv", header = TRUE) %>%
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

write.csv(Ref_data, file = "./tables/RefSeq_genomes_with_prophages.csv", row.names = FALSE)

write.csv(Vir_data, file = "./tables/Viromes_genomes_with_prophages.csv", row.names = FALSE)



#############
### PLOTS ###
#############

# Creates a plot of the distribution of all genomes in which at least one prophage was predicted
prophaged_genomes_data <- Ref_data %>%
  count(Phylum, sort = TRUE)

prophaged_genomes_plot <- dotchart(prophaged_genomes_data$n, labels = prophaged_genomes_data$Phylum, main = "Prophage distribution by phylum")

# Creates plots of the distribution of all genomes queried in the analysis by phylum
queried_genomes_data <- metadata %>%
  count(Phylum, sort = TRUE)

queried_genomes_plot <- dotchart(queried_genomes_data$n, labels = queried_genomes_data$Phylum, main = "All queried genomes by phylum")

# Creates a plot looking at distributions within the Proteobacteria
proteo_dist_data <- metadata %>%
  filter(Phylum == "Proteobacteria") %>%
  count(Class, sort = TRUE)

proteo_dist_plot <- dotchart(proteo_dist_data$n, labels = proteo_dist_data$Class, main = "Queried Proteobacteria distribution")

proteo_prophage_data <- Ref_data %>%
  filter(Phylum == "Proteobacteria") %>%
  count(Class, sort = TRUE)

proteo_prophage_plot <- dotchart(proteo_prophage_data$n, labels = proteo_prophage_data$Class, main = "Proteobacteria prophage distribution")
