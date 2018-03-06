library(tidyverse)

genome_id_phylum <- read.csv("./tables/formatted_genome_list.csv", header = TRUE) %>%
  select(Genome.ID, Phylum)

all_refseq_by_prophage <- read.csv("./tables/all_refseq_by_prophage.csv", header = TRUE) %>%
  
  # Adds phylum data to detected prophages so that we can look for phyla of interest
  left_join(genome_id_phylum, by = "Genome.ID") %>%
  
  # Selects phyla for which there were at least 10 genomes queried
  filter(Phylum %in% c("Firmicutes",
                       "Proteobacteria",
                       "Actinobacteria",
                       "Bacteroidetes",
                       "Cyanobacteria",
                       "Chloroflexi",
                       "Planctomycetes",
                       "Deferribacteres",
                       "Marinimicrobia",
                       "Thermotogae",
                       "Nitrospirae")
         ) %>%
  select(-Phylum)


all_viromes_by_prophage <- read.csv("./tables/all_viromes_by_prophage.csv", header = TRUE) %>%
  
  # Adds phylum data to detected prophages so that we can look for phyla of interest
  left_join(genome_id_phylum, by = "Genome.ID") %>%
  
  # Selects phyla for which there were at least 10 genomes queried
  filter(Phylum %in% c("Firmicutes",
                       "Proteobacteria",
                       "Actinobacteria",
                       "Bacteroidetes",
                       "Cyanobacteria",
                       "Chloroflexi",
                       "Planctomycetes",
                       "Deferribacteres",
                       "Marinimicrobia",
                       "Thermotogae",
                       "Nitrospirae")
  )