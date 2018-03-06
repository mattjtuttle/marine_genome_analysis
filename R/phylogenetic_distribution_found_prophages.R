################
# For RefSeq
################

# Calculating phylogenetic distributions
refseq_by_genome <-  read.csv("./tables/all_refseq_by_genome.csv", header = TRUE) %>%
  
  # Selects only taxonomic data
  select(Phylum,
         Class,
         Order,
         Family,
         Genus
  ) %>%
  
  # Chooses only phyla of interest
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
                       "Nitrospirae"
                       )
  )

# Calculates number of genomes per taxonomic group and creates a new table
refseq_prophages_taxonomy_alphabetical <- refseq_by_genome %>%
  
  # Calculates number of genomes queried per phylum
  group_by(Phylum) %>%
  mutate(
    Genomes.in.phylum = n()
  ) %>%
  ungroup() %>%
  
  # Calculates number of genomes queried per class
  group_by(Phylum, Class) %>%
  mutate(
    Genomes.in.class = n()
  ) %>%
  ungroup() %>%
  
  # Calculates number of genomes queried per class
  group_by(Phylum, Class, Order) %>%
  mutate(
    Genomes.in.order = n()
  ) %>%
  ungroup() %>%
  
  # Calculates number of genomes queried per class
  group_by(Phylum, Class, Order, Family) %>%
  mutate(
    Genomes.in.family = n()
  ) %>%
  ungroup() %>%
  
  # Calculates number of genomes queried per class
  group_by(Phylum, Class, Order, Family, Genus) %>%
  mutate(
    Genomes.in.genus = n()
  ) %>%
  ungroup() %>%
  
  # Condenses tables to unique rows
  distinct() %>%
  
  # Alphabatizes table
  arrange(Phylum, Class, Order, Family, Genus) %>%
  
  # Reorders columns
  select(Phylum,
         Genomes.in.phylum,
         Class,
         Genomes.in.class,
         Order,
         Genomes.in.order,
         Family,
         Genomes.in.family,
         Genus,
         Genomes.in.genus)

# Writes alphabetical table to a file
write.csv(refseq_prophages_taxonomy_alphabetical, file = "./tables/summaries/taxonomic_refseq_prophages_alphabetical.csv", row.names = FALSE)

# Reorganizes table based on numerical abundance
refseq_prophages_taxonomy_numerical <- refseq_prophages_taxonomy_alphabetical %>%
  arrange(-Genomes.in.phylum, -Genomes.in.class, -Genomes.in.order, -Genomes.in.family, -Genomes.in.genus)

# Writes numerical table to a file
write.csv(refseq_prophages_taxonomy_numerical, file = "./tables/summaries/taxonomic_refseq_prophages_numerical.csv", row.names = FALSE)





################
# For Viromes
################

# Calculating phylogenetic distributions
viromes_by_genome <-  read.csv("./tables/all_viromes_by_genome.csv", header = TRUE) %>%
  
  # Selects only taxonomic data
  select(Phylum,
         Class,
         Order,
         Family,
         Genus
  ) %>%
  
  # Chooses only phyla of interest
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
                       "Nitrospirae",
                       "Candidatus Marinimicrobia")
  )

# Calculates number of genomes per taxonomic group and creates a new table
viromes_prophages_taxonomy_alphabetical <- viromes_by_genome %>%
  
  # Calculates number of genomes queried per phylum
  group_by(Phylum) %>%
  mutate(
    Genomes.in.phylum = n()
  ) %>%
  ungroup() %>%
  
  # Calculates number of genomes queried per class
  group_by(Phylum, Class) %>%
  mutate(
    Genomes.in.class = n()
  ) %>%
  ungroup() %>%
  
  # Calculates number of genomes queried per class
  group_by(Phylum, Class, Order) %>%
  mutate(
    Genomes.in.order = n()
  ) %>%
  ungroup() %>%
  
  # Calculates number of genomes queried per class
  group_by(Phylum, Class, Order, Family) %>%
  mutate(
    Genomes.in.family = n()
  ) %>%
  ungroup() %>%
  
  # Calculates number of genomes queried per class
  group_by(Phylum, Class, Order, Family, Genus) %>%
  mutate(
    Genomes.in.genus = n()
  ) %>%
  ungroup() %>%
  
  # Condenses tables to unique rows
  distinct() %>%
  
  # Alphabatizes table
  arrange(Phylum, Class, Order, Family, Genus) %>%
  
  # Reorders columns
  select(Phylum,
         Genomes.in.phylum,
         Class,
         Genomes.in.class,
         Order,
         Genomes.in.order,
         Family,
         Genomes.in.family,
         Genus,
         Genomes.in.genus)

# Writes alphabetical table to a file
write.csv(viromes_prophages_taxonomy_alphabetical, file = "./tables/summaries/taxonomic_viromes_prophages_alphabetical.csv", row.names = FALSE)

# Reorganizes table based on numerical abundance
viromes_prophages_taxonomy_numerical <- viromes_prophages_taxonomy_alphabetical %>%
  arrange(-Genomes.in.phylum, -Genomes.in.class, -Genomes.in.order, -Genomes.in.family, -Genomes.in.genus)

# Writes numerical table to a file
write.csv(viromes_prophages_taxonomy_numerical, file = "./tables/summaries/taxonomic_viromes_prophages_numerical.csv", row.names = FALSE)

