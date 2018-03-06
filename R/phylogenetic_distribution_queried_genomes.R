# Calculating phylogenetic distributions
genome_list <-  read.csv("./tables/formatted_genome_list.csv", header = TRUE) %>%
  
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

# Changes names of phyla to remove "Candidatus"
genome_list$Phylum[genome_list$Phylum %in% "Candidatus Marinimicrobia"] <- "Marinimicrobia"

# Calculates number of genomes per taxonomic group and creates a new table
taxonomic_abundance_alphabetical <- genome_list %>%
  
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
write.csv(taxonomic_abundance_alphabetical, file = "./tables/summaries/taxonomic_abundance_alphabetacal.csv", row.names = FALSE)

# Reorganizes table based on numerical abundance
taxonomic_abundance_numerical <- taxonomic_abundance_alphabetical %>%
  arrange(-Genomes.in.phylum, -Genomes.in.class, -Genomes.in.order, -Genomes.in.family, -Genomes.in.genus)

# Writes numerical table to a file
write.csv(taxonomic_abundance_numerical, file = "./tables/summaries/taxonomic_abundance_numerical.csv", row.names = FALSE)

