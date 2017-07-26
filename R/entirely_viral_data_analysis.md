# Entirely viral sequence data analysis
Matt Tuttle  
July 11, 2017  



## Data information

The following data is comprised of sequences that VIRsorter detected as being solely viral (as opposed to integrated in a genome). Given that our original dataset of marine bacterial genomes was comprised of finished genomes as well as drafts genomes and SAGs, it is possible that these sequences annotated as being solely viral are in fact prophages since many of these genomes are highly fragmented and composed of many different contigs. According to the VirSorter paper (Roux, 2015), any contigs where 80% or more of a contig is identified as viral sequence, that contig is considered to be solely viral and not considered to be a prophage.

## Data import


```r
all_genome_data <- read.csv("../data/genome_list.csv", header = TRUE) %>%
  rename(Genome.ID = IMG.Genome.ID)

scaffold_data <- read.csv("../data/scaffold_list.csv", header = TRUE)

Ref_vir_only <- read.csv("../data/RefSeq_entirely_viral.csv", header = TRUE)

Vir_vir_only <- read.csv("../data/Viromes_entirely_viral.csv", header = TRUE)
```

## Detection of prophages at the prophage level

### Note

Prophage locations for these datasets cannot be calculated in the same way as the detected prophages since the program assumes that the entire scaffold that they prophage is found on is viral in nature. Based upon the headers found within the fasta file output, it also seems to assume that the scaffolds predicted as being entirely viral aare assumed to be circular in nature. However, I'm not sure if there is any additional information available to back this up.

RefSeq:


```r
# Makes changes to VirSorter output datatable to give more information about found prophages
Vir_refseq_by_prophage <- Ref_vir_only %>%
 
  # Removes prefix from Contig_id
  mutate(Contig_id = gsub("VIRSorter_", "", Contig_id)) %>%
  
  # Adds Scaffold.ID column to datatable by extracting it from the scaffold ID
  # Removing suffix and making string numeric allows this data to be combined with metadata
  mutate(Scaffold.ID = gsub("\\_.*", "", Contig_id)) %>%
  mutate(Scaffold.ID = as.numeric(Scaffold.ID)) %>%
  
  # Removes prefix from Fragment
  mutate(Fragment = gsub("VIRSorter_", "", Fragment)) %>%
  
  # Converts NAs in Nb.phage.hallmark.genes to zero
  mutate(Nb.phage.hallmark.genes = ifelse(is.na(Nb.phage.hallmark.genes), 0, Nb.phage.hallmark.genes)) %>%

  # Combines scaffold metadata with VirSorter output
  inner_join(scaffold_data, by = "Scaffold.ID") %>%
  
  # Calculates the number of prophages per scaffold
  mutate(Prophages.per.scaffold = 1) %>%
  group_by(Scaffold.ID) %>%
  mutate(Prophages.per.scaffold = sum(Prophages.per.scaffold)) %>%
  ungroup() %>%
  
  # Calculates the number of prophages per genome
  mutate(Prophages.per.genome = 1) %>%
  group_by(Genome.ID) %>%
  mutate(Prophages.per.genome = sum(Prophages.per.genome)) %>%
  ungroup()


# Calculates the total number of genomes found to contan at least 1 prophage
Vir_Refseq_genomes_with_prophage <- distinct(Vir_refseq_by_prophage, Genome.ID) %>%
  nrow()
  
# Writes the datatable to the tables folder
write.csv(Vir_refseq_by_prophage, file = "../tables/entirely_viral_refseq_by_prophage.csv", row.names = FALSE)

# Calculates the start and end points of prophage location and number of genes
# Since the entire fragment is annotated as being viral, the first gene and first bp are the starts of all these sequences
Vir_refseq_by_prophage_locations <- Vir_refseq_by_prophage %>%
  mutate(Gene.Start = 1,
         Gene.End = Nb.genes,
         bp.Start = 1,
         bp.End = Sequence.Length..bp.,
         Fragment.Length = Sequence.Length..bp.
         )

# Writes the datatable to the tables folder
write.csv(Vir_refseq_by_prophage_locations, file = "../tables/entirely_viral_refseq_by_prophage_locations.csv", row.names = FALSE)
```

Viromes:


```r
# Makes changes to VirSorter output datatable to give more information about found prophages
Vir_viromes_by_prophage <- Vir_vir_only %>%
 
  # Removes prefix from Contig_id
  mutate(Contig_id = gsub("VIRSorter_", "", Contig_id)) %>%
  
  # Adds Scaffold.ID column to datatable by extracting it from the scaffold ID
  # Removing suffix and making string numeric allows this data to be combined with metadata
  mutate(Scaffold.ID = gsub("\\_.*", "", Contig_id)) %>%
  mutate(Scaffold.ID = as.numeric(Scaffold.ID)) %>%
  
  # Removes prefix from Fragment
  mutate(Fragment = gsub("VIRSorter_", "", Fragment)) %>%
  
  # Converts NAs in Nb.phage.hallmark.genes to zero
  mutate(Nb.phage.hallmark.genes = ifelse(is.na(Nb.phage.hallmark.genes), 0, Nb.phage.hallmark.genes)) %>%

  # Combines scaffold metadata with VirSorter output
  inner_join(scaffold_data, by = "Scaffold.ID") %>%
  
  # Calculates the number of prophages per scaffold
  mutate(Prophages.per.scaffold = 1) %>%
  group_by(Scaffold.ID) %>%
  mutate(Prophages.per.scaffold = sum(Prophages.per.scaffold)) %>%
  ungroup() %>%
  
  # Calculates the number of prophages per genome
  mutate(Prophages.per.genome = 1) %>%
  group_by(Genome.ID) %>%
  mutate(Prophages.per.genome = sum(Prophages.per.genome)) %>%
  ungroup()

# Calculates the total number of genomes found to contan at least 1 prophage
Vir_Viromes_genomes_with_prophage <- distinct(Vir_viromes_by_prophage, Genome.ID) %>%
  nrow()
  
# Writes the datatable to the tables folder
write.csv(Vir_viromes_by_prophage, file = "../tables/entirely_viral_viromes_by_prophage.csv", row.names = FALSE)

# Calculates the start and end points of prophage location and number of genes
# Since the entire fragment is annotated as being viral, the first gene and first bp are the starts of all these sequences
Vir_viromes_by_prophage_locations <- Vir_viromes_by_prophage %>%
  mutate(Gene.Start = 1,
         Gene.End = Nb.genes,
         bp.Start = 1,
         bp.End = Sequence.Length..bp.,
         Fragment.Length = Sequence.Length..bp.
         )

# Writes the datatable to the tables folder
write.csv(Vir_viromes_by_prophage_locations, file = "../tables/entirely_viral_viromes_by_prophage_locations.csv", row.names = FALSE)
```

## Prophage detection by genome

RefSeq:


```r
# Collapses down VirSorter data to the genome level
Vir_refseq_by_genome <- Vir_refseq_by_prophage_locations %>%
  
  # Combines VirSorter data with genomic metadata
  inner_join(all_genome_data, by = "Genome.ID") %>%
  
  # Calculates the total size of 
  group_by(Genome.ID) %>%
  mutate(Total.Size.Prophages = sum(Fragment.Length)) %>%
  ungroup() %>%
  
  # Selects columns to keep in final table
  # Calculates some new columns as well
  transmute(Genome.ID = Genome.ID,
            NCBI.Taxon.ID = NCBI.Taxon.ID,
            Status = Status,
            Genome = Genome,
            Phylum = Phylum,
            Class = Class,
            Order = Order,
            Family = Family,
            Genus = Genus,
            Species = Species,
            Strain = Strain,
            Cultured = Cultured,
            Ecosystem = Ecosystem,
            Ecosystem.Category = Ecosystem.Category,
            Ecosystem.Type = Ecosystem.Type,
            Ecosystem.Subtype = Ecosystem.Subtype,
            Genome.Size = Genome.Size.....assembled,
            Gene.Count = Gene.Count.....assembled,
            Scaffold.Count = Scaffold.Count.....assembled,
            GOLD.Analysis.Project.Type = GOLD.Analysis.Project.Type,
            GOLD.Sequencing.Strategy = GOLD.Sequencing.Strategy,
            Prophages.per.genome = Prophages.per.genome,
            Total.Size.Prophages = Total.Size.Prophages,
            Genome.Size.Without.Prophage = Genome.Size - Total.Size.Prophages,
            Percentage.of.Genome.Prophage = 100 * (Total.Size.Prophages / Genome.Size)
            ) %>%
  
  # Only keeps a single row for each individual genome
  distinct(Genome.ID, .keep_all = TRUE)

# Writes the data to the tables folder
write.csv(Vir_refseq_by_genome, file = "../tables/entirely_viral_refseq_by_genome.csv", row.names = FALSE)
```

Viromes:


```r
# Collapses down VirSorter data to the genome level
Vir_viromes_by_genome <- Vir_viromes_by_prophage_locations %>%
  
  # Combines VirSorter data with genomic metadata
  inner_join(all_genome_data, by = "Genome.ID") %>%
  
  # Calculates the total size of 
  group_by(Genome.ID) %>%
  mutate(Total.Size.Prophages = sum(Fragment.Length)) %>%
  ungroup() %>%
  
  # Selects columns to keep in final table
  # Calculates some new columns as well
  transmute(Genome.ID = Genome.ID,
            NCBI.Taxon.ID = NCBI.Taxon.ID,
            Status = Status,
            Genome = Genome,
            Phylum = Phylum,
            Class = Class,
            Order = Order,
            Family = Family,
            Genus = Genus,
            Species = Species,
            Strain = Strain,
            Cultured = Cultured,
            Ecosystem = Ecosystem,
            Ecosystem.Category = Ecosystem.Category,
            Ecosystem.Type = Ecosystem.Type,
            Ecosystem.Subtype = Ecosystem.Subtype,
            Genome.Size = Genome.Size.....assembled,
            Gene.Count = Gene.Count.....assembled,
            Scaffold.Count = Scaffold.Count.....assembled,
            GOLD.Analysis.Project.Type = GOLD.Analysis.Project.Type,
            GOLD.Sequencing.Strategy = GOLD.Sequencing.Strategy,
            Prophages.per.genome = Prophages.per.genome,
            Total.Size.Prophages = Total.Size.Prophages,
            Genome.Size.Without.Prophage = Genome.Size - Total.Size.Prophages,
            Percentage.of.Genome.Prophage = 100 * (Total.Size.Prophages / Genome.Size)
            ) %>%
  
  # Only keeps a single row for each individual genome
  distinct(Genome.ID, .keep_all = TRUE)

# Writes the data to the tables folder
write.csv(Vir_viromes_by_genome, file = "../tables/entirely_viral_viromes_by_genome.csv", row.names = FALSE)
```
