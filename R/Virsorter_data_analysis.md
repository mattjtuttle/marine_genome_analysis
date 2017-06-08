# RefSeq data analysis
Matt Tuttle  
May 22, 2017  



## Data information

This document analyzes output data from a VirSorter analyses on marine bacterial genomes. 1778 genomes specified as belonging to a marine ecosystem type on JGI's IMG/M database were downloaded from JGI on May 15, 2017. These genomes comprise completed as well as permanent draft, draft, and single-cell amplified genomes. These genomes were then analyzed using VIRSorter version 1.0.3 in Cyverse's Discovery Environment. The analysis was performed twice, once using VirSorter's innate RefSeq database and a second time using the Viromes database.

One note to make regarding this dataset is that despite inclusion of whole as well as draft genomes, this dataset is apt to underestimate/underrepresent all of the bacterial diversity found within the ocean. This is due to the fact that not all genomes in the IMG/M database have complete metadata associated with the genome. For instance, our organism Sulfitobacter sp. CB2047 does not have an ecosystem type or habitat listed.

#### A note

A different dataset was originally run through VIRSorter that downloaded bacterial genomes based on habitats that stated that they were marine in nature. This dataset is not being used for further analysis due to the fact that subsequent viewing of the metadata for the genomes revealed that genomes from other aquatic environments were also somehow downloaded. These genomes were originally downloaded on April 28, 2017.

## Importing metadata and generating datatables

This code generates two datatables:
1. Metadata relating to genomes used in the Virsorter analyses
2. Metadata relating to scaffolds used in the Virsorter analyses


```r
# Imports genome metadata
all_genome_data <- read.csv("../data/genome_list.csv", header = TRUE) %>%
  rename(Genome.ID = IMG.Genome.ID)

# Only imports columns thought to be important for downstream analysis
select_genome_data <- all_genome_data %>%
    select(Taxon.oid = taxon_oid,
           Domain,
           Status,
           Genome.Name = Genome.Name...Sample.Name,
           Genome.ID,
           Phylum,
           Class,
           Order,
           Family,
           Genus,
           Species,
           NCBI.Taxon.ID,
           Strain, Cultured,
           Ecosystem,
           Ecosystem.Category,
           Ecosystem.Type,
           Ecosystem.Subtype,
           Genome.Size = Genome.Size.....assembled,
           Gene.Count = Gene.Count.....assembled,
           Scaffold.Count = Scaffold.Count.....assembled,
           GOLD.Analysis.Project.Type,
           GOLD.Sequencing.Strategy
           )
  
# Writes the curated genome metadata to the tables folder
write.csv(select_genome_data, file = "../tables/formatted_genome_list.csv", row.names = FALSE)

# Imports scaffold metadata
scaffold_data <- read.csv("../data/scaffold_list.csv", header = TRUE)
  
# Writes the curated scaffold metadata to the tables folder
write.csv(scaffold_data, file = "../tables/formatted_scaffold_list.csv", row.names = FALSE)
```

## Import of VirSorter output data


```r
# Imports VirSorter data output when using the RefSeq database
refseq_data <- read.csv("../data/Refseq_data.csv", header = TRUE)

# Imports VirSorter data output when using the Viromes database
viromes_data <- read.csv("../data/Viromes_data.csv", header = TRUE)
```

## Detection of prophages and calculating prophage genomic positions

For RefSeq data:


```r
# Makes changes to VirSorter output datatable to give more information about found prophages
refseq_by_prophage <- refseq_data %>%
 
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
Refseq_genomes_with_prophage <- distinct(refseq_by_prophage, Genome.ID) %>%
  nrow()
  
# Writes the datatable to the tables folder
write.csv(refseq_by_prophage, file = "../tables/refseq_by_prophage.csv", row.names = FALSE)

# Loads a datatable containing calculated locations of prophage fragments detected by VirSorter
# Makes Fragment column into a character vector so that it may be joined with other datatables
refseq_prophage_locations <- read.csv("../data/refseq_prophage_locations.csv", header = TRUE) %>%
  mutate(Fragment = as.character(Fragment))

# Creates a datatable that combines Refseq data at the prophage level with calculated prophage genomic locations
refseq_by_prophage_locations <- refseq_by_prophage %>%
  inner_join(refseq_prophage_locations, by = "Fragment")
  
# Some fragment lengths are negative in length due to how bp.Start/End were calculated for circular scaffolds
# Recalculates these negative value by taking into account the length of the contig that the prophage is on
refseq_by_prophage_locations$Fragment.Length <-ifelse(refseq_by_prophage_locations$Fragment.Length < 0, refseq_by_prophage_locations$Sequence.Length..bp. - refseq_by_prophage_locations$bp.Start + refseq_by_prophage_locations$bp.End, refseq_by_prophage_locations$Fragment.Length)
  
# Writes the combined datatable to the tables folder
write.csv(refseq_by_prophage_locations, file = "../tables/refseq_by_prophage_locations.csv", row.names = FALSE)
```

For Viromes data:


```r
# Makes changes to VirSorter output datatable to give more information about found prophages
viromes_by_prophage <- viromes_data %>%
 
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
Viromes_genomes_with_prophage <- distinct(viromes_by_prophage, Genome.ID) %>%
  nrow()
  
# Writes the datatable to the tables folder
write.csv(viromes_by_prophage, file = "../tables/viromes_by_prophage.csv", row.names = FALSE)

# Need to perform prophage location detection the same way as for the RefSeq data
#
#
#
#
#
#
#
```

## Prophage detection by genome

For RefSeq data:


```r
# Collapses down VirSorter data to the genome level
refseq_by_genome <- refseq_by_prophage_locations %>%
  
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
write.csv(refseq_by_genome, file = "../tables/refseq_by_genome.csv", row.names = FALSE)
```

For Viromes data:




## Prophage detection by species



## Plots distribution of taxa

