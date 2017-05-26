# Looks for SAGs that appear to have whole contigs as phage sequence
# Any contig determined to be 80% or more a prophage region is considered to be viral only sequence according to VirSOrter

library(tidyverse)

# Imports scaffold data
scaffold_data <- read.csv("./data/scaffold_list.csv", header = TRUE)

# Imports relavent genomic metadata
genomic_data <- read.csv("./data/genome_list.csv", header = TRUE) %>%
  select(Taxon.oid = taxon_oid,
         Domain,
         Status,
         Genome.Name = Genome.Name...Sample.Name,
         Genome.ID = IMG.Genome.ID,
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
         GOLD.Sequencing.Strategy
  )

# Imports VirSorter data
viral_seqs <- read.csv("./data/Refseq_phage_contigs_cat1-3.csv", header = TRUE) %>%
  mutate(Contig_id = gsub("VIRSorter_", "", Contig_id)) %>%
  mutate(Scaffold.ID = gsub("\\_.*", "", Contig_id)) %>%
  mutate(Scaffold.ID = as.numeric(Scaffold.ID)) %>%
  mutate(Fragment = gsub("VIRSorter_", "", Fragment)) %>%
  mutate(Nb.phage.hallmark.genes = ifelse(is.na(Nb.phage.hallmark.genes), 0, Nb.phage.hallmark.genes)) %>%
  inner_join(scaffold_data, by = "Scaffold.ID") %>%
  inner_join(genomic_data, by = "Genome.ID") %>%
  select(Genome,
         GOLD.Sequencing.Strategy) %>%
  View()