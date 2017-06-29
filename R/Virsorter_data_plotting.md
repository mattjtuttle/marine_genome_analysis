# RefSeq data plotting
Matt Tuttle  
`r Sys.Date()`  



## Data information

The following code imports processed data from the tables folder for using in making plots to visualize the data. To see how the data was processed, see `Virsorter_data_analysis.Rmd`.


```r
# Imports RefSeq data at the genome level from the tables folder
refseq_by_genome <- read.csv("../tables/refseq_by_genome.csv", header = TRUE)

# Imports Viromes data at the genome level from the tables folder
viromes_by_genome <- read.csv("../tables/viromes_by_genome.csv", header = TRUE)

# Imports genomic metadata
genome_metadata <- read.csv("../data/genome_list.csv", header = TRUE) %>%
  mutate(Genome.ID = IMG.Genome.ID)

# Imports RefSeq data at the prophage level
refseq_by_prophage <- read.csv("../tables/refseq_by_prophage.csv", header = TRUE)

# Imports Viromes data at the prophage level
viromes_by_prophage <- read.csv("../tables/viromes_by_prophage.csv", header = TRUE)
```

## Phylogenetic distribution of IMG dataset of marine bacterial genomes


```r
# 
metadata_by_phyla <- genome_metadata %>%
  count(Phylum, sort = TRUE) %>%
  transmute(Phylum,
            Dataset.Genomes = n)
```


## Phylogenetic distribution at the genome level

This plot compares the phylogenetic distribution of the RefSeq and Viromes datasets at the genome level. It creates a table which is saved to the tables folder as well as creates a plot to visualize differences between the two datasets.


```r
# Selects Genome IDs from all of the genomes for which prophages were predicted
refseq_IDs <- refseq_by_genome %>%
  select(Genome.ID)

viromes_IDs <- viromes_by_genome %>%
  select(Genome.ID)

# Adds selects metadata from genomes that are predicted to contain prophages
refseq_genome_distribution <- inner_join(refseq_IDs, genome_metadata, by = "Genome.ID") %>%
  select(Genome.Name = Genome.Name...Sample.Name,
         IMG.Genome.ID,
         Phylum,
         Class,
         Order,
         Family,
         Genus,
         Species,
         Strain)

viromes_genome_distribution <- inner_join(viromes_IDs, genome_metadata, by = "Genome.ID") %>%
  select(Genome.Name = Genome.Name...Sample.Name,
         IMG.Genome.ID,
         Phylum,
         Class,
         Order,
         Family,
         Genus,
         Species,
         Strain)

# Creates a table of Refseq data by phylum
refseq_genomes_by_phyla <- refseq_genome_distribution %>%
  count(Phylum, sort = TRUE) %>%
  transmute(Phylum,
            RefSeq.Genomes = n)

# Creates a table of Viromes data by phylum
viromes_genomes_by_phyla <- viromes_genome_distribution %>%
  count(Phylum, sort = TRUE) %>%
  transmute(Phylum,
            Viromes.Genomes = n)

# Combines the RefSeq and Viromes phyla tables with dataset metadata
genomes_by_phyla <- full_join(refseq_genomes_by_phyla, viromes_genomes_by_phyla, by = "Phylum") %>%
  full_join(metadata_by_phyla, by = "Phylum")
  
# Makes all NA values zeros
genomes_by_phyla[is.na(genomes_by_phyla)] <- 0

# Calculates percentage of bacterial phylum containing prophages
genomes_by_phyla <- genomes_by_phyla %>%
  mutate(RefSeq.Percent = 100 * (RefSeq.Genomes / Dataset.Genomes)) %>%
  mutate(Viromes.Percent = 100 * (Viromes.Genomes / Dataset.Genomes))


# Writes the table to the tables folder
write.csv(genomes_by_phyla, file = "../tables/genomes_by_phyla.csv", row.names = FALSE)


# Creates a plot of Refseq data by phylum
# Creates a plot of Viromes data by phylum
```

## Proteobacterial class distribution at the prophage level


```r
# Calculates number of genomes in each proteobacterial class for the overall dataset
proteo_metadata_by_class <- genome_metadata %>%
  filter(Phylum == "Proteobacteria") %>%
  count(Class, sort = TRUE) %>%
  transmute(Class,
            Dataset.Genomes = n)

# Calculates number of genomes in each proteobacterial class detected to contain a prophage via the refseq database
# Uses genome distributions from datasets in calculating
refseq_proteo_genomes_by_class <- refseq_genome_distribution %>%
  filter(Phylum == "Proteobacteria") %>%
  count(Class, sort = TRUE) %>%
  transmute(Class,
            RefSeq.Genomes = n)

# Calculates number of genomes in each proteobacterial class detected to contain a prophage via the viromes database
# Uses genome distributions from datasets in calculating
viromes_proteo_genomes_by_class <- viromes_genome_distribution %>%
  filter(Phylum == "Proteobacteria") %>%
  count(Class, sort = TRUE) %>%
  transmute(Class,
            Viromes.Genomes = n)

# Combines the RefSeq and Viromes prteobacterial class tables with dataset metadata
proteo_genomes_by_class <- full_join(refseq_proteo_genomes_by_class, viromes_proteo_genomes_by_class, by = "Class") %>%
  full_join(proteo_metadata_by_class, by = "Class")
  
# Makes all NA values zeros
proteo_genomes_by_class[is.na(proteo_genomes_by_class)] <- 0

# Calculates the percentage of proteobacterial classes containing prophages
proteo_genomes_by_class <- proteo_genomes_by_class %>%
  mutate(RefSeq.Percent = 100 * (RefSeq.Genomes / Dataset.Genomes)) %>%
  mutate(Viromes.Percent = 100 * (Viromes.Genomes / Dataset.Genomes))


# Writes the table to the tables folder
write.csv(proteo_genomes_by_class, file = "../tables/proteo_genomes_by_class.csv", row.names = FALSE)
```


## Phylogenetic distribution at the prophage level


```r
# Selects Genome IDs from each prophage that Virsorter predicted
prophages_refseq_IDs <- refseq_by_prophage %>%
  select(Genome.ID)

prophages_viromes_IDs <- viromes_by_prophage %>%
  select(Genome.ID)

# Adds selects metadata from genomes that are predicted to contain prophages
refseq_prophages_distribution <- inner_join(prophages_refseq_IDs, genome_metadata, by = "Genome.ID") %>%
  select(Genome.Name = Genome.Name...Sample.Name,
         IMG.Genome.ID,
         Phylum,
         Class,
         Order,
         Family,
         Genus,
         Species,
         Strain)

viromes_prophages_distribution <- inner_join(prophages_viromes_IDs, genome_metadata, by = "Genome.ID") %>%
  select(Genome.Name = Genome.Name...Sample.Name,
         IMG.Genome.ID,
         Phylum,
         Class,
         Order,
         Family,
         Genus,
         Species,
         Strain)

# Creates a table of Refseq data by phylum
refseq_prophages_by_phyla <- refseq_prophages_distribution %>%
  count(Phylum, sort = TRUE) %>%
  transmute(Phylum,
            RefSeq.Prophages = n)

# Creates a table of Viromes data by phylum
viromes_prophages_by_phyla <- viromes_prophages_distribution %>%
  count(Phylum, sort = TRUE) %>%
  transmute(Phylum,
            Viromes.Prophages = n)

# Combines the RefSeq and Viromes phyla tables
prophages_by_phyla <- full_join(refseq_prophages_by_phyla, viromes_prophages_by_phyla, by = "Phylum")

write.csv(prophages_by_phyla, file = "../tables/prophages_by_phyla.csv", row.names = FALSE)


# Creates a plot of Refseq data by phylum
# Creates a plot of Viromes data by phylum
```

