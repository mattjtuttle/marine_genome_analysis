# Data summary tables by phyla
Matt Tuttle  
July 11, 2017  



## Script information

The following scripts are for creating a data summary table that gives an overview of the overall dataset for the main body of the review paper.

## Data import


```r
# Imports the entire queried dataset
all_genome_data <- read.csv("../data/genome_list.csv", header = TRUE) %>%
  
  # Changes headers for some columns
  rename(Genome.ID = IMG.Genome.ID,
         Genome.Size = Genome.Size.....assembled,
         Scaffold.Count = Scaffold.Count.....assembled) %>%
  
  # Removes all SAGs that were not contamination screened
  filter(GOLD.Analysis.Project.Type != "Single Cell Analysis (unscreened)")
  
# Changes names of phyla to remove "Candidatus"
all_genome_data$Phylum[all_genome_data$Phylum %in% "Candidatus Marinimicrobia"] <- "Marinimicrobia"

# Imports various VirSorter outputs
all_refseq_by_genome_cat13 <- read.csv("../tables/all_refseq_by_genome.csv", header = TRUE) %>%
  filter(GOLD.Analysis.Project.Type != "Single Cell Analysis (unscreened)")

all_refseq_by_genome_cat12 <- read.csv("../tables/all_refseq_by_genome_cat12.csv", header = TRUE) %>%
  filter(GOLD.Analysis.Project.Type != "Single Cell Analysis (unscreened)")

all_viromes_by_genome_cat13 <- read.csv("../tables/all_viromes_by_genome.csv", header = TRUE) %>%
  filter(GOLD.Analysis.Project.Type != "Single Cell Analysis (unscreened)")

all_viromes_by_genome_cat12 <- read.csv("../tables/all_viromes_by_genome_cat12.csv", header = TRUE) %>%
  filter(GOLD.Analysis.Project.Type != "Single Cell Analysis (unscreened)")
```

## A function for the creation of summary tables


```r
# Calculates phyla that have at least ten genomes in the overall dataset
# Values calculated here used to create filter in below function
phyla_with_ten <- all_genome_data %>%
  group_by(Phylum) %>%
  filter(n() >= 10) %>%
  ungroup() %>%
  distinct(Phylum) %>%
  filter(Phylum != "unclassified") # Removes genomes where the Phylum is unclassified

# Creates a function that is able to make summary tables for VirSorter data at the genome level
create_summary_table <- function(df){
  
  table <- df %>%
    group_by(Phylum) %>%
    summarise(
      Number.of.genomes = n(),
      Percent.Complete.Genomes = 100 * (sum(Status == "Finished") / n()),
      Percent.Cultured.Genomes = 100 * (sum(Cultured == "Yes") / n()),
      SAGs.1 = sum(GOLD.Analysis.Project.Type == "Single Cell Analysis (screened)"),
      SAGs.2 = sum(GOLD.Analysis.Project.Type == "Single Cell Analysis (unscreened)"),
      Percent.SAGs = 100 * ( (SAGs.1 + SAGs.2) / n()),
      Percent.from.metagenomes = 100 * (sum(GOLD.Analysis.Project.Type == "Genome from Metagenome") / n()),
      Average.genome.size = mean(Genome.Size),
      Median.genome.size = median(Genome.Size),
      Min.genome.size = min(Genome.Size),
      Max.genome.size = max(Genome.Size),
      Average.Num.Contigs = mean(Scaffold.Count),
      Median.Num.Contigs = median(Scaffold.Count),
      Min.Num.Contigs = min(Scaffold.Count),
      Max.Num.Contigs = max(Scaffold.Count)
              ) %>%
    
    # Gets rid of useless columns used for calculations
    select(-SAGs.1, -SAGs.2) %>%
    
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
  return(table)
}
```

## Creation of summary datatables by phyla

For the overall dataset queried:


```r
# Creates table
queried_genomes_summary <- create_summary_table(all_genome_data)

# Writes it into a summaries subfolder of the tables folder
write.csv(queried_genomes_summary, file = "../tables/summaries/queried_genomes.csv", row.names = FALSE)
```


All of the following tables being created below represent a breakdown of these stats in the datasets of detected prophages for comparison with the above summary table of the overall dataset.

For RefSeq by genome (cat 1-3):


```r
# Creates table
all_refseq_by_genome_cat13_summary <- create_summary_table(all_refseq_by_genome_cat13)

# Writes it into a summaries subfolder of the tables folder
write.csv(all_refseq_by_genome_cat13_summary, file = "../tables/summaries/all_refseq_by_genome_summary.csv", row.names = FALSE)
```

For RefSeq by genome (cat 1-2):


```r
# Creates table
all_refseq_by_genome_cat12_summary <- create_summary_table(all_refseq_by_genome_cat12)

# Writes it into a summaries subfolder of the tables folder
write.csv(all_refseq_by_genome_cat12_summary, file = "../tables/summaries/all_refseq_by_genome_cat12_summary.csv", row.names = FALSE)
```

For Viromes by genome (cat 1-3):


```r
# Creates table
all_viromes_by_genome_cat13_summary <- create_summary_table(all_viromes_by_genome_cat13)

# Writes it into a summaries subfolder of the tables folder
write.csv(all_viromes_by_genome_cat13_summary, file = "../tables/summaries/all_viromes_by_genome_summary.csv", row.names = FALSE)
```

For Viromes by genome (cat 1-2):


```r
# Creates table
all_viromes_by_genome_cat12_summary <- create_summary_table(all_viromes_by_genome_cat12)

# Writes it into a summaries subfolder of the tables folder
write.csv(all_viromes_by_genome_cat12_summary, file = "../tables/summaries/all_viromes_by_genome_cat12_summary.csv", row.names = FALSE)
```


## Calculating statistics on the output datasets

The following code calculates some numbers that pertain to prophage detection of the overall output dataset:


```r
# Creates a function that calculates some basic overall stats on the VirSorter output datasets
# Only takes into account the phyla for which there are at least ten representatives in the overall dataset
calc_stats <- function(df){
  
  table <- df %>%
    filter(Phylum %in% c("Firmicutes",
                         "Proteobacteria",
                         "Actinobacteria",
                         "Bacteroidetes",
                         "Cyanobacteria",
                         "Chloroflexi",
                         "Planctomycetes",
                         "Verrucomicrobia",
                         "Deferribacteres",
                         "Gemmatimonadetes",
                         "Marinimicrobia",
                         "Thermotogae",
                         "Nitrospinae",
                         "Nitrospirae")
           ) %>%
    summarise(# Calulates the percent lysogens detected that were from cultivated isolates
              Percent.Cultured = 100 * (sum(Cultured == "Yes") / n()),
              
              # Calculates the percent of complete (finished) genomes containing prophage
              Percent.Compelete = 100 * (sum(Status == "Finished") / n())
              )
  
  return(table)
}
```

For RefSeq:


```r
# Calculates percent of genomes of isolated organisms where at least one prophage was detected
refseq_13_overall <- calc_stats(all_refseq_by_genome_cat13)
```

For Viromes:


```r
# Calculates percent of genomes of isolated organisms where at least one prophage was detected
viromes_13_overall <- calc_stats(all_viromes_by_genome_cat13)
```

