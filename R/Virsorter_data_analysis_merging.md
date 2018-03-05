# Merging VirSorter output data
Matt Tuttle  
July 11, 2017  



## Script information

The following scripts are for merging the category 1-3 (entirely viral) and category 4-6 (detected prophage) outputs of the RefSeq and Viromes VirSorter analyses. The two separate categories are being combined since our dataset contains SAGs and highly fragmented genome sequences which have the potential to hit as being entirely viral as opposed to prophages when using VirSorter.

## VirSorter output by prophage

RefSeq:


```r
# Imports entirely viral and prophage data
vir_refseq_by_prophage <- read.csv("../tables/entirely_viral_refseq_by_prophage.csv", header = TRUE)
pro_refseq_by_prophage <- read.csv("../tables/refseq_by_prophage.csv", header = TRUE)

# Joins the two dataset together
all_refseq_by_prophage <- bind_rows(pro_refseq_by_prophage, vir_refseq_by_prophage) %>%
  arrange(Category, Contig_id)
```

```
## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character
```

```r
# Writes the datatable to the tables folder
write.csv(all_refseq_by_prophage, file = "../tables/all_refseq_by_prophage.csv", row.names = FALSE)
```

Viromes:


```r
# Imports entirely viral and prophage data
vir_viromes_by_prophage <- read.csv("../tables/entirely_viral_viromes_by_prophage.csv", header = TRUE)
pro_viromes_by_prophage <- read.csv("../tables/viromes_by_prophage.csv", header = TRUE)

# Joins the two dataset together
all_viromes_by_prophage <- bind_rows(pro_viromes_by_prophage, vir_viromes_by_prophage) %>%
  arrange(Category, Contig_id)
```

```
## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character
```

```r
# Writes the datatable to the tables folder
write.csv(all_viromes_by_prophage, file = "../tables/all_viromes_by_prophage.csv", row.names = FALSE)
```


## VirSorter output by genome

RefSeq:


```r
# Imports entirely viral and prophage data
vir_refseq_by_genome <- read.csv("../tables/entirely_viral_refseq_by_genome.csv", header = TRUE)
pro_refseq_by_genome <- read.csv("../tables/refseq_by_genome.csv", header = TRUE)

# Joins the two dataset together
all_refseq_by_genome <- bind_rows(pro_refseq_by_genome, vir_refseq_by_genome) %>%
  arrange(Genome)
```

```
## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character
```

```r
# Writes the datatable to the tables folder
write.csv(all_refseq_by_genome, file = "../tables/all_refseq_by_genome.csv", row.names = FALSE)
```

Viromes:


```r
# Imports entirely viral and prophage data
vir_viromes_by_genome <- read.csv("../tables/entirely_viral_viromes_by_genome.csv", header = TRUE)
pro_viromes_by_genome <- read.csv("../tables/viromes_by_genome.csv", header = TRUE)

# Joins the two dataset together
all_viromes_by_genome <- bind_rows(pro_viromes_by_genome, vir_viromes_by_genome) %>%
  arrange(Genome)
```

```
## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character

## Warning in bind_rows_(x, .id): Unequal factor levels: coercing to character
```

```r
# Writes the datatable to the tables folder
write.csv(all_viromes_by_genome, file = "../tables/all_viromes_by_genome.csv", row.names = FALSE)
```

## Looking at categories of detected prophages

In analyzing the VirSorter output data, we would like to see how the phylogenetic distribution of prophages changes if only prophages detected with a confidence (category) of 1 and 2 are used as opposed to including 1, 2, and 3. According the VirSorter paper (Roux, 2015), if only categories 1+2 are looked at, the program is more precise. Mean that there is a greater number of correct predictions than if all categories (1+2+3) are included. However, inclusion of the prophages detected as category 3 will increase recall, by predicting more prophage regions at the expense of a highe rate of false positives.


```r
# Selects cat 1+2 data for refseq and viromes datasets
all_refseq_by_prophage_cat12 <- all_refseq_by_prophage %>%
  filter(Category != 3)

all_viromes_by_prophage_cat12 <- all_viromes_by_prophage %>%
  filter(Category != 3)


# Writes the datatables to the tables folder
write.csv(all_refseq_by_prophage_cat12, file = "../tables/all_refseq_by_prophage_cat12.csv", row.names = FALSE)

write.csv(all_viromes_by_prophage_cat12, file = "../tables/all_viromes_by_prophage_cat12.csv", row.names = FALSE)
```

Compressing cat 1+2 prophage data to the genome level for RefSeq data:


```r
# Loads genomic level data for compressing prphage level data
all_genome_data <- read.csv("../data/genome_list.csv", header = TRUE) %>%
  rename(Genome.ID = IMG.Genome.ID)


# Collapses down VirSorter data to the genome level
all_refseq_by_genome_cat12 <- all_refseq_by_prophage_cat12 %>%
  
  # Combines VirSorter data with genomic metadata
  inner_join(all_genome_data, by = "Genome.ID") %>%
  
  # Selects columns to keep in final table
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
            Prophages.per.genome = Prophages.per.genome
            ) %>%
  
  # Only keeps a single row for each individual genome
  distinct(Genome.ID, .keep_all = TRUE)

# Writes the data to the tables folder
write.csv(all_refseq_by_genome_cat12, file = "../tables/all_refseq_by_genome_cat12.csv", row.names = FALSE)
```

and Viromes data:


```r
# Collapses down VirSorter data to the genome level
all_viromes_by_genome_cat12 <- all_viromes_by_prophage_cat12 %>%
  
  # Combines VirSorter data with genomic metadata
  inner_join(all_genome_data, by = "Genome.ID") %>%
  
  # Selects columns to keep in final table
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
            Prophages.per.genome = Prophages.per.genome
            ) %>%
  
  # Only keeps a single row for each individual genome
  distinct(Genome.ID, .keep_all = TRUE)

# Writes the data to the tables folder
write.csv(all_viromes_by_genome_cat12, file = "../tables/all_viromes_by_genome_cat12.csv", row.names = FALSE)
```

