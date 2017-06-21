# This file looks at all of the single cell genomes for which there is percent coverage data avaailable

library(tidyverse)

percent_coverage <- all_genome_data %>%
  select(Genome.Name = Genome.Name...Sample.Name,
         JGI.Analysis.Project.Type,
         Genome.Completeness = Genome.Completeness..) %>%
  filter(JGI.Analysis.Project.Type != "Genome Analysis") %>%
  arrange(Genome.Completeness)

write.csv(percent_coverage, file = "../tables/percent_coverage.csv", row.names = FALSE)
