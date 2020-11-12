library(readxl)
library(reshape2)
library(readr)
library(here)
library(dplyr)
library(purrr)
library(tidyr)
library(broom)
library(magrittr)
detach("package:XML", unload = TRUE) #handles count()/n() conflict
detach("package:dplyr", unload = TRUE)
library(dplyr)

#read in sample meta data
pg_detection_data <- read_csv(here("data", "GSE68719", "pg_detection_data.csv"))

#perform differential expression          
expression <- read_tsv(here("data", "GSE68719", "GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz"))
dim(expression)
expression %<>% select(-EnsemblID)
expression %<>% melt(value.name = "expression", variable.name = "short_title") %>% as_tibble()
expression %<>% mutate(expression = log(expression+1))
expression %<>% rename(gene_symbol = symbol)

pg_detection_data_expression <- inner_join(pg_detection_data, expression)
pg_detection_data_expression %<>% group_by(gene_symbol, group) %>% do(fit_expression = lm(expression ~ pg_detected + Bases + age + pmi + rin, data = .))

#pg_detection_data_expression_fits = tidy(pg_detection_data_expression, fit_expression) #old code, depends on R libraries
#from https://stackoverflow.com/questions/62972843/retreiving-tidy-results-from-regression-by-group-with-broom
pg_detection_data_expression_fits <- pg_detection_data_expression %>% ungroup %>% 
  transmute(group, gene_symbol, HourCoef = map(fit_expression, tidy)) %>% 
  unnest(HourCoef)  #runs slow (less than 5 mins)

per_gene_pg_results <- pg_detection_data_expression_fits %>% filter(term == "pg_detectedTRUE") %>% 
  mutate(signed_log_p = -1 * sign(estimate) * log(p.value)) %>% group_by(group) %>%
  mutate(rank = rank(signed_log_p)) %>% arrange(signed_log_p) 

#filter for the genes used in the amino acid proportion tests
gene_universe <- read_csv(here("data", "gene_lists", "gene_universe.txt"), col_names = F)
per_gene_pg_results %<>% filter(gene_symbol %in% gene_universe$X1)

per_gene_pg_results %<>% group_by(group) %>% mutate(p.adjust = p.adjust(p.value, method = "fdr"))
per_gene_pg_results %>% filter(p.adjust < 0.05) %>% arrange(p.value)
per_gene_pg_results %<>% arrange(-signed_log_p)

per_gene_pg_results %>% filter(group == "C") %>% select(-term) %>% write_csv(here("results", "GSE68719", "genewise_pg_pd_controls.csv"))
per_gene_pg_results %>% filter(group == "P") %>% select(-term) %>% write_csv(here("results", "GSE68719", "genewise_pg_pd_cases.csv"))
