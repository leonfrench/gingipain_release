library(readxl)
library(reshape2)
library(readr)
library(here)
library(dplyr)
library(tidyr)
library(broom)
library(magrittr)
library(jsonlite)
library(RCurl)
library(ggplot2)
library(GEOquery)
detach("package:XML", unload = TRUE) #handles count()/n() conflict
detach("package:dplyr", unload = TRUE)
library(dplyr)

#data from "mRNA-Seq expression and MS3 proteomics profiling of human post-mortem BA9 brain tissue for Parkinson Disease and neurologically normal individuals"
GSE_from_GEO <- getGEO("GSE68719")

#SRP058181
geo_sample_table <- GSE_from_GEO$GSE68719_series_matrix.txt.gz@phenoData@data

geo_sample_table %<>% as_tibble()
geo_sample_table[1,] %>% as.data.frame()
geo_sample_table %<>% select(geo_accession, sex=`gender:ch1`, title, pmi = `post-mortem interval (pmi):ch1`, has_proteomics = `proteomics study:ch1`, age = `age at death:ch1`, rin=`rna integrity number (rin):ch1`)
geo_sample_table %>% select(title) %>% arrange(title) %>% .$title

sample_table <- read_csv(here("data", "GSE68719", paste0("SraRunTable.","GSE68719",".txt")))

colnames(sample_table)
sample_table[1,] %>% as.data.frame()
sample_table %<>% select(Run, geo_accession = `Sample Name`, Bases, Bytes)
sample_table <- inner_join(sample_table, geo_sample_table)

sample_table %<>% mutate(group = substring(title,1,1))
sample_table %<>% mutate(age = as.numeric(age))
sample_table %<>% mutate(rin = as.numeric(rin))
sample_table %<>% mutate(pmi = as.numeric(pmi))
sample_table %<>% mutate(short_title = gsub(" .*","", title))

cause_of_death_info <- read_xlsx(here("data", "GSE68719", "12920_2016_164_MOESM1_ESM_Subject_info.xlsx"))
colnames(cause_of_death_info)
cause_of_death_info %<>% select(short_title = `RNA-Seq Samples`, cause_of_death = `Cause of death`, Dementia, ApoE)
sample_table <- inner_join(sample_table, cause_of_death_info)

all_run_tax_tables <- NULL
#grab the data from sequence read archive
for (target_run in sample_table$Run) {
  print(target_run)
  x <- getURL(paste0("https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=", target_run))
  
  tax_table_text <- gsub(".*var oTaxAnalysisData = ", "", x)
  tax_table_text <- gsub(";.*", "", tax_table_text)
  tax_table_text <- gsub(", \n0\\]", "\\]", tax_table_text)
  tax_table <- fromJSON(txt=tax_table_text)
  tax_table <- flatten(tax_table, recursive = TRUE)
  tax_table %<>% as_tibble()
  print(paste0("Rows: ", dim(tax_table)[1]))
  #tax_table %>% arrange(-as.numeric(d.kbp)) %>% head(30) %>% as.data.frame
  tax_table %<>% mutate(Run = target_run)
  all_run_tax_tables <- bind_rows(all_run_tax_tables, tax_table)
}
all_run_tax_tables %<>% mutate(d.percent = as.numeric(d.percent))
all_run_tax_tables %<>% mutate(d.kbp = as.numeric(d.kbp))

all_run_tax_tables %>% select(d.name, Run) %>% count()
all_run_tax_tables %>% select(d.name, Run) %>% distinct() %>% count()

#complete dataset with zeroes
all_run_tax_tables %<>% select(-n,-p,-c) %>% tidyr::complete(Run, d.name, fill = list(d.percent = 0, d.kbp = 0))

all_run_tax_tables <- inner_join(all_run_tax_tables, sample_table)

all_run_tax_tables %>% filter(grepl("Porphyromonas", d.name), d.kbp > 0) %>% select(Run, d.name, group, rin, Bytes, Bases, d.kbp, sex,title) %>% distinct()
all_run_tax_tables %>% filter(grepl("Porphyromonas", d.name), d.kbp > 0) %>% select(d.name) %>% distinct()
all_run_tax_tables %>% filter(grepl("Porphyromonas endodontalis", d.name), d.kbp > 0) %>% select(Run, d.name, group, rin, Bytes, Bases, d.kbp, sex,title) %>% distinct()
all_run_tax_tables %>% select(d.name, Run) %>% count()
all_run_tax_tables %>% select(d.name, Run) %>% distinct() %>% count()

all_run_tax_tables %>% filter("Porphyromonas gingivalis" == d.name, d.kbp > 0, group == "C") %>% select(Run, group, rin, Bytes, Bases, d.kbp, sex, Dementia, cause_of_death, ApoE) %>% distinct()
all_run_tax_tables %>% filter("Porphyromonas gingivalis" == d.name, d.kbp == 0, group == "C") %>% select(Run, group, rin, Bytes, Bases, d.kbp, sex, Dementia, cause_of_death, ApoE) %>% distinct() %>% as.data.frame()

just_controls_for_demographics <- all_run_tax_tables %>% filter("Porphyromonas gingivalis" == d.name, group == "C") %>% select(Run, group, rin, Bytes, Bases, d.kbp, sex, Dementia, cause_of_death, ApoE, age, rin, pmi) %>% distinct()
just_controls_for_demographics %>% .$Run %>% unique() %>% length #44 controls
just_controls_for_demographics %>% .$Run %>% length #44 controls
just_controls_for_demographics %>% group_by(Run) %>% count() %>% filter(n > 1)

just_controls_for_demographics %<>% mutate(pg_detected = d.kbp > 0) 
just_controls_for_demographics %>% group_by(pg_detected) %>% count()
just_controls_for_demographics %>% group_by(pg_detected, cause_of_death) %>% count() %>% as.data.frame()
just_controls_for_demographics %>% filter(cause_of_death == "Pancreatic cancer")
just_controls_for_demographics %>% filter(pg_detected)
just_controls_for_demographics$Run %>% unique %>% length
just_controls_for_demographics %>% filter(ApoE != "N/A") %>% group_by(pg_detected, ApoE) %>% count()
range(just_controls_for_demographics$age)

#no hits for any red complex bacteria
all_run_tax_tables %>% filter(grepl("orsythia", d.name)) %>% select(Run, group, rin, Bytes, Bases, d.kbp, sex) %>% distinct()
all_run_tax_tables %>% filter(grepl("orsythus", d.name)) %>% select(Run, group, rin, Bytes, Bases, d.kbp, sex) %>% distinct()
all_run_tax_tables %>% filter(grepl("enticola", d.name)) %>% select(Run, group, rin, Bytes, Bases, d.kbp, sex) %>% distinct()
all_run_tax_tables %>% filter(grepl("alocis", d.name)) %>% select(Run, group, rin, Bytes, Bases, d.kbp, sex) %>% distinct()
all_run_tax_tables %>% filter(grepl("ilifactor", d.name)) %>% select(Run, group, rin, Bytes, Bases, d.kbp, sex) %>% distinct()

#forget about all other taxa
pg_detection_data <- all_run_tax_tables %>% filter("Porphyromonas gingivalis" == d.name)
pg_detection_data %<>% mutate(pg_detected = d.kbp > 0) 

pg_detection_data %>% group_by(group) %>% count()
pg_detection_data %>% group_by(group, pg_detected) %>% count()

#test for differences in rin, pmi, age
wilcox.test(rin ~ pg_detected, pg_detection_data)
wilcox.test(pmi ~ pg_detected, pg_detection_data)
wilcox.test(age ~ pg_detected, pg_detection_data)
wilcox.test(Bytes ~ pg_detected, pg_detection_data)
wilcox.test(Bases ~ pg_detected, pg_detection_data)

#in just controls - again strong relationship with bases
wilcox.test(rin ~ pg_detected, pg_detection_data %>% filter(group == "C"))
wilcox.test(pmi ~ pg_detected, pg_detection_data %>% filter(group == "C"))
wilcox.test(age ~ pg_detected, pg_detection_data %>% filter(group == "C"))
wilcox.test(Bytes ~ pg_detected, pg_detection_data %>% filter(group == "C"))
wilcox.test(Bases ~ pg_detected, pg_detection_data %>% filter(group == "C"))
t.test(Bases ~ pg_detected, pg_detection_data %>% filter(group == "C"))

ggplot(data= pg_detection_data, aes(y=Bases, x = pg_detected)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .3) + facet_wrap(.~group) + theme_bw()

write_csv(pg_detection_data, here("data", "GSE68719", "pg_detection_data.csv"))

