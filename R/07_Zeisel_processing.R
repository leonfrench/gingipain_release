library(homologene)
library(readr)
library(dplyr)
library(tidyr)
library(reshape2)
library(magrittr)
library(here)

Zeisel_file <- here("data", "Zeisel_etal_2018","l5_all.agg.tab")
if(!file.exists(Zeisel_file)){
  dir.create(here("data", "Zeisel_etal_2018"), recursive = T)
  download.file("https://storage.googleapis.com/linnarsson-lab-loom/l5_all.agg.tab", Zeisel_file)
}

table <- read_tsv(Zeisel_file, col_names = FALSE)
descriptions <- table[29, 9:273]

table %<>% filter(!is.na(X1) | X8 == 'ClusterName')
colnames(table) <- c(unlist(table[2,1:7]), unlist(table[1, 8:ncol(table)]))
table <- table[-c(1,2),]

# columns to drop
table <- table[!is.na(names(table))]
dropcols <- c('Accession', '_LogCV', '_LogMean', '_Selected', '_Total', '_Valid', 'ClusterName')
table %<>% select(-one_of(dropcols))

#to change all but the 'Gene' column name into numeric type
table[,2:length(colnames(table))] %<>% lapply(function(x) as.numeric(as.character(x)))

cell_type_info <- rbind(colnames(table[2:266]), descriptions)
cell_type_info <- t(cell_type_info)
colnames(cell_type_info) <- c('cluster_id', 'description')
cell_type_info <- as_tibble(cell_type_info)
#now melt
linnarsson <- as_tibble(reshape2::melt(table, id.vars='Gene', variable.name='cluster_id', value.name='expression'))
linnarsson %<>% mutate(log1Expression=log(1+expression))

# save tables for gene expression of Chat and AVP alone
linnarsson %>% filter(Gene == 'Chat') %>% write_csv(here('results/mouse', 'Zeisel_Chat_expression.csv'))
linnarsson %>% filter(Gene == 'Avp') %>% write_csv(here('results/mouse', 'Zeisel_Avp_expression.csv'))
linnarsson %>% filter(Gene == 'Vip') %>% write_csv(here('results/mouse', 'Zeisel_Vip_expression.csv'))

# zscore across genes
linnarsson %<>% 
  group_by(Gene) %>% 
  mutate(log1ExpZ = (log1Expression - mean(log1Expression)) / sd(log1Expression))

print(dim(linnarsson))
linnarsson %<>% filter(!is.na(log1ExpZ))       
print(dim(linnarsson))
linnarsson %<>% select(-expression, -log1Expression)

#some genes are duplicated - just average them
linnarsson %<>% group_by(Gene, cluster_id) %>% summarize(log1ExpZ = mean(log1ExpZ))
linnarsson %<>% group_by(cluster_id)
linnarsson %<>% mutate(log1ExpZRank = rank(log1ExpZ)) %>% select(-log1ExpZ)

linnarssonMatrixMouse <- dcast(linnarsson, Gene ~ cluster_id, value.var = "log1ExpZRank")

rownames(linnarssonMatrixMouse) <- linnarssonMatrixMouse$Gene
linnarssonMatrixMouse$Gene <- NULL

#repeat code again after filtering for mouse genes that can be reached via human symbols (should be refactored)
unique_genes_all <- rownames(linnarssonMatrixMouse)
unique_genes_human_reachable <- mouse2human(unique_genes_all)$mouseGene
linnarsson %<>% filter(Gene %in% unique_genes_human_reachable) 
linnarsson %<>% mutate(log1ExpZRank = rank(log1ExpZRank)) #rerank after filtering

linnarssonMatrixHumanReachable <- dcast(linnarsson, Gene ~ cluster_id, value.var = "log1ExpZRank")
rownames(linnarssonMatrixHumanReachable) <- linnarssonMatrixHumanReachable$Gene
linnarssonMatrixHumanReachable$Gene <- NULL

save(linnarssonMatrixMouse, linnarssonMatrixHumanReachable, file=here("data", "Zeisel_etal_2018","l5_all.agg.tab.processed.RData"))

