library(biomaRt)
library(reshape2)
library(readr)
library(here)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)
library(tmod)
source(here("R", "GeneSetBuilders.R"))
detach("package:XML", unload = TRUE) #handles count()/n() conflict
detach("package:dplyr", unload = TRUE)
library(dplyr)
library(conflicted)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("count", "dplyr")


#enrichment of postmortem brain results
per_gene_pg_results <- read_csv(here("results", "GSE68719", "genewise_pg_pd_controls.csv"))
per_gene_pg_results %>% filter(p.adjust < 0.05) %>% count()

#add names
human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
desc <-getBM(attributes=c("description", "ensembl_gene_id", "hgnc_symbol"), filters="hgnc_symbol", values=per_gene_pg_results$gene_symbol, mart=human) %>% as_tibble()
desc_slim <- desc %>% select(name = description, gene_symbol = hgnc_symbol) %>% distinct()
desc_slim %<>% mutate(name = gsub(" \\[.*\\]", "", name))
with_names <- right_join(desc_slim, per_gene_pg_results)
write_csv(with_names , here("results", "GSE68719", "genewise_pg_pd_controls.descriptions.csv"))
with_names
write_csv(with_names %>% head(10) %>% dplyr::select(name, symbol = gene_symbol, estimate, pFDR = p.adjust) %>% mutate_if(is.numeric, signif, digits=3),
          here("results","GSE68719", "genewise_pg_pd_controls.top10.csv"))
write_csv(with_names %>% tail(10) %>% dplyr::select(name, symbol = gene_symbol, estimate, pFDR = p.adjust) %>% mutate_if(is.numeric, signif, digits=3),
          here("results","GSE68719", "genewise_pg_pd_controls.bottom10.csv"))


per_gene_pg_results %>% arrange(p.value) %>% head(20)
#number of significant genes
dim(per_gene_pg_results)
per_gene_pg_results %>% filter(p.adjust < 0.05) %>% count()
hist(per_gene_pg_results$p.value)
per_gene_pg_results %>% filter(p.adjust < 0.05) %>% group_by(statistic < 0) %>% count()

#run cell-type estimates
cell_type_markers <- loadFileSets("Darmanis", convertToMouse = FALSE)

per_gene_pg_results %<>% arrange(-signed_log_p)
sortedGenes <- per_gene_pg_results$gene_symbol
filterGenes <- FALSE
result <- tbl_df(tmodUtest(c(sortedGenes), mset=cell_type_markers, qval = 1.01, filter = filterGenes))
#double because tmod doesn't
result %<>% mutate(P.Value=P.Value*2) 

#adjust again now that identical groups have been merged
result %<>% ungroup() %>% mutate(adj.P.Val=p.adjust(P.Value, method = "fdr")) %>% arrange(P.Value)
result

per_gene_pg_results %>% filter(gene_symbol %in% unlist(cell_type_markers$MODULES2GENES$Endothelial)) %>% as.data.frame()
per_gene_pg_results %>% filter(gene_symbol %in% unlist(cell_type_markers$MODULES2GENES$Neuron)) %>% as.data.frame()


#test relationship between rank and r + k ratio, not significant
rk_proportions_whole_genome <- read_csv(here("/results/residue_enrichment/rk_proportions_whole_genome.csv"))
just_srp_genes <- inner_join(just_srp_genes, rk_proportions_whole_genome)
cor.test(just_srp_genes$percent_target_amino_acids, just_srp_genes$rank, m='s')

joined_proportions <- inner_join(per_gene_pg_results, rk_proportions_whole_genome)
joined_proportions %>% group_by(p.adjust < 0.05, statistic < 0) %>% summarise(n=n(), proportion = mean(percent_target_amino_acids))
#test for difference between top and bottom hits
wilcox.test(percent_target_amino_acids ~ (statistic < 0), joined_proportions %>% filter(p.adjust < 0.05))
cor.test(joined_proportions$signed_log_p, joined_proportions$percent_target_amino_acids, m='s')

#not filtering for all genes in GO - does it matter? run all GO groups
if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { 
} else {
  geneSetsGO <- loadGOSets(sortedGenes)
}

#genes overlaping with SRP and in the top lists - overlap is week - only 11, signal isn't enriched at the top
per_gene_pg_results %>% filter(gene_symbol %in% geneSetsGO$MODULES2GENES$`GO:0006614`) %>% filter(p.adjust < 0.05) %>% group_by(statistic < 0) %>% count()
per_gene_pg_results %>% filter(p.adjust < 0.05) %>% group_by(statistic < 0) %>% count()
per_gene_pg_results %>% count()
per_gene_pg_results %>% filter(gene_symbol %in% geneSetsGO$MODULES2GENES$`GO:0006614`) %>% count()


filterGenes <- TRUE
result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSetsGO, qval = 1.01, filter = filterGenes))

result %<>% rowwise() %>% mutate(aspect = Ontology(ID))

#merge same GO groups using the genes in the set - these are sets that have the same set of genes
result %<>% rowwise() %>% mutate(genes = paste(sort(unlist(geneSetsGO$MODULES2GENES[ID])), collapse = " "))
result %<>% ungroup() %>% group_by(genes, N1) %>% arrange(Title) %>% 
  summarize(MainTitle = first(Title), ID=first(ID), all_IDs=paste(ID, collapse=","), 
            AUC = first(AUC), P.Value= first(P.Value), aspect= first(aspect), 
            otherNames = if_else(length(unique(Title)) > 1, paste(Title[2:length(Title)], collapse=", "), ""))
result %<>% ungroup() %>% dplyr::select(-genes)

#double p-value because tmod doesn't 
result %<>% mutate(P.Value=P.Value*2) 

#adjust again now that identical groups have been merged
result %<>% ungroup() %>% mutate(adj.P.Value=p.adjust(P.Value, method = "fdr")) %>% arrange(P.Value)
result %<>% dplyr::select(MainTitle, ID, geneCount = N1, AUC, P.Value, adj.P.Value, everything()) 

result %>% filter(adj.P.Value < 0.025) %>% group_by(AUC > 0.5) %>% count()
nrow(result)

write_csv(result, here("results", "GSE68719", "genewise_pg_pd_controls_GO_results.csv"))

write_csv(result %>% filter(AUC>0.5) %>% head(10) %>% dplyr::select(Title = MainTitle, ID, `Gene count` = geneCount, AUC, p = P.Value, pFDR = adj.P.Value) %>% mutate_if(is.numeric, signif, digits=3),
          here("results","GSE68719", "GO_enrichment_top_ten.csv"))
write_csv(result %>% filter(AUC<0.5) %>% head(10) %>% dplyr::select(Title = MainTitle, ID, `Gene count` = geneCount, AUC, p = P.Value, pFDR = adj.P.Value) %>% mutate_if(is.numeric, signif, digits=3),
          here("results","GSE68719", "GO_enrichment_bottom_ten.csv"))

#down-regulated GO groups
per_gene_pg_results %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0048167`)) %>% arrange(signed_log_p)
per_gene_pg_results %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0007156`)) %>% arrange(signed_log_p)
intersect(unlist(geneSetsGO$MODULES2GENES$`GO:0048167`), unlist(geneSetsGO$MODULES2GENES$`GO:0007156`))

source(here("R", "RunReviGO.R"))
#run up-regulated GO groups
go_results <- result %>% filter(AUC > 0.5, adj.P.Value < 0.05) %>% dplyr::select(ID, P.Value)
revigo_results <- run_revigo(go_results)
revigo_results %>% filter(dispensability == 0)

per_gene_pg_results %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:1990823`))

go_results <- result %>% filter(AUC < 0.5, adj.P.Value < 0.05) %>% dplyr::select(ID, P.Value)
revigo_results <- run_revigo(go_results)
revigo_results %>% filter(dispensability == 0)
per_gene_pg_results %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0050905`))
per_gene_pg_results %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0016052`))
#aminoacyl-tRNA ligase activity
per_gene_pg_results %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0004812`))
per_gene_pg_results %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0007156`))
per_gene_pg_results %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0005125`)) #cytokine activity


#compare top ten lists with r + k
top_ten_diff_exp <- as.character(result %>% filter(AUC > 0.5) %>% head(10) %>% pull(MainTitle))
top_ten_r_k <- read_csv(here("results", "residue_enrichment", "GO_enrichment_top_ten.csv")) %>% pull(Title)
intersect(top_ten_diff_exp, top_ten_r_k)
length(intersect(top_ten_diff_exp, top_ten_r_k))
setdiff(top_ten_diff_exp, top_ten_r_k)
setdiff(top_ten_r_k, top_ten_diff_exp)
