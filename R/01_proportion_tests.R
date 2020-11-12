library(ggplot2)
library(here)
library(readr)
library(stringr)
library(magrittr)
library(dplyr)
library(seqinr)
source(here("R", "GeneSetBuilders.R"))
detach("package:XML", unload = TRUE) #handles count()/n() conflict
detach("package:dplyr", unload = TRUE)
library(dplyr)

genecode_fasta <- read.fasta(file = here("data", "protein_gencode", "gencode.v32.pc_translations.fa.gz"))

name_table <- tibble(piped_name = names(genecode_fasta))
name_table %<>% dplyr::mutate(order = 1:n())

#extract gene_symbol from second last token
name_table %<>%
  mutate(
    splits = strsplit(piped_name, "[|]")
  ) %>% 
  rowwise() %>% 
  mutate(
    gene_symbol = splits[length(splits)-1]
  )


dir.create(here("results"), showWarnings = F)
dir.create(here("data", "gene_lists"), showWarnings = F)


#make sure arrangement is correct
genecode_fasta <- genecode_fasta[name_table$order]
length(genecode_fasta)

symbol_to_sequence <- data.frame(gene_symbol= name_table$gene_symbol, names=names(genecode_fasta), stringsAsFactors = F) %>% as_tibble()
symbol_to_sequence$protein_sequence = unlist(unname(lapply(genecode_fasta, c2s)))
symbol_to_sequence %<>% mutate(length = nchar(protein_sequence))
symbol_to_sequence %<>% filter(protein_sequence != "m") #clean up sequences set to just "m"
symbol_to_sequence %<>% select(-names) %>% distinct()

symbol_to_sequence$gene_symbol %>% unique() %>% length()
#remove sequences mapped to more than one gene/protein
sequence_to_gene_map <- symbol_to_sequence %>% group_by(protein_sequence) %>% summarize(n = n(), genes = list(gene_symbol)) %>% arrange(-n) %>% filter(n>1)
symbol_to_sequence %<>% filter(!protein_sequence %in% sequence_to_gene_map$protein_sequence)
symbol_to_sequence$gene_symbol %>% unique() %>% length()

symbol_to_sequence %>% filter(gene_symbol == "HIST1H2BL") %>% .$protein_sequence
symbol_to_sequence %>% filter(gene_symbol == "HIST1H2BM")%>% .$protein_sequence
symbol_to_sequence %>% filter(gene_symbol == "HIST2H2BE")%>% .$protein_sequence

gene_universe <- symbol_to_sequence$gene_symbol %>% unique()
gene_universe %>% length()

write_csv(gene_universe %>% as_tibble(), here("data", "gene_lists", "gene_universe.txt"), col_names = F)

#load the GO sets we've been using throughout
if (exists("geneSetsGO") && length(geneSetsGO$MODULES2GENES) > 1000 ) { 
} else {
  geneSetsGO <- loadGOSets(gene_universe)
}
filterGenes <- TRUE #work within the GO univerise
geneSetsGO

#run TMOD on r + k ranking
amino_acid_first_target <- "r"
amino_acid_second_target <- "k"
symbol_to_sequence %<>% mutate(r_count = str_count(protein_sequence, "r"), k_count = str_count(protein_sequence, "k"), 
                               target_pair_count = str_count(protein_sequence, amino_acid_first_target) + str_count(protein_sequence, amino_acid_second_target))
symbol_to_sequence %<>% mutate(percent_target_amino_acids = target_pair_count/length)

symbol_to_sequence_averaged <- symbol_to_sequence %>% group_by(gene_symbol) %>% summarize_if(is.numeric, mean)
symbol_to_sequence_averaged %<>% arrange(-percent_target_amino_acids)

dir.create(here("results","residue_enrichment"), showWarnings = F)
write_csv(symbol_to_sequence_averaged, here("results","residue_enrichment", "rk_proportions_whole_genome.csv"))

mean(symbol_to_sequence_averaged$percent_target_amino_acids)
symbol_to_sequence_averaged
proteins_used <- nrow(symbol_to_sequence_averaged)
symbol_to_sequence_averaged[round(proteins_used/40),]
symbol_to_sequence_averaged[round(proteins_used-proteins_used/40),]

sortedGenes <- symbol_to_sequence_averaged$gene_symbol

result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSetsGO, qval = 1.01, filter = filterGenes))
result %<>% rowwise() %>% mutate(aspect = Ontology(ID))
#merge same GO groups using the genes in the set - these are sets that have the same set of genes
result %<>% rowwise() %>% mutate(genes = paste(sort(unlist(geneSetsGO$MODULES2GENES[ID])), collapse = " "))
result %<>% ungroup() %>% group_by(genes, N1) %>% arrange(Title) %>% 
  summarize(MainTitle = dplyr::first(Title), ID=dplyr::first(ID), all_IDs=paste(ID, collapse=","), 
            AUC = dplyr::first(AUC), P.Value= dplyr::first(P.Value), aspect= dplyr::first(aspect), 
            otherNames = if_else(length(unique(Title)) > 1, paste(Title[2:length(Title)], collapse=", "), ""))
result %<>% ungroup() %>% dplyr::select(-genes)


#convert to one-sided p-value
result %<>% mutate(P.Value = if_else(AUC < 0.5, 1-P.Value, P.Value))
#readjust
result %<>% mutate(adj.P.Val = p.adjust(P.Value, method="fdr"))
result %<>% arrange(P.Value)
result %>% group_by(adj.P.Val < 0.05) %>% count()
result %>% group_by(adj.P.Val < 0.05) %>% filter(aspect == "BP") %>% count()


#write out table
dir.create(here("results","residue_enrichment"), showWarnings = F)
write_csv(result, here("results","residue_enrichment", "GO_enrichment.csv"))

write_csv(result %>% head(10) %>% dplyr::select(Title = MainTitle, ID, `Gene count` = N1, AUC, p = P.Value, pFDR = adj.P.Val) %>% mutate_if(is.numeric, signif, digits=3),
          here("results","residue_enrichment", "GO_enrichment_top_ten.csv"))

#write out SRP
symbol_to_sequence_averaged %>% group_by(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0006614`)) %>% summarize(mean = mean(percent_target_amino_acids))
write_csv(symbol_to_sequence_averaged %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0006614`)), 
          here("results","residue_enrichment", "SRP_stats.csv"))
write_csv(symbol_to_sequence_averaged %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0006614`)) %>% dplyr::select(gene_symbol), 
          here("data", "gene_lists", "SRP_list.txt"), col_names = F)

#write out nucleosome / histone -check polygenic
symbol_to_sequence_averaged %>% group_by(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0000786`)) %>% summarize(mean = mean(percent_target_amino_acids))
symbol_to_sequence_averaged %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0000786`))
write_csv(symbol_to_sequence_averaged %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0000786`)), 
          here("results","residue_enrichment", "nucleosome_stats.csv"))
write_csv(symbol_to_sequence_averaged %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0000786`)) %>% dplyr::select(gene_symbol), 
          here("data","gene_lists", "nucleosome_list.txt"), col_names = F)

symbol_to_sequence_averaged %>% group_by(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0098981`)) %>% summarize(mean = mean(percent_target_amino_acids))
write_csv(symbol_to_sequence_averaged %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0098981`)) %>% dplyr::select(gene_symbol), 
          here("data","gene_lists", "cholinergic_synapse_list.txt"), col_names = F)

symbol_to_sequence_averaged %>% filter(gene_symbol %in% unlist(geneSetsGO$MODULES2GENES$`GO:0030010`))  

#sumarize GO groups - optional
source(here("R", "RunReviGO.R"))
go_results <- result %>% filter(AUC > 0.5, adj.P.Val < 0.05) %>% dplyr::select(ID, P.Value)
revigo_results <- run_revigo(go_results)
revigo_results %>% filter(dispensability == 0)
revigo_results %>% filter(dispensability < 0.1)
revigo_results %>% filter(dispensability < 0.2)

#compare SRP to other pair combinations to show it's specific
amino_acid_alphabet <- c("a", "c", "d", "e", "f", "g", "h", "i", "k", "l", "m", "n", "p", "q", "r", "s", "t", "v", "w", "y")

all_pair_results <- as_tibble(NULL)
ptm <- Sys.time()
for(amino_acid_first_target in amino_acid_alphabet) {
  for(amino_acid_second_target in amino_acid_alphabet) {
    if (nrow (all_pair_results %>% filter(amino_acid_first == amino_acid_second_target & amino_acid_second == amino_acid_first_target)) == 1) {
      next
    }
    print(paste0(amino_acid_first_target, amino_acid_second_target))
    gene_stats_joined <- NULL
    symbol_to_sequence$target_pair_count <- NULL #clear for safety
    #count the amino acid hits
    symbol_to_sequence %<>% mutate(target_pair_count = str_count(protein_sequence, amino_acid_first_target) + str_count(protein_sequence, amino_acid_second_target))
    symbol_to_sequence %<>% mutate(percent_target_amino_acids = target_pair_count/length)
    
    #mean average percent across the many entries per gene
    symbol_to_sequence_averaged <- symbol_to_sequence %>% group_by(gene_symbol) %>% summarize(percent_target_amino_acids = mean(percent_target_amino_acids)) %>% arrange(-percent_target_amino_acids)
    
    symbol_to_sequence_averaged %<>% arrange(-percent_target_amino_acids)
    
    mean(symbol_to_sequence_averaged$percent_target_amino_acids)
    
    sortedGenes <- symbol_to_sequence_averaged$gene_symbol
    
    result <- tbl_df(tmodUtest(c(sortedGenes), mset=geneSetsGO, qval = 1.01, filter = filterGenes))
    

    #merge same GO groups using the genes in the set - these are sets that have the same set of genes
    #code duplication from above
    result %<>% rowwise() %>% mutate(genes = paste(sort(unlist(geneSetsGO$MODULES2GENES[ID])), collapse = " "))
    result %<>% ungroup() %>% group_by(genes, N1) %>% arrange(Title) %>% 
      summarize(MainTitle = dplyr::first(Title), ID=dplyr::first(ID), all_IDs=paste(ID, collapse=","), 
                AUC = dplyr::first(AUC), P.Value= dplyr::first(P.Value), 
                otherNames = if_else(length(unique(Title)) > 1, paste(Title[2:length(Title)], collapse=", "), ""))
    result %<>% ungroup() %>% dplyr::select(-genes)
    
    
    #convert to one-sided p-value
    result %<>% mutate(P.Value = if_else(AUC < 0.5, 1-P.Value, P.Value))
    #readjust
    result %<>% mutate(adj.P.Val = p.adjust(P.Value, method="fdr"))
    result %<>% mutate(rank = rank(P.Value))
    result %>% arrange(P.Value)
    
    #get rank of SRP
    SRP_rank <- result %>% filter(ID == "GO:0006614") %>% head(1) %>% .$rank
    SRP_adj.P.Value <- result %>% filter(ID == "GO:0006614") %>% head(1) %>% .$adj.P.Val
    SRP_AUC <- result %>% filter(ID == "GO:0006614") %>% head(1) %>% .$AUC
    #print(SRP_adj.P.Value)
    #print(SRP_rank)
    #print(SRP_AUC)
    
    #symbol_to_sequence_averaged %>% group_by(gene_symbol %in% geneSetsGO$MODULES2GENES$`GO:0006614`) %>% summarize(percent_target_amino_acids = mean(percent_target_amino_acids))
    
    single_row_result <- data.frame(amino_acid_first = amino_acid_first_target, amino_acid_second = amino_acid_second_target, SRP_AUC = SRP_AUC, SRP_rank = SRP_rank,SRP_adj.P.Value =SRP_adj.P.Value, stringsAsFactors = F)
    all_pair_results <- bind_rows(all_pair_results, single_row_result)
  }
}
Sys.time() - ptm

write_csv(all_pair_results, here("results", "residue_enrichment", "Amino_acid_pair_enrichment_GO0006614.csv"))
all_pair_results %>% arrange(-SRP_AUC) %>% filter(SRP_adj.P.Value < 0.05) %>% as.data.frame()
#level of single residues
all_pair_results %>% filter(amino_acid_second == amino_acid_first) %>% arrange(-SRP_AUC)

all_pair_results <- read_csv(here("results", "residue_enrichment", "Amino_acid_pair_enrichment_GO0006614.csv"))
all_pair_results_dup <- all_pair_results
all_pair_results_dup %<>% filter(amino_acid_second != amino_acid_first)
all_pair_results_dup$SRP_AUC = NA
all_pair_results_dup %<>% dplyr::rename(amino_acid_first = amino_acid_second, amino_acid_second = amino_acid_first)
all_pair_results <- bind_rows(all_pair_results, all_pair_results_dup)
all_pair_results %>% arrange(-SRP_AUC)

#amino acid combination plot
ggplot(data = all_pair_results, aes(x = amino_acid_first, y = amino_acid_second)) +
  geom_tile(aes(fill = SRP_AUC), colour='white') +
  ylab("Amino acid residue") + xlab("Amino acid residue") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0), limits = rev(unique(all_pair_results$amino_acid_first))) +
  scale_fill_distiller(palette = 'RdBu', limits = c(0,1), na.value = "white", name ="ER translocation\ngene enrichment\n(AUC)" ) +
  theme(legend.position = c(0.85, 0.8))+ coord_fixed()



