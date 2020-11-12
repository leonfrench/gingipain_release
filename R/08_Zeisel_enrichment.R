library(huxtable)
library(here)
library(readr)
library(magrittr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(readxl)
library(ggplot2)
library(ggrepel)
library(homologene)
library(pROC)
source(here('R', 'AUCFunction.R'))


apply_MWU <- function(column, targetIndices) {
  wilcox.test(column[targetIndices], column[!targetIndices], conf.int = F)$p.value
}

convert_genes <- function(input_genes) {
  mouse_genes <- human2mouse(input_genes)
  return(unique(mouse_genes$mouseGene))
}


get_enrichment <- function(df, gene_list) {
  # df is matrix of ranks, with gene_symbol + celltypes as columns
  forIndices <- as_tibble(df$gene_symbol)
  names(forIndices) <- 'gene_symbol'
  forIndices %<>% mutate(isTargetGene = gene_symbol %in% gene_list)
  
  targetIndices <- forIndices$isTargetGene
  print(paste0(
    sum(forIndices$isTargetGene),
    ' of ',
    length(gene_list),
    ' genes were found in the human reachable linnarson matrix'
  ))
  
  df %<>% select(-gene_symbol)
  AUROC <- map_df(df, auroc_analytic, targetIndices)
  MWU <- map_df(df, apply_MWU, targetIndices)
  
  enrichment <- gather(AUROC, key = 'celltype', value = 'AUC')
  pvals <- gather(MWU, key = 'celltype', value = 'pValue')
  enrichment %<>% inner_join(pvals, by = 'celltype') %>% mutate(rank = rank(desc(AUC))) %>% arrange(desc(AUC))
  #merge in celltype descriptions
  descriptions <-
    read_excel(here(
      'data',
      'Zeisel.TableS3.1-s2.0-S009286741830789X-mmc3.xlsx'
    )) %>% select(`Cluster name`, `Description`, `Probable location (manually assigned)`)
  enrichment %<>% left_join(descriptions, by = c('celltype' = 'Cluster name'))
  enrichment %<>% mutate(adjusted_P = signif(p.adjust(pValue, method = 'fdr'), digits = 3))
}

# load full dataset:
# --> linnarssonMatrixMouse // linnarssonMatrixHumanReachable

load(here('data', 'Zeisel_etal_2018', 'l5_all.agg.tab.processed.RData'), verbose = TRUE)

# we have a matrix of gene ranks by celltype
head(linnarssonMatrixHumanReachable)

# process matrix to be operating in same gene_universe
df <- as_tibble(linnarssonMatrixHumanReachable)
celltypes <- names(df) # store celltype names in order to select them for re-ranking after filtering for gene_universe
df$gene_symbol <- rownames(linnarssonMatrixHumanReachable)
gene_universe <- read_tsv(here('./data/processed_HBA/adult_brainarea_vs_genes_exp_default_donors_10021-9861-14380-15697-15496-12876.tsv')) %>% .$gene_symbol
mouse_gene_universe <- convert_genes(gene_universe)
# requires re-ranking the genes in the matrix after filtering some genes
df %<>% filter(gene_symbol %in% mouse_gene_universe) %>% modify_at(celltypes, rank)

# define gene list to test and convert symbols to mouse
# list loses 3 gene symbols in conversion from human to mouse
gene_list <- read_csv(here('data', 'gene_lists', 'SRP_list.txt'), col_names = 'gene_symbol')

pg_pd_genes <-
  read_csv(here(
    'results','GSE68719', 'genewise_pg_pd_controls.descriptions.csv'
  )) %>%
  filter((estimate > 0) & (p.adjust < 0.05)) %>%
  select(gene_symbol)

mouse_pg_pd_genes <- convert_genes(pg_pd_genes$gene_symbol)

mouse_pg_pd_genes %>% unique() %>% length()
mouse_gene_list <- convert_genes(gene_list$gene_symbol)
dim(gene_list)
length(mouse_gene_list)


srp_enrichment <- get_enrichment(df, mouse_gene_list)
srp_enrichment %<>% mutate(is_cholinergic = str_detect(Description, 'holine'))
pg_pd_enrichment <- get_enrichment(df, mouse_pg_pd_genes)
pg_pd_enrichment %>% filter((adjusted_P < 0.05) & (AUC > 0.5))

pg_pd_enrichment %>% filter( str_detect(Description, 'holine'))

pg_pd_enrichment %>% write_csv(here("results", "mouse", "Supplement table - PG_PD - All Zeisel clusters.csv"))
chat_enrichment <- get_enrichment(df, c('Chat'))
chat_enrichment %<>% mutate(is_cholinergic = str_detect(Description, 'holine'))



# Raw expression of Chat, Avp, Vip
chat_expression <- read_csv(here('results/mouse/Zeisel_Chat_expression.csv')) %>% 
  select(Chat_exp = expression, Cluster_ID = cluster_id)
avp_expression <- read_csv(here('results/mouse/Zeisel_Avp_expression.csv')) %>% 
  select(Avp_exp = expression, Cluster_ID = cluster_id)

dir.create(here("results/mouse/"), showWarnings = F, recursive = T)
supplement_table <- srp_enrichment %>% 
  select(Rank=rank, Cluster_ID=celltype, Name=Description, SRP_AUC=AUC, pValue, adjusted_P, `Probable location (manually assigned)`) %>% 
  inner_join(chat_expression, by='Cluster_ID') %>% 
  inner_join(avp_expression, by='Cluster_ID')

supplement_table %>% write_csv(here("results", "mouse", "Supplement table - All Zeisel clusters.csv"))
supplement_table %>% filter(str_detect(Name, 'holine')) %>% write_csv(here('results/mouse/Zeisel_cholinergic_cells.csv'))
supplement_table %>% filter(str_detect(Name, 'ypothal')) %>% write_csv(here('results/mouse/Zeisel_hypothalamic_cells.csv'))

final_table <- supplement_table %>% 
  select(Cluster_ID, Name, AUROC=SRP_AUC, adjusted_P)

ht <- hux(final_table, add_colnames = TRUE, autoformat = TRUE) %>% 
  set_right_padding(1) %>% 
  set_left_padding(1) %>% 
  set_bold(1, 1:4, TRUE) %>% 
  set_bottom_border(1, 1:4, 1) %>%
  set_right_border(everywhere, 1:3, 0.4) %>% 
  set_align(1:14, 3:4, 'right') %>%
  set_caption('SRP enrichment')

theme_plain(ht)
quick_html(ht, file = here("results", "mouse", "Table - Zeisel top clusters.html"), open=F)


##########################################################################################
# Plot SRP enrichment vs Chat enrichment
##########################################################################################

# define different subsets of the data to facilitate plotting
merged_enrichment <-
  inner_join(
    srp_enrichment %>% select(celltype, AUC_SRP = AUC, is_cholinergic, Description),
    chat_enrichment %>% select(celltype, AUC_CHAT = AUC),
    by = 'celltype'
  )

sub1 <- merged_enrichment %>% filter(AUC_SRP > 0.75)
sub2 <- merged_enrichment %>% filter((0.75 > AUC_SRP) & (AUC_SRP > 0.625))
sub3 <- merged_enrichment %>% filter(AUC_SRP < 0.625)

srp_chat_enrichment <-
  ggplot(merged_enrichment, aes(
    y = AUC_SRP,
    x = AUC_CHAT,
    label = paste(Description, celltype)
  )) +
  geom_point(color = if_else(merged_enrichment$is_cholinergic, 'red', 'black')) +
  geom_text_repel(
    data = sub1 %>% filter(is_cholinergic),
    xlim = 1.05,
    hjust = 1,
    force = 9,
    #box.padding = 0.15,
    direction = 'both',
    size = 3,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  geom_text_repel(
    data = sub2 %>% filter(is_cholinergic),
    xlim = 1.05,
    hjust = 1,
    size = 3,
    direction = 'both',
    nudge_y = -0.1,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  geom_text_repel(
    data = sub3 %>% filter(is_cholinergic),
    xlim = 1.05,
    hjust = 1,
    size = 3,
    direction = 'both',
    nudge_y = -0.2,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  scale_x_continuous(breaks = c(0.25, 0.50, 0.75, 1.00),
                     limits = c(NA, 1.5)) +
  theme_bw()

ggsave(filename = here('results','mouse', 'zeisel-srp_chat_enrichment.png'),
       width = 12, height = 9,
       dpi=300)


##########################################################################################
# Plot SRP enrichment vs Chat expression
##########################################################################################

# define different subsets of the data to facilitate plotting
merged_enrichment <-
  inner_join(
    srp_enrichment %>% select(celltype, AUC_SRP = AUC, is_cholinergic, Description),
    chat_expression %>% select(Cluster_ID, Chat_exp),
    by = c("celltype" = "Cluster_ID")
  )

merged_enrichment %<>% mutate(is_hypothalamic = str_detect(Description, 'ypothalam'))

sub1 <- merged_enrichment %>% filter(AUC_SRP > 0.75)
sub2 <- merged_enrichment %>% filter((0.75 > AUC_SRP) & (AUC_SRP > 0.625))
sub3 <- merged_enrichment %>% filter(AUC_SRP < 0.625)

srp_chat_enrichment <-
  ggplot(merged_enrichment, aes(
    y = AUC_SRP,
    x = Chat_exp,
    label = paste(Description, celltype)
  )) +
  geom_point(color = if_else(merged_enrichment$is_cholinergic, 'red', 'black')) +
  geom_text_repel(
    data = sub1 %>% filter(is_cholinergic),
    xlim = 1,
    hjust = 1,
    force = 9,
    #box.padding = 0.15,
    direction = 'both',
    size = 3,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  geom_text_repel(
    data = sub2 %>% filter(is_cholinergic),
    xlim = 1,
    hjust = 1,
    size = 3,
    direction = 'both',
    nudge_y = -0.1,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  geom_text_repel(
    data = sub3 %>% filter(is_cholinergic),
    xlim = 1,
    hjust = 1,
    size = 3,
    direction = 'both',
    nudge_y = -0.2,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1.00),
                limits = c(NA, 5000)) +
  theme_bw()

srp_chat_enrichment

ggsave(filename = here('results','mouse', 'zeisel-srp_enrichment_chat_exp.png'),
       width = 12, height = 9,
       dpi=300)

#########################
# Hypothalamic highlights
#########################

srp_chat_enrichment <-
  ggplot(merged_enrichment, aes(
    y = AUC_SRP,
    x = Chat_exp,
    label = paste(Description, celltype)
  )) +
  geom_point(color = if_else(merged_enrichment$is_hypothalamic, 'red', 'black')) +
  geom_text_repel(
    data = sub1 %>% filter(is_hypothalamic),
    xlim = 1,
    hjust = 1,
    force = 9,
    #box.padding = 0.15,
    direction = 'both',
    size = 3,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  geom_text_repel(
    data = sub2 %>% filter(is_hypothalamic),
    xlim = 1,
    hjust = 1,
    size = 3,
    direction = 'both',
    nudge_y = -0.1,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  geom_text_repel(
    data = sub3 %>% filter(is_hypothalamic),
    xlim = 1,
    hjust = 1,
    size = 3,
    direction = 'both',
    nudge_y = -0.2,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1.00),
                limits = c(NA, 5000)) +
  theme_bw()

srp_chat_enrichment

ggsave(filename = here('results','mouse', 'zeisel-srp_enrichment_chat_exp_hypothalamic_highlights.png'),
       width = 12, height = 9,
       dpi=300)

##########################################################################################
# Plot SRP enrichment vs Avp expression
##########################################################################################

# define different subsets of the data to facilitate plotting
merged_enrichment <-
  inner_join(
    srp_enrichment %>% select(celltype, AUC_SRP = AUC, is_cholinergic, Description),
    avp_expression %>% select(Cluster_ID, Avp_exp),
    by = c("celltype" = "Cluster_ID")
  )

merged_enrichment %<>% mutate(is_hypothalamic = str_detect(Description, 'ypothalam'))

sub1 <- merged_enrichment %>% filter(AUC_SRP > 0.75)
sub2 <- merged_enrichment %>% filter((0.75 > AUC_SRP) & (AUC_SRP > 0.625))
sub3 <- merged_enrichment %>% filter(AUC_SRP < 0.625)

Avp_chat_enrichment <-
  ggplot(merged_enrichment, aes(
    y = AUC_SRP,
    x = Avp_exp,
    label = paste(Description, celltype)
  )) +
  geom_point(color = if_else(merged_enrichment$is_cholinergic, 'red', 'black')) +
  geom_text_repel(
    data = sub1 %>% filter(is_cholinergic),
    xlim = 1,
    hjust = 1,
    force = 9,
    #box.padding = 0.15,
    direction = 'both',
    size = 3,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  geom_text_repel(
    data = sub2 %>% filter(is_cholinergic),
    xlim = 1,
    hjust = 1,
    size = 3,
    direction = 'both',
    nudge_y = -0.1,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  geom_text_repel(
    data = sub3 %>% filter(is_cholinergic),
    xlim = 1,
    hjust = 1,
    size = 3,
    direction = 'both',
    nudge_y = -0.2,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1.00),
                limits = c(NA, 5000)) +
  theme_bw()

Avp_chat_enrichment

ggsave(filename = here('results','mouse', 'zeisel-srp_enrichment_avp_exp-cholinergic_highlights.png'),
       width = 12, height = 9,
       dpi=300)

#########################
# Hypothalamic highlights
#########################

Avp_chat_enrichment <-
  ggplot(merged_enrichment, aes(
    y = AUC_SRP,
    x = Avp_exp,
    label = paste(Description, celltype)
  )) +
  geom_point(color = if_else(merged_enrichment$is_hypothalamic, 'red', 'black')) +
  geom_text_repel(
    data = sub1 %>% filter(is_hypothalamic),
    xlim = 1,
    hjust = 1,
    force = 9,
    #box.padding = 0.15,
    direction = 'both',
    size = 3,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  geom_text_repel(
    data = sub2 %>% filter(is_hypothalamic),
    xlim = 1,
    hjust = 1,
    size = 3,
    direction = 'both',
    nudge_y = -0.1,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  geom_text_repel(
    data = sub3 %>% filter(is_hypothalamic),
    xlim = 1,
    hjust = 1,
    size = 3,
    direction = 'both',
    nudge_y = -0.2,
    segment.size  = 0.2,
    segment.color = "grey50"
  ) +
  scale_x_log10(breaks = c(0.001, 0.01, 0.1, 1.00),
                limits = c(NA, 5000)) +
  theme_bw()

Avp_chat_enrichment

ggsave(filename = here('results','mouse', 'zeisel-srp_enrichment_avp_exp-hypothalamic_highlights.png'),
       width = 12, height = 9,
       dpi=300)

##########################################################################################

### AUC of AUCs
# 
srp_enrichment %<>% mutate(is_hypothalamic = str_detect(Description, 'ypothalam'))
head(srp_enrichment)
wilcox.test(data=srp_enrichment, rank ~ is_cholinergic)
wilcox.test(data=srp_enrichment, rank ~ is_hypothalamic)
t.test(data=srp_enrichment, AUC ~ is_cholinergic) #to see the means
t.test(data=srp_enrichment, AUC ~ is_hypothalamic) #to see the means
wilcox.test(data=srp_enrichment, AUC ~ is_cholinergic)
wilcox.test(data=srp_enrichment, AUC ~ is_hypothalamic)

srp_enrichment %>% group_by(is_cholinergic) %>% summarize(medianAUC = median(AUC))
srp_enrichment %>% group_by(is_hypothalamic) %>% summarize(medianAUC = median(AUC))

roc(srp_enrichment$is_cholinergic, srp_enrichment$AUC)
roc(is_cholinergic ~ AUC , srp_enrichment)
roc(is_hypothalamic ~ AUC , srp_enrichment)


###################################################################
# plots showing distribution of enrichment scores for all celltypes
# highlighting cholinergic celltypes

# raster
ggplot(srp_enrichment, x=rank, y=is_cholinergic) +
  geom_vline(data = filter(srp_enrichment, is_cholinergic == TRUE), aes(xintercept=rank)) +
  xlim(range(srp_enrichment$rank)) +
  coord_cartesian(expand = F) +
  theme(strip.background = element_blank(), strip.placement = "inside") 

#distplot with rug
ggplot(srp_enrichment, aes(x=AUC)) +
  geom_line(stat = 'density') +
  geom_rug(data=filter(srp_enrichment, is_cholinergic==TRUE))

#displot with lines for cholinergic cells
ggplot(srp_enrichment, aes(x=AUC)) +
  geom_line(stat = 'density') +
  geom_vline(data=filter(srp_enrichment, is_cholinergic==TRUE), aes(xintercept=AUC))
