library(here)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(cowplot)

#! This code requires running Allen_expression_for_region_specificty.ipynb in the root folder first
pg <- read_csv(here('/results/microarray/pg_pd_six_brains_regional_enrichment.csv')) %>% 
  select(brain_structure = X1, pg_AUC = AUROC)
srp <- read_csv(here('/results/microarray/srp_six_brains_regional_enrichment.csv')) %>% 
  select(brain_structure = X1, srp_AUC = AUROC)

####################
# Zeisel data
####################
pg_zeisel <- read_csv(here("/results", "mouse", "Supplement table - PG_PD - All Zeisel clusters.csv")) %>% 
  select(cluster_id = celltype, pg_AUC = AUC, Name = Description)
srp_zeisel <- read_csv(here("/results", "mouse", "Supplement table - All Zeisel clusters.csv")) %>% 
  select(cluster_id = Cluster_ID, srp_AUC = SRP_AUC, Name)
####################


srp_range <- range(c(srp$srp_AUC, srp_zeisel$srp_AUC))

chol_labels <- c('substantia innominata', 'basal nucleus of meynert', 'nucleus of the diagonal band', 
                 'bed  nucleus of stria terminalis', 'olfactory tubercle','islands of calleja')

hypo_labels <- c('anterior hypothalamic area', 'arcuate nucleus of the hypothalamus', 
                 'lateral hypothalamic area, anterior region', 'ventromedial hypothalamic nucleus',
                 'dorsomedial hypothalamic nucleus', 'pallidohypothalamic nucleus', 
                 'paraventricular nucleus of the hypothalamus', 'lateral hypothalamic area, tuberal region',
                 'posterior hypothalamic area', 'lateral hypothalamic area, mammillary region')


df <- left_join(pg, srp, by = 'brain_structure') %>%
  mutate(region = case_when(
    str_detect(brain_structure, 'ypothal') ~ "hypothalamic",
    str_detect(brain_structure, paste(chol_labels, collapse = '|')) ~ "substantia\ninnominata\nsubregions"
  )) %>% 
  replace_na(list(region = 'other')) 

# label only sub innom and anterior hypo
to_label <- c('substantia innominata', 'anterior hypothalamic area')
Palette1 <- c('blue', 'red', "grey")
df$region <- factor(df$region, levels = c( 'hypothalamic', 'substantia\ninnominata\nsubregions', 'other'))
levels(df$region)

#correlation label
corSummary <- cor.test(x = df$pg_AUC, y = df$srp_AUC, method = 'p')
correlationSummary <- paste0("\n  r = ", signif(corSummary$estimate, digits=2), "  \n  p-value = ", signif(corSummary$p.value, digits=2),"\n")

p1 <- ggplot(df, aes(x=pg_AUC, y=srp_AUC, color = region, label=brain_structure)) +
  theme_bw() + geom_blank() +
  geom_smooth(method = 'lm', se=F, colour='darkgrey') +
  geom_point(data=df %>% filter(region == 'other'), alpha=0.5) +
  geom_point(data=df %>% filter(region != 'other')) +
  scale_colour_manual(values=Palette1, name="Region") +
  geom_hline(yintercept = 0.5, color="darkgrey", size=.5) +
  geom_vline(xintercept = 0.5, color="darkgrey", size=.5) +
  geom_text_repel(data = df %>% filter(str_detect(brain_structure, paste(to_label, collapse='|'))), 
                  size = 3, segment.size  = 0.1, nudge_y = 0.05, color='black') +
  coord_fixed(ratio=1) +scale_y_continuous('ER translocation AUC', limits = srp_range) + scale_x_continuous('Pg detected AUC') +
  annotate("text", label=correlationSummary, x=-Inf, y=Inf, hjust=0, vjust=1, size = 3) +
  theme(legend.position = c(1, 0), legend.justification=c(1,0)) +
  theme(legend.box.background = element_rect(color="darkgrey",size=1)) 
p1

#ggsave(filename = here('/results', 'microarray', 'HBA-srp_auc_vs_pgpd_auc.pdf'), dpi=300)

####################
# Zeisel data
####################
pg_zeisel <- read_csv(here("/results", "mouse", "Supplement table - PG_PD - All Zeisel clusters.csv")) %>% 
  select(cluster_id = celltype, pg_AUC = AUC, Name = Description)
srp_zeisel <- read_csv(here("/results", "mouse", "Supplement table - All Zeisel clusters.csv")) %>% 
  select(cluster_id = Cluster_ID, srp_AUC = SRP_AUC, Name)

df_zeisel <- left_join(pg_zeisel, srp_zeisel, by = c('cluster_id', 'Name')) %>%
  mutate(celltype = case_when(
    str_detect(Name, 'ypothal') ~ "hypothalamic",
    str_detect(Name, 'holine') ~ "cholinergic"
  )) %>% 
  replace_na(list(celltype = 'other'))


corSummary <- cor.test(x = df_zeisel$pg_AUC, y = df_zeisel$srp_AUC, method = 'p')
correlationSummary <- paste0("\n  r = ", signif(corSummary$estimate, digits=2), "  \n  p-value = ", signif(corSummary$p.value, digits=2),"\n")
df_zeisel %<>% mutate(celltype = factor(celltype, levels = c("hypothalamic", "cholinergic", "other")))

p2 <- ggplot(df_zeisel, aes(x=pg_AUC, y=srp_AUC, color = celltype, label=cluster_id)) +
  theme_bw() + geom_blank() +
  geom_point(data=df_zeisel %>% filter(celltype == 'other'), alpha=0.5) +
  geom_smooth(method = 'lm', se=F, colour='darkgrey') +
  geom_point(data=df_zeisel %>% filter(celltype != 'other')) +
  scale_colour_manual(values=Palette1, name="Cell type") +
  geom_hline(yintercept = 0.5, color="darkgrey", size=.5) +
  geom_vline(xintercept = 0.5, color="darkgrey", size=.5) +
  coord_fixed(ratio=1) +scale_y_continuous('ER translocation AUC', limits = srp_range) + scale_x_continuous('Pg detected AUC') +
  #geom_text() + 
  annotate("text", label=correlationSummary, x=-Inf, y=Inf, hjust=0, vjust=1, size = 3) +
  theme(legend.position = c(1, 0), legend.justification=c(1,0)) +
  theme(legend.box.background = element_rect(color="darkgrey",size=1)) 

p2

plot_grid(p1, p2, labels = c('A', 'B'))

ggsave(filename = here('results', 'mouse', 'combined-srp_auc_vs_pgpd_auc.pdf'), height = 7, width = 8.25)

