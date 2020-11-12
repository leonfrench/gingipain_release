library(readr)
library(dplyr)
library(magrittr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(here)

#! This code requires running Allen_expression_for_region_specificty.ipynb in the root folder first

enrichment_table <- read_csv(here('data', 'processed_HBA', 'srp_pgpd_chat_avp_enrichment_table.csv'))
enrichment_table

chol_labels <- c('substantia innominata', 'basal nucleus of meynert', 'nucleus of the diagonal band', 
                 'bed  nucleus of stria terminalis', 'olfactory tubercle','islands of calleja')

hypo_labels <- c('anterior hypothalamic area', 'arcuate nucleus of the hypothalamus', 
                 'lateral hypothalamic area, anterior region', 'ventromedial hypothalamic nucleus',
                 'dorsomedial hypothalamic nucleus', 'pallidohypothalamic nucleus', 
                 'paraventricular nucleus of the hypothalamus', 'lateral hypothalamic area, tuberal region',
                 'posterior hypothalamic area', 'lateral hypothalamic area, mammillary region')


donor_info <- read_csv(here("/data/Donor_information.csv")) %>% select(donor = folder_name, id)
enrichment_table <- inner_join(enrichment_table, donor_info) %>% mutate(donor = paste(donor, id, sep="/"))

#setup the donor order
enrichment_table %<>% mutate(donor = factor(donor, levels= unique(c("9861/H0351.2001", "10021/H0351.2002", sort(unique(enrichment_table$donor))))))

#########################################################################################################
# Plot PG_PD_genes vs CHAT expression
#########################################################################################################
labels <- chol_labels
ggplot(enrichment_table, aes(y=pg_pd_AUC, x=CHAT, label=brain_structure)) +
  geom_point(data = enrichment_table %>% filter(!str_detect(brain_structure, paste(labels, collapse='|'))), color= 'darkgrey') +
  geom_point(data = enrichment_table %>% filter(str_detect(brain_structure, paste(labels, collapse='|'))), color= 'red') +
  geom_text_repel(data = enrichment_table %>% filter(str_detect(brain_structure, paste(labels, collapse='|'))),
                  #force = 10,
                  #box.padding = 1.8,
                  size = 2.5,
                  nudge_y = 0.05,
                  nudge_x = 0.1,
                  direction='y',
                  segment.size  = 0.2,
                  hjust = 0,
                  #xlim  = c(10,NA),
                  segment.color = "red") +
  facet_wrap(~ donor) + 
  scale_y_continuous(limits = c(0, 1)) + 
  xlim(min(enrichment_table$CHAT), max(enrichment_table$CHAT)+10) +
  ylab('PG_PD enriched genes (AUC)') + 
  xlab('CHAT gene expression (log2 intensity)') + 
  theme_bw()

ggsave(filename = here('results', 'microarray', 'HBA-PG_PD-CHAT_enrichment.pdf'), width=10, height=7, dpi=300)
#for PDF use 10x7


# same figure but highlighting the structures in hypothalamus
labels <- hypo_labels
ggplot(enrichment_table, aes(y=pg_pd_AUC, x=CHAT, label=brain_structure)) +
  geom_point(data = enrichment_table %>% filter(!str_detect(brain_structure, paste(labels, collapse='|'))), color= 'darkgrey') +
  geom_point(data = enrichment_table %>% filter(str_detect(brain_structure, paste(labels, collapse='|'))), color= 'red') +
  geom_text_repel(data = enrichment_table %>% filter(str_detect(brain_structure, paste(labels, collapse='|'))),
                  #force = 10,
                  #box.padding = 1.8,
                  size = 2.5,
                  nudge_y = 0.05,
                  nudge_x = 0.1,
                  #nudge_x = -0.5,
                  direction='y',
                  segment.size  = 0.2,
                  hjust = 0,
                  #xlim  = c(10,NA),
                  segment.color = "red") +
  facet_wrap(~ donor) + 
  scale_y_continuous(limits = c(0, 1)) + 
  xlim(min(enrichment_table$CHAT), max(enrichment_table$CHAT)+10) +
  ylab('PG_PD enriched genes (AUC)') + 
  xlab('CHAT gene expression (log2 intensity)') + 
  theme_bw()

ggsave(filename = here('results', 'microarray', 'HBA-PG_PD-CHAT_enrichment_hypothalamus.pdf'), width=10, height=7, dpi=300)

#########################################################################################################
# Plot PG_PD AUC vs AVP expression
#########################################################################################################
labels <- chol_labels
ggplot(enrichment_table, aes(y=pg_pd_AUC, x=AVP, label=brain_structure)) +
  geom_point(data = enrichment_table %>% filter(!str_detect(brain_structure, paste(labels, collapse='|'))), color= 'darkgrey') +
  geom_point(data = enrichment_table %>% filter(str_detect(brain_structure, paste(labels, collapse='|'))), color= 'red') +
  geom_text_repel(data = enrichment_table %>% filter(str_detect(brain_structure, paste(labels, collapse='|'))),
                  #force = 10,
                  #box.padding = 1.8,
                  size = 2.5,
                  #nudge_y = 0.6,
                  #nudge_x      = 0.15,
                  #nudge_x = -0.5,
                  direction='y',
                  segment.size  = 0.2,
                  hjust = 0,
                  #xlim  = c(10,NA),
                  segment.color = "red") +
  facet_wrap(~ donor) + 
  scale_y_continuous(limits = c(0, 1)) + 
  xlim(min(enrichment_table$AVP), max(enrichment_table$AVP)+10) +
  ylab('PG_PD enriched genes (AUC)') + 
  xlab('AVP gene expression (log2 intensity)') + 
  theme_bw()

ggsave(filename = here('results', 'microarray', 'HBA-PG_PD-AVP_enrichment.pdf'), width=10, height=7, dpi=300)
#for PDF use 10x7

# same figure but highlighting the structures in hypothalamus
labels <- hypo_labels
ggplot(enrichment_table, aes(y=pg_pd_AUC, x=AVP, label=brain_structure)) +
  geom_point(data = enrichment_table %>% filter(!str_detect(brain_structure, paste(labels, collapse='|'))), color= 'darkgrey') +
  geom_point(data = enrichment_table %>% filter(str_detect(brain_structure, paste(labels, collapse='|'))), color= 'red') +
  geom_text_repel(data = enrichment_table %>% filter(str_detect(brain_structure, paste(labels, collapse='|'))),
                  #force = 10,
                  #box.padding = 1.8,
                  size = 2.5,
                  #nudge_y = 0.6,
                  #nudge_x      = 0.15,
                  #nudge_x = -0.5,
                  direction='y',
                  segment.size  = 0.2,
                  hjust = 0,
                  #xlim  = c(10,NA),
                  segment.color = "red") +
  facet_wrap(~ donor) + 
  scale_y_continuous(limits = c(0, 1)) + 
  xlim(min(enrichment_table$AVP), max(enrichment_table$AVP)+10) +
  ylab('PG_PD enriched genes (AUC)') + 
  xlab('AVP gene expression (log2 intensity)') + 
  theme_bw()

ggsave(filename = here('results', 'microarray', 'HBA-PG_PD-AVP_enrichment_hypothalamus.pdf'), width=10, height=7, dpi=300)
#for PDF use 10x7

#########################################################################################################
#########################################################################################################

# correlations
cor.test(enrichment_table$CHAT, enrichment_table$srp_AUC, method = 's')

enrichment_table %>% 
  group_by(donor) %>% 
  summarise(correlation = cor(srp_AUC, CHAT, method='s'))

enrichment_table %>% 
  group_by(donor) %>% 
  summarise(meanSRP = mean(srp_AUC), meanCHAT = mean(CHAT))

#regions with both high colinergic and high SRP
means <- enrichment_table %>% 
  group_by(brain_structure) %>% 
  summarise(meanSRP = mean(srp_AUC), meanCHAT = mean(CHAT)) %>% arrange(-meanSRP)
means %>% mutate(SRPrank = rank(-meanSRP), CHATrank = rank(-meanCHAT)) %>% mutate(combined = SRPrank+CHATrank) %>% arrange(combined)
