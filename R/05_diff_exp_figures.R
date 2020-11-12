library(cowplot)
library(plotROC)
library(readxl)
library(reshape2)
library(readr)
library(here)
library(dplyr)
library(tidyr)
library(broom)
library(magrittr)
library(RCurl)
library(ggplot2)
detach("package:XML", unload = TRUE) #handles count()/n() conflict
detach("package:dplyr", unload = TRUE)
library(dplyr)

#read in sample meta data 
pg_detection_data <- read_csv(here("data", target_exp, "pg_detection_data.csv"))

#read in expression data from GEO
expression <- read_tsv(here("data", "GSE68719", "GSE68719_mlpd_PCG_DESeq2_norm_counts.txt.gz"))
target_gene <- "SRP9" 

single_gene <- expression %>% filter(symbol == target_gene) %>% select(-EnsemblID, -symbol)
single_gene <- data.frame(short_title = names(single_gene), expression = t(single_gene)) %>% as_tibble()
pg_detection_data_expression <- inner_join(pg_detection_data, single_gene, by = "short_title")

wilcox.test(expression ~ pg_detected, pg_detection_data_expression)
t.test(expression ~ pg_detected, pg_detection_data_expression)

wilcox.test(expression ~ pg_detected, pg_detection_data_expression %>% filter(group == "C"))

summary(lm(expression ~ pg_detected + group + Bases + age + pmi + rin, data = pg_detection_data_expression))
#just controls
summary(lm(log(1+expression) ~ pg_detected + Bases + age + pmi + rin, data = pg_detection_data_expression %>% filter(group == "C")))
pg_detection_data_expression %<>% filter(group== "C")
pg_detection_data_expression %<>% mutate(pg_detected_text = if_else(pg_detected, "detected", "not detected"))
pg_detection_data_expression %<>% mutate(pg_detected_text = factor(pg_detected_text, levels = c("not detected", "detected")))
SRP_plot <- ggplot(data= pg_detection_data_expression, aes(y= expression, x = pg_detected_text)) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = .5) + theme_bw() + 
  ylab(paste0(target_gene, " expression (counts)")) + xlab("P. gingivalis")
SRP_plot

#ROC curve for SRP genes
gene_wise_results <- read_csv(here("results","GSE68719","genewise_pg_pd_controls.csv"))
gene_wise_results %<>% select(gene_symbol, rank)
srp_genes <- read_csv(here("data","gene_lists","SRP_list.txt"), col_names = F) %>% pull(X1)

gene_wise_results %<>% mutate(present = gene_symbol %in% c("ITM2B","TMCO1","SRP9","EIF3M","RAB28","ASF1A","IMMP1L","MNAT1","ZNF267","MICU2"))

gene_wise_results %<>% mutate(present = gene_symbol %in% srp_genes)


(AUCPlot <- ggplot(gene_wise_results, aes(d = present, m = rank)) + ylab("") + 
    style_roc() + coord_cartesian(expand=F) +
    geom_roc(n.cuts=0, linetype = "solid") + 
    geom_abline(slope = 1, intercept = 0, colour = "grey", size = 0.2) +
    labs(color='Gene Group')  + 
    theme(strip.background = element_blank(), strip.placement = "inside", strip.text = element_blank()) +
    #theme(legend.key.width = unit(2, "line"), legend.position = c(1,0), legend.justification = c(1, 0), legend.background= element_rect(fill = "transparent", colour = "transparent"), plot.margin=unit(c(.5,.5,.5,.5),"cm")) +
    theme(legend.position = "none") 
)

plot_grid(SRP_plot, AUCPlot, labels = c("A", "B"), scale = 0.95)
#save as 9x5 pdf
