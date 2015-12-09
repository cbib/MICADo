# libs
library(ggplot2)
library(data.table)
library(dplyr)
library(gridExtra)
library(reshape2)
library(mosaic)
setwd("~/Documents/MICADo")

# sample labels 
pool_0_labels = fread("data/experimental_results/TP53/pool_0_groups.tsv")
known_snps_pos=c(417,664,588,768,841)

# MICADo

pool_0_results_micado = fread("data/tp53_analysis/summary/agg_unsupervised_micado_results_on_pool0_data.csv")
pool_0_results_micado=pool_0_labels %>% inner_join(pool_0_results_micado)
grouped_results_micado = pool_0_results_micado %>% group_by(sample_key, class) %>% summarise(n_significant_alterations = sum(z_score >= 10, na.rm = T)) %>% ungroup()
grouped_results_micado %>% filter(sample_key=="319_1")


# VarScan 
pool_0_results_varscan = fread("data/tp53_analysis/summary/agg_unsupervised_varscan_results_on_pool0_data.csv") %>% inner_join(pool_0_labels)
pool_0_results_varscan = pool_0_results_varscan %>% filter(!(start %in% known_snps_pos))
grouped_results_varscan = pool_0_results_varscan %>% group_by(sample_key,class) %>% summarise(n_significant_alterations=n())
grouped_results_varscan=grouped_results_varscan %>% mutate(classification = derivedFactor(
  "perfect" = n_significant_alterations==0 & class=="negative control",
  "perfect" = n_significant_alterations>=1 & class=="positive control",
  "fn" = n_significant_alterations<1 & class=="positive control",
  "fp" = n_significant_alterations>0 & class=="negative control",
  "Unk"
),.ordered=TRUE,.sort='given')

# GATK 
pool_0_results_gatk = fread("data/tp53_analysis/summary/agg_unsupervised_gatk_results_on_pool0_data.csv") %>% inner_join(pool_0_labels)
pool_0_results_gatk = pool_0_results_gatk %>% filter(!(start %in% known_snps_pos))
grouped_results_gatk = pool_0_results_gatk %>% group_by(sample_key,class) %>% summarise(n_significant_alterations=sum(n_caller_alterations))


grouped_results_all_caller=rbind(
  grouped_results_gatk %>% mutate(caller="gatk"),
  rbind(
    grouped_results_varscan %>% mutate(caller="varscan"),
    grouped_results %>% mutate(caller="micado"),fill=TRUE
  ),fill=TRUE
)


grouped_results_all_caller=rbind(
  grouped_results_gatk %>% mutate(caller="gatk"),
  rbind(
    grouped_results_varscan %>% mutate(caller="varscan"),
    grouped_results_micado %>% mutate(caller="micado"),fill=TRUE
  ),fill=TRUE
)

g=grouped_results_all_caller %>% ggplot(aes(x=sample_key,y=n_significant_alterations,fill=caller,label=sample_key))+geom_bar(stat="identity")+facet_grid(caller~class,scale="free")
g=g+scale_y_continuous(breaks=seq(0,15,2),limits=c(0,14))
g=g+theme(
  axis.title.y = element_text(face="bold", colour="#990000", size=11),
  axis.title.x = element_text(face="bold", colour="#990000", size=11),
  axis.text.x  = element_text(angle=90, vjust=0.5, size=11))+xlab("Sample")+ylab("Number of significant alterations")
g=g+scale_fill_discrete(guide=FALSE)
print(g)
# ggsave(g,filename = "manuscript/figures/tp53_all_callers_results.pdf",w=11.57,h=5.85)
