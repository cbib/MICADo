library(ggplot2)
library(data.table)
library(dplyr)
library(gridExtra)
library(reshape2)
library(mosaic)
setwd("~/Documents/MICADo")

# Aggregated results 
# obtained by running bin/tabulate_and_aggregate_xp_results.py

micado_result = fread("data/synthetic/summary/agg_micado_results_on_synthetic_data.csv") 
gatk_data = fread("data/synthetic/summary/agg_gatk_results_on_synthetic_data.csv") 
varscan_data = fread("data/synthetic/summary/agg_varscan_results_on_synthetic_data.csv") 

combined_results = rbind(micado_result,varscan_data ,gatk_data) %>% filter(n_reads!=1500)
table(micado_result$n_alterations)
table(combined_results$n_reads)
combined_results[is.na(combined_results)] <- 0



agg_results_discrete = combined_results %>% 
  filter(n_reads!=800,fraction_altered!=0.07,fraction_altered!=0.075,fraction_altered!=0.8) %>% 
  transform(result_class = as.character(interaction(ifelse(tp>0,"TP+","TP-"),
  ifelse(fp==0,"FP-","FP+"), 
  ifelse(fn==0,"FN-","FN+"), sep=" ")))

agg_results_discrete=agg_results_discrete %>% transform(result_class=ifelse(tp==n_alterations & fp==0 & fn==0,"Perfect",result_class))
agg_results_discrete$result_class = factor(agg_results_discrete$result_class,
                                           levels=c("TP- FP+ FN+", "TP- FP- FN+", "TP- FP+ FN-", 
                                                    "TP- FP- FN-", "TP+ FP+ FN+", "TP+ FP- FN+", 
                                                    "TP+ FP+ FN-", "TP+ FP- FN-","Perfect"))
agg_results_discrete$caller = factor(agg_results_discrete$caller,levels=c("gatk","micado","varscan"),labels=c("GATK","MICADo","VarScan"))
g=ggplot(agg_results_discrete,aes(x=factor(fraction_altered),fill=result_class))+geom_histogram(position="fill")+facet_grid(n_reads~caller)+scale_fill_brewer(palette=1,type="seq",name="Class")
g=g+scale_x_discrete(name="Fraction of altered reads (%)",labels=c("3.5%","4%","4.5%","5%","10%","50%","80%"))+theme(legend.position="top")+scale_y_continuous(name="Proportion of results",labels=c("0%","25%","50%","75%","100%"))

ggsave(plot=g,filename="manuscript/figures/performances_on_synthetic_data.pdf",w=11.57,h=7.85)

## Detailed results 
#
# data file obtained by runnig bin/post_process_results.py
micado_detailed_results=fread("data/synthetic/summary/micado_results_on_synthetic_data.csv")

#typing
micado_detailed_results$git_revision_hash=factor(micado_detailed_results$git_revision_hash)
micado_detailed_results$injected_alt_type=factor(micado_detailed_results$injected_alt_type)
micado_detailed_results$injected_hash=factor(micado_detailed_results$injected_hash)
micado_detailed_results$injected_len=factor(micado_detailed_results$injected_len)
micado_detailed_results$n_reads=factor(micado_detailed_results$n_reads)
micado_detailed_results$timestamp=factor(micado_detailed_results$timestamp)
micado_detailed_results$tool_sampler_alt_weight=factor(micado_detailed_results$tool_sampler_alt_weight)
micado_detailed_results$seed=factor(micado_detailed_results$seed)
micado_detailed_results$tool_sampler_fraction_altered=factor(micado_detailed_results$tool_sampler_fraction_altered)

micado_detailed_results=micado_detailed_results %>% filter(n_reads!=499) # 4 outliers
micado_detailed_results=micado_detailed_results %>% filter(injected_pos<=1000) # off target alignments to remove

micado_detailed_results=micado_detailed_results %>% mutate(class = derivedFactor(
  "tp" = is_match,
  "fn" = is.na(micado_hash),
  "fp" = is.na(injected_hash)
),.ordered=TRUE,.sort='given')

fraction_label = as.numeric(as.character(unique(micado_results$tool_sampler_fraction_altered)))*100

ggplot(micado_detailed_results,aes(x=tool_sampler_fraction_altered,fill=class))+geom_bar(position="fill")+
  facet_grid(n_reads~injected_alt_type,scale="free")+ 
  scale_x_discrete(labels=fraction_label)+
  xlab("Fraction of reads altered (%)")+
  ylab("Proportion")+
  scale_fill_manual(values=c("forestgreen","firebrick3"))

ggsave("manuscript/figures/match_by_alteration_type.pdf",w=11.57,h=5.85)



