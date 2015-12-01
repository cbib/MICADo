# libs
library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(gridExtra)
library(reshape2)
library(xtable)
library(hash)
library(mosaic)
setwd("~/Documents/MICADo")

# sample labels 
pool_0_labels = fread("data/experimental_results/TP53/pool_0_groups.tsv")
known_snps_pos=c(417,664,588,768,841)

# Known mutations 
known_mutations = fread("data/experimental_results/TP53/MutOnPool0.tsv")
setnames(known_mutations,c("patient_id","patient_biopsy","alt_description","fragment"))
# only keep rows with mutations 
known_mutations=known_mutations %>% filter(fragment!="/")

known_mutations=known_mutations%>% mutate(alt_type = derivedFactor(
  "I"=grepl("ins",alt_description,fixed=T),
  "X"=grepl(">",alt_description,fixed=T),
  "D"=grepl("del",alt_description,fixed=T)
))


# decompose start and end 
coordinates = stringr::str_split_fixed(stringr::str_extract_all(known_mutations$alt_description,"([0-9]+)_([0-9]+)|([0-9]+)"),"_",n=2)
content = stringr::str_split_fixed(stringr::str_extract_all(known_mutations$alt_description,"([ATCG]+)>?([ATCG]*)"),">",n=2)
known_mutations$start = as.numeric(as.character(coordinates[,1]))
known_mutations$end = as.numeric(as.character(coordinates[,2]))
known_mutations$ref=as.character(content[,1])
known_mutations$alt=as.character(content[,2])
known_mutations$sample=stringr::str_c(known_mutations$fragment,known_mutations$patient_id,known_mutations$patient_biopsy,sep="_")
known_mutations$sample_key=stringr::str_c(known_mutations$patient_id,known_mutations$patient_biopsy,sep="_")

known_mutations = known_mutations%>% mutate(uuid=stringr::str_c(sample,alt_type,start+2)) 
known_mutations %>% filter(sample=="N_276_1")

# function to allow for up to -5 or +5 start offset during matching 
augment_df_with_slack = function(orig_df){
  df_with_slack=orig_df %>% mutate(real_start=start)
  for(slack_offset in -5:+5){
    df_with_slack=rbind(df_with_slack,orig_df %>% mutate(real_start=start,start=start+slack_offset))
  }
  df_with_slack= df_with_slack %>% mutate(uuid=stringr::str_c(sample,alt_type,start+2)) 
  return(df_with_slack)
}

# return elements of src_df found in tgt_df, with some slack 
find_with_slack = function(src_df,tgt_df,negate=F){
  if (negate){
    src_df %>% filter(!(uuid %in% augment_df_with_slack(tgt_df)$uuid))
  }else{
    src_df %>% filter((uuid %in% augment_df_with_slack(tgt_df)$uuid))
  }
}
# build the uuid list 
# MICADo

pool_0_results_micado = fread("data/tp53_analysis/summary/agg_unsupervised_micado_results_on_pool0_data.csv")
pool_0_results_micado=pool_0_results_micado %>% mutate(uuid=stringr::str_c(sample,alt_type,start+2))
# 
# # weird bug on some deletion start pos 
# comparable_mutations = pool_0_results_micado  %>% 
#              mutate(start=start+2,micado_uid=1:nrow(pool_0_results_micado))  %>% 
#               left_join(known_mutations,by="sample")  %>% 
#               select(sample,fragment,contains("alt_type"),contains('start'),micado_uid,z_score,end.x,end.y)  %>% 
#               mutate(is_match=alt_type.x==alt_type.y & start.x==start.y)
# 
# # select the closest 
# comparable_mutations %>% head(30)
# comparable_mutations %>% filter(is_match==F) %>% arrange(-z_score) %>% head(30)
# comparable_mutations %>% filter(is_match==T) %>% arrange(-z_score) %>% head(30)
# represent each mutation by a comparable string hash 
micado_missed = known_mutations %>% filter(!(uuid %in% augment_df_with_slack(pool_0_results_micado)$uuid))
micado_fp = pool_0_results_micado %>% filter(!(uuid %in% augment_df_with_slack(known_mutations)$uuid)) %>% arrange(-z_score) %>% filter(!is.na(z_score))

# varscan

pool_0_results_varscan = fread("data/tp53_analysis/summary/agg_unsupervised_varscan_results_on_pool0_data.csv") %>% inner_join(pool_0_labels)
pool_0_results_varscan = pool_0_results_varscan %>% filter(!(start %in% known_snps_pos))
# manual correction for a particular tandem repeat in sample N_290_1/2, GCAGT is the ref seq, 
# GCAGT is also inserted => GCAGT two times in a row. Ref mutation indicates the rightmost, 
# VarScan the leftmost, we correct in favor of VarScan 
pool_0_results_varscan[pool_0_results_varscan$alt_sequence=="AGCAGT"]$start=pool_0_results_varscan[pool_0_results_varscan$alt_sequence=="AGCAGT"]$start+5
pool_0_results_varscan = pool_0_results_varscan %>% mutate(uuid=stringr::str_c(sample,alt_type,ifelse(alt_type=="I",start+alt_length-4,ifelse(alt_type=="D",start+1,start))))
varscan_missed = known_mutations %>% filter(!(uuid %in% augment_df_with_slack(pool_0_results_varscan)$uuid))
varscan_fp = pool_0_results_varscan %>% filter(!(uuid %in% augment_df_with_slack(known_mutations)$uuid)) 

# GATK 
pool_0_results_gatk = fread("data/tp53_analysis/summary/agg_unsupervised_gatk_results_on_pool0_data.csv") %>% inner_join(pool_0_labels)
pool_0_results_gatk = pool_0_results_gatk %>% filter(!(start %in% known_snps_pos))
pool_0_results_gatk = pool_0_results_gatk %>% mutate(uuid=stringr::str_c(sample,alt_type,ifelse(alt_type=="I",start+1,ifelse(alt_type=="D",start+nchar(alt_sequence)+2,start))))
gatk_missed = known_mutations %>% filter(!(uuid %in% augment_df_with_slack(pool_0_results_gatk)$uuid))
gatk_fp = pool_0_results_gatk %>% filter(!(uuid %in% augment_df_with_slack(known_mutations)$uuid))  %>% filter(!is.na(start))



######## Manuscript table 
# score and build table 

sample_category = hash(keys=stringr::str_extract(pool_0_labels$sample,"[0-9]+_[0-9]"),values=pool_0_labels$class)
stats_for_sample = function(a_sample){
  expected = known_mutations %>% filter(sample_key==a_sample)
  a_sample_micado  = pool_0_results_micado %>% filter(sample_key==a_sample,z_score>=35) %>% filter(!is.na(start))
  a_sample_micado_tp = find_with_slack(a_sample_micado,expected)
  a_sample_micado_fp = a_sample_micado %>% filter(sample_key==a_sample) %>% find_with_slack(expected,negate=T) %>% filter(!is.na(start))
  
  a_sample_varscan  = pool_0_results_varscan %>% filter(sample_key==a_sample) %>% filter(!is.na(start))
  a_sample_varscan_tp = find_with_slack(a_sample_varscan,expected)
  a_sample_varscan_fp = a_sample_varscan %>% filter(sample_key==a_sample) %>% find_with_slack(expected,negate=T) %>% filter(!is.na(start))
  
  a_sample_gatk  = pool_0_results_gatk %>% filter(sample_key==a_sample) %>% filter(!is.na(start))
  a_sample_gatk_tp = find_with_slack(a_sample_gatk,expected)
  a_sample_gatk_fp = a_sample_gatk %>% filter(sample_key==a_sample) %>% find_with_slack(expected,negate=T) %>% filter(!is.na(start))
  
  data.frame(a_sample=a_sample,
             a_sample_category = sample_category[[a_sample]],
             expected_mutations = nrow(expected),
             micado_tp_calls=nrow(a_sample_micado_tp),
             micado_fp_calls=nrow(a_sample_micado_fp),
             micado_tot_calls=nrow(a_sample_micado),
             varscan_tp_calls=nrow(a_sample_varscan_tp),
             varscan_fp_calls=nrow(a_sample_varscan_fp),
             varscan_tot_calls=nrow(a_sample_varscan),
             gatk_tp_calls=nrow(a_sample_gatk_tp),
             gatk_fp_calls=nrow(a_sample_gatk_fp),
             gatk_tot_calls=nrow(a_sample_gatk)
  )
}

stats_for_sample("169_1")
result_table = ldply(keys(sample_category),stats_for_sample)

# save as csv and as latex 
write.csv(result_table,file="manuscript/tables/tp53_pool_0_aggregate_by_caller.csv")
ms_table = result_table %>% arrange(a_sample_category,a_sample) %>% transmute(
  "EORTC ID"=a_sample,
  category=a_sample_category,
  "Exp. #"=expected_mutations,
  micado=paste(micado_tp_calls,micado_tot_calls,sep="/"),
  GATK=paste(gatk_tp_calls,gatk_tot_calls,sep="/"),
  VarScan=paste(varscan_tp_calls,varscan_tot_calls,sep="/")
  )
table_tex=print(xtable(ms_table))
fileConn<-file("manuscript/tables/tp53_pool_0_aggregate_by_caller.tex")
writeLines( table_tex, fileConn)
close(fileConn)

