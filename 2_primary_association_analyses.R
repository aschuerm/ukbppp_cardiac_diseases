### 1 - Import Libraries and data frames ####
  # 1a - Libraries ####

library(data.table)
library(boot)
library(broom)
library(dplyr)
library(survival)

  # 1b - Load dataframe and make new dataframes with time-varying covariates

df_inc <- fread("/medpop/esp2/aschuerm/ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz")
df_inc$Sex_numeric <- factor(df_inc$Sex_numeric)
df_inc$mergedrace <- factor(df_inc$mergedrace)
df_inc$ever_smoked <- factor(df_inc$ever_smoked)
df_inc$antihtnbase <- factor(df_inc$antihtnbase)
df_inc$cholmed <- factor(df_inc$cholmed)

df_cad <- df_inc
 df_cad$afib_fu <- ifelse(df_cad$afib_prev==1, 0, ifelse(df_cad$afib_inc==1, df_cad$afib_fu, NA))
 df_cad$hfail_fu <- ifelse(df_cad$hfail_prev==1, 0, ifelse(df_cad$hfail_inc==1, df_cad$hfail_fu, NA))
 df_cad$ao_sten_fu <- ifelse(df_cad$ao_sten_prev==1, 0, ifelse(df_cad$ao_sten_inc==1, df_cad$ao_sten_fu, NA))
 df_cad <- tmerge(data1=df_cad, data2=df_cad, id=id, tstop=cad_fu, event=event(cad_fu, cad_inc), afib=tdc(afib_fu), hfail=tdc(hfail_fu), ao_sten=tdc(ao_sten_fu))
df_afib <- df_inc
 df_afib$cad_fu <- ifelse(df_afib$cad_prev==1, 0, ifelse(df_afib$cad_inc==1, df_afib$cad_fu, NA))
 df_afib$hfail_fu <- ifelse(df_afib$hfail_prev==1, 0, ifelse(df_afib$hfail_inc==1, df_afib$hfail_fu, NA))
 df_afib$ao_sten_fu <- ifelse(df_afib$ao_sten_prev==1, 0, ifelse(df_afib$ao_sten_inc==1, df_afib$ao_sten_fu, NA))
 df_afib <- tmerge(data1=df_afib, data2=df_afib, id=id, tstop=afib_fu, event=event(afib_fu, afib_inc), cad=tdc(cad_fu), hfail=tdc(hfail_fu), ao_sten=tdc(ao_sten_fu))
df_hfail <- df_inc
 df_hfail$cad_fu <- ifelse(df_hfail$cad_prev==1, 0, ifelse(df_hfail$cad_inc==1, df_hfail$cad_fu, NA))
 df_hfail$afib_fu <- ifelse(df_hfail$afib_prev==1, 0, ifelse(df_hfail$afib_inc==1, df_hfail$afib_fu, NA))
 df_hfail$ao_sten_fu <- ifelse(df_hfail$ao_sten_prev==1, 0, ifelse(df_hfail$ao_sten_inc==1, df_hfail$ao_sten_fu, NA))
 df_hfail <- tmerge(data1=df_hfail, data2=df_hfail, id=id, tstop=hfail_fu, event=event(hfail_fu, hfail_inc), cad=tdc(cad_fu), afib=tdc(afib_fu), ao_sten=tdc(ao_sten_fu))
df_ao_sten <- df_inc
 df_ao_sten$cad_fu <- ifelse(df_ao_sten$cad_prev==1, 0, ifelse(df_ao_sten$cad_inc==1, df_ao_sten$cad_fu, NA))
 df_ao_sten$afib_fu <- ifelse(df_ao_sten$afib_prev==1, 0, ifelse(df_ao_sten$afib_inc==1, df_ao_sten$afib_fu, NA))
 df_ao_sten$hfail_fu <- ifelse(df_ao_sten$hfail_prev==1, 0, ifelse(df_ao_sten$hfail_inc==1, df_ao_sten$hfail_fu, NA))
 df_ao_sten <- tmerge(data1=df_ao_sten, data2=df_ao_sten, id=id, tstop=ao_sten_fu, event=event(ao_sten_fu, ao_sten_inc), cad=tdc(cad_fu), afib=tdc(afib_fu), hfail=tdc(hfail_fu))


### 5 - Create adjusted Cox models for each outcome across all proteins ####
  # 5a - Coronary artery disease adjusted model ####

cox_cad_adj <- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x) + age + age2 + 
                                          Sex_numeric + mergedrace + PC1 + PC2 + PC3 + 
                                          PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                          ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                          tchol_final + hdl_final + cholmed + 
                                          dm2_prev + tdi_log_final + creat_final + afib + hfail + ao_sten, data = df_cad)}, 
                                 exponentiate = TRUE, conf.int = TRUE)[1,]
result_cad_adj <- lapply(proteins, cox_cad_adj)
hr_cad_adj <- do.call(rbind, result_cad_adj)
hr_cad_adj$Outcome <- 'CAD'
hr_cad_adj$term <- proteins
rm(result_cad_adj, cox_cad_adj)

  # 5b - Atrial fibrillation adjusted model ####

cox_afib_adj <- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x) + age + age2 + 
                                         Sex_numeric + mergedrace + PC1 + PC2 + PC3 + 
                                         PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                         ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                         tchol_final + hdl_final + cholmed + 
                                         dm2_prev + tdi_log_final + creat_final + cad + hfail + ao_sten, data = df_afib)}, 
                                exponentiate = TRUE, conf.int = TRUE)[1,]
result_afib_adj <- lapply(proteins, cox_afib_adj)
hr_afib_adj <- do.call(rbind, result_afib_adj)
hr_afib_adj$Outcome <- 'Afib'
hr_afib_adj$term <- proteins
rm(result_afib_adj, cox_afib_adj)

# 5c - Heart failure adjusted model ####

cox_hfail_adj <- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x) + age + age2 + 
                                         Sex_numeric + mergedrace + PC1 + PC2 + PC3 + 
                                         PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                         ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                         tchol_final + hdl_final + cholmed + 
                                         dm2_prev + tdi_log_final + creat_final + cad + afib + ao_sten, data = df_hfail)}, 
                                exponentiate = TRUE, conf.int = TRUE)[1,]
result_hfail_adj <- lapply(proteins, cox_hfail_adj)
hr_hfail_adj <- do.call(rbind, result_hfail_adj)
hr_hfail_adj$Outcome <- 'HF'
hr_hfail_adj$term <- proteins
rm(result_hfail_adj, cox_hfail_adj)

# 5d - Aortic stenosis adjusted models ####

cox_ao_sten_adj <- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x) + age + age2 + 
                                         Sex_numeric + mergedrace + PC1 + PC2 + PC3 + 
                                         PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                         ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                         tchol_final + hdl_final + cholmed + 
                                         dm2_prev + tdi_log_final + creat_final + cad + afib + hfail, data = df_ao_sten)}, 
                                exponentiate = TRUE, conf.int = TRUE)[1,]
result_ao_sten_adj <- lapply(proteins, cox_ao_sten_adj)
hr_ao_sten_adj <- do.call(rbind, result_ao_sten_adj)
hr_ao_sten_adj$Outcome <- 'AS'
hr_ao_sten_adj$term <- proteins
rm(result_ao_sten_adj, cox_ao_sten_adj)


### 6 - Merge data frames, restructure and write to csv ####
  # 6a - Combine summary stats for base models of all outcomes to one data frame ####

adj_dfs <- list(hr_afib_adj, hr_cad_adj, hr_hfail_adj, hr_ao_sten_adj)
adj_df <- do.call("rbind", adj_dfs)
rm(hr_cad_base, hr_afib_adj, hr_hfail_adj, hr_ao_sten_adj, adj_dfs., 
   df_cad, df_afib, df_hfail, df_ao_sten)

  # 6b - Restructure data frame - remove unnecessary columns, rename, reorder ####

adj_df <- select(adj_df, -c('std.error', 'statistic'))
colnames(adj_df) <- c('Protein', 'HR', 'P_Value', 'CI_Lower', 'CI_Upper', 'Outcome')
adj_df <- adj_df[,c('Protein', 'Outcome', 'HR', 'CI_Lower', 'CI_Upper', 'P_Value')]

  # 6c - Write to csv ####

write.csv(adj_df, '/medpop/esp2/aschuerm/ukb_proteomics_cvd/output_files/1_adj_no_prev_timevar_model_hr.csv')

