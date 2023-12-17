### 1 - Import Libraries and data frames ####
  # 1a - Libraries ####

library(data.table)
library(boot)
library(broom)
library(dplyr)
library(survival)

  # 1b - Load dataframe ####

df <- fread(".../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz")

### 2 - Factorization and incident disease exclusion ####
  # 2a - Factorization of categorical variables ####
df$Sex_numeric <- factor(df$Sex_numeric)
df$mergedrace <- factor(df$mergedrace)
df$ever_smoked <- factor(df$ever_smoked)
df$antihtnbase <- factor(df$antihtnbase)
df$cholmed <- factor(df$cholmed)

  # 2b - Make a list of all proteins ####
proteins <- c(colnames(df)[which(colnames(df)=="CLIP2"):which(colnames(df)=="SCARB2")])



df_inc <- df %>% filter(!is.na(cad_inc) & !is.na(afib_inc) & !is.na(hfail_inc) 
                                                   & !is.na(ao_sten_inc))
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

 
### 5 - Create adjusted cox PH models for each outcome across all proteins including sex interaction term ####
  # 5a - CAD adjusted model ####

cox_cad_adj <- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x)*Sex_numeric + age + age2 + 
                                          mergedrace + PC1 + PC2 + PC3 + 
                                          PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                          ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                          tchol_final + hdl_final + cholmed + 
                                          dm2_prev + tdi_log_final + creat_final + afib + hfail + ao_sten, data = df_cad)}, 
                                 exponentiate = TRUE, conf.int = TRUE)[29,]
result_cad_adj <- lapply(proteins, cox_cad_adj)
hr_cad_adj <- do.call(rbind, result_cad_adj)
hr_cad_adj$Outcome <- 'CAD'
hr_cad_adj$term <- proteins
rm(result_cad_adj, cox_cad_adj)


  # 5b - Atrial fibrillation adjusted model ####

cox_afib_adj <- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x)*Sex_numeric + age + age2 + 
                                         mergedrace + PC1 + PC2 + PC3 + 
                                         PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                         ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                         tchol_final + hdl_final + cholmed + 
                                         dm2_prev + tdi_log_final + creat_final + cad + hfail + ao_sten, data = df_afib)}, 
                                exponentiate = TRUE, conf.int = TRUE)[29,]
result_afib_adj <- lapply(proteins, cox_afib_adj)
hr_afib_adj <- do.call(rbind, result_afib_adj)
hr_afib_adj$Outcome <- 'Afib'
hr_afib_adj$term <- proteins
rm(result_afib_adj, cox_afib_adj)


  # 5c - Heart failure adjusted model ####

cox_hfail_adj <- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x)*Sex_numeric + age + age2 + 
                                         mergedrace + PC1 + PC2 + PC3 + 
                                         PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                         ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                         tchol_final + hdl_final + cholmed + 
                                         dm2_prev + tdi_log_final + creat_final + cad + afib + ao_sten, data = df_hfail)}, 
                                exponentiate = TRUE, conf.int = TRUE)[29,]
result_hfail_adj <- lapply(proteins, cox_hfail_adj)
hr_hfail_adj <- do.call(rbind, result_hfail_adj)
hr_hfail_adj$Outcome <- 'HF'
hr_hfail_adj$term <- proteins
rm(result_hfail_adj, cox_hfail_adj)


  # 5d - Aortic stenosis adjusted models ####

cox_ao_sten_adj <- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x)*Sex_numeric + age + age2 + 
                                         mergedrace + PC1 + PC2 + PC3 + 
                                         PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                         ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                         tchol_final + hdl_final + cholmed + 
                                         dm2_prev + tdi_log_final + creat_final + cad + afib + hfail, data = df_ao_sten)}, 
                                exponentiate = TRUE, conf.int = TRUE)[29,]
result_ao_sten_adj <- lapply(proteins, cox_ao_sten_adj)
hr_ao_sten_adj <- do.call(rbind, result_ao_sten_adj)
hr_ao_sten_adj$Outcome <- 'AS'
hr_ao_sten_adj$term <- proteins
rm(result_ao_sten_adj, cox_ao_sten_adj)


#### 7 - Run sex-stratified Cox models for significant proteins ####
  # 7a - Subgroup into male and female dateframes for each outcome ####

df_cad_m <- df_cad %>% filter(df_cad$Sex_numeric ==1)
df_cad_f <- df_cad %>% filter(df_cad$Sex_numeric ==0)
df_afib_m <- df_afib %>% filter(df_afib$Sex_numeric ==1)
df_afib_f <- df_afib %>% filter(df_afib$Sex_numeric ==0)
df_hfail_m <- df_hfail %>% filter(df_hfail$Sex_numeric ==1)
df_hfail_f <- df_hfail %>% filter(df_hfail$Sex_numeric ==0)
df_as_m <- df_ao_sten %>% filter(df_ao_sten$Sex_numeric ==1)
df_as_f <- df_ao_sten %>% filter(df_ao_sten$Sex_numeric ==0)


  # 7c - Run CAD models stratified by male/female ####

cox_cad_adj_m <- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x) + Sex_numeric + age + age2 + 
                                         mergedrace + PC1 + PC2 + PC3 + 
                                         PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                         ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                         tchol_final + hdl_final + cholmed + 
                                         dm2_prev + tdi_log_final + creat_final + afib + hfail + ao_sten, data = df_cad_m)}, 
                                exponentiate = TRUE, conf.int = TRUE)[1,]
result_cad_adj_m <- lapply(proteins, cox_cad_adj_m)
hr_cad_adj_m <- do.call(rbind, result_cad_adj_m)
hr_cad_adj_m$Outcome <- 'CAD'
hr_cad_adj_m$term <- proteins
colnames(hr_cad_adj_m) <- c('Protein', 'HR_M', 'SE_M', 'Stat_M', 'P_Val_M', 'CI_Low_M', 'CI_High_M', 'Outcome')
rm(result_cad_adj_m, cox_cad_adj_m)

cox_cad_adj_f<- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x) + Sex_numeric + age + age2 + 
                                           mergedrace + PC1 + PC2 + PC3 + 
                                           PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                           ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                           tchol_final + hdl_final + cholmed + 
                                           dm2_prev + tdi_log_final + creat_final + afib + hfail + ao_sten, data = df_cad_f)}, 
                                  exponentiate = TRUE, conf.int = TRUE)[1,]
result_cad_adj_f <- lapply(proteins, cox_cad_adj_f)
hr_cad_adj_f <- do.call(rbind, result_cad_adj_f)
hr_cad_adj_f$Outcome <- 'CAD'
hr_cad_adj_f$term <- proteins
colnames(hr_cad_adj_f) <- c('Protein', 'HR_F', 'SE_F', 'Stat_F', 'P_Val_F', 'CI_Low_F', 'CI_High_F', 'Outcome')
rm(result_cad_adj_f, cox_cad_adj_f)

cad_strat <- cbind(hr_cad_adj_m, hr_cad_adj_f)
cad_strat <- cbind(cad_strat, hr_cad_adj)


  # 7d - Run AF models stratified by male/female ####

cox_afib_adj_m <- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x) + Sex_numeric + age + age2 + 
                                           mergedrace + PC1 + PC2 + PC3 + 
                                           PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                           ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                           tchol_final + hdl_final + cholmed + 
                                           dm2_prev + tdi_log_final + creat_final + cad + hfail + ao_sten, data = df_afib_m)}, 
                                  exponentiate = TRUE, conf.int = TRUE)[1,]
result_afib_adj_m <- lapply(proteins, cox_afib_adj_m)
hr_afib_adj_m <- do.call(rbind, result_afib_adj_m)
hr_afib_adj_m$Outcome <- 'AF'
hr_afib_adj_m$term <- proteins
colnames(hr_afib_adj_m) <- c('Protein', 'HR_M', 'SE_M', 'Stat_M', 'P_Val_M', 'CI_Low_M', 'CI_High_M', 'Outcome')
rm(result_afib_adj_m, cox_afib_adj_m)

cox_afib_adj_f<- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x) + Sex_numeric + age + age2 + 
                                          mergedrace + PC1 + PC2 + PC3 + 
                                          PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                          ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                          tchol_final + hdl_final + cholmed + 
                                          dm2_prev + tdi_log_final + creat_final + cad + hfail + ao_sten, data = df_afib_f)}, 
                                 exponentiate = TRUE, conf.int = TRUE)[1,]
result_afib_adj_f <- lapply(proteins, cox_afib_adj_f)
hr_afib_adj_f <- do.call(rbind, result_afib_adj_f)
hr_afib_adj_f$Outcome <- 'AF'
hr_afib_adj_f$term <- proteins
colnames(hr_afib_adj_f) <- c('Protein', 'HR_F', 'SE_F', 'Stat_F', 'P_Val_F', 'CI_Low_F', 'CI_High_F', 'Outcome')
rm(result_afib_adj_f, cox_afib_adj_f)

afib_strat <- cbind(hr_afib_adj_m, hr_afib_adj_f)
afib_strat <- cbind(afib_strat, hr_afib_adj)


  # 7e - Run HF Models stratified by male/female ####

cox_hfail_adj_m <- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x) + Sex_numeric + age + age2 + 
                                            mergedrace + PC1 + PC2 + PC3 + 
                                            PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                            ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                            tchol_final + hdl_final + cholmed + 
                                            dm2_prev + tdi_log_final + creat_final + cad + afib + ao_sten, data = df_hfail_m)}, 
                                   exponentiate = TRUE, conf.int = TRUE)[1,]
result_hfail_adj_m <- lapply(proteins, cox_hfail_adj_m)
hr_hfail_adj_m <- do.call(rbind, result_hfail_adj_m)
hr_hfail_adj_m$Outcome <- 'HF'
hr_hfail_adj_m$term <- proteins
colnames(hr_hfail_adj_m) <- c('Protein', 'HR_M', 'SE_M', 'Stat_M', 'P_Val_M', 'CI_Low_M', 'CI_High_M', 'Outcome')
rm(result_hfail_adj_m, cox_hfail_adj_m)

cox_hfail_adj_f<- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x) + Sex_numeric + age + age2 + 
                                           mergedrace + PC1 + PC2 + PC3 + 
                                           PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                           ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                           tchol_final + hdl_final + cholmed + 
                                           dm2_prev + tdi_log_final + creat_final + cad + afib + ao_sten, data = df_hfail_f)}, 
                                  exponentiate = TRUE, conf.int = TRUE)[1,]
result_hfail_adj_f <- lapply(proteins, cox_hfail_adj_f)
hr_hfail_adj_f <- do.call(rbind, result_hfail_adj_f)
hr_hfail_adj_f$Outcome <- 'HF'
hr_hfail_adj_f$term <- proteins
colnames(hr_hfail_adj_f) <- c('Protein', 'HR_F', 'SE_F', 'Stat_F', 'P_Val_F', 'CI_Low_F', 'CI_High_F', 'Outcome')
rm(result_hfail_adj_f, cox_hfail_adj_f)

hfail_strat <- cbind(hr_hfail_adj_m, hr_hfail_adj_f)
hfail_strat <- cbind(hfail_strat, hr_hfail_adj)


  # 7f - Run AS models stratified by male/female ####

cox_as_adj_m <- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x) + Sex_numeric + age + age2 + 
                                             mergedrace + PC1 + PC2 + PC3 + 
                                             PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                             ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                             tchol_final + hdl_final + cholmed + 
                                             dm2_prev + tdi_log_final + creat_final + cad + afib + hfail, data = df_as_m)}, 
                                    exponentiate = TRUE, conf.int = TRUE)[1,]
result_as_adj_m <- lapply(proteins, cox_as_adj_m)
hr_as_adj_m <- do.call(rbind, result_as_adj_m)
hr_as_adj_m$Outcome <- 'AS'
hr_as_adj_m$term <- proteins
colnames(hr_as_adj_m) <- c('Protein', 'HR_M', 'SE_M', 'Stat_M', 'P_Val_M', 'CI_Low_M', 'CI_High_M', 'Outcome')
rm(result_as_adj_m, cox_as_adj_m)

cox_as_adj_f<- function(x) tidy({coxph(Surv(tstart, tstop, event) ~ get(x) + Sex_numeric + age + age2 + 
                                            mergedrace + PC1 + PC2 + PC3 + 
                                            PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                            ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                            tchol_final + hdl_final + cholmed + 
                                            dm2_prev + tdi_log_final + creat_final + cad + afib + hfail, data = df_as_f)}, 
                                   exponentiate = TRUE, conf.int = TRUE)[1,]
result_as_adj_f <- lapply(proteins, cox_as_adj_f)
hr_as_adj_f <- do.call(rbind, result_as_adj_f)
hr_as_adj_f$Outcome <- 'AS'
hr_as_adj_f$term <- proteins
colnames(hr_as_adj_f) <- c('Protein', 'HR_F', 'SE_F', 'Stat_F', 'P_Val_F', 'CI_Low_F', 'CI_High_F', 'Outcome')
rm(result_as_adj_f, cox_as_adj_f)

as_strat <- cbind(hr_as_adj_m, hr_as_adj_f)
as_strat <- cbind(as_strat, hr_ao_sten_adj)


### 8 - Combine and restructure dataframe ####

strat_full <- rbind(cad_strat, afib_strat)
strat_full <- rbind(strat_full, hfail_strat)
strat_full <- rbind(strat_full, as_strat)
strat_full <- strat_full[,-c(9,16,17,20,24)]
colnames(strat_full) <- c(colnames(strat_full)[-(15:19)], "HR_Int", "SE_Int", "P_Val_Int", "CI_Low_Int", "CI_High_Int")
write.csv(strat_full, '.../ukb_proteomics_cvd/output_files/4_sex_strat_full.csv')
