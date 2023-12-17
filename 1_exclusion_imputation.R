### 1 - Libraries and data.frames ####
  # 1a - Libraries ####

  library(data.table)
  library(ukbtools)
  library(dplyr)
  library(R.utils)


  # 1b - Data.frame ####

  a <- fread(".../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_noexcl_noimput.tsv.gz")
  
  
### 2 - Exclude biomarkers + proteins ####
  
  index <- c(which(colnames(a)=="id"), which(colnames(a)=="CLIP2"):which(colnames(a)=="NPM1")) 
  b <- a[, ..index]
  rm(index)
  
  # dim(b)          
  # [1] 52705  1464 

  
  # 2a - Exclude biomarkers with >10% missingness ####
  
  b <- b %>% select(where(~mean(is.na(.)) <= 0.1))
  
  # dim(b)                                                                                                
  # [1] 52705  1460                                                                                       

  
  # 2b - Exclude individuals with >10% biomarker missingness ####

  b <- b[rowMeans(is.na(b[,-1])) <= 0.1,]
  
  # dim(b)                                                                                                
  # [1] 48674  1460                                                                                       

  
  a <- a[a$id %in% b$id, -c("PCOLCE","TACSTD2","CTSS","NPM1")]
  rm(b)
  

  # 2c - Exclude those without baseline covariates (missing self-reported race or principal components) ####
  
  a <- a[!is.na(a$race),]
  # dim(a)                                                                                                
  # [1] 48442  1563                                                                                       

  a <- a[!is.na(a$PC1),]
  # dim(a)                                                                                                
  # [1] 48094  1564                                                                                       
  
  
  # 2d - Exclude relatives ####

  b <- fread('..../UKBiobank/v3data/ukb7089_rel_s488363.dat') # this is the UKB document that contains info on genetically inferred kinship
  set.seed(1)
  c <- ukb_gen_samples_to_remove(b, a$id, cutoff = 0.0884)                    
  
  a <- a[!(a$id %in% c), ] 			                                              
  rm(b, c)
  # dim(a)                                                                                                
  # [1] 47665  1564                                                                                       
  
  
  # 2e - Exclude participants with prevalent disease ####

  a <- a[a$cad_prev==0,]
  # dim(a)                                                                                                
  # [1] 45258  1564                                                                                       

  a <- a[a$hfail_prev==0,]
  # dim(a)                                                                                                
  # [1] 45053  1564                                                                                       

  a <- a[a$afib_prev==0,]
  # dim(a)                                                                                                
  # [1] 44372  1564                                                                                       

  a <- a[a$ao_sten_prev==0,]
  # dim(a)                                                                                                
  # [1] 44313  1564                                                                                       

  notimputed <- a
  
    
### 3 - Imputation of covariates ####

  library(RNOmni)
  b <- a[!is.na(a$tdi), c("id", "tdi")]
  b$tdi_log <- scale(RankNorm(b$tdi))
  a <- merge(a, b[,c("id", "tdi_log")], by="id", all.x=T)
  # dim(a)                                                                                                
  # [1] 44313  1565                                                                                       
  
  
  # 3a - Make missingness variables  ####

  a$BMI_missing <- ifelse(is.na(a$BMI), 1, 0)
  a$tdi_log_missing <- ifelse(is.na(a$tdi_log), 1, 0)
  a$alc_int_num_missing <- ifelse(is.na(a$alc_int_num), 1, 0)
  a$SBP_missing <- ifelse(is.na(a$SBP), 1, 0)
  a$DBP_missing <- ifelse(is.na(a$DBP), 1, 0)
  a$tchol_missing <- ifelse(is.na(a$tchol), 1, 0)
  a$ldl_missing <- ifelse(is.na(a$ldl), 1, 0)
  a$hdl_missing <- ifelse(is.na(a$hdl), 1, 0)
  a$tg_missing <- ifelse(is.na(a$tg), 1, 0)
  a$creat_missing <- ifelse(is.na(a$creat), 1, 0)

  
  # 3b - Construct linear regression models to predict the variables that are to be imputed  ####

  BMI_predict <- lm(BMI ~ Sex_numeric + age + mergedrace + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = a)
  tdi_log_predict <- lm(tdi_log ~ Sex_numeric + age + mergedrace + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = a)
  alc_int_num_predict <- lm(alc_int_num ~ Sex_numeric + age + mergedrace + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = a)
  SBP_predict <- lm(SBP ~ Sex_numeric + age + mergedrace + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = a)
  DBP_predict <- lm(DBP ~ Sex_numeric + age + mergedrace + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = a)
  tchol_predict <- lm(tchol ~ Sex_numeric + age + mergedrace + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = a)
  ldl_predict <- lm(ldl ~ Sex_numeric + age + mergedrace + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = a)
  hdl_predict <- lm(hdl ~ Sex_numeric + age + mergedrace + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = a)
  tg_predict <- lm(tg ~ Sex_numeric + age + mergedrace + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = a)
  creat_predict <- lm(creat ~ Sex_numeric + age + mergedrace + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data = a)

  
  # 3c - Use the predict() function to predict variables based on the previous linear regression models  ####

  a$BMI_predicted <- predict(BMI_predict, a[,c("Sex_numeric", "age", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])
  a$tdi_log_predicted <- predict(tdi_log_predict, a[,c("Sex_numeric", "age", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])
  a$alc_int_num_predicted <- predict(alc_int_num_predict, a[,c("Sex_numeric", "age", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])
  a$SBP_predicted <- predict(SBP_predict, a[,c("Sex_numeric", "age", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])
  a$DBP_predicted <- predict(DBP_predict, a[,c("Sex_numeric", "age", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])
  a$tchol_predicted <- predict(tchol_predict, a[,c("Sex_numeric", "age", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])
  a$ldl_predicted <- predict(ldl_predict, a[,c("Sex_numeric", "age", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])
  a$hdl_predicted <- predict(hdl_predict, a[,c("Sex_numeric", "age", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])
  a$tg_predicted <- predict(tg_predict, a[,c("Sex_numeric", "age", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])
  a$creat_predicted <- predict(creat_predict, a[,c("Sex_numeric", "age", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")])

  
  # 3d - Impute missing covariates  ####

  a$BMI_final <- ifelse(is.na(a$BMI), a$BMI_predicted, a$BMI)
  a$tdi_log_final <- ifelse(is.na(a$tdi_log), a$tdi_log_predicted, a$tdi_log)
  a$alc_int_num_final <- ifelse(is.na(a$alc_int_num), a$alc_int_num_predicted, a$alc_int_num)
  a$SBP_final <- ifelse(is.na(a$SBP), a$SBP_predicted, a$SBP)
  a$DBP_final <- ifelse(is.na(a$DBP), a$DBP_predicted, a$DBP)
  a$tchol_final <- ifelse(is.na(a$tchol), a$tchol_predicted, a$tchol)
  a$ldl_final <- ifelse(is.na(a$ldl), a$ldl_predicted, a$ldl)
  a$hdl_final <- ifelse(is.na(a$hdl), a$hdl_predicted, a$hdl)
  a$tg_final <- ifelse(is.na(a$tg), a$tg_predicted, a$tg)
  a$creat_final <- ifelse(is.na(a$creat), a$creat_predicted, a$creat)

  
### 4 - Imputation and rank-base INT of biomarkers ####
  
  library(impute)
  
  proteomics_imputed <- impute.knn(as.matrix(a[,which(colnames(a)=="CLIP2"):which(colnames(a)=="SCARB2")]), rowmax = 0.1, colmax = 0.1, k = 10, rng.seed=1)
  a[,which(colnames(a)=="CLIP2"):which(colnames(a)=="SCARB2")] <- as.data.frame(proteomics_imputed$data)
  # dim(a)
  # [1] 44313  1595
  
  for (i in (which(colnames(a)=="CLIP2"):which(colnames(a)=="SCARB2"))) {
   a[,colnames(a)[i]] <- scale(as.numeric(as.matrix(a[,..i])))
  }
  # dim(a)
  # [1] 44313  1595
  
    
### 5 - Add plate / well variables  ####
  
  platewell <- fread('..../UK_Biobank/baskets/ukb672544/ukb672544.tab.gz', select=c("f.eid", "f.30901.0.0", "f.30902.0.0"))
  colnames(platewell) <- c("id", "plate", "well")
  a <- merge(a, platewell, by="id", all.x=T, all.y=F)
  rm(platewell)
  
  
### 6 - Write table  ####
  
  write.table(a, ".../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv", sep="\t", row.names=F, col.names=T, quote=F)
  gzip(".../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv")

  
### 7 - Quantify missingness  ####

  full <- fread(".../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_noexcl_noimput.tsv.gz")
  notimputed <- merge(notimputed, full[,c("id", "PCOLCE", "CTSS", "NPM1", "TACSTD2")], by="id", all.x=T, all.y=F)
  
  df_impcount <- data.frame(colname=names(colSums(is.na(notimputed))), missing_count=colSums(is.na(notimputed)))
  rownames(df_impcount) <- NULL
  df_impcount$missing_proportion <- df_impcount$missing_count / nrow(notimputed)
  write.csv(df_impcount, ".../ukb_proteomics_cvd/input_files/missingness_basefile.csv", row.names=F)
