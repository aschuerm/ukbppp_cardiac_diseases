### 1 - Load Libraries and upload input dataframe ####
  library(glmnet)
  library(tidyverse)
  library(dplyr)
  library(data.table)
  library(pROC)
  library(forcats)

  df <- fread('.../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz')

### 2 - Preprocess Dataframe to create predcitors and outcome ####
  # 2a - Shuffle Dataframe and split into training and test set ####
    set.seed(1234)
    df <- df[sample(1:nrow(df)), ]
    set.seed(1)
    df$id <- 1:nrow(df)
    train <- df %>% dplyr::sample_frac(0.80)
    test  <- dplyr::anti_join(df, train, by = 'id')
  
    
  # 2b - Split training and test sets into predictors and labels) ####
    
    df_res_allmodels <- data.frame(model=NA)[-1,]
    
    for (i in c("cad", "hfail", "afib", "ao_sten")) {
    
    df_pr_train <- train[,106:1564]
      df_pr_train <- data.matrix(df_pr_train)
      df_pr_test <- test[,106:1564]
      df_pr_test <- data.matrix(df_pr_test)
    df_cl_train <- train[,c(3,4,7,39,41,44,53,1586,1589,1591,1593, 1595)]
      df_cl_train <- data.matrix(df_cl_train)
      df_cl_test <- test[,c(3,4,7,39,41,44,53,1586,1589,1591,1593, 1595)]
      df_cl_test <- data.matrix(df_cl_test)
    df_comb_train <- train[,c(3,4,7,39,41,44,53,1586,1589,1591,1593, 1595, 106:1564)]
      df_comb_train <- data.matrix(df_comb_train)
      df_comb_test <- test[,c(3,4,7,39,41,44,53,1586,1589,1591,1593, 1595, 106:1564)]
      df_comb_test <- data.matrix(df_comb_test)
    index <- which(colnames(train)==paste0(i, "_inc"))
    y_train <- train[,..index]
      y_train <- data.matrix(y_train)
      y_test <- test[,..index]
      y_test <- data.matrix(y_test)
  
      
### 3 - Cross-validation and model fitting  ####
  # 3a - Perform cross validation (and visualize AUC/coeff) ####
      
    set.seed(1234)
    cv.lasso_pr <- cv.glmnet(df_pr_train, y_train, alpha = 1, family = "binomial", nfolds =10, type.measure = 'auc', keep = TRUE)
    set.seed(1234)
    cv.lasso_cl <- cv.glmnet(df_cl_train, y_train, alpha = 1, family = "binomial", nfolds =10, type.measure = 'auc', keep = TRUE)
    set.seed(1234)
    cv.lasso_comb <- cv.glmnet(df_comb_train, y_train, alpha = 1, family = "binomial", nfolds =10, type.measure = 'auc', keep = TRUE)
    
    pdf(paste0(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/1_lambda_featureselection_plot_", i, "_pr.pdf"))
    plot(cv.lasso_pr)
    dev.off()
    pdf(paste0(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/1_lambda_featureselection_plot_", i, "_cl.pdf"))
    plot(cv.lasso_cl)
    dev.off()
    pdf(paste0(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/1_lambda_featureselection_plot_", i, "_comb.pdf"))
    plot(cv.lasso_comb)
    dev.off()
    
    # Visualize coef of lambda 1SE (value of lambda that gives most regularized model such that CV error is within 1 SD of min)
    
    df_pr <- data.frame(lambda=cv.lasso_pr$lambda, cvm=cv.lasso_pr$cvm, cvsd=cv.lasso_pr$cvsd, 
                        cvup=cv.lasso_pr$cvup, cvlo=cv.lasso_pr$cvlo, nzero=cv.lasso_pr$nzero)
    df_cl <- data.frame(lambda=cv.lasso_cl$lambda, cvm=cv.lasso_cl$cvm, cvsd=cv.lasso_cl$cvsd, 
                        cvup=cv.lasso_cl$cvup, cvlo=cv.lasso_cl$cvlo, nzero=cv.lasso_cl$nzero)
    df_comb <- data.frame(lambda=cv.lasso_comb$lambda, cvm=cv.lasso_comb$cvm, cvsd=cv.lasso_comb$cvsd, 
                        cvup=cv.lasso_comb$cvup, cvlo=cv.lasso_comb$cvlo, nzero=cv.lasso_comb$nzero)
    
    write.csv(df_pr, paste0(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/1_lambda_featureselection_df_", i, "_pr.csv"))
    write.csv(df_cl, paste0(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/1_lambda_featureselection_df_", i, "_cl.csv"))
    write.csv(df_comb, paste0(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/1_lambda_featureselection_df_", i, "_comb.csv"))
    rm(df_pr, df_cl, df_comb)
    
    write.csv(as.matrix(coef(cv.lasso_pr, s= "lambda.1se")), paste0(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_", i, "_pr.csv"))
    write.csv(as.matrix(coef(cv.lasso_cl, s= "lambda.min")), paste0(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_", i, "_cl.csv"))
    write.csv(as.matrix(coef(cv.lasso_comb, s= "lambda.1se")), paste0(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_", i, "_comb.csv"))
  
  
  # 3b - Refit model with optimal lambda value ####
    
    model_pr <- glmnet(df_pr_train, y_train, alpha = 1, family = "binomial", lambda = cv.lasso_pr$lambda.1se)
    model_cl <- glmnet(df_cl_train, y_train, alpha = 1, family = "binomial", lambda = cv.lasso_cl$lambda.min)
    model_comb <- glmnet(df_comb_train, y_train, alpha = 1, family = "binomial", lambda = cv.lasso_comb$lambda.1se)
    
    
    # Make predictions 
    
    roc_pr <- roc.glmnet(model_pr, newx=df_pr_test, newy=as.numeric(y_test))
      roc_pr_predict <- predict(model_pr, newx=df_pr_test, s=cv.lasso_pr$lambda.1se, type="response")
      auc_model_pr <- roc(response=as.numeric(y_test[,1]), predictor=as.numeric(roc_pr_predict), levels=c(0,1), ci=T)$ci[1:3]
      auc_roc_pr <- data.frame(model="Proteins", outc=i, auc=auc_model_pr[2], auc_lci=auc_model_pr[1], auc_uci=auc_model_pr[3])
    
    roc_cl <- roc.glmnet(model_cl, newx=df_cl_test, newy=as.numeric(y_test))
      roc_cl_predict <- predict(model_cl, newx=df_cl_test, s=cv.lasso_cl$lambda.min, type="response")
      auc_model_cl <- roc(response=as.numeric(y_test[,1]), predictor=as.numeric(roc_cl_predict), levels=c(0,1), ci=T)$ci[1:3]
      auc_roc_cl <- data.frame(model="Clinical", outc=i, auc=auc_model_cl[2], auc_lci=auc_model_cl[1], auc_uci=auc_model_cl[3])
    
    roc_comb <- roc.glmnet(model_comb, newx=df_comb_test, newy=as.numeric(y_test))
      roc_comb_predict <- predict(model_comb, newx=df_comb_test, s=cv.lasso_comb$lambda.1se, type="response")
      auc_model_comb <- roc(response=as.numeric(y_test[,1]), predictor=as.numeric(roc_comb_predict), levels=c(0,1), ci=T)$ci[1:3]
      auc_roc_comb <- data.frame(model="Combination", outc=i, auc=auc_model_comb[2], auc_lci=auc_model_comb[1], auc_uci=auc_model_comb[3])
  
    write.csv(roc_pr, paste0(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_", i, "_pr.csv"))
      write.csv(roc_cl, paste0(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_", i, "_cl.csv"))
      write.csv(roc_comb, paste0(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_", i, "_comb.csv"))
  
    df_res_allmodels <- rbind(df_res_allmodels, auc_roc_pr, auc_roc_cl, auc_roc_comb)
      write.csv(df_res_allmodels, ".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_roc_full.csv")
  
    rm(roc_pr, roc_cl, roc_comb, auc_roc_pr, auc_roc_cl, auc_roc_comb)
    
    }
  
    
