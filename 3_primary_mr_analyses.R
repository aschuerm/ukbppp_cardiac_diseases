 library(synapser) 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)
 
 setwd("/broad/hptmp/aschuerm/synapser")
  synLogin(email='XXX', authToken="XXX") # This is to log in to the Synapse platform, needed to gain access to the protein sum stats

 sumstats_info <- fread(".../ukb_proteomics_cvd/input_files/2_olink_protein_map_mr.txt")
 primary_findings <- fread(".../ukb_proteomics_cvd/output_files/1_adj_no_prev_timevar_model_hr.csv")
  cad_listformr <- fread(".../ukb_proteomics_cvd/output_files/1_primary_listformr/cad_listformr.txt")
  hf_listformr <- fread(".../ukb_proteomics_cvd/output_files/1_primary_listformr/hf_listformr.txt")
  af_listformr <- fread(".../ukb_proteomics_cvd/output_files/1_primary_listformr/af_listformr.txt")
  as_listformr <- fread(".../ukb_proteomics_cvd/output_files/1_primary_listformr/as_listformr.txt")
 sumstats_info <- sumstats_info[sumstats_info$Assay %in% cad_listformr$Protein | sumstats_info$Assay %in% hf_listformr$Protein | sumstats_info$Assay %in% af_listformr$Protein | sumstats_info$Assay %in% as_listformr$Protein,]
  sumstats_info <- sumstats_info[!((sumstats_info$Assay=="IL6" & sumstats_info$Panel!="Cardiometabolic")|(sumstats_info$Assay=="TNF" & sumstats_info$Panel!="Oncology")),]
 
 cad <- fread("https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_CHD.gz") 
  cad$"#chrom" <- ifelse(cad$"#chrom"==23, "X", cad$"#chrom") 
 hf <- fread("https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_HEARTFAIL.gz")               
  hf$"#chrom" <- ifelse(hf$"#chrom"==23, "X", hf$"#chrom")                                                                   
 af <- fread("https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_AF.gz")               
  af$"#chrom" <- ifelse(af$"#chrom"==23, "X", af$"#chrom")                                                                   
 as <- fread("https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_CAVS_OPERATED.gz")              
  as$"#chrom" <- ifelse(as$"#chrom"==23, "X", as$"#chrom")                                                                  

 ref_rsid <- fread(".../tools/hg38_common_chrpos_X.txt")
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
   
 for (i in sumstats_info$Code){

   syn_code <- synGet(entity=i, downloadLocation = paste(getwd(), "sumstat_prot", sep="/"))
   untar(paste(syn_code$path), list=F, exdir=paste(syn_code$cacheDir))
   chrom_u <- fread(paste0(syn_code$cacheDir, "/", gsub(".tar", "", sumstats_info[sumstats_info$Code==i,]$Docname[1]), "/", 
                  "discovery_chr", sumstats_info[sumstats_info$Code==i,]$chr[1], "_", sumstats_info[sumstats_info$Code==i,]$UKBPPP_ProteinID[1], 
                  ":", sumstats_info[sumstats_info$Code==i,]$Panel[1], ".gz"))
   
   chrom_u <- chrom_u[chrom_u$GENPOS > (sumstats_info[sumstats_info$Code==i,]$gene_start[1] - 200000) &  
                    chrom_u$GENPOS < (sumstats_info[sumstats_info$Code==i,]$gene_end[1] + 200000), ]
   chrom_u$P <- 10^-chrom_u$LOG10P
   
  for (outcome_get in c("cad", "hf", "af", "as")){
    outcome <- get(outcome_get)
    if (!(sumstats_info[sumstats_info$Code==i,]$Assay[1] %in% get(paste0(outcome_get, "_listformr"))$Protein)) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1], " and ", outcome_get))
    } else {
   
   for (pval_thresh in c(5e-6)){
   chrom <- chrom_u[chrom_u$P<pval_thresh,]                                                                                     
   chrom$CHROM <- ifelse(chrom$CHROM==23, "X", chrom$CHROM)                                                                     
   
   if (is.null(chrom) || nrow(chrom) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_overlap <- outcome[outcome$`#chrom` == sumstats_info[sumstats_info$Code==i,]$chr[1],]                              
     outcome_overlap <- outcome_overlap[outcome_overlap$pos %in% chrom$GENPOS,]
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")] 
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste(outcome_get)
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     chrom_overlap <- chrom[chrom$GENPOS %in% outcome_overlap$pos.outcome,]
     chrom_overlap_2 <- chrom_overlap
     chrom_overlap_2$BETA <- chrom_overlap_2$BETA * -1
     chrom_overlap_2$A1FREQ  <- 1 - chrom_overlap_2$A1FREQ
     colnames(chrom_overlap_2)[colnames(chrom_overlap_2) %in% c("ALLELE0", "ALLELE1")] <- c("ALLELE1", "ALLELE0")
     chrom_overlap <- rbind(chrom_overlap, chrom_overlap_2)
     chrom_overlap$ID <- paste(chrom_overlap$CHROM, chrom_overlap$GENPOS, chrom_overlap$ALLELE1, chrom_overlap$ALLELE0, sep=":")
     chrom_overlap$phen <- sumstats_info[sumstats_info$Code==i,]$Assay[1]
     chrom_overlap <- as.data.frame(chrom_overlap)
     chrom_overlap <- format_data(chrom_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="BETA", se_col="SE", eaf_col="A1FREQ",
                                effect_allele_col="ALLELE1", other_allele_col="ALLELE0", pval_col="LOG10P", chr_col="CHROM", samplesize_col="N", pos_col="GENPOS", log_pval=T)
     rm(chrom_overlap_2, chrom)
     dat_u <- harmonise_data(exposure_dat=chrom_overlap, outcome_dat=outcome_overlap) 
     
     if (sumstats_info[sumstats_info$Code==i,]$chr[1]=="X"){  
       dat_u <- merge(dat_u, ref_rsid[,c("V1", "V2", "V3")], by.x="pos.exposure", by.y="V2", all.x=T) 
       colnames(dat_u)[colnames(dat_u) %in% c("V1", "V3")] <- c("#chrom", "rsids")
     } else {
      dat_u <- merge(dat_u, outcome_rsid, by.x="pos.exposure", by.y="pos", all.x=T)
     }
     
     colnames(dat_u)[colnames(dat_u) %in% c("SNP", "rsids")] <- c("pos_id", "SNP") 
     dat_u <- dat_u[order(dat_u$pval.exposure),]
     dat_u <- dat_u[!duplicated(dat_u$SNP),]
     dat_u <- dat_u[!duplicated(dat_u$pos_id),]
     
   for (rsq_thresh in c(0.1)){
     clump <- ld_clump(dplyr::tibble(rsid=dat_u$SNP, pval=dat_u$pval.exposure, id=dat_u$id.exposure),
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq_thresh,
                      bfile = ".../tools/g1000_eur")
     dat <- dat_u[dat_u$SNP %in% clump$rsid,]
     rm(chrom_overlap, outcome_overlap, outcome_rsid, clump)
     
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
    if (nrow(dat)==1) {                                                             
       dat$marker_ld <- paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_")
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outcome_get), pvalthreshold=pval_thresh, 
                              rsqthreshold=rsq_thresh, nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval)
      } else if (nrow(dat)==2) { 
       ld <- ld_matrix(dat$SNP, bfile=".../tools/g1000_eur", plink_bin=genetics.binaRies::get_plink_binary())
       rownames(ld) <- gsub("\\_.*", "", rownames(ld))
       dat_dupl <- dat[!(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        colnames(dat_dupl)[which(colnames(dat_dupl) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
        dat_dupl$beta.exposure <- dat_dupl$beta.exposure*-1
        dat_dupl$beta.outcome <- dat_dupl$beta.outcome*-1
        dat_dupl$eaf.exposure <- 1-dat_dupl$eaf.exposure
        dat_dupl$eaf.outcome <- 1-dat_dupl$eaf.outcome
        dat <- dat[(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        dat <- rbind(dat, dat_dupl)
        dat <- dat[ order(match(dat$SNP, rownames(ld))), ]
        dat$marker_ld <- paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_")
        dat2 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome,
                                                 correlation = ld)
       output_mr_ivw_corr <- MendelianRandomization::mr_ivw(dat2, correl = TRUE)
       results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outcome_get), pvalthreshold=pval_thresh, 
                              rsqthreshold=rsq_thresh, nsnp=output_mr_ivw_corr@SNPs, 
                             method="Inverse variance weighted (correlation inc)", b=output_mr_ivw_corr@Estimate, 
                             se=output_mr_ivw_corr@StdError, pval=output_mr_ivw_corr@Pvalue)
      } else {
       ld <- ld_matrix(dat$SNP, bfile=".../tools/g1000_eur", plink_bin=genetics.binaRies::get_plink_binary())
       rownames(ld) <- gsub("\\_.*", "", rownames(ld))
       dat_dupl <- dat[!(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        colnames(dat_dupl)[which(colnames(dat_dupl) %in% c("effect_allele.exposure", "other_allele.exposure", "effect_allele.outcome", "other_allele.outcome"))] <- c("other_allele.exposure", "effect_allele.exposure", "other_allele.outcome", "effect_allele.outcome")
        dat_dupl$beta.exposure <- dat_dupl$beta.exposure*-1
        dat_dupl$beta.outcome <- dat_dupl$beta.outcome*-1
        dat_dupl$eaf.exposure <- 1-dat_dupl$eaf.exposure
        dat_dupl$eaf.outcome <- 1-dat_dupl$eaf.outcome
        dat <- dat[(paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_") %in% colnames(ld)),]
        dat <- rbind(dat, dat_dupl)
        dat <- dat[ order(match(dat$SNP, rownames(ld))), ]
        dat$marker_ld <- paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_")
       dat2 <- MendelianRandomization::mr_input(bx = dat$beta.exposure, bxse = dat$se.exposure,
                                                 by = dat$beta.outcome, byse = dat$se.outcome,
                                                 correlation = ld)
       output_mr_ivw_corr <- MendelianRandomization::mr_ivw(dat2, correl = TRUE)
       output_mr_egger_corr <- MendelianRandomization::mr_egger(dat2, correl = TRUE)
       results1 <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outcome_get), pvalthreshold=pval_thresh, 
                              rsqthreshold=rsq_thresh, nsnp=output_mr_ivw_corr@SNPs, 
                             method="Inverse variance weighted (correlation inc)", b=output_mr_ivw_corr@Estimate, 
                             se=output_mr_ivw_corr@StdError, pval=output_mr_ivw_corr@Pvalue)
       results2 <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outcome_get), pvalthreshold=pval_thresh, 
                              rsqthreshold=rsq_thresh, nsnp=output_mr_egger_corr@SNPs, 
                             method="Egger (correlation inc)", b=output_mr_egger_corr@Estimate, 
                             se=output_mr_egger_corr@StdError.Est, pval=output_mr_egger_corr@Pvalue.Est)
       results3 <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outcome_get), pvalthreshold=pval_thresh, 
                              rsqthreshold=rsq_thresh, nsnp=output_mr_egger_corr@SNPs, 
                             method="Egger intercept (correlation inc)", b=output_mr_egger_corr@Intercept, 
                             se=output_mr_egger_corr@StdError.Int, pval=output_mr_egger_corr@Pvalue.Int)
       results <- rbind(results1, results2, results3)
       rm(results1, results2, results3)
       
      }
      }
     
   if (is.null(results) || nrow(results) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
      df_sum <- rbind(df_sum, results)
      df_instr <- rbind(df_instr, dat[,-c(34)])
      rm(dat, results, ld, dat2, output_mr_ivw_corr, output_mr_egger_corr)
   }
   }
   }
     }
   }
    }
  }
      
   write.csv(df_sum, ".../tmp/outcomes/primaryanalyses_loop_corr.csv") 
   write.csv(df_instr, ".../tmp/outcomes/primaryanalyses_loop_corr_instr.csv") 
   print(paste(sumstats_info[sumstats_info$Code==i,]$Assay[1], "done"))
   unlink(paste0(gsub("\\\\[^\\\\]*$", "", syn_code$cacheDir)), recursive=T)
}
