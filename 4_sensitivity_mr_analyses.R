 library(synapser) 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)
 
 setwd(".../tmp/synapser")
  synLogin(email='XXX', authToken="XXX") # This is to log in to the Synapse platform, needed to gain access to the protein sum stats

 listformr <- fread(".../ukb_proteomics_cvd/output_files/3_MR_signfindings.csv")
 sumstats_info <- fread(".../ukb_proteomics_cvd/input_files/2_olink_protein_map_mr.txt")                   # Information on the sum stats / protein codes etc.; functions as a linker file
 sumstats_info <- sumstats_info[sumstats_info$Assay %in% listformr$exp, ]
  sumstats_info <- sumstats_info[!((sumstats_info$Assay=="IL6" & sumstats_info$Panel!="Cardiometabolic")|(sumstats_info$Assay=="TNF" & sumstats_info$Panel!="Oncology")),]
 
 cad <- fread("https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_CHD.gz")                     # This is to download your outcome summary statistics of interest
  cad$"#chrom" <- ifelse(cad$"#chrom"==23, "X", cad$"#chrom")                                                                   # We also transfer the 23rd chromosome to "X" for consistency across sum stats
  cad_listformr <- listformr[listformr$outc=="CAD",]
 hf <- fread("https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_HEARTFAIL.gz")               
  hf$"#chrom" <- ifelse(hf$"#chrom"==23, "X", hf$"#chrom")                                                                   
  hf_listformr <- listformr[listformr$outc=="HF",]
 af <- fread("https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_AF.gz")               
  af$"#chrom" <- ifelse(af$"#chrom"==23, "X", af$"#chrom")                                                                   
  af_listformr <- listformr[listformr$outc=="AF",]
 as <- fread("https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_CAVS_OPERATED.gz")              
  as$"#chrom" <- ifelse(as$"#chrom"==23, "X", as$"#chrom")                                                                  
  as_listformr <- listformr[listformr$outc=="AS",]
 
 ref_rsid <- fread(".../tools/hg38_common_chrpos_X.txt")                                                   # This is a linker file to match RSIDs with chromosome positions for the X chromosome
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
   
 for (i in sumstats_info$Code){

   syn_code <- synGet(entity=i, downloadLocation = paste(getwd(), "sumstat_prot", sep="/")) # Downloading the summary statistics for the protein of interest
   untar(paste(syn_code$path), list=F, exdir=paste(syn_code$cacheDir))
   chrom_u <- fread(paste0(syn_code$cacheDir, "/", gsub(".tar", "", sumstats_info[sumstats_info$Code==i,]$Docname[1]), "/", 
                  "discovery_chr", sumstats_info[sumstats_info$Code==i,]$chr[1], "_", sumstats_info[sumstats_info$Code==i,]$UKBPPP_ProteinID[1], 
                  ":", sumstats_info[sumstats_info$Code==i,]$Panel[1], ".gz"))
   
   chrom_u <- chrom_u[chrom_u$GENPOS > (sumstats_info[sumstats_info$Code==i,]$gene_start[1] - 200000) &                              # Selecting the cis region only (here defined as 1Mb before or after the protein-encoding region)
   chrom_u$GENPOS < (sumstats_info[sumstats_info$Code==i,]$gene_end[1] + 200000), ]
   chrom_u$P <- 10^-chrom_u$LOG10P
   
  for (outcome_get in c("cad", "hf", "af", "as")){
    outcome <- get(outcome_get)
    if (!(sumstats_info[sumstats_info$Code==i,]$Assay[1] %in% get(paste0(outcome_get, "_listformr"))$exp)) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1], " and ", outcome_get))
    } else {
   
   for (pval_thresh in c(5e-4, 5e-6, 5e-8)){
   chrom <- chrom_u[chrom_u$P<pval_thresh,]                                                                                                # Selecting "region-wide" significant cis-pQTs
   chrom$CHROM <- ifelse(chrom$CHROM==23, "X", chrom$CHROM)                                                                     # Transferring 23rd chromosome to "X" for consistency
   
   if (is.null(chrom) || nrow(chrom) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_overlap <- outcome[outcome$`#chrom` == sumstats_info[sumstats_info$Code==i,]$chr[1],]                                      # Only selecting the chromosome of interest to speed stuff up downstream from here
     outcome_overlap <- outcome_overlap[outcome_overlap$pos %in% chrom$GENPOS,]                                                 # Only selecting the variants that are overlapping
     
   if (is.null(outcome_overlap) || nrow(outcome_overlap) == 0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
   } else {
     outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              # These next couple of lines of code just wrangle the data so that the "TwoSampleMR" package can read everything and do its magic
     outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
     outcome_overlap$phen <- paste(outcome_get)
     outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
     outcome_overlap <- as.data.frame(outcome_overlap)
     outcome_overlap <- format_data(outcome_overlap, type="outcome", phenotype_col="phen", snp_col="id", beta_col="beta", se_col="sebeta", eaf_col="af_alt",
                                effect_allele_col="alt", other_allele_col="ref", pval_col="pval", chr_col="#chrom", pos_col="pos")
     
     chrom_overlap <- chrom[chrom$GENPOS %in% outcome_overlap$pos.outcome,]                                                     # Again, we just take the overlapping variants (now in the other direction)
     chrom_overlap_2 <- chrom_overlap                                                                                           # Because the order of effect allele and other allele is random, we make a second dataframe with the opposite order of these alleles to optimize matching between sum stats
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
     dat_u <- harmonise_data(exposure_dat=chrom_overlap, outcome_dat=outcome_overlap)                                             # This is where the matching happens
     
     if (sumstats_info[sumstats_info$Code==i,]$chr[1]=="X"){                                                                    # This little if-else-statement just makes sure that you get the appropriate RSIDs for each variant; because the X-chromosome requires an additional file, this one is in a separate loop
       dat_u <- merge(dat_u, ref_rsid[,c("V1", "V2", "V3")], by.x="pos.exposure", by.y="V2", all.x=T) 
       colnames(dat_u)[colnames(dat_u) %in% c("V1", "V3")] <- c("#chrom", "rsids")
     } else {
      dat_u <- merge(dat_u, outcome_rsid, by.x="pos.exposure", by.y="pos", all.x=T)
     }
     
     colnames(dat_u)[colnames(dat_u) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  # We make sure that our reference column (which should have the name "SNP") is the RSID column
     dat_u <- dat_u[order(dat_u$pval.exposure),]                                                                                      # We make sure there are no duplicate SNPs (e.g., SNPs with the same position but other alleles [this messes the MR itself up])
     dat_u <- dat_u[!duplicated(dat_u$SNP),]
     
   for (rsq_thresh in c(0.001, 0.01, 0.1, 0.2, 0.4)){
     clump <- ld_clump(dplyr::tibble(rsid=dat_u$SNP, pval=dat_u$pval.exposure, id=dat_u$id.exposure),                                 # Clumping (i.e., excluding the variants that are correlated with each other); you'll need the 1000G LD reference file for this
                      plink_bin = genetics.binaRies::get_plink_binary(), clump_kb=10000, clump_r2 = rsq_thresh,
                      bfile = ".../tools/g1000_eur")
     dat <- dat_u[dat_u$SNP %in% clump$rsid,]
     rm(chrom_overlap, outcome_overlap, outcome_rsid, clump)
     
     
     # Note: this particular script uses the IVW method adjusted for between-variant correlation. This is not standard, but is a good method to use when using a lenient R2 threshold such as the one we use (0.1) when using proteins as the exposure.
     
     if (nrow(dat[dat$mr_keep,])==0) {
        print(paste0("Skipping ", sumstats_info[sumstats_info$Code==i,]$Assay[1]))
       results <- NULL
      } else {
    
    if (nrow(dat)==1) {                                                                                                        # This is where the magic happens: if you have 1 variant, you use the Wald ratio as your method
       dat$marker_ld <- paste(dat$SNP, dat$effect_allele.exposure, dat$other_allele.exposure, sep="_")
       results_mr <- mr(dat, method_list=c("mr_wald_ratio"))
       results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outcome_get), pvalthreshold=pval_thresh, 
                              rsqthreshold=rsq_thresh, nsnp=results_mr$nsnp, method=results_mr$method, b=results_mr$b, 
                           se=results_mr$se, pval=results_mr$pval)
      } else if (nrow(dat)==2) {                                                                                                # If you have 2 variants, you can use the classic IVW method but not the MR-Egger method
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
      } else {                                                                                                                  # If you have more than 2 variants, you can do anything (including IVW and MR-Egger)
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
      output_mr_ivw_corr <- try(MendelianRandomization::mr_ivw(dat2, correl = TRUE))
      
      if (inherits(output_mr_ivw_corr, "try-error")) {results <- NULL} else {
      
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
      
   write.csv(df_sum, ".../tmp/outcomes/sensitivityanalyses_loop_corr.csv") 
   write.csv(df_instr, ".../tmp/outcomes/sensitivityanalyses_loop_corr_instr.csv") 
   print(paste(sumstats_info[sumstats_info$Code==i,]$Assay[1], "done"))
   unlink(paste0(gsub("\\\\[^\\\\]*$", "", syn_code$cacheDir)), recursive=T)
}
 
 
