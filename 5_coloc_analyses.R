 library(synapser) 
 library(coloc) 
 library(R.utils)
 library(TwoSampleMR)
 library(tidyr)
 library(dplyr)
 library(data.table)
 library(ieugwasr)
 library(genetics.binaRies)
 library(MendelianRandomization)
 
 synLogin()  # This is to log in to the Synapse platform, you need your own login data

 listformr <- fread(".../ukb_proteomics_cvd/output_files/3_MR_signfindings.csv")
  listformr <- listformr[,-"V1"]
  sens <- fread(".../ukb_proteomics_cvd/input_files/proteins_prs/sensitivityanalyses_posttwosamplesens_postallele_postallelemv.csv")
  listformr <- merge(listformr, sens[,c("exp", "outc", "all_sens")], by=c("exp", "outc"), all.x=T, all.y=F)
 sumstats_info <- fread(".../ukb_proteomics_cvd/input_files/2_olink_protein_map_mr.txt")                   
 sumstats_info <- sumstats_info[sumstats_info$Assay %in% listformr$exp, ]
 
 cad <- fread("https://storage.googleapis.com/finngen-public-data-r9/summary_stats/finngen_R9_I9_CHD.gz")                     
  cad$"#chrom" <- ifelse(cad$"#chrom"==23, "X", cad$"#chrom")                                                                 
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
 
 ref_rsid <- fread(".../tools/hg38_common_chrpos_X.txt")                                                  
 df_sum <- data.frame(exp=NA, outc=NA, nsnp=NA, method=NA, b=NA, se=NA, pval=NA)[-1,]
 df_instr <- data.frame(pos.exposure=NA, pos_id=NA, effect_allele.exposure=NA, other_allele.exposure=NA, effect_allele.outcome=NA, 
                       other_allele.outcome=NA, beta.exposure=NA, beta.outcome=NA, eaf.exposure=NA, eaf.outcome=NA, remove=NA, 
                       palindromic=NA, ambiguous=NA, id.outcome=NA, chr.outcome=NA, pos.outcome=NA, pval.outcome=NA, se.outcome=NA,
                       outcome=NA, mr_keep.outcome=NA, pval_origin.outcome=NA, chr.exposure=NA, samplesize.exposure=NA, se.exposure=NA,
                       pval.exposure=NA, exposure=NA, pval=NA, mr_keep.exposure=NA, pval_origin.exposure=NA, id.exposure=NA,
                       action=NA, mr_keep=NA, samplesize.outcome=NA, SNP=NA)[-1,]

 
   
 for (i in sumstats_info$Code[98:nrow(sumstats_info)]){

   syn_code <- synGet(entity=i, downloadLocation = paste(getwd(), "sumstat_prot", sep="/")) 
   untar(paste(syn_code$path), list=F, exdir=paste(syn_code$cacheDir))
   chrom <- fread(paste0(syn_code$cacheDir, "/", gsub(".tar", "", sumstats_info[sumstats_info$Code==i,]$Docname[1]), "/", 
                  "discovery_chr", sumstats_info[sumstats_info$Code==i,]$chr[1], "_", sumstats_info[sumstats_info$Code==i,]$UKBPPP_ProteinID[1], 
                  ":", sumstats_info[sumstats_info$Code==i,]$Panel[1], ".gz"))
   
   chrom <- chrom[chrom$GENPOS > (sumstats_info[sumstats_info$Code==i,]$gene_start[1] - 200000) &                        
                    chrom$GENPOS < (sumstats_info[sumstats_info$Code==i,]$gene_end[1] + 200000), ]
   chrom$P <- 10^-chrom$LOG10P
   chrom$CHROM <- ifelse(chrom$CHROM==23, "X", chrom$CHROM)                                                              
   
   
  for (outcome_get in c("cad", "hf", "af", "as")){
    
    outcome <- get(outcome_get)
    outcome_overlap <- outcome[outcome$`#chrom` == sumstats_info[sumstats_info$Code==i,]$chr[1],]                               
    outcome_overlap <- outcome_overlap[outcome_overlap$pos %in% chrom$GENPOS,]                                                 
    outcome_rsid <- outcome_overlap[,c("#chrom", "pos", "rsids")]                                                              
    outcome_rsid <- outcome_rsid %>% mutate(rsids = strsplit(as.character(rsids), ",")) %>% unnest(rsids)
    outcome_overlap$phen <- paste(outcome_get)
    outcome_overlap$id <- paste(outcome_overlap$`#chrom`, outcome_overlap$pos, outcome_overlap$alt, outcome_overlap$ref, sep=":")
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
    chrom_overlap <- format_data(chrom_overlap, type="exposure", phenotype_col="phen", snp_col="ID", beta_col="BETA", se_col="SE", eaf_col="A1FREQ",
                              effect_allele_col="ALLELE1", other_allele_col="ALLELE0", pval_col="LOG10P", chr_col="CHROM", samplesize_col="N", pos_col="GENPOS", log_pval=T)
    rm(chrom_overlap_2)
    dat <- harmonise_data(exposure_dat=chrom_overlap, outcome_dat=outcome_overlap)                                             
    
     if (sumstats_info[sumstats_info$Code==i,]$chr[1]=="X"){                                                                 
       dat <- merge(dat, ref_rsid[,c("V1", "V2", "V3")], by.x="pos.exposure", by.y="V2", all.x=T) 
       colnames(dat)[colnames(dat) %in% c("V1", "V3")] <- c("#chrom", "rsids")
     } else {
      dat <- merge(dat, outcome_rsid, by.x="pos.exposure", by.y="pos", all.x=T)
     }
     
    colnames(dat)[colnames(dat) %in% c("SNP", "rsids")] <- c("pos_id", "SNP")                                                  
    dat <- dat[order(dat$pval.exposure),]                                                                                      
    dat <- dat[!duplicated(dat$SNP),]
    dat <- dat[!duplicated(dat$pos_id),]
    ld <- ld_matrix(dat$SNP, bfile=".../tools/g1000_eur", plink_bin=genetics.binaRies::get_plink_binary())
    rownames(ld) <- gsub("_.*", "", rownames(ld))
    colnames(ld) <- gsub("_.*", "", colnames(ld))
    dat <- dat[dat$SNP %in% rownames(ld),]
    dat <- dat[match(rownames(ld), dat$SNP),]
   
    dat_exp <- list(beta=dat$beta.exposure, varbeta=dat$se.exposure^2, snp=dat$SNP, position=dat$pos.exposure, type="quant", sdY=1, LD=ld, N=unique(dat$samplesize.exposure))
    dat_outc <- list(beta=dat$beta.outcome, varbeta=dat$se.outcome^2, snp=dat$SNP, position=dat$pos.outcome, type="cc", LD=ld, N=ifelse(outcome_get=="cad", 377277, ifelse(outcome_get=="hf", 377277, ifelse(outcome_get=="af", 237690, 377277))), s=ifelse(outcome_get=="cad", 0.115347609, ifelse(outcome_get=="hf", 0.072371229, ifelse(outcome_get=="af", 0.192544911, 0.024260689))))
    
    coloc_outc <- coloc.abf(dataset1=dat_exp, dataset2=dat_outc)

    results <- data.frame(exp=sumstats_info[sumstats_info$Code==i,]$Assay[1], outc=paste(outcome_get), nsnps=coloc_outc$summary["nsnps"],
                      pp_h0=coloc_outc$summary["PP.H0.abf"], pp_h1=coloc_outc$summary["PP.H1.abf"], pp_h2=coloc_outc$summary["PP.H2.abf"],
                      pp_h3=coloc_outc$summary["PP.H3.abf"], pp_h4=coloc_outc$summary["PP.H4.abf"])
    df_sum <- rbind(df_sum, results)
    rm(results, outcome, outcome_overlap, outcome_rsid, chrom_overlap, dat, dat_exp, dat_outc, coloc_outc)

    
  }
    write.csv(df_sum, ".../tmp/outcomes/coloc_loop.csv") 
    print(paste(sumstats_info[sumstats_info$Code==i,]$Assay[1], "done"))
    unlink(paste0(gsub("\\\\[^\\\\]*$", "", syn_code$cacheDir)), recursive=T)
 }
  
 
