### 0 - Libraries and dataframes ####  

  library(ggVennDiagram)
  library(data.table)
  library(sf)
  library(ggplot2)
  library(compareGroups)
  library(pROC)
  library(forcats)
  library(ggpubr)
  library(survival)
  library(survminer)
  library(ggtext)
  library(tidyverse)
  library(dplyr)
 

  ### 0 - Proteins ####

  proteins <- fread(".../ukb_proteomics_cvd/input_files/tested_proteins_suppltable.csv")
    hpa <- fread(".../ukb_proteomics_cvd/input_files/proteinatlas.tsv")
    hpa <- hpa[hpa$`Secretome location`=="Secreted to blood",]
    hpa <- hpa[,c("Gene", "Secretome location")]
  proteins$olink <- proteins$"Protein (abbreviation)"
    proteins <- separate(proteins, olink, into = c("before", "after"), sep = "_")
    proteins$after <- ifelse(is.na(proteins$after), "Nothing", proteins$after)
    proteins$after <- ifelse(proteins$before=="NT-proBNP", "NPPB", proteins$after)
    proteins$blood_secr <- ifelse(proteins$before %in% hpa$Gene, "Yes", ifelse(proteins$after %in% hpa$Gene, "Yes", "No"))
    # write.csv(proteins, ".../ukb_proteomics_cvd/output_files/0_tested_proteins_hpa.csv", row.names=F)

  df <- fread('.../ukb_proteomics_cvd/output_files/1_adj_no_prev_timevar_model_hr.csv')[,-1]
    df_unique <- df[df$P_Value<(0.05/(1459*4)),]
    df_unique <- data.frame(Protein=unique(df_unique$Protein))
    df_unique$Protein <- ifelse(df_unique$Protein=="NTproBNP", "NT-proBNP", df_unique$Protein)
    df_unique <- merge(df_unique, proteins[,c("Protein (abbreviation)", "blood_secr")], by.x="Protein", by.y="Protein (abbreviation)", all.x=F, all.y=F)
    df_unique <- df_unique[!duplicated(df_unique$Protein),]
    df_unique$group <- "sign"
    
    df_unique_non <- data.frame(Protein=unique(df$Protein[!(df$Protein%in%df_unique$Protein)]))
    df_unique_non <- merge(df_unique_non, proteins[,c("Protein (abbreviation)", "blood_secr")], by.x="Protein", by.y="Protein (abbreviation)", all.x=F, all.y=F)
    df_unique_non <- df_unique_non[!duplicated(df_unique_non$Protein),]
    df_unique_non$group <- "nonsign"
    
    df_together <- rbind(df_unique, df_unique_non)
      summary(as.factor(df_unique$blood_secr))
        sum(df_unique$blood_secr=="Yes")/nrow(df_unique)
      summary(as.factor(df_unique_non$blood_secr))
        sum(df_unique_non$blood_secr=="Yes")/nrow(df_unique_non)
    chisq.test(df_together$blood_secr, df_together$group)
    

  ### 1 - Baseline table ####  
    # 1.A - Baseline ####  

  df <- fread(".../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz")
  df$SmokingStatusv2_nona <- ifelse(is.na(df$SmokingStatusv2), "Never", df$SmokingStatusv2)
  
  table_df <- compareGroups(y = df$Sex, data = df[,c(1:105, 1598)], method = 1, alpha = 0.05, min.dis = 5, max.ylev = 5, include.label = TRUE, Q1 = 0.25, Q3 = 0.75, simplify = TRUE, ref = 1, fact.ratio = 1, ref.y = 1)
    pvals_df  <- getResults(table_df , "p.overall")
    export_table_df <- createTable(table_df, hide.no = 0, type = 2, show.ratio = FALSE, show.n = FALSE, show.descr = TRUE, show.all = TRUE, show.p.overall = TRUE, show.p.trend = FALSE, show.p.mul = TRUE)
    export2xls(export_table_df, file = ".../ukb_proteomics_cvd/output_files/0_baseline_table.xlsx")
  table_df <- compareGroups(y = df$Sex, data = df[,c(1:105, 1598)], method = NA, alpha = 0.05, min.dis = 5, max.ylev = 5, include.label = TRUE, Q1 = 0.25, Q3 = 0.75, simplify = TRUE, ref = 1, fact.ratio = 1, ref.y = 1)
    pvals_df  <- getResults(table_df , "p.overall")
    export_table_df <- createTable(table_df, hide.no = 0, type = 2, show.ratio = FALSE, show.n = FALSE, show.descr = TRUE, show.all = TRUE, show.p.overall = TRUE, show.p.trend = FALSE, show.p.mul = TRUE)
    export2xls(export_table_df, file = ".../ukb_proteomics_cvd/output_files/0_baseline_table_nonnormal.xlsx")
    rm(table_df, pvals_df, export_table_df)
  
    
    # 1.B - Cumulative incidence ####  

  library(data.table)
  library(boot)
  library(broom)
  library(dplyr)
  library(survival)
  library(survminer)

  a <- fread(".../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd_incdis.tsv.gz")

  surv_object <- Surv(time = a$fu, event = a$inc)
  fit1 <- survfit(surv_object ~ dis, data=a)
  pl <- ggsurvplot(fit1, pval=F, fun='event', color = "dis", ylim=c(0,0.1), legend=c(0.19,0.76), palette=c("#3C5488", "#397C7D", "#36A471", "#33CC66"), 
                     censor=F, risk.table=T, ylab="Cumulative incidence", xlab="Follow-up time (years)", legend.title = "Outcome:",
                     ggtheme = theme_classic(), risk.table.fontsize=3, risk.table.y.text=F)
  
  tiff(".../ukb_proteomics_cvd/output_files/1_inc_outcomes.tiff", width=2000, height=1800, res=300)
  pl
  dev.off()
  
  df_source <- data.frame(time=fit1$time, n.risk=fit1$cumhaz, fit1$n.censor, 
                          strata=c(rep("Aortic stenosis", fit1$strata[1]), rep("Atrial fibrillation", fit1$strata[2]), rep("Coronary artery disease", fit1$strata[3]), rep("Heart failure", fit1$strata[4])))
  write.csv(df_source, ".../ukb_proteomics_cvd/output_files/source_data/ext_data_fig_1.csv", row.names=F)
  
  
  
    # 1.C - Correlation matrix ####

  df <- fread(".../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz")
    index <- which(colnames(df)=="CLIP2"):which(colnames(df)=="SCARB2")
    df_prot <- df[,..index]
    cor_prot <- cor(df_prot, method="pearson")
    # write.csv(cor_prot, ".../ukb_proteomics_cvd/output_files/0_correlation_proteins.csv", row.names=F)
    
    
    # 1.D - Correlation matrix visual ####
    
  library(pheatmap)
    
  corr_prot <- fread(".../ukb_proteomics_cvd/output_files/0_correlation_proteins.csv")
    heatmap <- pheatmap(corr_prot, show_colnames=F, show_rownames=F,)
    
    tiff(".../ukb_proteomics_cvd/output_files/0_corr_matrix.tiff", width=2500, height=2500, res=300)
    heatmap
    dev.off()

    
  ### 2 - Venn ####  
  
  timevar_excl <- fread('.../ukb_proteomics_cvd/output_files/1_adj_no_prev_timevar_model_hr.csv')

  timevar_excl_list <- list(AF=timevar_excl$Protein[timevar_excl$Outcome=="Afib" & timevar_excl$P_Value<(0.05/(1459*4))], 
                    AS=timevar_excl$Protein[timevar_excl$Outcome=="AS" & timevar_excl$P_Value<(0.05/(1459*4))], 
                    CAD=timevar_excl$Protein[timevar_excl$Outcome=="CAD" & timevar_excl$P_Value<(0.05/(1459*4))], 
                    HF=timevar_excl$Protein[timevar_excl$Outcome=="HF" & timevar_excl$P_Value<(0.05/(1459*4))])
  venn <- Venn(timevar_excl_list)
    d <- process_data(venn)
    d2 <- process_data(venn)
    d2@region <- st_polygonize(d@setEdge)
  tiff('C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/0_venn_assoc_timevar_excl.tiff', width=3000, height=2100, res=300)
  ggplot() +
    geom_sf(aes(fill = name), data = venn_region(d2)) +
    geom_sf(aes(color = name), data = venn_setedge(d)) +
    geom_sf_text(aes(label = name), data = venn_setlabel(d)) +
    geom_sf_text(aes(label = count), data = venn_region(d)) +
    scale_color_manual(values = alpha(c("#DC0000B2", "navy", "darkgreen", "grey"), .1), name="Outcome:") +
    scale_fill_manual(values = alpha(c("#DC0000B2", "navy", "darkgreen", "grey"), .1), name="Outcome:") +
    theme_void() 
  dev.off()
  
  timevar_excl[order(timevar_excl$P_Value),][timevar_excl[order(timevar_excl$P_Value),]$Outcome=="CAD",][1:15, ]
    timevar_excl[order(timevar_excl$P_Value),][timevar_excl[order(timevar_excl$P_Value),]$Outcome=="HF",][1:15, ]
    timevar_excl[order(timevar_excl$P_Value),][timevar_excl[order(timevar_excl$P_Value),]$Outcome=="Afib",][1:15, ]
    timevar_excl[order(timevar_excl$P_Value),][timevar_excl[order(timevar_excl$P_Value),]$Outcome=="AS",][1:15, ]

  write.csv(as.data.frame(d2$regionData)[,c(2,4)], ".../ukb_proteomics_cvd/output_files/source_data/ext_data_fig_3.csv", row.names=F)
  

  ### 3 - Manhattan ####

    # 3.A - Overall ####
    
  library(data.table)
  library(tidyverse)
  library(ggtext)
  library(gtools)
  library(ggrepel)
  library(ggbreak)
    
  df <- fread('.../ukb_proteomics_cvd/output_files/1_adj_no_prev_timevar_model_hr.csv')[,-1]
  link <- fread(".../ukb_proteomics_cvd/input_files/2_olink_protein_map_mr.txt", select=c("Assay", "chr", "gene_start"))
  link <- distinct(link, Assay, .keep_all = T)
  df <- merge(df, link, by.x="Protein", by.y="Assay", all.x=T)
  df$Protein <- sub("_", " / ", df$Protein)
  df$Protein <- sub("NTproBNP", "NT-proBNP", df$Protein)
  df$chr <- factor(df$chr, levels=c(1:22, "X"))
  df$p_min <- ifelse(df$HR>1, -log10(df$P_Value), log10(df$P_Value))
  
  data_cum <- df %>% 
                group_by(chr) %>% 
                summarise(max_bp = max(gene_start)) %>%
                mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) %>% 
                select(chr, bp_add)

  df <- df %>% 
                inner_join(data_cum, by = "chr") %>% 
                mutate(bp_cum = gene_start + bp_add) 
  axis_set <- df %>% 
                group_by(chr) %>% 
                summarize(center = mean(bp_cum))
  ylim <- df %>% 
                filter(P_Value == min(P_Value)) %>% 
                mutate(ylim = abs(floor(log10(P_Value))) + 2) %>% 
                pull(ylim)
  sig <- 0.5/5836

  
  
    # 3.B - CAD ####

  df_sub <- df[df$Outcome=="CAD",]
  df_sub$col <- factor(ifelse(df_sub$P_Value<sig & df_sub$Protein%in%df[df$Outcome %in% c("HF", "Afib", "AS")&df$P_Value<sig,]$Protein,
                        "A", ifelse(df_sub$P_Value<sig, "B", df_sub$chr)), levels=c(1:23, "A", "B"))
  df_sub$Protein_lab <- ifelse(df_sub$P_Value<5e-15, df_sub$Protein, "")
  df_sub_cad <- df_sub
  plot_cad <- ggplot(df_sub, aes(x = bp_cum, y = p_min, 
                                    color = col)) +
    geom_hline(yintercept = -log10(sig), color = "#DC0000B2", linetype = "solid") + 
    geom_hline(yintercept = 0, color = "black", linetype = "solid") + 
    geom_hline(yintercept = log10(sig), color = "#DC0000B2", linetype = "solid") + 
    geom_point(alpha = 1, size=1.5) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(-15, 65)) +
     scale_color_manual(values = c(rep(c("grey90", "grey80"), 
                                    11), "grey90", "#3C5488", "#6BB56B")) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = "-log<sub>10</sub>(P)") + 
    theme_classic() +
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    ) +
    geom_text_repel(aes(label=Protein_lab), size=2)

  
    # 3.C - HF ####

  df_sub <- df[df$Outcome=="HF",]
  df_sub$col <- factor(ifelse(df_sub$P_Value<sig & df_sub$Protein%in%df[df$Outcome %in% c("CAD", "Afib", "AS")&df$P_Value<sig,]$Protein,
                        "A", ifelse(df_sub$P_Value<sig, "B", df_sub$chr)), levels=c(1:23, "A", "B"))
  df_sub$Protein_lab <- ifelse(df_sub$P_Value<5e-25, df_sub$Protein, "")
  df_sub_hf <- df_sub
  plot_hf <- ggplot(df_sub, aes(x = bp_cum, y = p_min, 
                                    color = col)) +
    geom_hline(yintercept = -log10(sig), color = "#DC0000B2", linetype = "solid") + 
    geom_hline(yintercept = 0, color = "black", linetype = "solid") + 
    geom_hline(yintercept = log10(sig), color = "#DC0000B2", linetype = "solid") + 
    geom_point(alpha = 1, size=1.5) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(-15, 65)) +
     scale_color_manual(values = c(rep(c("grey90", "grey80"), 
                                    11), "grey90", "#3C5488", "#6BB56B")) +
   scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = "-log<sub>10</sub>(P)") + 
    theme_classic() +
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    ) +
    geom_text_repel(aes(label=Protein_lab), size=2)


    # 3.D - AF ####

  df_sub <- df[df$Outcome=="Afib",]
  df_sub$col <- factor(ifelse(df_sub$P_Value<sig & df_sub$Protein%in%df[df$Outcome %in% c("HF", "CAD", "AS")&df$P_Value<sig,]$Protein,
                        "A", ifelse(df_sub$P_Value<sig, "B", df_sub$chr)), levels=c(1:23, "A", "B"))
  df_sub$Protein_lab <- ifelse(df_sub$P_Value<5e-15, df_sub$Protein, "")
  df_sub_af <- df_sub
  plot_af <- ggplot(df_sub, aes(x = bp_cum, y = p_min, 
                                    color = col)) +
    geom_hline(yintercept = -log10(sig), color = "#DC0000B2", linetype = "solid") + 
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +
    geom_hline(yintercept = log10(sig), color = "#DC0000B2", linetype = "solid") + 
    geom_point(alpha = 1, size=1.5) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(-15, 65)) +
     scale_color_manual(values = c(rep(c("grey90", "grey80"), 
                                    11), "grey90", "#3C5488", "#6BB56B")) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = "-log<sub>10</sub>(P)") + 
    theme_classic() +
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    ) +
    geom_text_repel(aes(label=Protein_lab), size=2)


  plot_af_add1 <- ggplot(df_sub, aes(x = bp_cum, y = p_min, 
                                    color = col)) +
    geom_hline(yintercept = -log10(sig), color = "#DC0000B2", linetype = "solid") + 
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +
    geom_hline(yintercept = log10(sig), color = "#DC0000B2", linetype = "solid") + 
    geom_point(alpha = 1, size=1.5) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(170,180), breaks=c(170,180)) +
     scale_color_manual(values = c(rep(c("grey90", "grey80"), 
                                    11), "grey90", "#3C5488", "#6BB56B")) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = "-log<sub>10</sub>(P)") + 
    theme_classic() +
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    ) +
    geom_text_repel(aes(label=Protein_lab), size=2)


  plot_af_add2 <- ggplot(df_sub, aes(x = bp_cum, y = p_min, 
                                    color = col)) +
    geom_hline(yintercept = -log10(sig), color = "#DC0000B2", linetype = "solid") + 
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +
    geom_hline(yintercept = log10(sig), color = "#DC0000B2", linetype = "solid") + 
    geom_point(alpha = 1, size=1.5) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(90,100), breaks=c(90,100)) +
     scale_color_manual(values = c(rep(c("grey90", "grey80"), 
                                    11), "grey90", "#3C5488", "#6BB56B")) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = "-log<sub>10</sub>(P)") + 
    theme_classic() +
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    ) +
    geom_text_repel(aes(label=Protein_lab), size=2)


    # 3.E - AS ####

  df_sub <- df[df$Outcome=="AS",]
  df_sub$col <- factor(ifelse(df_sub$P_Value<sig & df_sub$Protein%in%df[df$Outcome %in% c("HF", "Afib", "CAD")&df$P_Value<sig,]$Protein,
                        "A", ifelse(df_sub$P_Value<sig, "B", df_sub$chr)), levels=c(1:23, "A", "B"))
  df_sub$Protein_lab <- ifelse(df_sub$P_Value<1e-5, df_sub$Protein, "")
  df_sub_as <- df_sub
  plot_as <- ggplot(df_sub, aes(x = bp_cum, y = p_min, 
                                    color = col)) +
    geom_hline(yintercept = -log10(sig), color = "#DC0000", linetype = "solid") + 
    geom_hline(yintercept = 0, color = "black", linetype = "solid") +
    geom_hline(yintercept = log10(sig), color = "#DC0000", linetype = "solid") + 
    geom_point(alpha = 1, size=1.5) +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(-15, 65)) +
     scale_color_manual(values = c(rep(c("grey90", "grey80"), 
                                    11), "grey90", "#3C5488", "#6BB56B")) +
    scale_size_continuous(range = c(0.5,3)) +
    labs(x = NULL, 
         y = "-log<sub>10</sub>(P)") + 
    theme_classic() +
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    ) +
    geom_text_repel(aes(label=Protein_lab), size=2)


    # 3.F - PPT ####
  
  df_sub <- rbind(df_sub_cad, df_sub_af, df_sub_hf, df_sub_as)
  write.csv(df_sub[,c("Protein", "Outcome", "p_min", "bp_cum", "col")], ".../ukb_proteomics_cvd/output_files/source_data/fig_2.csv", row.names=F)
  
  library(rvg)
  library(officer)
  
    Comp_1_dml <- dml(code = print(plot_cad, newpage = FALSE))
    Comp_2_dml <- dml(code = print(plot_hf, newpage = FALSE))
    Comp_3_dml <- dml(code = print(plot_af, newpage = FALSE))
    Comp_4_dml <- dml(code = print(plot_af_add1, newpage = FALSE))
    Comp_5_dml <- dml(code = print(plot_af_add2, newpage = FALSE))
    Comp_6_dml <- dml(code = print(plot_as, newpage = FALSE))
    
    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_1_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")
    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_2_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")
    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_3_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")
    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_4_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")
    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_5_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")
    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_6_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")

    rm(plot_af, plot_af_add1, plot_af_add2, plot_as, plot_cad, plot_hf, 
       Comp_1_dml, Comp_2_dml, Comp_3_dml, Comp_4_dml, Comp_5_dml, Comp_6_dml, df, df_sub, link, axis_set, sig, ylim, data_cum)
    
    
    # 3.G - Violin plots top hits ####
    
    df <- fread('.../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz')
    
    plot_cad <- ggplot(df, aes(x=as.factor(cad_inc), y=GDF15, fill=as.factor(cad_inc))) +
      geom_violin(trim=F) +
      geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
      scale_fill_manual(values=c("#9EAAC4","#3C5488"), labels=c("Controls", "Cases"), name="") +
      scale_x_discrete(labels=c("Controls", "Cases")) +
      scale_y_continuous(limits=c(-5,8)) +
      theme_classic() +
      labs(title="Coronary artery disease", x="Group", y = "Circulating GDF15 (SD)")
    plot_hfail <- ggplot(df, aes(x=as.factor(hfail_inc), y=WFDC2, fill=as.factor(hfail_inc))) +
      geom_violin(trim=F) +
      geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
      scale_fill_manual(values=c("#9EAAC4","#3C5488"), labels=c("Controls", "Cases"), name="") +
      scale_x_discrete(labels=c("Controls", "Cases")) +
      scale_y_continuous(limits=c(-5,8)) +
      theme_classic() +
      labs(title="Heart failure", x="Group", y = "Circulating WFDC2 (SD)")
    plot_afib <- ggplot(df, aes(x=as.factor(afib_inc), y=NTproBNP, fill=as.factor(afib_inc))) +
      geom_violin(trim=F) +
      geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
      scale_fill_manual(values=c("#9EAAC4","#3C5488"), labels=c("Controls", "Cases"), name="") +
      scale_x_discrete(labels=c("Controls", "Cases")) +
      scale_y_continuous(limits=c(-5,8)) +
      theme_classic() +
      labs(title="Atrial fibrillation", x="Group", y = "Circulating NT-proBNP (SD)")
    plot_ao_sten <- ggplot(df, aes(x=as.factor(ao_sten_inc), y=GDF15, fill=as.factor(ao_sten_inc))) +
      geom_violin(trim=F) +
      geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
      scale_fill_manual(values=c("#9EAAC4","#3C5488"), labels=c("Controls", "Cases"), name="") +
      scale_x_discrete(labels=c("Controls", "Cases")) +
      scale_y_continuous(limits=c(-5,8)) +
      theme_classic() +
      labs(title="Aortic stenosis", x="Group", y = "Circulating GDF15 (SD)")

      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_violins.tiff", width=5300, height=1500, res=300)
        ggarrange(plot_cad, plot_hfail, plot_afib, plot_ao_sten, common.legend=T, ncol=4, nrow=1, legend="bottom")
        dev.off()

    
    # 3.H - Distributions top hits ####
    
    df <- fread('.../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz')
    
    plot_cad <- ggplot(df, aes(x=GDF15, fill=as.factor(cad_inc))) +
      geom_density(alpha=0.4) +
      scale_fill_manual(values=c("#9EAAC4","#3C5488"), labels=c("Controls", "Cases"), name="") +
      theme_classic() +
      labs(title="Coronary artery disease", x="Circulating GDF15 (SD)", y = "Density")
    plot_hfail <- ggplot(df, aes(x=WFDC2, fill=as.factor(hfail_inc))) +
      geom_density(alpha=0.4) +
      scale_fill_manual(values=c("#9EAAC4","#3C5488"), labels=c("Controls", "Cases"), name="") +
      theme_classic() +
      labs(title="Heart failure", x = "Circulating WFDC2 (SD)", y="Density")
    plot_afib <- ggplot(df, aes(x=NTproBNP, fill=as.factor(afib_inc))) +
      geom_density(alpha=0.4) +
      scale_fill_manual(values=c("#9EAAC4","#3C5488"), labels=c("Controls", "Cases"), name="") +
      theme_classic() +
      labs(title="Atrial fibrillation", x = "Circulating NT-proBNP (SD)", y="Density")
    plot_ao_sten <- ggplot(df, aes(x=GDF15, fill=as.factor(ao_sten_inc))) +
      geom_density(alpha=0.4) +
      scale_fill_manual(values=c("#9EAAC4","#3C5488"), labels=c("Controls", "Cases"), name="") +
      theme_classic() +
      labs(title="Aortic stenosis", x = "Circulating GDF15 (SD)", y="Density")

      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_distributions.tiff", width=2*3000, height=2*800, res=2*300)
        ggarrange(plot_cad, plot_hfail, plot_afib, plot_ao_sten, common.legend=T, ncol=4, nrow=1, legend="bottom")
        dev.off()
      ggsave(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_distributions.svg", ggarrange(plot_cad, plot_hfail, plot_afib, plot_ao_sten, common.legend=T, ncol=4, nrow=1, legend="bottom"), width=12, height=3.2)
        
    plot_hfail_1 <- ggplot(df, aes(x=TNFRSF10B, fill=as.factor(hfail_inc))) +
      geom_density(alpha=0.4) +
      scale_x_continuous(limits=c(-5,5)) +
      geom_vline(xintercept = quantile(df$TNFRSF10B[df$hfail_inc==0], 0.95)) +
      scale_fill_manual(values=c("#9EAAC4","#3C5488"), labels=c("Controls", "Cases"), name="") +
      theme_classic() +
      labs(title="Heart failure", x = "Circulating TNFRSF10B (SD)", y="Density")
    plot_hfail_2 <- ggplot(df, aes(x=IGFBP4, fill=as.factor(hfail_inc))) +
      geom_density(alpha=0.4) +
      scale_x_continuous(limits=c(-5,5)) +
      geom_vline(xintercept = quantile(df$IGFBP4[df$hfail_inc==0], 0.95)) +
      scale_fill_manual(values=c("#9EAAC4","#3C5488"), labels=c("Controls", "Cases"), name="") +
      theme_classic() +
      labs(title="Heart failure", x = "Circulating IGFBP4 (SD)", y=" ")

      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_distributions_reb_1.tiff", width=2*1500, height=2*800, res=2*300)
        ggarrange(plot_hfail_1, plot_hfail_2, common.legend=T, ncol=2, nrow=1, legend="bottom")
        dev.off()
        
    quantile(df$TNFRSF10B[df$hfail_inc==0], 0.95) 
      percentile <- ecdf(df$TNFRSF10B[df$hfail_inc==1])
      1 - percentile(quantile(df$TNFRSF10B[df$hfail_inc==0], 0.95))
    quantile(df$IGFBP4[df$hfail_inc==0], 0.95)
      percentile <- ecdf(df$IGFBP4[df$hfail_inc==1])
      1 - percentile(quantile(df$IGFBP4[df$hfail_inc==0], 0.95))
    
    
    plot_cad_1 <- ggplot(df, aes(x=TNFRSF10B, fill=as.factor(cad_inc))) +
      geom_density(alpha=0.4) +
      scale_x_continuous(limits=c(-5,5)) +
      scale_fill_manual(values=c("#9EAAC4","#3C5488"), labels=c("Controls", "Cases"), name="") +
      geom_vline(xintercept = quantile(df$TNFRSF10B[df$cad_inc==0], 0.95)) +
      theme_classic() +
      labs(title="Coronary artery disease", x="Circulating TNFRSF10B (SD)", y = "Density")
   
      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_distributions_reb_2.tiff", width=2*750, height=2*800, res=2*300)
        ggarrange(plot_cad_1, common.legend=T, ncol=1, nrow=1, legend="bottom")
        dev.off()
        
    quantile(df$TNFRSF10B[df$cad_inc==0], 0.95) 
      percentile <- ecdf(df$TNFRSF10B[df$cad_inc==1])
      1 - percentile(quantile(df$TNFRSF10B[df$cad_inc==0], 0.95))
    
  ### 4 - GO enrichment ####
  
  library(data.table)
  library(tidyverse)
  library(ggtext)
  library(gtools)
  library(Hmisc)
  library(readxl)


  df_go <- read_excel(".../ukb_proteomics_cvd/output_files/2_GO.xlsx")
  df_go$Domain <- factor(df_go$Domain, levels=c("Biological process", "Molecular function", "Cellular component"))
  
  df_go$Ont_clean <- sapply( strwrap(capitalize(tolower(sub("\\s*\\(GO:[^)]+\\)", "", df_go$Ont))), 25, simplify=FALSE), paste, collapse="\n" )
  df_go$alpha <- ifelse(df_go$FDR_P<0.05, 0.9, 0.8)
  
  cad <- ggplot(data=df_go[df_go$Disease=="Coronary artery disease",], aes(x=fct_reorder(Ont_clean, desc(Number)), 
                         y=-log(P), fill=Domain)) +
    geom_col(aes(alpha=alpha)) +
    scale_alpha_continuous(guide=FALSE, range=c(0.5, 0.9)) + 
    scale_fill_manual(values=c("#F39B7F", "#B97E64", "#7E6148")) + 
    scale_y_continuous(limits = c(0,20)) +
    theme_classic() +
    labs(y="-log<sub>10</sub>(P)", x="") +
    coord_flip() +
    theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
  af <- ggplot(data=df_go[df_go$Disease=="Atrial fibrillation",], aes(x=fct_reorder(Ont_clean, desc(Number)), 
                         y=-log(P), fill=Domain)) +
    geom_col(aes(alpha=alpha)) +
    scale_alpha_continuous(guide=FALSE, range=c(0.5, 0.9)) + 
    scale_fill_manual(values=c("#F39B7F", "#B97E64", "#7E6148")) + 
    theme_classic() +
    scale_y_continuous(limits = c(0,20)) +
    labs(y="-log<sub>10</sub>(P)", x="") +
    coord_flip() +
    theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
  hf <- ggplot(data=df_go[df_go$Disease=="Heart failure",], aes(x=fct_reorder(Ont_clean, desc(Number)), 
                         y=-log(P), fill=Domain)) +
    geom_col(aes(alpha=alpha)) +
    scale_alpha_continuous(guide=FALSE, range=c(0.5, 0.9)) + 
    scale_fill_manual(values=c("#F39B7F", "#B97E64", "#7E6148")) + 
    theme_classic() +
    scale_y_continuous(limits = c(0,20)) +
    labs(y="-log<sub>10</sub>(P)", x="") +
    coord_flip() +
    theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
  as <- ggplot(data=df_go[df_go$Disease=="Aortic stenosis",], aes(x=fct_reorder(Ont_clean, desc(Number)), 
                         y=-log(P), fill=Domain)) +
    geom_col(aes(alpha=alpha)) +
    scale_alpha_continuous(guide=FALSE, range=c(0.5, 0.9)) + 
    scale_fill_manual(values=c("#F39B7F", "#B97E64", "#7E6148")) + 
    scale_y_continuous(limits = c(0,20)) +
    theme_classic() +
    labs(y="-log<sub>10</sub>(P)", x="") +
    coord_flip() +
    theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
  
  
  library(rvg)
  library(officer)
  
    Comp_1_dml <- dml(code = print(cad, newpage = FALSE))
    Comp_2_dml <- dml(code = print(hf, newpage = FALSE))
    Comp_3_dml <- dml(code = print(af, newpage = FALSE))
    Comp_4_dml <- dml(code = print(as, newpage = FALSE))

    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_1_dml, ph_location(left = 1, top = 1, width = 4, height = 7)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")
    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_2_dml, ph_location(left = 1, top = 1, width = 4, height = 7)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")
    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_3_dml, ph_location(left = 1, top = 1, width = 4, height = 7)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")
    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_4_dml, ph_location(left = 1, top = 1, width = 4, height = 7)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")

    rm(cad, af, hf, as, Comp_1_dml, Comp_2_dml, Comp_3_dml, Comp_4_dml)

  write.csv(df_go[,c(1,2,4,5)], ".../ukb_proteomics_cvd/output_files/source_data/ext_data_fig_4.csv", row.names=F)
  

  ### 5 - Sex-specific analyses ####

  library(ggtext)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
    
  df <- fread(".../ukb_proteomics_cvd/output_files/4_sex_strat_full.csv")
  
  df$diff_M_F <- log(df$HR_M)-log(df$HR_F)
  df <- df[order(df$diff_M_F),]
  df$num_all <- nrow(df):1
  df$num_peroutc <- (nrow(df)/4)+1-ave(df$num_all, df$Outcome, FUN = seq_along)
  
  # df_bonf <- df[df$P_Val_Int < 0.05/5836, ]
  # df_sign_ss <- df[(df$P_Val_F>=0.05&df$P_Val_M<0.05/(2*5836)) | (df$P_Val_M>=0.05&df$P_Val_F<0.05/(2*5836)), ]
  # df_sign_s_inv <- df[(df$P_Val_F<0.05&df$P_Val_M<0.05&(log(df$HR_M)*log(df$HR_F))<0), ]
  
  # plist <- c(df$P_Val_M, df$P_Val_F)
  # plist <- data.frame(Sex=c(rep("M", 5836), rep("F", 5836)), P=plist, P_adj=p.adjust(plist, "BH"))
  # df$P_Val_M_Adj <- plist$P_adj[1:5836]
  # df$P_Val_F_Adj <- plist$P_adj[5837:11672]
  df$col <- ifelse(((df$P_Val_F>0.05 | log(df$HR_F)*log(df$HR_M)<0) & df$P_Val_M<0.05/(5836) & df$P_Val_Int<0.05), "Male-specific biomarker", 
                          ifelse((df$P_Val_M>0.05 | log(df$HR_F)*log(df$HR_M)<0) & df$P_Val_F<0.05/(5836) & df$P_Val_Int<0.05, "Female-specific biomarker", "Other"))
  df$col_2 <- ifelse(((df$P_Val_F > 0.05 | log(df$HR_F)*log(df$HR_M)<0) & df$P_Val_M<0.05/(5836) & df$P_Val_Int<0.05), "Male-specific biomarker", 
                          ifelse((df$P_Val_M > 0.05 | log(df$HR_F)*log(df$HR_M)<0) & df$P_Val_F<0.05/(5836) & df$P_Val_Int<0.05, "Female-specific biomarker",
                          ifelse((df$P_Val_M < 0.05/(5836)) & df$P_Val_Int<0.05, "Stronger in males", 
                          ifelse((df$P_Val_F < 0.05/(5836)) & df$P_Val_Int<0.05, "Stronger in females", "Other"))))
  df$col_3 <- ifelse((((df$P_Val_F>0.05 | log(df$HR_F)*log(df$HR_M)<0) & df$P_Val_M<0.05/(5836)) | 
                        ((df$P_Val_M>0.05 | log(df$HR_F)*log(df$HR_M)<0) & df$P_Val_F<0.05/(5836))) & df$P_Val_Int<0.05 & df$diff_M_F>0, "Sex-specific biomarker (higher in males)", 
                          ifelse((((df$P_Val_F>0.05 | log(df$HR_F)*log(df$HR_M)<0) & df$P_Val_M<0.05/(5836)) | 
                                  ((df$P_Val_M>0.05 | log(df$HR_F)*log(df$HR_M)<0) & df$P_Val_F<0.05/(5836))) & df$P_Val_Int<0.05 & df$diff_M_F<0, "Sex-specific biomarker (higher in females)", 
                          ifelse((df$diff_M_F>0 & df$P_Val_Int<0.05 & (df$P_Val_M<0.05/(5836)|df$P_Val_F<0.05/(5836))), "Higher in males", 
                          ifelse(df$diff_M_F<0 & df$P_Val_Int<0.05 & (df$P_Val_M<0.05/(5836)|df$P_Val_F<0.05/(5836)), "Higher in females", "Other"))))

  write.csv(df[,c("Protein", "Outcome", "diff_M_F", "num_peroutc", "col_3")], ".../ukb_proteomics_cvd/output_files/source_data/fig_4_a.csv", row.names=F)
  
   a <- ggplot(df[df$Outcome=="CAD",], aes(x=num_peroutc, y=diff_M_F, color=col_3)) +
    geom_segment( aes(x=num_peroutc, xend=num_peroutc, y=0, yend=diff_M_F), color="grey90") +
    geom_point(data=df[df$Outcome=="CAD" & df$col_3=="Other",], color="grey90", size=2) +
    geom_point(data=df[df$Outcome=="CAD" & df$col_3!="Other",], aes(color=col_3), size=2) +
    geom_point(data=df[df$Outcome=="CAD" & df$col!="Other",], size=2) +
    scale_color_manual(values=c("#FDB7C7", "#BAB1D1", "grey90", "#fa6f8f", "#7563a3"), name="") +
    xlab("") +
    ylab("") +
    scale_y_continuous(limits = c(-0.45, 0.45), breaks=c(-0.4, -0.2, 0, 0.2, 0.4)) +
    theme_classic() +
    theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(), legend.position="none") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    ggtitle("Coronary artery disease") +
    geom_label_repel(aes(label=ifelse(df$col[df$Outcome=="CAD"]=="Other", "", df$Protein[df$Outcome=="CAD"])), max.overlaps=300, force=5, size=3, show.legend = F)
  b <- ggplot(df[df$Outcome=="HF",], aes(x=num_peroutc, y=diff_M_F, color=col_3)) +
    geom_segment( aes(x=num_peroutc, xend=num_peroutc, y=0, yend=diff_M_F), color="grey90") +
    geom_point(data=df[df$Outcome=="HF" & df$col_3=="Other",], color="grey90", size=2) +
    geom_point(data=df[df$Outcome=="HF" & df$col_3!="Other",], aes(color=col_3), size=2) +
    geom_point(data=df[df$Outcome=="HF" & df$col!="Other",], size=2) +
    scale_color_manual(values=c("#FDB7C7", "#BAB1D1", "grey90", "#fa6f8f", "#7563a3"), name="") +
    xlab("") +
    ylab("") +
    scale_y_continuous(limits = c(-0.45, 0.45), breaks=c(-0.4, -0.2, 0, 0.2, 0.4)) +
    theme_classic() +
    theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(), legend.position="none") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    ggtitle("Heart failure") +
    geom_label_repel(aes(label=ifelse(df$col[df$Outcome=="HF"]=="Other", "", df$Protein[df$Outcome=="HF"])), max.overlaps=300, force=2, size=3, show.legend = F)
  c <- ggplot(df[df$Outcome=="AF",], aes(x=num_peroutc, y=diff_M_F, color=col_3)) +
    geom_segment( aes(x=num_peroutc, xend=num_peroutc, y=0, yend=diff_M_F), color="grey90") +
    geom_point(data=df[df$Outcome=="AF" & df$col_3=="Other",], color="grey90", size=2) +
    geom_point(data=df[df$Outcome=="AF" & df$col_3!="Other",], aes(color=col_3), size=2) +
    geom_point(data=df[df$Outcome=="AF" & df$col!="Other",], size=2) +
    scale_color_manual(values=c("#FDB7C7", "#BAB1D1", "grey90", "#fa6f8f", "#7563a3"), name="") +
    xlab("") +
    ylab("") +
    scale_y_continuous(limits = c(-0.45, 0.45), breaks=c(-0.4, -0.2, 0, 0.2, 0.4)) +
    theme_classic() +
    theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(), legend.position="none") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    ggtitle("Atrial fibrillation") +
    geom_label_repel(aes(label=ifelse(df$col[df$Outcome=="AF"]=="Other", "", df$Protein[df$Outcome=="AF"])), max.overlaps=300, force=5, size=3, show.legend = F)
  d <- ggplot(df[df$Outcome=="AS",], aes(x=num_peroutc, y=diff_M_F, color=col_3)) +
    geom_segment( aes(x=num_peroutc, xend=num_peroutc, y=0, yend=diff_M_F), color="grey90") +
    geom_point(data=df[df$Outcome=="AS" & df$col_3=="Other",], color="grey90", size=2) +
    geom_point(data=df[df$Outcome=="AS" & df$col_3!="Other",], aes(color=col_3), size=2) +
    geom_point(data=df[df$Outcome=="AS" & df$col!="Other",], size=2) +
    scale_color_manual(values=c("#FDB7C7", "grey90", "#fa6f8f"), name="") +
    xlab("") +
    ylab("") +
    scale_y_continuous(limits = c(-0.45, 0.45), breaks=c(-0.4, -0.2, 0, 0.2, 0.4)) +
    theme_classic() +
    theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(), legend.position="none") +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    ggtitle("Aortic stenosis") +
    geom_label_repel(aes(label=ifelse(df$col[df$Outcome=="AS"]=="Other", "", df$Protein[df$Outcome=="AS"])), max.overlaps=300, force=5, size=3, show.legend = F)

  library(ggpubr)
  tiff(".../ukb_proteomics_cvd/output_files/4_diff_sexes_alloutcomes.tiff", width=3100, height=2200, res=300)
  ggarrange(a, b, c, d, common.legend=T, ncol=2, nrow=2)
  dev.off()
  # tiff(".../ukb_proteomics_cvd/output_files/4_diff_sexes_alloutcomes_fdr.tiff", width=3100, height=2200, res=300)
  # ggarrange(a, b, c, d, common.legend=T, ncol=2, nrow=2)
  # dev.off()
  
  
  library(rvg)
  library(officer)
  
    Comp_1_dml <- dml(code = print(ggarrange(a, b, c, d, common.legend=T, ncol=2, nrow=2), newpage = FALSE))

    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_1_dml, ph_location(left = 1, top = 1, width = 3100/300, height = 2200/300)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")

  
    # 5.A - Interaction - top 5 ####

  df <- df[order(df$P_Val_Int),]
  df$num_peroutc <- ave(df$P_Val_Int, df$Outcome, FUN = seq_along)
  df <- df[df$num_peroutc %in% 1:5,]
  
  df_long <- data.frame(protein=rep(df$Protein, 2), sex=c(rep("F", 20), rep("M", 20)), hr=c(df$HR_F, df$HR_M), lci=c(df$CI_Low_F, df$CI_Low_M),
                        uci=c(df$CI_High_F, df$CI_High_M), pval=c(df$P_Val_F, df$P_Val_M), rank=rep(df$num_peroutc, 2), outcome=rep(df$Outcome, 2))
  df_long <- df_long[order( df_long$outcome, df_long$sex, df_long$rank),]

  write.csv(df_long, ".../ukb_proteomics_cvd/output_files/source_data/fig_4_b.csv", row.names=F)
  
  a <- ggplot(data = df_long[df_long$outcome=="CAD",], aes(x = hr, y = as.factor(rank), color = sex, xmin=lci, xmax=uci)) +
      geom_pointrange(size=0.3, shape=20, position = position_dodge2(width = 0.65, padding = 0.5)) +
      scale_color_manual(values = c("#fa6f8f", "#7563a3"), labels=c("Female", "Male")) +
      labs(color = "") +
      ylab("") +
      xlab("HR (95% CI)") +
      ggtitle("Coronary artery disease") +
      scale_y_discrete(limits=rev, labels=rev(df_long[df_long$outcome=="CAD",]$protein)) +
      scale_x_continuous(trans="log", limits=c(0.7, 2.1), breaks=c(0.8, 1.2, 1.6, 2.0)) +
      geom_vline(xintercept = 1, linetype=2) +
      theme_classic() + theme(legend.position = "none")
  b <- ggplot(data = df_long[df_long$outcome=="HF",], aes(x = hr, y = as.factor(rank), color = sex, xmin=lci, xmax=uci)) +
      geom_pointrange(size=0.3, shape=20, position = position_dodge2(width = 0.65, padding = 0.5)) +
      scale_color_manual(values = c("#fa6f8f", "#7563a3"), labels=c("Female", "Male")) +
      labs(color = "") +
      ylab("") +
      xlab("HR (95% CI)") +
      ggtitle("Heart failure") +
      scale_y_discrete(limits=rev, labels=rev(df_long[df_long$outcome=="HF",]$protein)) +
      scale_x_continuous(trans="log", limits=c(0.7, 2.1), breaks=c(0.8, 1.2, 1.6, 2.0)) +
      geom_vline(xintercept = 1, linetype=2) +
      theme_classic() + theme(legend.position = "none")
  c <- ggplot(data = df_long[df_long$outcome=="AF",], aes(x = hr, y = as.factor(rank), color = sex, xmin=lci, xmax=uci)) +
      geom_pointrange(size=0.3, shape=20, position = position_dodge2(width = 0.65, padding = 0.5)) +
      scale_color_manual(values = c("#fa6f8f", "#7563a3"), labels=c("Female", "Male")) +
      labs(color = "") +
      ylab("") +
      xlab("HR (95% CI)") +
      ggtitle("Atrial fibrillation") +
      scale_y_discrete(limits=rev, labels=rev(df_long[df_long$outcome=="AF",]$protein)) +
      scale_x_continuous(trans="log", limits=c(0.7, 2.1), breaks=c(0.8, 1.2, 1.6, 2.0)) +
      geom_vline(xintercept = 1, linetype=2) +
      theme_classic() + theme(legend.position = "none")
  d <- ggplot(data = df_long[df_long$outcome=="AS",], aes(x = hr, y = as.factor(rank), color = sex, xmin=lci, xmax=uci)) +
      geom_pointrange(size=0.3, shape=20, position = position_dodge2(width = 0.65, padding = 0.5)) +
      scale_color_manual(values = c("#fa6f8f", "#7563a3"), labels=c("Female", "Male")) +
      labs(color = "") +
      ylab("") +
      xlab("HR (95% CI)") +
      ggtitle("Aortic stenosis") +
      scale_y_discrete(limits=rev, labels=rev(df_long[df_long$outcome=="AS",]$protein)) +
      scale_x_continuous(trans="log", limits=c(0.7, 2.1), breaks=c(0.8, 1.2, 1.6, 2.0)) +
      geom_vline(xintercept = 1, linetype=2) +
      theme_classic() + theme(legend.position = "none")

  library(ggpubr)
  tiff("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/4_diff_sexes_hrs.tiff", width=1650, height=1200, res=300)
  ggarrange(a,b,c,d, common.legend=T)
  dev.off()
  
  
  library(rvg)
  library(officer)
  
    Comp_1_dml <- dml(code = print(ggarrange(a,b,c,d, common.legend=T), newpage = FALSE))

    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_1_dml, ph_location(left = 1, top = 1, width = 1650/300, height = 1200/300)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")

  rm(list=ls())
  
  
    # 5.B - Correlation between sexes ####

  df <- fread(".../ukb_proteomics_cvd/output_files/4_sex_strat_full.csv")
  df$diff_M_F <- log(df$HR_M)-log(df$HR_F)
  df$col <- ifelse((df$P_Val_F>=0.05 & df$P_Val_M<0.05/(5836) & df$P_Val_Int<0.05), "Male-specific biomarker", 
                          ifelse(df$P_Val_M>=0.05 & df$P_Val_F<0.05/(5836) & df$P_Val_Int<0.05, "Female-specific biomarker", "Other"))
  df$sign <- ifelse(df$P_Val_M<0.05/(5836) | df$P_Val_F<0.05/(5836), T, F)
  
  cor(log(df$HR_M), log(df$HR_F), method="pearson")  # r=0.64
  cor(log(df$HR_M[df$Outcome=="CAD"]), log(df$HR_F[df$Outcome=="CAD"]), method="pearson")  # r=0.71
  cor(log(df$HR_M[df$Outcome=="HF"]), log(df$HR_F[df$Outcome=="HF"]), method="pearson")  # r=0.79
  cor(log(df$HR_M[df$Outcome=="AF"]), log(df$HR_F[df$Outcome=="AF"]), method="pearson")  # r=0.67
  cor(log(df$HR_M[df$Outcome=="AS"]), log(df$HR_F[df$Outcome=="AS"]), method="pearson")  # r=0.47
  
  df$Outcome <- factor(ifelse(df$Outcome=="CAD", "Coronary artery disease", 
                       ifelse(df$Outcome=="HF", "Heart failure", 
                       ifelse(df$Outcome=="AF", "Atrial fibrillation", 
                       ifelse(df$Outcome=="AS", "Aortic stenosis", NA)))), 
                       levels=c("Coronary artery disease", "Heart failure", "Atrial fibrillation", "Aortic stenosis"))
  
  write.csv(df[,c("Protein", "HR_M", "HR_F", "Outcome", "col")], ".../ukb_proteomics_cvd/output_files/source_data/ext_dat_fig_5.csv", row.names=F)
  
  a <- ggplot(data = df, aes(x = log(HR_M), y = log(HR_F), color = diff_M_F)) +
      geom_abline(intercept = 0, slope = 1, color="black") +
      geom_point(alpha=0.4) +
      ylab("log(HR) for females") +
      xlab("log(HR) for males") +
      facet_wrap(~Outcome) +
      scale_color_continuous(low="#7563a3", high="#fa6f8f") +
      theme_classic() + theme(legend.position = "none") + 
      theme(panel.grid.major = element_blank(),
        strip.text.x = element_text(hjust = 0, margin=margin(l=0)),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.text = element_text(size=12)) +
      geom_label_repel(aes(label=ifelse(df$col=="Other", "", df$Protein)), max.overlaps=300, force=5, size=3)
  
  tiff(".../ukb_proteomics_cvd/output_files/4_diff_sexes_correlation.tiff", width=3000, height=2500, res=300)
  a
  dev.off()  


### 6 - MR ####
  
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(tidyr)

  instr <- fread("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/MR_loop/2_correlatedinstruments/primaryanalyses_loop_corr_instr.csv")
  res <- fread("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/MR_loop/2_correlatedinstruments/primaryanalyses_loop_corr.csv")
  fstat <- fread("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/MR_loop/2_correlatedinstruments/rsq_fstat_combined_estimated_and_allelescore.csv")
    res <- res[!(res$exp %in% fstat$exposure[fstat$f<10]),]
    res <- res[res$method=="Inverse variance weighted (correlation inc)" | res$method=="Wald ratio",]
    res <- res[!duplicated( res[,c("exp", "outc")]),]
    res$outc <- factor(toupper(res$outc), levels=c("CAD", "HF", "AF", "AS"))
    res$col <- ifelse(res$pval<0.05 & res$b<0, "Protective", ifelse(res$pval<0.05 & res$b>0, "Detrimental", "Other"))
    # write.csv(res, ".../ukb_proteomics_cvd/output_files/3_MR_allfindings.csv")
    # res_sign <- res[res$col!="Other",]
    # write.csv(res_sign, ".../ukb_proteomics_cvd/output_files/3_MR_signfindings.csv")
    # rm(res_sign)
    rm(fstat, instr)
  
  sign_prim <- fread('.../ukb_proteomics_cvd/output_files/1_adj_no_prev_timevar_model_hr.csv')
    sign_prim <- sign_prim[sign_prim$P_Value<0.05/5836,]
    sign_prim$Outcome <- ifelse(sign_prim$Outcome == "Afib", "AF", sign_prim$Outcome)
    sign_prim$beta_prim_obs <- log(sign_prim$HR)
    res <- merge(res, sign_prim[,c("Protein", "Outcome", "beta_prim_obs")], by.x=c("exp", "outc"), by.y=c("Protein", "Outcome"), all.x=T)
    rm(sign_prim)
    
  sens <- fread("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/MR_loop/2_correlatedinstruments/sensitivityanalyses_loop_corr.csv")
    sens$outc <- toupper(sens$outc)
    sens <- merge(sens, res[res$pval<0.05, c("exp", "outc")], by=c("exp", "outc"), all.x=F, all.y=F)
    sens_1 <- sens[sens$method=="Inverse variance weighted (correlation inc)" | sens$method=="Wald ratio",]
    sens_2 <- sens[sens$method=="Egger (correlation inc)" & sens$pvalthreshold==5e-06 & sens$rsqthreshold==0.1,]
    sens_3 <- sens[sens$method=="Egger intercept (correlation inc)" & sens$pvalthreshold==5e-06 & sens$rsqthreshold==0.1,]
    sens <- rbind(sens_1, sens_2)
    sens <- sens[!duplicated( sens[,c("exp", "outc","pvalthreshold", "rsqthreshold", "method")]),]
    sens_3 <- merge(sens_3, sens_2[,c("exp", "outc", "pval")], all.x=T, all.y=T, by=c("exp", "outc"))
    sens_3 <- sens_3[!duplicated(sens_3[,c("exp", "outc","pvalthreshold", "rsqthreshold", "method")]),]
    rm(sens_1, sens_2)

  sens_4 <- fread("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/MR_loop/2_correlatedinstruments/sensitivityanalysis_allelescore_oneproteinpermodel.csv")
    colnames(sens_4) <- c("exp", "outc", "b", "se", "pval") 
    sens_4$V1 <- NA
    sens_4$pvalthreshold <- "Allele score"
    sens_4$rsqthreshold <- "Allele score"
    sens_4$nsnp <- 0
    sens_4$method <- "Allele score"
    sens_4$outc <- toupper(ifelse(sens_4$outc=="hfail", "hf", ifelse(sens_4$outc=="afib", "af", ifelse(sens_4$outc=="ao_sten", "as", sens_4$outc))))
    sens <- rbind(sens, sens_4)

  sens$method <- ifelse(sens$method=="Wald ratio", "Wald", ifelse(sens$method=="Inverse variance weighted (correlation inc)", "IVW", 
                        ifelse(sens$method=="Egger (correlation inc)", "Egg", 
                        ifelse(sens$method=="Allele score", "All", NA))))
    sens_w <- sens[,-c("V1", "nsnp", "se", "pval")] %>% pivot_wider(names_from=c(pvalthreshold, rsqthreshold, method), values_from=b)
    # sens_w2 <- sens
    # sens_w2$method <- ifelse(sens_w2$method %in% c("Wald", "IVW"), "Prim", sens_w2$method)
    # sens_w2 <- sens_w2[,-c("V1")] %>% pivot_wider(names_from=c(pvalthreshold, rsqthreshold, method), values_from=c(nsnp, b, se, pval))
    # write.csv(sens_w2, ".../ukb_proteomics_cvd/output_files/3_MR_sens_formatted.csv")
    sens_w <- as.data.frame(sens_w)
    sens_w$prim_analysis <- ifelse(is.na(sens_w$`5e-06_0.1_IVW`), sens_w$`5e-06_0.1_Wald`, sens_w$`5e-06_0.1_IVW`)
    sens_w <- sens_w[,-c(which(colnames(sens_w)=="5e-06_0.1_IVW"), which(colnames(sens_w)=="5e-06_0.1_Wald"))]
    
    directional_consistency_df <- sens_w
    directional_consistency_df$outc <- toupper(directional_consistency_df$outc)
    for (col in colnames(sens_w)[-c(which(colnames(sens_w)=="exp"), which(colnames(sens_w)=="outc"))]) {
      for (i in 1:nrow(sens_w)) {
        directional_consistency_df[i, col] <- sign(sens_w[i, col]) == sign(sens_w[i, "prim_analysis"])
      }
    }
    for (i in 1:nrow(directional_consistency_df)) {
      directional_consistency_df[i, "non_na_count"] <- sum(!is.na(directional_consistency_df[i, -c(which(colnames(directional_consistency_df)=="exp"), which(colnames(directional_consistency_df)=="outc"), 
                                                                                                   which(colnames(directional_consistency_df)=="prim_analysis"), which(colnames(directional_consistency_df)=="5e-06_0.1_Egg"))]))
      directional_consistency_df[i, "sum_same_columns"] <- sum(directional_consistency_df[i, -c(which(colnames(directional_consistency_df)=="exp"), which(colnames(directional_consistency_df)=="outc"), 
                                                                                                which(colnames(directional_consistency_df)=="prim_analysis"), which(colnames(directional_consistency_df)=="5e-06_0.1_Egg"), 
                                                                                                which(colnames(directional_consistency_df)=="non_na_count"))], na.rm=T)
      directional_consistency_df[i, "ratio"] <- directional_consistency_df[i, "sum_same_columns"] / directional_consistency_df[i, "non_na_count"]
    }
    directional_consistency_df$all_sens <- ifelse(directional_consistency_df$ratio==1 & (is.na(directional_consistency_df$"5e-06_0.1_Egg")|(directional_consistency_df$"5e-06_0.1_Egg"==1)),
                                                  "Yes", "No")
    rm(sens)                                              
    
    res <- merge(res, directional_consistency_df[,c("exp", "outc", "all_sens")], by=c("exp", "outc"), all.x=T)
      res$exp_outc <- paste0(res$exp, "_", res$outc)
      sens_3$exp_outc <- paste0(sens_3$exp, "_", sens_3$outc)
      res$egger_int_viol <- ifelse(res$exp_outc %in% sens_3$exp_outc[sens_3$pval.x>0.05], "No violation",
                                   ifelse(res$exp_outc %in% sens_3$exp_outc[sens_3$pval.x<0.05 & sens_3$pval.y<0.05], "No violation",
                                   ifelse(res$exp_outc %in% sens_3$exp_outc[sens_3$pval.x<0.05], "Violation", NA
                                   )))
      res$all_sens <- ifelse(is.na(res$egger_int_viol), res$all_sens, ifelse(res$egger_int_viol=="Violation", "No", res$all_sens))
      res$all_sens <- ifelse(is.na(res$all_sens), "Nonsignificant", res$all_sens)
      
    # write.csv(res, "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/MR_loop/2_correlatedinstruments/sensitivityanalyses_posttwosamplesens_postallele.csv")
      
  all_mv <- fread("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/MR_loop/2_correlatedinstruments/regression_prs_proteins_output_multiv_analyses.csv")
    colnames(all_mv) <- c("exp", "outc", "b_allelescore_mv", "se", "pval")
    all_mv$exp <- gsub("prs_", "", all_mv$exp)
    all_mv$outc <- ifelse(all_mv$outc=="cad", "CAD", ifelse(all_mv$outc=="hfail", "HF", ifelse(all_mv$outc=="afib", "AF", "AS")))
  
  res <- merge(res, all_mv[,c("exp", "outc", "b_allelescore_mv")], by=c("exp", "outc"), all.x=T, all.y=F)
    res$sens_robust_allelescore_mv <- ifelse(res$b_allelescore_mv*res$b > 0, T, F)

  res$all_sens <- ifelse(is.na(res$sens_robust_allelescore_mv), res$all_sens, ifelse(res$sens_robust_allelescore_mv==F, "No", res$all_sens))
    res$directional_obs <- ifelse(res$b*res$beta_prim_obs>0, "Consistent", "Not consistent")
    res$group <- ifelse(res$directional_obs == "Consistent" & res$all_sens == "Yes", "Robust / consistent with\nobservational findings",
                          ifelse(res$directional_obs == "Not consistent" & res$all_sens == "Yes", "Robust / not consistent\nwith observational findings",
                          ifelse(res$directional_obs == "Consistent" & res$all_sens == "No", "Not robust / consistent with\nobservational findings",
                          ifelse(res$directional_obs == "Not consistent" & res$all_sens == "No", "Not robust / not consistent\nwith observational findings",
                          "Nonsignificant"))))
    # write.csv(res, "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/MR_loop/2_correlatedinstruments/sensitivityanalyses_posttwosamplesens_postallele_postallelemv.csv")
    # write.csv(res[res$all_sens=="Yes", c("exp", "outc")], ".../ukb_proteomics_cvd/output_files/3_MR_robust_findings.csv", row.names=F)
    res$exp <- ifelse(res$exp=="NTproBNP", "NT-proBNP", res$exp)
  rm(directional_consistency_df, fstat, instr, sens_3, sens_4)
      
  write.csv(res[,c("exp", "outc", "b", "pval", "group")], ".../ukb_proteomics_cvd/output_files/source_data/fig_3.csv", row.names=F)
  
  a <- ggplot(res[res$outc=="CAD",], aes(b, -log10(pval), color=group)) +
      scale_color_manual(values=c("grey90", "#FFF1B9", "#ECB5A2", "#E8B900", "red3")) +
      labs(y="-log<sub>10</sub>(*P*)", x="") +
      geom_point()+
      theme_classic() +
      scale_y_continuous(limits=c(0,4.5)) +
      scale_x_continuous(limits=c(-1.2, 1.2)) +
      guides(fill=guide_legend(override.aes = aes(label = NA), ncol=2, nrow=2)) +
      ggtitle("Coronary artery disease") +
      theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(), legend.position="none") +
      geom_text_repel(aes(label=ifelse(res[res$outc=="CAD",]$all_sens == "Yes" & (res[res$outc=="CAD",]$pval<0.05), 
                                       res[res$outc=="CAD",]$exp, "")), size=3, force=25, max.overlaps=100, show.legend = FALSE)
  b <- ggplot(res[res$outc=="HF",], aes(b, -log10(pval), color=group)) +
      scale_color_manual(values=c("grey90", "#FFF1B9", "#ECB5A2", "#E8B900", "red3"), name="") +
      labs(y="", x="") +
      geom_point()+
      theme_classic() +
      scale_y_continuous(limits=c(0,4.5)) +
      scale_x_continuous(limits=c(-1.2, 1.2)) +
      guides(fill=guide_legend(override.aes = aes(label = NA), ncol=2, nrow=2)) +
      ggtitle("Heart failure") +
      theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),legend.position = "bottom") +
      geom_text_repel(aes(label=ifelse(res[res$outc=="HF",]$all_sens == "Yes" & (res[res$outc=="HF",]$pval<0.05), 
                                       res[res$outc=="HF",]$exp, "")), size=3, force=25, max.overlaps=100, show.legend = FALSE)
  c <- ggplot(res[res$outc=="AF",], aes(b, -log10(pval), color=group)) +
      scale_color_manual(values=c("grey90", "#FFF1B9", "#ECB5A2", "#E8B900", "red3"), name="") +
      labs(y="-log<sub>10</sub>(*P*)", x="log(OR)") +
      geom_point()+
      theme_classic() + 
      scale_y_continuous(limits=c(0,4.5)) +
      scale_x_continuous(limits=c(-1.2, 1.2)) +
      guides(fill=guide_legend(override.aes = aes(label = NA), ncol=2, nrow=2)) +
      ggtitle("Atrial fibrillation") +
      theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(), legend.position="none") +
      geom_text_repel(aes(label=ifelse(res[res$outc=="AF",]$all_sens == "Yes" & res[res$outc=="AF",]$pval<0.05, 
                                       res[res$outc=="AF",]$exp, "")), size=3, force=25, max.overlaps=100, show.legend = FALSE)
  d <- ggplot(res[res$outc=="AS",], aes(b, -log10(pval), color=group)) +
      scale_color_manual(values=c("grey90", "#E8B900"), name="Outcome") +
      labs(y="", x="log(OR)") +
      geom_point()+
      theme_classic() +
      guides(fill=guide_legend(override.aes = aes(label = NA), ncol=2, nrow=2)) +
      scale_y_continuous(limits=c(0,4.5)) +
      scale_x_continuous(limits=c(-1.2, 1.2)) +
      ggtitle("Aortic stenosis") +
      theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(), legend.position="none") +
      geom_text_repel(aes(label=ifelse(res[res$outc=="AS",]$all_sens == "Yes" & res[res$outc=="AS",]$pval<0.05, 
                                       res[res$outc=="AS",]$exp, "")), size=3, force=25, max.overlaps=100, show.legend = FALSE)
  legend <- cowplot::get_plot_component(b, 'guide-box-bottom', return_all = TRUE)
  b <- b + theme(legend.position="none")
  
  library(ggpubr)
  tiff(".../ukb_proteomics_cvd/output_files/3_MR.tiff", width=6200, height=4400, res=600)
  ggarrange(ggarrange(a, b, c, d, common.legend=F, ncol=2, nrow=2), cowplot::plot_grid(legend), nrow=2, ncol=1, heights=c(1,0.1))
  dev.off()

  pdf(".../ukb_proteomics_cvd/output_files/Fig. 3.pdf", width=6200/600, height=4400/600)
  ggarrange(ggarrange(a, b, c, d, common.legend=F, ncol=2, nrow=2), cowplot::plot_grid(legend), nrow=2, ncol=1, heights=c(1,0.1))
  dev.off()
  
  
    # 6.A - Druggability MR findings ####

  library(data.table)
  
  mr <- fread(".../ukb_proteomics_cvd/output_files/3_MR_robust_findings.csv")
  mr$outc <- ifelse(mr$outc=="CAD", "Coronary artery disease", 
                      ifelse(mr$outc=="HF", "Heart failure", 
                      ifelse(mr$outc=="AF", "Atrial fibrillation", "Aortic stenosis"
                             )))
  drug <- fread(".../ukb_proteomics_cvd/input_files/finan_druggability.csv")
  
  mr <- merge(mr, drug, by.x="exp", by.y="hgnc_names", all.x=T, all.y=F)
  mr$exp_outc <- paste0(mr$exp, "_", mr$outc)

  obs_est <- fread(".../ukb_proteomics_cvd/output_files/3_druggability_robust_obs_genet/obs_estimates.csv")
    obs_est$exp <- ifelse(obs_est$exp=="NT-proBNP", "NTproBNP", obs_est$exp)
  mr_est <- fread(".../ukb_proteomics_cvd/output_files/3_druggability_robust_obs_genet/primary_mr_estimates.csv")
  obs_est$exp_outc <- paste0(obs_est$exp, "_", obs_est$outc)
  mr_est$exp_outc <- paste0(mr_est$exp, "_", mr_est$outc)
  obs_est <- obs_est[obs_est$exp_outc %in% mr$exp_outc,]
  mr_est <- mr_est[mr_est$exp_outc %in% mr$exp_outc,]
  colnames(obs_est) <- c("exp", "outc", "beta_obs", "exp_outc")
  colnames(mr_est) <- c("exp", "outc", "beta_mr", "exp_outc")
  est <- merge(obs_est[,c("beta_obs", "exp_outc")], mr_est[,c("beta_mr", "exp_outc")], by.x=c("exp_outc"), by.y=c("exp_outc"))
  est$cons <- ifelse(est$beta_obs*est$beta_mr>0, "Yes", "No")
  
  mr <- merge(mr, est[,c("exp_outc", "cons")], by="exp_outc")
  
  write.csv(mr, ".../ukb_proteomics_cvd/output_files/3_MR_robust_findings_druggability.csv", row.names=F)
  
  rm(drug, est, mr, mr_est, obs_est)
  
  
  
  ### 7 - Coloc ####
  
  library(data.table)
  mr <- fread(".../ukb_proteomics_cvd/output_files/3_MR_robust_findings.csv")
    coloc <- fread("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/MR_loop/2_correlatedinstruments/coloc_loop.csv")
    coloc$outc <- toupper(coloc$outc)
    coloc <- merge(coloc, mr, by=c("exp", "outc"))
    rm(mr)

  write.csv(coloc, ".../ukb_proteomics_cvd/output_files/3_coloc_robust_findings.csv", row.names=F)
      
  
  ### 8 - Prediction ####
    # 8.A - UKB - original ####
      # 8.A.1 - Regression coefficients / weights ####
        # 8.A.1a - CAD ####
        
      cad_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_cl.csv")
          cad_feat_cl <- cad_feat_cl[cad_feat_cl$s1!=0 & cad_feat_cl$V1!="(Intercept)",]
          cad_feat_cl$sd <- c(1, 8.2, 1, 1, 1, 1, 1, 4.7, 19.6, 43.8, 14.9, 0.19)
          cad_feat_cl$s1_weighted <- abs(cad_feat_cl$s1*cad_feat_cl$sd)
          cad_feat_cl$name_form <- c("Sex (male vs. female)", "Age (per SD)", "Smoking (ever vs. never)", "Cholesterol-lowering medication use (yes vs. no)", 
                                     "Race/ethnicity (non-White vs. White)", "Antihypertensive medication use (yes vs. no)", "Type 2 diabetes mellitus (yes vs. no)", 
                                     "BMI (per SD)", "Systolic blood pressure (per SD)", "Total cholesterol (per SD)", "HDL cholesterol (per SD)", "Creatinine (per SD)")
          cad_feat_cl <- cad_feat_cl[order(cad_feat_cl$s1_weighted),]
          cad_feat_cl$name_form <- sapply( strwrap(cad_feat_cl$name_form, 22, simplify=FALSE), paste, collapse="\n" )
      cad_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_pr.csv")
          cad_feat_pr <- cad_feat_pr[cad_feat_pr$s1!=0 & cad_feat_pr$V1!="(Intercept)",]
          cad_feat_pr$s1_abs <- abs(cad_feat_pr$s1)
          cad_feat_pr <- cad_feat_pr[order(cad_feat_pr$s1_abs),]
          cad_feat_pr$V1 <- ifelse(cad_feat_pr$V1=="NTproBNP", "NT-proBNP", cad_feat_pr$V1)
          
      cadplot_cl <- ggplot(data=cad_feat_cl, aes(x=fct_reorder(name_form, s1_weighted), y=s1_weighted)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("CAD") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
      cadplot_pr <- ggplot(data=cad_feat_pr, aes(x=fct_reorder(V1, s1_abs), y=s1_abs)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("Coronary artery disease") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),axis.text.y = element_text(size = 5))
      rm(cad_feat_cl)
       
        
        # 8.A.1b - HF ####
        
      hfail_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_cl.csv")
          hfail_feat_cl <- hfail_feat_cl[hfail_feat_cl$s1!=0 & hfail_feat_cl$V1!="(Intercept)",]
          hfail_feat_cl$sd <- c(1, 8.2, 1, 1, 1, 1, 4.7, 19.6, 43.8, 0.19)
          hfail_feat_cl$s1_weighted <- abs(hfail_feat_cl$s1*hfail_feat_cl$sd)
          hfail_feat_cl$name_form <- c("Sex (male vs. female)", "Age (per SD)", "Smoking (ever vs. never)", "Cholesterol-lowering medication use (yes vs. no)", 
                                     "Antihypertensive medication use (yes vs. no)", "Type 2 diabetes mellitus (yes vs. no)", 
                                     "BMI (per SD)", "Systolic blood pressure (per SD)", "Total cholesterol (per SD)", "Creatinine (per SD)")
          hfail_feat_cl <- hfail_feat_cl[order(hfail_feat_cl$s1_weighted),]
          hfail_feat_cl$name_form <- sapply( strwrap(hfail_feat_cl$name_form, 22, simplify=FALSE), paste, collapse="\n" )
      hfail_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_pr.csv")
          hfail_feat_pr <- hfail_feat_pr[hfail_feat_pr$s1!=0 & hfail_feat_pr$V1!="(Intercept)",]
          hfail_feat_pr$s1_abs <- abs(hfail_feat_pr$s1)
          hfail_feat_pr <- hfail_feat_pr[order(hfail_feat_pr$s1_abs),]
          hfail_feat_pr$V1 <- ifelse(hfail_feat_pr$V1=="NTproBNP", "NT-proBNP", hfail_feat_pr$V1)
          
      hfailplot_cl <- ggplot(data=hfail_feat_cl, aes(x=fct_reorder(name_form, s1_weighted), y=s1_weighted)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("HF") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
      hfailplot_pr <- ggplot(data=hfail_feat_pr, aes(x=fct_reorder(V1, s1_abs), y=s1_abs)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("Heart failure") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),axis.text.y = element_text(size = 5))
      rm(hfail_feat_cl)
                      
        
        # 8.A.1c - AF ####
        
      afib_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_cl.csv")
          afib_feat_cl <- afib_feat_cl[afib_feat_cl$s1!=0 & afib_feat_cl$V1!="(Intercept)",]
          afib_feat_cl$sd <- c(1, 8.2, 1, 1, 1, 1, 1, 4.7, 19.6, 43.8, 14.9, 0.19)
          afib_feat_cl$s1_weighted <- abs(afib_feat_cl$s1*afib_feat_cl$sd)
          afib_feat_cl$name_form <- c("Sex (male vs. female)", "Age (per SD)", "Smoking (ever vs. never)", "Cholesterol-lowering medication use (yes vs. no)", 
                                     "Race/ethnicity (non-White vs. White)", "Antihypertensive medication use (yes vs. no)", "Type 2 diabetes mellitus (yes vs. no)", 
                                     "BMI (per SD)", "Systolic blood pressure (per SD)", "Total cholesterol (per SD)", "HDL cholesterol (per SD)", "Creatinine (per SD)")
          afib_feat_cl <- afib_feat_cl[order(afib_feat_cl$s1_weighted),]
          afib_feat_cl$name_form <- sapply( strwrap(afib_feat_cl$name_form, 22, simplify=FALSE), paste, collapse="\n" )
      afib_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_pr.csv")
          afib_feat_pr <- afib_feat_pr[afib_feat_pr$s1!=0 & afib_feat_pr$V1!="(Intercept)",]
          afib_feat_pr$s1_abs <- abs(afib_feat_pr$s1)
          afib_feat_pr <- afib_feat_pr[order(afib_feat_pr$s1_abs),]
          afib_feat_pr$V1 <- ifelse(afib_feat_pr$V1=="NTproBNP", "NT-proBNP", afib_feat_pr$V1)
         
      afibplot_cl <- ggplot(data=afib_feat_cl, aes(x=fct_reorder(name_form, s1_weighted), y=s1_weighted)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("AF") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
      afibplot_pr <- ggplot(data=afib_feat_pr, aes(x=fct_reorder(V1, s1_abs), y=s1_abs)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("Atrial fibrillation") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),axis.text.y = element_text(size = 5))
      rm(afib_feat_cl)
       
      
        # 8.A.1d - AS ####
        
      ao_sten_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_cl.csv")
          ao_sten_feat_cl <- ao_sten_feat_cl[ao_sten_feat_cl$s1!=0 & ao_sten_feat_cl$V1!="(Intercept)",]
          ao_sten_feat_cl$sd <- c(1, 8.2, 1, 1, 1, 1, 4.7, 19.6, 43.8, 14.9, 0.19)
          ao_sten_feat_cl$s1_weighted <- abs(ao_sten_feat_cl$s1*ao_sten_feat_cl$sd)
          ao_sten_feat_cl$name_form <- c("Sex (male vs. female)", "Age (per SD)", "Smoking (ever vs. never)", "Cholesterol-lowering medication use (yes vs. no)", 
                                     "Race/ethnicity (non-White vs. White)", "Antihypertensive medication use (yes vs. no)", 
                                     "BMI (per SD)", "Systolic blood pressure (per SD)", "Total cholesterol (per SD)", "HDL cholesterol (per SD)", "Creatinine (per SD)")
          ao_sten_feat_cl <- ao_sten_feat_cl[order(ao_sten_feat_cl$s1_weighted),]
          ao_sten_feat_cl$name_form <- sapply( strwrap(ao_sten_feat_cl$name_form, 22, simplify=FALSE), paste, collapse="\n" )
      ao_sten_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_pr.csv")
          ao_sten_feat_pr <- ao_sten_feat_pr[ao_sten_feat_pr$s1!=0 & ao_sten_feat_pr$V1!="(Intercept)",]
          ao_sten_feat_pr$s1_abs <- abs(ao_sten_feat_pr$s1)
          ao_sten_feat_pr <- ao_sten_feat_pr[order(ao_sten_feat_pr$s1_abs),]
          ao_sten_feat_pr$V1 <- ifelse(ao_sten_feat_pr$V1=="NTproBNP", "NT-proBNP", ao_sten_feat_pr$V1)
          
      ao_stenplot_cl <- ggplot(data=ao_sten_feat_cl, aes(x=fct_reorder(name_form, s1_weighted), y=s1_weighted)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("AS") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
      ao_stenplot_pr <- ggplot(data=ao_sten_feat_pr, aes(x=fct_reorder(V1, s1_abs), y=s1_abs)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("Aortic stenosis") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),axis.text.y = element_text(size = 5))
      rm(ao_sten_feat_cl)
       
      
        # 8.A.1e - Pics ####
    
      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_all_cl.tiff", width=5300, height=1800, res=300)
        ggarrange(cadplot_cl, hfailplot_cl, afibplot_cl, ao_stenplot_cl, ncol=4, nrow=1)
        dev.off()
      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_all_pr.tiff", width=5300, height=1800, res=300)
        ggarrange(cadplot_pr, hfailplot_pr, afibplot_pr, ao_stenplot_pr, ncol=4, nrow=1)
        dev.off()
      rm(cadplot_cl, hfailplot_cl, afibplot_cl, ao_stenplot_cl, cadplot_pr, hfailplot_pr, afibplot_pr, ao_stenplot_pr)
    
      cad_feat_pr$outc <- "Coronary artery disease"
      hfail_feat_pr$outc <- "Heart failure"
      afib_feat_pr$outc <- "Atrial fibrillation"
      ao_sten_feat_pr$outc <- "Aortic stenosis"
        
      write.csv(rbind(cad_feat_pr, hfail_feat_pr, afib_feat_pr, ao_sten_feat_pr), ".../ukb_proteomics_cvd/output_files/source_data/ext_dat_fig_6.csv", row.names=F)
      
      
      # 8.A.2 - ROC ####
    
      roc_cad_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_cad_cl.csv")
        roc_cad_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_cad_pr.csv")
        roc_cad_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_cad_comb.csv")
        roc_cad_cl$model <- "Clinical parameters"
        roc_cad_pr$model <- "Proteins"
        roc_cad_comb$model <- "Combined"
        roc_cad <- rbind(roc_cad_cl, roc_cad_pr, roc_cad_comb)
        rm(roc_cad_cl, roc_cad_pr, roc_cad_comb)
        
      roc_hfail_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_hfail_cl.csv")
        roc_hfail_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_hfail_pr.csv")
        roc_hfail_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_hfail_comb.csv")
        roc_hfail_cl$model <- "Clinical parameters"
        roc_hfail_pr$model <- "Proteins"
        roc_hfail_comb$model <- "Combined"
        roc_hfail <- rbind(roc_hfail_cl, roc_hfail_pr, roc_hfail_comb)
        rm(roc_hfail_cl, roc_hfail_pr, roc_hfail_comb)
      
      roc_afib_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_afib_cl.csv")
        roc_afib_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_afib_pr.csv")
        roc_afib_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_afib_comb.csv")
        roc_afib_cl$model <- "Clinical parameters"
        roc_afib_pr$model <- "Proteins"
        roc_afib_comb$model <- "Combined"
        roc_afib <- rbind(roc_afib_cl, roc_afib_pr, roc_afib_comb)
        rm(roc_afib_cl, roc_afib_pr, roc_afib_comb)
      
      roc_ao_sten_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_ao_sten_cl.csv")
        roc_ao_sten_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_ao_sten_pr.csv")
        roc_ao_sten_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/3_roc_rawdata_ao_sten_comb.csv")
        roc_ao_sten_cl$model <- "Clinical parameters"
        roc_ao_sten_pr$model <- "Proteins"
        roc_ao_sten_comb$model <- "Combined"
        roc_ao_sten <- rbind(roc_ao_sten_cl, roc_ao_sten_pr, roc_ao_sten_comb)
        rm(roc_ao_sten_cl, roc_ao_sten_pr, roc_ao_sten_comb)
      
      
      rocplot_cad <- ggplot(roc_cad, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "True positive rate") +
                        ggtitle("Coronary artery disease") +
                        coord_equal() +
                        theme_classic()
      rocplot_hfail <- ggplot(roc_hfail, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "") +
                        ggtitle("Heart failure") +
                        coord_equal() +
                        theme_classic()
      rocplot_af <- ggplot(roc_afib, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "") +
                        ggtitle("Atrial fibrillation") +
                        coord_equal() +
                        theme_classic()
      rocplot_ao_sten <- ggplot(roc_ao_sten, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "") +
                        ggtitle("Aortic stenosis") +
                        coord_equal() +
                        theme_classic()
      
      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_rocs.tiff", width=4500, height=1300, res=300)
        ggarrange(rocplot_cad, rocplot_hfail, rocplot_af, rocplot_ao_sten, common.legend=T, ncol=4, nrow=1, legend="bottom")
        dev.off()
      fig5_d <- ggarrange(rocplot_cad, rocplot_hfail, rocplot_af, rocplot_ao_sten, common.legend=T, ncol=4, nrow=1, legend="bottom")

      roc_cad$outc <- "Coronary artery disease"
      roc_hfail$outc <- "Heart failure"
      roc_afib$outc <- "Atrial fibrillation"
      roc_ao_sten$outc <- "Aortic stenosis"
        
      write.csv(rbind(roc_cad, roc_hfail, roc_afib, roc_ao_sten), ".../ukb_proteomics_cvd/output_files/source_data/fig_5_a.csv", row.names=F)
 
        
        # 8.A.2a - DR / FPR (crude) ####
        
      roc_sum <- data.frame(outcome=NA, model=NA, FPR=NA, TPR=NA)[-1,]
      
      for (j in c("afib", "cad", "hfail", "ao_sten")) {
        for (k in c("Clinical parameters", "Proteins", "Combined")) {
          for (i in c(0.01, 0.05, 0.1)) {
      
          roc_df <- get(paste0("roc_", j))
          roc_df <- roc_df[roc_df$model == k,]
          roc_df <- roc_df[roc_df$FPR < i,]
          
          roc_sum_row <- data.frame(outcome=j, model=k, FPR=i, TPR=roc_df$TPR[nrow(roc_df)])
          roc_sum <- rbind(roc_sum, roc_sum_row)
          
          }
        }
      }
        
      write.csv(roc_sum, ".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_fpr_tpr_models.csv", row.names=F)
      
      rocs <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_roc_full.csv")
      rm(rocplot_af, rocplot_cad, rocplot_hfail, rocplot_ao_sten, roc_afib, roc_cad, roc_hfail, roc_ao_sten, rocs, roc_sum, roc_sum_row)
      
      
       # 8.A.3 - KM ####
  
    df <- fread('.../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz')
      set.seed(1234)
      df <- df[sample(1:nrow(df)), ]
      set.seed(1)
      df$id <- 1:nrow(df)
      train <- df %>% dplyr::sample_frac(0.80)
      test  <- dplyr::anti_join(df, train, by = 'id')
      rm(df)
      
    cad_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_pr.csv")
      hfail_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_pr.csv")
      afib_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_pr.csv")
      ao_sten_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_pr.csv")
    test$cad_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(cad_feat_pr$s1))
      test$hfail_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(hfail_feat_pr$s1))
      test$afib_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(afib_feat_pr$s1))
      test$ao_sten_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(ao_sten_feat_pr$s1))
    test$cad_protriskscore_bin <- exp(test$cad_protriskscore)/(1+exp(test$cad_protriskscore))
      test$hfail_protriskscore_bin <- exp(test$hfail_protriskscore)/(1+exp(test$hfail_protriskscore))
      test$afib_protriskscore_bin <- exp(test$afib_protriskscore)/(1+exp(test$afib_protriskscore))
      test$ao_sten_protriskscore_bin <- exp(test$ao_sten_protriskscore)/(1+exp(test$ao_sten_protriskscore))
    test$cad_protriskscore_quintile <- factor(ntile(test$cad_protriskscore_bin, 5), levels=c(5:1))
      test$hfail_protriskscore_quintile <- factor(ntile(test$hfail_protriskscore_bin, 5), levels=c(5:1))
      test$afib_protriskscore_quintile <- factor(ntile(test$afib_protriskscore_bin, 5), levels=c(5:1))
      test$ao_sten_protriskscore_quintile <- factor(ntile(test$ao_sten_protriskscore_bin, 5), levels=c(5:1))
      
    surv_object <- Surv(time = test$cad_fu, event = test$cad_inc)
    fit1 <- survfit(surv_object ~ cad_protriskscore_quintile, data=test)
    survplot_cad <- ggsurvplot(fit1, pval=F, fun='event', color = "cad_protriskscore_quintile", ylim=c(0,0.22), legend=c(0.21,0.76), palette=c("#A35C14", "#E08529", "#FF9933", "#FFBB77","#FFDDBB"), 
                       censor=F, risk.table=F, ylab="Cumulative incidence", xlab="Follow-up time (years)", legend.title = "Protein score percentile:", title="Coronary artery disease",
                       legend.labs=c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"),ggtheme = theme_classic())$plot
    df_source_1 <- data.frame(time=fit1$time, n.risk=fit1$cumhaz, fit1$n.censor, outc="Coronary artery disease",
                        strata=c(rep("80-100%", fit1$strata[1]), rep("60-80%", fit1$strata[2]), rep("40-60%", fit1$strata[3]), rep("20-40%", fit1$strata[4]), rep("0-20%", fit1$strata[5])))

    surv_object <- Surv(time = test$hfail_fu, event = test$hfail_inc)
    fit1 <- survfit(surv_object ~ hfail_protriskscore_quintile, data=test)
    survplot_hfail <- ggsurvplot(fit1, pval=F, fun='event', color = "hfail_protriskscore_quintile", ylim=c(0,0.22), legend=c(0.2,0.76), palette=c("#A35C14", "#E08529", "#FF9933", "#FFBB77","#FFDDBB"), 
                       censor=F, risk.table=F, ylab=" ", xlab="Follow-up time (years)", legend.title = "Protein score percentile:", title="Heart failure",
                       legend.labs=c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"),ggtheme = theme_classic())$plot
    df_source_2 <- data.frame(time=fit1$time, n.risk=fit1$cumhaz, fit1$n.censor, outc="Heart failure",
                        strata=c(rep("80-100%", fit1$strata[1]), rep("60-80%", fit1$strata[2]), rep("40-60%", fit1$strata[3]), rep("20-40%", fit1$strata[4]), rep("0-20%", fit1$strata[5])))
    
    surv_object <- Surv(time = test$afib_fu, event = test$afib_inc)
    fit1 <- survfit(surv_object ~ afib_protriskscore_quintile, data=test)
    survplot_afib <- ggsurvplot(fit1, pval=F, fun='event', color = "afib_protriskscore_quintile", ylim=c(0,0.22), legend=c(0.2,0.76), palette=c("#A35C14", "#E08529", "#FF9933", "#FFBB77","#FFDDBB"), 
                       censor=F, risk.table=F, ylab=" ", xlab="Follow-up time (years)", legend.title = "Protein score percentile:", title="Atrial fibrillation",
                       legend.labs=c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"),ggtheme = theme_classic())$plot
    df_source_3 <- data.frame(time=fit1$time, n.risk=fit1$cumhaz, fit1$n.censor, outc="Atrial fibrillation",
                        strata=c(rep("80-100%", fit1$strata[1]), rep("60-80%", fit1$strata[2]), rep("40-60%", fit1$strata[3]), rep("20-40%", fit1$strata[4]), rep("0-20%", fit1$strata[5])))
    
    surv_object <- Surv(time = test$ao_sten_fu, event = test$ao_sten_inc)
    fit1 <- survfit(surv_object ~ ao_sten_protriskscore_quintile, data=test)
    survplot_ao_sten <- ggsurvplot(fit1, pval=F, fun='event', color = "ao_sten_protriskscore_quintile", ylim=c(0,0.22), legend=c(0.2,0.76), palette=c("#A35C14", "#E08529", "#FF9933", "#FFBB77","#FFDDBB"), 
                       censor=F, risk.table=F, ylab=" ", xlab="Follow-up time (years)", legend.title = "Protein score percentile:", title="Aortic stenosis",
                       legend.labs=c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"),ggtheme = theme_classic())$plot
    df_source_4 <- data.frame(time=fit1$time, n.risk=fit1$cumhaz, fit1$n.censor, outc="Aortic stenosis",
                        strata=c(rep("80-100%", fit1$strata[1]), rep("60-80%", fit1$strata[2]), rep("40-60%", fit1$strata[3]), rep("20-40%", fit1$strata[4]), rep("0-20%", fit1$strata[5])))
    
            
    
    tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_km_all.tiff", width=4500, height=1000, res=300)
    ggarrange(survplot_cad, survplot_hfail, survplot_afib, survplot_ao_sten, nrow=1, ncol=4, common.legend=T, legend="bottom")
    dev.off()
    
    fig5_b <- ggarrange(survplot_cad, survplot_hfail, survplot_afib, survplot_ao_sten, nrow=1, ncol=4, common.legend=T, legend="bottom")
    
    write.csv(rbind(df_source_1, df_source_2, df_source_3, df_source_4), ".../ukb_proteomics_cvd/output_files/source_data/fig_5_b.csv", row.names=F)

    
    
        # 8.A.3a - Comparison AUCs ####
    
    df <- fread('.../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz')
      set.seed(1234)
      df <- df[sample(1:nrow(df)), ]
      set.seed(1)
      df$id <- 1:nrow(df)
      train <- df %>% dplyr::sample_frac(0.80)
      test  <- dplyr::anti_join(df, train, by = 'id')
      rm(df)
      
    df_pr_test <- test[,106:1564]
      test <- test[,-(106:1564)]
      test <- cbind(test, df_pr_test)
      rm(df_pr_test)
      
    cad_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_pr.csv")
      hfail_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_pr.csv")
      afib_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_pr.csv")
      ao_sten_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_pr.csv")
    cad_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_cl.csv")
      hfail_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_cl.csv")
      afib_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_cl.csv")
      ao_sten_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_cl.csv")
    cad_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_comb.csv")
      hfail_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_comb.csv")
      afib_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_comb.csv")
      ao_sten_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_comb.csv")
    test$cad_protriskscore <- as.numeric(cbind(1,as.matrix(test[,139:1597])) %*% as.numeric(cad_feat_pr$s1))
      test$hfail_protriskscore <- as.numeric(cbind(1,as.matrix(test[,139:1597])) %*% as.numeric(hfail_feat_pr$s1))
      test$afib_protriskscore <- as.numeric(cbind(1,as.matrix(test[,139:1597])) %*% as.numeric(afib_feat_pr$s1))
      test$ao_sten_protriskscore <- as.numeric(cbind(1,as.matrix(test[,139:1597])) %*% as.numeric(ao_sten_feat_pr$s1))
        test$cad_protriskscore_bin <- exp(test$cad_protriskscore)/(1+exp(test$cad_protriskscore))
      test$hfail_protriskscore_bin <- exp(test$hfail_protriskscore)/(1+exp(test$hfail_protriskscore))
      test$afib_protriskscore_bin <- exp(test$afib_protriskscore)/(1+exp(test$afib_protriskscore))
      test$ao_sten_protriskscore_bin <- exp(test$ao_sten_protriskscore)/(1+exp(test$ao_sten_protriskscore))
        test$cad_protriskscore_quintile <- factor(ntile(test$cad_protriskscore_bin, 5), levels=c(5:1))
      test$hfail_protriskscore_quintile <- factor(ntile(test$hfail_protriskscore_bin, 5), levels=c(5:1))
      test$afib_protriskscore_quintile <- factor(ntile(test$afib_protriskscore_bin, 5), levels=c(5:1))
      test$ao_sten_protriskscore_quintile <- factor(ntile(test$ao_sten_protriskscore_bin, 5), levels=c(5:1))
    test$cad_cl_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 127, 130, 132, 134, 136)])) %*% as.numeric(cad_feat_cl$s1))
      test$hfail_cl_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 127, 130, 132, 134, 136)])) %*% as.numeric(hfail_feat_cl$s1))
      test$afib_cl_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 127, 130, 132, 134, 136)])) %*% as.numeric(afib_feat_cl$s1))
      test$ao_sten_cl_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 127, 130, 132, 134, 136)])) %*% as.numeric(ao_sten_feat_cl$s1))
        test$cad_cl_score_bin <- exp(test$cad_cl_score)/(1+exp(test$cad_cl_score))
      test$hfail_cl_score_bin <- exp(test$hfail_cl_score)/(1+exp(test$hfail_cl_score))
      test$afib_cl_score_bin <- exp(test$afib_cl_score)/(1+exp(test$afib_cl_score))
      test$ao_sten_cl_score_bin <- exp(test$ao_sten_cl_score)/(1+exp(test$ao_sten_cl_score))
        test$cad_cl_score_quintile <- factor(ntile(test$cad_cl_score_bin, 5), levels=c(5:1))
      test$hfail_cl_score_quintile <- factor(ntile(test$hfail_cl_score_bin, 5), levels=c(5:1))
      test$afib_cl_score_quintile <- factor(ntile(test$afib_cl_score_bin, 5), levels=c(5:1))
      test$ao_sten_cl_score_quintile <- factor(ntile(test$ao_sten_cl_score_bin, 5), levels=c(5:1))
    test$cad_comb_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 127, 130, 132, 134, 136, 139:1597)])) %*% as.numeric(cad_feat_comb$s1))
      test$hfail_comb_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 127, 130, 132, 134, 136, 139:1597)])) %*% as.numeric(hfail_feat_comb$s1))
      test$afib_comb_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 127, 130, 132, 134, 136, 139:1597)])) %*% as.numeric(afib_feat_comb$s1))
      test$ao_sten_comb_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 127, 130, 132, 134, 136, 139:1597)])) %*% as.numeric(ao_sten_feat_comb$s1))
        test$cad_comb_score_bin <- exp(test$cad_comb_score)/(1+exp(test$cad_comb_score))
      test$hfail_comb_score_bin <- exp(test$hfail_comb_score)/(1+exp(test$hfail_comb_score))
      test$afib_comb_score_bin <- exp(test$afib_comb_score)/(1+exp(test$afib_comb_score))
      test$ao_sten_comb_score_bin <- exp(test$ao_sten_comb_score)/(1+exp(test$ao_sten_comb_score))
        test$cad_comb_score_quintile <- factor(ntile(test$cad_comb_score_bin, 5), levels=c(5:1))
      test$hfail_comb_score_quintile <- factor(ntile(test$hfail_comb_score_bin, 5), levels=c(5:1))
      test$afib_comb_score_quintile <- factor(ntile(test$afib_comb_score_bin, 5), levels=c(5:1))
      test$ao_sten_comb_score_quintile <- factor(ntile(test$ao_sten_comb_score_bin, 5), levels=c(5:1))
      
    library(pROC)
    roc.test(roc(response=as.numeric(test$cad_inc), predictor=as.numeric(test$cad_cl_score_bin), levels=c(0,1), ci=T),
             roc(response=as.numeric(test$cad_inc), predictor=as.numeric(test$cad_comb_score_bin), levels=c(0,1), ci=T))
    roc.test(roc(response=as.numeric(test$hfail_inc), predictor=as.numeric(test$hfail_cl_score_bin), levels=c(0,1), ci=T),
             roc(response=as.numeric(test$hfail_inc), predictor=as.numeric(test$hfail_comb_score_bin), levels=c(0,1), ci=T))
    roc.test(roc(response=as.numeric(test$afib_inc), predictor=as.numeric(test$afib_cl_score_bin), levels=c(0,1), ci=T),
             roc(response=as.numeric(test$afib_inc), predictor=as.numeric(test$afib_comb_score_bin), levels=c(0,1), ci=T))
    roc.test(roc(response=as.numeric(test$ao_sten_inc), predictor=as.numeric(test$ao_sten_cl_score_bin), levels=c(0,1), ci=T),
             roc(response=as.numeric(test$ao_sten_inc), predictor=as.numeric(test$ao_sten_comb_score_bin), levels=c(0,1), ci=T))


    rm(list=ls()[-c(which(ls()=="test"), which(ls()=="train"))])

    

    
        # 8.A.3b - Evaluation of scores in multivariable-adjusted Cox models ####

    library(survival)   
    library(broom)   
    
    df <- test
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
    
    tidy(coxph(Surv(tstart, tstop, event) ~ Sex_numeric + age + age2 + 
                                        mergedrace + PC1 + PC2 + PC3 + 
                                        PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                        ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                        tchol_final + hdl_final + cholmed + 
                                        dm2_prev + tdi_log_final + creat_final + afib + hfail + ao_sten + cad_protriskscore, data = df_cad), exponentiate = TRUE, conf.int = TRUE)[28,]
    tidy(coxph(Surv(tstart, tstop, event) ~ Sex_numeric + age + age2 + 
                                        mergedrace + PC1 + PC2 + PC3 + 
                                        PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                        ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                        tchol_final + hdl_final + cholmed + 
                                        dm2_prev + tdi_log_final + creat_final + afib + cad + ao_sten + hfail_protriskscore, data = df_hfail), exponentiate = TRUE, conf.int = TRUE)[28,]
    tidy(coxph(Surv(tstart, tstop, event) ~ Sex_numeric + age + age2 + 
                                        mergedrace + PC1 + PC2 + PC3 + 
                                        PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                        ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                        tchol_final + hdl_final + cholmed + 
                                        dm2_prev + tdi_log_final + creat_final + cad + hfail + ao_sten + afib_protriskscore, data = df_afib), exponentiate = TRUE, conf.int = TRUE)[28,]
    tidy(coxph(Surv(tstart, tstop, event) ~ Sex_numeric + age + age2 + 
                                        mergedrace + PC1 + PC2 + PC3 + 
                                        PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
                                        ever_smoked + BMI_final + SBP_final + antihtnbase + 
                                        tchol_final + hdl_final + cholmed + 
                                        dm2_prev + tdi_log_final + creat_final + afib + hfail + cad + ao_sten_protriskscore, data = df_ao_sten), exponentiate = TRUE, conf.int = TRUE)[28,]

      
    
      # 8.A.4 - Quantile plots ####
    
    library(data.table)
    library(dplyr)
    
    df <- fread('.../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz')
      set.seed(1234)
      df <- df[sample(1:nrow(df)), ]
      set.seed(1)
      df$id <- 1:nrow(df)
      train <- df %>% dplyr::sample_frac(0.80)
      test  <- dplyr::anti_join(df, train, by = 'id')
      rm(df)
      
    cad_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_pr.csv")
      hfail_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_pr.csv")
      afib_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_pr.csv")
      ao_sten_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_pr.csv")
    test$cad_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(cad_feat_pr$s1))
      test$hfail_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(hfail_feat_pr$s1))
      test$afib_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(afib_feat_pr$s1))
      test$ao_sten_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(ao_sten_feat_pr$s1))
    test$cad_protriskscore_bin <- exp(test$cad_protriskscore)/(1+exp(test$cad_protriskscore))
      test$hfail_protriskscore_bin <- exp(test$hfail_protriskscore)/(1+exp(test$hfail_protriskscore))
      test$afib_protriskscore_bin <- exp(test$afib_protriskscore)/(1+exp(test$afib_protriskscore))
      test$ao_sten_protriskscore_bin <- exp(test$ao_sten_protriskscore)/(1+exp(test$ao_sten_protriskscore))
    test$cad_protriskscore_quintile <- factor(ntile(test$cad_protriskscore_bin, 5), levels=c(5:1))
      test$hfail_protriskscore_quintile <- factor(ntile(test$hfail_protriskscore_bin, 5), levels=c(5:1))
      test$afib_protriskscore_quintile <- factor(ntile(test$afib_protriskscore_bin, 5), levels=c(5:1))
      test$ao_sten_protriskscore_quintile <- factor(ntile(test$ao_sten_protriskscore_bin, 5), levels=c(5:1))
    test$cad_protriskscore_quantile <- factor(ntile(test$cad_protriskscore_bin, 10), levels=c(10:1))
      test$hfail_protriskscore_quantile <- factor(ntile(test$hfail_protriskscore_bin, 10), levels=c(10:1))
      test$afib_protriskscore_quantile <- factor(ntile(test$afib_protriskscore_bin, 10), levels=c(10:1))
      test$ao_sten_protriskscore_quantile <- factor(ntile(test$ao_sten_protriskscore_bin, 10), levels=c(10:1))
    
    library(fmsb)
    library(ggplot2)
      
    df_cad <- data.frame(quant=NA, inc_r=NA, inc_r_lci=NA, inc_r_uci=NA)[-1,]
    for (i in 1:10) {
      df_tmp <- test[test$cad_protriskscore_quantile==i,]
      n <- nrow(df_tmp)
      inc <- sum(df_tmp$cad_inc)
      fu <- sum(df_tmp$cad_fu)
      inc_r <- IRCI(inc, fu, conf.level=0.95)$IR
      inc_r_lci <- IRCI(inc, fu, conf.level=0.95)$IRL
      inc_r_uci <- IRCI(inc, fu, conf.level=0.95)$IRU
      row <- data.frame(quant=i, inc_r=inc_r, inc_r_lci=inc_r_lci, inc_r_uci=inc_r_uci)
      df_cad <- rbind(df_cad, row)
    }
    rm(n, inc, fu, inc_r, inc_r_lci, inc_r_uci, row)
    cad_plot <- ggplot(data=df_cad, aes(x=quant*10-5, y=1000*inc_r, color=quant*10-5)) +
                  geom_point() +
                  scale_color_gradient2(low="#FFDDBB", mid="#FF9933", high="#A35C14", midpoint = 50, limits=c(0,100))  + 
                  scale_x_continuous(limits=c(0,100), name="Protein score percentile") +
                  scale_y_continuous(name="Incidence rate per 1,000 years", limits = c(0.05, 32), trans="log2") +
                  ggtitle("Coronary artery disease") +
                  theme_classic() +
                  theme(legend.position = "none")
      
    df_hfail <- data.frame(quant=NA, inc_r=NA, inc_r_lci=NA, inc_r_uci=NA)[-1,]
    for (i in 1:10) {
      df_tmp <- test[test$hfail_protriskscore_quantile==i,]
      n <- nrow(df_tmp)
      inc <- sum(df_tmp$hfail_inc)
      fu <- sum(df_tmp$hfail_fu)
      inc_r <- IRCI(inc, fu, conf.level=0.95)$IR
      inc_r_lci <- IRCI(inc, fu, conf.level=0.95)$IRL
      inc_r_uci <- IRCI(inc, fu, conf.level=0.95)$IRU
      row <- data.frame(quant=i, inc_r=inc_r, inc_r_lci=inc_r_lci, inc_r_uci=inc_r_uci)
      df_hfail <- rbind(df_hfail, row)
    }
    rm(n, inc, fu, inc_r, inc_r_lci, inc_r_uci, row)
    hfail_plot <- ggplot(data=df_hfail, aes(x=quant*10-5, y=1000*inc_r, color=quant*10-5)) +
                  geom_point() +
                  scale_color_gradient2(low="#FFDDBB", mid="#FF9933", high="#A35C14", midpoint = 50, limits=c(0,100))  + 
                  scale_x_continuous(limits=c(0,100), name="Protein score percentile") +
                  scale_y_continuous(name="", limits = c(0.05, 32), trans="log2") +
                  ggtitle("Heart failure") +
                  theme_classic() +
                  theme(legend.position = "none")
      
   df_afib <- data.frame(quant=NA, inc_r=NA, inc_r_lci=NA, inc_r_uci=NA)[-1,]
    for (i in 1:10) {
      df_tmp <- test[test$afib_protriskscore_quantile==i,]
      n <- nrow(df_tmp)
      inc <- sum(df_tmp$afib_inc)
      fu <- sum(df_tmp$afib_fu)
      inc_r <- IRCI(inc, fu, conf.level=0.95)$IR
      inc_r_lci <- IRCI(inc, fu, conf.level=0.95)$IRL
      inc_r_uci <- IRCI(inc, fu, conf.level=0.95)$IRU
      row <- data.frame(quant=i, inc_r=inc_r, inc_r_lci=inc_r_lci, inc_r_uci=inc_r_uci)
      df_afib <- rbind(df_afib, row)
    }
    rm(n, inc, fu, inc_r, inc_r_lci, inc_r_uci, row)
    afib_plot <- ggplot(data=df_afib, aes(x=quant*10-5, y=1000*inc_r, color=quant*10-5)) +
                  geom_point() +
                  scale_color_gradient2(low="#FFDDBB", mid="#FF9933", high="#A35C14", midpoint = 50, limits=c(0,100))  + 
                  scale_x_continuous(limits=c(0,100), name="Protein score percentile") +
                  scale_y_continuous(name="", limits = c(0.05, 32), trans="log2") +
                  ggtitle("Atrial fibrillation") +
                  theme_classic() +
                  theme(legend.position = "none")
      
   df_ao_sten <- data.frame(quant=NA, inc_r=NA, inc_r_lci=NA, inc_r_uci=NA)[-1,]
    for (i in 1:10) {
      df_tmp <- test[test$ao_sten_protriskscore_quantile==i,]
      n <- nrow(df_tmp)
      inc <- sum(df_tmp$ao_sten_inc)
      fu <- sum(df_tmp$ao_sten_fu)
      inc_r <- IRCI(inc, fu, conf.level=0.95)$IR
      inc_r_lci <- IRCI(inc, fu, conf.level=0.95)$IRL
      inc_r_uci <- IRCI(inc, fu, conf.level=0.95)$IRU
      row <- data.frame(quant=i, inc_r=inc_r, inc_r_lci=inc_r_lci, inc_r_uci=inc_r_uci)
      df_ao_sten <- rbind(df_ao_sten, row)
    }
    rm(n, inc, fu, inc_r, inc_r_lci, inc_r_uci, row)
    df_ao_sten <- df_ao_sten[df_ao_sten$inc_r>0,]
    ao_sten_plot <- ggplot(data=df_ao_sten, aes(x=quant*10-5, y=1000*inc_r, color=quant*10-5)) +
                  geom_point() +
                  scale_color_gradient2(low="#FFDDBB", mid="#FF9933", high="#A35C14", midpoint = 50, limits=c(0,100))  + 
                  scale_x_continuous(limits=c(0,100), name="Protein score percentile") +
                  scale_y_continuous(name="", limits = c(0.05, 32), trans="log2") +
                  ggtitle("Aortic stenosis") +
                  theme_classic() +
                  theme(legend.position = "none")

    
    library(ggpubr)
    
    tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_quantileplots_all.tiff", width=4500, height=1000, res=300)
    ggarrange(cad_plot, hfail_plot, afib_plot, ao_sten_plot, nrow=1, ncol=4)
    dev.off()
    
    fig5_c <- ggarrange(cad_plot, hfail_plot, afib_plot, ao_sten_plot, nrow=1, ncol=4)
    
      df_cad$outc <- "Coronary artery disease"
      df_hfail$outc <- "Heart failure"
      df_afib$outc <- "Atrial fibrillation"
      df_ao_sten$outc <- "Aortic stenosis"
        
      write.csv(rbind(df_cad, df_hfail, df_afib, df_ao_sten), ".../ukb_proteomics_cvd/output_files/source_data/fig_5_c.csv", row.names=F)

    
      # 8.A.5 - Violin plots and FPR-TPR tables ####
    
    library(data.table)
    library(dplyr)
    library(pROC)
    
    df <- fread('.../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz')
      set.seed(1234)
      df <- df[sample(1:nrow(df)), ]
      set.seed(1)
      df$id <- 1:nrow(df)
      train <- df %>% dplyr::sample_frac(0.80)
      test  <- dplyr::anti_join(df, train, by = 'id')
      rm(df)
      
    cad_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_pr.csv")
      hfail_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_pr.csv")
      afib_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_pr.csv")
      ao_sten_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_pr.csv")
        cad_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_cl.csv")
      hfail_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_cl.csv")
      afib_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_cl.csv")
      ao_sten_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_cl.csv")
        cad_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_comb.csv")
      hfail_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_comb.csv")
      afib_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_comb.csv")
      ao_sten_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_comb.csv")
        cad_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_pr.csv")
      hfail_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_pr.csv")
      afib_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_pr.csv")
      ao_sten_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_pr.csv")
    test$cad_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(cad_feat_pr$s1))
      test$hfail_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(hfail_feat_pr$s1))
      test$afib_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(afib_feat_pr$s1))
      test$ao_sten_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(ao_sten_feat_pr$s1))
        test$cad_protriskscore_bin <- exp(test$cad_protriskscore)/(1+exp(test$cad_protriskscore))
      test$hfail_protriskscore_bin <- exp(test$hfail_protriskscore)/(1+exp(test$hfail_protriskscore))
      test$afib_protriskscore_bin <- exp(test$afib_protriskscore)/(1+exp(test$afib_protriskscore))
      test$ao_sten_protriskscore_bin <- exp(test$ao_sten_protriskscore)/(1+exp(test$ao_sten_protriskscore))
    test$cad_cl_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595)])) %*% as.numeric(cad_feat_cl$s1))
      test$hfail_cl_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595)])) %*% as.numeric(hfail_feat_cl$s1))
      test$afib_cl_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595)])) %*% as.numeric(afib_feat_cl$s1))
      test$ao_sten_cl_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595)])) %*% as.numeric(ao_sten_feat_cl$s1))
        test$cad_cl_score_bin <- exp(test$cad_cl_score)/(1+exp(test$cad_cl_score))
      test$hfail_cl_score_bin <- exp(test$hfail_cl_score)/(1+exp(test$hfail_cl_score))
      test$afib_cl_score_bin <- exp(test$afib_cl_score)/(1+exp(test$afib_cl_score))
      test$ao_sten_cl_score_bin <- exp(test$ao_sten_cl_score)/(1+exp(test$ao_sten_cl_score))
        test$cad_cl_score_quintile <- factor(ntile(test$cad_cl_score_bin, 5), levels=c(5:1))
      test$hfail_cl_score_quintile <- factor(ntile(test$hfail_cl_score_bin, 5), levels=c(5:1))
      test$afib_cl_score_quintile <- factor(ntile(test$afib_cl_score_bin, 5), levels=c(5:1))
      test$ao_sten_cl_score_quintile <- factor(ntile(test$ao_sten_cl_score_bin, 5), levels=c(5:1))
    test$cad_comb_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595, 106:1564)])) %*% as.numeric(cad_feat_comb$s1))
      test$hfail_comb_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595, 106:1564)])) %*% as.numeric(hfail_feat_comb$s1))
      test$afib_comb_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595, 106:1564)])) %*% as.numeric(afib_feat_comb$s1))
      test$ao_sten_comb_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595, 106:1564)])) %*% as.numeric(ao_sten_feat_comb$s1))
        test$cad_comb_score_bin <- exp(test$cad_comb_score)/(1+exp(test$cad_comb_score))
      test$hfail_comb_score_bin <- exp(test$hfail_comb_score)/(1+exp(test$hfail_comb_score))
      test$afib_comb_score_bin <- exp(test$afib_comb_score)/(1+exp(test$afib_comb_score))
      test$ao_sten_comb_score_bin <- exp(test$ao_sten_comb_score)/(1+exp(test$ao_sten_comb_score))
        test$cad_comb_score_quintile <- factor(ntile(test$cad_comb_score_bin, 5), levels=c(5:1))
      test$hfail_comb_score_quintile <- factor(ntile(test$hfail_comb_score_bin, 5), levels=c(5:1))
      test$afib_comb_score_quintile <- factor(ntile(test$afib_comb_score_bin, 5), levels=c(5:1))
      test$ao_sten_comb_score_quintile <- factor(ntile(test$ao_sten_comb_score_bin, 5), levels=c(5:1))


    plot_cad <- ggplot(test, aes(x=as.factor(cad_inc), y=scale(cad_protriskscore), fill=as.factor(cad_inc))) +
      geom_violin(trim=F) +
      geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
      scale_fill_manual(values=c("#FFDDBB","#A35C14"), labels=c("Controls", "Cases"), name="") +
      scale_x_discrete(labels=c("Controls", "Cases")) +
      scale_y_continuous(limits=c(-5,8)) +
      theme_classic() +
      labs(title="Coronary artery disease", x="Group", y = "Protein score")
    plot_hfail <- ggplot(test, aes(x=as.factor(hfail_inc), y=scale(hfail_protriskscore), fill=as.factor(hfail_inc))) +
      geom_violin(trim=F) +
      geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
      scale_fill_manual(values=c("#FFDDBB","#A35C14"), labels=c("Controls", "Cases"), name="") +
      scale_x_discrete(labels=c("Controls", "Cases")) +
      scale_y_continuous(limits=c(-5,8)) +
      theme_classic() +
      labs(title="Heart failure", x="Group", y = "")
    plot_afib <- ggplot(test, aes(x=as.factor(afib_inc), y=scale(afib_protriskscore), fill=as.factor(afib_inc))) +
      geom_violin(trim=F) +
      geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
      scale_fill_manual(values=c("#FFDDBB","#A35C14"), labels=c("Controls", "Cases"), name="") +
      scale_x_discrete(labels=c("Controls", "Cases")) +
      scale_y_continuous(limits=c(-5,8)) +
      theme_classic() +
      labs(title="Atrial fibrillation", x="Group", y = "")
    plot_ao_sten <- ggplot(test, aes(x=as.factor(ao_sten_inc), y=scale(ao_sten_protriskscore), fill=as.factor(ao_sten_inc))) +
      geom_violin(trim=F) +
      geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
      scale_fill_manual(values=c("#FFDDBB","#A35C14"), labels=c("Controls", "Cases"), name="") +
      scale_x_discrete(labels=c("Controls", "Cases")) +
      scale_y_continuous(limits=c(-5,8)) +
      theme_classic() +
      labs(title="Aortic stenosis", x="Group", y = "")

      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_violins.tiff", width=5300, height=1500, res=300)
        ggarrange(plot_cad, plot_hfail, plot_afib, plot_ao_sten, common.legend=T, ncol=4, nrow=1, legend="bottom")
        dev.off()
       
        
        
        
        
        
    plot_cad <- ggplot(test, aes(x=scale(cad_protriskscore), fill=as.factor(cad_inc))) +
      geom_density(alpha=0.4) +
      scale_fill_manual(values=c("#FFDDBB","#A35C14"), labels=c("Controls", "Cases"), name="") +
      geom_vline(xintercept = quantile(scale(test$cad_protriskscore[test$cad_inc==0]), 0.95) ) +
      theme_classic() +
      labs(title="Coronary artery disease", x="Protein score (SD)", y = "Density")
    plot_hfail <- ggplot(test, aes(x=scale(hfail_protriskscore), fill=as.factor(hfail_inc))) +
      geom_density(alpha=0.4) +
      scale_fill_manual(values=c("#FFDDBB","#A35C14"), labels=c("Controls", "Cases"), name="") +
      geom_vline(xintercept = quantile(scale(test$hfail_protriskscore[test$hfail_inc==0]), 0.95) ) +
      theme_classic() +
      labs(title="Heart failure", x = "Protein score (SD)", y=" ")
    plot_afib <- ggplot(test, aes(x=scale(afib_protriskscore), fill=as.factor(afib_inc))) +
      geom_density(alpha=0.4) +
      scale_fill_manual(values=c("#FFDDBB","#A35C14"), labels=c("Controls", "Cases"), name="") +
      geom_vline(xintercept = quantile(scale(test$afib_protriskscore[test$afib_inc==0]), 0.95) ) +
      theme_classic() +
      labs(title="Atrial fibrillation", x = "Protein score (SD)", y=" ")
    plot_ao_sten <- ggplot(test, aes(x=scale(ao_sten_protriskscore), fill=as.factor(ao_sten_inc))) +
      geom_density(alpha=0.4) +
      scale_fill_manual(values=c("#FFDDBB","#A35C14"), labels=c("Controls", "Cases"), name="") +
      geom_vline(xintercept = quantile(scale(test$ao_sten_protriskscore[test$ao_sten_inc==0]), 0.95) ) +
      theme_classic() +
      labs(title="Aortic stenosis", x = "Protein score (SD)", y=" ")

      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_distributions.tiff", width=4500, height=1000, res=300)
        ggarrange(plot_cad, plot_hfail, plot_afib, plot_ao_sten, common.legend=T, ncol=4, nrow=1, legend="bottom")
        dev.off()

    fig5_a <- ggarrange(plot_cad, plot_hfail, plot_afib, plot_ao_sten, common.legend=T, ncol=4, nrow=1, legend="bottom")
    
    cad <- ggplot_build(plot_cad)$data[[1]]
      cad$outc <- "Coronary artery disease"
    hfail <- ggplot_build(plot_hfail)$data[[1]]
      hfail$outc <- "Heart failure"
    afib <- ggplot_build(plot_afib)$data[[1]]
      afib$outc <- "Heart failure"
    ao_sten <- ggplot_build(plot_ao_sten)$data[[1]]
      ao_sten$outc <- "Aortic stenosis"
    
    df_source <- rbind(cad, hfail, afib, ao_sten)
    write.csv(df_source[,c("fill", "y", "x", "outc")], ".../ukb_proteomics_cvd/output_files/source_data/fig_5_d.csv", row.names=F)
    
    
    Comp_1_dml <- dml(code = print(ggarrange(fig5_a, fig5_b, fig5_c, fig5_d, nrow=4, ncol=1, heights = c(1000, 1000, 1000, 1300)), newpage = FALSE))
    
    read_pptx("C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_1_dml, ph_location(left = 1, top = 1, width = 4500/300, height = 4300/300)) %>%
      print(target = "C:/Users/artsc/Documents/- Onderzoek/MCH/13 - Proteomics CVD UKB/Figures.pptx")
    


    test$cad_protriskscore_norm <- scale(test$cad_protriskscore)      
    test$hfail_protriskscore_norm <- scale(test$hfail_protriskscore)      
    test$afib_protriskscore_norm <- scale(test$afib_protriskscore)      
    test$ao_sten_protriskscore_norm <- scale(test$ao_sten_protriskscore)      
      
    quantile(test$cad_protriskscore_norm[test$cad_inc==0], 0.95) 
      percentile <- ecdf(test$cad_protriskscore_norm[test$cad_inc==1])
      1 - percentile(quantile(test$cad_protriskscore_norm[test$cad_inc==0], 0.95))
    quantile(test$hfail_protriskscore_norm[test$hfail_inc==0], 0.95)
      percentile <- ecdf(test$hfail_protriskscore_norm[test$hfail_inc==1])
      1 - percentile(quantile(test$hfail_protriskscore_norm[test$hfail_inc==0], 0.95))
    quantile(test$afib_protriskscore_norm[test$afib_inc==0], 0.95) 
      percentile <- ecdf(test$afib_protriskscore_norm[test$afib_inc==1])
      1 - percentile(quantile(test$afib_protriskscore_norm[test$afib_inc==0], 0.95))
    quantile(test$ao_sten_protriskscore_norm[test$ao_sten_inc==0], 0.95) 
      percentile <- ecdf(test$ao_sten_protriskscore_norm[test$ao_sten_inc==1])
      1 - percentile(quantile(test$ao_sten_protriskscore_norm[test$ao_sten_inc==0], 0.95))

        
        
        
        
    scores <- c("_cl_score_bin", "_protriskscore_bin", "_comb_score_bin")
    outcomes <- c("afib", "cad", "hfail", "ao_sten")
    df_res <- data.frame(Protein=NA, Outcome=NA, FPR=NA, TPR=NA)[-1,]
    
    for (j in outcomes) {
      for (i in scores) {
        for (k in c(0.01, 0.05, 0.1)) {
      roc_df <- data.frame(TPR=roc(response=as.numeric(test[[paste0(j, "_inc")]]), predictor=as.numeric(test[[paste0(j, i)]]), levels=c(0,1), ci=T)$sensitivities, 
                                FPR=1-roc(response=as.numeric(test[[paste0(j, "_inc")]]), predictor=as.numeric(test[[paste0(j, i)]]), levels=c(0,1), ci=T)$specificities)
      roc_df <- roc_df[roc_df$FPR < k,]
      df_res_row <- data.frame(Protein=i, Outcome=j, FPR=k, TPR=roc_df$TPR[1])
      df_res <- rbind(df_res, df_res_row)
       }
      }
    }

    df_res$LR <- df_res$TPR / df_res$FPR
    df_res$inc <- ifelse(df_res$Outcome == "afib", (sum(test$afib_inc)+sum(train$afib_inc)) / (length(test$afib_inc)+length(train$afib_inc)), 
                        ifelse(df_res$Outcome == "cad", (sum(test$cad_inc)+sum(train$cad_inc)) / (length(test$cad_inc)+length(train$cad_inc)), 
                        ifelse(df_res$Outcome == "hfail", (sum(test$hfail_inc)+sum(train$hfail_inc)) / (length(test$hfail_inc)+length(train$hfail_inc)), 
                         (sum(test$ao_sten_inc)+sum(train$ao_sten_inc)) / (length(test$ao_sten_inc)+length(train$ao_sten_inc)))))
    df_res$odds <- df_res$inc / (1-df_res$inc)
    df_res$odds_ifpositive <- df_res$odds * df_res$LR
    df_res$prob_ifpositive <- df_res$odds_ifpositive / (1+df_res$odds_ifpositive)
    # write.csv(df_res, ".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_riskscores_probabilities_odds.csv", row.names=F)
    
        
      # 8.A.6 - DR / FPR (approximated) ####
  
    library(data.table)
    library(dplyr)
    library(pROC)
    
    df <- fread('.../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz')
      set.seed(1234)
      df <- df[sample(1:nrow(df)), ]
      set.seed(1)
      df$id <- 1:nrow(df)
      train <- df %>% dplyr::sample_frac(0.80)
      test  <- dplyr::anti_join(df, train, by = 'id')
      rm(df)
      
    cad_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_pr.csv")
      hfail_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_pr.csv")
      afib_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_pr.csv")
      ao_sten_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_pr.csv")
        cad_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_cl.csv")
      hfail_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_cl.csv")
      afib_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_cl.csv")
      ao_sten_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_cl.csv")
        cad_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_comb.csv")
      hfail_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_comb.csv")
      afib_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_comb.csv")
      ao_sten_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_comb.csv")
        cad_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_pr.csv")
      hfail_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_pr.csv")
      afib_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_pr.csv")
      ao_sten_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_pr.csv")
    test$cad_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(cad_feat_pr$s1))
      test$hfail_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(hfail_feat_pr$s1))
      test$afib_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(afib_feat_pr$s1))
      test$ao_sten_protriskscore <- as.numeric(cbind(1,as.matrix(test[,106:1564])) %*% as.numeric(ao_sten_feat_pr$s1))
        test$cad_protriskscore_bin <- exp(test$cad_protriskscore)/(1+exp(test$cad_protriskscore))
      test$hfail_protriskscore_bin <- exp(test$hfail_protriskscore)/(1+exp(test$hfail_protriskscore))
      test$afib_protriskscore_bin <- exp(test$afib_protriskscore)/(1+exp(test$afib_protriskscore))
      test$ao_sten_protriskscore_bin <- exp(test$ao_sten_protriskscore)/(1+exp(test$ao_sten_protriskscore))
    test$cad_cl_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595)])) %*% as.numeric(cad_feat_cl$s1))
      test$hfail_cl_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595)])) %*% as.numeric(hfail_feat_cl$s1))
      test$afib_cl_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595)])) %*% as.numeric(afib_feat_cl$s1))
      test$ao_sten_cl_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595)])) %*% as.numeric(ao_sten_feat_cl$s1))
        test$cad_cl_score_bin <- exp(test$cad_cl_score)/(1+exp(test$cad_cl_score))
      test$hfail_cl_score_bin <- exp(test$hfail_cl_score)/(1+exp(test$hfail_cl_score))
      test$afib_cl_score_bin <- exp(test$afib_cl_score)/(1+exp(test$afib_cl_score))
      test$ao_sten_cl_score_bin <- exp(test$ao_sten_cl_score)/(1+exp(test$ao_sten_cl_score))
        test$cad_cl_score_quintile <- factor(ntile(test$cad_cl_score_bin, 5), levels=c(5:1))
      test$hfail_cl_score_quintile <- factor(ntile(test$hfail_cl_score_bin, 5), levels=c(5:1))
      test$afib_cl_score_quintile <- factor(ntile(test$afib_cl_score_bin, 5), levels=c(5:1))
      test$ao_sten_cl_score_quintile <- factor(ntile(test$ao_sten_cl_score_bin, 5), levels=c(5:1))
    test$cad_comb_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595, 106:1564)])) %*% as.numeric(cad_feat_comb$s1))
      test$hfail_comb_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595, 106:1564)])) %*% as.numeric(hfail_feat_comb$s1))
      test$afib_comb_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595, 106:1564)])) %*% as.numeric(afib_feat_comb$s1))
      test$ao_sten_comb_score <- as.numeric(cbind(1,as.matrix(test[,c(3, 4, 7, 39, 41, 44, 53, 1586, 1589, 1591, 1593, 1595, 106:1564)])) %*% as.numeric(ao_sten_feat_comb$s1))
        test$cad_comb_score_bin <- exp(test$cad_comb_score)/(1+exp(test$cad_comb_score))
      test$hfail_comb_score_bin <- exp(test$hfail_comb_score)/(1+exp(test$hfail_comb_score))
      test$afib_comb_score_bin <- exp(test$afib_comb_score)/(1+exp(test$afib_comb_score))
      test$ao_sten_comb_score_bin <- exp(test$ao_sten_comb_score)/(1+exp(test$ao_sten_comb_score))
        test$cad_comb_score_quintile <- factor(ntile(test$cad_comb_score_bin, 5), levels=c(5:1))
      test$hfail_comb_score_quintile <- factor(ntile(test$hfail_comb_score_bin, 5), levels=c(5:1))
      test$afib_comb_score_quintile <- factor(ntile(test$afib_comb_score_bin, 5), levels=c(5:1))
      test$ao_sten_comb_score_quintile <- factor(ntile(test$ao_sten_comb_score_bin, 5), levels=c(5:1))

    rm(afib_feat_cl, afib_feat_comb, afib_feat_pr, cad_feat_cl, cad_feat_comb, cad_feat_pr, 
       hfail_feat_cl, hfail_feat_comb, hfail_feat_pr, ao_sten_feat_cl, ao_sten_feat_comb, ao_sten_feat_pr)
      
    calculate_detection_rate <- function(mean_cases, mean_controls, sd_cases, sd_controls, FPR) {
      z_controls <- qnorm(1 - FPR)
      z_shift <- (mean_cases - mean_controls) / sd_controls
      z_affected <- z_controls - z_shift
      TPR_2 <- pnorm(z_affected, mean = 0, sd = sd_cases / sd_controls)
      TPR <- 1 - TPR_2
      return(TPR)
    }

    df_res <- data.frame(Protein=NA, Outcome=NA, Cases_N=NA, Cases_Mean=NA, Cases_SD=NA, Controls_N=NA, Controls_Mean=NA, Controls_SD=NA)[-1,]
    for (j in c("cad", "afib", "hfail", "ao_sten")) {
      for (i in c("protriskscore", "cl_score", "comb_score")) {
      
      cases_n <- sum(test[[paste0(j, "_inc")]]==1)
      cases_mean <- mean(test[[paste0(j, "_", i)]][test[[paste0(j, "_inc")]]==1])
      cases_sd <- sd(test[[paste0(j, "_", i)]][test[[paste0(j, "_inc")]]==1])
      controls_n <- sum(test[[paste0(j, "_inc")]]==0)
      controls_mean <- mean(test[[paste0(j, "_", i)]][test[[paste0(j, "_inc")]]==0])
      controls_sd <- sd(test[[paste0(j, "_", i)]][test[[paste0(j, "_inc")]]==0])
      dr1 <- calculate_detection_rate(cases_mean, controls_mean, cases_sd, controls_sd, 0.01)
      dr5 <- calculate_detection_rate(cases_mean, controls_mean, cases_sd, controls_sd, 0.05)
      dr10 <- calculate_detection_rate(cases_mean, controls_mean, cases_sd, controls_sd, 0.1)
      
      df_res_row <- data.frame(Protein=i, Outcome=j, Cases_N=cases_n, Cases_Mean=cases_mean, Cases_SD=cases_sd, 
                               Controls_N=controls_n, Controls_Mean=controls_mean, Controls_SD=controls_sd, 
                               DR_1=dr1, DR_5=dr5, DR_10=dr10)
      df_res <- rbind(df_res, df_res_row)
      
      }
    }
    
    write.csv(df_res, ".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_fpr_tpr_models_approximated.csv", row.names=F)
    rm(df_res, df_res_row, test, train, cases_mean, cases_n, cases_sd, controls_mean, controls_n, controls_sd, dr1, dr5, dr10, i, j, calculate_detection_rate())
    

    # 8.B - UKB - overlap WHI ####
        # 8.B.1a - CAD ####
        
      cad_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_formula_cad_cl.csv")
          cad_feat_cl <- cad_feat_cl[cad_feat_cl$s1!=0 & cad_feat_cl$V1!="(Intercept)",]
          cad_feat_cl$sd <- c(1, 8.2, 1, 1, 1, 1, 1, 4.7, 19.6, 43.8, 14.9, 0.19)
          cad_feat_cl$s1_weighted <- abs(cad_feat_cl$s1*cad_feat_cl$sd)
          cad_feat_cl$name_form <- c("Sex (male vs. female)", "Age (per SD)", "Smoking (ever vs. never)", "Cholesterol-lowering medication use (yes vs. no)", 
                                     "Race/ethnicity (non-White vs. White)", "Antihypertensive medication use (yes vs. no)", "Type 2 diabetes mellitus (yes vs. no)", 
                                     "BMI (per SD)", "Systolic blood pressure (per SD)", "Total cholesterol (per SD)", "HDL cholesterol (per SD)", "Creatinine (per SD)")
          cad_feat_cl <- cad_feat_cl[order(cad_feat_cl$s1_weighted),]
          cad_feat_cl$name_form <- sapply( strwrap(cad_feat_cl$name_form, 22, simplify=FALSE), paste, collapse="\n" )
      cad_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_formula_cad_pr.csv")
          cad_feat_pr <- cad_feat_pr[cad_feat_pr$s1!=0 & cad_feat_pr$V1!="(Intercept)",]
          cad_feat_pr$s1_abs <- abs(cad_feat_pr$s1)
          cad_feat_pr <- cad_feat_pr[order(cad_feat_pr$s1_abs),]
          
      cadplot_cl <- ggplot(data=cad_feat_cl, aes(x=fct_reorder(name_form, s1_weighted), y=s1_weighted)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("CAD") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
      cadplot_pr <- ggplot(data=cad_feat_pr, aes(x=fct_reorder(V1, s1_abs), y=s1_abs)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("CAD") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),axis.text.y = element_text(size = 5))
      rm(cad_feat_cl, cad_feat_pr)
       
        
        # 8.B.1b - HF ####
        
      hfail_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_formula_hfail_cl.csv")
          hfail_feat_cl <- hfail_feat_cl[hfail_feat_cl$s1!=0 & hfail_feat_cl$V1!="(Intercept)",]
          hfail_feat_cl$sd <- c(1, 8.2, 1, 1, 1, 1, 4.7, 19.6, 43.8, 0.19)
          hfail_feat_cl$s1_weighted <- abs(hfail_feat_cl$s1*hfail_feat_cl$sd)
          hfail_feat_cl$name_form <- c("Sex (male vs. female)", "Age (per SD)", "Smoking (ever vs. never)", "Cholesterol-lowering medication use (yes vs. no)", 
                                     "Antihypertensive medication use (yes vs. no)", "Type 2 diabetes mellitus (yes vs. no)", 
                                     "BMI (per SD)", "Systolic blood pressure (per SD)", "Total cholesterol (per SD)", "Creatinine (per SD)")
          hfail_feat_cl <- hfail_feat_cl[order(hfail_feat_cl$s1_weighted),]
          hfail_feat_cl$name_form <- sapply( strwrap(hfail_feat_cl$name_form, 22, simplify=FALSE), paste, collapse="\n" )
      hfail_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_formula_hfail_pr.csv")
          hfail_feat_pr <- hfail_feat_pr[hfail_feat_pr$s1!=0 & hfail_feat_pr$V1!="(Intercept)",]
          hfail_feat_pr$s1_abs <- abs(hfail_feat_pr$s1)
          hfail_feat_pr <- hfail_feat_pr[order(hfail_feat_pr$s1_abs),]
          
      hfailplot_cl <- ggplot(data=hfail_feat_cl, aes(x=fct_reorder(name_form, s1_weighted), y=s1_weighted)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("HF") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
      hfailplot_pr <- ggplot(data=hfail_feat_pr, aes(x=fct_reorder(V1, s1_abs), y=s1_abs)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("HF") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),axis.text.y = element_text(size = 5))
      rm(hfail_feat_cl, hfail_feat_pr)
                      
        
        # 8.B.1c - AF ####
        
      afib_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_formula_afib_cl.csv")
          afib_feat_cl <- afib_feat_cl[afib_feat_cl$s1!=0 & afib_feat_cl$V1!="(Intercept)",]
          afib_feat_cl$sd <- c(1, 8.2, 1, 1, 1, 1, 1, 4.7, 19.6, 43.8, 14.9, 0.19)
          afib_feat_cl$s1_weighted <- abs(afib_feat_cl$s1*afib_feat_cl$sd)
          afib_feat_cl$name_form <- c("Sex (male vs. female)", "Age (per SD)", "Smoking (ever vs. never)", "Cholesterol-lowering medication use (yes vs. no)", 
                                     "Race/ethnicity (non-White vs. White)", "Antihypertensive medication use (yes vs. no)", "Type 2 diabetes mellitus (yes vs. no)", 
                                     "BMI (per SD)", "Systolic blood pressure (per SD)", "Total cholesterol (per SD)", "HDL cholesterol (per SD)", "Creatinine (per SD)")
          afib_feat_cl <- afib_feat_cl[order(afib_feat_cl$s1_weighted),]
          afib_feat_cl$name_form <- sapply( strwrap(afib_feat_cl$name_form, 22, simplify=FALSE), paste, collapse="\n" )
      afib_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_formula_afib_pr.csv")
          afib_feat_pr <- afib_feat_pr[afib_feat_pr$s1!=0 & afib_feat_pr$V1!="(Intercept)",]
          afib_feat_pr$s1_abs <- abs(afib_feat_pr$s1)
          afib_feat_pr <- afib_feat_pr[order(afib_feat_pr$s1_abs),]
          
      afibplot_cl <- ggplot(data=afib_feat_cl, aes(x=fct_reorder(name_form, s1_weighted), y=s1_weighted)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("AF") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
      afibplot_pr <- ggplot(data=afib_feat_pr, aes(x=fct_reorder(V1, s1_abs), y=s1_abs)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("AF") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),axis.text.y = element_text(size = 5))
      rm(afib_feat_cl, afib_feat_pr)
       
      
        # 8.B.1d - AS ####
        
      ao_sten_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_formula_ao_sten_cl.csv")
          ao_sten_feat_cl <- ao_sten_feat_cl[ao_sten_feat_cl$s1!=0 & ao_sten_feat_cl$V1!="(Intercept)",]
          ao_sten_feat_cl$sd <- c(1, 8.2, 1, 1, 1, 1, 4.7, 19.6, 43.8, 14.9, 0.19)
          ao_sten_feat_cl$s1_weighted <- abs(ao_sten_feat_cl$s1*ao_sten_feat_cl$sd)
          ao_sten_feat_cl$name_form <- c("Sex (male vs. female)", "Age (per SD)", "Smoking (ever vs. never)", "Cholesterol-lowering medication use (yes vs. no)", 
                                     "Race/ethnicity (non-White vs. White)", "Antihypertensive medication use (yes vs. no)", 
                                     "BMI (per SD)", "Systolic blood pressure (per SD)", "Total cholesterol (per SD)", "HDL cholesterol (per SD)", "Creatinine (per SD)")
          ao_sten_feat_cl <- ao_sten_feat_cl[order(ao_sten_feat_cl$s1_weighted),]
          ao_sten_feat_cl$name_form <- sapply( strwrap(ao_sten_feat_cl$name_form, 22, simplify=FALSE), paste, collapse="\n" )
      ao_sten_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_formula_ao_sten_pr.csv")
          ao_sten_feat_pr <- ao_sten_feat_pr[ao_sten_feat_pr$s1!=0 & ao_sten_feat_pr$V1!="(Intercept)",]
          ao_sten_feat_pr$s1_abs <- abs(ao_sten_feat_pr$s1)
          ao_sten_feat_pr <- ao_sten_feat_pr[order(ao_sten_feat_pr$s1_abs),]
          
      ao_stenplot_cl <- ggplot(data=ao_sten_feat_cl, aes(x=fct_reorder(name_form, s1_weighted), y=s1_weighted)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("AS") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
      ao_stenplot_pr <- ggplot(data=ao_sten_feat_pr, aes(x=fct_reorder(V1, s1_abs), y=s1_abs)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("AS") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),axis.text.y = element_text(size = 5))
      rm(ao_sten_feat_cl, ao_sten_feat_pr)
       
      
        # 8.B.1e - Pics ####
    
      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_all_pr.tiff", width=5300, height=1800, res=300)
        ggarrange(cadplot_pr, hfailplot_pr, afibplot_pr, ao_stenplot_pr, ncol=4, nrow=1)
        dev.off()
      rm(cadplot_cl, hfailplot_cl, afibplot_cl, ao_stenplot_cl, cadplot_pr, hfailplot_pr, afibplot_pr, ao_stenplot_pr)
    
      

      # 8.B.2 - ROC ####
    
      roc_cad_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/3_roc_rawdata_cad_cl.csv")
        roc_cad_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/3_roc_rawdata_cad_pr.csv")
        roc_cad_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/3_roc_rawdata_cad_comb.csv")
        roc_cad_cl$model <- "Clinical parameters"
        roc_cad_pr$model <- "Proteins"
        roc_cad_comb$model <- "Combined"
        roc_cad <- rbind(roc_cad_cl, roc_cad_pr, roc_cad_comb)
        rm(roc_cad_cl, roc_cad_pr, roc_cad_comb)
        
      roc_hfail_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/3_roc_rawdata_hfail_cl.csv")
        roc_hfail_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/3_roc_rawdata_hfail_pr.csv")
        roc_hfail_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/3_roc_rawdata_hfail_comb.csv")
        roc_hfail_cl$model <- "Clinical parameters"
        roc_hfail_pr$model <- "Proteins"
        roc_hfail_comb$model <- "Combined"
        roc_hfail <- rbind(roc_hfail_cl, roc_hfail_pr, roc_hfail_comb)
        rm(roc_hfail_cl, roc_hfail_pr, roc_hfail_comb)
      
      roc_afib_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/3_roc_rawdata_afib_cl.csv")
        roc_afib_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/3_roc_rawdata_afib_pr.csv")
        roc_afib_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/3_roc_rawdata_afib_comb.csv")
        roc_afib_cl$model <- "Clinical parameters"
        roc_afib_pr$model <- "Proteins"
        roc_afib_comb$model <- "Combined"
        roc_afib <- rbind(roc_afib_cl, roc_afib_pr, roc_afib_comb)
        rm(roc_afib_cl, roc_afib_pr, roc_afib_comb)
      
      roc_ao_sten_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/3_roc_rawdata_ao_sten_cl.csv")
        roc_ao_sten_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/3_roc_rawdata_ao_sten_pr.csv")
        roc_ao_sten_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/3_roc_rawdata_ao_sten_comb.csv")
        roc_ao_sten_cl$model <- "Clinical parameters"
        roc_ao_sten_pr$model <- "Proteins"
        roc_ao_sten_comb$model <- "Combined"
        roc_ao_sten <- rbind(roc_ao_sten_cl, roc_ao_sten_pr, roc_ao_sten_comb)
        rm(roc_ao_sten_cl, roc_ao_sten_pr, roc_ao_sten_comb)
      
      
      rocplot_cad <- ggplot(roc_cad, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "True positive rate") +
                        ggtitle("Coronary artery disease") +
                        coord_equal() +
                        theme_classic()
      rocplot_hfail <- ggplot(roc_hfail, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "") +
                        ggtitle("Heart failure") +
                        coord_equal() +
                        theme_classic()
      rocplot_af <- ggplot(roc_afib, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "") +
                        ggtitle("Atrial fibrillation") +
                        coord_equal() +
                        theme_classic()
      rocplot_ao_sten <- ggplot(roc_ao_sten, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "") +
                        ggtitle("Aortic stenosis") +
                        coord_equal() +
                        theme_classic()
      
      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/4_rocs.tiff", width=5300, height=1500, res=300)
        ggarrange(rocplot_cad, rocplot_hfail, rocplot_af, rocplot_ao_sten, common.legend=T, ncol=4, nrow=1)
        dev.off()
      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/4_rocs_cadhfaf.tiff", width=3975, height=1500, res=300)
        ggarrange(rocplot_cad, rocplot_hfail, rocplot_af, common.legend=T, ncol=3, nrow=1)
        dev.off()
  
      rocs <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/4_roc_full.csv")
      rm(rocplot_af, rocplot_cad, rocplot_hfail, rocplot_ao_sten, roc_afib, roc_cad, roc_hfail, roc_ao_sten, rocs)
      
      
      
      # 8.B.3 - KM ####
  
    df <- fread('.../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz')
      set.seed(1234)
      df <- df[sample(1:nrow(df)), ]
      set.seed(1)
      df$id <- 1:nrow(df)
      train <- df %>% dplyr::sample_frac(0.80)
      test  <- dplyr::anti_join(df, train, by = 'id')
      rm(df)
      
    sumstats_info_whi <- fread(".../ukb_proteomics_cvd/input_files/WHI_proteomics_protein_map_2021-09-07.csv")
      names(sumstats_info_whi)[names(sumstats_info_whi) == 'Uniprot ID'] <- 'UniprotID'
      sumstats_info_whi <- sumstats_info_whi %>% mutate(UniprotID = strsplit(as.character(UniprotID), ",")) %>% unnest(UniprotID)
      sumstats_info_ukb <- fread(".../ukb_proteomics_cvd/input_files/2_olink_protein_map_mr.txt")
      sumstats_info_ukb <- sumstats_info_ukb %>% mutate(UniProt = strsplit(as.character(UniProt), "_")) %>% unnest(UniProt)
      sumstats_info_ukb <- sumstats_info_ukb[sumstats_info_ukb$UniProt %in% sumstats_info_whi$UniprotID | sumstats_info_ukb$Assay=="NTproBNP",]
    
    df_pr_test <- test[,106:1564]
      index <- which(colnames(df_pr_test) %in% sumstats_info_ukb$Assay)
      df_pr_test_overl <- df_pr_test[,..index]
      test <- test[,-(106:1564)]
      test <- cbind(test, df_pr_test_overl)
      rm(index, df_pre_test_overl, sumstats_info_whi, sumstats_info_ukb, df_pr_test, df_pr_test_overl)
      
    cad_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_formula_cad_pr.csv")
      hfail_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_formula_hfail_pr.csv")
      afib_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_formula_afib_pr.csv")
      ao_sten_feat_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/2_features_formula_ao_sten_pr.csv")
    test$cad_protriskscore <- as.numeric(cbind(1,as.matrix(test[,139:656])) %*% as.numeric(cad_feat_pr$s1))
      test$hfail_protriskscore <- as.numeric(cbind(1,as.matrix(test[,139:656])) %*% as.numeric(hfail_feat_pr$s1))
      test$afib_protriskscore <- as.numeric(cbind(1,as.matrix(test[,139:656])) %*% as.numeric(afib_feat_pr$s1))
      test$ao_sten_protriskscore <- as.numeric(cbind(1,as.matrix(test[,139:656])) %*% as.numeric(ao_sten_feat_pr$s1))
    test$cad_protriskscore_bin <- exp(test$cad_protriskscore)/(1+exp(test$cad_protriskscore))
      test$hfail_protriskscore_bin <- exp(test$hfail_protriskscore)/(1+exp(test$hfail_protriskscore))
      test$afib_protriskscore_bin <- exp(test$afib_protriskscore)/(1+exp(test$afib_protriskscore))
      test$ao_sten_protriskscore_bin <- exp(test$ao_sten_protriskscore)/(1+exp(test$ao_sten_protriskscore))
    test$cad_protriskscore_quintile <- factor(ntile(test$cad_protriskscore_bin, 5), levels=c(5:1))
      test$hfail_protriskscore_quintile <- factor(ntile(test$hfail_protriskscore_bin, 5), levels=c(5:1))
      test$afib_protriskscore_quintile <- factor(ntile(test$afib_protriskscore_bin, 5), levels=c(5:1))
      test$ao_sten_protriskscore_quintile <- factor(ntile(test$ao_sten_protriskscore_bin, 5), levels=c(5:1))
      
    surv_object <- Surv(time = test$cad_fu, event = test$cad_inc)
    fit1 <- survfit(surv_object ~ cad_protriskscore_quintile, data=test)
    survplot_cad <- ggsurvplot(fit1, pval=F, fun='event', color = "cad_protriskscore_quintile", ylim=c(0,0.22), legend=c(0.21,0.76), palette=c("#A35C14", "#E08529", "#FF9933", "#FFBB77","#FFDDBB"), 
                       censor=F, risk.table=F, ylab="Cumulative incidence of CAD", xlab="Follow-up time (years)", legend.title = "CAD protein score\npercentile:",
                       legend.labs=c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"),ggtheme = theme_classic())$plot
    surv_object <- Surv(time = test$hfail_fu, event = test$hfail_inc)
    fit1 <- survfit(surv_object ~ hfail_protriskscore_quintile, data=test)
    survplot_hfail <- ggsurvplot(fit1, pval=F, fun='event', color = "hfail_protriskscore_quintile", ylim=c(0,0.22), legend=c(0.2,0.76), palette=c("#A35C14", "#E08529", "#FF9933", "#FFBB77","#FFDDBB"), 
                       censor=F, risk.table=F, ylab="Cumulative incidence of HF", xlab="Follow-up time (years)", legend.title = "HF protein score\npercentile:",
                       legend.labs=c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"),ggtheme = theme_classic())$plot
    surv_object <- Surv(time = test$afib_fu, event = test$afib_inc)
    fit1 <- survfit(surv_object ~ afib_protriskscore_quintile, data=test)
    survplot_afib <- ggsurvplot(fit1, pval=F, fun='event', color = "afib_protriskscore_quintile", ylim=c(0,0.22), legend=c(0.2,0.76), palette=c("#A35C14", "#E08529", "#FF9933", "#FFBB77","#FFDDBB"), 
                       censor=F, risk.table=F, ylab="Cumulative incidence of AF", xlab="Follow-up time (years)", legend.title = "AF protein score\npercentile:",
                       legend.labs=c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"),ggtheme = theme_classic())$plot
    surv_object <- Surv(time = test$ao_sten_fu, event = test$ao_sten_inc)
    fit1 <- survfit(surv_object ~ ao_sten_protriskscore_quintile, data=test)
    survplot_ao_sten <- ggsurvplot(fit1, pval=F, fun='event', color = "ao_sten_protriskscore_quintile", ylim=c(0,0.22), legend=c(0.2,0.76), palette=c("#A35C14", "#E08529", "#FF9933", "#FFBB77","#FFDDBB"), 
                       censor=F, risk.table=F, ylab="Cumulative incidence of AS", xlab="Follow-up time (years)", legend.title = "AS protein score\npercentile:",
                       legend.labs=c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"),ggtheme = theme_classic())$plot
            
    
    tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0_overlap_w_whi/4_km_all.tiff", width=6000, height=1300, res=300)
    ggarrange(survplot_cad, survplot_hfail, survplot_afib, survplot_ao_sten, nrow=1, ncol=4)
    dev.off()
    
      
    # 8.C - UKB - NT-proBNP ####
      # 8.C.1 - ROC - only NT-proBNP ####
    
    library(pROC)
    
    df <- fread('.../ukb_proteomics_cvd/input_files/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz')
      set.seed(1234)
      df <- df[sample(1:nrow(df)), ]
      set.seed(1)
      df$id <- 1:nrow(df)
      train <- df %>% dplyr::sample_frac(0.80)
      test  <- dplyr::anti_join(df, train, by = 'id')
      rm(df)
 
      auc_model_pr_hfail <- pROC::roc(response=as.numeric(test$hfail_inc), predictor=as.numeric(test$NTproBNP), levels=c(0,1))
      auc_model_pr_afib <- pROC::roc(response=as.numeric(test$afib_inc), predictor=as.numeric(test$NTproBNP), levels=c(0,1))

      roc_hfail_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0b_onlyntprobnp/3_roc_rawdata_hfail_cl.csv")
        roc_hfail_pr <- data.frame(V1=1:length(auc_model_pr_hfail$specificities), FPR=rev(1-auc_model_pr_hfail$specificities), TPR=rev(auc_model_pr_hfail$sensitivities))
        roc_hfail_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0b_onlyntprobnp/3_roc_rawdata_hfail_comb.csv")
        roc_hfail_cl$model <- "Clinical parameters"
        roc_hfail_pr$model <- "Proteins"
        roc_hfail_comb$model <- "Combined"
        roc_hfail_ntprobnp <- rbind(roc_hfail_cl, roc_hfail_pr, roc_hfail_comb)
        rm(roc_hfail_cl, roc_hfail_pr, roc_hfail_comb)
      
      roc_afib_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0b_onlyntprobnp/3_roc_rawdata_afib_cl.csv")
        roc_afib_pr <- data.frame(V1=1:length(auc_model_pr_afib$specificities), FPR=rev(1-auc_model_pr_afib$specificities), TPR=rev(auc_model_pr_afib$sensitivities))
        roc_afib_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0b_onlyntprobnp/3_roc_rawdata_afib_comb.csv")
        roc_afib_cl$model <- "Clinical parameters"
        roc_afib_pr$model <- "Proteins"
        roc_afib_comb$model <- "Combined"
        roc_afib_ntprobnp <- rbind(roc_afib_cl, roc_afib_pr, roc_afib_comb)
        rm(roc_afib_cl, roc_afib_pr, roc_afib_comb)

      
      rocplot_hfail_ntprobnp <- ggplot(roc_hfail_ntprobnp, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "True positive rate") +
                        ggtitle("NT-proBNP") +
                        coord_equal() +
                        theme_classic()
      rocplot_af_ntprobnp <- ggplot(roc_afib_ntprobnp, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "True positive rate") +
                        ggtitle("NT-proBNP") +
                        coord_equal() +
                        theme_classic()

        
      # 8.C.2 - ROC - no BNP or NT-proBNP ####

      roc_hfail_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0b_nobnp/3_roc_rawdata_hfail_cl.csv")
        roc_hfail_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0b_nobnp/3_roc_rawdata_hfail_pr.csv")
        roc_hfail_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0b_nobnp/3_roc_rawdata_hfail_comb.csv")
        roc_hfail_cl$model <- "Clinical parameters"
        roc_hfail_pr$model <- "Proteins"
        roc_hfail_comb$model <- "Combined"
        roc_hfail_nobnp <- rbind(roc_hfail_cl, roc_hfail_pr, roc_hfail_comb)
        rm(roc_hfail_cl, roc_hfail_pr, roc_hfail_comb)
      
      roc_afib_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0b_nobnp/3_roc_rawdata_afib_cl.csv")
        roc_afib_pr <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0b_nobnp/3_roc_rawdata_afib_pr.csv")
        roc_afib_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/0b_nobnp/3_roc_rawdata_afib_comb.csv")
        roc_afib_cl$model <- "Clinical parameters"
        roc_afib_pr$model <- "Proteins"
        roc_afib_comb$model <- "Combined"
        roc_afib_nobnp <- rbind(roc_afib_cl, roc_afib_pr, roc_afib_comb)
        rm(roc_afib_cl, roc_afib_pr, roc_afib_comb)

      
      rocplot_hfail_nobnp <- ggplot(roc_hfail_nobnp, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "") +
                        ggtitle("All proteins except NT-proBNP and NPPB") +
                        coord_equal() +
                        theme_classic()
      rocplot_af_nobnp <- ggplot(roc_afib_nobnp, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "") +
                        ggtitle("All proteins except NT-proBNP and NPPB") +
                        coord_equal() +
                        theme_classic()

      
      tiff(".../ukb_proteomics_cvd/output_files/7_ntprobnp_hfail.tiff", width=5300/2, height=1500, res=300)
        ggarrange(rocplot_hfail_ntprobnp, rocplot_hfail_nobnp, common.legend=T, legend="bottom", ncol=2, nrow=1)
        dev.off()
      tiff(".../ukb_proteomics_cvd/output_files/7_ntprobnp_afib.tiff", width=5300/2, height=1500, res=300)
        ggarrange(rocplot_af_ntprobnp, rocplot_af_nobnp, common.legend=T, legend="bottom", ncol=2, nrow=1)
        dev.off()
       
      roc_hfail_nobnp$group <- "No NT-proBNP or NPPB"
      roc_hfail_ntprobnp$group <- "NT-proBNP"
      roc_afib_nobnp$outc <- "No NT-proBNP or NPPB"
      roc_afib_ntprobnp$outc <- "NT-proBNP"
        
      write.csv(rbind(roc_afib_ntprobnp, roc_afib_nobnp), ".../ukb_proteomics_cvd/output_files/source_data/ext_data_fig_9.csv", row.names=F)
      write.csv(rbind(roc_hfail_ntprobnp, roc_hfail_nobnp), ".../ukb_proteomics_cvd/output_files/source_data/ext_data_fig_10.csv", row.names=F)
       
      
