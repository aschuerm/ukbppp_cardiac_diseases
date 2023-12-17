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
 

  ### 1 - Baseline ####  
    # 1.A - Between-group comparisons (e.g., Table 1) ####  

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
  tiff('.../ukb_proteomics_cvd/output_file/0_venn_assoc_timevar_excl.tiff', width=3000, height=2100, res=300)
  ggplot() +
    geom_sf(aes(fill = name), data = venn_region(d2)) +
    geom_sf(aes(color = name), data = venn_setedge(d)) +
    geom_sf_text(aes(label = name), data = venn_setlabel(d)) +
    geom_sf_text(aes(label = count), data = venn_region(d)) +
    scale_color_manual(values = alpha(c("#DC0000B2", "navy", "darkgreen", "grey"), .1), name="Outcome:") +
    scale_fill_manual(values = alpha(c("#DC0000B2", "navy", "darkgreen", "grey"), .1), name="Outcome:") +
    theme_void() + labs(title="Time-varying covariates and exclusion of all prevalent cases:")
  dev.off()


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
  sig <- 0.05/5836

    # 3.B - CAD ####

  df_sub <- df[df$Outcome=="CAD",]
  df_sub$col <- factor(ifelse(df_sub$P_Value<sig & df_sub$Protein%in%df[df$Outcome %in% c("HF", "Afib", "AS")&df$P_Value<sig,]$Protein,
                        "A", ifelse(df_sub$P_Value<sig, "B", df_sub$chr)), levels=c(1:23, "A", "B"))
  df_sub$Protein_lab <- ifelse(df_sub$P_Value<5e-15, df_sub$Protein, "")
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
  
  library(rvg)
  library(officer)
  
    Comp_1_dml <- dml(code = print(plot_cad, newpage = FALSE))
    Comp_2_dml <- dml(code = print(plot_hf, newpage = FALSE))
    Comp_3_dml <- dml(code = print(plot_af, newpage = FALSE))
    Comp_4_dml <- dml(code = print(plot_af_add1, newpage = FALSE))
    Comp_5_dml <- dml(code = print(plot_af_add2, newpage = FALSE))
    Comp_6_dml <- dml(code = print(plot_as, newpage = FALSE))
    
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_1_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_2_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_3_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_4_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_5_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_6_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")

    rm(plot_af, plot_af_add1, plot_af_add2, plot_as, plot_cad, plot_hf, 
       Comp_1_dml, Comp_2_dml, Comp_3_dml, Comp_4_dml, Comp_5_dml, Comp_6_dml, df, df_sub, link, axis_set, sig, ylim, data_cum)
    
    
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

    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_1_dml, ph_location(left = 1, top = 1, width = 4, height = 7)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_2_dml, ph_location(left = 1, top = 1, width = 4, height = 7)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_3_dml, ph_location(left = 1, top = 1, width = 4, height = 7)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_4_dml, ph_location(left = 1, top = 1, width = 4, height = 7)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")

    rm(cad, af, hf, as, Comp_1_dml, Comp_2_dml, Comp_3_dml, Comp_4_dml)


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
  
    Comp_1_dml <- dml(code = print(a, newpage = FALSE))
    Comp_2_dml <- dml(code = print(b, newpage = FALSE))
    Comp_3_dml <- dml(code = print(c, newpage = FALSE))
    Comp_4_dml <- dml(code = print(d, newpage = FALSE))

    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_1_dml, ph_location(left = 1, top = 1, width = 6, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_2_dml, ph_location(left = 1, top = 1, width = 6, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_3_dml, ph_location(left = 1, top = 1, width = 6, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_4_dml, ph_location(left = 1, top = 1, width = 6, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")

  
  
    # 5.A - Interaction - top 5 ####

  df <- df[order(df$P_Val_Int),]
  df$num_peroutc <- ave(df$P_Val_Int, df$Outcome, FUN = seq_along)
  df <- df[df$num_peroutc %in% 1:5,]
  
  df_long <- data.frame(protein=rep(df$Protein, 2), sex=c(rep("F", 20), rep("M", 20)), hr=c(df$HR_F, df$HR_M), lci=c(df$CI_Low_F, df$CI_Low_M),
                        uci=c(df$CI_High_F, df$CI_High_M), pval=c(df$P_Val_F, df$P_Val_M), rank=rep(df$num_peroutc, 2), outcome=rep(df$Outcome, 2))
  df_long <- df_long[order( df_long$outcome, df_long$sex, df_long$rank),]

  
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
  tiff(".../ukb_proteomics_cvd/output_files/4_diff_sexes_hrs.tiff", width=1650, height=1200, res=300)
  ggarrange(a,b,c,d, common.legend=T)
  dev.off()
  
  library(rvg)
  library(officer)
  
    Comp_1_dml <- dml(code = print(a, newpage = FALSE))
    Comp_2_dml <- dml(code = print(b, newpage = FALSE))
    Comp_3_dml <- dml(code = print(c, newpage = FALSE))
    Comp_4_dml <- dml(code = print(d, newpage = FALSE))

    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_1_dml, ph_location(left = 1, top = 1, width = 3, height = 2)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_2_dml, ph_location(left = 1, top = 1, width = 3, height = 2)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_3_dml, ph_location(left = 1, top = 1, width = 3, height = 2)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_4_dml, ph_location(left = 1, top = 1, width = 3, height = 2)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")

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

  instr <- fread("../tmp/outcomes/primaryanalyses_loop_corr_instr.csv")
  res <- fread("../tmp/outcomes/primaryanalyses_loop_corr.csv")
  fstat <- fread("../tmp/outcomes/rsq_fstat_combined_estimated_and_allelescore.csv")
    res <- res[!(res$exp %in% fstat$exposure[fstat$f<10]),]
    res <- res[res$method=="Inverse variance weighted (correlation inc)" | res$method=="Wald ratio",]
    res_dupl <- res[,c("exp", "outc")]
    res <- res[!duplicated( res[,c("exp", "outc")]),]
      rm(res_dupl)
    res$outc <- factor(toupper(res$outc), levels=c("CAD", "HF", "AF", "AS"))
    res$col <- ifelse(res$pval<0.05 & res$b<0, "Protective", ifelse(res$pval<0.05 & res$b>0, "Detrimental", "Other"))
    # write.csv(res, ".../ukb_proteomics_cvd/output_files/3_MR_allfindings.csv")
    # res_sign <- res[res$col!="Other",]
    # write.csv(res_sign, ".../ukb_proteomics_cvd/output_files/3_MR_signfindings.csv")
    # rm(res_sign)
  
  sign_prim <- fread('.../ukb_proteomics_cvd/output_files/1_adj_no_prev_timevar_model_hr.csv')
    sign_prim <- sign_prim[sign_prim$P_Value<0.05/5836,]
    sign_prim$Outcome <- ifelse(sign_prim$Outcome == "Afib", "AF", sign_prim$Outcome)
    sign_prim$beta_prim_obs <- log(sign_prim$HR)
    res <- merge(res, sign_prim[,c("Protein", "Outcome", "beta_prim_obs")], by.x=c("exp", "outc"), by.y=c("Protein", "Outcome"), all.x=T)
    rm(sign_prim)
    
  sens <- fread("../tmp/outcomes/sensitivityanalyses_loop_corr.csv")
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

  sens_4 <- fread("../tmp/outcomes/sensitivityanalysis_allelescore_oneproteinpermodel.csv")
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
      
    # write.csv(res, "../tmp/outcomes/sensitivityanalyses_posttwosamplesens_postallele.csv")
      
  all_mv <- fread("../tmp/outcomes/regression_prs_proteins_output_multiv_analyses.csv")
    colnames(all_mv) <- c("exp", "outc", "b_allelescore_mv", "se", "pval")
    all_mv$exp <- gsub("prs_", "", all_mv$exp)
    all_mv$exp <- gsub("HLA_E", "HLA-E", all_mv$exp)
    all_mv$outc <- ifelse(all_mv$outc=="cad", "CAD", ifelse(all_mv$outc=="hfail", "HF", ifelse(all_mv$outc=="afib", "AF", "AS")))
  
  res <- merge(res, all_mv[,c("exp", "outc", "b_allelescore_mv")], by=c("exp", "outc"), all.x=T, all.y=F)
    res$sens_robust_allelescore_mv <- ifelse(res$b_allelescore_mv*res$b > 0, T, F)

  res$all_sens <- ifelse(is.na(res$sens_robust_allelescore_mv), res$all_sens, ifelse(res$sens_robust_allelescore_mv==F, "No", res$all_sens))
    res$directional_obs <- ifelse(res$b*res$beta_prim_obs>0, "Consistent", "Not consistent")
    res$group <- ifelse(res$directional_obs == "Consistent" & res$all_sens == "Yes", "Robust / consistent with observational findings",
                          ifelse(res$directional_obs == "Not consistent" & res$all_sens == "Yes", "Robust / not consistent with observational findings",
                          ifelse(res$directional_obs == "Consistent" & res$all_sens == "No", "Not robust / consistent with observational findings",
                          ifelse(res$directional_obs == "Not consistent" & res$all_sens == "No", "Not robust / not consistent with observational findings",
                          "Nonsignificant"))))
    # write.csv(res, "../tmp/outcomes/sensitivityanalyses_posttwosamplesens_postallele_postallelemv.csv")
  rm(directional_consistency_df, fstat, instr, sens_3, sens_4)
      
      
      
  a <- ggplot(res[res$outc=="CAD",], aes(b, -log10(pval), color=group)) +
      scale_color_manual(values=c("grey90", "#FFF1B9", "#ECB5A2", "#E8B900", "red3"), name="") +
      labs(y="-log<sub>10</sub>(P)", x="log(OR)") +
      geom_point()+
      theme_classic() +
      scale_y_continuous(limits=c(0,10.5)) +
      scale_x_continuous(limits=c(-0.75, 0.75)) +
      guides(fill=guide_legend(override.aes = aes(label = NA), ncol=2, nrow=2)) +
      ggtitle("Coronary artery disease") +
      theme(axis.title.y = element_markdown(),axis.title.x = element_markdown()) +
      geom_text_repel(aes(label=ifelse(res[res$outc=="CAD",]$all_sens == "Yes" & (res[res$outc=="CAD",]$pval<0.001|abs(res[res$outc=="CAD",]$b)>0.3), res[res$outc=="CAD",]$exp, "")), size=3, force=10, max.overlaps=100)
  a_bis <- ggplot(res[res$outc=="CAD",], aes(b, -log10(pval), color=group)) +
      scale_color_manual(values=c("grey90", "#FFF1B9", "#ECB5A2", "#E8B900", "red3"), name="") +
      labs(y="-log<sub>10</sub>(P)", x="log(OR)") +
      geom_point()+
      theme_classic() +
      scale_y_continuous(limits=c(12,13)) +
      scale_x_continuous(limits=c(-0.75, 0.75)) +
      guides(fill=guide_legend(override.aes = aes(label = NA), ncol=2, nrow=2)) +
      ggtitle("Coronary artery disease") +
      theme(axis.title.y = element_markdown(),axis.title.x = element_markdown()) +
      geom_text_repel(aes(label=ifelse(res[res$outc=="CAD",]$all_sens == "Yes" & (res[res$outc=="CAD",]$pval<0.001|abs(res[res$outc=="CAD",]$b)>0.3), res[res$outc=="CAD",]$exp, "")), size=3, force=10, max.overlaps=100)
  b <- ggplot(res[res$outc=="HF",], aes(b, -log10(pval), color=group)) +
      scale_color_manual(values=c("grey90", "#FFF1B9", "#ECB5A2", "#E8B900", "red3"), name="Outcome") +
      labs(y="-log<sub>10</sub>(P)", x="log(OR)") +
      geom_point()+
      theme_classic() +
      scale_y_continuous(limits=c(0,10.5)) +
      scale_x_continuous(limits=c(-0.75, 0.75)) +
      guides(fill=guide_legend(override.aes = aes(label = NA), ncol=2, nrow=2)) +
      ggtitle("Heart failure") +
      theme(axis.title.y = element_markdown(),axis.title.x = element_markdown()) +
      geom_text_repel(aes(label=ifelse(res[res$outc=="HF",]$all_sens == "Yes" & (res[res$outc=="HF",]$pval<0.001|abs(res[res$outc=="HF",]$b)>0.3), res[res$outc=="HF",]$exp, "")), size=3, force=10, max.overlaps=100)
  b_bis <- ggplot(res[res$outc=="HF",], aes(b, -log10(pval), color=group)) +
      scale_color_manual(values=c("grey90", "#FFF1B9", "#ECB5A2", "#E8B900", "red3"), name="Outcome") +
      labs(y="-log<sub>10</sub>(P)", x="log(OR)") +
      geom_point()+
      theme_classic() +
      scale_y_continuous(limits=c(15,16)) +
      scale_x_continuous(limits=c(-0.75, 0.75)) +
      guides(fill=guide_legend(override.aes = aes(label = NA), ncol=2, nrow=2)) +
      ggtitle("Heart failure") +
      theme(axis.title.y = element_markdown(),axis.title.x = element_markdown()) +
      geom_text_repel(aes(label=ifelse(res[res$outc=="HF",]$all_sens == "Yes" & (res[res$outc=="HF",]$pval<0.001|abs(res[res$outc=="HF",]$b)>0.3), res[res$outc=="HF",]$exp, "")), size=3, force=10, max.overlaps=100)
  c <- ggplot(res[res$outc=="AF",], aes(b, -log10(pval), color=group)) +
      scale_color_manual(values=c("grey90", "#FFF1B9", "#ECB5A2", "#E8B900", "red3"), name="Outcome") +
      labs(y="-log<sub>10</sub>(P)", x="log(OR)") +
      geom_point()+
      theme_classic() + 
      scale_y_continuous(limits=c(0,10.5)) +
      scale_x_continuous(limits=c(-0.75, 0.75)) +
      guides(fill=guide_legend(override.aes = aes(label = NA), ncol=2, nrow=2)) +
      ggtitle("Atrial fibrillation") +
      theme(axis.title.y = element_markdown(),axis.title.x = element_markdown()) +
      geom_text_repel(aes(label=ifelse(res[res$outc=="AF",]$all_sens == "Yes" & res[res$outc=="AF",]$pval<0.01, res[res$outc=="AF",]$exp, "")), size=3, force=10, max.overlaps=100)
  d <- ggplot(res[res$outc=="AS",], aes(b, -log10(pval), color=group)) +
      scale_color_manual(values=c("grey90", "#FFF1B9", "#E8B900", "red3"), name="Outcome") +
      labs(y="-log<sub>10</sub>(P)", x="log(OR)") +
      geom_point()+
      theme_classic() +
      scale_y_continuous(limits=c(0,10.5)) +
      guides(fill=guide_legend(override.aes = aes(label = NA), ncol=2, nrow=2)) +
      scale_x_continuous(limits=c(-0.75, 0.75)) +
      ggtitle("Aortic stenosis") +
      theme(axis.title.y = element_markdown(),axis.title.x = element_markdown()) +
      geom_text_repel(aes(label=ifelse(res[res$outc=="AS",]$all_sens == "Yes" & res[res$outc=="AS",]$pval<0.05, res[res$outc=="AS",]$exp, "")), size=3, force=10, max.overlaps=100)

  library(rvg)
  library(officer)
  
    Comp_1_dml <- dml(code = print(a, newpage = FALSE))
    Comp_2_dml <- dml(code = print(a_bis, newpage = FALSE))
    Comp_3_dml <- dml(code = print(b, newpage = FALSE))
    Comp_4_dml <- dml(code = print(b_bis, newpage = FALSE))
    Comp_5_dml <- dml(code = print(c, newpage = FALSE))
    Comp_6_dml <- dml(code = print(d, newpage = FALSE))
    
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_1_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_2_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_3_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_4_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_5_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")
    read_pptx(".../ukb_proteomics_cvd/output_file/Figures.pptx") %>%
      add_slide() %>%
      ph_with(Comp_6_dml, ph_location(left = 1, top = 1, width = 8, height = 4)) %>%
      print(target = ".../ukb_proteomics_cvd/output_file/Figures.pptx")

  rm(list=ls())

  
    # 6.A - Create suppl table with sensitivity analyses ####

  
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggtext)
  library(tidyr)

  instr <- fread(".../tmp/outcomes//primaryanalyses_loop_corr_instr.csv")
  res <- fread(".../tmp/outcomes//primaryanalyses_loop_corr.csv")
    res <- res[res$method=="Inverse variance weighted (correlation inc)" | res$method=="Wald ratio",]
    res_dupl <- res[,c("exp", "outc")]
    res <- res[!duplicated( res[,c("exp", "outc")]),]
      rm(res_dupl)
    res$outc <- factor(toupper(res$outc), levels=c("CAD", "HF", "AF", "AS"))
    res$col <- ifelse(res$pval<0.05 & res$b<0, "Protective", ifelse(res$pval<0.05 & res$b>0, "Detrimental", "Other"))

  sign_prim <- fread('.../ukb_proteomics_cvd/output_files/1_adj_no_prev_timevar_model_hr.csv')
    sign_prim <- sign_prim[sign_prim$P_Value<0.05/5836,]
    sign_prim$Outcome <- ifelse(sign_prim$Outcome == "Afib", "AF", sign_prim$Outcome)
    sign_prim$beta_prim_obs <- log(sign_prim$HR)
    res <- merge(res, sign_prim[,c("Protein", "Outcome", "beta_prim_obs")], by.x=c("exp", "outc"), by.y=c("Protein", "Outcome"), all.x=T)
    rm(sign_prim)
    
  sens <- fread(".../tmp/outcomes//sensitivityanalyses_loop_corr.csv")
    sens$outc <- toupper(sens$outc)
    sens <- merge(sens, res[res$pval<0.05, c("exp", "outc")], by=c("exp", "outc"), all.x=F, all.y=F)
    sens_1 <- sens[sens$method=="Inverse variance weighted (correlation inc)" | sens$method=="Wald ratio",]
    sens_2 <- sens[sens$method=="Egger (correlation inc)" & sens$pvalthreshold==5e-06 & sens$rsqthreshold==0.1,]
    sens_3 <- sens[sens$method=="Egger intercept (correlation inc)" & sens$pvalthreshold==5e-06 & sens$rsqthreshold==0.1,]
    sens <- rbind(sens_1, sens_2, sens_3)
    sens <- sens[!duplicated( sens[,c("exp", "outc","pvalthreshold", "rsqthreshold", "method")]),]
    rm(sens_1, sens_2, sens_3)
  
  sens$method <- ifelse(sens$method=="Wald ratio", "IVW/Wald", ifelse(sens$method=="Inverse variance weighted (correlation inc)", "IVW/Wald", 
                        ifelse(sens$method=="Egger intercept (correlation inc)", "Egg_Int",
                        ifelse(sens$method=="Egger (correlation inc)", "Egg", NA))))
    sens$or <- exp(sens$b)
    sens$b_ci <- paste(round(sens$b+qnorm(0.025)*sens$se, 3), "to", round(sens$b+qnorm(0.975)*sens$se, 3))
    sens$or_ci <- paste(round(exp(sens$b+qnorm(0.025)*sens$se), 3), "to", round(exp(sens$b+qnorm(0.975)*sens$se), 3))
    sens_w <- sens[,-c("V1", "se")] %>% pivot_wider(names_from=c(pvalthreshold, rsqthreshold, method), values_from=c(nsnp, b, b_ci, or, or_ci, pval))
    # write.csv(sens_w, ".../ukb_proteomics_cvd/output_files/3_MR_sens_formatted.csv")
  
  
  
    # 6.B - Coloc table ####
    
    
  library(data.table)

  res <- fread(".../tmp/outcomes/sensitivityanalyses_posttwosamplesens_postallele_postallelemv.csv")
    res <- res[res$method=="Inverse variance weighted (correlation inc)" | res$method=="Wald ratio",]
    res <- res[res$all_sens=="Yes",]
    res <- res[,c("exp", "outc")]
  coloc <- fread(".../tmp/outcomes/coloc_loop.csv")
    coloc$outc <- toupper(coloc$outc)
  coloc <- merge(res, coloc, all.x=F, all.y=F, by=c("exp", "outc"))
    
  write.csv(coloc,".../tmp/outcomes/coloc_proteins_surviving_sensanalyses.csv")

  
  

  
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
                        ggtitle("Coronary artery disease") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
      cadplot_pr <- ggplot(data=cad_feat_pr, aes(x=fct_reorder(V1, s1_abs), y=s1_abs)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("Coronary artery disease") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),axis.text.y = element_text(size = 5))
      rm(cad_feat_cl, cad_feat_pr)
       
        
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
                        ggtitle("Heart failure") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
      hfailplot_pr <- ggplot(data=hfail_feat_pr, aes(x=fct_reorder(V1, s1_abs), y=s1_abs)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("Heart failure") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),axis.text.y = element_text(size = 5))
      rm(hfail_feat_cl, hfail_feat_pr)
                      
        
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
                        ggtitle("Atrial fibrillation") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
      afibplot_pr <- ggplot(data=afib_feat_pr, aes(x=fct_reorder(V1, s1_abs), y=s1_abs)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("Atrial fibrillation") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),axis.text.y = element_text(size = 5))
      rm(afib_feat_cl, afib_feat_pr)
       
      
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
                        ggtitle("Aortic stenosis") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown())
      ao_stenplot_pr <- ggplot(data=ao_sten_feat_pr, aes(x=fct_reorder(V1, s1_abs), y=s1_abs)) +
                        geom_col(color="black", fill="black") +
                        theme_classic() +
                        labs(x="", y="Regression coefficient in prediction model") +
                        ggtitle("Aortic stenosis") +
                        coord_flip() +
                        theme(axis.title.y = element_markdown(),axis.title.x = element_markdown(),axis.text.y = element_text(size = 5))
      rm(ao_sten_feat_cl, ao_sten_feat_pr)
       
      
        # 8.A.1e - Pics ####
    
      library(forcats)
      library(ggplot2)
      library(ggpubr)
      
      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_all_cl.tiff", width=5300, height=1800, res=300)
        ggarrange(cadplot_cl, hfailplot_cl, afibplot_cl, ao_stenplot_cl, ncol=4, nrow=1)
        dev.off()
      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_all_pr.tiff", width=5300, height=1800, res=300)
        ggarrange(cadplot_pr, hfailplot_pr, afibplot_pr, ao_stenplot_pr, ncol=4, nrow=1)
        dev.off()
      rm(cadplot_cl, hfailplot_cl, afibplot_cl, ao_stenplot_cl, cadplot_pr, hfailplot_pr, afibplot_pr, ao_stenplot_pr)
    
      
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
                        theme_classic()
      rocplot_hfail <- ggplot(roc_hfail, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "") +
                        ggtitle("Heart failure") +
                        theme_classic()
      rocplot_af <- ggplot(roc_afib, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "") +
                        ggtitle("Atrial fibrillation") +
                        theme_classic()
      rocplot_ao_sten <- ggplot(roc_ao_sten, aes(x=FPR, y=TPR, color=model)) +
                        geom_line() +
                        geom_abline(linetype = "dotted", color = "grey90") +
                        scale_color_manual(values=c("black", "#1EC000", "#C2EEB9"), name="") +
                        labs(x = "False positive rate", y = "") +
                        ggtitle("Aortic stenosis") +
                        theme_classic()
      
      tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_rocs.tiff", width=4000, height=1000, res=300)
        ggarrange(rocplot_cad, rocplot_hfail, rocplot_af, rocplot_ao_sten, common.legend=T, ncol=4, nrow=1, legend="bottom")
        dev.off()
  
      rocs <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_roc_full.csv")
      rm(rocplot_af, rocplot_cad, rocplot_hfail, rocplot_ao_sten, roc_afib, roc_cad, roc_hfail, roc_ao_sten, rocs)
      
      
      
      
      # 8.A.3 - KM ####
  
    library(survival)
    library(ggpubr)
    library(survminer)
      
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
    # cad_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_cl.csv")
    #   hfail_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_cl.csv")
    #   afib_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_cl.csv")
    #   ao_sten_feat_cl <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_cl.csv")
    # cad_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_cad_comb.csv")
    #   hfail_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_hfail_comb.csv")
    #   afib_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_afib_comb.csv")
    #   ao_sten_feat_comb <- fread(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_ao_sten_comb.csv")
    # cad_table <- merge(cad_feat_pr, cad_feat_cl, by="V1", all.x=T, all.y=T)
    #   cad_table <- merge(cad_table, cad_feat_comb, by="V1", all.x=T, all.y=T)
    #   hfail_table <- merge(hfail_feat_pr, hfail_feat_cl, by="V1", all.x=T, all.y=T)
    #   hfail_table <- merge(hfail_table, hfail_feat_comb, by="V1", all.x=T, all.y=T)
    #   afib_table <- merge(afib_feat_pr, afib_feat_cl, by="V1", all.x=T, all.y=T)
    #   afib_table <- merge(afib_table, afib_feat_comb, by="V1", all.x=T, all.y=T)
    #   ao_sten_table <- merge(ao_sten_feat_pr, ao_sten_feat_cl, by="V1", all.x=T, all.y=T)
    #   ao_sten_table <- merge(ao_sten_table, ao_sten_feat_comb, by="V1", all.x=T, all.y=T)
    # table <- data.frame(V1=c(cad_feat_pr$V1, cad_feat_cl$V1[-1]))
    #   table <- merge(table, cad_table, by="V1", all.x=T, all.y=T)
    #   table <- merge(table, hfail_table, by="V1", all.x=T, all.y=T)
    #   table <- merge(table, afib_table, by="V1", all.x=T, all.y=T)
    #   colnames(table) <- c("V1", "cad_pr", "cad_cl", "cad_comb", "hfail_pr", "hfail_cl", "hfail_comb", "afib_pr", "afib_cl", "afib_comb")
    #   table <- merge(table, ao_sten_table, by="V1", all.x=T, all.y=T)
    #   colnames(table) <- c("Protein", "cad_pr", "cad_cl", "cad_comb", "hfail_pr", "hfail_cl", "hfail_comb", "afib_pr", "afib_cl", "afib_comb", "ao_sten_pr", "ao_sten_cl", "ao_sten_comb")
    #   write.csv(table, ".../ukb_proteomics_cvd/output_files/5_modeling_lasso/2_features_formula_all.csv", row.names=F)
    #   rm(list=ls())
      
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
                       censor=F, risk.table=F, ylab="Cumulative incidence", xlab="Follow-up time (years)", legend.title = "Protein score\npercentile:",
                       legend.labs=c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"),ggtheme = theme_classic())$plot
    surv_object <- Surv(time = test$hfail_fu, event = test$hfail_inc)
    fit1 <- survfit(surv_object ~ hfail_protriskscore_quintile, data=test)
    survplot_hfail <- ggsurvplot(fit1, pval=F, fun='event', color = "hfail_protriskscore_quintile", ylim=c(0,0.22), legend=c(0.2,0.76), palette=c("#A35C14", "#E08529", "#FF9933", "#FFBB77","#FFDDBB"), 
                       censor=F, risk.table=F, ylab="", xlab="Follow-up time (years)", legend.title = "Protein score\npercentile:",
                       legend.labs=c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"),ggtheme = theme_classic())$plot
    surv_object <- Surv(time = test$afib_fu, event = test$afib_inc)
    fit1 <- survfit(surv_object ~ afib_protriskscore_quintile, data=test)
    survplot_afib <- ggsurvplot(fit1, pval=F, fun='event', color = "afib_protriskscore_quintile", ylim=c(0,0.22), legend=c(0.2,0.76), palette=c("#A35C14", "#E08529", "#FF9933", "#FFBB77","#FFDDBB"), 
                       censor=F, risk.table=F, ylab="", xlab="Follow-up time (years)", legend.title = "Protein score\npercentile:",
                       legend.labs=c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"),ggtheme = theme_classic())$plot
    surv_object <- Surv(time = test$ao_sten_fu, event = test$ao_sten_inc)
    fit1 <- survfit(surv_object ~ ao_sten_protriskscore_quintile, data=test)
    survplot_ao_sten <- ggsurvplot(fit1, pval=F, fun='event', color = "ao_sten_protriskscore_quintile", ylim=c(0,0.22), legend=c(0.2,0.76), palette=c("#A35C14", "#E08529", "#FF9933", "#FFBB77","#FFDDBB"), 
                       censor=F, risk.table=F, ylab="", xlab="Follow-up time (years)", legend.title = "Protein score\npercentile:",
                       legend.labs=c("80-100%", "60-80%", "40-60%", "20-40%", "0-20%"),ggtheme = theme_classic())$plot
            
    tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_km_all.tiff", width=4000, height=1000, res=300)
    ggarrange(survplot_cad, survplot_hfail, survplot_afib, survplot_ao_sten, nrow=1, ncol=4, common.legend=T, legend = "bottom")
    dev.off()
    
    
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
    cad_plot <- ggplot(data=df_cad, aes(x=quant*10-5, y=1000*inc_r, ymin=1000*inc_r_lci, ymax=1000*inc_r_uci, color=quant*10-5)) +
                  geom_pointrange() +
                  scale_color_gradient2(low="#FFDDBB", mid="#FF9933", high="#A35C14", midpoint = 50, limits=c(0,100))  + 
                  scale_x_continuous(limits=c(0,100), name="Protein score percentile") +
                  scale_y_continuous(name="Incidence rate per 1,000 years (95% CI)", limits = c(-1, 27)) +
                  ggtitle("Coronary artery disease") +
                  theme_classic() +
                  theme(legend.position = "none")
      
    df_hfail <- data.frame(quant=NA, inc_r=NA, inc_r_lci=NA, inc_r_uci=NA)[-1,]
    for (i in 1:10) {7
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
    hfail_plot <- ggplot(data=df_hfail, aes(x=quant*10-5, y=1000*inc_r, ymin=1000*inc_r_lci, ymax=1000*inc_r_uci, color=quant*10-5)) +
                  geom_pointrange() +
                  scale_color_gradient2(low="#FFDDBB", mid="#FF9933", high="#A35C14", midpoint = 50, limits=c(0,100))  + 
                  scale_x_continuous(limits=c(0,100), name="Protein score percentile") +
                  scale_y_continuous(name="", limits = c(-1, 27)) +
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
    afib_plot <- ggplot(data=df_afib, aes(x=quant*10-5, y=1000*inc_r, ymin=1000*inc_r_lci, ymax=1000*inc_r_uci, color=quant*10-5)) +
                  geom_pointrange() +
                  scale_color_gradient2(low="#FFDDBB", mid="#FF9933", high="#A35C14", midpoint = 50, limits=c(0,100))  + 
                  scale_x_continuous(limits=c(0,100), name="Protein score percentile") +
                  scale_y_continuous(name="", limits = c(-1, 27)) +
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
    ao_sten_plot <- ggplot(data=df_ao_sten, aes(x=quant*10-5, y=1000*inc_r, ymin=1000*inc_r_lci, ymax=1000*inc_r_uci, color=quant*10-5)) +
                  geom_pointrange() +
                  scale_color_gradient2(low="#FFDDBB", mid="#FF9933", high="#A35C14", midpoint = 50, limits=c(0,100))  + 
                  scale_x_continuous(limits=c(0,100), name="Protein score percentile") +
                  scale_y_continuous(name="", limits = c(-1, 27)) +
                  ggtitle("Aortic stenosis") +
                  theme_classic() +
                  theme(legend.position = "none")

    
    library(ggpubr)
    
    tiff(".../ukb_proteomics_cvd/output_files/5_modeling_lasso/4_quantileplots_all.tiff", width=4000, height=1000, res=300)
    ggarrange(cad_plot, hfail_plot, afib_plot, ao_sten_plot, nrow=1, ncol=4)
    dev.off()