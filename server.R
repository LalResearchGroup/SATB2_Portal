################################################
################# SATB2 PORTAL ##################
################################################
################## LIBRARIES ################
#setwd("C:/Users/tobia/OneDrive/Studium Biologie/Praktikum/SATB2_Portal/")

library(shiny)
library(plotly)
library(readxl)
library(readr)
library(r3dmol)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(tidyverse)
library(seqinr)
library(bio3d)
library(RColorBrewer)
library(plyr)
library(ggwordcloud)

############## FUNCTIONS, STYLE ############

numextract <- function(string) {
    as.numeric(str_extract(string, "\\-*\\d+\\.*\\d*"))
}

convert_aa <- function(aa){
  
  return(ifelse(aa == "*","Stop",aaa(aa)))
  
}

log_fun <- function(x) { 
    x = case_when(
        as.numeric(x) > 0 ~ log10(x), #log10 transformation to make the visuals clean 
        as.numeric(x) == 0 ~ 0,
        as.numeric(x) < 0 ~ log10(x * -1) * -1, 
        TRUE ~ -9999)
    return(x)
    
}

plotly_font <- list(
    family = "sans-serif",
    size = 15)

goodbye <- c("zoomIn2d", "zoomOut2d", "hoverClosestGeo", "hoverClosestGl2d",
             "hoverClosestCartesian","lasso2d","select2d","resetScale2d",
             "hoverCompareCartesian", "hoverClosestPie","toggleSpikelines")

# coloring for colorblindness 
lolliplot_fill_scheme <-  #("Missense"="#D55E00","PTV"="#0072B2","Control" ="#000000", "PER" = "red"))+
  c(
    "N-terminus" = "#CCCCCC",
    "Domain-Linker1" = "#CCCCCC",
    "Domain-Linker2" = "#CCCCCC",
    "Domain-Linker3" = "#CCCCCC",
    "Domain-Linker4" = "#CCCCCC",
    "C-terminus" = "#CCCCCC",
    "ULD Domain" = "#6699CC",
    "CUTL Domain" = "#FF9999",
    "CUT1 Domain" = "#99CC33",
    "CUT2 Domain" = "#FFFF00",
    "HOX" = "#FFCCCC",
    "no" = "#FFFFFF",
    "Missense" = "#D55E00",
    "Synonymous" = "#D55E00",
    "PTV" = "#0072B2",
    "Other" = "#0072B2"
  )



######Functions######

#Basic Information 
Phenotype_fac_1.fun <- function(phenotype,color_sel){
  
  plot_ly(Patient_data.df %>% 
            dplyr::rename(phenotype_sel = phenotype) %>% 
            select(phenotype_sel) %>%
            arrange(phenotype_sel) %>% 
            filter(!is.na(phenotype)) %>% 
            mutate(phenotype_sel = factor(phenotype_sel, levels = c("None","1 to 10","10 to 50","Greater than 50"))) %>% 
            group_by(phenotype_sel) %>% 
            dplyr::summarise(n = n()) %>% 
            mutate(n_gene = sum(n),
                   prop = n/n_gene*100)%>% 
            assign("save",.,envir = .GlobalEnv), 
          x = ~ phenotype_sel, 
          y = ~ round(prop, digits = 2), 
          color = ~ phenotype_sel, 
          colors = color_sel,
          type = "bar", 
          hoverinfo = "text", showlegend = FALSE,
          text= ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
    layout(title="", 
           font=plotly_font,
           xaxis = list(title="",showline = T, tickangle = 45),
           yaxis = list(title="Share of individuals (%)",showline = T),
           margin = list(b = 160)) %>%
    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  
  
}


Phenotype_fac_2.fun <- function(select_gene,phenotype,color_sel){
  
  plot_ly(Patient_data.df %>% 
            dplyr::rename(phenotype_fac = phenotype) %>% 
            filter(phenotype_fac != "NA") %>% 
            mutate(phenotype_fac =ifelse(phenotype_fac == "Yes"," Yes",phenotype_fac)) %>% 
            select(phenotype_fac) %>%
            arrange(phenotype_fac) %>% 
            mutate(phenotype_sel = factor(phenotype_fac, levels = unique(phenotype_fac))) %>% 
            group_by(phenotype_fac) %>% 
            dplyr::summarise(n = n()) %>% 
            mutate(n_gene = sum(n),
                   prop = n/n_gene*100)%>% 
            assign("save",.,envir = .GlobalEnv),
          x = ~ phenotype_fac, 
          y = ~ round(prop, digits = 2), 
          color = ~ phenotype_fac, 
          colors = color_sel,
          type = "bar", 
          hoverinfo = "text", showlegend = FALSE,
          text= ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
    layout(title="", 
           font=plotly_font,
           xaxis = list(title="",showline = T, tickangle = 45),
           yaxis = list(title="Share of individuals (%)",showline = T),
           margin = list(b = 160)) %>%
    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  
  
}

Phenotype_fac_3.fun <- function(phenotype,color_sel){
  
  plot_ly(Patient_data.df %>% 
            dplyr::rename(phenotype_fac = phenotype) %>% 
            filter(phenotype_fac != "NA") %>% 
            mutate(phenotype_fac =ifelse(phenotype_fac == "Yes"," Yes",phenotype_fac)) %>% 
            select(phenotype_fac) %>%
            arrange(phenotype_fac) %>% 
            mutate(phenotype_sel = factor(phenotype_fac, levels = c("Yes", "No"))) %>% 
            group_by(phenotype_fac) %>% 
            dplyr::summarise(n = n()) %>% 
            mutate(n_gene = sum(n),
                   prop = n/n_gene*100)%>% 
            assign("save",.,envir = .GlobalEnv),
          x = ~ phenotype_fac, 
          y = ~ round(prop, digits = 2), 
          color = ~ phenotype_fac, 
          colors = color_sel,
          type = "bar", 
          hoverinfo = "text", showlegend = FALSE,
          text= ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
    layout(title="", 
           font=plotly_font,
           xaxis = list(title="",showline = T, tickangle = 45),
           yaxis = list(title="Share of individuals (%)",showline = T),
           margin = list(b = 160)) %>%
    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  
  
}

Phenotype_num_1.fun <- function(select_gene,phenotype,colors_sel,y_axis_title,detailed_hov,scale,min_filt){
  
  if(detailed_hov == T){
  out <- plot_ly(Patient_data.df %>% 
            dplyr::rename(onset = phenotype) %>% 
            filter(onset > min_filt) %>% ##introduced a filter what small numbers are allowed 
            mutate(p_variant = paste0("p.",AA_ref,AA_pos,AA_alt)) %>% 
            filter(!is.na(onset)) %>% 
            ungroup() %>% 
            assign("save",.,envir = .GlobalEnv),
          y = ~onset, type = "box",x = "",
          boxpoints = "all", jitter = 0.3,
          pointpos = 0, hoverinfo = "text", showlegend = FALSE,
          text= ~paste0(round(onset, digits = 2), " ",scale,", ", Protein)) %>% ##sclae accounts for "months" or "years"
    layout(font=plotly_font,  
           title="",
           xaxis = list(title="", tickangle = 45, showline = T),
           yaxis = list(
             title = y_axis_title,
             tickmode = "array",
             showline = T
           )) %>%
    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  }else{ #This ensures that we have an option to remove the exact variant (in the age distribution)
  out <- plot_ly(Patient_data.df %>% 
                   dplyr::rename(Onset_days = phenotype) %>% 
                   mutate(p_variant = paste0("p.",AA_ref,AA_pos,AA_alt)) %>% 
                   filter(Onset_days > min_filt) %>% 
                   filter(!is.na(Onset_days)) %>% 
                   ungroup() %>% 
                   assign("save",.,envir = .GlobalEnv),
                 y = ~Onset_days, type = "box",x = "",
                 boxpoints = "all", jitter = 0.3,
                 pointpos = 0, hoverinfo = "text", showlegend = FALSE,
                 text= ~paste0(round(Onset_days, digits = 2), " ",scale)) %>% 
    layout(font=plotly_font,  
           title="",
           xaxis = list(title="", tickangle = 45, showline = T),
           yaxis = list(
             title = y_axis_title,
             tickmode = "array",
             showline = T
           )) %>%
    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  }
  
  return(out)
}


Onset_days_2.fun <- function(select_gene,phenotype,colors_sel){
  plot_ly(Patient_data.df %>% 
            filter(Gene == select_gene,
                   !is.na(Sz_onset)) %>% 
            rename(phenotype_fac = phenotype) %>% 
            filter(Vartype != "Synonymous",
                   Vartype != "Complex",
                   Vartype != "Missense mosaic") %>% 
            mutate(ID = paste0(phenotype_fac,ifelse(Vartype == "Missense","Missense","PTV"))) %>% 
            arrange(ID) %>% 
            mutate(ID = factor(ID, levels = unique(ID))) %>% 
            group_by(ID) %>% 
            summarise(n = n()) %>% 
            mutate(n_gene = sum(n))%>% 
            assign("save",.,envir = .GlobalEnv),
          x = ~ ID, 
          y = ~ n, 
          color = ~ ID, 
          colors = colors_sel,
          type = "bar", 
          hoverinfo = "text", showlegend = FALSE,
          text= ~ paste0(round(n/n_gene, digits = 2), " (" ,n," individuals)")) %>% 
    layout(title="", 
           font=plotly_font,
           xaxis = list(title="",tickangle=45,showline = T),
           yaxis = list(title="Number of patients",showline = T)) %>%
    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  
}

Phenotype_wordcl.fun <- function(Patient_data.df){
  
  Patient_filtered.df <- Patient_data.df %>% 
    select(26:28,31:38) 
  wc_input <- colSums(is.na(Patient_filtered.df))
  
  word_cloud_input.df<- tibble(word = names(wc_input) %>% str_replace_all(.,"_"," "),freq = wc_input) %>% 
    mutate(angle = 90 * sample(c(0, 1), nrow(.), replace = TRUE, prob = c(30, 70)))
  
  set.seed(24)
  gg_out <- ggplot(word_cloud_input.df, aes(
    label = word, size = freq,
    angle = angle
  )) +
    geom_text_wordcloud_area() +
    scale_size_area(max_size = 12) +
    theme_minimal()
 
  return(gg_out) 
}

basic_onset_legend <- function(){
  legend <- data.frame(x=c(1,3), y=c(3,3), text=c("Missense", "PTV       "))
  plot <- ggplot(legend, aes(x=x, y=y, color=text))+
    geom_point(size = 8)+
    scale_color_manual(values = c("Missense"="#BEBADA","PTV       "="#FDB462"))+
    ylim(c(2.9,3.1))+
    xlim(c(0.8,5.2))+
    theme_void()+
    geom_text(aes(label=text), hjust=-0.2, color="black", size =5)+
    theme(legend.position = "none")
  
  
  return(plot)
}

#variant analysis<

#### Variant Analysis- Variant selection ####
get_var_conseq_cDNA <- function(cDNA_input,gene_input, all_exchanges.df){
  
  cDNA_sel <- str_remove(cDNA_input,"c\\.") %>% tolower()
  
  cDNA_pos_sel <- str_extract_all(cDNA_sel,"[0-9]+") %>% unlist() %>% as.numeric()
  
  cDNA_alt_sel <- str_split(cDNA_sel,">|dup|del|ins",simplify = T) %>% .[,length(.)]
  
  cDNA_ref_sel <- case_when(str_detect(cDNA_sel,"ins") ~ str_split(cDNA_sel,"ins",simplify = T) %>% .[1] %>% str_sub(.,-1,-1),
                            str_detect(cDNA_sel,"del") ~"del",
                            str_detect(cDNA_sel,">") ~ str_split(cDNA_sel,">",simplify = T) %>% .[1] %>% str_sub(.,-1,-1),
                            str_detect(cDNA_sel,"dup") ~ str_split(cDNA_sel,"dup",simplify = T) %>% .[1] %>% str_sub(.,-1,-1),
                            TRUE~"error")
  
  variant_conseq <- case_when(length(cDNA_pos_sel) == 1 & str_detect(cDNA_sel,"ins") & nchar(cDNA_alt_sel) %%3 == 0 ~ "Inframe Insertion",
                              length(cDNA_pos_sel) == 1 & str_detect(cDNA_sel,"ins") & nchar(cDNA_alt_sel) %%3 != 0 ~ "Frameshift Insertion",
                              length(cDNA_pos_sel) == 1 & str_detect(cDNA_sel,"del") & nchar(cDNA_alt_sel) %%3 == 0 ~ "Inframe Deletion",
                              length(cDNA_pos_sel) == 1 & str_detect(cDNA_sel,"del") & nchar(cDNA_alt_sel) %%3 != 0 ~ "Frameshift Deletion",
                              length(cDNA_pos_sel) == 1 & str_detect(cDNA_sel,"dup") & nchar(cDNA_alt_sel) %%3 == 0 ~ "Inframe Duplication",
                              length(cDNA_pos_sel) == 1 & str_detect(cDNA_sel,"dup") & nchar(cDNA_alt_sel) %%3 != 0 ~ "Frameshift Duplication",
                              length(cDNA_pos_sel) == 2 & str_detect(cDNA_sel,"\\+|\\-") & cDNA_pos_sel[2] <=5~ "splice-site",
                              length(cDNA_pos_sel) == 2 & str_detect(cDNA_sel,"\\+|\\-") & cDNA_pos_sel[2] >5~ "intronic",
                              length(cDNA_pos_sel) == 1 & cDNA_ref_sel != "error"~ "Missense",
                              TRUE ~ "No valid input")
  
  variant_conseq <- ifelse(cDNA_ref_sel == "error" | 
                             (str_detect(cDNA_sel,"ins") & str_detect(cDNA_sel,"dup")) |
                             (str_detect(cDNA_sel,"ins") & str_detect(cDNA_sel,"del")) |
                             (str_detect(cDNA_sel,"dup") & str_detect(cDNA_sel,"del")), "No valid input", variant_conseq)
  
  if(variant_conseq == "Missense"){
    
    #if cDNA_pos_error is empty this cDNA_pos does not exist in this gene 
    cDNA_pos_true.df <- all_exchanges.df %>%
      filter(Gene == gene_input,
             cDNA_pos == cDNA_pos_sel)
    
    
    variant_conseq_miss <- all_exchanges.df %>%
      filter(Gene == gene_input,
             cDNA_pos == cDNA_pos_sel,
             Allele == cDNA_ref_sel %>% toupper(),
             cDNA_alt == cDNA_alt_sel %>% toupper()) %>% 
      .$Vartype
    
    
    variant_conseq <- case_when(nrow(cDNA_pos_true.df) == 0~"This is not a valid cDNA position",
                                identical(variant_conseq_miss, character(0)) & cDNA_ref_sel != cDNA_alt_sel~paste("No valid nucleotide exchange. The reference nuleotide at this position is",cDNA_pos_true.df$Allele[1]),
                                cDNA_ref_sel == cDNA_alt_sel~ "Reference and alternative nucelotide must not be the same.")
    
    variant_conseq <- ifelse(!is.na(variant_conseq), variant_conseq,variant_conseq_miss)
    
  }
  
  if(variant_conseq == "splice-site"){
    
    if(cDNA_pos_sel[1]< 6){
      
      variant_conseq <- "5'UTR variant"
      
    } else if (cDNA_pos_sel[1]> max(master.df %>% filter(Gene == gene_input) %>% .$cDNA_pos,na.rm = T)){
      
      variant_conseq <- "3'UTR variant"
    } else{
      
      splice_master.df <- master.df %>% 
        filter(Gene ==gene_input,
               cDNA_pos > cDNA_pos_sel[1]-6,
               cDNA_pos < cDNA_pos_sel[1]+6)
      
      if(max(splice_master.df$Genomic_pos)-min(splice_master.df$Genomic_pos) <12){
        
        variant_conseq <- "Selected splice variant is not close to a splice site"
        
      } 
      
    }  
    
  }
  
  return(variant_conseq)
}


#Protein


get_var_conseq_Protein <- function(protein_input,gene_input, all_exchanges.df){
  
  protein_sel <- str_remove(protein_input,"p\\.") %>% tolower()
  
  protein_pos_sel <- str_extract_all(protein_sel,"[0-9]+") %>% unlist() %>% as.numeric()
  
  
  #alternative
  protein_alt_sel <- str_split(protein_sel,"[0-9]+",simplify = T) %>% .[2]
  
  
  if(is.na(protein_alt_sel)){
    protein_alt_sel <- NA
  } else if(protein_alt_sel == "X"| protein_alt_sel == "x"|protein_alt_sel == "*"| tolower(protein_alt_sel) == "stop"){
    ##search for stop codoons
    protein_alt_sel <- "Stop"
    
  }else if (nchar(protein_alt_sel) == 3){
    
    protein_alt_sel <- paste0(str_sub(protein_alt_sel,1,1) %>% toupper(),str_sub(protein_alt_sel,2,3) %>% tolower()) %>% a() 
    
  }else{
    protein_alt_sel <- protein_alt_sel %>% toupper() 
  }
  
  protein_ref_sel <- str_split(protein_sel,"[0-9]+",simplify = T) %>% .[1]
  #reference
  
  if(is.na(protein_ref_sel)){
    protein_ref_sel <- NA
  }else if(protein_ref_sel == "X"| protein_ref_sel == "x"|protein_ref_sel == "*"| tolower(protein_ref_sel) == "stop"){
    
    protein_ref_sel <- "Stop"
    
  }else if (nchar(protein_ref_sel) == 3){
    
    protein_ref_sel <- paste0(str_sub(protein_ref_sel,1,1) %>% toupper(),str_sub(protein_ref_sel,2,3) %>% tolower()) %>% a() 
    
  }else{
    protein_ref_sel <- protein_ref_sel %>% toupper() 
  }
  
  
  #generate variant consequence 
  if(is.na(protein_ref_sel) | is.na(protein_alt_sel) | length(protein_pos_sel)>1 | length(protein_pos_sel) == 0 | (nchar(protein_ref_sel)  != 1 & protein_ref_sel != "Stop") | (nchar(protein_alt_sel)  != 1 & protein_alt_sel != "Stop")){
    
    variant_conseq <- "No valid input. Only single aminoacid exchanges are accepted."
    
  }else{
    
    protein_ref_sel <- ifelse(protein_ref_sel != "Stop", aaa(protein_ref_sel),"Stop")
    protein_alt_sel <- ifelse(protein_alt_sel != "Stop", aaa(protein_alt_sel),"Stop")
    
    #check if aminoacid position exist in this gene
    AA_pos_true.df <- all_exchanges.df %>% 
      filter(Gene == gene_input,
             AA_pos == protein_pos_sel)
    
    #check if selected reference aminoacid is correct 
    true_ref_aminoacid_check.df <- all_exchanges.df %>% 
      filter(Gene == gene_input,
             AA_pos == protein_pos_sel,
             AA_ref == protein_ref_sel)
    
    true_ref_aminoacid.df <- all_exchanges.df %>% 
      filter(Gene == gene_input,
             AA_pos == protein_pos_sel)
    
    
    variant_conseq <- case_when(nrow(AA_pos_true.df) == 0~ "Amino acid position does not exist in this gene",
                                nrow(true_ref_aminoacid_check.df) == 0~ paste0("Wrong reference amino acid. The reference amino acid is:",true_ref_aminoacid.df$AA_ref[1]),
                                protein_alt_sel == protein_ref_sel ~"Synonymous",
                                protein_alt_sel == "Stop"~ "Stop-gain",
                                protein_pos_sel ==1~ "Stop-loss",
                                (protein_alt_sel %in% (all_exchanges.df$AA_alt %>% unique() %>% .[.!= "Stop"] ))~ "Missense",
                                TRUE~"No valid alternative amino acid.")
    
    
    
  }
  
  
  return(variant_conseq)
}


extract_gnomad_features <- function(Control_data.df,selected.df,variable,label){

  if(label == "exchange"){
    
    control_int.df <- Control_data.df %>% filter(Gene == selected.df$Gene[1], AA_pos == selected.df$AA_pos[1], AA_alt == selected.df$AA_alt[1])
    
    control_int2.df <- Control_data.df %>% filter(Gene == selected.df$Gene[1], AA_pos == selected.df$AA_pos[1]) 
    
    if(variable == "Allele count"){
      
      out <- ifelse(nrow(control_int.df) == 0,0,control_int.df$Allele_count)
      
    }else if(variable == "Allele freq"){
      
      out <- ifelse(nrow(control_int.df) == 0,0,control_int.df$Allele_freq)
      
    }
  }else if(label == "position"){
    
    control_int.df <- Control_data.df %>% filter(Gene == selected.df$Gene[1], AA_pos == selected.df$AA_pos[1]) 
    
    if(variable == "Allele count"){
      
      out <- ifelse(nrow(control_int.df) == 0,0,sum(control_int.df$Allele_count))
      
    }else if(variable == "Allele freq"){
      
      out <- ifelse(nrow(control_int.df) == 0,0,sum(control_int.df$Allele_freq))
      
    }
    
  }
  
  return(out)
  
  
}

##map paraz score 
map_paraz <- function(data_scores.df){
  plot <- ggplot(data_scores.df , aes(x = AA_pos, y = Paraz_score, fill = group)) +
    geom_bar(stat = "identity") +
    theme_bw() + theme(panel.border = element_blank(),
                       legend.position = "none",
                       legend.title = element_blank()) +
    scale_fill_manual(values = c("gray","indianred"))+
    coord_cartesian(ylim = c((min(data_scores.df$Paraz_score,na.rm = T)-0.25),
                             max(data_scores.df$Paraz_score,na.rm = T)+0.25),
                    expand = FALSE) +
    labs(title="Paralog conservation",
         y="Parazscore",
         x = paste0("Amino acid sequence"))
  
  return(plot)
}


##map mtrscore 
map_mtr <- function(data_scores.df){
  
  mtr_threshold = data_scores.df$MTR_score %>% sort() %>% .[round(length(.)/4)] #lowest 25% of values. This threshold is dervied from the publication of the score 
  mtr_threshold_hard = data_scores.df$MTR_score %>% sort() %>% .[round(length(.)/20)] #lowest 5% of values. This threshold is dervied from the publication of the score 
  
  plot <- ggplot(data_scores.df, aes(x = AA_pos, y = MTR_score, colour=(MTR_score<mtr_threshold))) +
    #geom_line(aes(group=group)) +
    geom_point(aes(x = AA_pos2), color = "purple")+
    #scale_color_manual(values = c("gray","indianred1"),
    #                   labels = c("Tolerant region", "Intolerant region")) +
    geom_hline(aes(yintercept=mtr_threshold), colour="orange", linetype="dashed")  +
    geom_hline(aes(yintercept=mtr_threshold_hard), colour="indianred1", linetype="dashed")  +
    coord_cartesian(expand = FALSE) +
    theme_bw() + theme(panel.border = element_blank(),
                       legend.position = "right",
                       legend.direction = "vertical") +
    labs(title="Missense Tolerance Ratio",
         y="MTR",
         x = paste0("Amino acid sequence"),
         colour = "Patient Variants")
  
  return(plot)
}



#Research 3d mapping
map_var_3d <- function(data,Gene_sel,gnomad_bool,pdb_sel,structure_coordinates,sub_color_i){

  variant.df <- data %>%
    filter(Gene == Gene_sel) %>%
    filter(Vartype == "Missense") %>%
    mutate(label = "pathogenic") %>%
    group_by(AA_pos,AA_ref,Gene) %>%
    summarise(n_occ = n()) %>%
    select(AA_pos,AA_ref,n_occ,Gene)
  
  gnomad.df <- Control_data.df %>%
    filter(Gene == Gene_sel) %>%
    group_by(AA_pos,AA_ref,Gene) %>%
    summarise(n_occ = n()) %>%
    select(AA_pos,AA_ref,n_occ,Gene)
  
  structure.df <- read_delim(structure_coordinates,delim = "\t") %>%
    mutate(Aminoacid = aaa(Aminoacid)) %>%
    select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
  
  
  variant.df <- variant.df %>%
    left_join(structure.df,by = c("AA_pos" = "Uniprot_position","AA_ref" = "Aminoacid","Gene" = "gene")) %>%
    mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not
    filter(struc_cov == "yes")
  
  gnomad.df <- gnomad.df %>%
    left_join(structure.df,by = c("AA_pos" = "Uniprot_position","AA_ref" = "Aminoacid","Gene" = "gene")) %>%
    mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not
    filter(struc_cov == "yes")
  
  sub_color <- c("#e2f970","#6fbbf7","#ee6c71","#ffbc5a","#bf73cc")
  sub_scale <- c(1.2,0.8)
  struc_color <- "wheat"
  
  rot = 270
  rot_axis = "x"
  spin_axis = "vy"
  
  #Specify yourself- color of the cartoon per subunit
  subunit_color <- c("wheat","white") #
  
  print(variant.df)
  print("now gnomad")
  print(gnomad.df)
  #Model for the protein complex
  
  modelo <- r3dmol(
    viewer_spec = m_viewer_spec(
      cartoonQuality = 10,
      lowerZoomLimit = 50,
      upperZoomLimit = 1000
    )
  )
  
  modelo <- modelo %>% m_add_model(data = pdb_sel, format = "pdb")
  
  # Zoom to encompass the whole scene
  modelo <- modelo %>% m_zoom_to() %>%
    # Set color o cartoon representation
    m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
    # Set subunit colors
    m_set_style(
      sel = m_sel(chain = c("A")),
      style = m_style_cartoon(color = subunit_color[2])
    ) %>%
    # visualize variants grin1
    m_set_style(
      sel = m_sel(resi = variant.df$Position_in_structure,
                  atom = "CA",
                  chain = c("A")),
      style = m_style_sphere(colorScheme = NULL,
                             color = sub_color[sub_color_i],
                             scale = sub_scale[1])
    )
  
  if  (gnomad_bool == TRUE) {
    
    modelo <- modelo %>% m_set_style(
      sel = m_sel(resi = gnomad.df$Position_in_structure,
                  atom = "CA",
                  chain = c("A")),
      style = m_style_sphere(colorScheme = NULL,
                             color = "#333333",
                             scale = sub_scale[2]))
    
  }

  return(modelo) 
}     

## substitute for seqinr function
three_to_one_aa <- function(sequence){
  
  sapply(sequence, function(x){
    
    ifelse(is.na(x),NA,seqinr::aaa(x))
    
  }) %>% 
    as.vector() %>% 
    return()
  
}
############## DATA ############
#Load domain data 
Domain_data.df <- read_delim("data/Domain.txt", delim = "\t") %>% 
  select(Aln_pos,Domain,Domain_color)

#Load all possible exchanges 
all_exchanges.df <- read_delim("data/master_table_exchanges.txt",delim = "\t") %>% 
  left_join(Domain_data.df %>% distinct(Aln_pos,Domain,Domain_color)) %>% 
  mutate(AA_ref = aaa(AA_ref),
         AA_alt = convert_aa(AA_alt),
         Vartype = ifelse(Vartype == "Stop Gained","Nonsense", Vartype))

#Load master table 
master.df <- read_delim("data/master_table.txt", delim = "\t") 

#load domains
Domain_gene.df <- master.df %>% 
  left_join(Domain_data.df %>% distinct(Aln_pos,Domain,Domain_color)) %>% 
  distinct(Gene,AA_pos,Domain,Domain_color) %>% 
  ungroup()


# #PER-2D
# per_family2d.df <- read_delim("data/per2d_genewise.txt", delim = "\t")
# 
# per2d_for_var_analysis.df <- read_delim("data/per2d_genewise.txt", delim = "\t")
# 
# per2d.df <- read_delim("data/per2d_genewise.txt", delim = "\t") %>% 
#   mutate(Hotzone_2D = ifelse(per == "PER","PER","No-PER")) %>% 
#   select(Hotzone_2D,AA_pos,pvalue,odds) %>% 
#   rename(pvalue_per3d = "pvalue",
#          odds_per3d = "odds")
# #PER-3D
# per3d.df <- read_delim("data/per3d.txt", delim = "\t") %>% 
#   select(PER3D,AA_pos,pvalue,odds) %>% 
#   rename(pvalue_per2d  = "pvalue",
#          odds_per2d  = "odds")

# #Load functional data
# Functional_data.df <- read_delim("data/Functional_data.txt", delim = "\t") %>% 
#   rename(uptake = "GABA uptake (vs wt)",
#          surface_exp = "Surface expression (vs wt)",
#          total_exp = "Total expression (vs wt)",
#          relative_update_surface_exp = "Relative uptake/surface expression",
#          relative_surface_exp_tot_exp = "Relative surface expression/total expression") %>% 
#   filter(!is.na(AA_pos),
#          AA_alt != "X") %>% 
#   select(AA_pos,AA_alt,uptake,surface_exp,total_exp,relative_update_surface_exp,relative_surface_exp_tot_exp)

#Load patient and control data 
Patient_data.df <- read_delim("data/SATB2_Patient_variants_v1.txt", delim = "\t") %>% 
  select(-Transcript) %>% 
  mutate(AA_pos = as.numeric(AA_pos)) %>% 
  #this ensures that splice cite mutation are shown in the genotype overview plot
  mutate(cDNA_pos = str_extract(Original_cDNA_change, "[0-9]+") %>% as.numeric()) %>% 
  left_join(master.df %>% 
              select(AA_pos,cDNA_pos) %>% 
              dplyr::rename(AA_pos_splice_site = "AA_pos")) %>% 
  mutate(AA_pos = ifelse(Vartype == "splice site",AA_pos_splice_site,AA_pos),
         AA_alt = ifelse(!is.na(AA_alt) & AA_alt == "STOP","Stop",
                         ifelse(!is.na(AA_alt),AA_alt,NA))) %>% 
  ##
  left_join(master.df %>% distinct(Transcript,Gene,AA_pos), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) %>% 
  left_join(Domain_gene.df %>% distinct(Domain,Gene,AA_pos,Domain_color), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) %>%  
  mutate(AA_ref = ifelse(!is.na(AA_ref),AA_ref,"XXX"), #%>% aaa(), ##warnings due to none matching aminoacids are fine 
         #AA_alt_complex = ifelse(Vartype == "Missense",AA_alt,AA_alt_complex),
         cDNA = ifelse(!is.na(cDNA_pos), paste0("c.",cDNA_pos,cDNA_ref,">",cDNA_alt), "Not available"),
         #Protein = paste0("p.",AA_ref,AA_pos,AA_alt_complex),
         Protein = Original_AA_change)

Patient_data_missense_only.df <- Patient_data.df %>% 
  filter(Vartype == "Missense") #%>% 
 # mutate(AA_alt = three_to_one_aa(AA_alt))

Control_data.df <- read_delim("data/gnomad_variants.txt", delim = "\t") %>% 
  mutate(AA_ref = aaa(AA_ref)) %>% #,
         #AA_alt = convert_aa(AA_alt)) %>%
  # left_join(per3d.df) %>% 
  # left_join(per2d.df) %>% 
  left_join(Domain_gene.df %>% distinct(Domain,Gene,AA_pos,Domain_color), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) 

#Load Scores
#Paraz/MTR
# paraz_mtr.df <- read_delim("data/mtr_paraz_slc6a1.txt",delim = "\t")

#hotzones3D on structure 
#PER3D_struc.df <- read_delim("data/pdb/6j8e_varburden.txt", delim = "\t")

#Load ClinVar data for Variant Analysis
clinvar.df <- read_delim("data/Clinvar_links_SATB2.txt", delim = "\t")

##### Variables #####
  #Basic Information
  basic_gene1 = "SATB2"
  basic_phenotype_fac = "Epilepsy"
  basic_phenotype_num = "Sz_onset"
  
  
  phenotype_name1 <- "Epilepsy"
  phenotype_name2 <- "Seizure onset (months)"
  phenotype_name3 <- "Seizure onset (months)"
  
  basic_phenotype_colors <- RColorBrewer::brewer.pal(20,"Set3") # Warning occurs at the moment, can be ignored for now 
  
  ##add here the colors for all phenotypes 
  basic_phenotype_colors <- c("Female"="#F97EE5","Male"="#0404B1",
                                  " Yes" ="#F97EE5", "No"= "#0404B1", "Yes" = "#F97EE5",
                                  "None" = "#9191fd" ,"1 to 10"= "#4747fb" ,"10 to 50" = "#0404b1","Greater than 50" ="#02026a")
  
  basic_phenotype_colors_ID <- c("black","grey","#994714")
  
  basic_phenotype_colors_yes_no <-  c("No" = "grey",
                                      " Yes" = "#EFC56F", 
                                      "Yes" = "#EFC56F")
  basic_phenotype_colors_speech <- c("10 to 50" = "grey40",
                                  "Greater than 50" = "grey20",
                                  "None" = "grey80",
                                  "1 to 10" = "grey60")
  colors_gene1 <-  c("#BEBADA","#FDB462","#BEBADA","#FDB462")
 
##### Variant Analysis variable #####
  variant_title1 <- "Epilepsy"
  variant_title2 <- "Autism"
  variant_title3 <- "Seizure onset"

#####Research variable #####
  Gene_colors <-  c("Patient"="red", 
                    "Control" = "blue",
                    "Other" = "#333333")
  
  research_phenotype1_title <- "Number of patient variants per unit"
  research_phenotype2_title <- "Epilepsy"
  research_phenotype3_title <- "Seizure onset"
  research_phenotype4_title <- "Autism"
  research_phenotype5_title <- "Cognitive development"
  research_phenotype6_title <- "Epilepsy syndrome classification"
  
  research_functional1_title <- "Variants with molecular function assessment"
  #3d_mapping Genotype interface 
  pdb_sel_gene1 = "data/pdb/SCN1A_model.pdb1"
  structure_coordinates_gene1 <- "data/pdb/7dtd_structure_coordinates.txt"
  
  pdb_sel_gene2 = "data/pdb/SCN2A_model.pdb1"
  structure_coordinates_gene2 <- "data/pdb/6j8e_structure_coordinates.txt"
  
  pdb_sel_gene3 = "data/pdb/SCN2A_model.pdb1"
  structure_coordinates_gene3 <- "data/pdb/SCN3A_6j8e_structure_coordinates.txt"
  
  pdb_sel_gene4 = "data/pdb/SCN2A_model.pdb1"
  structure_coordinates_gene4 <- "data/pdb/SCN8A_6j8e_structure_coordinates.txt"

############### SERVER ###############
shinyServer(function(input, output, session) {
  
  # output$disclaimer <- renderValueBox({
  #   valueBox(
  #     value = tags$p("Pre-release version",style = "font-size:60%"),
  #     div("Website only online during ASHG 2021. New final version will be released during AES 2021 in December."),color = "red"
  #   )
  #   
  # })

  ############### PANEL DECISION Welcome page###############
  observeEvent(input$infoBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "infoTab")
  })
  
  observeEvent(input$familyBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "familyTab")
  })
  
  observeEvent(input$variantBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "variantTab")
  })
  
  observeEvent(input$researchBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "researchTab")
  })
  
  observeEvent(input$registryBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "registryTab")
  })
  
  ##### BASIC INFORMATION #####
  #can be replicated for each additional gene 
  #Factor Gene1 
  
  output$Phenotype_fac1_1 <- renderPlotly({
    Phenotype_fac_2.fun(basic_gene1, "Sex",basic_phenotype_colors)
  })
  
  output$Phenotype_fac3_1 <- renderPlotly({
    Phenotype_fac_1.fun("Total_speech",basic_phenotype_colors)
  })
  
  output$Phenotype_fac4_1 <- renderPlotly({
    Phenotype_fac_3.fun("Cleft_palate",basic_phenotype_colors)
  })
  
  output$Phenotype_fac4_2 <- renderPlotly({
    Phenotype_fac_3.fun("Low_BMD",basic_phenotype_colors)
  })
  
  output$Phenotype_fac5 <- renderPlotly({
    Phenotype_wordcl.fun(Patient_data.df)
  })
  
  
  
  # Numeric Gene 1
  
  output$Phenotype_num1_2 <- renderPlotly({
    Phenotype_num_1.fun(basic_gene1, "Age_years", colors_gene1, "Age (years)",F,"years",0)
  })
  
  output$Phenotype_num2_1 <- renderPlotly({
    Phenotype_num_1.fun(basic_gene1, "Age_first_word_months", colors_gene1, "Age (months)",T,"months",5)
  })
  
  output$Phenotype_num2_2 <- renderPlotly({
    Phenotype_num_1.fun(basic_gene1, "Age_walk_months", colors_gene1, "Age (months)",T, "months",5)
  })
  

  output$basic_legend_plot1 <- renderPlot({
    basic_onset_legend()
  })
  
  
  ##### FOR FAMILIES #####
  
  # nothing to be calculated here # 
  
  ##### VARIANT ANALYSIS ##### 
  
  #Updates after gene change 
  
  output$value <- renderText({ input$cDNA })
  
  output$vartype <- renderText({ 6 })
  ##Updates after Gene selection
  observeEvent("SATB2", {
    
    variant_conseq <- get_var_conseq_cDNA(input$cDNA,"SATB2",all_exchanges.df)
    
    output$value <- renderText({ variant_conseq })
    
  })
  
  
  #Updates after cDNA change 
  observeEvent(input$cDNA, {
    
    variant_conseq <- get_var_conseq_cDNA(input$cDNA,"SATB2",all_exchanges.df)
    
    output$value <- renderText({ variant_conseq })
    
    output$consequenceBox <- renderValueBox({
      valueBox(
        value = tags$p("Variant consequence",style = "font-size:70%"), 
        subtitle = tags$p(variant_conseq,style = "font-size:150%;font-color:black"),
        color = "olive"
      )
    })
    
  })
  
  
  observeEvent(input$Protein, {
    
    variant_conseq <- get_var_conseq_Protein(input$Protein,"SATB2",all_exchanges.df)
    
    output$value <- renderText({ variant_conseq })
    
    output$consequenceBox <- renderValueBox({
      valueBox(
        value = tags$p("Variant consequence",style = "font-size:70%"), 
        subtitle = tags$p(variant_conseq,style = "font-size:150%;font-color:black"),
        color = "olive"
      )
    })
    
  })
  
  ## Once the search button is pressed 
  
  varFilterInput <- reactiveValues(data=NULL)
  
  varFilterInputClinvar <- reactiveValues(data=NULL)
  
  observeEvent(input$search_var_c, {
    
    variant_conseq <- get_var_conseq_cDNA(input$cDNA,"SATB2",all_exchanges.df)
    
    output$value <- renderText({ variant_conseq })
    
    output$vartype <- reactive({
      g <- ifelse(variant_conseq == "Missense","Missense",
             ifelse(variant_conseq %in% c("Frameshift Insertion", "Frameshift Deletion","Frameshift Duplication", "splice-site", "Nonsense"),"Nonsense","Other"))
      
      print(g)
    })
    
    outputOptions(output, 'vartype', suspendWhenHidden = FALSE)
    
    output$consequenceBox <- renderValueBox({
      valueBox(
        value = tags$p("Variant consequence",style = "font-size:70%"), 
        subtitle = tags$p(variant_conseq,style = "font-size:150%;font-color:black"),
        color = "olive"
      )
    })
    #generate datafile for subsequent visualz 
    
    if(variant_conseq == "Missense"){
      
      cDNA_sel <- str_remove(input$cDNA,"c\\.") %>% tolower()
      
      cDNA_pos_sel <- str_extract_all(cDNA_sel,"[0-9]+") %>% unlist() %>% as.numeric()
      
      cDNA_ref_sel <- str_split(cDNA_sel,">",simplify = T) %>% .[1] %>% str_sub(.,-1,-1) %>% toupper()
      
      cDNA_alt_sel <- str_split(cDNA_sel,">",simplify = T) %>% .[2] %>% str_sub(.,-1,-1) %>% toupper()
      
      
      varFilterInput$data <- all_exchanges.df %>% 
        filter(Gene=="SATB2", 
               cDNA_pos == cDNA_pos_sel,
               Allele == cDNA_ref_sel,
               cDNA_alt == cDNA_alt_sel)
      
      #For clinvar significance
      varFilterInputClinvar$data <- clinvar.df %>% 
        filter(Gene=="SATB2", 
               AA_pos == varFilterInput$data$AA_pos[1],
               AA_ref == varFilterInput$data$AA_ref[1],
               AA_alt == varFilterInput$data$AA_alt[1])
    }else{
      
      print("test")
      
      varFilterInput$data <-tibble(cDNA_input = input$cDNA, Gene = "SATB2")
      
    }
    
  })
  
  observeEvent(input$search_var_p, {
    
    variant_conseq <- get_var_conseq_Protein(input$Protein,"SATB2",all_exchanges.df)
    
    output$value <- renderText({ variant_conseq })
    
    output$vartype <- renderText({
      ifelse(variant_conseq == "Missense","Missense",
             ifelse(variant_conseq %in% c("Frameshift Insertion", "Frameshift Deletion","Frameshift Duplication", "splice-site", "Stop-gain","Nonsense"),"Nonsense","Other"))
    })
    
    output$consequenceBox <- renderValueBox({
      valueBox(
        value = tags$p("Variant consequence",style = "font-size:70%"), 
        subtitle = tags$p(variant_conseq,style = "font-size:150%;font-color:black"),
        color = "olive"
      )
    })
    
    #generate datafile for subsequent visualz 
    
    if(variant_conseq == "Missense"){
      
      protein_sel <- str_remove(input$Protein,"p\\.") %>% tolower()
      
      protein_pos_sel <- str_extract_all(protein_sel,"[0-9]+") %>% unlist() %>% as.numeric()
      
      protein_ref_sel <- str_split(protein_sel,"[0-9]+",simplify = T) %>% .[1]
      protein_ref_sel <- ifelse(nchar(protein_ref_sel) ==3, 
                                paste0(str_sub(protein_ref_sel,1,1) %>% toupper(),str_sub(protein_ref_sel,2,3) %>% tolower()),
                                aaa(toupper(protein_ref_sel)))
      
      protein_alt_sel <- str_split(protein_sel,"[0-9]+",simplify = T) %>% .[2]
      protein_alt_sel <- ifelse(nchar(protein_alt_sel) ==3, 
                                paste0(str_sub(protein_alt_sel,1,1) %>% toupper(),str_sub(protein_alt_sel,2,3) %>% tolower()),
                                aaa(toupper(protein_alt_sel)))
      
      varFilterInput$data <- all_exchanges.df %>% 
        filter(Gene=="SATB2", 
               AA_pos == protein_pos_sel,
               AA_ref == protein_ref_sel,
               AA_alt == protein_alt_sel) 
      
      #For clinvar significance
      varFilterInputClinvar$data <- clinvar.df %>% 
        filter(Gene=="SATB2", 
               AA_pos == varFilterInput$data$AA_pos[1],
               AA_ref == varFilterInput$data$AA_ref[1],
               AA_alt == varFilterInput$data$AA_alt[1])
      
    }else{
      varFilterInput$data <- tibble(input = input$Protein)
    }
    
  })
  
  #Updates after gene change 
  
  # observeEvent("SATB2", {
  # 
  #   filtered_cDNA_pos <- all_exchanges.df %>% filter(cDNA_pos == input$search_cDNA_pos, Gene == "SATB2")
  #   
  # 
  #   updatePickerInput(session = session, inputId = "search_Allele",
  #                     choices = unique(filtered_cDNA_pos$Allele))
  # 
  #   updatePickerInput(session = session, inputId = "search_cDNA_alt",
  #                     choices = if(all(is.na(filtered_cDNA_pos$cDNA_alt))){NA}else{unique(filtered_cDNA_pos$cDNA_alt[!is.na(filtered_cDNA_pos$cDNA_alt)])},
  #                     selected = filtered_cDNA_pos$cDNA_alt[1])
  # 
  #   updatePickerInput(session = session, inputId = "search_AA_alt",
  #                     choices = all_exchanges.df %>% filter(Gene == "SATB2", AA_pos == filtered_cDNA_pos$AA_pos[1]) %>% .$AA_alt %>% unique(),
  #                     selected = filtered_cDNA_pos %>% filter(cDNA_alt == cDNA_alt[1]) %>% .$AA_alt)
  # 
  #   updatePickerInput(session = session, inputId = "search_AA_ref",
  #                     choices = filtered_cDNA_pos %>% filter(cDNA_alt == cDNA_alt[1]) %>% .$AA_ref,
  #                     selected = filtered_cDNA_pos %>% filter(cDNA_alt == cDNA_alt[1]) %>% .$AA_ref)
  # })
  #   
  # #Updates after DNA change 
  # 
  # observeEvent(input$search_cDNA_pos, {
  # 
  #   filtered_cDNA_pos <- all_exchanges.df %>% filter(cDNA_pos == input$search_cDNA_pos, Gene == "SATB2")
  #   filtered_cDNA_pos_aa <- all_exchanges.df %>% filter(AA_pos == filtered_cDNA_pos$AA_pos, Gene == "SATB2")
  # 
  #   updatePickerInput(session = session, inputId = "search_Allele",
  #                     choices = unique(filtered_cDNA_pos$Allele))
  # 
  #   updatePickerInput(session = session, inputId = "search_cDNA_alt",
  #                     choices = if(all(is.na(filtered_cDNA_pos$cDNA_alt))){NA}else{unique(filtered_cDNA_pos$cDNA_alt[!is.na(filtered_cDNA_pos$cDNA_alt)])},
  #                     selected = "A")
  # 
  #   updateNumericInputIcon(session = session, inputId = "search_AA_pos",
  #                     value = if(all(is.na(filtered_cDNA_pos$AA_pos))){NA}else{unique(filtered_cDNA_pos$AA_pos[!is.na(filtered_cDNA_pos$AA_pos)])})
  # 
  #   updatePickerInput(session = session, inputId = "search_AA_ref",
  #                     choices = unique(filtered_cDNA_pos_aa$AA_ref))
  # 
  #   updatePickerInput(session = session, inputId = "search_AA_alt",
  #                     choices = unique(filtered_cDNA_pos_aa$AA_alt))
  # 
  # })
  # 
  # observeEvent(input$search_cDNA_alt, {
  # 
  #   filtered_cDNA_pos <- all_exchanges.df %>% filter(cDNA_pos == input$search_cDNA_pos, Gene == "SATB2", cDNA_alt == input$search_cDNA_alt)
  # 
  #   updatePickerInput(session = session, inputId = "search_AA_alt",
  #                     selected = unique(filtered_cDNA_pos$AA_alt))
  # 
  # })
  # 
  # #Updates after AA change
  # observeEvent(input$search_AA_alt, {
  #   filtered_AA_alt <- all_exchanges.df %>% filter(AA_pos == input$search_AA_pos, Gene == "SATB2", AA_alt == input$search_AA_alt)
  # 
  #   updatePickerInput(session = session,
  #                     inputId = c("get_var_type"),
  #                     choices = c(unique(filtered_AA_alt$Vartype)),
  #                     selected = c(unique(filtered_AA_alt$Vartype)))
  # })
  # 
  # observeEvent(input$search_AA_pos, {
  # 
  #   filtered_AA_pos <- all_exchanges.df %>% filter(AA_pos == input$search_AA_pos, Gene == "SATB2")
  # 
  #   updatePickerInput(session = session, inputId = "search_AA_ref",
  #                     choices = sort(unique(filtered_AA_pos$AA_ref)))
  # 
  #   updatePickerInput(session = session, inputId = "search_AA_alt",
  #                     choices = unique(filtered_AA_pos$AA_alt))
  #   
  # })
  # 
  # ## Actions once the search button is pressed 
  # 
  # varFilterInput <- reactiveValues(data=NULL)
  # 
  # varFilterInputClinvar <- reactiveValues(data=NULL)
  # 
  # observeEvent(input$search_var_c, {
  #   varFilterInput$data <- all_exchanges.df %>% filter(Gene=="SATB2") %>% filter(cDNA_pos==input$search_cDNA_pos) %>%
  #     filter(Allele==input$search_Allele) %>% filter(cDNA_alt==input$search_cDNA_alt)
  #   
  #   varFilterInputClinvar$data <- clinvar.df %>% filter(Gene=="SATB2") %>% filter(AA_pos==input$search_AA_pos) %>%
  #     filter(AA_ref==input$search_AA_ref) %>% filter(AA_alt==input$search_AA_alt)
  #   
  # })
  # 
  # observeEvent(input$search_var_p, {
  #   varFilterInput$data <- all_exchanges.df %>% filter(Gene=="SATB2") %>% filter(AA_pos==input$search_AA_pos) %>%
  #     filter(AA_ref==input$search_AA_ref) %>% filter(AA_alt==input$search_AA_alt)
  #   
  #   varFilterInputClinvar$data <- clinvar.df %>% filter(Gene=="SATB2") %>% filter(AA_pos==input$search_AA_pos) %>%
  #     filter(AA_ref==input$search_AA_ref) %>% filter(AA_alt==input$search_AA_alt)
  #   
  # })
  # 
  # # update cDNA input according to user search input , not 100% precise as many option may be possible 
  # 
  #  
  # observeEvent(input$search_var_p, {
  #   updateNumericInputIcon(
  #     session = session,
  #     inputId = c("search_cDNA_pos"),
  #     value = c(unique(varFilterInput$data$cDNA_pos)[1])
  #   )
  # 
  #   updatePickerInput(
  #     session = session,
  #     inputId = c("search_Allele"),
  #     choices = c(unique(varFilterInput$data$Allele)[1])
  #   )
  # 
  #   updatePickerInput(
  #     session = session,
  #     inputId = c("search_cDNA_alt"),
  #     selected = c(unique(varFilterInput$data$cDNA_alt))
  #                  
  #   )
  #    
  # })
  
  
  ##ClinVar variant interpretation
  output$ClinVarbox <- renderValueBox({
    
    clinvar_input <- ifelse(nrow(varFilterInputClinvar$data)>0,varFilterInputClinvar$data$ClinicalSignificance[1],"No ClinVar variants available")
    
    clinvar_link <-  ifelse(nrow(varFilterInputClinvar$data)>0,varFilterInputClinvar$data$Link[1],"https://www.ncbi.nlm.nih.gov/clinvar/")
    
    valueBox(
      value = tags$p(clinvar_input,style = "font-size:60%"),
      div(HTML(paste0("To access variant click <a style=color:black;  href=\"",clinvar_link,"\">here</a>"))),
      # div("To access variant click ", shiny::a("here", 
      #                                 href=clinvar_link, 
      #                                 target="_blank",
      #                                 <a style=color:white;
      #                                 )),
      color = "red"
    )
  })
  
##### Variant Information #####
  
  output$geneBox1 <- renderValueBox({
    
    valueBox(
      value = tags$p(paste0("Gene: ",unique(varFilterInput$data$Gene)),style = "font-size:50%"), 
      div("Transcript: ",unique(varFilterInput$data$Transcript),br(),br(),""), icon = icon("dna"),
      color = "purple"
    )
  })
  
  output$geneBox2 <- renderValueBox({
    
    if(!is.null(varFilterInput$data)){
      
      p_old <- unique(varFilterInput$data$AA_ref)  
      p_old_aa <- seqinr::a(p_old)
      p_new <- unique(varFilterInput$data$AA_alt) 
      p_new_aa <- ifelse(p_new != "Stop",seqinr::a(p_new),"*")
      p_domain <- Domain_gene.df %>% filter(Gene == varFilterInput$data$Gene[1], AA_pos == varFilterInput$data$AA_pos[1]) %>% .$Domain %>% unique()
      p_pos <- unique(varFilterInput$data$AA_pos)
      
    }else{
      
      p_old <- ""
      p_old_aa <- ""
      p_new <- ""
      p_new_aa <- ""
      p_domain <- ""
      p_pos <- ""
      
    }
    
    valueBox(
      value = tags$p(paste0("Domain: ", p_domain),style = "font-size:50%"),  
      div("Amino Acid Position: ",p_pos, br(),
          paste0("Amino Acid Change: ", p_old, " (",p_old_aa, ") "), "-" ,
          paste0(p_new, " (",ifelse(p_new != "Stop",p_new_aa,"*"),")")), icon = icon("dna"),
      color = "light-blue"
    )
  })
  
  output$geneBox3 <- renderValueBox({
    
    if(is.null(varFilterInput$data)){
      
      valueBox(
        value = tags$p(paste0("Control variants: "),style = "font-size:50%"),  
        div(paste0("gnomAD allele count at same position: "),br(),
                                                                  paste0("gnomAD allele frequency at same substitution: ")), icon = icon("dna"),
        color = "green"
      )
      
    }else{
    
      gnomad_count_exchange <- round(extract_gnomad_features(Control_data.df,varFilterInput$data,"Allele count","exchange"),2)
      
      gnomad_freq_exchange <- extract_gnomad_features(Control_data.df,varFilterInput$data,"Allele freq","exchange") 
      
      if(gnomad_freq_exchange !=0){
         gnomad_freq_exchange <- paste0(str_sub(gnomad_freq_exchange,1,3),str_sub(gnomad_freq_exchange,-3,-1))
      }
      
      gnomad_count_position <- extract_gnomad_features(Control_data.df,varFilterInput$data,"Allele count","position")
      
      valueBox(
        value = tags$p(paste0("Control variants: ",  gnomad_count_exchange),style = "font-size:50%"),  
        div(paste0("gnomAD allele count at same position: ", gnomad_count_position),br(),
                                                                  paste0("gnomAD allele frequency at same substitution ", gnomad_freq_exchange)), icon = icon("dna"),
        color = "green"
      )
    }
  })
  
  # Table with patient information ####
  output$patientTable <- DT::renderDataTable({
    
    print(varFilterInput$data)
    
    validate(
      need(!plyr::empty(varFilterInput$data),
           "There is no data that matches your filters.")) 
   
    datatable(Patient_data.df %>% 
                filter(AA_pos==varFilterInput$data$AA_pos[1]) %>%
                filter(AA_alt == varFilterInput$data$AA_alt[1]) %>% 
                       #!is.na(Published_in)) %>%
      select(Domain, Original_cDNA_change, Original_AA_change, Origin,Cleft_palate,Low_BMD,Abnormal_brainMRI,Age_walk_months,Age_first_word_months,Total_speech,Dental_issues,Clinical_seizures,Behavior_anomalies,Sleep_problems,Published_in),
    colnames = c("Domain","cDNA level","Protein level","Origin","CP","Low BMD","Abnl MRI","Walk at (months)", "Talk at (months)","Speech (words)","Abnl Teeth","Seizures", "Abnl behaviour","Abnl sleep","Link"),
    options = list(dom = 't', scrollY = TRUE), escape=FALSE)
  })
  
  output$patientTable_nonsense <- DT::renderDataTable({
    
    #print(varFilterInput$data)
    
    validate(
      need(!plyr::empty(varFilterInput$data),
           "There is no data that matches your filters.")) 

    
    if(!is.null(varFilterInput$data$input)){
    
      AA_ref_sel <- str_sub(varFilterInput$data$input,3,5)
      AA_pos_sel <- str_extract(varFilterInput$data$input,"[0-9]+")
      
      table_filt.df <- Patient_data.df %>% 
        filter(AA_pos == AA_pos_sel,
               AA_ref == AA_ref_sel,
               AA_alt == "Stop") 
    }else{
      
      table_filt.df <- Patient_data.df %>% 
        filter(Original_cDNA_change == varFilterInput$data$cDNA_input)
      
    }
    
    datatable(table_filt.df %>% 
                #!is.na(Published_in)) %>%
                select(Domain, Original_cDNA_change, Original_AA_change, Origin,Cleft_palate,Low_BMD,Abnormal_brainMRI,Age_walk_months,Age_first_word_months,Total_speech,Dental_issues,Clinical_seizures,Behavior_anomalies,Sleep_problems,Published_in),
              colnames = c("Domain","cDNA level","Protein level","Origin","CP","Low BMD","Abnl MRI","Walk at (months)", "Talk at (months)","Speech (words)","Abnl Teeth","Seizures", "Abnl behaviour","Abnl sleep","Link"),
              options = list(dom = 't', scrollY = TRUE), escape=FALSE)
  })
  
  observe_helpers(withMathJax = TRUE)
  
  output$comparePlot <- renderPlotly({
    
    validate(need(
      !plyr::empty(varFilterInput$data),
      "There is no data that matches your filters."
    ))
    

    z <- Patient_data.df %>%
      filter(case_when(
        input$compareButtons =="Variant Type" ~ Vartype ==varFilterInput$data$Vartype,
        input$compareButtons =="Protein region" ~ Domain==varFilterInput$data$Domain,
        #input$compareButtons == "Functional Consequence" ~ CALL == varFilterInput$data$CALL,
        input$compareButtons =="Amino Acid Position" ~ AA_pos==varFilterInput$data$AA_pos)) 

    
    validate(need(
      !plyr::empty(z),
      "There is no data that matches your filters."
    ))

    selected_data.df <- z %>% 
      dplyr::rename(phenotype_sel = Cleft_palate) %>% 
      select(phenotype_sel) %>%
      arrange(phenotype_sel) %>% 
      filter(!is.na(phenotype_sel)) %>% 
      group_by(phenotype_sel) %>% 
      dplyr::summarise(n = n()) %>% 
      assign("save",.,envir = .GlobalEnv)
    
    all_data.df <- Patient_data.df %>% 
      dplyr::rename(phenotype_sel = Cleft_palate) %>% 
      select(phenotype_sel) %>%
      arrange(phenotype_sel) %>% 
      filter(!is.na(phenotype_sel)) %>% 
      group_by(phenotype_sel) %>% 
      dplyr::summarise(n = n()) 
    
    n_tot_y <- ifelse(length(selected_data.df$n[which(selected_data.df == "Yes")]) != 0,selected_data.df$n[which(selected_data.df == "Yes")],0)
    
    n_tot_n <- ifelse(length(selected_data.df$n[which(selected_data.df == "No")]) != 0,selected_data.df$n[which(selected_data.df == "No")],0)
    
    plot1_input.df <- tibble(phenotype_sel = c("Yes","No"), n_tot = c(n_tot_y,n_tot_n), per = c(n_tot_y/all_data.df$n[which(all_data.df$phenotype_sel == "Yes")]*100,n_tot_n/all_data.df$n[which(all_data.df$phenotype_sel == "No")]*100))
    
    plotty1 <- plot_ly(data = plot1_input.df, 
                       x = ~ phenotype_sel, 
                       y = ~ round(per, digits = 2), 
                       color = ~ phenotype_sel, 
                       colors = basic_phenotype_colors,
                       type = "bar", 
                       hoverinfo = "text", showlegend = FALSE,
                       text= ~ paste0(round(n_tot, digits = 2), " (" ,n_tot," individuals)")) %>% 
      layout(title="", 
             font=plotly_font,
             xaxis = list(title="",showline = T, tickangle = 45),
             yaxis = list(title="Share of individuals (%)",showline = T),
             margin = list(b = 160)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  
                       
    
    selected_data.df <- z %>% 
      dplyr::rename(phenotype_sel = Total_speech) %>% 
      select(phenotype_sel) %>%
      arrange(phenotype_sel) %>% 
      filter(!is.na(phenotype_sel)) %>% 
      mutate(phenotype_sel = factor(phenotype_sel, levels = c("None","1 to 10","10 to 50","Greater than 50"))) %>% 
      group_by(phenotype_sel) %>% 
      dplyr::summarise(n = n()) %>% 
      assign("save",.,envir = .GlobalEnv)
    
    all_data.df <- Patient_data.df %>% 
      dplyr::rename(phenotype_sel = Total_speech) %>% 
      select(phenotype_sel) %>%
      arrange(phenotype_sel) %>% 
      filter(!is.na(phenotype_sel)) %>% 
      mutate(phenotype_sel = factor(phenotype_sel, levels = c("None","1 to 10","10 to 50","Greater than 50"))) %>% 
      group_by(phenotype_sel) %>% 
      dplyr::summarise(n = n())
    
    n_tot_none <- ifelse(length(selected_data.df$n[which(selected_data.df == "None")]) != 0,selected_data.df$n[which(selected_data.df == "None")],0)
    
    n_tot_1_10 <- ifelse(length(selected_data.df$n[which(selected_data.df == "1 to 10")]) != 0,selected_data.df$n[which(selected_data.df == "1 to 10")],0)
    
    n_tot_10_50 <- ifelse(length(selected_data.df$n[which(selected_data.df == "10 to 50")]) != 0,selected_data.df$n[which(selected_data.df == "10 to 50")],0)
    
    n_tot_g50 <- ifelse(length(selected_data.df$n[which(selected_data.df == "Greater than 50")]) != 0,selected_data.df$n[which(selected_data.df == "Greater than 50")],0)
    
    plot2_input.df <- tibble(phenotype_sel = c(" None","1 to 10","10 to 50","Greater than 50"), 
                             n_tot = c(n_tot_none,n_tot_1_10,n_tot_10_50,n_tot_g50), 
                             per = c(n_tot_none/all_data.df$n[which(all_data.df$phenotype_sel == "None")]*100,
                                     n_tot_1_10/all_data.df$n[which(all_data.df$phenotype_sel == "1 to 10")]*100,
                                     n_tot_10_50/all_data.df$n[which(all_data.df$phenotype_sel == "10 to 50")]*100,
                                     n_tot_g50/all_data.df$n[which(all_data.df$phenotype_sel == "Greater than 50")]*100))
    
    plotty2 <- plot_ly(data = plot2_input.df, 
                       x = ~ phenotype_sel, 
                       y = ~ round(per, digits = 2), 
                       color = ~ phenotype_sel, 
                       colors = basic_phenotype_colors,
                       type = "bar", 
                       hoverinfo = "text", showlegend = FALSE,
                       text= ~ paste0(round(n_tot, digits = 2), " (" ,n_tot," individuals)")) %>% 
      layout(title="", 
             font=plotly_font,
             xaxis = list(title="",showline = T, tickangle = 45),
             yaxis = list(title="Share of individuals (%)",showline = T),
             margin = list(b = 160)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  
    
    
    plotty3 <- plot_ly(data =z %>% 
                         dplyr::rename(Onset_days = Age_first_word_months) %>% 
                         mutate(p_variant = paste0("p.",AA_ref,AA_pos,AA_alt)) %>% 
                         filter(!is.na(Onset_days)) %>% 
                         ungroup() %>% 
                         assign("save",.,envir = .GlobalEnv),
                       y = ~Onset_days, type = "box",x = "",
                       boxpoints = "all", jitter = 0.3,
                       pointpos = 0, hoverinfo = "text", showlegend = FALSE,
                       fillcolor = "rgba(31,119,180,0.5)",
                       boxpoints = 'all', marker = list (color = "#1f77b4"), line = list (color = "#1f77b4"),
                       text= ~paste0(round(Onset_days, digits = 2), " months, ", Protein)) %>% 
      layout(font=plotly_font,  
             title="",
             xaxis = list(title="", tickangle = 45, showline = T),
             yaxis = list(
               title = "Age first words (months)",
               tickmode = "array",
               showline = T
             ),
             margin = list(b = 160)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  
    subplot(plotty1,plotty2,plotty3, nrows = 1,titleY = T,margin = 0.05) %>% 
      layout(annotations = list(
        list(x = 0.1 , y = 1.1, text = "Cleft Palate", showarrow = F, xref='paper', yref='paper'),
        list(x = 0.5 , y = 1.1, text = "Total words spoken", showarrow = F, xref='paper', yref='paper'),
        list(x = 0.9 , y = 1.1, text = "Age at first word (months)", showarrow = F, xref='paper', yref='paper'))
      )
  })
  
  output$compareTable <- DT::renderDataTable({
    
    validate(need(
      !plyr::empty(varFilterInput$data),
      "There is no data that matches your filters."
    ))
    
    z <- Patient_data.df %>%
      filter(case_when(
        input$compareButtons =="Variant Type" ~ Vartype ==varFilterInput$data$Vartype,
        input$compareButtons =="Protein region" ~ Domain==varFilterInput$data$Domain,
        #input$compareButtons == "Functional Consequence" ~ CALL == varFilterInput$data$CALL,
        input$compareButtons =="Amino Acid Position" ~ AA_pos==varFilterInput$data$AA_pos)) 
    
    
    datatable(z %>% 
                select(Domain, Original_cDNA_change, Original_AA_change, Origin,Cleft_palate,Low_BMD,Abnormal_brainMRI,Age_walk_months,Age_first_word_months,Total_speech,Dental_issues,Clinical_seizures,Behavior_anomalies,Sleep_problems,Published_in), 
              extensions = "Buttons", 
              colnames = c("Domain","cDNA level","Protein level","Origin","CP","Low BMD","Abnl MRI","Walk at (months)", "Talk at (months)","Speech (words)","Abnl Teeth","Seizures", "Abnl behaviour","Abnl sleep","Link"),
              options = list(dom = 'Brtip',
                             buttons = c('csv', 'excel'), pageLength=100, scrollY = "350px"), escape = FALSE)
    
  })
  

  output$Var_analyis_paraz <- renderPlotly({
    
    validate(need(
      nrow(varFilterInput$data) >0,
      "There is no data that matches your filters."
    ))
    
    
    paraz_mtr_input.df <- paraz_mtr.df %>% 
      filter(Gene == varFilterInput$data$Gene[1]) %>% 
      left_join(Patient_data_missense_only.df %>% filter(Gene == varFilterInput$data$Gene[1], Vartype == "Missense") %>% mutate(p = 1) %>% select(p, AA_pos)) %>% 
      left_join(Control_data.df %>% filter(Gene == varFilterInput$data$Gene[1]) %>% mutate(g = 1) %>% select(g, AA_pos)) %>% 
      replace(is.na(.),0) %>% 
      assign("save",.,envir = .GlobalEnv) %>% 
      filter(p == 1) %>% 
      mutate(g = "pathogenic") %>% 
      rbind(save %>% filter(g == 1) %>% mutate(g = "control"))
    
    col1 <- "darkred"
    col2 <- "darkblue"
    
    plot <- plot_ly(colors=c(col1, col2)) %>% 
      add_boxplot(data = paraz_mtr_input.df %>% filter(Gene == varFilterInput$data$Gene[1], g == "pathogenic"),
                  y = ~Paraz_score, type = "box", x= 0,
                  color = I(col1)) %>% 
      add_trace(y = ~ paraz_mtr_input.df %>% filter(Gene == varFilterInput$data$Gene[1], g == "control") %>%  .$Paraz_score,
                type = "box", x = 1,
                color = I(col2)) %>% 
      add_trace(data = paraz_mtr.df%>%
                  filter(Gene == varFilterInput$data$Gene[1],
                         AA_pos == varFilterInput$data$AA_pos[1]) ,
                y = ~Paraz_score, x = 0.5,
                marker = list(color = "black",
                              size = 15),
                name = "Paraz-score",
                mode = "markers",
                type = "scatter") %>%
      add_segments(x = c(-0.5), 
                   xend = 1.5, 
                   y = 0, 
                   yend = 0,
                   line = list(color = 'gray', width = 4, dash = 'dot')) %>% 
      layout(yaxis = list(title = "Paraz-score"),
             xaxis = list(title = varFilterInput$data$Gene[1],
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    return(plot)
    
  })
  
  
  output$paraz_legend <- renderPlot({
    
    legend <- data.frame(x=c(1,4,7), y=c(3,3,3), text=c("Pathogenic variants", "Control variants", "Selected variant"))
    plot <- ggplot(legend, aes(x=x, y=y, color=text))+
      geom_point(size = 8)+
      scale_color_manual(values = c("darkblue","darkred","black"))+
      ylim(c(2.9,3.1))+
      xlim(c(0.8,9))+
      theme_void()+
      geom_text(aes(label=text), hjust=-0.2, color="black", size =5)+
      theme(legend.position = "none")
    
    
    return(plot)
    
  })
  
  
  output$Var_analyis_mtr <- renderPlotly({
    
    validate(need(
      nrow(varFilterInput$data) >0,
      "There is no data that matches your filters."
    ))
    
    
    paraz_mtr_gene.df <- paraz_mtr.df %>% filter(Gene == varFilterInput$data$Gene[1])
    
    mtr_25 <- paraz_mtr_gene.df %>% .$MTR_score %>% sort %>% .[floor(length(.)/4)]
    mtr_5 <- paraz_mtr_gene.df %>% .$MTR_score %>% sort %>% .[floor(length(.)/20)]
    
    
    paraz_mtr_input.df <- paraz_mtr.df %>% 
      filter(Gene == varFilterInput$data$Gene[1]) %>% 
      left_join(Patient_data_missense_only.df %>% filter(Gene == varFilterInput$data$Gene[1], Vartype == "Missense") %>% mutate(p = 1) %>% select(p, AA_pos)) %>% 
      left_join(Control_data.df %>% filter(Gene == varFilterInput$data$Gene[1]) %>% mutate(g = 1) %>% select(g, AA_pos)) %>% 
      replace(is.na(.),0) %>% 
      assign("save",.,envir = .GlobalEnv) %>% 
      filter(p == 1) %>% 
      mutate(g = "pathogenic") %>% 
      rbind(save %>% filter(g == 1) %>% mutate(g = "control"))
    
    col1 <- "darkred"
    col2 <- "darkblue"
    
    plot_ly(colors=c(col1, col2)) %>% 
      add_boxplot(data = paraz_mtr_input.df %>% filter(Gene == varFilterInput$data$Gene[1], g == "pathogenic"),
                  y = ~MTR_score, type = "box", x= 0,
                  color = I(col1)) %>% 
      add_trace(y = ~ paraz_mtr_input.df %>% filter(Gene == varFilterInput$data$Gene[1], g == "control") %>%  .$MTR_score,
                type = "box", x = 1,
                color = I(col2)) %>% 
      add_trace(data = paraz_mtr_gene.df%>% 
                  filter(AA_pos == varFilterInput$data$AA_pos[1]), 
                x = 0.5,
                y = ~MTR_score, x = "",
                marker = list(color = "black",
                              size = 15),
                name = "MTR-score",
                mode = "markers",
                type = "scatter") %>%  
      add_segments(x = c(-0.5), 
                   xend = 1.5, 
                   y = mtr_25, 
                   yend = mtr_25,
                   line = list(color = '#ffb366', width = 4, dash = 'dot')) %>% 
      add_segments(x = c(-0.5), 
                   xend = 1.5, 
                   name = "Lowest 5% MTR-score",
                   y = mtr_5, 
                   yend = mtr_5,
                   line = list(color = '#e60000', width = 4, dash = 'dot')) %>% 
      add_text(y = mtr_25-0.05,
               x = 0.15,
               text = "Lowest 25% MTR-score",
               textposition = "middle right") %>% 
      add_text(y = mtr_5-0.05,
               x = 0.15,
               text = "Lowest 5% MTR-score",
               textposition = "middle right") %>% 
      layout(yaxis = list(title = "MTR-score"),
             xaxis = list(title = varFilterInput$data$Gene[1],
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE) 
    
    
  })
  
  output$Var_analyis_per <- renderPlotly({
    
    validate(need(
      nrow(varFilterInput$data) >0,
      "There is no data that matches your filters."
    ))
    
    per_input.df <-per2d_for_var_analysis.df %>% 
      filter(Gene == varFilterInput$data$Gene[1]) %>% 
      left_join(Patient_data_missense_only.df %>% filter(Gene == varFilterInput$data$Gene[1], Vartype == "Missense") %>% mutate(p = 1) %>% select(p, AA_pos)) %>% 
      left_join(Control_data.df %>% filter(Gene == varFilterInput$data$Gene[1]) %>% mutate(g = 1) %>% select(g, AA_pos)) %>% 
      replace(is.na(.),0) %>% 
      assign("save",.,envir = .GlobalEnv) %>% 
      filter(p == 1) %>% 
      mutate(g = "pathogenic") %>% 
      rbind(save %>% filter(g == 1) %>% mutate(g = "control")) %>% 
      mutate(odds = log2(odds))
    
    col1 <- "darkred"
    col2 <- "darkblue"
    
    plot <- plot_ly(colors=c(col1, col2)) %>% 
      add_boxplot(data = per_input.df %>% filter(g == "pathogenic"),
                  y = ~odds, type = "box", x= 0,
                  color = I(col1)) %>% 
      add_trace(y = ~ per_input.df %>% filter(g == "control") %>%  .$odds,
                type = "box", x = 1,
                color = I(col2)) %>% 
      add_trace(data = per2d_for_var_analysis.df %>%
                  filter(AA_pos == varFilterInput$data$AA_pos[1],
                         Gene == varFilterInput$data$Gene[1]),
                y = ~odds, x = 0.5,
                marker = list(color = "black",
                              size = 15),
                name = "Fold enrichment of pathogenic variants (log2)",
                mode = "markers",
                type = "scatter") %>% 
      add_segments(x = c(-0.5), 
                   xend = 1.5, 
                   y = 0, 
                   yend = 0,
                   line = list(color = 'gray', width = 4, dash = 'dot')) %>% 
      layout(yaxis = list(title = "Fold enrichment of pathogenic variants (log2)"),
             xaxis = list(title = varFilterInput$data$Gene[1],
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    return(plot)
    
  })
  
  #####Research #####
  # Filter for subset of variants
  res_mod <- callModule(
    module = selectizeGroupServer,
    id = "research-filters",
    data = Patient_data.df %>% filter(Vartype != "Intronic"),
    vars = c("Vartype",  "AA_alt", "Domain","Clinical_seizures","Cleft_palate","Total_speech", "Abnormal_brainMRI")
  )
  
  output$filtered_n <- renderText({
    x <- nrow(res_mod())
    x <- paste("n = ", x)
    return(x)
  })
  
  # Table with displayed variants
  output$subsetTable <- DT::renderDataTable({
    req(res_mod())

     z <- res_mod() 
     
    patient_table <- datatable(z %>% 
                                 select(Transcript,Gene, Domain, cDNA, Protein, Phenotype, Onset_days,functional_effect, Published_in) , 
                               extensions = "Buttons",
                               colnames = c("Transcript","Gene", "Domain"," cDNA", "Protein", phenotype_name1,phenotype_name2, "Functional Consequence","Source"), 
                               escape=FALSE,
                               options = list(dom = 'Brtip',buttons = c('csv', 'excel'), pageLength=300, scrollY = "350px"))

    annotation_table <-  datatable(z %>% 
                                     select(Transcript,Gene, Domain, cDNA, Protein, Phenotype, Onset_days,functional_effect, Published_in), 
                                   extensions = "Buttons",
                                   colnames = c("Transcript","Gene", "Domain"," cDNA", "Protein", phenotype_name1,phenotype_name2, "Functional Consequence","Source"), 
                                   escape=FALSE,
                                   options = list(dom = 'Brtip',buttons = c('csv', 'excel'), pageLength=300, scrollY = "350px"))

    if (input$patientFunSwitch == FALSE) {
      return(patient_table)
    } else {
      return(annotation_table)
    }
  })

  ### Genotype Interface

  output$Genotype_overview_plot <- renderPlotly({
    

    # 2D lolliplot with ´SATB2 variants
    
    research_genotype_domain.df <- all_exchanges.df %>% 
      distinct(Domain,Gene,AA_pos,Domain_color) %>% 
      group_by(Gene,Domain,Domain_color) %>% 
      dplyr::summarise(start = min(AA_pos),
                       end = max(AA_pos))  
    

    g <- ggplot(data=all_exchanges.df %>% 
                  distinct(AA_pos,Gene,Domain,Domain_color) %>% 
                  left_join(res_mod() %>% 
                              filter(!is.na(AA_pos)) %>% 
                              select(Gene,Protein,AA_pos,Vartype,Original_cDNA_change) %>% 
                              mutate(Protein = ifelse(Vartype == "splice site",Original_cDNA_change,Protein)) %>% 
                              group_by(Protein,AA_pos,Gene,Vartype) %>% 
                              dplyr::summarise(var_count = n()) %>% 
                              ungroup() %>% 
                              mutate(Protein_count =  paste0(" ",Protein,", Variant count ",var_count)) %>% 
                              group_by(Gene,AA_pos,Vartype) %>% 
                              dplyr::summarise(Protein_final = paste(Protein_count, collapse = ";"))) %>% 
                  mutate(Vartype = ifelse(Vartype %in% c("Missense","PTV"),Vartype,
                                          ifelse(!is.na(Vartype),"Other",NA)))) +
      geom_segment(aes(x=AA_pos, xend=AA_pos, y=2, yend=ifelse(Vartype=="Missense", 7,8)), colour="black")+
      geom_point(aes(x=AA_pos, y=ifelse(Vartype=="Missense", 7,8), color=Vartype, text=Protein_final))+
      geom_rect(data=research_genotype_domain.df, aes(xmin=start, xmax=end, ymin=2, ymax=2.5, fill=Domain, text=Domain))+
      theme_classic()+
      ylim(c(1.7,10))+
      labs( x= "Amino acid sequence")+
      scale_color_manual(values = lolliplot_fill_scheme)+
      scale_fill_manual(values = lolliplot_fill_scheme)+
      #facet_grid(Gene ~ .)+
      theme(
        text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none")
    
   
    if (input$gnomad_m == TRUE) {
      g <- g + geom_point(data=Control_data.df ,
                          size=2, color = "black", aes(x=AA_pos, y=2, alpha=1, text=paste0("Position: ",AA_pos,", Allele count: ", Allele_count)))
    }
    

    g <- ggplotly(g, tooltip = "text") %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  %>%
      layout(title="",font=plotly_font,
             xaxis = list(title = "Amino acid sequence")
      )

  })
  
  output$Genotype_legend_plot <- renderPlot({
    
    legend <- data.frame(x=c(1,11,21), y=c(1, 1, 1), text=c("Missense", "   Null", "Control"))
    plot <- ggplot(legend, aes(x=x, y=y, color=text))+
      geom_point(size = 6)+
      scale_color_manual(values = c("Missense"="#D55E00","   Null"="#0072B2","Control" ="#000000"))+
      ylim(c(0,2))+
      xlim(c(0,40))+
      theme_void()+
      geom_text(aes(label=text), hjust=-0.4, color="black")+
      theme(legend.position = "none")
    
    return(plot)
    
  })
  
  
  output$threeDmolGene_all <- renderR3dmol({
    
    validate(need(
      nrow(res_mod() %>% filter(Vartype == "Missense")) >0,
      "There is no data that matches your filters."
    ))
    
    
    variant.df <- res_mod() %>%
      filter(Vartype == "Missense") %>%
      mutate(label = "pathogenic") %>%
      group_by(AA_pos,AA_ref,Gene) %>%
      dplyr::summarise(n_occ = n()) %>%
      select(AA_pos,AA_ref,n_occ,Gene)
    
    gnomad.df <- Control_data.df %>%
      group_by(AA_pos,AA_ref,Gene) %>%
      filter(Domain %in% unique(res_mod() %>% .$Domain)) %>% 
      dplyr::summarise(n_occ = n()) %>%
      select(AA_pos,AA_ref,n_occ,Gene)
    
    structure.df <- read_delim("data/pdb/SATB2_AF.txt",delim = "\t") %>%
      mutate(Aminoacid = aaa(Aminoacid)) %>%
      select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
    
    
    variant.df <- variant.df %>%
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position","AA_ref" = "Aminoacid","Gene" = "gene")) %>%
      mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not
      filter(struc_cov == "yes") %>% 
      distinct(Position_in_structure,Gene) %>% 
      group_by(Position_in_structure) %>% 
      dplyr::summarise(var_mut = ifelse(n() >1,"mutiple",Gene))
    
    gnomad.df <- gnomad.df %>%
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position","AA_ref" = "Aminoacid","Gene" = "gene")) %>%
      mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not
      filter(struc_cov == "yes") %>% 
      distinct(Position_in_structure,Gene) %>% 
      group_by(Position_in_structure) %>% 
      dplyr::summarise(var_mut = ifelse(n() >1,"mutiple",Gene))
    
    
    sub_color <- c("red","black")
    sub_scale <- c(1.2,0.8)
    struc_color <- "white"
    
    rot = 270
    rot_axis = "x"
    spin_axis = "vy"
    
    #Specify yourself- color of the cartoon per subunit
    subunit_color <- c("wheat","white") #first color for GRIN1 second or GRIN2A
    
    #Model for the protein complex
    
    modelo <- r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 50,
        upperZoomLimit = 1000
      )
    )
    
    modelo <- modelo %>% m_add_model(data = "data/pdb/SATB2.pdb", format = "pdb")
    
    # Zoom to encompass the whole scene
    modelo <-modelo %>% m_zoom_to() %>%
      # Set color o cartoon representation
      m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
      # Set subunit colors
      m_set_style(
        sel = m_sel(chain = c("A")),
        style = m_style_cartoon(color = subunit_color[2])
      ) %>%
      # visualize variants all
      m_set_style(
        sel = m_sel(resi = variant.df %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[1],
                               scale = sub_scale[1])
      ) 
    
    
    
    if  (input$gnomad_m == TRUE) {
      
      modelo <- modelo %>% m_set_style(
        sel = m_sel(resi = gnomad.df$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[2],
                               scale = sub_scale[2]))
      
    }
    return(modelo)
    
  })
  
  
  output$compareTableResearch <- DT::renderDataTable({
    req(res_mod())
    
    z <- res_mod()
    
    datatable(z %>% 
                select(Domain, Original_cDNA_change, Original_AA_change, Origin,Cleft_palate,Low_BMD,Abnormal_brainMRI,Age_walk_months,Age_first_word_months,Total_speech,Dental_issues,Clinical_seizures,Behavior_anomalies,Sleep_problems,Published_in), 
              extensions = "Buttons", 
              colnames = c("Domain","cDNA level","Protein level","Origin","CP","Low BMD","Abnl MRI","Walk at (months)", "Talk at (months)","Speech (words)","Abnl Teeth","Seizures", "Abnl behaviour","Abnl sleep","Link"),
              options = list(dom = 'Brtip',
                             buttons = c('csv', 'excel'), pageLength=100, scrollY = "350px"), escape = FALSE)
    
  })
  
  output$compareTableResearch2 <- DT::renderDataTable({
    req(res_mod())
    
    z <- res_mod()
    
    datatable(z %>% 
                select(Domain, Original_cDNA_change, Original_AA_change, Origin,Cleft_palate,Low_BMD,Abnormal_brainMRI,Age_walk_months,Age_first_word_months,Total_speech,Dental_issues,Clinical_seizures,Behavior_anomalies,Sleep_problems,Published_in), 
              extensions = "Buttons", 
              colnames = c("Domain","cDNA level","Protein level","Origin","CP","Low BMD","Abnl MRI","Walk at (months)", "Talk at (months)","Speech (words)","Abnl Teeth","Seizures", "Abnl behaviour","Abnl sleep","Link"),
              options = list(dom = 'Brtip',
                             buttons = c('csv', 'excel'), pageLength=100, scrollY = "350px"), escape = FALSE)
    
  })

  ### Phenotype Interface ####

  output$research_phenotype1 <- renderPlotly({

    plot <- plot_ly(res_mod() %>% 
                      dplyr::rename(phenotype_fac = Domain) %>% 
                      filter(phenotype_fac != "NA") %>% 
                      mutate(phenotype_fac =ifelse(phenotype_fac == "Yes"," Yes",phenotype_fac)) %>% 
                      select(phenotype_fac) %>%
                      arrange(phenotype_fac) %>% 
                      mutate(phenotype_fac = factor(phenotype_fac, 
                                                    levels = c("N-terminus","ULD Domain","Domain-Linker1","CUTL Domain","Domain-Linker2","CUT1 Domain","Domain-Linker3","CUT2 Domain","Domain-Linker4", "HOX","C-terminus"))) %>% 
                      group_by(phenotype_fac) %>% 
                      dplyr::summarise(n = n()) %>% 
                      mutate(n_gene = sum(n))%>% 
                      assign("save",.,envir = .GlobalEnv),
                    x = ~ phenotype_fac, 
                    y = ~ round(n, digits = 2), 
                    color = ~ phenotype_fac, 
                    colors = lolliplot_fill_scheme,
                    type = "bar", 
                    hoverinfo = "text", showlegend = FALSE,
                    text= ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
      layout(title = "Patient distribution according to variant location", 
             font=plotly_font,
             xaxis = list(title="",showline = T, tickangle = 45),
             yaxis = list(title="N of individuals",showline = T),
             margin = list(b = 160, t =100)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE) 

    return(plot)

  })

  output$research_phenotype2 <- renderPlotly({
    plot <- plot_ly(res_mod() %>% 
                     dplyr::rename(phenotype_fac = Total_speech) %>% 
                     filter(phenotype_fac != "NA") %>% 
                      filter(!is.na(phenotype_fac)) %>% 
                     mutate(phenotype_fac = factor(phenotype_fac, levels = c("None","1 to 10","10 to 50","Greater than 50"))) %>% 
                     group_by(phenotype_fac) %>% 
                     dplyr::summarise(n = n()) %>% 
                      mutate(n_gene = sum(n),
                             prop = n/n_gene*100)%>% 
                     assign("save",.,envir = .GlobalEnv),
                   x = ~ phenotype_fac, 
                   y = ~ prop, 
                   color = ~ phenotype_fac, 
                   colors = basic_phenotype_colors,
                   type = "bar", 
                   hoverinfo = "text", showlegend = FALSE,
                   text = ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
      layout(title = "Spoken words", 
             font=plotly_font,
             xaxis = list(title = "",showline = T, tickangle = 45),
             yaxis = list(title = "Share of individuals (%)", showline = T),
             margin = list(b = 160, t =100)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  

    return(plot)


  })

  output$research_phenotype3 <- renderPlotly({

    plot <- plot_ly(res_mod() %>% 
                      dplyr::rename(phenotype_fac = Clinical_seizures) %>% 
                      filter(phenotype_fac != "NA") %>% 
                      filter(!is.na(phenotype_fac)) %>% 
                      mutate(phenotype_fac =factor(phenotype_fac, levels = c("Yes","No"))) %>% 
                      group_by(phenotype_fac) %>% 
                      dplyr::summarise(n = n()) %>% 
                      mutate(n_gene = sum(n),
                             prop = n/n_gene*100)%>% 
                      assign("save",.,envir = .GlobalEnv),
                    x = ~ phenotype_fac, 
                    y = ~ prop, 
                    color = ~ phenotype_fac, 
                    colors = basic_phenotype_colors,
                    type = "bar", 
                    hoverinfo = "text", showlegend = FALSE,
                    text= ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
      layout(title = "Clinical seizures", 
             font = plotly_font,
             xaxis = list(title= "", showline = T, tickangle = 45),
             yaxis = list(title = "Share of individuals (%)", showline = T),
             margin = list(b = 160, t =100)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  
    

    return(plot)

  })
  
  output$research_phenotype4 <- renderPlotly({
    plot <- plot_ly(res_mod() %>% 
                      dplyr::rename(phenotype_fac = Dental_issues) %>% 
                      filter(phenotype_fac != "NA") %>% 
                      filter(!is.na(phenotype_fac)) %>% 
                      mutate(phenotype_fac =factor(phenotype_fac, levels = c("Yes","No"))) %>% 
                      group_by(phenotype_fac) %>% 
                      dplyr::summarise(n = n()) %>% 
                      mutate(n_gene = sum(n),
                             prop = n/n_gene*100)%>% 
                      assign("save",.,envir = .GlobalEnv),
                    x = ~ phenotype_fac, 
                    y = ~ prop, 
                    color = ~ phenotype_fac, 
                    colors = basic_phenotype_colors,
                    type = "bar", 
                    hoverinfo = "text", showlegend = FALSE,
                    text= ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
      layout(title = "Abnormal teeth", 
             font = plotly_font,
             xaxis = list(title= "", showline = T, tickangle = 45),
             yaxis = list(title = "Share of individuals (%)", showline = T),
             margin = list(b = 160, t =100)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  
    
    
    return(plot)
    
    
  })
  
  output$research_phenotype5 <- renderPlotly({
    plot <- plot_ly(res_mod() %>% 
                      dplyr::rename(phenotype_fac = Behavior_anomalies) %>% 
                      filter(phenotype_fac != "NA") %>% 
                      filter(!is.na(phenotype_fac)) %>% 
                      mutate(phenotype_fac =factor(phenotype_fac, levels = c("Yes","No"))) %>% 
                      group_by(phenotype_fac) %>% 
                      dplyr::summarise(n = n()) %>% 
                      mutate(n_gene = sum(n),
                             prop = n/n_gene*100)%>% 
                      assign("save",.,envir = .GlobalEnv),
                    x = ~ phenotype_fac, 
                    y = ~ prop, 
                    color = ~ phenotype_fac, 
                    colors = basic_phenotype_colors,
                    type = "bar", 
                    hoverinfo = "text", showlegend = FALSE,
                    text= ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
      layout(title = "Abnormal behavior", 
             font = plotly_font,
             xaxis = list(title= "", showline = T, tickangle = 45),
             yaxis = list(title = "Share of individuals (%)", showline = T),
             margin = list(b = 160, t =100)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  
    
    
    return(plot)
    
    
  })
  
  output$research_phenotype6 <- renderPlotly({
    plot <- plot_ly(res_mod() %>% 
                      dplyr::rename(phenotype_fac = Sleep_problems) %>% 
                      filter(phenotype_fac != "NA") %>% 
                      filter(!is.na(phenotype_fac)) %>% 
                      mutate(phenotype_fac =factor(phenotype_fac, levels = c("Yes","No"))) %>% 
                      group_by(phenotype_fac) %>% 
                      dplyr::summarise(n = n()) %>% 
                      mutate(n_gene = sum(n),
                             prop = n/n_gene*100)%>%  
                      assign("save",.,envir = .GlobalEnv),
                    x = ~ phenotype_fac, 
                    y = ~ prop, 
                    color = ~ phenotype_fac, 
                    colors = basic_phenotype_colors,
                    type = "bar", 
                    hoverinfo = "text", showlegend = FALSE,
                    text= ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
      layout(title = "Abnormal sleep", 
             font = plotly_font,
             xaxis = list(title= "", showline = T, tickangle = 45),
             yaxis = list(title = "Share of individuals (%)", showline = T),
             margin = list(b = 160, t =100)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  
    
    
    return(plot)
    
    
  })

  ### Functional Interface

  output$research_functional1 <- renderPlotly({
    
    data.df <- res_mod() %>% 
      rename(func_effect = "uptake") %>% 
      filter(!is.na(func_effect)) %>% 
      distinct(func_effect,AA_pos,AA_alt)
    
    validate(need(
      nrow(data.df) >0,
      "There is no data that matches your filters."
    ))
    
    
    plot <- plot_ly() %>% 
      add_boxplot(data = data.df ,
                  y = ~func_effect, type = "box", x= 0, boxpoints = "all",
                  color = I("gray")) %>% 
      layout(yaxis = list(title = "GABA uptake realtive to WT"),
             title = "GABA uptake rate",
             xaxis = list(title = "",
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F,
             margin = list(t = 60)) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    

    return(plot)
  })
  
  output$research_functional2 <- renderPlotly({
    
    data.df <- res_mod() %>% 
      rename(func_effect = "surface_exp") %>% 
      filter(!is.na(func_effect)) %>% 
      distinct(func_effect,AA_pos,AA_alt)
    
    validate(need(
      nrow(data.df) >0,
      "There is no data that matches your filters."
    ))
    
    
    plot <- plot_ly() %>% 
      add_boxplot(data = data.df ,
                  y = ~func_effect, type = "box", x= 0, boxpoints = "all",
                  color = I("gray")) %>% 
      layout(yaxis = list(title = "Surface expression realtive to WT"),
             title = "Surface expression",
             titlefont = list(size = 20),
             xaxis = list(title = "",
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F,
             margin = list(t = 60)) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    return(plot)
  })
  
  output$research_functional3 <- renderPlotly({
    
    data.df <- res_mod() %>% 
      rename(func_effect = "total_exp") %>% 
      filter(!is.na(func_effect)) %>% 
      distinct(func_effect,AA_pos,AA_alt)
    
    validate(need(
      nrow(data.df) >0,
      "There is no data that matches your filters."
    ))
    
    
    plot <- plot_ly() %>% 
      add_boxplot(data = data.df ,
                  y = ~func_effect, type = "box", x= 0, boxpoints = "all",
                  color = I("gray")) %>% 
      layout(yaxis = list(title = "Total expression realtive to WT"),
             title = "Total expression",
             titlefont = list(size = 20),
             xaxis = list(title = "",
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F,
             margin = list(t = 60)) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    return(plot)
  })
  
  output$research_functional4 <- renderPlotly({
    
    data.df <- res_mod() %>% 
      rename(func_effect = "relative_update_surface_exp") %>% 
      filter(!is.na(func_effect)) %>% 
      distinct(func_effect,AA_pos,AA_alt)
    
    validate(need(
      nrow(data.df) >0,
      "There is no data that matches your filters."
    ))
    
    
    plot <- plot_ly() %>% 
      add_boxplot(data = data.df ,
                  y = ~func_effect, type = "box", x= 0, boxpoints = "all",
                  color = I("gray")) %>% 
      layout(yaxis = list(title = "Relative uptake to surface expression"),
             title = "Relative uptake to surface expression",
             titlefont = list(size = 20),
             xaxis = list(title = "",
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F,
             margin = list(t = 60)) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    return(plot)
  })
  
  output$research_functional5 <- renderPlotly({
    
    data.df <- res_mod() %>% 
      rename(func_effect = "relative_surface_exp_tot_exp") %>% 
      filter(!is.na(func_effect)) %>% 
      distinct(func_effect,AA_pos,AA_alt)
    
    validate(need(
      nrow(data.df) >0,
      "There is no data that matches your filters."
    ))
    
    
    plot <- plot_ly() %>% 
      add_boxplot(data = data.df ,
                  y = ~func_effect, type = "box", x= 0, boxpoints = "all",
                  color = I("gray")) %>% 
      layout(yaxis = list(title = "Relative surface expression\n to total expression"),
             title = "Relative surface expression to total expression",
             titlefont = list(size = 20),
             xaxis = list(title = "",
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F,
             margin = list(t = 60)) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    return(plot)
  })
  
  output$functional_structure_legend_plot <- renderPlot({
    
    
    legend <- data.frame(x=c(1,7,14,21,28), y=c(1, 1, 1,1,1), text=c("LoF", "GoF", "Mixed", "Complex","STW"))
    plot <- ggplot(legend, aes(x=x, y=y, color=text))+
      geom_point(size = 6)+
      scale_color_manual(values = c("LoF"="#cc0000","GoF"="#66ffff","Mixed"="#ff8000","Complex"="#6600cc","STW"="#c0c0c0"))+
      ylim(c(0,2))+
      xlim(c(0,32))+
      theme_void()+
      geom_text(aes(label=text), hjust=-0.4, color="black")+
      theme(legend.position = "none")
    
    
    return(plot)
  })
  
  
  output$threeDmolfunctional <- renderR3dmol({
    
    
    selection_data.df <- res_mod() 
    
    if(selection_data.df$Gene %>% unique() %>% length()== 4){
      
      genes_selected <-master.df$Gene %>% unique()
    }else{
      genes_selected <- selection_data.df$Gene %>% unique()
    }
    
    
    variant.df <- master.df %>% 
      left_join(Functional_data.df %>% distinct(Gene,AA_pos,functional_effect), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) %>% 
      left_join(all_exchanges.df %>% distinct(Domain,Gene,AA_pos), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) %>% 
      mutate(AA_ref = aaa(AA_ref)) %>% 
      filter(Domain %in% selection_data.df$Domain,
             AA_ref %in% selection_data.df$AA_ref[!is.na(selection_data.df$AA_ref)],
             Gene %in% genes_selected,
             functional_effect %in% selection_data.df$functional_effect) %>% 
      filter(!is.na(functional_effect)) %>% 
      distinct(Aln_pos,Domain,functional_effect) %>% 
      group_by(Aln_pos) %>% 
      summarise(functional_effect = ifelse(unique(functional_effect) %>% length() ==1,functional_effect,"complex")) %>% 
      left_join(master.df %>% filter(Gene == "SCN2A") %>% distinct(AA_pos,AA_ref,Aln_pos,Gene)) 
    
    
    structure.df <- read_delim("data/pdb/6j8e_structure_coordinates.txt",delim = "\t") %>%
      select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
    
    
    variant.df <- variant.df %>%
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position","AA_ref" = "Aminoacid","Gene" = "gene")) %>%
      mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>% 
      filter(struc_cov == "yes") 
    
    
    sub_color <- c("#cc0000","#66ffff", "#ff8000","#6600cc","#c0c0c0")
    sub_scale <- c(1.2,0.8)
    struc_color <- "white"
    
    rot = 270
    rot_axis = "x"
    spin_axis = "vy"
    
    #Specify yourself- color of the cartoon per subunit
    subunit_color <- c("wheat","white") #first color for GRIN1 second or GRIN2A
    
    #Model for the protein complex
    
    modelo <- r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 50,
        upperZoomLimit = 1000
      )
    )
    
    modelo <- modelo %>% m_add_model(data = "data/pdb/SCN2A_model.pdb1", format = "pdb")
    
    # Zoom to encompass the whole scene
    modelo <-modelo %>% m_zoom_to() %>%
      # Set color o cartoon representation
      m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
      # Set subunit colors
      m_set_style(
        sel = m_sel(chain = c("A")),
        style = m_style_cartoon(color = subunit_color[2])
      ) %>%
      # visualize variants all
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(functional_effect == "LoF") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[1],
                               scale = sub_scale[1])
      ) %>% 
      #visualize variants SCN1A
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(functional_effect == "GoF") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[2],
                               scale = sub_scale[1])
      ) %>% 
      #visualize variants SCN2A
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(functional_effect == "Mixed") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[3],
                               scale = sub_scale[1])
      )  %>% 
      # visualize variants SCN3A
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(functional_effect == "complex") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[4],
                               scale = sub_scale[1])
      ) %>% 
      # visualize variants SCN8A
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(functional_effect == "STW") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[5],
                               scale = sub_scale[1])
      ) 
    
    
    return(modelo)
    
  })
  
})
