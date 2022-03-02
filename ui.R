################################################
################# SATB2 PORTAL ##################
################################################
################## LIBRARIES ################

library(shiny)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(plotly)
library(tippy)
library(r3dmol)
library(DT)
library(readr)
library(tidyverse)
library(vembedr)
library(RColorBrewer)
library(shinyhelper)
library(plyr)
library(rsconnect)

# CSS/STYLE/INFO #
landing_panel <- "color: #333333;
      height: 200px;
      width: 260px"

spinner_color <- "#2196f3"

sub_style <- "color:gray; font-style: italic; font-size: 14px;"

#Layout colors 

schema_color <- "#85C1E9"

schema_color_strong <- "#2E86C1"

schema_color_light <- "gainsboro"

##### General Variable #####

contact_us <- "mailto:lald@ccf.org; stefana4@ccf.org; yazarate@uams.edu"
contact_yuri <- "mailto:yazarate@uams.edu"
contact_katie <- "mailto:kbbosanko@uams.edu"

##### Landing Page Variable #####
#General
landing_portal_title =  "SATB2 Portal"

landing_bannername = "banner01_italic.png"

#Tabs 

landing_tab1 = p(HTML(paste0("Features and genetics of ", em("SATB2"), "-Asssociated syndrome")))
landing_tab2 = p("Resources for families and caregivers\n     ",
                 br(),
                 br())
               #  "    ")
landing_tab3 = p("Information on individual variant interpretation")
landing_tab4 = p("Comprehensive genotype-phenotype correlations")
landing_tab5 = p(HTML(paste0("Information on the ", em("SATB2"), " Registry and multidisciplinary ", em("SATB2"), " clinic")))

#####Basic Information Variable #####
#General 
basic_title = HTML(paste0(em("SATB2"),"-Associated syndrome"))

##Genes 
#Gene 1
gene1 = "SATB2"
gene1_info_header = div(HTML(paste0("Information about ",em("SATB2"))),style = "font-size:17px;")

basic_text_title_1 <- HTML(paste0(em("SATB2"),"-Associated syndrome (Glass syndrome)"))

quick_facts_png <- "SATB2_quick_facts.png"

#Visualisations
basic_visual_title = "Clinical Information of our patients"

basic_phenotype_title1 = "Population breakdown"
basic_phenotype_title2 = "Speech and walking"
basic_phenotype_title3 = "Verbal output" 
basic_phenotype_title4 = "Clefte plate & Low bone density" 
basic_phenotype_title5 = "Other phenotypes"

basic_plot_title_fac1_1 = "Sex"
basic_plot_title_fac1_2 = "Age Distribution"

basic_plot_title_fac2_1 = "Age at first word"
basic_plot_title_fac2_2 = "Age at first step"

basic_plot_title_fac3_1 = "Total words spoken"

basic_plot_title_fac4_1 = "Cleft palate"
basic_plot_title_fac4_2 = "Low bone density"

basic_plot_title_fac3 = "Epilepsy syndromes"

basic_plot_title_num1 = ""
basic_plot_title_num2 = ""

basic_abbreviations1 = ""

basic_abbreviations2 = "PTV: Protein truncating variants"

basic_abbreviations3 = "EMAS = Epilepsy with myoclonic-atonic seizures; GGE = Genetic generalized epilepsy;
NAFE = Non-acquired focal epilepsy; EOAE = Early onset absence epilepsy; DEE = Developmental epileptic encephalopathy; CAE = Childhood absence epilepsy; TLE = Temporal lobe epilepsy; LGS = Lennox-Gastaut syndrome"

basic_abbreviations4 = "ID = Intellectual disability; DD = Developmental delay"

#tippy
epi_syndrome_tippy <- h5(HTML(paste0("<ul><li>EMAS = Epilepsy with myoclonic-atonic seizures</li>")),
                         HTML(paste0("<li>GGE = Genetic generalized epilepsy</li>")), 
                         HTML(paste0("<li>NAFE = Non-acquired focal epilepsy</li>")), 
                         HTML(paste0("<li>EOAE = Early onset absence epilepsy</li>")), 
                         HTML(paste0("<li>DEE = Developmental epileptic encephalopathy</li>")), 
                         HTML(paste0("<li>CAE = Childhood absence epilepsy</li>")), 
                         HTML(paste0("<li>TLE = Temporal lobe epilepsy</li>")), 
                         HTML(paste0("<li>LGS = Lennox-Gastaut syndrome</li>")), 
                         HTML(paste0("<li>No SZ = No seizures</li>")), 
                         HTML(paste0("<li>UG = Unclassified general</li>")), align = "left")

dd_id_tippy <- h5(HTML(paste0("<ul><li>DD = Developmental delay</li>")),
                  HTML(paste0("<li>ID = Intellectual disability</li>")), align = "left")

Gene1_basic_text <- p(HTML(paste0("",em("SATB2"),"-Associated syndrome (SAS, Glass syndrome) is an autosomal dominant multisystem disorder characterized by significant
                                  neurodevelopmental compromise with limited to absent speech, intellectual disability of variable severity, behavioral issues, skeletal 
                                  anomalies, and craniofacial abnormalities. More details can be found on ", 
                                  shiny::a("Genereviews.", href = "https://www.ncbi.nlm.nih.gov/books/NBK458647/"))))


#####Families Variable #####  
family_title <- HTML(paste0("What is ",em("SATB2"),"Associated syndrome"))

family_video_url <- "https://www.youtube.com/watch?v=kCcfPM5s3iM"

detailed_video_url <- "https://www.youtube.com/watch?v=IwpdBfGuxus"

family_video_sub_description <- HTML(paste0("Short video that describes the basics of ", em("SATB2"), "-Associated syndrome"))

detailed_video_sub_description <- HTML(paste0("Detailed video from the 2019 International ",em("SATB2")," Family meeting with Dr. Zarate"))

family_video_description_title <- HTML(paste0("What is ",em("SATB2"),"-Associated syndrome"))

family_video_description_content <- p(em("SATB2"),"-Associated syndrome (SAS) is a rare genetic condition. While knowledge and awareness on SAS continues to improve, most people have not heard of it. The core features of SAS can be remembered by the following acronym using the name of the ", em("SATB2"), " gene:", 
                                      br(),
                                      br(),
                                      div("S: Severe speech anomalies"), 
                                      div("A: Abnormalities of the palate"),
                                      div("T: Teeth anomalies"),
                                      div("B: Behavioral issues with or without Bone or Brain MRI anomalies"),
                                      div("2: 2 years of age or earlier is usually the age of onset"),
                                      br(),
                                      br(),
                                      ("In these videos you'll get an overview of the main features of SAS. You have the option of a simple and quick version or a more detailed version with more information."))



##### Variant Analysis variable #####

master.df <- read_delim("data/master_table.txt", delim = "\t") 

var_possible_genes_title <- "Select gene"

var_possible_genes <- "SATB2"

var_possible_phenotype <- c("Epilepsy + ID", "Epilepsy without ID",
                            "Autism", "Other")

var_patient_info_title <- p("Individuals with the same variant in the SATB2 Portal",
                            tippy(icon("question-circle"),
                                  tooltip = h5(HTML(paste0("Variant collection:")),
                                               HTML(paste0("<ul><li>Zarate et al.,????</li>")),
                                               HTML(paste0("<li>Internal variant database (unpublished)</li>")), align = "left"),
                                                           animation = "scale", theme = "light"))
var_patient_info_abb <- "Abnl: abnormal; BMD: bone mineral density; CP: cleft palate"


var_paralog_info_abb <- ""


variants_compareplots_abb <- "DEE: Developmental and Epileptic Encephalopathies; BNFS: Benign Familial Neonatal Seizures; EOEE: Early-Onset Epileptic Encephalopathies; ASD: Autism Spectrum Disorder "

#####Research Variable ####
research_pheno_abb <- "DD: Developmental delay, ID: Intellectual disability"
research_pheno_abb2 <- "EMAS = Epilepsy with myoclonic-atonic seizures; GGE = Genetic generalized epilepsy; UG = Unclassified general;
NAFE = Non-acquired focal epilepsy; EOAE = Early onset absence epilepsy; DEE = Developmental epileptic encephalopathy; CAE = Childhood absence epilepsy; TLE = Temporal lobe epilepsy; LGS = Lennox-Gastaut syndrome; No SZ = No seizures; NA = Not available"

research_geno_transcripts <- p("The following transcript was used:",em("SATB2"),": ENST00000287766",
                               br(),
                               "H: Helix")

#Datasets required for research tab 
Patient_data.df <- read_delim("data/SATB2_Patient_variants_v1.txt", delim = "\t") 
# %>%
  # rename(Sz_onset = "Age at seizure onset (months)",
  #        Epilepsy = "Epilepsy") %>% 
  # rename(Autism = "Autistic traits",
  #        Epilepsy_syndrome = "Epilepsy Syndrome Classification",
  #        ID_after_sz_onset = "Cognitive Level AFTER Seizure Onset")

Patient_data_missense.df <- read_delim("data/SATB2_Patient_variants_v1.txt", delim = "\t") %>% 
  mutate(Sex = as.factor(Sex)) %>% 
  filter(Vartype == "Missense")

##### About Variable #####
about_terms_of_use <- p("All data here are publicly for the benefit of researchers, clinicians, and caregivers. You can freely explore the data, and we encourage the use and publication of results generated from these data. However, we encourage you to contact us before embarking on analyses to check if your proposed analysis overlaps with work currently underway by our team. Further, we request that you actively acknowledge and give attribution to the SATB2 Portal project, and link back to the relevant page, wherever possible. All users of our data agree to not attempt to reidentify participants. Our data set has been subjected to extensive quality control, but maybe imperfect so errors may remain.
If you spot mistakes or have suggestions for SATB2 Portal improvements:")

about_data_generation <- "Information for this portal is obtained by careful review of the literature and extensive family-provided data curation. Standardized data collection forms are filled and aggregated in all instances. While we procure obtaining medical records for all individuals enrolled in the registry, in some instances, we can’t verify all information provided to be accurate. For cases reported in the literature, data is obtained from the original publication as reported by the authors."


## Tutorials variables
tut_basic_info <- "https://youtu.be/q9Usxkj_wBs"

tut_family <- "https://youtu.be/FdVgO6bWaWc"

tut_variant_analysis <- "https://youtu.be/gCPjOzyH5u8"

tut_research <- "https://youtu.be/i2w8Z2qH9kQ"



################# UI #################
shinyUI(
    ##### LANDING PAGE #####
    ui <- 
      navbarPage(#fluid = FALSE, 
        windowTitle = "SATB2 Portal",  
        id = "TabDisplay",
        theme = "mytheme2.css",
        title = p(icon("dna"), landing_portal_title), 
        header = tagList(useShinydashboard()),
        
        tabPanel(
          tags$style(HTML(paste0("
                      .navbar-default .navbar-brand {
                      background-color: ",schema_color_light,";
                      color: black
                    }"))),
          tags$style(HTML(paste0("
                      .navbar-default {
                      background-color: ",schema_color_light,";
                    }"))),
          
          title = "Welcome",
          value = "welcomeTab",
          tabName = "welcomeTab",
          
          div(style = "font-size:100%", 
              fluidRow(
                column(12, 
                   style = "background-color: white; color: #676767",
                   align = "center",
                   br(), br(),
                   img(src = landing_bannername, width = "100%")   #insert your banner, saved in the www-folder 
            )
          ),
          
          div(style = 'background-color: #F8FCFE',
              fluidRow(
                
                # Basic Information Tab
                div(width = "100%", 
                    column(2, offset = 1, align = "center",
                           br(),
                           panel(width = 12, 
                                 status = "success",
                                 heading = "",
                                 h2(tags$i(class = "fa fa-dna", style = "color:#676767")),
                                 br(),
                                 div(landing_tab1,  #c
                                     br(), br(),
                                     actionBttn(
                                       inputId = "infoBtn",
                                       label = p("Basic Information",
                                                 br(),
                                                 br()),
                                       color = "success",
                                       block = TRUE,
                                       size = "md",
                                       style = "stretch")),
                                 style =  "background-color: #f3faf4;",   #set color of the tab
                           ), 
                    ),
                    
                    #Family tab 
                    column(2,align = "center",
                           br(),
                           div(panel(
                             status = "warning",
                             heading = "",
                             h2(tags$i(class = "fa fa-child", style = "color:#676767")),
                             br(),
                             div(landing_tab2,
                                 br(), 
                                 br(),
                                 actionBttn(
                                   inputId = "familyBtn",
                                   label = p("Families",
                                             br(),
                                             br()),
                                   color = "warning",
                                   block = TRUE,
                                   style = "stretch"
                                 )),
                             style = "background-color: #fff8ef;",   #set color of the tab
                           ))
                    ),
                    
                    # Variant Analysis Tab
                    column(2, align = "center",
                           br(),
                           div(panel(
                             status = "info",
                             heading = "",
                             h2(tags$i(class = "fa fa-code-branch", style = "color:#676767")),
                             br(),
                             div(landing_tab3,br(), 
                                 br(),
                                 actionBttn(
                                   inputId = "variantBtn",
                                   label = p("Variant Analysis",
                                            br(),
                                            br()),
                                   color = "royal",
                                   block = TRUE,
                                   style = "stretch"
                                 )),
                             style = "background-color: #f9f1fa;",    #set color of the tab
                           ))
                    ),
                    
                    # Research Tab
                    column(2, align = "center",
                           br(),
                           div(panel(
                             status = "danger",
                             heading = "",
                             h2(tags$i(class = "fa fa-microscope", style = "color:#676767")),
                             br(),
                             div(landing_tab4,br(), 
                                 br(),
                                 actionBttn(
                                   inputId = "researchBtn",
                                   label = p("Genotype-Phenotype",
                                             br(),
                                             br()),
                                   color = "danger",
                                   block = TRUE,
                                   style = "stretch"
                                 )),
                             style = "background-color: #fdf0f1;",
                             
                           ))
                    ),
                   # Registry Tab
                    column(2, align = "center", 
                           br(),
                        div(panel(
                            status = "primary",
                            heading = "",
                            h2(tags$i(class = "fa fa-laptop-code", style = "color:#676767")),
                            br(),
                            div(landing_tab5,br(),  #c
                            br(),
                            actionBttn(
                                #class = "btn-primary",
                                inputId = "registryBtn",
                                label = p(em("SATB2 ", "Registry and Clinic")),   #c
                                color = "primary",
                                block = TRUE,
                                style = "stretch"
                            )), 
                            style = "background-color: #f1f8fe;",   #set color of the tab
                        )) , br(), br(), br(), br(), br(), br(), br(), br(), br()
                    )
                )),
              fluidRow(
                column(4, offset = 2, align = "center",
                       p("Visit our other Portals:", style = "font-size:14px;", align = "center", style = "font-size:14px;"),
                       p(shiny::a("SCN-Portal", href="http://scn-portal.broadinstitute.org/", target = '_blank'), "& ",
                         shiny::a("GRIN-Portal", href="http://grin-portal.broadinstitute.org/", target = '_blank'),
                         style = "font-size:12px;")
                ),
                column(4, align = "center",
                       
                       p("You want to join the project or provide feedback?", align = "center",style = "font-size:14px;"),
                       p("Please ", shiny::a("contact us!", href=contact_us), style = "font-size:12px;", align = "center"),
                       br()
                ),
                column(12, align = "center",offset = 4,
                      fluidRow(
                        div(width = "100%",valueBoxOutput("disclaimer")))
               )))
            )), # end tab Panel

        ##### BASIC INFORMATION #####
        tabPanel(title = "Basic Information", value = "infoTab",

                 fluidRow(
                   column(10, offset = 1,
                     panel(heading = basic_title,
                           status = "success",
                       fluidRow(
                         column(4,
                           box(width = 12,
                             title = p("History of ",em(gene1), " research"),
                             timelineBlock(
                                 reversed = FALSE,
                                 width = 12,

                                 #Add as many blocks as desired for new publications
                                 timelineLabel(1989, color = "teal"),
                                 timelineItem(
                                   title = div(strong("First report of an individual with a chromosome 2 deletion including the",em("SATB2")," gene and the origin of the term 'Glass syndrome'")),
                                   icon = icon("dna"),
                                   color = "olive",
                                   #ADD link to the publication
                                   time = shiny::a("Glass et al.", href="https://pubmed.ncbi.nlm.nih.gov/2918541/", target = '_blank'),
                                   border = FALSE,
                                 ),
                                 timelineLabel(2003, color = "teal"),
                                 timelineItem(
                                     title = div(em("SATB2"), " gene discovery"),
                                     icon = icon("user"),
                                     color = "aqua",
                                     #ADD link to the publication
                                     time = shiny::a("FitzPatrick et al.", href="https://pubmed.ncbi.nlm.nih.gov/12915443/", target = '_blank'),
                                     border = FALSE,
                                 ),
                                 timelineLabel(2017, color = "teal"),
                                 timelineItem(
                                   title = HTML(paste0("Delineation of the phenotype and functional studies")),
                                   icon = icon("user"),
                                   color = "aqua",
                                   #ADD link to the publication
                                   time = shiny::a("Bengani et al.", href="https://pubmed.ncbi.nlm.nih.gov/28151491/", target = '_blank'),
                                   border = FALSE,
                                 ),
                                 timelineLabel(2018, color = "teal"),
                                 timelineItem(
                                   title = HTML(paste0("Delineation of the phenotype and natural history of SAS")),
                                   icon = icon("user"),
                                   color = "aqua",
                                   #ADD link to the publication
                                   time = shiny::a("Zarate et al.", href="https://pubmed.ncbi.nlm.nih.gov/29436146/", target = '_blank'),
                                   border = FALSE,
                                 ),
                                 timelineLabel(2019, color = "teal"),
                                 timelineItem(
                                   title = div("Comprehensive mutation spectrum overview"),
                                   icon = icon("user"),
                                   color = "aqua",
                                   #ADD link to the publication
                                   time = shiny::a("Zarate et al.", href="https://pubmed.ncbi.nlm.nih.gov/31021519/", target = '_blank'),
                                   border = FALSE,
                                 ),
                             )
                           )),
                           column(8, # image removed and here, to increase column width 4->8. Change back to 4, when adding an image again
                            box(title=basic_text_title_1, width = 12,
                                  Gene1_basic_text)),
                           # column(4, align = "center", # image to the right of basic info text
                           #   box(title =  HTML(paste0(em("SATB2") ," main characteristics")), width = 12, align = "center",
                           #       img(src = quick_facts_png, width = "100%"))),

                           column(8,
                             box(title = p("Clinical information from the SAS Registry and Literature for individuals with point mutations (310 individuals):",
                                           tippy(icon("question-circle"),
                                 tooltip = h5(HTML(paste0("The ", strong(em("SATB2")), strong(" mutation database "), "is an aggregate of the following sources, as of October 17, 2021:")),
                                              br(),br(),
                                              HTML(paste0("<ul><li>Zarate et al., 2019 (n = 49; PMID:31021519)</li>")),
                                              HTML(paste0("<li>Zarate et al., 2018 (n = 46; PMID:29436146)</li>")),
                                              HTML(paste0("<li>Bengani et al., 2017 (n = 11; PMID:28151491)</li>")),
                                              HTML(paste0("<li>ClinVar (n = 44)</li>")),
                                              HTML(paste0("<li>Internal variant database (n = 123; unpublished)</li>")),
                                              HTML(paste0("<li>Others</li>")), align = "left"),
                                 animation = "scale",
                                 theme = "light")),
                                 width = 12,

                              tabsetPanel(
                              #Phenotype summary
                                tabPanel(title = basic_phenotype_title1,
                                         br(),
                                      column(6,
                                         div(p("n = " ,Patient_data.df %>% filter(Gene == gene1,!is.na(Sex)) %>% nrow(), ""), style = "font-size:15px;color:black;")),
                                      column(6,
                                             div(p("n = " ,Patient_data.df %>% filter(Gene == gene1,!is.na(Age_years)) %>% nrow(), ""),style = "font-size:15px;color:black;")),
                                 fluidRow(align = "center",
                                      column(6,
                                            box(title=div(basic_plot_title_fac1_1, style = "font-size: 15px"), width=12, plotlyOutput(outputId = "Phenotype_fac1_1"))),
                                     column(6,
                                           box(title=div(basic_plot_title_fac1_2, style = "font-size: 15px"), width=12, plotlyOutput(outputId = "Phenotype_num1_2"))),
                                     column(12,
                                            align="center",
                                            p(basic_abbreviations1, style=sub_style)))),
                                
                                ##Tab 2
                                tabPanel(title = basic_phenotype_title2,
                                         br(),
                                         column(6,
                                                div(p("n = " ,Patient_data.df %>% filter(Gene == gene1,!is.na(Age_first_word_months)) %>% nrow(), ""),style = "font-size:15px;color:black;")),
                                         column(6,
                                                div(p("n =" ,Patient_data.df %>% filter(Gene == gene1,!is.na(Age_walk_months)) %>% nrow(), ""),style = "font-size:15px;color:black;")),
                                         fluidRow(align = "center",
                                                  column(6,
                                                         box(title=div(basic_plot_title_fac2_1, style = "font-size: 15px"), width=12, plotlyOutput(outputId = "Phenotype_num2_1"))),
                                                  column(6,
                                                         box(title=div(basic_plot_title_fac2_2, style = "font-size: 15px"), width=12, plotlyOutput(outputId = "Phenotype_num2_2"))),
                                                  column(12,
                                                         align="center",
                                                         p(basic_abbreviations1, style=sub_style)))),
                                
                                #Tab3 
                                tabPanel(title = basic_phenotype_title3,
                                         br(),
                                         column(6,
                                                div(p("n =" ,Patient_data.df %>% filter(Gene == gene1,!is.na(Total_speech)) %>% nrow(), ""),style = "font-size:15px;color:black;")),
                                         column(12,
                                                ""),
                                         fluidRow(align = "center",
                                                  column(6,
                                                         box(title=div(basic_plot_title_fac3_1, style = "font-size: 15px"), width=12, plotlyOutput(outputId = "Phenotype_fac3_1"))),
                                                  column(12,
                                                         align="center",
                                                         p(basic_abbreviations1, style=sub_style)))),
                              #Tab4 
                                tabPanel(title = basic_phenotype_title4,
                                         br(),
                                         column(6,
                                                div(p("n = " ,Patient_data.df %>% filter(Gene == gene1,!is.na(Cleft_palate)) %>% nrow(), ""),style = "font-size:15px;color:black;")),
                                         column(6,
                                                div(p("n = " ,Patient_data.df %>% filter(Gene == gene1,!is.na(Low_BMD)) %>% nrow(), ""),style = "font-size:15px;color:black;")),
                                         fluidRow(align = "center",
                                                  column(6,
                                                         box(title=div(basic_plot_title_fac4_1, style = "font-size: 15px"), width=12, plotlyOutput(outputId = "Phenotype_fac4_1"))),
                                                  column(6,
                                                         box(title=div(basic_plot_title_fac4_2, style = "font-size: 15px"), width=12, plotlyOutput(outputId = "Phenotype_fac4_2"))),
                                                  column(12,
                                                         align="center",
                                                         p(basic_abbreviations1, style=sub_style)))),
                              #Tab5
                              tabPanel(title = basic_phenotype_title5,
                                       br(),
                                       fluidRow(align = "center",
                                                column(6, box(title=div("Further Phenotypes", style = "font-size: 15px"), width=12, align="center",
                                                              img(src = "worldcloud_phenotype.png", width = "100%"),
                                                              footer = div("Frequently described phenotypes associaed with SATB2."))),
                                                column(12,
                                                       align="center",
                                                       p(basic_abbreviations1, style=sub_style))))
                             ))))
        )))), #end basic info

##### FAMILIES #####

tabPanel(title = "Families", value = "familyTab",
         tags$style(HTML(paste0("
                      .box.box-solid.box-danger>.box-header {
                      background-color: ",schema_color_light,";
                      border: black
                      }"))),
         tags$style(HTML(paste0("
                      .box.box-solid.box-danger {
                      border: 1px solid black
                      }"))),
         fluidRow(
           column(10, offset = 1,
                  panel(heading = "Families", status = "warning",
                        tabsetPanel(
                          # First tab for videos
                          tabPanel(title = HTML(paste0("What is ",em("SATB2"),"-Associated syndrome?")),
                                   fluidRow(
                                     column(5, offset = 1,
                                            br(),
                                            br(),
                                            box(title= HTML(paste0("What is ",em("SATB2"),"-Associated syndrome?")),
                                                width = 12,
                                                embed_url(family_video_url) %>% use_rounded() %>% use_align("center") %>% use_bs_responsive(),
                                                footer = p(div(align="center", family_video_sub_description, style=sub_style))
                                            ),
                                            box(title= HTML(paste0(em("SATB2"), "-Associated syndrome Webinar Fall 2019")),
                                                width = 12,
                                                embed_url(detailed_video_url) %>% use_rounded() %>% use_align("center") %>% use_bs_responsive(),
                                                footer = p(div(align = "center", detailed_video_sub_description, style = sub_style))
                                            )),
                                     column(5, offset=0, br(), align = "justify",  br(), br(), br(),
                                            box(title= div(family_video_description_title, style = "color:black;"),
                                                width = 12,
                                                status = "danger", solidHeader = TRUE,
                                                family_video_description_content
                                            ))
                                   )),
                          
                          # Second tab for Logos linked to family organizations
                          tabPanel(title = "Family Organizations",
                                   panel(status="default", heading = "Family Organizations", width = 12,
                                         fluidRow(
                                           column(12, align = "center",
                                                  offset = 0,
                                                  div(style = "display: inline-block; width: 40%;", 
                                                      width = 12, shiny::a(
                                                        img(src = "SATB2_LOGO_1.jpg", 
                                                            width = 300, 
                                                            height = 300),
                                                        href = "https://satb2gene.org/",
                                                        target = '_blank')),
                                                  div(style = "display: inline-block; width: 40%;", 
                                                      width = 12, shiny::a(
                                                        img(src = "SATB2_LOGO_Australia.jpg", 
                                                            width = 300, 
                                                            height = 300),
                                                        href = "https://satb2.org.au/",
                                                        target = '_blank')),
                                                  div(style = "display: inline-block; width: 40%;", 
                                                      width = 12, shiny::a(
                                                        img(src = "SATB2_LOGO_esp.png", 
                                                            width = 600, 
                                                            height = 600),
                                                        href = "http://satb2.es/",
                                                        target = '_blank')),
                                                  div(style = "display: inline-block; width: 40%",
                                                      width = 12, shiny::a(
                                                        img(src = "SATB2_LOGO_2_Europe.jpg", 
                                                            width = 300, 
                                                            height = 300),
                                                        href = "https://www.satb2europe.org/", 
                                                        target = '_blank'))
                                           ),
                                           column(12, align = "center",
                                                  br(),
                                                  br(),
                                                  div("Please", 
                                                      shiny::a("contact us", 
                                                               href="mailto:yazarate@uams.edu") ,"if you would like your organization to be listed here."))
                                         )
                                   ),
                                   # panel(status="default",heading = "Related Organizations", width = 12,
                                   #   fluidRow(
                                   #     column(3, align = "center", offset = 1,
                                   #            div(width = 12, shiny::a(
                                   #               img(src = "aes.png", width = "70%"),
                                   #               href = "https://www.aesnet.org/",
                                   #               target = '_blank'))
                                   #     ),
                                   #     column(3, align = "center",
                                   #            div(width = 12, shiny::a(
                                   #               img(src = "nord.png", width = "70%"),
                                   #               href = "https://rarediseases.org/",
                                   #               target = '_blank'))
                                   #     ),
                                   #     column(3, align = "center",
                                   #            div(width = 12, shiny::a(
                                   #               img(src = "epi_foundation.png", width = "70%"),
                                   #               href = "https://www.epilepsy.com/",
                                   #               target = '_blank'))
                                   #     ),
                                   #     column(3, align = "center", offset = 1,
                                   #            div(width = 12, shiny::a(
                                   #              img(src = "combined.png", width = "70%"),
                                   #              href = "https://combinedbrain.org/",
                                   #              target = '_blank'))
                                   #     ),
                                   #     column(3, align = "center",
                                   #            div(width = 12, shiny::a(
                                   #               img(src = "global_genes.png", width = "70%"),
                                   #               href = "https://globalgenes.org/",
                                   #               target = '_blank'))
                                   #     ),
                                   #     column(3, align = "center",
                                   #            div(width = 12, shiny::a(
                                   #              img(src = "ren.png", width = "70%"),
                                   #              href = "https://www.rareepilepsynetwork.org/",
                                   #              target = '_blank'))
                                   #     ),
                                   #     column(3, align = "center",
                                   #            div(width = 12, shiny::a(
                                   #               img(src = "epi_council.png", width = "70%"),
                                   #               href = "https://www.epilepsyleadershipcouncil.org/",
                                   #               target = '_blank'))
                                   #     ),
                                   #     column(12, align = "center",
                                   #            div("...and many more! Please", 
                                   #                shiny::a("contact us", 
                                   #                         href="mailto:GRIN@medizin.uni-leipzig.de") ,"if you would like your organization to be listed here."))
                                   #   ))
                          ),
                          # Third tab for Resources
                          tabPanel(title = "Resources",
                                   panel(status="default", heading = "Helpful documents with detailed information on important topics", width = 12,
                                         # h2("pdf", tags$iframe(style = "height:400px; width:100%;scrolling = yes", src = "03Background_info.pdf")),
                                         fluidRow(
                                           column(10, align = "left",
                                                  br(),
                                                  div(style = "display: inline-block; width: 40%;", 
                                                      width = 12, shiny::a(
                                                        img(src = "pdf_icon.png", 
                                                            width = 16, 
                                                            height = 20),
                                                        href = "03Background_info.pdf",
                                                        target = '_blank'), "Background Information"),
                                                  div(style = "display: inline-block; width: 40%;", 
                                                      width = 12, shiny::a(
                                                        img(src = "pdf_icon.png", 
                                                            width = 16, 
                                                            height = 20),
                                                        href = "04SATB2_Quick_Facts.pdf",
                                                        target = '_blank'), "Quick facts"),
                                                  br(), br(),
                                                  div(style = "display: inline-block; width: 40%;", 
                                                      width = 12, shiny::a(
                                                        img(src = "pdf_icon.png", 
                                                            width = 16, 
                                                            height = 20),
                                                        href = "05SATB2_Information_Sheet_Families.pdf",
                                                        target = '_blank'), "Information Sheet for Families"),
                                                  div(style = "display: inline-block; width: 40%;", 
                                                      width = 12, shiny::a(
                                                        img(src = "pdf_icon.png", 
                                                            width = 16, 
                                                            height = 20),
                                                        href = "06SATB2_Information_sheet_Provider.pdf",
                                                        target = '_blank'), "Information Sheet for Providers"),
                                                  br(), br(),
                                                  div(style = "display: inline-block; width: 40%;", 
                                                      width = 12, shiny::a(
                                                        img(src = "pdf_icon.png", 
                                                            width = 16, 
                                                            height = 20),
                                                        href = "07SATB2_Infographic_Behavior.pdf",
                                                        target = '_blank'), "Infographic Behavior"),
                                                  div(style = "display: inline-block; width: 40%;", 
                                                      width = 12, shiny::a(
                                                        img(src = "pdf_icon.png", 
                                                            width = 16, 
                                                            height = 20),
                                                        href = "08SATB2_Infographic_Dental.pdf",
                                                        target = '_blank'), "Infographic Dental"),
                                                  br(), br(),
                                                  div(style = "display: inline-block; width: 40%;", 
                                                      width = 12, shiny::a(
                                                        img(src = "pdf_icon.png", 
                                                            width = 16, 
                                                            height = 20),
                                                        href = "09SATB2_Infographic_Speech.pdf",
                                                        target = '_blank'), "Infographic Speech"),
                                                  div(style = "display: inline-block; width: 40%;", 
                                                      width = 12, shiny::a(
                                                        img(src = "pdf_icon.png", 
                                                            width = 16, 
                                                            height = 20),
                                                        href = "11Final_Dental_letter_for_SATB2_professionals.pdf",
                                                        target = '_blank'), "Dental letter for ", em("SATB2 "), "Professionals"),
                                                  br(), br(),
                                                  div(style = "display: inline-block; width: 40%;", 
                                                      width = 12, shiny::a(
                                                        img(src = "pdf_icon.png", 
                                                            width = 16, 
                                                            height = 20),
                                                        href = "10Final_Speech_letter_for_SATB2_Families_and_Professionals.pdf",
                                                        target = '_blank'), "Speech letter for ", em("SATB2 "), "Families and Professionals")
                                                  #div(h2(tags$i(class = "far fa-file-pdf", style = "color:#676767")), "Background Information", shiny::a("PDF", href = "03Background_info.pdf", download = NA, target = "_blank")),
                                                  #div("Quick facts", shiny::a("PDF", href = "04SATB2_Quick_Facts.pdf")),
                                                  #div("Information Sheet for Families", shiny::a("PDF", href = "05SATB2_Information_Sheet_Families.pdf")),
                                                  #div("Information Sheet for Providers", shiny::a("PDF", href = "06SATB2_Information_sheet_Provider.pdf")),
                                                  #div("Infographic Behavior", shiny::a("PDF", href = "07SATB2_Infographic_Behavior.pdf")),
                                                  #div("Infographic Dental", shiny::a("PDF", href = "08SATB2_Infographic_Dental.pdf")),
                                                  #div("Infographic Speech", shiny::a("PDF", href = "09SATB2_Infographic_Speech.pdf")),
                                                  #div("Final Speech letter for SATB2 Families and Professionals", shiny::a("PDF", href = "10Final_Speech_letter_for_SATB2_Families_and_Professionals.pdf")),
                                                  #div("Final Dental letter for SATB2 Professionals", shiny::a("PDF", href = "11Final_Dental_letter_for_SATB2_professionals.pdf")),
                                                  # p(HTML(paste0("<li> Background Information PDF</li>"))),
                                                  # p(HTML(paste0("<li> Quick facts PDF</li>"))),
                                                  # p(HTML(paste0("<li> Information Sheet for Families PDF</li>"))),
                                                  # p(HTML(paste0("<li> Information Sheet for Providers PDF</li>"))),
                                                  # p(HTML(paste0("<li> Infographic Behavior PDF</li>"))),
                                                  # p(HTML(paste0("<li> Infographic Dental PDF</li>"))),
                                                  # p(HTML("<li> Infographic Speech</li>")),
                                                  # p(HTML("<li> Final Speech letter for SATB2 Families and Professionals</li>")),
                                                  # p(HTML("<li> Final Dental letter for SATB2 Professionals</li>"))
                                                  ),
                                           ),
                                         
                                         # fluidRow(
                                         #   column(12, align = "center",
                                         #          offset = 0,
                                         #          div(style = "display: inline-block; width: 40%",
                                         #              width = 12, shiny::a(
                                         #                img(src = "SATB2_LOGO_1.jpg", width = 300, height = 300),
                                         #                href = "03Background_info.pdf",
                                         #                target = '_blank')),
                                         #          div(style = "display: inline-block; width: 40%",
                                         #              width = 12, shiny::a(
                                         #                img(src = "SATB2_LOGO_2_Europe.jpg", width = 300, height = 300),
                                         #                href = "https://www.satb2europe.org/",
                                         #                target = '_blank'))
                                         #   ))
                                   ))
                        ))))
         
), # end families tab

        # ##### VARIANT ANALYSIS #####

        tabPanel(title = "Variant Analysis", value = "variantTab",

         fluidRow(
           column(10, offset = 1,
                  panel(heading = "Analyse your variants", status = "info",
           fluidRow(
             column(12,style='padding:30px;',
               fluidRow(
                 panel(status="default", heading = "Enter variant",
                   column(3,
                     pickerInput(
                         inputId = 'var_gene',
                         label =  h5(strong(var_possible_genes_title)),
                         width = "100%",
                         choices = var_possible_genes,
                         options = list(`style` = "default")
                     )
                   ),
                   column(3,
                     numericInputIcon(
                         inputId = "search_cDNA_pos",
                         label = h5(strong("cDNA Position")),
                         min = 1,
                         max = 10000,
                         value = 1166,
                         width = "100%"
                     ),
                     pickerInput(
                         inputId = "search_Allele",
                         label = "Ref",
                         choices = c("G", "A", "C", "T"),
                         selected = "G"
                     ),
                     pickerInput(
                         inputId = "search_cDNA_alt",
                         label = "Alt",
                         choices = c("G", "A", "C", "T", "null"),
                         selected = "A"
                     ),
                     actionButton(inputId = "search_var_c", label = "Search")),
                   column(3,
                     numericInputIcon(
                         inputId = "search_AA_pos",
                         label = h5(strong("Amino Acid Position")),
                         min = 1,
                         max = 10000,
                         value = 389,
                         width = "100%"
                     ),
                     pickerInput(
                         inputId = "search_AA_ref",
                         label = "Ref",
                         choices = sort(unique(master.df$AA_ref)),
                         selected = "Arg"
                     ),
                     pickerInput(
                         inputId = "search_AA_alt",
                         label = "Alt",
                         choices = sort(unique(master.df$AA_alt)),
                         selected = "Cys"
                     ),
                     actionButton(inputId = "search_var_p", label = "Search")),
                   column(2,
                     pickerInput(
                       inputId = "get_var_type",
                       label = h5(strong("Variant Consequence")),
                       choices = "Missense",
                       selected = "Missense"))
           )))),
           fluidRow(
             column(12,style='padding:30px;',
               fluidRow(
                 panel(status="default",
                       heading = "Variant Information",
                   fluidRow(
                     div(width = "100%",valueBoxOutput("geneBox1")),
                     div(width = "100%",valueBoxOutput("geneBox2")),
                     div(width = "100%",valueBoxOutput("geneBox3"))
                   ),
                   ##Clinical significance according to ClinVar
                   div(h4("Clinical significance according to ClinVar",
                          tippy(icon("question-circle"),
                                tooltip = h6(HTML(paste0("ClinVar version: 14th July 2021")),
                                             align = "left"),
                                animation = "scale",
                                theme = "light"))),
                   fluidRow(
                     #tags$head(tags$style(HTML("a {color: black}"))),
                     #column(4,
                     div(width = "100%",valueBoxOutput("ClinVarbox")))
                 ),
           ))),
           fluidRow(
             column(12,style='padding:30px;',
               fluidRow(
                 panel(status="default",
                       heading = "Patient information",
                       box(title = var_patient_info_title, width = 12,
                       DT::dataTableOutput(outputId = "patientTable"),
                       br(), p(var_patient_info_abb, style=sub_style, align = "center")))#,
             ))),
           fluidRow(
             column(12,style='padding:30px;',
              fluidRow(
                panel(status="default",
                      heading = "Custom variant analysis",
                  tabsetPanel(
                    tabPanel("Comparative Information",
                             br(),
                             h4("Compare the selected variant with other similar variants."),
                             br(),
                             radioGroupButtons(inputId = "compareButtons",
                                               label = "Variants with the same:",
                                               choices = c("Amino Acid Position", "Protein region", "Variant Type"),
                                               justified = TRUE,
                                               status = "default",
                                               checkIcon = list(yes = icon("ok",lib = "glyphicon"))
                             ),
                    br(),
                    div(width = "100%", plotlyOutput("comparePlot")),
                    br(),
                    br(),
                    DT::dataTableOutput(outputId = "compareTable")
                  )#,
#                   tabPanel("In silico scores",
#                            br(),
#                            #column(5,plotOutput(outputId = "paraz_legend", height = 30, width = 450)),
#                            #column(12,div(width = "100%", plotlyOutput("paraz_legend"))),
#                            column(4,
#                                   "Paralog conservation score",
#                                   align = "center",
#                                   div(width = "100%", plotlyOutput("Var_analyis_paraz"))
#                            ),
#                            column(4,
#                                   "Missense constraint score",
#                                   align = "center",
#                                   div(width = "100%", plotlyOutput("Var_analyis_mtr"))
#                            ),
#                            column(4,
#                                   "Pathogenic variant enrichment",
#                                   align = "center",
#                                   div(width = "100%", plotlyOutput("Var_analyis_per"))
#                            ),
#                            column(12,plotOutput(outputId = "paraz_legend", height = 30, width = 600))

                  )
  ))
             ))
                  )))
         ), # end variant analysis tab

##### RESEARCH #####

tabPanel(title = "Genotype-Phenotype", value = "researchTab", # title name of top menu
         fluidRow(
           column(10, offset = 1,
                  panel(heading = "Analyse your variants", status = "danger",
                        fluidRow(
                          column(12,style='padding:15px;',
                                 panel(heading="Filter Registry", status="default",
                                       fluidRow(
                                         column(12,
                                                selectizeGroupUI(
                                                  id = "research-filters",
                                                  btn_label = "Reset filters",
                                                  params = list(
                                                    varcons = list(
                                                      inputId = "Vartype",
                                                      title = p(strong("Variant Type")),
                                                      placeholder="all",
                                                      multiple = TRUE,
                                                      choices = unique(Patient_data.df$Vartype)
                                                    ),
                                                    aachange = list(
                                                      inputId = "AA_alt",
                                                      title = p(strong("Amino Acid Change")),
                                                      placeholder="all",
                                                      multiple = TRUE,
                                                      choices = unique(Patient_data.df$AA_alt)
                                                    ),
                                                    domain = list(
                                                      inputId = "Domain",
                                                      title = p(strong("Protein region")),
                                                      placeholder="all",
                                                      multiple = TRUE,
                                                      choices = unique(Patient_data.df$Domain)
                                                    ),
                                                    Clinical_seizures = list(
                                                      inputId = "Clinical_seizures",
                                                      title = p(strong("Seizures")),
                                                      placeholder="all",
                                                      choices = unique(Patient_data.df$Clinical_seizures)
                                                    ),
                                                    Cleft_palate = list(
                                                      inputId = "Cleft_palate",
                                                      title = p(strong("Cleft palate")),
                                                      placeholder="all",
                                                      choices = unique(Patient_data.df$Cleft_palate)
                                                    ),
                                                    Total_speech = list(
                                                      inputId = "Total_speech",
                                                      title = p(strong("Total speech")),
                                                      placeholder="all",
                                                      choices = unique(Patient_data.df$Total_speech)
                                                    ),
                                                    Abnormal_brainMRI= list(
                                                      inputId = "Abnormal_brainMRI",
                                                      title = p(strong("Abnormal Brain MRI"),
                                                                placeholder="all",
                                                                choices = unique(Patient_data.df$Abnormal_brainMRI),
                                                                multiple = TRUE
                                                      ))
                                                  )),
                                         ))))),
                        fluidRow(
                          column(12,style='padding:15px;',
                                 panel(status = "default",heading = "Custom variant exploration",
                                       div(textOutput(outputId = "filtered_n"),
                                           br(),
                                           "  "),
                                       tabsetPanel(
                                         tabPanel(
                                           "Genotype Interface",
                                           fluidRow(
                                             column(
                                               12,
                                               br(),
                                               p(
                                                 strong("Selected variants are displayed in 2D (lolliplot) and 3D (protein structure).")
                                               ),
                                               fluidRow(
                                                 column(7,
                                                        materialSwitch(
                                                          inputId = "gnomad_m",
                                                          label = "Reference population missense variants (gnomAD)",
                                                          status = "primary",
                                                          right = T,
                                                          inline = T
                                                        ),
                                                 )),
                                               fluidRow(
                                                 column(5,plotOutput(outputId = "Genotype_legend_plot", height = 30, width = 450) ),
                                               ),
                                               fluidRow(
                                                 column(7,
                                                        addSpinner(plotlyOutput(outputId = "Genotype_overview_plot"), color =
                                                                     spinner_color),
                                                        br(),
                                                        p(research_geno_transcripts , align="center", style=sub_style),
                                                        br()
                                                 ),
                                                 column(5,
                                                        addSpinner(color = spinner_color,
                                                                   r3dmolOutput(
                                                                     outputId = "threeDmolGene_all",
                                                                     width = "100%",
                                                                     height = "400px"
                                                                   )),
                                                        div("UniProt:",
                                                            align="center",
                                                            style=sub_style,
                                                            shiny::a("Q9UPW6", href = "https://www.uniprot.org/uniprot/Q9UPW6", target = "_blank"),
                                                            br(),
                                                            div("Predicted structure", align = "center", style=sub_style,
                                                                shiny::a("Alpha fold", href = "https://www.alphafold.ebi.ac.uk/entry/Q9UPW6", target = "_blank"))))
                                               )),
                                             fluidRow(
                                               column(12, 
                                                      panel(heading = strong("Display variants"), status = "default",
                                                            DT::dataTableOutput(outputId = "compareTableResearch")
                                                      ))
                                             )
                                           )),
                                         
                                         tabPanel("Phenotype Interface",
                                                  br(),
                                                  fluidRow(
                                                    column(12, align = "justify", plotlyOutput("research_phenotype1"))),
                                                  fluidRow(
                                                    column(6,align="justify", plotlyOutput("research_phenotype2")),
                                                    column(6,align="justify", plotlyOutput("research_phenotype3"))
                                                  ),
                                                  fluidRow(
                                                    column(6,align="justify", plotlyOutput("research_phenotype4")),
                                                    column(6,align="justify", plotlyOutput("research_phenotype5"))
                                                  ),
                                                  fluidRow(
                                                    column(6,align="justify", plotlyOutput("research_phenotype6"))
                                                  ),
                                                  fluidRow(
                                                    column(12, 
                                                           panel(heading = strong("Display variants"), status = "default",
                                                                 DT::dataTableOutput(outputId = "compareTableResearch2")
                                                                )
                                                           )
                                                          )
                                                  )
           ))))
                  )))
      ), #end research tab

##### SATB2 REGISTRY #####

tabPanel(title = p(HTML(paste0(em("SATB2"), " Registry"))), value = "registryTab", # change font size -> style = "font-size:16px;", 
         fluidRow(column(10, offset=1,
                         panel(heading = p(HTML(paste0("How to enter the ", em("SATB2"), " Registry"))), status="primary",
                               fluidRow(
                                 panel(column(8,
                                              #style = 'padding:30px;', 
                                              panel(heading = p(HTML(paste0(em("SATB2"), " Registry information"))),
                                                      "This large worldwide registry for SATB2-Associated syndrome (SAS) is led by Dr. Yuri Zarate and Katie Bosanko. 
                                                      After you contact the team, the participation consent is sent and reviewed for questions. Next, a link to the
                                                      database with SAS-specific questions is sent to the family to complete. The registry team collects medical records 
                                                      from providers and provides a high level of data curation. The registry has provided a great deal of information 
                                                      that is published in the literature as it is analyzed.", 
                                                    br(), br(), 
                                                    p(HTML(paste0("Arkansas Children’s Hospital hosts the only dedicated multidisciplinary clinic for the ", em("SATB2"), "-Associated 
                                                    syndrome in the United States. This clinic is hosted quarterly and participants can anticipate undergoing a comprehensive evaluation
                                                    by our experienced ", em("SATB2"), " team. More information at ", em("SATB2"), " International Clinic (", 
                                                                  shiny::a("archildrens.org", 
                                                                           href = "https://www.archildrens.org/programs-and-services/satb2-international-clinic?&journey=symptoms", 
                                                                           target = "_blank"), ").")))
                                 )),
                                 column(4,align="center", img(src="SAS_registry_LOGO.jpg", width="95%"), 
                                        actionBttn(
                                          inputId="link2",
                                          label = shiny::a("Participate", href =  contact_katie, style="color: white;", target =
                                                             "_blank"),
                                          style = "simple",
                                          color = "warning",
                                          icon = icon("arrow-right"),
                                          size = "md",
                                          block = FALSE,
                                          no_outline = TRUE
                                        )))),
                         )))),
                               
                               # fluidRow(
                               #   panel(
                               #     column(8,
                               #            panel(
                               #              heading = "If you are from North or South America or Australia",
                               #              "please contact Jenifer Sargent (study coordinator) at:
                               #                Jenifer.Sargent@cuanschutz.edu. This chapter of the study is led by neurologists Drs. Tim Benke and Kristen Park
                               #                (University of Colorado) and Dr. Jennifer Bain (", shiny::a("Simons Foundation/Columbia University", href = "https://www.simonssearchlight.org/"
                               #                ),")",
                               #              
                               #              br(),
                               #              br(),
                               #              "1. Email Jenifer Sargent your de-identified genetic testing results to make sure you qualify.",
                               #              br(),
                               #              "2. After we review your genetic testing results, Jenifer will send you a consent form; sign and return to her.",
                               #              br(),
                               #              "3. Jenifer will send you the clinical questionnaire. We are updating the questionnaire and want you to fill out the latest version, even if you've already done it before.",
                               #              br(),
                               #              "4. Fill out the questionnaire; you can ask your doctor to help fill it out.",
                               #              br(),
                               #              "5. Return the filled-out questionnaire to Jenifer; Jenifer will contact CFERV to conduct functional studies for novel variants.",
                               #              br(),
                               #              "6. In the near future, Jenifer will email you for yearly updates. If you haven't heard from us, then please contact us again."
                               #            ),
                               #            panel(
                               #              p("Both chapters of the registry will be merged, and patients only need to be enrolled at Colorado or Leipzig.
                               #                If you have already enrolled but have additional novel information, please get in contact and update your entry."
                               #              )
                               #            )
                               #     ),
                               #     column(4,
                               #            align="center", img(src="grin-world-1.PNG", width="95%"),
                               #            actionBttn(
                               #              inputId="link1",
                               #              label = shiny::a("Participate", 
                               #                               href="mailto:Jenifer.Sargent@cuanschutz.edu?subject=GRIN%20Registry&body=Dear%20Jenifer,%0D%0A%0D%0AWe%20would%20like%20to%20participate%20in%20the%20GRIN%20Registry.%20Please%20let%20us%20know%20how%20we%20can%20contribute.", 
                               #                               style="color:white;"
                               #              ),
                               #              icon = icon("arrow-right"),
                               #              style = "simple",
                               #              color = "warning",
                               #              size = "md",
                               #              block = F,
                               #              no_outline = TRUE
                               #            )))),
                               
                         #))),
         
#), # end GRIN Registry tab

##### ABOUT #####
tabPanel(title = "About",
         value = "aboutTab",
         fluidRow(
           column(10, offset = 1,
                  panel(heading = "About", status = "primary",
                        tabsetPanel(
                          tabPanel("General Information",
                                   panel(heading ="Portal version", status = "default",
                                         p(strong("This is the alpha version of the  SATB2 Portal")),
                                         p(HTML(paste0("The SATB2 Portal is the result of a collaborative effort to disseminate quality relevant data generated during the study of ", em("SATB2"), "-Associated syndrome."))),
                                         p("The goals of this project are: "),
                                         br(),
                                         fluidRow(
                                           column(6,
                                                  p(HTML(paste0("<ul><li> Summarize clinical and molecular data interactively</li>"))),
                                                  p(HTML(paste0("<li> Provide information on ", em("SATB2")," -Associated syndrome</li>"))),
                                                  p(HTML(paste0("<li> Support and promote research on ", em("SATB2"),"-Associated syndrome</li>"))),
                                                  p(HTML(paste0("<li> Facilitate ongoing recruitment to the ", em("SATB2")," -Associated syndrome registry</li>")))),
                                           column(6,
                                                  p(HTML(paste0("<ul><li> Aid in ", em("SATB2"), " variant interpretation and classification</li>"))),
                                                  p(HTML(paste0("<li> Visualizing data from the global ", em("SATB2")," -Associated syndrome registry</li>"))),
                                                  p(HTML("<li> Promote collaboration between clinicians and other researchers</li>")))),
                                         br(),
                                         footer = div("The SATB2 Portal is an ongoing project of the scientific community. Interested collaborators are invited to reach out to join the project.")),
                                   panel(heading = "Teams and People", status = "default",
                                         p("The current version of the SATB2 Portal has been developed by an international team of researchers and clinicians:"),
                                         fluidRow(
                                           column(10,
                                                  panel(heading = "Team Leaders",
                                                        p(strong("Yuri Zarate"), "(Little Rock, US): Clinical & genetic data"),
                                                        p(strong("Dennis Lal"), "(Cleveland, US): General concept, web development, bioinformatics, video production"),
                                                        # p(strong("Kimberly Goodspeed"), "(Dallas, US): Clinical & genetic data"),
                                                        # p(strong("Katty (Jing-Qiong) Kang"), "(Vanderbilt, US): Molecular data"),
                                                  )),
                                           column(12, p("")),
                                           column(3,
                                                  panel(style="height: 230px;",heading = "Clinical & Genetic Data",
                                                        div(style="height: 100%;",
                                                            p(strong("Yuri Zarate")))
                                                        # )),
                                                        # 
                                                        # column(2,
                                                        #   panel(style="height: 230px;",heading = "Molecular Data",
                                                        #         div(style="height: 100%",p(strong("Yuri Zarate")))
                                                  )),
                                           column(3,
                                                  panel(style="height: 230px;",heading = "Web Development",
                                                        div(style="height: 100%",p(strong("Arthur Stefanski")),
                                                            p(strong("Tobias Brünger")),
                                                            p("Eduardo Perez-Palma"),
                                                            p("Marie Macnee"),
                                                            p("Chiara Klöckner"))
                                                  )),
                                           column(3,
                                                  panel(style="height: 230px;",heading = "Bioinformatics",
                                                        div(style="height: 100%",
                                                            p(strong("Tobias Brünger")),
                                                            p("Eduardo Perez-Palma"),
                                                            p("Marie Mcnee"),
                                                            p("Patrick May"),
                                                            p("Chiara Klöckner")) #,
                                                        # p("Johannes Lemke"))
                                                        # )),
                                                        # column(2,
                                                        #   panel(style="height: 230px;", heading = "Video",
                                                        #         div(style="height: 100%",p(strong("Arthur Stefanski")),
                                                        #         p("Amber Freed"))
                                                  ))
                                         )),
                                   panel(heading = "Imprint",
                                         status = "default",
                                         p("We object to any commercial use and disclosure of data."),
                                         p(strong("Copyright and use:"), "The authors grants you the right of use to make a private copy for personal purposes.
                        However, you are not entitled to change and/or pass on the materials or publish them yourself.
                        Upon request, you will receive free information about the personal data stored about you.
                        To do so, please contact the administrator."),
                                         p(strong("No liability:"), "The contents of this web project have been carefully checked and created to
                        the best of our knowledge. But for the information presented here is no claim to completeness,
                        timeliness, quality and accuracy. No responsibility can be accepted for any damage caused by reliance
                        on or use of the contents of this website."))
                          ),
                          tabPanel("Terms and Data Information",
                                   panel(heading = "Terms of Use", status = "default",
                                         p(about_terms_of_use,
                                           shiny::a(href=contact_us, "Contact us"),
                                           "that we can improve.")),
                                   panel(heading = "Data generation information", status = "default",
                                         about_data_generation)
                          ),
                          tabPanel(title = "Tutorials", value = "tutorialTab", #### Tutorials
                                   fluidRow(
                                     column(10, offset = 1,
                                            panel(heading = "Tutorials (EXAMPLE)", status = "primary",
                                                  div(align = "center",
                                                      tabsetPanel(
                                                        tabPanel("Basic Information",
                                                                 column(8, offset = 2,
                                                                        box(title= "",
                                                                            align = "right",
                                                                            width = 12,
                                                                            embed_url(tut_basic_info) %>%
                                                                              use_rounded() %>% use_align("center") %>% use_bs_responsive()
                                                                        ))
                                                        ),
                                                        tabPanel("Families",
                                                                 column(8, offset = 2,
                                                                        box(title= "",
                                                                            align = "right",
                                                                            width = 12,
                                                                            embed_url(tut_family) %>%
                                                                              use_rounded() %>% use_align("center") %>% use_bs_responsive()
                                                                        ))
                                                        ),
                                                        tabPanel("Variant Analysis",
                                                                 column(8, offset = 2,
                                                                        box(title= "",
                                                                            align = "right",
                                                                            width = 12,
                                                                            embed_url(tut_variant_analysis) %>%
                                                                              use_rounded() %>% use_align("center") %>% use_bs_responsive()
                                                                        ))
                                                        ),
                                                        tabPanel("Genotype-Phenotype",
                                                                 column(8, offset = 2,
                                                                        box(title= "",
                                                                            align = "right",
                                                                            width = 12,
                                                                            embed_url(tut_research) %>%
                                                                              use_rounded() %>% use_align("center") %>% use_bs_responsive()
                                                                        )))
                                                      )))))
                          ))))))
)) # end ui