#Setup working directory
setwd("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/analysis_data")

#Installs and initializes necessary packages
if (!require(BiocManager)){
  install.packages(BiocManager)
}
if (!require(TCGAbiolinks)){
  install.packages(BiocManager)
}
if (!require(survival)){
  install.packages(BiocManager)
}
if (!require(survminer)){
  install.packages(BiocManager)
}
library(BiocManager)
library(TCGAbiolinks)
library(survival)
library(survminer)

#Re-queries, reads, and obtains clinical data
clinical_query <- GDCquery(project = "TCGA-BRCA",
                           data.category = "Clinical",
                           file.type = "xml")
clinical <- read.csv("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/analysis_data/brca_clinical_data.csv")
clinical.drug <- GDCprepare_clinic(query = clinical_query,
                                  clinical.info = "drug")
clinical.rad <- GDCprepare_clinic(query = clinical_query,
                                  clinical.info = "radiation")

#Checks for relative low N/A proportion in selected variables
sum(is.na(clinical$lymph_node_examined_count))
# Number of N/A: 139
sum(!is.na(clinical$lymph_node_examined_count))
# Number of valid values: 1035
sum(is.na(clinical$race_list))
# Number of N/A: 0
sum(!is.na(clinical$race_list))
# Number of valid values: 1174

#Filters unwanted rows (blank or N/A data) from the dataset
filter_mask1 = ifelse(is.numeric(clinical$lymph_node_examined_count), TRUE, FALSE)
clinical_filtered1 = clinical[filter_mask1, ]
filter_mask2 = ifelse(!is.na(clinical_filtered1$race_list), TRUE, FALSE)
clinical_filtered2 = clinical_filtered1[filter_mask2, ]
filter_mask3 = ifelse(clinical_filtered2$race_list != "", TRUE, FALSE)
clinical_filtered3 = clinical_filtered2[filter_mask3, ]

#Plots lymph_node_exampined_count vs. race_list
plot(factor(clinical_filtered3$race_list), clinical_filtered3$lymph_node_examined_count, xlab = "Race", ylab = "Lymph Node Examined Count")
#ggplot(data = clinical_filtered2, aes(x = race_list, y = lymph_node_examined_count))

#Creates a survival plot for race groups
clinical_filtered3$survival_time = ifelse(is.na(clinical_filtered3$days_to_death), clinical_filtered3$days_to_last_followup, clinical_filtered3$days_to_death)
clinical_filtered3$death_event = ifelse(clinical_filtered3$vital_status == "Dead", TRUE, FALSE)

surv_object_age <- Surv(time = clinical_filtered3$survival_time,
                        event = clinical_filtered3$death_event)

race_fit <- surv_fit(surv_object_age ~ clinical_filtered3$race_list,
                    data = clinical_filtered3 )

survplot_race = ggsurvplot(race_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")

KM_plot_race = survplot_race$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10))

KM_plot_race

#Creates a survival plot for lymph node examined counts by making a new grouped lymph node variable
#Stratification of lymph node variable
clinical_filtered3$lymphgroups = ifelse(clinical_filtered3$lymph_node_examined_count < 9, "0-8 Lymph Nodes Examined",
  ifelse(clinical_filtered3$lymph_node_examined_count < 18, "9-17 Lymph Nodes Examined",
  ifelse(clinical_filtered3$lymph_node_examined_count < 27, "18-26 Lymph Nodes Examined",
  ifelse(clinical_filtered3$lymph_node_examined_count < 36, "27-35 Lymph Nodes Examined", "36-44 Lymph Nodes Examined"
))))

#Graphs the survival plot for lymph node groups
clinical_filtered3$survival_time = ifelse(is.na(clinical_filtered3$days_to_death), clinical_filtered3$days_to_last_followup, clinical_filtered3$days_to_death)
clinical_filtered3$death_event = ifelse(clinical_filtered3$vital_status == "Dead", TRUE, FALSE)

surv_object_age <- Surv(time = clinical_filtered3$survival_time,
                        event = clinical_filtered3$death_event)

lymph_fit <- surv_fit(surv_object_age ~ clinical_filtered3$lymphgroups,
                     data = clinical_filtered3 )

survplot_lymph = ggsurvplot(lymph_fit, 
                           pval=TRUE, 
                           ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                           legend = "right")

KM_plot_lymph = survplot_lymph$plot + 
  theme_bw() +  #
  theme(axis.title = element_text(size=16), 
        axis.text = element_text(size=12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=10))

KM_plot_lymph

