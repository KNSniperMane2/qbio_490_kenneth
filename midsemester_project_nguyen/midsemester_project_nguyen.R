---
  title: "midsemester_project_nguyen"
output: html_notebook
---

#directory initialization
mkdir("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/midsemester_project_nguyen")
dir.create("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/midsemester_project_nguyen/outputs")
setwd("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/midsemester_project_nguyen/outputs")

#package/library initialization
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if(!require(TCGAbiolinks)) {BiocManager::install("TCGAbiolinks")}
if (!require(maftools)) {BiocManager::install("maftools")}
if (!require(ggplot2)) {install.packages(ggplot2)}
if (!require(SummarizedExperiment)) {install.packages(SummarizedExperiment)}
if (!require(survival)) {install.packages("survival")}
if (!require(survminer)) {install.packages("survminer")}
library(BiocManager)
library(maftools)
library(TCGAbiolinks)
library(ggplot2)
library(SummarizedExperiment)
library(survival)
library(survminer)

#querying
clinical_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
GDCdownload(clinical_query) 
clinical <- GDCprepare_clinic(query = clinical_query, clinical.info = "patient")
clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")
clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")
colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", # we only have access to somatic mutations which are open access
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(maf_query)
maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)
rna_query <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(rna_query)
rna_se <- GDCprepare(rna_query)
load("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/clinical_maf_rna.RData")

#loading necessary files (querying alternative due to query issues, can skip chunk)
#clinical_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
#GDCdownload(clinical_query) 
#clinical <- GDCprepare_clinic(query = clinical_query, clinical.info = "patient")
#colnames(clinical)[ colnames(clinical) == "bcr_patient_barcode" ] <- "Tumor_Sample_Barcode"
#clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")
#clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")
#maf_object <- read.maf(maf = "C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/analysis_data/brca_maf_data.csv_maftools.maf", clinicalData = clinical, isTCGA = TRUE)
#load("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/clinical_maf_rna.RData")
#rna_clinical <- read.csv("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/analysis_data/brca_rna_clincial_data.csv", header = TRUE)
#rna_counts <- read.csv("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/analysis_data/brca_rna_count_data.csv", header = TRUE)
#rna_genes <- read.csv("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/analysis_data/brca_rna_gene_data.csv", header = TRUE)

###Gene Counts Plots
#initializes rna_clinical
age_mask <-  is.na(rna_se@colData$age_at_index)
rna_clinical <-  rna_se@colData[!age_mask, ]
rna_clinical <- as.data.frame(rna_clinical)
treatments_mask <- ifelse(colnames(rna_clinical) == "treatments" | colnames(rna_clinical) == "disease_type" | colnames(rna_clinical) == "primary_site", FALSE, TRUE)
rna_clinical <- rna_clinical[, treatments_mask] 
#initializes rna_counts
rna_counts <- rna_se@assays@data$unstranded[, !age_mask]
#initializes rna_genes
rna_genes <- rna_se@rowRanges@elementMetadata
rna_genes <- as.data.frame(rna_genes)
#sets row names for rna_clinical
row.names(rna_clinical) = rna_clinical$barcode
#makes new age column
rna_clinical$age_category = ifelse(rna_clinical$age_at_index < 50, "young", "old")
#sets row names for rna_genes
row.names(rna_genes) = rna_genes$gene_id
#sets row/column names for rna_genes
colnames(rna_counts) = rna_clinical$barcode
row.names(rna_counts) = rna_genes$gene_id
#subsets out samples with normal tissue type
unique(rna_clinical$definition)
sample_mask <- ifelse(rna_clinical$definition == "Solid Tissue Normal", FALSE, TRUE)
rna_clinical <- rna_clinical[sample_mask, ]
rna_counts <- rna_counts[, sample_mask]
#analyzes 5 year survival, accounting for NA values and lack of death data
five_yr_death <- ifelse(rna_clinical$days_to_death == "NA", NA, ifelse(rna_clinical$days_to_death >= (365.25 * 5), TRUE, FALSE))
five_yr_death_and_followup <- ifelse(is.na(five_yr_death), ifelse(rna_clinical$days_to_last_follow_up > 365.25*5, TRUE, FALSE), five_yr_death)
rna_clinical$five_year_surv <- five_yr_death_and_followup
#masks out all but a desired gene
TP53_mask <- ifelse(rna_genes$gene_name == "TP53", TRUE, FALSE) # create mask based on the gene name (your mask should be true when the name is 'geneA' and false for any other gene name)
MKI67_mask <- ifelse(rna_genes$gene_name == "MKI67", TRUE, FALSE)
#gets gene ensemble IDs
TP53_ensembl <- rna_genes$gene_id[TP53_mask] # think about which column we should apply out mask to
MKI67_ensembl <- rna_genes$gene_id[MKI67_mask]
#creates list with counts for each gene
TP53_counts <- as.numeric(rna_counts[TP53_mask, ])
MKI67_counts <- as.numeric(rna_counts[MKI67_mask, ])
#gets summary of counts
summary(TP53_counts)
summary(MKI67_counts)

###plots MKI67 counts vs. TP53 counts
jpeg("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/midsemester_project_nguyen/outputs/TP53_MKI67_counts.jpg")
par(mar=c(4,4,4,4))
plot(TP53_counts, MKI67_counts, xlab = "TP53 Counts", ylab = "MKI67 Counts", main = "MKI67 vs. TP53 counts")
abline(lm(TP53_counts ~ MKI67_counts))
dev.off()

###plots 5 year survival for each gene (TP53 and PTGS2)
jpeg("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/midsemester_project_nguyen/outputs/TP53_five_year_surv.jpg")
par(mar=c(4,4,4,4))
boxplot(TP53_counts ~ rna_clinical$five_year_surv, xlab = "5 Year Survival", ylab = "TP53 Counts")
dev.off()
jpeg("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/midsemester_project_nguyen/outputs/MKI67_five_year_surv.jpg")
par(mar=c(4,4,4,4))
boxplot(MKI67_counts ~ rna_clinical$five_year_surv, xlab = "5 Year Survival", ylab = "MKI67 Counts")
dev.off()

###Draftsman plot
jpeg("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/midsemester_project_nguyen/outputs/Draftsman.jpg")
PTGS2_mask <- ifelse(rna_genes$gene_name == "PTGS2", TRUE, FALSE)
TTN_mask <- ifelse(rna_genes$gene_name == "TTN", TRUE, FALSE)
PTGS2_counts <- as.numeric(rna_counts[PTGS2_mask, ])
TTN_counts <- as.numeric(rna_counts[TTN_mask, ])
geneABCD_counts <- data.frame(TP53_counts, TTN_counts, MKI67_counts, PTGS2_counts)
colnames(geneABCD_counts) <- c("TP53", "TTN", "MKI67", "PTGS2")
cols <- character(nrow(rna_clinical)) 
cols[rna_clinical$five_year_surv == TRUE] <- "blue"
cols[rna_clinical$five_year_surv == FALSE] <- "red"
pairs(geneABCD_counts, col = cols, lower.panel=NULL)
dev.off()

###Oncoplot from maf data
jpeg("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/midsemester_project_nguyen/outputs/oncoplot.jpg")
oncoplot(maf = maf_object, genes = c("TP53", "MKI67"))
dev.off()

###KM plot
jpeg("C:/Users/kenne/OneDrive/Desktop/QBIO490/qbio_490_kenneth/midsemester_project_nguyen/outputs/KMplot.jpg")
filter_mask = ifelse((clinical$radiation_therapy == "NO" | clinical$radiation_therapy == "YES"), TRUE, FALSE)
clinical_filtered = clinical[filter_mask, ]
clinical_filtered$survival_time = ifelse(is.na(clinical_filtered$days_to_death), clinical_filtered$days_to_last_followup, clinical_filtered$days_to_death)
clinical_filtered$death_event = ifelse(clinical_filtered$vital_status == "Dead", TRUE, FALSE)
surv_object_rad <- Surv(time = clinical_filtered$survival_time,
                        event = clinical_filtered$death_event)
#creates a fit object for the KM plot
rad_fit <- surv_fit(surv_object_rad ~ clinical_filtered$radiation_therapy, data = clinical_filtered)
#formats KM plot
survplot_rad = ggsurvplot(rad_fit, 
                          pval=TRUE, 
                          ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                          legend = "right")
KM_plot_rad = survplot_rad$plot + 
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))
KM_plot_rad
dev.off()

###clinical.rad radiation data
summary(as.numeric(levels(clinical.rad_cleaned$radiation_dosage))[clinical.rad$radiation_dosage])


