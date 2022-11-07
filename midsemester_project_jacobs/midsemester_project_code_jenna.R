
dir.create("/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/midsemester_project_jacobs/outputs")
knitr::opts_knit$set(root.dir = normalizePath("/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/midsemester_project_jacobs/outputs"))
 
if(!require(BiocManager)){ install.packages("BiocManager")
}
library(BiocManager)    
BiocManager::install("maftools")
library(maftools) 

if(!require(survival)) {BiocManager::install("survival")
}
library(survival) 

if(!require(TCGAbiolinks)) {BiocManager::install("TCGAbiolinks")
}
library(TCGAbiolinks)   

if(!require(survminer)) {BiocManager::install("survminer")
}
library(survminer)   

if(!require(ggplot2)) {BiocManager::install("ggplot2")
}
library(ggplot2)          
         
clinical_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
GDCdownload(clinical_query)
clinical <- GDCprepare_clinic(clinical_query, clinical.info = "patient")
#querying the clinical data
clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")


ethnicity_mask <- ifelse(clinical$ethnicity == "HISPANIC OR LATINO", TRUE, FALSE)
ethnicity_NA_mask <- ifelse(clinical$ethnicity == "", F, T)

ethnicity_cleaned_clinical <-clinical[ethnicity_NA_mask, ] #removes the NA values

colnames(clinical)[colnames(clinical) == "bcr_patient_barcode"] <- "Tumor_Sample_Barcode"
#changing the name of one of the columns in clinical data in order to use maf
#View(CPR_cleaned_clinical)

#creating a survival time column, and death event column
ethnicity_cleaned_clinical$survival_time <-ifelse(is.na(ethnicity_cleaned_clinical$days_to_death), ethnicity_cleaned_clinical$survival_time <-ethnicity_cleaned_clinical$days_to_last_followup, ethnicity_cleaned_clinical$survival_time <-ethnicity_cleaned_clinical$days_to_death)
ethnicity_cleaned_clinical$death_event <- ifelse(ethnicity_cleaned_clinical$vital_status == "Alive", ethnicity_cleaned_clinical$days_to_death <-FALSE, ethnicity_cleaned_clinical$death_event <- TRUE)

#creating survival and fit object
ethnicity_surv_object <- Surv(time = ethnicity_cleaned_clinical$survival_time, event = ethnicity_cleaned_clinical$death_event)
ethnicity_fit <- survfit(ethnicity_surv_object ~ ethnicity_cleaned_clinical$ethnicity, data = ethnicity_cleaned_clinical)

#creating the survival plot
survplot_ethnicity <- ggsurvplot(ethnicity_fit, pval = TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

#saving the survival plot
jpeg("/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/midsemester_project_jacobs/outputs/ethnicity_KM_Plot.jpg")
KM_Plot_ethnicity = survplot_ethnicity$plot + theme_bw() + theme(axis.title = element_text(size = 20), axis.text = element_text(size=16), legend.title = element_text(size = 14))
KM_Plot_ethnicity
dev.off()



#download and query maf 
maf_query <- GDCquery(
  project = "TCGA-BRCA", 
  data.category = "Simple Nucleotide Variation", 
  access = "open",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(maf_query) 
maf <- GDCprepare(maf_query)
maf_object <- read.maf(maf = maf, 
                       clinicalData = clinical,
                       isTCGA = TRUE)

#mask for only patients that are hispanic or latino, and create barcodes
His_Lat_clinical <- ethnicity_cleaned_clinical[ethnicity_mask, ]
HL_barcodes <- ethnicity_cleaned_clinical[ethnicity_cleaned_clinical$ethnicity == "HISPANIC OR LATINO", 1]
HL_maf <- subsetMaf(maf = maf_object, tsb = HL_barcodes)

#mask for patients that are not hispanic or latino, and create barcodes
NOT_His_Lat_clinical <- ethnicity_cleaned_clinical[!ethnicity_mask, ]
NotHL_barcodes <- ethnicity_cleaned_clinical[ethnicity_cleaned_clinical$ethnicity == "NOT HISPANIC OR LATINO", 1]
NotHL_maf <- subsetMaf(maf = maf_object, tsb = NotHL_barcodes)

#par(mar=c(2,2,2,2))
#jpeg("/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/midsemester_project_jacobs/outputs/ethnicity_coonco.jpg")

#coOncoplot(m1 = positive_maf,
           #m2 = negative_maf,
           #m1Name = 'Hispanic or Latino',
           #m2Name = 'Not Hispanic or Latino'
           #)


jpeg("/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/midsemester_project_jacobs/outputs/his_lat_onco.jpg")
vector <- c("TP53", "PIK3CA", "TTN", "CDH1", "MAP3K1", "MUC16", "GATA3")
oncoplot(maf = HL_maf, genes = vector)
dev.off()

jpeg("/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/midsemester_project_jacobs/outputs/Not_his_lat_onco.jpg")
vector <- c("TP53", "PIK3CA", "TTN", "CDH1", "MAP3K1", "MUC16", "GATA3")
oncoplot(maf = NotHL_maf, genes = vector)
dev.off()

View(clinical.rad)

#creating and saving the lollipop plots
jpeg("/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/midsemester_project_jacobs/outputs/ethnicity_lollipop_PIK3CA.jpg")
lollipopPlot2(m1 = HL_maf, m2 = NotHL_maf, m1_name = 'Hispanic or Latino',
              m2_name = 'Not Hispanic or Latino', gene = "PIK3CA")
dev.off()
jpeg("/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/midsemester_project_jacobs/outputs/ethnicity_lollipop_BRCA1.jpg")
lollipopPlot2(m1 = HL_maf, m2 = NotHL_maf, m1_name = 'Hispanic or Latino',
              m2_name = 'Not Hispanic or Latino', gene = "BRCA1")
dev.off()

radiation_mask <- ifelse(clinical.rad$radiation_dosage == "", F, T)
boxplot_clinical<- clinical.rad[radiation_mask, ]

jpeg("/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/midsemester_project_jacobs/outputs/boxplot_clinical_radiation.jpg")
boxplot(boxplot_clinical$radiation_dosage, xlab = "Radiation Dosage")
dev.off()



#Saving CSV files of created datasets
write.csv(ethnicity_cleaned_clinical, "/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/midsemester_project_jacobs/outputs/ethnicity_cleaned_clinical.csv")
write.csv(clinical.rad, "/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/midsemester_project_jacobs/outputs/clinical.rad.csv")
write.csv(boxplot_clinical, "/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/midsemester_project_jacobs/outputs/boxplot_clinical.csv")
