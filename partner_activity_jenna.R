knitr::opts_knit$setwd(root.dir = normalizePath("/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/analysis_data")) 
clinical <- read.csv("/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/brca_clinical_data.csv")
# load all the packages and files required 
library(BiocManager)
library(TCGAbiolinks)
library(ggplot2)
GDCdownload(clinical_query)
clinical_query <- GDCquery(project="TCGA-BRCA", data.category="Clinical", file.type = "xml")
clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")
clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")

# Task1
# create a clean data frame with no NA
percent_mask <- ifelse(clinical$progesterone_receptor_level_cell_percent_category=="", F, T)
percent_cleaned <- clinical[percent_mask,]
# change the progesterone receptor level data from strings to integers
percent_cleaned$progesterone_receptor_level_cell_percent_category <- gsub("%","",as.character(percent_cleaned$progesterone_receptor_level_cell_percent_category))
percent_cleaned$progesterone_receptor_level_cell_percent_category <- gsub("<","",as.character(percent_cleaned$progesterone_receptor_level_cell_percent_category))
percent_cleaned$progesterone_receptor_level_cell_percent_category <- gsub("-","",as.character(percent_cleaned$progesterone_receptor_level_cell_percent_category))
percent_cleaned$progesterone_receptor_level_cell_percent_category <- strtoi(percent_cleaned$progesterone_receptor_level_cell_percent_category)
percent_cleaned$progesterone_receptor_level_cell_percent_category <- ifelse(percent_cleaned$progesterone_receptor_level_cell_percent_category==10, 10.1, percent_cleaned$progesterone_receptor_level_cell_percent_category%/%100)
# use boolean indexing to get catagorical variable "status"
status <- percent_cleaned$breast_carcinoma_progesterone_receptor_status
positive <- percent_cleaned[status=="Positive",]
negative <- percent_cleaned[status=="Negative",]

# plot a histogram with both positive and negative status vs. progesterone receptor level
jpeg(file="/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/histogram.jpg")
hist(positive$progesterone_receptor_level_cell_percent_category, 
     freq=F, xlim=c(0,100), ylim=c(0,0.12),
     col=rgb(0,0,1,1/4), main="Histogram of Progesterone Receptor Level Cell Percent",
     xlab="Percentage", cex.main=0.8)
hist(negative$progesterone_receptor_level_cell_percent_category, 
     freq=F, col=rgb(1,0,0,1/4), add=T)

labels <- c("positive", "negative")
legend("topright", legend=labels, cex=0.8, inset=0.01, pch=15, col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))
dev.off()

# Task 2
library(survival)
library(survminer)
# create survival and fit object
# make the continuous variable catagorical
percent_cat <- ifelse(percent_cleaned$progesterone_receptor_level_cell_percent_category<=33, "low", 
                      ifelse(percent_cleaned$progesterone_receptor_level_cell_percent_category<=66, "medium", "high"))
surv_object_percent <- Surv(time = percent_cleaned$survival_time,
                            event = percent_cleaned$death_event)
percent_fit <- surv_fit( surv_object_percent ~ percent_cat,
                         data = percent_cleaned)
# plot the KM plot
survplot_percent = ggsurvplot(percent_fit, 
                              pval=TRUE, 
                              ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                              legend = "right")
KM_plot_percent = survplot_percent$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=10), 
        axis.text = element_text(size=6),
        legend.title = element_text(size=10),
        legend.text = element_text(size=7))

jpeg(file="/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/KM1.jpg")
KM_plot_percent
dev.off()

# Task 3
# create survival and fit object
surv_object_status <- Surv(time = percent_cleaned$survival_time,
                           event = percent_cleaned$death_event)
status_fit <- surv_fit( surv_object_percent ~ status,
                        data = percent_cleaned)
# plot the KM plot
survplot_status = ggsurvplot(status_fit, 
                             pval=TRUE, 
                             ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), 
                             legend = "right")
KM_plot_status = survplot_status$plot + 
  theme_bw() +  
  theme(axis.title = element_text(size=10), 
        axis.text = element_text(size=6),
        legend.title = element_text(size=10),
        legend.text = element_text(size=7))
jpeg(file="/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/KM2.jpg")
KM_plot_status
dev.off()

# save all the figures and data
write.csv(percent_cleaned, "/Users/JennaJacobs/Documents/USC/qbio490/qbio_490_jenna/analysis_data/percent_cleaned.csv", row.names = FALSE)

