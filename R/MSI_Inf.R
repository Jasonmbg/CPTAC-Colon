## For MSI status perform prediction for the CPTAC cohort with additional tools to evaluate the overall concordance within the ACCC-CRC cohort

# We start with the saved and normalized expression matrix; the first approach would be based on the RNA-Seq data to infer MSI status

library(PreMSIm)
library(tidyverse)

dr_here()

# this to provide the full path to the normalized gene expression data to reproduce the MSI inference
pt <- here("Data/CPTAC.Colon.Norm.Filt.Exp.ACCCPaper.MSI.tsv")

# also the same object in rda format: CPTAC.Colon.Norm.Filt.Exp.ACCCPaper.MSI.rda

input_cptac_colon = data_pre(pt, type = "Symbol")

Out.MSI.CPTAC.Colon <- msi_pre(input_cptac_colon) # Prediction results  ("1" indicates MSI-high, and "0" MSI-low/microsatellite stability).

# write_tsv(Out.MSI.CPTAC.Colon, file="MSIStatus.PreMSIm.Predict.CPTAC.Colon.NormFilt.tsv") 

# continue to merge the files;

clinical.cptac.final <- read_tsv(here("Data","CPTAC_clinical_final_mod_INFO.tsv"))
table(clinical.cptac.final$`MSI Status`)

combo.MSI.eval.dat.CPTAC <- inner_join(clinical.cptac.final, Out.MSI.CPTAC.Colon, 
by=c("Tumor_Sample_Barcode"="Sample")) %>% 
dplyr::select(Tumor_Sample_Barcode, `MSI Status`, MSI_status)

# write_tsv(combo.MSI.eval.dat.CPTAC, file="CPTAC.Colon.PreMSIm.Eval.MSI.tsv")

############ Create confusion matrix for CPTAC MSI inference #########################
# extra respective sources of explaining more the terms and relative metrics
# library(caret); -> confusionMatrix
# https://www.statology.org/confusion-matrix-in-r/

premsim <- as.character(combo.MSI.eval.dat.CPTAC$MSI_status)
premsim.refined <- recode(premsim, "0"="MSS", "1"="MSI-H")
combo.MSI.eval.dat.CPTAC$MSI_status <- premsim.refined

cm <- table(combo.MSI.eval.dat.CPTAC$MSI_status, combo.MSI.eval.dat.CPTAC$`MSI Status`)

accuracy <- sum(cm[1], cm[4]) / sum(cm[1:4]) 
precision <- cm[4] / sum(cm[4], cm[2])
sensitivity <- cm[4] / sum(cm[4], cm[3])
fscore <- (2 * (sensitivity * precision))/(sensitivity + precision)
specificity <- cm[1] / sum(cm[1], cm[2])

######################################################################################