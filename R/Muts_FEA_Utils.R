## TRY TO IMPLEMENT A MUTATIONAL ANALYSIS ON THE MAF COHORT, BASED ON THE TOP-SCORING GENES;
## THE ULTIMATE GOAL IS TO WORK WITH maftools R PACKAGE, AIMING TO MATCH THE UPDATED COORDINATES TO THE MATCHED SAMPLES

library(here)
library(tidyverse)
library(clusterProfiler)
library(ReactomePA)
library(data.table)
library(DOSE)
library(ggplot2)

dr_here()

######################### start with the MSI samples ##########################################

csv.files <- list.files(here("Mut_Freq/MSI_CPTAC"),pattern = "csv")

DT <- rbindlist(sapply(here("Mut_Freq","MSI_CPTAC",csv.files),fread,
                       simplify=FALSE,
                       select=c("chrom","pos", "ref_base", "alt_base", "achange", 
                                "all_mappings", "transcript", "hugo", "so", "cchange",
                                "Final_Score", "samples"),
                       verbose=getOption("datatable.verbose", TRUE)),
                use.names= TRUE,
                idcol="file")  %>%
  as_tibble()  %>%
  dplyr::filter(Final_Score >= 0.5) %>%
  dplyr::select(chrom, pos, ref_base, alt_base, samples, hugo, so, Final_Score) %>% 
  dplyr::filter(so %in% c("frameshift_elongation","frameshift_truncation","inframe_deletion", 
  "missense_variant", "splice_site_variant","start_lost","stop_gained","stop_lost",
  "inframe_insertion"))

# to break alterations occurring in multiple samples into separate lines
DT_2 <- DT %>% separate_rows(samples, sep=";") 

group_size_MSI <- length(unique(DT_2$samples))
group_size_MSI # check

msi_cptac_dat <- DT_2 %>% select(hugo, samples) %>% distinct(hugo,samples) %>% 
group_by(hugo) %>% dplyr::count(hugo,sort=T,name="MSI_sample_Freq") # make the count

freq_genes_MSI_CPTAC_dat <- msi_cptac_dat %>% filter(MSI_sample_Freq>=4)
freq_genes_MSI_CPTAC <- msi_cptac_dat %>% filter(MSI_sample_Freq>=4) %>% pull(hugo)

############################ also regarding MSS samples #####################################

csv.files <- list.files(here("Mut_Freq/MSS_CPTAC"),pattern = "csv")

DT <- rbindlist(sapply(here("Mut_Freq","MSS_CPTAC",csv.files),fread,
                       simplify=FALSE,
                       select=c("chrom","pos", "ref_base", "alt_base", "achange", "all_mappings", 
                                "transcript", "hugo", "so", "cchange",
                                "Final_Score", "samples"),
                       verbose=getOption("datatable.verbose", TRUE)),
                use.names= TRUE,
                idcol="file")  %>%
  as_tibble()  %>%
  dplyr::filter(Final_Score >= 0.5) %>%
  mutate(
    file = file %>% str_remove("Finaloutput.Variant.Ranked.OpenCRAVAT.ACCC.CRC.SNVs.") %>% 
    str_remove("Finaloutput.Variant.Ranked.OpenCRAVAT.ACCC.CRC.InDels."))  %>%
  dplyr::select(chrom, pos, ref_base, alt_base, samples, hugo, so, Final_Score) %>% 
  dplyr::filter(so %in% c("frameshift_elongation","frameshift_truncation","inframe_deletion", 
  "missense_variant", "splice_site_variant","start_lost","stop_gained","stop_lost",
  "inframe_insertion"))

DT_2 <- DT %>% separate_rows(samples, sep=";")

group_size_MSS <- length(unique(DT_2$samples)) # check
group_size_MSS

mss_cptac_dat <- DT_2 %>% select(hugo, samples) %>% distinct(hugo,samples) %>% 
group_by(hugo) %>% dplyr::count(hugo,sort=T, name="MSS_sample_Freq")

freq_genes_MSS_CPTAC_dat <- mss_cptac_dat %>% filter(MSS_sample_Freq>=4)
freq_genes_MSS_CPTAC <- mss_cptac_dat %>% filter(MSS_sample_Freq>=4) %>% pull(hugo)

################## step of the customized process ##################################

m1Name = "CPTAC_MSI"
m2Name = "CPTAC_MSS"

sampleSummary = data.table::data.table(Cohort = c(m1Name, 
m2Name), SampleSize = c(group_size_MSI, group_size_MSS))

uniqueGenes = unique(c(freq_genes_MSI_CPTAC, freq_genes_MSS_CPTAC)) # or union

MSI.CPTAC.comGenes = intersect(freq_genes_MSI_CPTAC,uniqueGenes)
MSS.CPTAC.comGenes = intersect(freq_genes_MSS_CPTAC, uniqueGenes)

cptac.msi.mss.merged = full_join(freq_genes_MSI_CPTAC_dat,freq_genes_MSS_CPTAC_dat)

cptac.msi.mss.merged[is.na(cptac.msi.mss.merged)] = 0

# modified version of the function mafCompare from maftools R package
# credits also here: https://rdrr.io/bioc/maftools/man/mafCompare.html

fisherTable = lapply(seq_len(nrow(cptac.msi.mss.merged)), function(i) {
  gene = cptac.msi.mss.merged[i, 1]$hugo
  m1Mut = cptac.msi.mss.merged[i, 2]$MSI_sample_Freq
  m2Mut = cptac.msi.mss.merged[i, 3]$MSS_sample_Freq
  ft_mat = matrix(c(m1Mut, group_size_MSI - m1Mut, m2Mut, 
                    group_size_MSS - m2Mut), byrow = TRUE, nrow = 2)
  if (length(which(x = ft_mat == 0)) > 0) {
    ft_mat = ft_mat + 1
  }
  xf = fisher.test(ft_mat, conf.int = TRUE, conf.level = 0.95)
  pval = xf$p.value
  or = xf$estimate
  ci.up = xf$conf.int[2]
  ci.low = xf$conf.int[1]
  tdat = data.table::data.table(Hugo_Symbol = gene, m1Mut, 
                                m2Mut, pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
  tdat
})

fisherTable = data.table::rbindlist(l = fisherTable, use.names = TRUE, 
                                    fill = TRUE)
fisherTable = fisherTable[order(pval)]
fisherTable[, `:=`(adjPval, p.adjust(p = pval, method = "fdr"))]
colnames(fisherTable)[2:3] = c(m1Name, m2Name)

CPTAC_Fisher_out_hg38_scored <- fisherTable %>% filter(adjPval <= 0.1) %>% 
mutate(Freq_Ratio_MSI_CPTAC=CPTAC_MSI/14) %>% mutate(Freq_Ratio_MSS_CPTAC=CPTAC_MSS/65)   

# write_tsv(CPTAC_Fisher_out_hg38_scored, file="CPTAC.FisherAdaptation.MSIvsMSS.hg38.ScoredGenes.Direct.79samples.tsv")

##########################################################################################
############ continue with one additional filter for mutational frequency ################

hg38.out.sel <- CPTAC_Fisher_out_hg38_scored %>% filter(Freq_Ratio_MSI_CPTAC > Freq_Ratio_MSS_CPTAC)

# write_tsv(hg38.out.sel, file="CPTAC.Colon.hg38.Diff.Mut.Freq.Scored.AdjPVal.FreqRatio.tsv")

hg38.out.sel.genes <- hg38.out.sel %>% pull("Hugo_Symbol")

# continue with relative functional analysis of the selected genes:

# load the saved object of expressed genes in the CPTAC cohort
exp.genes.cptac <- read_tsv(here("Data","ExpressedGenes.CPTAC.RNAseq.ACCC.Paper.MSI.tsv"))
exp.genes.input.cptac <- as.character(exp.genes.cptac$ExpressedGenes.CPTAC.CRC.SelSamples.RNAseq)

background_genes.id <- bitr(exp.genes.input.cptac, fromType="SYMBOL", toType="ENTREZID", 
OrgDb="org.Hs.eg.db")

background_genes_vector <- background_genes.id$ENTREZID

eg.cptac.freq <- bitr(hg38.out.sel.genes, fromType="SYMBOL", toType="ENTREZID", 
OrgDb="org.Hs.eg.db")

CPTAC_enrich_freq <- enrichPathway(gene=eg.cptac.freq$ENTREZID, pvalueCutoff = 0.05, 
                                   readable=TRUE, universe = background_genes_vector,
                                   minGSSize = 10, qvalueCutoff = 0.1)


CPTAC_enrich_freq_2 <- mutate(CPTAC_enrich_freq, richFactor = Count / as.numeric
                              (sub("/\\d+", "", BgRatio)))

CPTAC_enrich_dt <- CPTAC_enrich_freq_2[]

# write_tsv(CPTAC_enrich_dt, file="CPTAC.MSIvsMSS.Signif.Freq.Scored.Reactome.clusterProfiler.ACCCPaper.BackgroundSet.tsv")

############################################################################################
########## *load and re-create the same plot with only selected terms* #####################
############################################################################################

common_paths <- read_tsv(here("Data","Common_Mechanisms_TCGA_CPTAC.txt"))
common_paths_out <- common_paths %>% pull("Common_Terms")

cptac_mod_plot <- CPTAC_enrich_dt %>% dplyr::select(ID,Description,p.adjust,Count,richFactor) %>% 
filter(ID%in%common_paths_out)

ggplot(cptac_mod_plot,
       aes(richFactor, Description, p.adjust)) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2",
                                  "#7e62a3"),
                        trans = "log10",
                        guide=guide_colorbar(reverse=TRUE,
                                             order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_dose(12) +
  xlab("Rich Factor") +
  ylab(NULL) + 
  ggtitle("CPTAC MSI_vs_MSS")

#############################################################################################
