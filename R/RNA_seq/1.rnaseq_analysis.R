######################################################
# 1.rnaseq_analysis_larger_cohort_v2.R
# created on Sept 15, 2021
# lucile.neyton@ucsf.edu

# This script aims at analysing RNA-Seq data for dex vs on-dex patients
######################################################

rm(list = ls())
setwd("/Users/lucileneyton/Box Sync/COMET/data_analysis/")

# load libraries
library(magmaR)
library(DESeq2)
library(stringr)
library(ggrepel)
library(msigdbr)
library(fgsea)
library(mixOmics)
library(DescTools)
library(ggforestplot)

# data files
alt_data_path <- "/Users/lucileneyton/Box Sync/COMET/pratik_lca/data/"
data_path <- "/Users/lucileneyton/Box Sync/COMET/data/larger_cohort/"
results_path <- "/Users/lucileneyton/Box Sync/COMET/results/larger_cohort/"

# reset plot settings
theme_set(theme_grey())

#########################
# DATA LOADING
#########################
# metadata
metadata_df <- readRDS(paste0(data_path, "processed/df_metadata.RDS"))

# annotation data
annot_data <- read.table(paste0(data_path, "../raw/Homo_sapiens.GRCh38.103.gtf"), sep="\t")
annot_data_gene <- annot_data[annot_data$V3 == "gene", ]
annot_data_gene$gene_id <- sapply(annot_data_gene$V9, function(x) strsplit(strsplit(x, ";")[[1]][1], " ")[[1]][2])
annot_data_gene$biotype <- sapply(annot_data_gene$V9, function(x) strsplit(strsplit(x, ";")[[1]][5], " ")[[1]][3])

# keep only protein-coding genes
annot_data_pcgene <- annot_data_gene[annot_data_gene$biotype == "protein_coding", ]

# collection times and dex ttt info
dex_tp_data <- readRDS(paste0(data_path, "processed/dex_tp_data_unfiltered_with_impacc_status.RDS"))

# cnt data
cnt_files <- c(paste0(paste0(data_path,"../processed/"), "2022-03-07.D0pbmc.salmon.merged.gene_counts.tsv"),
               paste0(paste0(data_path,"../processed/"), "2022-06-13.D0pbmc.salmon.merged.gene_counts.tsv"),
               paste0(paste0(data_path,"../processed/"), "2022-06-21.D0_4_eta.salmon.merged.gene_counts.tsv"),
               paste0(paste0(data_path,"../processed/"), "2022-02-07.D0eta.salmon.merged.gene_counts.tsv"))
cnt_data_list <- lapply(cnt_files, function(x) read.table(x, sep="\t", header = T, row.names = 1))
cnt_data <- do.call(cbind, cnt_data_list)

# redcap data to get collection times
# read data from the latest pull
date_ <- readLines(paste(data_path, "raw/last_pull.txt", sep=""))
files_list <- list.files(paste0(data_path, "raw/"))
redcap_data_file <- files_list[grepl("COMETPatientAndSampl_DATA_", files_list) & grepl(date_, files_list)]
redcap_data <- read.csv(paste(data_path, paste("raw/", redcap_data_file, sep=""), sep=""), colClasses=c("record_id"="character"))

redcap_data_D0 <- redcap_data[redcap_data$redcap_event_name=="d0_arm_1", ]
rownames(redcap_data_D0) <- redcap_data_D0$record_id

# for dex and non-dex patients
dex_tp_data_all <- readRDS(paste0(data_path, "processed/dex_tp_data_unfiltered_with_impacc_status.RDS"))

# intubation data
vfd_data_filt <- readRDS(paste0(data_path, "processed/vfd_data_filt.RDS"))

# cytokine data
cytok_mofa <- readRDS(paste0(data_path, "processed/cytok_data.RDS"))

#########################
# DATA PRE-PROCESSING
#########################
# protein-coding genes only
cnt_data <- cnt_data[rownames(cnt_data) %in% annot_data_pcgene$gene_id, ]

# save gene mapping
gene_names <- cnt_data[, "gene_name", drop = F]

# filtered data
dex_tp_data$overlap[(is.na(dex_tp_data$overlap)) & (dex_tp_data$dex)] <- FALSE
dex_tp_data_overlap <- dex_tp_data[(!(dex_tp_data$dex)) | (dex_tp_data$overlap), ]

dex_tp_data_overlap$dex <- factor(dex_tp_data_overlap$dex)
dex_tp_data_overlap$dex <- relevel(dex_tp_data_overlap$dex, ref = "FALSE")

# EA sc
dex_tp_data_all_easc <- dex_tp_data_overlap[dex_tp_data_overlap$data_type=="EA_scRNASeq", ]
dex_tp_data_all_easc$sample <- dex_tp_data_all_easc$timepoint

dex_tp_data_all_easc$intub <- rep(TRUE, nrow(dex_tp_data_all_easc))

# earliest time point (max 4)
dex_tp_data_all_easc <- dex_tp_data_all_easc[grepl(dex_tp_data_all_easc$timepoint, pattern = "*D0$|*D1$|*D2$|*D3$|*D4$"), ]
dex_tp_data_all_easc <- dex_tp_data_all_easc[order(dex_tp_data_all_easc$timepoint), ]
dex_tp_data_all_easc <- dex_tp_data_all_easc[!duplicated(dex_tp_data_all_easc$patient), ]
saveRDS(dex_tp_data_all_easc[, c("sample", "dex", "intub")], paste0(results_path, "sceta_data_for_dems.RDS"))

# time difference bw dex first dose and sample
dex_tp_data_overlap$dt <- as.Date(dex_tp_data_overlap$collection_date) - dex_tp_data_overlap$earliest_dex

# differences using collection and treatment times
dex_tp_data_overlap$collection_date_char <- as.POSIXlt(dex_tp_data_overlap$collection_date_char, format="%Y-%m-%d %H:%M")
dex_tp_data_overlap$earliest_dex_char <- as.POSIXlt(dex_tp_data_overlap$earliest_dex_char, format="%Y-%m-%d %H:%M")
dex_tp_data_overlap$dt_char <- dex_tp_data_overlap$collection_date_char - dex_tp_data_overlap$earliest_dex_char

# filter
dex_tp_data_overlap_filt <- dex_tp_data_overlap[dex_tp_data_overlap$data_type %in% c("EA_RNASeq", "PBMC_RNASeq"), ]

# format
dex_tp_data_overlap_filt$patient_short <- dex_tp_data_overlap_filt$patient
dex_tp_data_overlap_filt$patient_short <- sapply(dex_tp_data_overlap_filt$patient_short, function(x) strsplit(x, "-")[[1]][2])

dex_tp_data_overlap_filt$timepoint_short <- dex_tp_data_overlap_filt$timepoint
dex_tp_data_overlap_filt$timepoint_short <- sapply(dex_tp_data_overlap_filt$timepoint_short, function(x) strsplit(x, "-")[[1]][3])

dex_tp_data_overlap_filt$data_type_short <- as.character(dex_tp_data_overlap_filt$data_type)
dex_tp_data_overlap_filt$data_type_short[dex_tp_data_overlap_filt$data_type_short == "EA_RNASeq"] <- "ETA"
dex_tp_data_overlap_filt$data_type_short[dex_tp_data_overlap_filt$data_type_short == "PBMC_RNASeq"] <- "PBMC"

dex_tp_data_overlap_filt$rnaseq_samp <- paste(dex_tp_data_overlap_filt$patient_short, dex_tp_data_overlap_filt$timepoint_short, sep=".")
dex_tp_data_overlap_filt$rnaseq_samp <- paste0(dex_tp_data_overlap_filt$rnaseq_samp, dex_tp_data_overlap_filt$data_type_short)

rsq_vect <- sapply(colnames(cnt_data), function(x) strsplit(x, "\\.")[[1]][3])
cnt_data_complete_sample_names <- colnames(cnt_data)
colnames(cnt_data) <- sapply(colnames(cnt_data), function(x) strsplit(x, "1.R")[[1]][1])
names(rsq_vect) <- colnames(cnt_data)

dex_tp_data_overlap_filt_seq <- dex_tp_data_overlap_filt[dex_tp_data_overlap_filt$rnaseq_samp %in% colnames(cnt_data), ]

metadata_df$patient <- rownames(metadata_df)
dex_tp_data_overlap_filt_seq <- merge(dex_tp_data_overlap_filt_seq, metadata_df[, !(colnames(metadata_df) %in% c("timepoint", "dex"))], all.x=T, by=c("patient"))

# day 0 only
cnt_data_complete_sample_names <- cnt_data_complete_sample_names[grepl("D0", colnames(cnt_data))]
d0_cnt_data <- cnt_data[, grepl("D0", colnames(cnt_data))]

d0_4_cnt_data <- cnt_data[, grepl("*D0ETA|*D1ETA|*D2ETA|*D3ETA|*D4ETA", colnames(cnt_data))]

d0_eta_cnt_data <- d0_4_cnt_data[, grepl("ETA", colnames(d0_4_cnt_data))]

cnt_data_complete_sample_names <- cnt_data_complete_sample_names[grepl("PBMC", colnames(d0_cnt_data))]
d0_pbmc_cnt_data <- d0_cnt_data[, grepl("PBMC", colnames(d0_cnt_data))]

eta_uniq_patients <- unique(sapply(colnames(d0_eta_cnt_data), function(x) strsplit(x, "\\.")[[1]][1]))
eta_to_keep <- c()
for (eta_uniq_patient in eta_uniq_patients){
  eta_uniq_patient_cols <- colnames(d0_eta_cnt_data)[grepl(paste0(eta_uniq_patient, "\\."), colnames(d0_eta_cnt_data))]
  eta_to_keep <- c(eta_to_keep, sort(eta_uniq_patient_cols)[1])
}

d0_eta_cnt_data <- d0_eta_cnt_data[, eta_to_keep]

d0_eta_dex_tp_data <- dex_tp_data_overlap_filt_seq[dex_tp_data_overlap_filt_seq$rnaseq_samp %in% colnames(d0_eta_cnt_data), ]
d0_pbmc_dex_tp_data <- dex_tp_data_overlap_filt_seq[dex_tp_data_overlap_filt_seq$rnaseq_samp %in% colnames(d0_pbmc_cnt_data), ]

# only keep pre-selected samples
d0_eta_cnt_data <- d0_eta_cnt_data[, colnames(d0_eta_cnt_data) %in% d0_eta_dex_tp_data$rnaseq_samp]

cnt_data_complete_sample_names <- cnt_data_complete_sample_names[colnames(d0_pbmc_cnt_data) %in% d0_pbmc_dex_tp_data$rnaseq_samp]
d0_pbmc_cnt_data <- d0_pbmc_cnt_data[, colnames(d0_pbmc_cnt_data) %in% d0_pbmc_dex_tp_data$rnaseq_samp]

rownames(d0_eta_dex_tp_data) <- d0_eta_dex_tp_data$rnaseq_samp
d0_eta_dex_tp_data <- d0_eta_dex_tp_data[colnames(d0_eta_cnt_data), ]
rownames(d0_pbmc_dex_tp_data) <- d0_pbmc_dex_tp_data$rnaseq_samp
d0_pbmc_dex_tp_data <- d0_pbmc_dex_tp_data[colnames(d0_pbmc_cnt_data), ]

d0_eta_dex_tp_data$rsq <- rsq_vect[rownames(d0_eta_dex_tp_data)]
d0_pbmc_dex_tp_data$rsq <- rsq_vect[rownames(d0_pbmc_dex_tp_data)]

d0_eta_dex_tp_data$collection_date <- redcap_data_D0[as.character(d0_eta_dex_tp_data$record_id), "dos_esc_ta"]
d0_pbmc_dex_tp_data$collection_date <- redcap_data_D0[as.character(d0_pbmc_dex_tp_data$record_id), "process_toc_4"]

# all ETA ventilated (TA collected from vent device)
d0_eta_dex_tp_data$vent_started_before_D0[!d0_eta_dex_tp_data$vent_status] <- FALSE
d0_eta_dex_tp_data$vent_started_before_D0[(d0_eta_dex_tp_data$vent_status) & (d0_eta_dex_tp_data$first_intub < d0_eta_dex_tp_data$collection_date)] <- TRUE
d0_eta_dex_tp_data$vent_started_before_D0[(d0_eta_dex_tp_data$vent_status) & (d0_eta_dex_tp_data$first_intub >= d0_eta_dex_tp_data$collection_date)] <- FALSE

d0_pbmc_dex_tp_data$vent_started_before_D0[!d0_pbmc_dex_tp_data$vent_status] <- FALSE
d0_pbmc_dex_tp_data$vent_started_before_D0[(d0_pbmc_dex_tp_data$vent_status) & (d0_pbmc_dex_tp_data$first_intub < d0_pbmc_dex_tp_data$collection_date)] <- TRUE
d0_pbmc_dex_tp_data$vent_started_before_D0[(d0_pbmc_dex_tp_data$vent_status) & (d0_pbmc_dex_tp_data$first_intub >= d0_pbmc_dex_tp_data$collection_date)] <- FALSE

# check intubation status
# agreement between cytokine and rna seq (same tube used)
dex_tp_data_all <- dex_tp_data_all[dex_tp_data_all$data_type %in% c("PBMC_RNASeq", "Cytokine"), ]

intub_vect <- c()
for (patient_ in rownames(metadata_df)){
  patient_longid <- as.numeric(strsplit(patient_, "-HS")[[1]][2]) + 1000
  
  if (paste0(patient_, "-D0") %in% dex_tp_data_all$timepoint){
    patient_col <- unique(dex_tp_data_all[dex_tp_data_all$patient==patient_ & grepl("-D0", dex_tp_data_all$timepoint), "collection_date"])
    patient_col <- as.Date(patient_col)
    patient_vfd <- vfd_data_filt[as.character(patient_longid), ]
    
    patient_intubs <- list()
    if (!is.na(patient_vfd$vent_start_date)){
      patient_intubs[["1"]] <- c(patient_vfd$vent_start_date, patient_vfd$vent_end_date)
      if (!is.na(patient_vfd$reintubation1)){
        patient_intubs[["2"]] <- c(patient_vfd$reintubation1, patient_vfd$reextubation1)
        if (!is.na(patient_vfd$reintubation2)){
          patient_intubs[["3"]] <- c(patient_vfd$reintubation2, patient_vfd$reextubation2)
        }
      }
    }
    
    intub_status <- F
    for (patient_intub in patient_intubs){
      intub_start <- patient_intub[1]
      intub_end <- patient_intub[2]
      
      if (patient_col>intub_start & patient_col<=intub_end){
        intub_status <- T
      }
    }
  }else{
    intub_status <- NA
  }
  
  intub_vect <- c(intub_vect, intub_status)
}

metadata_df$intub <- intub_vect

# check WB
#vfd_data <- read.csv(paste0(data_path, "raw/COMET_VFD_DATA.csv")) # REDCap export
vfd_data <- read.csv("/Users/lucileneyton/OneDrive - University of California, San Francisco/UCSF/COMET/data/raw/COMET_VFD_DATA.csv") 

dex_tp_data_wb <- readRDS(paste0(data_path, "processed/dex_tp_data_unfiltered_with_impacc_status.RDS"))
dex_tp_data_wb <- dex_tp_data_wb[dex_tp_data_wb$data_type %in% c("WB_scRNASeq"), ]

dex_tp_data_wb$overlap[(is.na(dex_tp_data_wb$overlap)) & (dex_tp_data_wb$dex)] <- FALSE
dex_tp_data_wb <- dex_tp_data_wb[(!(dex_tp_data_wb$dex)) | (dex_tp_data_wb$overlap), ]

intub_vect_wb <- c()
for (patient_ in dex_tp_data_wb$patient){
  patient_longid <- as.numeric(strsplit(patient_, "-HS")[[1]][2]) + 1000
  
  if (paste0(patient_, "-D0") %in% dex_tp_data_wb$timepoint){
    patient_col <- unique(dex_tp_data_wb[dex_tp_data_wb$patient==patient_ & grepl("-D0", dex_tp_data_wb$timepoint), "collection_date"])
    patient_col <- as.Date(patient_col)
    patient_vfd <- vfd_data[vfd_data$record_id==as.character(patient_longid), ]
    
    patient_intubs <- list()
    if (!is.na(patient_vfd$vent_start_date)){
      patient_intubs[["1"]] <- c(patient_vfd$vent_start_date, patient_vfd$vent_end_date)
      if (!is.na(patient_vfd$reintubation1)){
        patient_intubs[["2"]] <- c(patient_vfd$reintubation1, patient_vfd$reextubation1)
        if (!is.na(patient_vfd$reintubation2)){
          patient_intubs[["3"]] <- c(patient_vfd$reintubation2, patient_vfd$reextubation2)
        }
      }
    }
    
    intub_status <- F
    for (patient_intub in patient_intubs){
      intub_start <- patient_intub[1]
      intub_end <- patient_intub[2]
      
      if (intub_start!=""){
        if (patient_col>intub_start & patient_col<=intub_end){
          intub_status <- T
        }
      }
      
    }
  }else{
    intub_status <- NA
  }
  
  intub_vect_wb <- c(intub_vect_wb, intub_status)
}

dex_tp_data_wb$intub_wb <- intub_vect_wb
dex_tp_data_wb$intub <- intub_vect_wb

dex_tp_data_wb$sample <- dex_tp_data_wb$timepoint
dex_tp_data_wb <- dex_tp_data_wb[dex_tp_data_wb$intub, ]
saveRDS(dex_tp_data_wb[, c("sample", "dex", "intub")], paste0(results_path, "scwb_data_for_dems.RDS"))

# save list of included patients
wb_vent_patients <- dex_tp_data_wb[, c("patient", "dex"), drop=F]
write.csv(wb_vent_patients, row.names = F,
          paste0(results_path, "wb_vent_patients.csv"))

#########################
# DGEA
#########################
# PBMC
# remove outlier(s) -> stranded libraries
outliers_ <- c("HS51.D0PBMC", "HS53.D0PBMC", "HS130.D0PBMC", "HS149.D0PBMC", 
               "HS29.D0PBMC", "HS33.D0PBMC", "HS37.D0PBMC", "HS3.D0PBMC",
               "HS41.D0PBMC")

cnt_data_complete_sample_names <- cnt_data_complete_sample_names[!(colnames(d0_pbmc_cnt_data) %in% outliers_)]
d0_pbmc_cnt_data <- d0_pbmc_cnt_data[, !(colnames(d0_pbmc_cnt_data) %in% outliers_)]
d0_pbmc_dex_tp_data <- d0_pbmc_dex_tp_data[!(d0_pbmc_dex_tp_data$rnaseq_samp %in% outliers_), ]

# build DESeq2 object
d0_pbmc_dex_tp_data$age_scaled <- scale(d0_pbmc_dex_tp_data$age)
d0_pbmc_dex_tp_data$sex_at_birth <- as.factor(d0_pbmc_dex_tp_data$sex_at_birth)

d0_pbmc_dex_tp_data$intub <- metadata_df[d0_pbmc_dex_tp_data$patient, "intub"]
d0_pbmc_dex_tp_data$intub <- as.factor(d0_pbmc_dex_tp_data$intub)

d0_pbmc_dex_tp_data$sample <- d0_pbmc_dex_tp_data$timepoint

# use only the 500 most variable genes
d0_pbmc_dex_tp_data$intub <- factor(d0_pbmc_dex_tp_data$intub, levels=c(T, F))

# intub
saveRDS(d0_pbmc_dex_tp_data[d0_pbmc_dex_tp_data$intub==TRUE, c("sample", "dex", "intub")], paste0(results_path, "bulkpbmc_data_for_dems.RDS"))
pbmc_intub <- rownames(d0_pbmc_dex_tp_data[d0_pbmc_dex_tp_data$intub==TRUE, ])

cnt_data_complete_sample_names <- cnt_data_complete_sample_names[colnames(d0_pbmc_cnt_data) %in% pbmc_intub]
d0_pbmc_cnt_data_intub <- d0_pbmc_cnt_data[, colnames(d0_pbmc_cnt_data) %in% pbmc_intub]
d0_pbmc_dex_tp_data_intub <- d0_pbmc_dex_tp_data[d0_pbmc_dex_tp_data$rnaseq_samp %in% pbmc_intub, ]

write.csv(cnt_data_complete_sample_names, paste0(data_path, "processed/cnt_data_complete_sample_names.csv"))

# build DESeq2 object
d0_pbmc_dex_tp_data_intub$age_scaled <- scale(d0_pbmc_dex_tp_data_intub$age)
d0_pbmc_dex_tp_data_intub$sex_at_birth <- as.factor(d0_pbmc_dex_tp_data_intub$sex_at_birth)

dds_pbmc_intub <- DESeqDataSetFromMatrix(countData = round(d0_pbmc_cnt_data_intub), 
                                         colData = d0_pbmc_dex_tp_data_intub, design = ~ dex + age_scaled + sex_at_birth)

# save list of included patients
pbmc_vent_patients <- d0_pbmc_dex_tp_data_intub[, "dex", drop=F]
rownames(pbmc_vent_patients) <- paste0("MVIR1-", sapply(rownames(pbmc_vent_patients), function(x) str_split(x, ".D0PBMC")[[1]][1]))
write.csv(pbmc_vent_patients,
          paste0(results_path, "pbmc_vent_patients.csv"))

# run DESeq
dds_pbmc_intub <- DESeq(dds_pbmc_intub)

# extract results
res_pbmc_intub <- lfcShrink(dds_pbmc_intub, coef=2, type="apeglm")

# sort the genes from lowest to highest given adjusted p-values
res_pbmc_intub <- res_pbmc_intub[order(res_pbmc_intub$padj, decreasing = F), ]

# replace NA values with 1s
res_pbmc_intub$padj[is.na(res_pbmc_intub$padj)] <- 1

res_pbmc_intub$gene_name <- gene_names[rownames(res_pbmc_intub), ]
dim(res_pbmc_intub[res_pbmc_intub$padj < 0.1, ])

res_pbmc_intub$sig <- res_pbmc_intub$padj < 0.1
res_pbmc_intub <- as.data.frame(res_pbmc_intub)

write.csv(res_pbmc_intub, paste0(results_path, "pbmc_intub_dgea.csv"))

label_data <- res_pbmc_intub[1:min(c(50, nrow(res_pbmc_intub[res_pbmc_intub$sig, ]))), ]

pdf(paste0(results_path, "pbmc_intub_volcano_plot.pdf"), width = 7, height = 7)
p <- ggplot(res_pbmc_intub, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(col = sig)) +
  scale_color_manual(values = c("black", "red")) +
  theme(legend.position = "none") +
  geom_text_repel(data = label_data,
                  aes(label = gene_name),
                  max.overlaps=20, size=5) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value")  +
  theme_bw() + theme(text = element_text(size = 18))
print(p)
dev.off()

# low log2fc genes should be black
res_pbmc_intub_filt <- res_pbmc_intub
res_pbmc_intub_filt[res_pbmc_intub_filt$sig & abs(res_pbmc_intub_filt$log2FoldChange) < 0.1, "sig"] <- FALSE
label_data_filt <- res_pbmc_intub_filt[res_pbmc_intub_filt$sig, ][1:50, ]

pdf(paste0(results_path, "pbmc_intub_filt_volcano_plot.pdf"), width = 7, height = 7)
p <- ggplot(res_pbmc_intub_filt, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(col = sig)) +
  scale_color_manual(values = c("black", "red")) +
  theme(legend.position = "none") +
  geom_text_repel(data = label_data_filt,
                  aes(label = gene_name),
                  max.overlaps=20, size=5) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value")  +
  theme_bw() + theme(text = element_text(size = 18))
print(p)
dev.off()

# transform data for PCA plot
vsd_pbmc_intub <- vst(dds_pbmc_intub, blind = FALSE)

#########################
# GSEA
#########################
# load gene set DB
m_dbs <- list(reactome=msigdbr(species="Homo sapiens", "C2", "CP:REACTOME"))

res_pbmc_intub <- lfcShrink(dds_pbmc_intub, coef=2, type="apeglm")
res_pbmc_intub <- res_pbmc_intub[order(res_pbmc_intub$padj, decreasing = F), ]

# replace NA values with 1s
res_pbmc_intub$padj[is.na(res_pbmc_intub$padj)] <- 1

res_pbmc_intub$gene_name <- gene_names[rownames(res_pbmc_intub), ]

res_pbmc_intub$sig <- res_pbmc_intub$padj < 0.1
res_pbmc_intub <- as.data.frame(res_pbmc_intub)

for (m_db_id in names(m_dbs)){
  m_db <- m_dbs[[m_db_id]]    
  
  m_db_list <- lapply(unique(m_db$gs_name), function(x) unname(unlist(m_db[m_db$gs_name == x, "gene_symbol"])))
  names(m_db_list) <- unique(m_db$gs_name)
  
  # sort genes given log fc values
  res_pbmc_intub_sorted <- res_pbmc_intub[order(res_pbmc_intub$log2FoldChange, decreasing = T), ]
  gene_ranks_pbmc_intub <- res_pbmc_intub_sorted$log2FoldChange
  names(gene_ranks_pbmc_intub) <- res_pbmc_intub_sorted$gene_name
  gene_ranks_pbmc_intub <- gene_ranks_pbmc_intub[!is.na(gene_ranks_pbmc_intub)]
  
  # GSEA
  set.seed(123)
  gsea_res_pbmc_intub <- fgsea(m_db_list, stats=gene_ranks_pbmc_intub, nPermSimple=10000)
  
  # filter
  gsea_res_pbmc_intub$padj[is.na(gsea_res_pbmc_intub$padj)] <- 1
  gsea_res_pbmc_intub_df <- data.frame(gsea_res_pbmc_intub[gsea_res_pbmc_intub$padj < 0.1, ])
  gsea_res_pbmc_intub_df$leadingEdge <- as.character(gsea_res_pbmc_intub_df$leadingEdge)
  gsea_res_pbmc_intub_df$pathway <- as.factor(gsea_res_pbmc_intub_df$pathway)
  
  write.csv(gsea_res_pbmc_intub_df, paste0(results_path, paste0(m_db_id, "_gsea_pbmc_intub.csv")))
  
  if (dim(gsea_res_pbmc_intub_df)[1] > 0){
    # dot plot
    gsea_res_pbmc_df_pos <- gsea_res_pbmc_intub_df[gsea_res_pbmc_intub_df$NES > 0, ]
    gsea_res_pbmc_df_neg <- gsea_res_pbmc_intub_df[gsea_res_pbmc_intub_df$NES < 0, ]
    
    gsea_res_pbmc_df_pos <- gsea_res_pbmc_df_pos[order(gsea_res_pbmc_df_pos$padj, decreasing = F), ]
    gsea_res_pbmc_df_neg <- gsea_res_pbmc_df_neg[order(gsea_res_pbmc_df_neg$padj, decreasing = F), ]
    
    gsea_res_pbmc_df_pos_top10 <- as.character(gsea_res_pbmc_df_pos[1:(min(nrow(gsea_res_pbmc_df_pos), 10)), "pathway"])
    gsea_res_pbmc_df_neg_top10 <- as.character(gsea_res_pbmc_df_neg[1:(min(nrow(gsea_res_pbmc_df_neg), 10)), "pathway"])
    
    gsea_res_pbmc_intub_df <- gsea_res_pbmc_intub_df[order(gsea_res_pbmc_intub_df$NES, decreasing=T), ]
    gsea_res_pbmc_intub_df_filt <- gsea_res_pbmc_intub_df[gsea_res_pbmc_intub_df$pathway %in% c(gsea_res_pbmc_df_pos_top10, gsea_res_pbmc_df_neg_top10), ]
    
    if (max(sapply(as.character(gsea_res_pbmc_intub_df_filt$pathway), nchar)) > 100){
      width_ <- 14
    }else{
      if (max(sapply(as.character(gsea_res_pbmc_intub_df_filt$pathway), nchar)) > 75){
        width_ <- 12
      }else{
        width_ <- 10
      }
    }
    
    pdf(paste0(results_path, paste0(m_db_id, "_dotplot_gsea_pbmc_intub.pdf")), width = width_*1.2, height = 8)
    print(ggplot(gsea_res_pbmc_intub_df_filt, aes(x= NES, y=pathway)) +
            geom_point(aes(colour=NES), size=3) +
            scale_colour_gradient2(
              low = "#388ecc",
              mid = "white",
              high = "#f68b33",
              midpoint = 0,
              aesthetics = "colour"
            ) + 
            theme_bw() + theme(text = element_text(size = 18)) +
            xlab("NES") +
            ylab(paste0(m_db_id, " gene set")) +
            scale_y_discrete(limits = rev(gsea_res_pbmc_intub_df_filt$pathway)))
    dev.off()
  }
}

#########################
# DIABLO - WITHOUT ETA
#########################
metadata_df_pbmc_cytok <- metadata_df[!is.na(metadata_df$intub), ]
vent_status <- T

# add cytokine data
cytok_mofa <- readRDS(paste0(data_path, "processed/cytok_data.RDS"))

vsd_pbmc_intub_df <- assay(vsd_pbmc_intub)

colnames(vsd_pbmc_intub_df) <- str_replace(colnames(vsd_pbmc_intub_df), "PBMC", "")
colnames(vsd_pbmc_intub_df) <- str_replace(colnames(vsd_pbmc_intub_df), "\\.", "-")

vsd_pbmc_df_filt <- vsd_pbmc_intub_df

# filter by intubation status
intub_samples <- rownames(metadata_df_pbmc_cytok[metadata_df_pbmc_cytok$intub==vent_status, ])
cytok_mofa <- cytok_mofa[rownames(cytok_mofa) %in% intub_samples, ]

rownames(cytok_mofa) <- paste0(str_replace(rownames(cytok_mofa), "MVIR1-", ""), "-D0")

# samples in common only
vsd_pbmc_df_filt <- vsd_pbmc_df_filt[, intersect(colnames(vsd_pbmc_df_filt), rownames(cytok_mofa))]
cytok_mofa <- cytok_mofa[intersect(colnames(vsd_pbmc_df_filt), rownames(cytok_mofa)), ]

# scale data
cytok_mofa <- scale(cytok_mofa)
cytok_mofa <- t(cytok_mofa)

# scale and select features
vsd_pbmc_top500 <- t(vsd_pbmc_df_filt)[,colnames(t(vsd_pbmc_df_filt))[order(colVars(t(vsd_pbmc_df_filt)), decreasing = T)][1:2000]]

# replace ENSG with gene symbols
colnames(vsd_pbmc_top500) <- gene_names[colnames(vsd_pbmc_top500), "gene_name"]

# create list of data modalities
pbmc_mofa <- scale(vsd_pbmc_top500)
pbmc_mofa <- t(pbmc_mofa)

# samples in common
mofa_samples <- Reduce(union, list(colnames(pbmc_mofa), colnames(cytok_mofa)))

# add NA columns for samples with missing data type
data_list <- list(PBMC = as.matrix(pbmc_mofa),
                  cytok = as.matrix(cytok_mofa))

# add meta data
mofa_metadata <- metadata_df[paste0("MVIR1-", str_replace(mofa_samples, "-D0", "")), ]
mofa_metadata$sample <- str_replace(rownames(mofa_metadata), "MVIR1-", "")
mofa_metadata$sample <- paste0(mofa_metadata$sample, "-D0")

mofa_metadata$sample <- paste0(mofa_metadata$patient, "-D0")
mofa_metadata$intub[is.na(mofa_metadata$intub)] <- T # ETA only

# cytok
# filter by intubation status
cytok_mofa <- readRDS(paste0(data_path, "processed/cytok_data.RDS"))
cytok_mofa <- cytok_mofa[rownames(cytok_mofa) %in% intub_samples, ]

# more than 10% missing values
cyto_drop <- colnames(cytok_mofa)[colSums(is.na(cytok_mofa))/dim(cytok_mofa)[1]>0.1]
cytok_mofa <- cytok_mofa[, !(colnames(cytok_mofa) %in% cyto_drop)]

# drop rows with NA values
samp_drop <- rownames(cytok_mofa)[rowSums(is.na(cytok_mofa))>0]

# filtered set
cytok_mofa <- cytok_mofa[!(rownames(cytok_mofa) %in% samp_drop), ]

# scale data
rownames(cytok_mofa) <- paste0(str_replace(rownames(cytok_mofa), "MVIR1-", ""), "-D0")
cytok_mofa <- scale(cytok_mofa)
cytok_mofa <- t(cytok_mofa)

# PBMC
# bulk RNA
vsd_pbmc_df_filt <- vsd_pbmc_intub_df

# scale and select features
vsd_pbmc_top500 <- t(vsd_pbmc_df_filt)[,colnames(t(vsd_pbmc_df_filt))[order(colVars(t(vsd_pbmc_df_filt)), decreasing = T)][1:2000]]

# replace ENSG with gene symbols
colnames(vsd_pbmc_top500) <- gene_names[colnames(vsd_pbmc_top500), "gene_name"]

# create list of data modalities
pbmc_mofa <- scale(vsd_pbmc_top500)
pbmc_mofa <- t(pbmc_mofa)

# samples in common
# 20 and not 21 because one sample had missing cytokine values
diablo_samples <- Reduce(intersect, list(colnames(pbmc_mofa), colnames(cytok_mofa)))

data_list <- list(PBMC = t(pbmc_mofa[, diablo_samples]),
                  cytok = t(cytok_mofa[, diablo_samples]))

diablo_samples_md <- sapply(diablo_samples, function(x) paste0("MVIR1-", str_replace(x, "-D0", "")))
dex_vect <- mofa_metadata[diablo_samples_md, "dex"]
names(dex_vect) <- diablo_samples

saveRDS(data.frame(sample=paste0("MVIR1-", diablo_samples), dex=dex_vect, intub=mofa_metadata[diablo_samples_md, "intub"]), paste0(results_path, "diablonoeta_data_for_dems.RDS"))

# balance model discriminative ability
val_ <- 0.5
design <- matrix(val_, ncol = length(data_list), nrow = length(data_list), 
                 dimnames = list(names(data_list), names(data_list)))
diag(design) <- 0
res_diablo <- block.splsda(X= data_list, Y = dex_vect, ncomp=2, scale=F, design = design)

dex_tp_data_overlap_with_dt_char <- readRDS(paste0(results_path, "dex_tp_data_overlap_with_dt_char.RDS"))
rownames(dex_tp_data_overlap_with_dt_char) <- dex_tp_data_overlap_with_dt_char$timepoint

# sorted by weight on component 1 (cytokine)
sorted_ids_cytok1 <- sort(res_diablo$variates$cytok[, "comp1"])
sorted_ids_cytok1_times <- dex_tp_data_overlap_with_dt_char[paste0("MVIR1-", names(sorted_ids_cytok1)), c("timepoint", "dt", "dt_char", "dex")]
sorted_ids_cytok1_times <- cbind(sorted_ids_cytok1_times, comp1=sort(res_diablo$variates$cytok[, "comp1"]))
sorted_ids_cytok1_times_filt <- sorted_ids_cytok1_times[sorted_ids_cytok1_times$dex, ]

# sorted by weight on component 2 (cytokine)
sorted_ids_cytok2 <- sort(res_diablo$variates$cytok[, "comp2"])
sorted_ids_cytok2_times <- dex_tp_data_overlap_with_dt_char[paste0("MVIR1-", names(sorted_ids_cytok2)), c("timepoint", "dt", "dt_char", "dex")]
sorted_ids_cytok2_times <- cbind(sorted_ids_cytok2_times, comp2=sort(res_diablo$variates$cytok[, "comp2"]))
sorted_ids_cytok2_times_filt <- sorted_ids_cytok2_times[sorted_ids_cytok2_times$dex, ]

# sorted by weight on component 1 (pbmc)
sorted_ids_pbmc1 <- sort(res_diablo$variates$PBMC[, "comp1"])
sorted_ids_pbmc1_times <- dex_tp_data_overlap_with_dt_char[paste0("MVIR1-", names(sorted_ids_pbmc1)), c("timepoint", "dt", "dt_char", "dex")]
sorted_ids_pbmc1_times <- cbind(sorted_ids_pbmc1_times, comp1=sort(res_diablo$variates$PBMC[, "comp1"]))
sorted_ids_pbmc1_times_filt <- sorted_ids_pbmc1_times[sorted_ids_pbmc1_times$dex, ]

# sorted by weight on component 2 (pbmc)
sorted_ids_pbmc2 <- sort(res_diablo$variates$PBMC[, "comp2"])
sorted_ids_pbmc2_times <- dex_tp_data_overlap_with_dt_char[paste0("MVIR1-", names(sorted_ids_pbmc2)), c("timepoint", "dt", "dt_char", "dex")]
sorted_ids_pbmc2_times <- cbind(sorted_ids_pbmc2_times, comp2=sort(res_diablo$variates$PBMC[, "comp2"]))
sorted_ids_pbmc2_times_filt <- sorted_ids_pbmc2_times[sorted_ids_pbmc2_times$dex, ]

pdf(paste0(results_path, paste0(paste0("noeta_all_indiv_diablo_", paste(val_, vent_status, sep="_")), ".pdf")), width = 8, height = 4)
plotIndiv(res_diablo, legend=TRUE, ind.names = F, group=dex_vect, pch=rep(1, length(dex_vect)), 
          size.title=18, size.subtitle=18, size.xlabel=18, size.ylabel=18, 
          size.axis=18, size.legend=18, size.legend.title=18)
dev.off()

pdf(paste0(results_path, paste0(paste0("noeta_all_loadings_1_diablo_", paste(val_, vent_status, sep="_")), ".pdf")), width = 8, height = 4)
plotLoadings(res_diablo, comp = 1, contrib = "max", legend = FALSE, method="median", 
             size.name  =1, xlab="Loadings")
dev.off()

pdf(paste0(results_path, paste0(paste0("noeta_all_loadings_2_diablo_", paste(val_, vent_status, sep="_")), ".pdf")), width = 8, height = 4)
plotLoadings(res_diablo, comp = 2, contrib = "max", legend = FALSE, method="median")
dev.off()

pdf(paste0(results_path, paste0(paste0("noeta_all_aucroc_cytok_1_diablo_", paste(val_, vent_status, sep="_")), ".pdf")), width = 7, height = 5)
auroc(res_diablo, roc.block = "cytok", roc.comp = 1)
dev.off()

pdf(paste0(results_path, paste0(paste0("noeta_all_aucroc_cytok_2_diablo_", paste(val_, vent_status, sep="_")), ".pdf")), width = 7, height = 5)
auroc(res_diablo, roc.block = "cytok", roc.comp = 2)
dev.off()

pdf(paste0(results_path, paste0(paste0("noeta_all_aucroc_pbmc_1_diablo_", paste(val_, vent_status, sep="_")), ".pdf")), width = 7, height = 5)
auroc(res_diablo, roc.block = "PBMC", roc.comp = 1)
dev.off()

pdf(paste0(results_path, paste0(paste0("noeta_all_aucroc_pbmc_2_diablo_", paste(val_, vent_status, sep="_")), ".pdf")), width = 7, height = 5)
auroc(res_diablo, roc.block = "PBMC", roc.comp = 2)
dev.off()

# GSEA
pbmc_1 <- res_diablo$loadings$PBMC[, "comp1", drop=F]
pbmc_2 <- res_diablo$loadings$PBMC[, "comp2", drop=F]
pbmc_1 <- pbmc_1[order(pbmc_1, decreasing = T), ]
pbmc_2 <- pbmc_2[order(pbmc_2, decreasing = T), ]

for (m_db_id in names(m_dbs)){
  print(m_db_id)
  m_db <- m_dbs[[m_db_id]]
  m_db_list <- lapply(unique(m_db$gs_name), function(x) unname(unlist(m_db[m_db$gs_name == x, "gene_symbol"])))
  names(m_db_list) <- unique(m_db$gs_name)
  
  set.seed(123)
  pbmc_1_gsea <- fgsea(m_db_list, stats=pbmc_1, nPermSimple=10000)
  pbmc_2_gsea <- fgsea(m_db_list, stats=pbmc_2, nPermSimple=10000)
  
  pbmc_1_gsea[pbmc_1_gsea$padj<0.1, "pathway"]$pathway
  pbmc_2_gsea[pbmc_2_gsea$padj<0.1, "pathway"]$pathway
  
  if (length(pbmc_1_gsea[pbmc_1_gsea$padj<0.1, "pathway"]$pathway)>0){
    pbmc_1_gsea_df <- data.frame(pbmc_1_gsea)
    pbmc_1_gsea_df$leadingEdge <- as.character(pbmc_1_gsea_df$leadingEdge)
    write.csv(pbmc_1_gsea_df, paste0(results_path, paste0(m_db_id, paste0(paste0("_noeta_gsea_pbmc_1_", paste(val_, vent_status, sep="_")), ".csv"))))
    
    # dot plots
    # PBMC 1
    pbmc_1_gsea_df_pos <- pbmc_1_gsea_df[pbmc_1_gsea_df$NES > 0 & pbmc_1_gsea_df$padj < 0.1, ]
    pbmc_1_gsea_df_neg <- pbmc_1_gsea_df[pbmc_1_gsea_df$NES < 0 & pbmc_1_gsea_df$padj < 0.1, ]
    
    pbmc_1_gsea_df_pos <- pbmc_1_gsea_df_pos[order(pbmc_1_gsea_df_pos$padj, decreasing = F), ]
    pbmc_1_gsea_df_neg <- pbmc_1_gsea_df_neg[order(pbmc_1_gsea_df_neg$padj, decreasing = F), ]
    
    pbmc_1_gsea_df_pos_top10 <- as.character(pbmc_1_gsea_df_pos[1:(min(nrow(pbmc_1_gsea_df_pos), 10)), "pathway"])
    pbmc_1_gsea_df_neg_top10 <- as.character(pbmc_1_gsea_df_neg[1:(min(nrow(pbmc_1_gsea_df_neg), 10)), "pathway"])
    
    pbmc_1_gsea_df <- pbmc_1_gsea_df[order(pbmc_1_gsea_df$NES, decreasing=T), ]
    pbmc_1_gsea_df_filt <- pbmc_1_gsea_df[pbmc_1_gsea_df$pathway %in% c(pbmc_1_gsea_df_pos_top10, pbmc_1_gsea_df_neg_top10), ]
    
    # rename pathways
    pbmc_1_gsea_df_filt$pathway <- str_wrap(sapply(pbmc_1_gsea_df_filt$pathway, function(x) StrCap(tolower(str_replace(str_replace_all(x, "_", " "), "REACTOME ", "")), "title")), 50, exdent = 2)

    pdf(paste0(results_path, paste0(m_db_id, paste0(paste0("_noeta_dotplot_gsea_pbmc_1_", paste(val_, vent_status, sep="_")), ".pdf"))), 
        width = 9, height = 8)
    print(ggplot(pbmc_1_gsea_df_filt, aes(x= NES, y=pathway)) +
            geom_point(aes(colour=NES), size=3) +
            scale_colour_gradient2(
              low = "#388ecc",
              mid = "white",
              high = "#f68b33",
              midpoint = 0,
              aesthetics = "colour"
            ) + 
            theme_bw() + theme(text = element_text(size = 18)) +
            xlab("NES") +
            ylab(paste0(m_db_id, " gene set")) +
            scale_y_discrete(limits = rev(pbmc_1_gsea_df_filt$pathway))+
            geom_vline(xintercept=0, linetype="dashed") +
            geom_stripes(odd = "#33333333", even = "#00000000"))
    dev.off()
    
    # NES sorted
    pbmc_1_gsea_df_pos <- pbmc_1_gsea_df_pos[order(pbmc_1_gsea_df_pos$NES, decreasing = T), ]
    pbmc_1_gsea_df_neg <- pbmc_1_gsea_df_neg[order(pbmc_1_gsea_df_neg$NES, decreasing = F), ]
    
    pbmc_1_gsea_df_pos_top10 <- as.character(pbmc_1_gsea_df_pos[1:(min(nrow(pbmc_1_gsea_df_pos), 10)), "pathway"])
    pbmc_1_gsea_df_neg_top10 <- as.character(pbmc_1_gsea_df_neg[1:(min(nrow(pbmc_1_gsea_df_neg), 10)), "pathway"])
    
    pbmc_1_gsea_df <- pbmc_1_gsea_df[order(pbmc_1_gsea_df$NES, decreasing=T), ]
    pbmc_1_gsea_df_filt <- pbmc_1_gsea_df[pbmc_1_gsea_df$pathway %in% c(pbmc_1_gsea_df_pos_top10, pbmc_1_gsea_df_neg_top10), ]
    
    pdf(paste0(results_path, paste0(m_db_id, paste0(paste0("_noeta_dotplot_gsea_NES_sorted_pbmc_1_", paste(val_, vent_status, sep="_")), ".pdf"))), width = 10, height = 6)
    print(ggplot(pbmc_1_gsea_df_filt, aes(x= NES, y=pathway)) +
            geom_point(aes(colour=NES), size=3) +
            scale_colour_gradient2(
              low = "#388ecc",
              mid = "white",
              high = "#f68b33",
              midpoint = 0,
              aesthetics = "colour"
            ) + 
            theme_bw() +
            xlab("NES") +
            ylab(paste0(m_db_id, " gene set")) +
            scale_y_discrete(limits = rev(pbmc_1_gsea_df_filt$pathway)))
    dev.off()
  }
  
  if (length(pbmc_2_gsea[pbmc_2_gsea$padj<0.1, "pathway"]$pathway)>0){
    pbmc_2_gsea_df <- data.frame(pbmc_2_gsea)
    pbmc_2_gsea_df$leadingEdge <- as.character(pbmc_2_gsea_df$leadingEdge)
    write.csv(pbmc_2_gsea_df, paste0(results_path, paste0(m_db_id, paste0(paste0("_noeta_gsea_pbmc_2_", paste(val_, vent_status, sep="_")), ".csv"))))
    
    # PBMC 2
    pbmc_2_gsea_df_pos <- pbmc_2_gsea_df[pbmc_2_gsea_df$NES > 0 & pbmc_2_gsea_df$padj < 0.1, ]
    pbmc_2_gsea_df_neg <- pbmc_2_gsea_df[pbmc_2_gsea_df$NES < 0 & pbmc_2_gsea_df$padj < 0.1, ]
    
    pbmc_2_gsea_df_pos <- pbmc_2_gsea_df_pos[order(pbmc_2_gsea_df_pos$padj, decreasing = F), ]
    pbmc_2_gsea_df_neg <- pbmc_2_gsea_df_neg[order(pbmc_2_gsea_df_neg$padj, decreasing = F), ]
    
    pbmc_2_gsea_df_pos_top10 <- as.character(pbmc_2_gsea_df_pos[1:(min(nrow(pbmc_2_gsea_df_pos), 10)), "pathway"])
    pbmc_2_gsea_df_neg_top10 <- as.character(pbmc_2_gsea_df_neg[1:(min(nrow(pbmc_2_gsea_df_neg), 10)), "pathway"])
    
    pbmc_2_gsea_df <- pbmc_2_gsea_df[order(pbmc_2_gsea_df$NES, decreasing=T), ]
    pbmc_2_gsea_df_filt <- pbmc_2_gsea_df[pbmc_2_gsea_df$pathway %in% c(pbmc_2_gsea_df_pos_top10, pbmc_2_gsea_df_neg_top10), ]
    
    pdf(paste0(results_path, paste0(m_db_id, paste0(paste0("_noeta_dotplot_gsea_pbmc_2_", paste(val_, vent_status, sep="_")), ".pdf"))), width = 10, height = 6)
    print(ggplot(pbmc_2_gsea_df_filt, aes(x= NES, y=pathway)) +
            geom_point(aes(colour=NES), size=3) +
            scale_colour_gradient2(
              low = "#388ecc",
              mid = "white",
              high = "#f68b33",
              midpoint = 0,
              aesthetics = "colour"
            ) + 
            theme_bw() +
            xlab("NES") +
            ylab(paste0(m_db_id, " gene set")) +
            scale_y_discrete(limits = rev(pbmc_2_gsea_df_filt$pathway)))
    dev.off()
  }
}











