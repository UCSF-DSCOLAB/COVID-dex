######################################################
# 1.cytok_analysis-larger_cohort_v2.R
# created on August 5 2021
# lucile.neyton@ucsf.edu

# This script aims at analysing cytokine data for dex vs on-dex patients
######################################################

rm(list = ls())
setwd("/Users/lucileneyton/Box Sync/COMET/data_analysis/")

# libraries
library(magmaR)
library(stringr)
library(ggrepel)
library(scales)
source('/Users/lucileneyton/Box Sync/tools/R_tools/distrib_utils.R')

# data files
data_path <- "/Users/lucileneyton/Box Sync/COMET/data/larger_cohort/"
results_path <- "/Users/lucileneyton/Box Sync/COMET/results/larger_cohort/"

#########################
# FUNCTIONS
#########################
# figs only
plot_distrib_continuous_nop <- function(meta_data, var_name, class_name,
                                        diff_test=c("param", "non-param"), 
                                        clust_folder, logscale=F, logscale_val=NULL) {
  #check parameters
  diff_test <- match.arg(diff_test)
  
  # extract data of interest
  var_data <- meta_data[, c(var_name, class_name)]
  colnames(var_data)[1] <- "var"
  colnames(var_data)[2] <- "cluster_label"
  
  # colors
  col_levels <- levels(var_data$cluster_label)
  colors_ <- hue_pal()(length(col_levels))
  # sort alphabetically
  names(colors_) <- sort(col_levels)
  
  print(col_levels)
  print(colors_)
  
  # generate a boxplot
  pdf(paste(clust_folder, paste(var_name, "boxplot.pdf", sep = "_"), sep = ""),
      width = 5, height = 5)
  
  if (length(col_levels)==2){
    print(ggplot(var_data, aes(x = cluster_label, y = var,
                               fill = cluster_label)) +
            ggtitle(paste(var_name, "distribution")) +
            geom_boxplot() +
            geom_point() +
            labs(x = class_name) +
            scale_fill_manual(values=c("dex" = "#f68b33", "nondex" = "#388ecc"), breaks=col_levels) + 
            theme_bw()
    )
  }else{
    print(ggplot(var_data, aes(x = cluster_label, y = var,
                               fill = cluster_label)) +
            ggtitle(paste(var_name, "distribution")) +
            geom_boxplot() +
            geom_point() +
            labs(x = class_name) +
            theme_bw()
    )
  }
  
  dev.off()
}

#########################
# DATA LOADING
#########################
# analyte data
# data loading
# download analyte tables from COMET DL
# prod <- magmaRset()
# analyte_data <- retrieve(target = prod,
#                      projectName = "mvir1",
#                      modelName = "analyte",
#                      recordNames = "all",
#                      attributeNames = "all",
#                      filter = "")

# save data to fix version
# date_ <- Sys.Date()
# saveRDS(analyte_data, paste(data_path, paste(paste("raw/analyte_data_", date_, sep=""), ".RDS", sep=""), sep = ""))

#date_ <- "2022-05-31"
date_ <- "2022-08-24"
analyte_data <- readRDS(paste(data_path, paste(paste("raw/analyte_data_", date_, sep=""), ".RDS", sep=""), sep = ""))

# metadata
metadata_df <- readRDS(paste0(data_path, "processed/df_metadata.RDS"))

# collection times and dex ttt info
dex_tp_data <- readRDS(paste0(data_path, "processed/dex_tp_data_unfiltered_with_impacc_status.RDS"))

# intubation data
vfd_data_filt <- readRDS(paste0(data_path, "processed/vfd_data_filt.RDS"))

# store unfiltered data
analyte_data_full <- analyte_data

#########################
# DATA PRE-PROCESSING
#########################
# add patient identifiers
analyte_data <- analyte_data[grepl("-D0PL1-CTK1", analyte_data$immunoassay), ]
analyte_data$patient <- str_replace(analyte_data$immunoassay, "-D0PL1-CTK1", "")

# filtered data
dex_tp_data$overlap[(is.na(dex_tp_data$overlap)) & (dex_tp_data$dex)] <- FALSE
dex_tp_data_overlap <- dex_tp_data[(!(dex_tp_data$dex)) | (dex_tp_data$overlap), ]
dex_tp_data_overlap <- dex_tp_data_overlap[dex_tp_data_overlap$data_type == "Cytokine", ]

# dex vs non-dex
table(dex_tp_data_overlap$dex)

# difference between sample collection day and first dex dose for dex
dex_tp_data_overlap$dt <- as.Date(dex_tp_data_overlap$collection_date) - dex_tp_data_overlap$earliest_dex
table(dex_tp_data_overlap$dt)

# differences using collection and treatment times
dex_tp_data_overlap$collection_date_char <- as.POSIXlt(dex_tp_data_overlap$collection_date_char, format="%Y-%m-%d %H:%M")
dex_tp_data_overlap$earliest_dex_char <- as.POSIXlt(dex_tp_data_overlap$earliest_dex_char, format="%Y-%m-%d %H:%M")
dex_tp_data_overlap$dt_char <- dex_tp_data_overlap$collection_date_char - dex_tp_data_overlap$earliest_dex_char

# filter analyte data to only keep selected patients
analyte_data <- analyte_data[analyte_data$patient %in% dex_tp_data_overlap$patient, ]
analyte_data$timepoint <- sapply(analyte_data$immunoassay, function(x) str_replace(x, "PL1-CTK1", ""))
dex_tp_data_overlap <- dex_tp_data_overlap[dex_tp_data_overlap$timepoint %in% analyte_data$timepoint, ]

# TNF-R1 reported as TNF RI and TNF R1
analyte_data$analyte_name[analyte_data$analyte_name == "TNF RI"] <- "TNF R1"
analyte_data <- analyte_data[(analyte_data$analyte_name != "TNF R1") | ((analyte_data$analyte_name == "TNF R1") & (!is.na(analyte_data$value))), ]

# reformat to have one column per analyte
analyte_data <- dcast(analyte_data, patient ~ analyte_name)
rownames(analyte_data) <- analyte_data$patient
analyte_data <- analyte_data[, colnames(analyte_data)!="patient"]

# store cyto variable names
cyto_vars <- colnames(analyte_data)

# add metadata
analyte_data <- cbind(analyte_data, metadata_df[rownames(analyte_data), ])

# list dex and non-dex patients
dex_patients <- dex_tp_data_overlap[dex_tp_data_overlap$dex, "patient"]
nondex_patients <- dex_tp_data_overlap[!dex_tp_data_overlap$dex, "patient"]

# list less than 24hrs and more than 24 hours patients
h24_patients <- dex_tp_data_overlap[(!is.na(dex_tp_data_overlap$dt_char)) & dex_tp_data_overlap$dt_char<=24	, "patient"]
h24plus_patients <- dex_tp_data_overlap[(!is.na(dex_tp_data_overlap$dt_char)) & dex_tp_data_overlap$dt_char>24, "patient"]

saveRDS(dex_tp_data_overlap, paste0(results_path, "dex_tp_data_overlap_with_dt_char.RDS"))

# dex timing
analyte_data$dex_hrs[rownames(analyte_data) %in% nondex_patients] <- NA
analyte_data$dex_hrs[rownames(analyte_data) %in% h24_patients] <- "H24"
analyte_data$dex_hrs[rownames(analyte_data) %in% h24plus_patients] <- "H24+"

# add dex status
analyte_data$dex <- rownames(analyte_data) %in% dex_patients

# were patients intubated when samples were collected?
# filter for D0 cytokine samples
dex_tp_data_all <- dex_tp_data[dex_tp_data$data_type=="Cytokine" & grepl("D0", dex_tp_data$timepoint), ]

intub_vect <- c()
for (patient_ in rownames(analyte_data)){
  patient_longid <- as.numeric(strsplit(patient_, "-HS")[[1]][2]) + 1000
  
  patient_col <- dex_tp_data_all[dex_tp_data_all$patient==patient_, "collection_date"]
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
  
  intub_vect <- c(intub_vect, intub_status)
}

analyte_data$intub <- intub_vect

analyte_data$sample <- paste0(rownames(analyte_data), "-D0")

# new vars
non_cyto_vars <- colnames(analyte_data)[!(colnames(analyte_data) %in% cyto_vars)]

# only intubated
analyte_data <- analyte_data[analyte_data$intub, ]
saveRDS(analyte_data[, c("sample", "dex", "intub")], paste0(results_path, "cytok_data_for_dems.RDS"))

#########################
# DATA ANALYSIS - LINEAR MODELS
#########################
# for each cytokine
# dex vs non-dex
# non-normally distributed
wilcox_tests_res <- c()
log2fc_res <- c()
for (cyto_ in cyto_vars){
  wilcox_tests_res[cyto_] <- wilcox.test(analyte_data[rownames(analyte_data) %in% dex_patients , cyto_], analyte_data[rownames(analyte_data) %in% nondex_patients , cyto_])[['p.value']]
  log2fc_res[cyto_] <- log2(mean(analyte_data[rownames(analyte_data) %in% dex_patients , cyto_], na.rm=T)/mean(analyte_data[rownames(analyte_data) %in% nondex_patients , cyto_], na.rm=T))
}

wilcox_tests_res_cor <- p.adjust(wilcox_tests_res, "BH")

wilcox_tests_res_cor[wilcox_tests_res_cor<0.1]

# volcano plot
res_cytok_df <- data.frame(cytok=cyto_vars, log2fc= log2fc_res, padj= wilcox_tests_res_cor)
res_cytok_df$sig <- res_cytok_df$padj < 0.1

pdf(paste0(results_path, "cytok_volcano_plot.pdf"), width = 7, height = 7)
p <- ggplot(res_cytok_df, aes(log2fc, -log10(padj))) +
  geom_point(aes(col = sig)) +
  scale_color_manual(values = c("black", "red")) +
  theme(legend.position = "none") +
  geom_text_repel(data = res_cytok_df,
                  aes(label = cytok),
                  max.overlaps=20) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value") +
  theme_bw()
print(p)
dev.off()

# boxplots
# dex status and time point
analyte_data$cat3[rownames(analyte_data) %in% intersect(dex_patients, h24_patients)] <- "dex_h24"
analyte_data$cat3[rownames(analyte_data) %in% intersect(dex_patients, h24plus_patients)] <- "dex_h24plus"
analyte_data$cat3[rownames(analyte_data) %in% nondex_patients] <- "nondex"

analyte_data$cat2[rownames(analyte_data) %in% dex_patients] <- "dex"
analyte_data$cat2[rownames(analyte_data) %in% nondex_patients] <- "nondex"

analyte_data$cat2 <- factor(analyte_data$cat2, levels=c("nondex", "dex"))
analyte_data$cat3 <- factor(analyte_data$cat3, levels=c("nondex", "dex_h24", "dex_h24plus"))

for (cyto_ in cyto_vars){
  if (cyto_ %in% cyto_vars){
    analyte_data_tmp <- analyte_data
    analyte_data_tmp[, cyto_] <- log2(analyte_data_tmp[, cyto_] + 1)
    plot_distrib_continuous_nop(analyte_data_tmp, cyto_, "cat2", "non-param", paste0(results_path, "log2_cat2_"))
  }
}

nonna_samples <- rownames(analyte_data_tmp[!is.na(analyte_data_tmp$cat3), ])
for (cyto_ in cyto_vars){
  if (cyto_ %in% cyto_vars){
    analyte_data_tmp <- analyte_data
    analyte_data_tmp[, cyto_] <- log2(analyte_data_tmp[, cyto_] + 1)
    analyte_data_tmp <- analyte_data_tmp[nonna_samples, ]
    plot_distrib_continuous_nop(analyte_data_tmp, cyto_, "cat3", "non-param", paste0(results_path, "log2_cat3_"))
  }
}

#########################
# SAVE CYTOK DATA
#########################
saveRDS(analyte_data[, cyto_vars], paste0(data_path, "processed/cytok_data.RDS"))





