suppressMessages({
  library(magmaR)
  library(tidyverse)
  library(RColorBrewer)
  library(cowplot)
  library(ggpubr)
  library(pheatmap)
  library(pals)
})

make_sample_timepoint_plot_v2 <- function(data, data_type="PBMC_RNASeq", redcap_data_mini, collection_date_column = "dos", return_df=F) {
  ea_rnaseq_patients = data %>% filter(data_type == !!data_type) %>% pull(patient) %>% unique
  ea_rnaseq_timepoints = data %>% filter(data_type == !!data_type) %>% pull(timepoint) %>% unique
  
  # Extract Dex start-end dates for the selected patient ids above
  ea_rnaseq_dates = redcap_data_mini %>% 
    filter( grepl("enrollment", redcap_event_name) & patient_id %in% ea_rnaseq_patients ) %>% 
    select(patient_id, meds_hosp_cs_dex_start, meds_hosp_cs_dex_end, meds_hosp_cs_po_dex_start, meds_hosp_cs_po_dex_end, hosp_dt, vent_start_date, vent_end_date, date_discharged) %>%
    ## Add dex status
    #mutate( dex_status = ifelse( (meds_hosp_cs_dex_start == "" | meds_hosp_cs_dex_end == "") & (meds_hosp_cs_po_dex_start == "" | meds_hosp_cs_po_dex_start == ""), FALSE, TRUE)  ) %>%
    inner_join( data %>% select(patient, dex) %>% unique(), by=c("patient_id"="patient") ) %>%
    inner_join( redcap_data_mini %>% filter(.data[[collection_date_column]] != "") %>% select(patient_id, !!collection_date_column, redcap_event_name) ) %>%
    mutate( tp = ifelse(redcap_event_name == "d1_arm_1b", "1", 
                        ifelse(redcap_event_name == "d1_arm_1", "-1", 
                               gsub("d", "", str_match(redcap_event_name, "d\\d+")[,1]) )) ) %>%
    mutate( dex_start = as.Date(meds_hosp_cs_dex_start) - as.Date(hosp_dt),
            dex_end = as.Date(meds_hosp_cs_dex_end) - as.Date(hosp_dt),
            #dex_p_start = as.Date(meds_hosp_cs_po_dex_start) - as.Date(hosp_dt),
            #dex_p_end = as.Date(meds_hosp_cs_po_dex_end) - as.Date(hosp_dt),
            sampling = as.Date(.data[[collection_date_column]]) - as.Date(hosp_dt),
            vent_start = as.Date(vent_start_date) - as.Date(hosp_dt),
            vent_end = as.Date(vent_end_date) - as.Date(hosp_dt),
            discharge = as.Date(date_discharged) - as.Date(hosp_dt)
    )
  
  ea_rnaseq_dates = ea_rnaseq_dates %>% filter( paste0(patient_id, "-D", tp) %in% ea_rnaseq_timepoints )
  
  ea_rnaseq_dates = ea_rnaseq_dates %>% 
    mutate( patient_id = factor(patient_id, levels = ea_rnaseq_dates %>% select(patient_id, hosp_dt) %>% unique() %>% arrange(as.Date(hosp_dt)) %>% pull(patient_id) ),
            tp = factor(tp, levels= as.character(sort(unique(as.numeric(tp)))) ))
  #return(ea_rnaseq_dates)
  dex_plot_df = ea_rnaseq_dates %>% 
    select( patient_id, starts_with("dex_")) %>%
    unique() %>%
    # Adding 1 to the "end"s, if the start and end dates are same.
    mutate( dex_end = ifelse(dex_start == dex_end, dex_end + 1, dex_end),
            #dex_p_end = ifelse(dex_p_start == dex_p_end, dex_p_end + 1, dex_p_end),
            dex_end = as.difftime(dex_end, units = "days"),
            #dex_p_end = as.difftime(dex_p_end, units = "days") 
    ) %>%
    pivot_longer( cols=-patient_id) %>%
    filter( ! is.na(value)) %>%
    dplyr::rename( dex_type = name, dex_starts_ends = value ) %>%
    mutate(  dex_type = gsub("_start|_end", "", dex_type) )
  
  
  vent_plot_df = ea_rnaseq_dates %>% 
    select( patient_id, vent_start, vent_end) %>%
    unique() %>%
    pivot_longer( cols=-patient_id) %>%
    filter( ! is.na(value)) %>%
    dplyr::rename( vent_type = name, vent_starts_ends = value )
  
  ea_rnaseq_dates = rbind( 
    cbind(data.frame(matrix(ncol = ncol(ea_rnaseq_dates)-1, nrow=nrow(dex_plot_df)) ) %>% setNames( grep("patient_id",colnames(ea_rnaseq_dates), value = T, invert = T) ), dex_plot_df),
    
    ea_rnaseq_dates = ea_rnaseq_dates %>% mutate( dex_type=NA, dex_starts_ends=NA))
  
  
  ea_rnaseq_dates = rbind( 
    cbind(data.frame(matrix(ncol = ncol(ea_rnaseq_dates)-1, nrow=nrow(vent_plot_df)) ) %>% setNames( grep("patient_id",colnames(ea_rnaseq_dates), value = T, invert = T) ), vent_plot_df),
    
    ea_rnaseq_dates = ea_rnaseq_dates %>% mutate( vent_type=NA, vent_starts_ends=NA))
  
  
  #pdf("~/../Desktop/Rplot.pdf", height = 11, width = 13)
  p = ggplot( ea_rnaseq_dates, aes(x=sampling, y=patient_id, group=patient_id, color=tp) ) + 
    geom_line(aes(x=vent_starts_ends, group=paste(patient_id) ), size=8, color="lightgrey" ) +
    geom_line(aes(x=dex_starts_ends, group=paste(patient_id,dex_type) ), size=4, color="lightpink" ) +
    #geom_point(aes(x=dex_start), pch=0, color="black", size=3) + 
    #geom_point(aes(x=dex_end), pch=15, color="black", size=3, alpha=0.2) + 
    #geom_point(aes(x=dex_p_start), pch=0, color="blue", size=3) + 
    #geom_point(aes(x=dex_p_end), pch=15, color="blue", size=3, alpha=0.2) + 
    #geom_point(aes(x=vent_start), pch=5, color="black", size=2) + 
    #geom_point(aes(x=vent_end), pch=18, color="black", size=4, alpha=0.2) + 
    geom_point(aes(x=-25,color=factor(as.character(as.numeric(dex)))), size=3) +
    geom_line(color="black") + geom_point(aes(shape=tp, fill=tp),size=3, color="black") + 
    theme_classic() + 
    scale_shape_manual(values=c(21,22,23,24)) +
    xlab("Days since hosp.") +
    ylab("Patients (ordered by hosp.)") +
    scale_fill_manual(values = pals::kelly(10)[-1]) +
    scale_color_manual( values = c(brewer.paired(12),alphabet2(26)) %>% setNames(NULL) )
  #dev.off()
  print(p)
  
  if(return_df)
    return(ea_rnaseq_dates)
}

## Load data from Lucile and from RedCap
data = readRDS("dex_tp_data_unfiltered_with_impacc_status_20220527.RDS")
redcap_data = read.csv("COMETPatientAndSampl_DATA_2022-04-25_1215_screening_enrollConcent_demog_covidtest_admission_meds_smplColl_ards.csv", fileEncoding = "UTF-8-BOM")
studyTracker = read.csv("COMETPatientAndSampl_DATA_2021-10-01_1455_studyTracker.csv", fileEncoding = "UTF-8-BOM")
vent_discharge_data = studyTracker %>% select(record_id, vent_start_date, vent_end_date, contains("tubation"), date_discharged) %>% filter( vent_start_date != "" | date_discharged != "" )

## Download patient data from Data Library
TOKEN = "GET_TOKEN_FROM_JANUS"
prod <- magmaRset(token = TOKEN)

patient <- retrieve(
  target = prod,
  projectName = "mvir1",
  modelName = "patient",
  recordNames = "all",
  attributeNames = "all",
  filter = "")
patient <- patient %>% mutate( death = ifelse(age_at_death == 0, FALSE, TRUE) )
patient$comet_id_2 = sapply(patient$comet_id, function(x){ val=unlist(str_split(x,"\\."))[2]; str_pad(val, width=4, side = "right", pad = "0") } )




## Combine all of the data loaded above
t = inner_join( data %>% select(patient, timepoint, dex), patient %>% select(name, comet_id_2), by=c("patient"="name") ) %>% unique() %>% mutate(timepoint_d = gsub("MVIR1-HS\\d+-","",timepoint) )

redcap_data$patient_id = patient$name[ match(redcap_data$record_id, patient$comet_id_2) ]

redcap_data_mini = redcap_data %>% select(c("patient_id","record_id", "dos", "date_swab", "dos_esc_ta","redcap_event_name", "meds_hosp_cs_dex_start", "meds_hosp_cs_dex_end", "meds_hosp_cs_po_dex_start", "meds_hosp_cs_po_dex_end", "meds_hosp_dex_total_start", "meds_hosp_dex_total_end", "inpatient_medications_complete", "hosp_dt")) %>%
  filter( grepl("arm_1", redcap_event_name) ) %>%
  left_join( vent_discharge_data, by="record_id")

# Copy the meds_hosp_cs_po_dex_end dates to meds_hosp_cs_dex_end for HS233
pati_enroll_row_idx = which(redcap_data_mini$patient_id == "MVIR1-HS233" & grepl("enrollment", redcap_data_mini$redcap_event_name))
redcap_data_mini$meds_hosp_cs_dex_start[ pati_enroll_row_idx ] = redcap_data_mini$meds_hosp_cs_po_dex_start[ pati_enroll_row_idx ]
redcap_data_mini$meds_hosp_cs_dex_end[ pati_enroll_row_idx ] = redcap_data_mini$meds_hosp_cs_po_dex_end[ pati_enroll_row_idx ]

##### ETA scRNA sample plot
# Final list.
ea_scRNA_ids_aartik = c("MVIR1-HS72-D4","MVIR1-HS414-D2","MVIR1-HS1-D1","MVIR1-HS2-D1","MVIR1-HS24-D0","MVIR1-HS50-D0","MVIR1-HS89-D1","MVIR1-HS107-D0","MVIR1-HS161-D4","MVIR1-HS340-D1","MVIR1-HS357-D2","MVIR1-HS389-D2","MVIR1-HS398-D4","MVIR1-HS403-D2","MVIR1-HS410-D1","MVIR1-HS415-D4", "MVIR1-HS471-D2")

# Add vent start and end dates for HS471 based on Lucile's med. record review.
redcap_data_mini$vent_start_date[ redcap_data_mini$patient_id == "MVIR1-HS471" ] = "2021-07-12"
redcap_data_mini$vent_end_date[ redcap_data_mini$patient_id == "MVIR1-HS471" ] = "2021-08-10"

tdat = data %>% filter(timepoint %in% ea_scRNA_ids_aartik)
pdf("eta_scrna_seq.pdf", height = 5.5, width = 6.5)
d = make_sample_timepoint_plot_v2(tdat, "EA_scRNASeq", redcap_data_mini, "dos_esc_ta", return_df = T)
dev.off()

# Export data points for Nat Comms.
dplyr::select(d, all_of(c("dex", "dex_type", "dex_starts_ends", "vent_starts_ends", "tp", "patient_id", "sampling"))) %>% 
  write.csv("eta_scrna_seq_plotdata.csv", row.names = F)

##### WB scRNA sample plot
# Make a sample plot for the selected samples.
tdat = read.csv("wb_vent_patients_final.csv") %>%
  setNames(c("patient","dex")) %>%
  mutate(data_type="WB_scRNASeq", timepoint = paste0(patient, "-D0"))

pdf("wb_scrna_seq.pdf", height = 5.5, width = 6.5)
d = make_sample_timepoint_plot_v2(tdat, "WB_scRNASeq", redcap_data_mini, "dos", return_df = T)
dev.off()

# Export data points for Nat Comms.
dplyr::select(d, all_of(c("dex", "dex_type", "dex_starts_ends", "vent_starts_ends", "tp", "patient_id", "sampling"))) %>% 
  write.csv("wb_scrna_seq_plotdata.csv", row.names = F)


