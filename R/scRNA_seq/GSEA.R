library(tidyverse)
library(BiocParallel)
library(fgsea)
MulticoreParam(6)
rm(list = ls())
#This script runs within a folder containing subfolders with DGEA results: 
#/COMETETA_DGEA
#/COMETWB_DGEA
#/LiaoBAL_DGEA
#/SinhaWB_DGEA
#Each subfolder contains csv files from Ravi's MAST analysis 

setwd("~/Research/COVIDDex/blueprintresults/") 
filenames <- dir(pattern = "_main_dgea", recursive = T)

gsea <- function(deseq.res, #A DESeq result object
                 gmtfile ="~/Research/data/shared/c2.cp.reactome.v7.1.symbols.gmt", #A GMT file containing gene signatures of interest. These can be downloaded from MSigDB
                 collapseCutoff = 0.1, 
                 collapse = FALSE,#Collapse similar pathways
                 rankstat = "l2f", 
                 nPerm = 10000, 
                 minSize = 25
){
  stopifnot(rankstat %in% c("stat", "l2f"))
  
  if(rankstat == "stat"){
    cat("Ranking genes by Wald statistic...\n")
    ranks <- deseq.res %>%
      as.data.frame %>%
      # rownames_to_column("ensembl_gene_id") %>%
      # left_join(ensembltohgnc) %>%
      dplyr::select(stat, hgnc_symbol) %>%
      na.omit() %>% 
      distinct() %>% 
      group_by(hgnc_symbol) %>% 
      summarize(stat=mean(stat)) %>% 
      arrange(desc(stat)) %>%
      deframe()
  }
  
  if(rankstat == "l2f"){
    cat("Ranking genes by log2-fold change...\n")
    ranks <- deseq.res %>%
      as.data.frame %>%
      # rownames_to_column("ensembl_gene_id") %>%
      # left_join(ensembltohgnc) %>%
      dplyr::select(log2FoldChange, hgnc_symbol) %>%
      na.omit() %>% 
      distinct() %>% 
      group_by(hgnc_symbol) %>% 
      summarize(log2FoldChange=mean(log2FoldChange)) %>% 
      arrange(desc(log2FoldChange)) %>%
      deframe()
  }
  
  pathways <<- gmtPathways(gmtfile)
  
  cat("Running fgsea...\n")
  fgseaRes <- fgseaMultilevel(pathways=pathways, stats=ranks,eps = 0, nPermSimple = nPerm, minSize = minSize)
  
  if(collapse == TRUE){
    cat("Collapsing pathways...\n")
    collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < collapseCutoff], 
                                          pathways, ranks)
    
    fgseaRes <- fgseaRes[pathway %in% collapsedPathways$mainPathways]
  }
  
  return(fgseaRes)
  
}

for(file in filenames){
  cellname <- file %>%
    str_remove("_main_dgea.csv") %>%
    str_remove("_DGEA\\/") %>%
    str_remove("\\+ ") %>%
    str_remove(" ") %>%
    str_remove("-")
  print(cellname)
  tmp <- read_csv(file) %>% 
    dplyr::rename(hgnc_symbol = 1, log2FoldChange = avg_log2FC, pvalue = 2, fdr = 4)
  assign(cellname, tmp)
  tmp <- gsea(tmp, gmtfile = "~/Research/data/shared/c5.go.bp.v7.5.1.symbols.gmt.txt")
  assign(paste(cellname, "GSEA", sep = "."), tmp)
}

gseaobjs <- ls(pattern = "*.GSEA")

gsea.merged <- NULL
for(obj in gseaobjs) {
  tmp <- get(obj)
  tmp <- cbind(tmp, obj) 
  gsea.merged <- rbind(gsea.merged,tmp)
}

saveRDS(gsea.merged, "~/Research/COVIDDex/merged.gsea.C5.Rds")

padj_cutoff <- 0.1
gsea.merged %>%
  mutate(celltype = str_remove(obj, "ETA") %>%
           str_remove("WB") %>%
           str_remove("BAL") %>%
           str_remove("COMET|Liao|Sinha") %>%
           str_remove(".GSEA"), 
         study = str_extract(obj, "COMETWB|COMETETA|LiaoBAL|SinhaWB"), 
         compartment = case_when(
           str_detect(study, "WB") ~ "Blood", 
           TRUE ~ "Lung"
         )) %>% 
  dplyr::select(celltype, study, compartment, pathway, NES, padj) %>%
  filter(!str_detect(celltype, "CLP|DC|Epithe|Macroph|Bcells|Megakaryo|Plasma|Tcells")) %>% 
  mutate(pathway = str_remove(pathway, "REACTOME_") %>% str_replace_all("_", " ")) %>%
  pivot_wider(id_cols = c(pathway, celltype), names_from = study, values_from = c(NES, padj)) -> 
  gsea.wide 

gsea.wide %>%
  filter(
    padj_COMETETA < padj_cutoff, padj_LiaoBAL < padj_cutoff, 
    padj_SinhaWB < padj_cutoff, padj_COMETWB < padj_cutoff
  ) %>%
  filter(
    sign(NES_COMETETA) == sign(NES_LiaoBAL), 
    sign(NES_COMETWB) == sign(NES_SinhaWB),
    sign(NES_COMETETA) != sign(NES_COMETWB)
  ) %>% write_csv("Pathways that differ between blood and lung.csv")

pdf(height = 15, width = 10, file = "MSigDB pathways that are concordant between lung and blood.pdf")
gsea.wide %>%
  filter(
    padj_COMETETA < padj_cutoff, padj_LiaoBAL < padj_cutoff, 
    padj_SinhaWB < padj_cutoff, padj_COMETWB < padj_cutoff
  ) %>%
  filter(
    sign(NES_COMETETA) == sign(NES_LiaoBAL), 
    sign(NES_COMETWB) == sign(NES_SinhaWB),
    sign(NES_COMETETA) == sign(NES_COMETWB)
  ) %>%
  mutate(ETAscore = NES_COMETETA) %>%
  pivot_longer(cols = matches("NES"), names_to = "cohort", values_to = "NES") %>%
  mutate(compartment = case_when(
    str_detect(cohort, "WB") ~ "Blood", 
    TRUE ~ "Lung"
  ))%>%
  mutate(cohort = case_when(
    str_detect(cohort, "COMET") ~ "COMET", 
    cohort == "NES_LiaoBAL" ~ "Liao", 
    cohort == "NES_SinhaWB" ~ "Sinha"
  )) %>%
  ggplot(aes(x = cohort, fill = NES, y = reorder(pathway, ETAscore))) + 
  facet_grid(celltype~compartment, scale = "free" , space = "free") + 
  geom_tile()+ 
  theme_classic() + 
  scale_fill_viridis_c(option = "A") + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 15), 
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust =1), 
        panel.border = element_rect(size = 2, fill = NA), 
        axis.line = element_blank()) + 
  labs(x = "Cohort", y = "MSigDB pathway", fill = "NES") 
dev.off()

pdf(height = 30, width = 10, file = "MSigDB pathways that are discordant between lung and blood.pdf")
gsea.wide %>%
  filter(
    padj_COMETETA < padj_cutoff, padj_LiaoBAL < padj_cutoff, 
    padj_SinhaWB < padj_cutoff, padj_COMETWB < padj_cutoff
  ) %>%
  filter(
    sign(NES_COMETETA) == sign(NES_LiaoBAL), 
    sign(NES_COMETWB) == sign(NES_SinhaWB),
    sign(NES_COMETETA) != sign(NES_COMETWB)
  ) %>%
  mutate(ETAscore = NES_COMETETA) %>%
  pivot_longer(cols = matches("NES"), names_to = "cohort", values_to = "NES") %>%
  mutate(compartment = case_when(
    str_detect(cohort, "WB") ~ "Blood", 
    TRUE ~ "Lung"
  ))%>%
  mutate(cohort = case_when(
    str_detect(cohort, "COMET") ~ "COMET", 
    cohort == "NES_LiaoBAL" ~ "Liao", 
    cohort == "NES_SinhaWB" ~ "Sinha"
  )) %>%
  ggplot(aes(x = cohort, fill = NES, y = reorder(pathway, ETAscore))) + 
  facet_grid(celltype~compartment, scale = "free" , space = "free") + 
  geom_tile()+ 
  theme_classic() + 
  scale_fill_viridis_c(option = "A") + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 15), 
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust =1), 
        panel.border = element_rect(size = 2, fill = NA), 
        axis.line = element_blank()) + 
  labs(x = "Cohort", y = "MSigDB pathway", fill = "NES") 
dev.off()

pdf(height = 200, width = 10, file = "MSigDB pathways that are reproducible in blood.pdf")
gsea.wide %>%
  filter(
    padj_SinhaWB < padj_cutoff, padj_COMETWB < padj_cutoff
  ) %>%
  filter(
    sign(NES_COMETWB) == sign(NES_SinhaWB),
  ) %>%
  mutate(ETAscore = NES_COMETETA) %>%
  pivot_longer(cols = matches("NES"), names_to = "cohort", values_to = "NES") %>%
  mutate(compartment = case_when(
    str_detect(cohort, "WB") ~ "Blood", 
    TRUE ~ "Lung"
  ))%>%
  mutate(cohort = case_when(
    str_detect(cohort, "COMET") ~ "COMET", 
    cohort == "NES_LiaoBAL" ~ "Liao", 
    cohort == "NES_SinhaWB" ~ "Sinha"
  )) %>%
  ggplot(aes(x = cohort, fill = NES, y = reorder(pathway, ETAscore))) + 
  facet_grid(celltype~compartment, scale = "free" , space = "free") + 
  geom_tile()+ 
  theme_classic() + 
  scale_fill_viridis_c(option = "A") + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 15), 
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust =1), 
        panel.border = element_rect(size = 2, fill = NA), 
        axis.line = element_blank()) + 
  labs(x = "Cohort", y = "MSigDB pathway", fill = "NES") 
dev.off()

pdf(height = 200, width = 10, file = "MSigDB pathways that are reproducible in lung.pdf")
gsea.wide %>%
  filter(
    padj_COMETETA < padj_cutoff, padj_LiaoBAL < padj_cutoff
  ) %>%
  filter(
    sign(NES_COMETETA) == sign(NES_LiaoBAL), 
  ) %>%
  mutate(ETAscore = NES_COMETETA) %>%
  pivot_longer(cols = matches("NES"), names_to = "cohort", values_to = "NES") %>%
  mutate(compartment = case_when(
    str_detect(cohort, "WB") ~ "Blood", 
    TRUE ~ "Lung"
  ))%>%
  mutate(cohort = case_when(
    str_detect(cohort, "COMET") ~ "COMET", 
    cohort == "NES_LiaoBAL" ~ "Liao", 
    cohort == "NES_SinhaWB" ~ "Sinha"
  )) %>%
  ggplot(aes(x = cohort, fill = NES, y = reorder(pathway, ETAscore))) + 
  facet_grid(celltype~compartment, scale = "free" , space = "free") + 
  geom_tile()+ 
  theme_classic() + 
  scale_fill_viridis_c(option = "A") + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 15), 
        axis.title = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 45, hjust = 1, vjust =1), 
        panel.border = element_rect(size = 2, fill = NA), 
        axis.line = element_blank()) + 
  labs(x = "Cohort", y = "MSigDB pathway", fill = "NES") 
dev.off()


gsea.wide <- gsea.wide %>% 
  mutate(
    sig_lung = 
      case_when(
        padj_COMETETA < padj_cutoff & padj_LiaoBAL < padj_cutoff ~ "Both",
        padj_COMETETA < padj_cutoff  ~ "COMET",
        padj_LiaoBAL < padj_cutoff ~ "Liao",
        TRUE ~ "None"
      ),
    sig_blood = 
      case_when(
        padj_COMETWB < padj_cutoff & padj_SinhaWB < padj_cutoff ~ "Both",
        padj_COMETWB < padj_cutoff  ~ "COMET",
        padj_SinhaWB < padj_cutoff ~ "Sinha",
        TRUE ~ "None"
      ), 
    alpha_lung = 
      case_when(
        sig_lung == "Both" ~ 1, 
        sig_lung == "None" ~ 0.1, 
        TRUE ~ 0.2
      ),
    alpha_blood = 
      case_when(
        sig_blood == "Both" ~ 1, 
        sig_blood == "None" ~ 0.1, 
        TRUE ~ 0.2
      )
  )

pdf(file = "COVID Dex - Pathways by cell and compartment.pdf", height = 10, width = 10)
gsea.wide %>%
  ggplot(aes(x = NES_COMETETA, y = NES_LiaoBAL, color = sig_lung, alpha = alpha_lung)) + 
  geom_point() + 
  theme_classic() + 
  facet_wrap(~celltype)

gsea.wide %>%
  ggplot(aes(x = NES_COMETWB, y = NES_SinhaWB, color = sig_blood, alpha = alpha_blood)) + 
  geom_point() + 
  theme_classic() + 
  facet_wrap(~celltype)
dev.off()

