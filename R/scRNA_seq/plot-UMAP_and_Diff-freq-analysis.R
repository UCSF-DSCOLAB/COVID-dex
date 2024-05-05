library(Seurat)
library(tidyverse)
library(ggrepel)
library(pals)
eta_rdata = "../ETA_scRNA/scRNA_analysis/out_blueprint_20220927/covid.combined_2022-09-28.Rdata" 
wb_rdata = "../WB_scRNA/scRNA_analysis/out_blueprint_20220927/covid.combined_2022-10-04.Rdata" 
sinha_rdata = "../natMed_scRNA/scRNA_analysis/out_blueprint_20220927/covid.combined_2022-10-04.Rdata" 
liao_rdata = "../natMed_liao_scRNA/scRNA_analysis/out_blueprint_20220927/covid.combined_2022-10-04.Rdata"

load(eta_rdata, verbose=T)
covid.combined.eta = covid.combined
load(wb_rdata, verbose=T)
covid.combined.wb = covid.combined
load(sinha_rdata, verbose=T)
covid.combined.sinha = covid.combined
load(liao_rdata, verbose=T)
covid.combined.liao = covid.combined

redo_major_labels = function(sobj) {
  celltype_change = read.csv("celltype_label_conversion_v2.csv")
  sobj$major.label = sobj$major.label.uncleaned
  sobj$major.label = ifelse( sobj$major.label %in% names(which(table(sobj$major.label) > 100)), sobj$major.label, "other")
  sobj$major.label = celltype_change$major.label2[ match(sobj$major.label, celltype_change$major.label) ]
  if(sum(is.na(sobj$major.label)) > 0) {
   stop("NAs were introduced in the cleaned up major.labels. Fix the issue and rerun.")
  }
  return(sobj)
}

covid.combined.eta = redo_major_labels(covid.combined.eta)
covid.combined.wb = redo_major_labels(covid.combined.wb)
covid.combined.sinha = redo_major_labels(covid.combined.sinha)
covid.combined.liao = redo_major_labels(covid.combined.liao)


# Make UMAPs
ctypes = unique(c(covid.combined.eta$major.label, covid.combined.wb$major.label))
cols = c("other"=kelly(1), kelly(length(ctypes))[-1] %>% setNames(ctypes[ctypes != "other"]) )
pdf("final_umaps_v2.pdf", width=8, height=7)
DimPlot(covid.combined.eta, pt.size = 0.1, group.by = "major.label", label = T) + scale_color_manual(values = cols) + ggtitle("ETA")
DimPlot(covid.combined.eta, pt.size = 0.1, group.by = "dex", label = T) + ggtitle("ETA") + scale_color_manual(values = c("TRUE"="#F68B33","FALSE"="#388ECD"))
DimPlot(covid.combined.wb, pt.size = 0.1, group.by = "major.label", label = T) + scale_color_manual(values = cols) + ggtitle("WB")
DimPlot(covid.combined.wb, pt.size = 0.1, group.by = "dex", label = T) + ggtitle("WB") + scale_color_manual(values = c("TRUE"="#F68B33","FALSE"="#388ECD"))
dev.off()

FetchData(covid.combined.wb, c("UMAP_1","UMAP_2","major.label","dex")) %>% write.csv("Fig3CE.csv")
FetchData(covid.combined.eta, c("UMAP_1","UMAP_2","major.label","dex")) %>% write.csv("Fig3DF.csv")





### Differential frequency analysis
removeNeuts = TRUE
# ETA
covid.combined_sub <- subset(covid.combined.eta, cells=Cells(covid.combined.eta)[ covid.combined.eta$major.label %in% names(which(table(covid.combined.eta$major.label) > 100))] )
if(removeNeuts)
  covid.combined_sub = subset(covid.combined_sub, cells=Cells(covid.combined_sub)[covid.combined_sub$major.label != "Neutrophils"] )
freq <- as.data.frame(table(covid.combined_sub@meta.data[ ,c("patient", "major.label")] )) %>%
  inner_join( (covid.combined_sub@meta.data %>% select(patient, dex) %>% unique()), by="patient"  )

norm_freq <- freq %>% group_by(patient) %>%
  mutate(Freq = Freq/sum(Freq)*100)

de_freq_eta = data.frame()
for(cell_t in unique(norm_freq$major.label)) {
  norm_freq_ct <- norm_freq %>% filter(major.label == cell_t)
  dex_grp = na.omit(norm_freq_ct$Freq[norm_freq_ct$dex])
  nondex_grp = na.omit(norm_freq_ct$Freq[!norm_freq_ct$dex])
  wt = wilcox.test(dex_grp, nondex_grp)
  de_freq_eta = rbind(de_freq_eta, data.frame(cell_t = cell_t, pval = wt$p.value, log2fc=log2( (mean(dex_grp)+0.01)/(mean(nondex_grp)+0.01)) ))
}


# Liao BAL
covid.combined_sub <- subset(covid.combined.liao, cells=Cells(covid.combined.liao)[ covid.combined.liao$major.label %in% names(which(table(covid.combined.liao$major.label) > 100))] )
if(removeNeuts)
  covid.combined_sub = subset(covid.combined_sub, cells=Cells(covid.combined_sub)[covid.combined_sub$major.label != "Neutrophils"] )
freq <- as.data.frame(table(covid.combined_sub@meta.data[ ,c("patient", "major.label")] )) %>%
  inner_join( (covid.combined_sub@meta.data %>% select(patient, dex) %>% unique()), by="patient"  )

norm_freq <- freq %>% group_by(patient) %>%
  mutate(Freq = Freq/sum(Freq)*100)

de_freq_liao = data.frame()
for(cell_t in unique(norm_freq$major.label)) {
  norm_freq_ct <- norm_freq %>% filter(major.label == cell_t)
  dex_grp = na.omit(norm_freq_ct$Freq[norm_freq_ct$dex])
  nondex_grp = na.omit(norm_freq_ct$Freq[!norm_freq_ct$dex])
  wt = wilcox.test(dex_grp, nondex_grp)
  de_freq_liao = rbind(de_freq_liao, data.frame(cell_t = cell_t, pval = wt$p.value, log2fc=log2( (mean(dex_grp)+0.01)/(mean(nondex_grp)+0.01)) ))
}



# WB
covid.combined_sub <- subset(covid.combined.wb, cells=Cells(covid.combined.wb)[ covid.combined.wb$major.label %in% names(which(table(covid.combined.wb$major.label) > 100))] )
if(removeNeuts)
  covid.combined_sub = subset(covid.combined_sub, cells=Cells(covid.combined_sub)[covid.combined_sub$major.label != "Neutrophils"] )
freq <- as.data.frame(table(covid.combined_sub@meta.data[ ,c("patient", "major.label")] )) %>%
  inner_join( (covid.combined_sub@meta.data %>% select(patient, dex) %>% unique()), by="patient"  )

# Remove neutrophils from the freq analysis
freq <- freq %>% filter(! major.label %in% c("Neutrophils","Megakaryocytes"))

norm_freq <- freq %>% group_by(patient) %>%
  mutate(Freq = Freq/sum(Freq)*100)

de_freq_wb = data.frame()
for(cell_t in unique(norm_freq$major.label)) {
  norm_freq_ct <- norm_freq %>% filter(major.label == cell_t)
  dex_grp = na.omit(norm_freq_ct$Freq[norm_freq_ct$dex])
  nondex_grp = na.omit(norm_freq_ct$Freq[!norm_freq_ct$dex])
  wt = wilcox.test(dex_grp, nondex_grp)
  de_freq_wb = rbind(de_freq_wb, data.frame(cell_t = cell_t, pval = wt$p.value, log2fc=log2( (mean(dex_grp)+0.01)/(mean(nondex_grp)+0.01)) ))
}



# Sinha
covid.combined_sub <- subset(covid.combined.sinha, cells=Cells(covid.combined.sinha)[ covid.combined.sinha$major.label %in% names(which(table(covid.combined.sinha$major.label) > 100))] )
if(removeNeuts)
  covid.combined_sub = subset(covid.combined_sub, cells=Cells(covid.combined_sub)[covid.combined_sub$major.label != "Neutrophils"] )
freq <- as.data.frame(table(covid.combined_sub@meta.data[ ,c("patient", "major.label")] )) %>%
  inner_join( (covid.combined_sub@meta.data %>% select(patient, dex) %>% unique()), by="patient"  )

# Remove neutrophils from the freq analysis
freq <- freq %>% filter(! major.label %in% c("Neutrophils","Megakaryocytes"))

norm_freq <- freq %>% group_by(patient) %>%
  mutate(Freq = Freq/sum(Freq)*100)

de_freq_sinha = data.frame()
for(cell_t in unique(norm_freq$major.label)) {
  norm_freq_ct <- norm_freq %>% filter(major.label == cell_t)
  dex_grp = na.omit(norm_freq_ct$Freq[norm_freq_ct$dex])
  nondex_grp = na.omit(norm_freq_ct$Freq[!norm_freq_ct$dex])
  wt = wilcox.test(dex_grp, nondex_grp)
  de_freq_sinha = rbind(de_freq_sinha, data.frame(cell_t = cell_t, pval = wt$p.value, log2fc=log2( (mean(dex_grp)+0.01)/(mean(nondex_grp)+0.01)) ))
}






de_freq = rbind(de_freq_eta %>% mutate(tis="eta"), de_freq_wb %>% mutate(tis="wb"), de_freq_sinha %>% mutate(tis="sinha"), de_freq_liao %>% mutate(tis="liaoBAL"))
de_freq = de_freq %>% mutate(cell_t = factor(cell_t,
                             #levels=unique(de_freq$cell_t[order(de_freq$log2fc)])
                             levels=rev(c(de_freq[de_freq$tis=="wb",]$cell_t[order(de_freq[de_freq$tis=="wb",]$log2fc)], "Macrophages","DC"))
                        )) %>%
                      mutate(tis = factor(tis, levels = c("eta","liaoBAL","wb","sinha")))

# Remove "other"
de_freq = dplyr::filter(de_freq, cell_t != "other")

pdf("diff_frequency_analysis.pdf", height=6, width=7)

ggplot(de_freq, aes(log2fc, cell_t, group=tis, fill=tis, shape=tis, size=-log10(pval))) +
  geom_vline(xintercept = 0, linetype=2, alpha=0.5) +
  geom_vline(xintercept = -0.5, linetype=2, alpha=0.5) +
  geom_vline(xintercept = 0.5, linetype=2, alpha=0.5) +
  geom_point(color = "black", alpha=0.7) +
  scale_fill_manual(values = pals::kelly(6)[c(4:5, 4:5)] %>% setNames(c("wb","eta","sinha","liaoBAL")) ) +
  scale_shape_manual(values = c(21, 21, 23, 25) %>% setNames(c("wb","eta","sinha","liaoBAL")) ) +
  geom_point(data = de_freq %>% filter(pval < 0.1), color = "black", shape = 0, aes(size = -log10(pval) + 5) ) +
  theme_classic() +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  xlab("log2FC(dex/non-dex)")

dev.off()

