library(tidyverse)
library(ggrepel)
eta_dir = "../../data/dge_files/" # "../ETA_scRNA/scRNA_analysis/out_blueprint_20220622/"
wb_dir = "../../data/dge_files/" # "../WB_scRNA/scRNA_analysis/out_blueprint_20220531/


# fc_thresh.x is abs(log2FC) > 0.5
lfcCutoff = 0.5

eta_files = list.files(eta_dir, pattern="cometETA_covidDex-vs-covidNondex.*_main_dgea.csv")
wb_files = list.files(wb_dir, pattern="cometWB_covidDex-vs-covidNondex.*_main_dgea.csv")
eta_files = str_remove(eta_files, "cometETA_")
wb_files = str_remove(wb_files, "cometWB_")
cmn_files = intersect(eta_files, wb_files)


cross_tissue_counts = data.frame()
all_dge_res = data.frame()

pdf("dge_results.pdf", width=6, height=6)
for(f in cmn_files) {
  cell_t = gsub("_main_dgea.csv","",f)
  eta = read.csv(file.path(eta_dir, paste0("cometETA_",f)), row.names=1)
  wb = read.csv(file.path(wb_dir, paste0("cometWB_",f)), row.names=1)
  col = pals::kelly()

  # For cross_tissue_counts, the Both contains genes that are significantly DE in the same direction. The features that are significantly DE but in different directions are counted towards ETA, which needs to be fixed....
  tmp_df = merge(eta,wb,by="row.names") %>%
      mutate(fc_thresh.x = abs(avg_log2FC.x) > lfcCutoff, fc_thresh.y = abs(avg_log2FC.y) > lfcCutoff) %>%
      filter( !(is.na(avg_log2FC.x) | is.na(avg_log2FC.y)) ) %>%
      mutate(
        sign.x = sign(avg_log2FC.x),
        sign.y = sign(avg_log2FC.y),
        padj_thresh.x = p_val_adj.x < 0.1,
        padj_thresh.y = p_val_adj.y < 0.1,
        new_sig.x = padj_thresh.x & fc_thresh.x,
        new_sig.y = padj_thresh.y & fc_thresh.y,
        eta_group=ifelse(new_sig.x & new_sig.y & sign.x==sign.y, "Both", ifelse(new_sig.x, "ETA", NA ) ),
        wb_group=ifelse(new_sig.x & new_sig.y & sign.x==sign.y, "Both", ifelse(new_sig.y, "WB", NA ) ),
        )
  counts = data.frame(Both = sum(tmp_df$eta_group=="Both", na.rm=T), ETA = sum(tmp_df$eta_group=="ETA", na.rm=T), WB = sum(tmp_df$wb_group=="WB", na.rm=T), row.names=cell_t)
  #counts = table(tmp_df$group)
  #counts = as.data.frame(t(as.matrix(table(tmp_df$group))), row.names=cell_t)
  if(nrow(cross_tissue_counts) == 0) {
    cross_tissue_counts = counts
  } else {
    cross_tissue_counts = bind_rows(cross_tissue_counts, counts )
  }
  tmp_df2 = tmp_df %>%
    mutate( sig=ifelse(new_sig.x & new_sig.y, "Both", ifelse(new_sig.x, "ETA", ifelse(new_sig.y, "WB", "NS") ) ) ) %>%
    select(Row.names,avg_log2FC.x, avg_log2FC.y, sig)

  #tmp_cor = cor.test(tmp_df2$avg_log2FC.x, tmp_df2$avg_log2FC.y, method = "spearman")

  #pval_text = paste0("r: ", round(tmp_cor$estimate,2), "; pval: ", format.pval(tmp_cor$p.value,digits=1) )

  all_dge_res = rbind(all_dge_res, tmp_df2 %>% mutate(ctype = gsub("covidDex-vs-covidNondex_","",cell_t)))
  p = ggplot(tmp_df2, aes(avg_log2FC.y, avg_log2FC.x, color=sig, size=sig)) +
      geom_point() +
      scale_color_manual(values=c("ETA"=col[6],"WB"=col[8],"Both"=col[7],"NS"="grey75")) +
      scale_size_manual(values=c("ETA"=1,"WB"=1,"Both"=1.5,"NS"=0.5)) +
      theme_classic() +
      geom_vline(xintercept=0, linetype=2) +
      geom_hline(yintercept=0, linetype=2) +
      geom_text_repel(data=. %>% filter(sig %in% c("Both","ETA","WB")), aes(label=Row.names), size=2, max.overlaps=100, color="black") +
      xlab("log2FC(Dex/Nondex) WB") + ylab("log2FC(Dex/Nondex) ETA") +
      ggtitle(cell_t) +
      ggpubr::stat_cor(aes(color="black"), method = "spearman", size=4, digits=1)

  print(p)
}
dev.off()

all_dge_res = rename(all_dge_res, log2FC_dex_over_nondex_WB = avg_log2FC.y, log2FC_dex_over_nondex_ETA = avg_log2FC.x )

write.csv(all_dge_res, "Fig3GH_FigS4.csv")

write.csv(cross_tissue_counts, file="cross_tissue_DEG_counts.csv")

