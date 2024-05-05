library(CellChat)
library(patchwork)
library(tidyverse)
library(pheatmap)
library(ggpubr)

rdss = list.files(".", pattern="cellchat_WB.*.rds", full.names=T)
objnames = gsub(".rds","",gsub("cellchat_","",rdss)) %>% setNames(rdss)

object.list = list()
for(f in rdss) {
  ccsub <- readRDS(f)
  ccsub <- computeCommunProbPathway(ccsub)
  ccsub <- netAnalysis_computeCentrality(ccsub, slot.name = "netP")
  object.list[[objnames[f]]] = ccsub
}

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix=T)
cellchat_wb <- cellchat
object.list_wb <- object.list


rdss = list.files(".", pattern="cellchat_SINHA.*.rds", full.names=T)
objnames = gsub(".rds","",gsub("cellchat_","",rdss)) %>% setNames(rdss)

object.list = list()
for(f in rdss) {
  ccsub <- readRDS(f)
  ccsub <- computeCommunProbPathway(ccsub)
  ccsub <- netAnalysis_computeCentrality(ccsub, slot.name = "netP")
  object.list[[objnames[f]]] = ccsub
}

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix=T)
cellchat_sinha <- cellchat
object.list_sinha <- object.list


rdss = list.files(".", pattern="cellchat_ETA.*.rds", full.names=T)
objnames = gsub(".rds","",gsub("cellchat_","",rdss)) %>% setNames(rdss)

object.list = list()
for(f in rdss) {
  ccsub <- readRDS(f)
  ccsub <- computeCommunProbPathway(ccsub)
  ccsub <- netAnalysis_computeCentrality(ccsub, slot.name = "netP")
  object.list[[objnames[f]]] = ccsub
}

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix=T)
cellchat_eta <- cellchat
object.list_eta <- object.list



rdss = list.files(".", pattern="cellchat_BAL.*.rds", full.names=T)
objnames = gsub(".rds","",gsub("cellchat_","",rdss)) %>% setNames(rdss)

object.list = list()
for(f in rdss) {
  ccsub <- readRDS(f)
  ccsub <- computeCommunProbPathway(ccsub)
  ccsub <- netAnalysis_computeCentrality(ccsub, slot.name = "netP")
  object.list[[objnames[f]]] = ccsub
}

cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix=T)
cellchat_bal <- cellchat
object.list_bal <- object.list


celltypes = unique(c(
  unique(cellchat_eta@meta$major.label.o),
  unique(cellchat_wb@meta$major.label.o),
  unique(cellchat_sinha@meta$major.label.o),
  unique(cellchat_bal@meta$major.label.o)
  ))

celltypes_cols = pals::kelly(12)[3:12] %>% setNames(celltypes)

celltypes_cols2 = c()
for(s in c("WB","WBH","ETA","SINHA","BAL","BALH")) {
  celltypes_cols2 = c(
    celltypes_cols2,
    pals::kelly(11)[-1] %>% setNames(paste0(s,":",celltypes))
    )
}


path = c("COLLAGEN","ANNEXIN","ICAM","CD99","ITGB2","MIF","SELPLG","MHC-II")
path_wb = c("COLLAGEN","ANNEXIN","ICAM","CD99","ITGB2","MIF","SELPLG","MHC-II")
path_eta = c("SELPLG","MHC-II","GALECTIN","ANNEXIN","ITGB2","CD86","ICAM","IL16","CXCL","LCK")






library(gridExtra)
# Significant pathway heatmaps
all_wb = mergeCellChat(append(object.list_wb, object.list_sinha), add.names=names(append(object.list_wb, object.list_sinha)), cell.prefix=T)

make_infoFlow_heatmap = function(myobj, show_rownames=TRUE, silent=FALSE, return_data=FALSE) {
  dc = length(levels(myobj@meta$datasets))
  # Select pathways that are significant (padj < 0.1 and log2FC > 1) in at least one comparison
  paths = c(); 
  for(i in 1:(dc-1)){ 
    for(j in (i+1):dc) { 
      path_tmp = rankNet(myobj, mode = "comparison", stacked = F, do.stat = TRUE, comparison=c(i,j))$data %>% 
        mutate(cont=contribution.scaled, padjs = p.adjust(pvalues, "BH")) %>% 
        group_by(name, padjs) %>% 
        summarise(cont_lfc = abs(log2(cont[1]+0.1)-log2(cont[2]+0.1)) ) %>% 
        dplyr::filter(padjs < 0.1 & cont_lfc > 1) %>% 
        pull(name) %>% as.vector()
      paths=unique(c(paths,path_tmp)) 
    } 
  }
  plot = rankNet(myobj, mode = "comparison", stacked = F, do.stat = TRUE, comparison=1:dc)
  tmp_df = plot$data %>% 
    filter( name %in% paths) %>% 
    select(name, contribution.scaled, group) %>% 
    pivot_wider(names_from=group, values_from=contribution.scaled) %>% 
    column_to_rownames("name")
  p = log2(tmp_df+0.1) %>% 
    pheatmap(border_color=NA, show_rownames=show_rownames, silent=silent)
  if(return_data) {
    return(log2(tmp_df+0.1))
  }
  return(p)
}
pdf("significant_pathway_heatmaps.pdf", width=6, height=6)
p1 = make_infoFlow_heatmap(all_wb, T, T)
p2 = make_infoFlow_heatmap(cellchat_eta, T, T)
grid.arrange(grobs=list(p1[[4]],p2[[4]]), ncol=2, widths=c(2,1.5))
data_p1 = make_infoFlow_heatmap(all_wb, T, T, T)
data_p2 = make_infoFlow_heatmap(cellchat_eta, T, T, T)
dev.off()

write.csv(data_p2, "../../Fig5A.csv")
write.csv(data_p1, "../../Fig5D.csv")




# Cell communication plots.
make_cell_communication_plots <- function(object.list, path) {
  for(i in 1:length(path)) {
    pathways.show <- path[i] 
    print(pathways.show)
    object.list.tmp = object.list[ unlist(lapply(object.list, function(x) any(x@netP$pathways == pathways.show) )) ]
    celltypes = unique(unlist(lapply(object.list.tmp, function(x) {unique(x@meta$major.label.o)} )))
    weight.max <- getMaxWeight(object.list.tmp, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
    par(mfrow = c(2,2), xpd=TRUE)
    sig.features = c()
    for (j in 1:length(object.list.tmp)) {
      tmp_obj = object.list.tmp[[j]]
      tmp_obj = updateClusterLabels(tmp_obj, new.cluster.name=gsub("^[^:]+:","",levels(tmp_obj@idents)), new.order=celltypes)
      netVisual_aggregate(tmp_obj, color.use=celltypes_cols[celltypes], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list.tmp)[j]))
      sig.feat = extractEnrichedLR(object.list.tmp[[j]], signaling = pathways.show, geneLR.return=T)$geneLR
      print(sig.feat)
      print(length(sig.feat))
      if(length(sig.features) == 0) {
        sig.features = sig.feat
      } else {
        sig.features = unique(c(sig.features, sig.feat))
      }
    }

    print(sig.features)
    print(plotGeneExpression(
      object = mergeCellChat(object.list, add.names = names(object.list), cell.prefix=T), 
      features = sig.features, signaling = pathways.show, split.by = "datasets", colors.ggplot = T, group.by = "major.label.o"))
  }
}


pdf("select_cellchat_pathway_figs_v1_wb.pdf", height = 12, width = 12)
make_cell_communication_plots(object.list_wb, path_wb)
dev.off()
pdf("select_cellchat_pathway_figs_v1_eta.pdf", height = 12, width = 12)
make_cell_communication_plots(object.list_eta, path_eta)
dev.off()




make_information_flow_plot <- function(cellchat, pathw) {

  if( length(levels(cellchat@meta$datasets)) == 3 ) {

    t12 = rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison=c(1:2))
    t23 = rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison=c(2,3))
    t13 = rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison=c(1,3))
    t3 = rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison=c(1:3))

    # Collect pairwise comparison results
    df1 = rbind( 
      t13$data[ t13$data$name %in% pathw, ] %>% mutate(comp="NonDex_Healthy"),
      t12$data[ t12$data$name %in% pathw, ] %>% mutate(comp="NonDex_Dex"),
      t23$data[ t23$data$name %in% pathw, ] %>% mutate(comp="Dex_Healthy")
      )

    df1$pvalues = p.adjust(df1$pvalues, "BH")

    # Get max information flow for each pathway
    y.position = t3$data %>% group_by(name) %>% summarise(y.position = max(contribution.scaled))

    # Gather pvalues for each pairwise comparison, add significance stars for pvalues and add y.position for the star labels.
    df2 =
      df1 %>% dplyr::select(name, comp, pvalues) %>% unique() %>%
      inner_join(y.position, by="name") %>%
      mutate(group1 = gsub("_.+$","",comp), 
        group2 = gsub("^.+_","",comp), 
        p = pvalues, 
        p.signif = case_when(p < 0.00001 ~ "****", p < 0.0001 ~ "***", p < 0.001 ~ "**", p < 0.05 ~ "*", TRUE ~ ""),
        y.position = case_when(comp == "NonDex_Healthy" ~ (y.position + 7), comp == "NonDex_Dex" ~ (y.position + 3), comp == "Dex_Healthy" ~ (y.position + 5))
        )

    # Gather information flow for all groups
    df3 = t3$data %>% filter(name %in% pathw) %>% 
      mutate(group = case_when( grepl("FALSE",group) ~ "NonDex", grepl("TRUE",group) ~ "Dex", TRUE ~ "Healthy" ),
      group = factor(group, levels=c("NonDex","Dex","Healthy")),
      name = factor(name, levels=pathw) )

    # Make plot
    p = df3 %>% ggplot(aes(group, contribution.scaled, fill=group)) + geom_bar(stat="identity") + facet_wrap(~name, nrow=1) + theme_classic() + theme(axis.text.x=element_blank())
    p = p + stat_pvalue_manual(df2 %>% filter(p.signif != "") %>% mutate(group="Dex"), label="p.signif", tip.length=0.01)
    return(p)
  }

  if( length(levels(cellchat@meta$datasets)) == 2 ) {
    t12 = rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, comparison=c(1:2))

    # Collect pairwise comparison results
    df1 = t12$data[ t12$data$name %in% pathw, ] %>% mutate(comp="NonDex_Dex")    

    df1$pvalues = p.adjust(df1$pvalues, "BH")

    # Get max information flow for each pathway
    y.position = t12$data %>% group_by(name) %>% summarise(y.position = max(contribution.scaled))

    # Gather pvalues for each pairwise comparison, add significance stars for pvalues and add y.position for the star labels.
    df2 =
      df1 %>% dplyr::select(name, comp, pvalues) %>% unique() %>%
      inner_join(y.position, by="name") %>%
      mutate(group1 = gsub("_.+$","",comp), 
        group2 = gsub("^.+_","",comp), 
        p = pvalues, 
        p.signif = case_when(p < 0.00001 ~ "****", p < 0.0001 ~ "***", p < 0.001 ~ "**", p < 0.05 ~ "*", TRUE ~ ""),
        y.position = case_when(comp == "NonDex_Dex" ~ (y.position + 1))
        )

    # Gather information flow for all groups
    df3 = t12$data %>% filter(name %in% pathw) %>% 
      mutate(group = case_when( grepl("FALSE",group) ~ "NonDex", grepl("TRUE",group) ~ "Dex", TRUE ~ "Healthy" ),
      group = factor(group, levels=c("NonDex","Dex","Healthy")),
      name = factor(name, levels=pathw) )

    # Make plot
    p = df3 %>% ggplot(aes(group, contribution.scaled, fill=group)) + geom_bar(stat="identity") + facet_wrap(~name, nrow=1) + theme_classic() + theme(axis.text.x=element_blank())
    p = p + stat_pvalue_manual(df2 %>% filter(p.signif != "") %>% mutate(group="Dex"), label="p.signif", tip.length=0.01) + ylab("Information flow")
    return(p)
  }
}


pdf("select_cellchat_pathway_figs_v2_wb.pdf", height = 5, width = 7)
p = make_information_flow_plot(cellchat_wb, path_wb)
dev.off()
p$data %>% dplyr::select(name, contribution.scaled, group) %>% write.csv("../../Fig5E.csv")
#pdf("select_cellchat_pathway_figs_v2_eta.pdf", height = 5, width = 7)
#make_information_flow_plot(cellchat_eta, path_eta)
#dev.off()


