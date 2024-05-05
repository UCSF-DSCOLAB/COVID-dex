library(CellChat)
library(patchwork)


eta_rdata = "../../ETA_scRNA/scRNA_analysis/out_blueprint_20220927/covid.combined_2022-09-28.Rdata" 
wb_rdata = "../../../WB_scRNA/scRNA_analysis/out_blueprint_20220927/covid.combined_2022-10-04.Rdata"
wb_h_rdata = "../../../WB_scRNA/scRNA_analysis/healthy_out_blueprint_20220927/covid.combined_2022-10-12.Rdata"
sinha_rdata = "../../../natMed_scRNA/scRNA_analysis/out_blueprint_20220927/covid.combined_2022-10-04.Rdata" 
bal_rdata = "../../../natMed_liao_scRNA/scRNA_analysis/out_blueprint_20220927/covid.combined_2022-10-04.Rdata"
bal_h_rdata = "../../../natMed_liao_scRNA/scRNA_analysis/healthy_out_blueprint_20220927/covid.combined_2022-10-12.Rdata"


load(eta_rdata, verbose=T) 
eta = covid.combined
load(wb_rdata, verbose=T) 
wb = covid.combined
load(wb_h_rdata, verbose=T)
wb_h = covid.combined
load(sinha_rdata, verbose=T)
sinha = covid.combined
load(bal_rdata, verbose=T)
bal = covid.combined
load(bal_h_rdata, verbose=T)
bal_h = covid.combined


wb$major.label.o = wb$major.label
eta$major.label.o = eta$major.label
wb_h$major.label.o = wb_h$major.label
sinha$major.label.o = sinha$major.label
bal_h$major.label.o = bal_h$major.label
bal$major.label.o = bal$major.label

wb$major.label = paste0("WB:", wb$major.label)
eta$major.label = paste0("ETA:", eta$major.label)
wb_h$major.label = paste0("WBH:", wb_h$major.label)
sinha$major.label = paste0("SINHA:", sinha$major.label)
bal$major.label = paste0("BAL:", bal$major.label)
bal_h$major.label = paste0("BALH:", bal_h$major.label)

gt100_wb = names(which(table(wb$major.label) > 100))
wb = subset(wb, cells=Cells(wb)[wb$major.label %in% gt100_wb] )
gt100_eta = names(which(table(eta$major.label) > 100))
eta = subset(eta, cells=Cells(eta)[eta$major.label %in% gt100_eta] )
gt100_wb_h = names(which(table(wb_h$major.label) > 100))
wb_h = subset(wb_h, cells=Cells(wb_h)[wb_h$major.label %in% gt100_wb_h] )
gt100_sinha = names(which(table(sinha$major.label) > 100))
sinha = subset(sinha, cells=Cells(sinha)[sinha$major.label %in% gt100_sinha] )
gt100_bal = names(which(table(bal$major.label) >= 20))
bal = subset(bal, cells=Cells(bal)[bal$major.label %in% gt100_bal] )
gt100_bal_h = names(which(table(bal_h$major.label) >= 20))
bal_h = subset(bal_h, cells=Cells(bal_h)[bal_h$major.label %in% gt100_bal_h] )


#seurat_obj = merge(wb, eta)
wb_dex_nondex_obj_list = Seurat::SplitObject(wb, split.by="dex")
names(wb_dex_nondex_obj_list) = paste0("WB:", names(wb_dex_nondex_obj_list))
eta_dex_nondex_obj_list = Seurat::SplitObject(eta, split.by="dex")
names(eta_dex_nondex_obj_list) = paste0("ETA:", names(eta_dex_nondex_obj_list))
sinha_dex_nondex_obj_list = Seurat::SplitObject(sinha, split.by="dex")
names(sinha_dex_nondex_obj_list) = paste0("SINHA:", names(sinha_dex_nondex_obj_list))
bal_dex_nondex_obj_list = Seurat::SplitObject(bal, split.by="dex")
names(bal_dex_nondex_obj_list) = paste0("BAL:", names(bal_dex_nondex_obj_list))
healthy_list = list( "WBH" = wb_h, "BALH" = bal_h )
all_objs = append( 
		append(
			append(wb_dex_nondex_obj_list, eta_dex_nondex_obj_list),
			append(sinha_dex_nondex_obj_list, bal_dex_nondex_obj_list)
		),
	    healthy_list)



## Select the common cell types separately for each tissue
objs_wb = all_objs[ grep("WB|SINHA", names(all_objs), value=T) ]
objs_eta = all_objs[ grep("ETA|BAL", names(all_objs), value=T) ]
common_cell_wb = c()
for(obj in objs_wb) {
  if(length(common_cell_wb) == 0) {
    common_cell_wb = unique(obj$major.label.o)
  }
  else {
    common_cell_wb = intersect( common_cell_wb,  unique(obj$major.label.o))
  }
}

common_cell_eta = c()
for(obj in objs_eta) {
  if(length(common_cell_eta) == 0) {
    common_cell_eta = unique(obj$major.label.o)
  }
  else {
    common_cell_eta = intersect( common_cell_eta,  unique(obj$major.label.o))
  }
}

# Remove "other"
common_cell_wb = common_cell_wb[ common_cell_wb != "other" ]
common_cell_eta = common_cell_eta[ common_cell_eta != "other" ]

print(common_cell_wb)
print(common_cell_eta)


for(n in names(objs_wb)) {
  obj = objs_wb[[n]]
  objs_wb[[n]] = subset(obj, cells=Cells(obj)[obj$major.label.o %in% common_cell_wb] )
}

for(n in names(objs_eta)) {
  obj = objs_eta[[n]]
  objs_eta[[n]] = subset(obj, cells=Cells(obj)[obj$major.label.o %in% common_cell_eta] )
}

all_objs = append(objs_wb, objs_eta)





#set.seed(1234)
#seurat_obj = merge(subset(wb, cells=sample(Cells(wb), 10000)), subset(eta, cells=sample(Cells(eta), 10000)))
for(o in names(all_objs)) {
  seurat_obj_sub <- all_objs[[o]]
  data.input <- GetAssayData(seurat_obj_sub, assay = "RNA", slot = "data") # normalized data matrix
  Idents(seurat_obj_sub) <- seurat_obj_sub$major.label
  labels <- Idents(seurat_obj_sub)
  meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
  meta <- cbind( meta, seurat_obj_sub@meta.data)

  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")


  CellChatDB <- CellChatDB.human
  showDatabaseCategory(CellChatDB)

  cellchat@DB <- CellChatDB
  # subset the expression data of signaling genes for saving computation cost
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  future::plan("multiprocess", workers = 4) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # project gene expression data onto PPI network (optional)
  #cellchat <- projectData(cellchat, PPI.human)


  options(future.globals.maxSize= 1091289600)
  cellchat <- computeCommunProb(cellchat)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  cellchat <- filterCommunication(cellchat, min.cells = 10)


  cellchat <- aggregateNet(cellchat)

  saveRDS(cellchat, paste0("cellchat_",o,".rds"))
}


