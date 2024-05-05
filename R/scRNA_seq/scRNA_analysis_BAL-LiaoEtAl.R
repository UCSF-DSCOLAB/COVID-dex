

### Clean-up the environment
rm(list = ls())
gc()
#Sys.setenv(R_MAX_VSIZE = "100Gb")

### Load packages
library(tidyverse)
library(Seurat)
library(BiocGenerics)
library(patchwork)
library(hdf5r)
library(future)
library(ggpubr)
library(ggrepel)
library(pheatmap)
library(pals)
library(SingleR)
library(BiocParallel)
library(MAST)
library(SingleCellExperiment)
#library(doParallel)
#library(foreach)
library(ComplexHeatmap)
library(DoubletFinder)


######################################### Functions:
make_breaks <- function(axismax) {
  breaks_gt100 = c( 100, 200, 300, seq(400,1000,200), seq(2000,10000,1000), 20000, 30000 )
  breaks_lt100 = c(0:10, seq(15, 100, 5))
  if(axismax > 100) {
    breaks = breaks_gt100[ breaks_gt100 < axismax ]
  } else {
    breaks = breaks_lt100[ breaks_lt100 < axismax ]
  }
  return(breaks)
}

scatterhist <- function(x, y, sobj, params) {
  x_upper = as.numeric(params[paste0(x,".upper")])
  x_lower = as.numeric(params[paste0(x,".lower")])
  y_upper = as.numeric(params[paste0(y,".upper")])
  y_lower = as.numeric(params[paste0(y,".lower")])
  xy_data = FetchData(sobj, vars=c(x, y))
  filter_cells = xy_data[,1] <= x_upper &
    xy_data[,1] >= x_lower &
    xy_data[,2] <= y_upper &
    xy_data[,2] >= y_lower
  filter_cell_percent = round( sum(filter_cells) / ncol(sobj) ,3) * 100
  xmax = max(sobj@meta.data[,x])
  ymax = max(sobj@meta.data[,y])
  p = ggplot(sobj@meta.data, aes(x=!!sym(x), y=!!sym(y))) +
    #geom_point(size=0.1,alpha=0.1) +
    geom_hex(bins=100) +
    scale_fill_distiller(palette = "RdYlBu") +
    theme_bw() +
    geom_vline(xintercept = x_upper) +
    geom_vline(xintercept = x_lower) +
    geom_hline(yintercept = y_upper) +
    geom_hline(yintercept = y_lower) +
    geom_rect(aes(xmin=x_lower, xmax=x_upper, ymin=y_lower, ymax=y_upper), color="red", alpha=0) +
    geom_text(data = data.frame(x=x_lower, y=y_upper,text=paste0(filter_cell_percent, "%")), aes(x=x,y=y,label=text), hjust=0, vjust=1, color="darkred", size=8) +
    #scale_x_continuous(trans='log2') +
    #  scale_y_continuous(trans='log2')
    #scale_x_sqrt() +
    #scale_y_sqrt() +
    scale_x_continuous(breaks = make_breaks(xmax), trans = scales::sqrt_trans(), sec.axis = dup_axis()) +
    scale_y_continuous(breaks = make_breaks(ymax), trans = scales::sqrt_trans(), sec.axis = dup_axis()) +
    theme(axis.text.x = element_text(angle=90), panel.grid.minor = element_blank())
  return(p)
}


runDoubletFinder = function(sngObj) {
  doubletF_stat = list()
  # If the number of cells is smaller than npcs, adjust the npcs.
  npcs = 35  # Setting number of PCs/dimensions to use for UMAP/Clustering.
  if(length(Cells(sngObj)) <= npcs) npcs = length(Cells(sngObj)) - 1
  print_message("\tThe number of PCs used for DoubletFinder analysis: ", npcs)
  doubletF_stat[["nPCs"]] = npcs
  
  
  print_message("\tRunning Seurat pipeline to make sample(freemuxlet cluster)-specific UMAP/clustering for the use in DoubletFinder")
  sngObj <- SCTransform(sngObj, verbose = FALSE)
  sngObj <- RunPCA(sngObj, npcs=npcs, verbose = FALSE)
  sngObj <- RunUMAP(sngObj, dims = 1:npcs, verbose = FALSE)
  sngObj <- FindNeighbors(sngObj, dims = 1:npcs, verbose = FALSE)
  sngObj <- FindClusters(sngObj, verbose = FALSE)
  
  print_message("\tEstimating pK")
  sweep.res.list <- paramSweep_v3(sngObj, PCs = 1:npcs, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK = as.numeric(as.vector(bcmvn[ bcmvn$BCmetric == max(bcmvn$BCmetric) ,]$pK))
  print_message("\tChosen pK: ", pK)
  doubletF_stat[["pK"]] = pK
  
  # The modeling of homotypic doublets requires celltype annotation. In place of the celltype annotation, I am using the cluster identity as celltype labels (as suggested by the authors of doubletFinder) - since clusters should correspond to different celltypes anyway.
  annotations = sngObj@active.ident
  homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- sngObj@meta.data$ClusteringResults
  effectiveDblRate = sngObj@misc$scStat$effDblRate
  nExp_poi <- ceiling( effectiveDblRate * ncol(sngObj))
  nExp_poi.adj <- ceiling(nExp_poi*(1-homotypic.prop))
  
  print_message("\tEstimated no. of heterotypic doublets:", nExp_poi.adj)
  
  print_message("\tRunning DoubletFinder")
  sngObj <- doubletFinder_v3(sngObj, PCs = 1:npcs, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = TRUE)
  ind = grep("^DF", colnames(sngObj@meta.data))[1]
  dfCalls = sngObj@meta.data[,ind]
  names(dfCalls) = row.names(sngObj@meta.data)
  doubletF_stat[["nDFSng"]] = sum(dfCalls == "Singlet")
  doubletF_stat[["nDFDbl"]] = sum(dfCalls == "Doublet")
  return(list(dfCalls=dfCalls, doubletF_stat=doubletF_stat))
}

find_doublets <- function(sObj, freemuxlet=FALSE, multi_run=FALSE, num_cells_loaded = 15000) {
  effectiveIntraDblRate = 0 # Initialize
  
  # Adapted from Arjun's get_10x_multiplet_rate.R
  df <- read.table("/krummellab/data1/rpatel5/data/10x_multiplet_rate.tsv", header=T)
  y <- lm(multiplet_rate_pct ~ num_cells_loaded, data=df)
  effectiveDblRate = stats::predict(y, new=data.frame(num_cells_loaded=num_cells_loaded)) / 100
  
  names(effectiveDblRate) = NULL
  sObj@misc$scStat$effDblRate = effectiveDblRate
  print_message("Effective intra-sample doublet rate: ", round(effectiveDblRate * 100, 4), "%")
  
  
  # Now that the effective intra-sample doublet rate has been calculated, let's run doubletFinder. If freemuxlet isn't available, the DoubletFinder will be run on all cells as a single run. If freemuxlet is available and if the "multi" mode is off, DoubletFinder is run on all cells from all freemuxlet clusters as a single run, if "multi" is on, DoubletFinder is run on cells from each freemuxlet cluster separately.
  # If multiple samples were pooled in a well and "multi" mode is on, we will first separate cells based on freemuxlet clusters (i.e. make one seurat object per freemuxlet cluster) and run doubletFinder on each cluster. We will retain the doubletFinder calls and discard these seurat sub-objects.
  doubletF_stat = list()
  doubletFinderCalls = c()
  
  sng.cells = rownames(sObj@meta.data)
  
  # Subset sObj to include only singlets
  print_message("\tNumber of cells to use for DoubletFinder: ", length(sng.cells))
  doubletF_stat[["nCells"]] = length(sng.cells)
  sngObj = subset(sObj, cells=sng.cells)
  
  # Run DoubletFinder on single-cells
  dfCalls_stat = runDoubletFinder(sngObj)
  doubletF_stat = dfCalls_stat[[ "doubletF_stat" ]]
  doubletFinderCalls = dfCalls_stat[[ "dfCalls" ]]
  
  sObj@misc$scStat$doubletF_stat = doubletF_stat
  sObj = AddMetaData(sObj, metadata=data.frame(doubletFinderCalls))
  print("Number of doublets")
  print(table(sObj$doubletFinderCalls, useNA="ifany"))
  return(sObj)
}

print_message <- function(...) {
  cat("[", format(Sys.time()), "] ", ..., "\n", sep="")
}

######################################### 



### INPUT::: set absolute paths to the following directories.
dataset = "ET" # ET or WB

data_dir = "natMed_liao_scRNA/scRNA_analysis/data/"
out_dir = "natMed_liao_scRNA/scRNA_analysis/out_blueprint_20220927/"
gsea_dir = "gsea_pathways"


### Set parallel options
#plan("multicore", workers = 8)
#options(future.globals.maxSize = 5000 * 1024^2)


### Set global environment
if(! dir.exists(out_dir))
  dir.create(out_dir)
setwd(out_dir)
today <- as.Date(Sys.time(), format="YYYY-MM-DD")


# If the preprocessed and doublet-removed data already exist, skip to QC.
sobj_list_Rds = list.files(out_dir, pattern = "covid_sobj.list_.*.Rds", full.names = T)
if(length(sobj_list_Rds) == 1) {
  print("Found the sobj list Rds, skipping the preprocessing and doublet finding. Going directly to QC.")
}

# If the processed data already exist, skip to differential freq and DGE analysis.
processed_Rdata = list.files(out_dir, pattern = "covid.combined_.*.Rdata", full.names = T)
if(length(processed_Rdata) == 1) {
  print("Found the processed Rdata, skipping the preprocessing and annotation. Going directly to frequency and DGE analysis.")
}

if( length(sobj_list_Rds) == 0 ) {
  ### Read-in the H5 files
  filenames <- dir(data_dir, full.names = T, pattern = ".*.h5$")
  # Select files for eta_samples
  #timepoint_from_filenames <- as.vector(str_match(filenames, "MVIR1-HS\\d+-D\\d+"))
  
  #filenames <- filenames[ timepoint_from_filenames %in% eta_samples ]
  #libnames <- filenames %>% str_extract("MVIR1-HS[:digit:]+-[^-]+-[^-]+")
  #samplename <- filenames[1] %>% str_extract("MVIR1-HS[:digit:]+-D\\d+")
  samplename <- filenames[1] %>% str_extract("S\\d+")
  covid_sobj <-  Read10X_h5(filenames[1]) %>% CreateSeuratObject(project =samplename, min.features = 100, min.cells = 3 )
  #covid_sobj <- find_doublets(covid_sobj, num_cells_loaded = 15000)
  
  for (file in filenames[2:length(filenames)]) {
    samplename <- file %>% str_extract("S\\d+")
    print(paste("Loading", samplename, "H5 file..."))
    tmp <- Read10X_h5(file) %>% CreateSeuratObject(project = samplename, min.features = 100, min.cells = 3)
    #tmp <- find_doublets(tmp, num_cells_loaded = 15000)
    covid_sobj  <- merge(covid_sobj, y = tmp)
    rm(tmp)
    gc()
  }
  
  saveRDS(covid_sobj, file = paste0("covid_sobj_",today,".Rds"))

  covid_sobj$patient <- covid_sobj$orig.ident
  covid_sobj$dex <- ifelse( covid_sobj$patient %in% paste0("S", 2:5) , TRUE, FALSE )
  
  if(FALSE) {
  # Drop doublets.
  sng.cells <- Cells(covid_sobj)[is.na(covid_sobj$doubletFinderCalls) | covid_sobj$doubletFinderCalls == "Singlet"]
  covid_sobj <- subset(covid_sobj, cells=sng.cells)
  
  ### Add patient data
  covid_sobj$patient <- gsub("-[^-]+$", "", covid_sobj$orig.ident)
  covid_sobj$timepoint <- covid_sobj$orig.ident
  # Add scaled age
  select_patients <- rownames(patient_data) %in% unique(covid_sobj$patient)
  patient_data$age_scaled = NA
  patient_data[select_patients,]$age_scaled <- scale( patient_data$age[ select_patients ] )
  covid_sobj$age_scaled <- patient_data[covid_sobj$patient, ]$age_scaled
  # Add dex status
  covid_sobj$dex <- patient_data_wOverlapStatus[covid_sobj$timepoint, ]$overlap
  # Add sex
  covid_sobj$sex_at_birth <- patient_data[covid_sobj$patient, ]$sex_at_birth
  # Add days_bw_sampling_ventstart
  # If days_bw_sampling_ventstart is missing for any patient, use the value of -100. This is fine since I am not using this in DGE. But if you decide to use it for anything, make sure to find the correct days_bw_sampling_ventstart for HS471-D2.
  patient_data_wOverlapStatus$days_bw_sampling_ventstart[ is.na(patient_data_wOverlapStatus$days_bw_sampling_ventstart) ] = -100
  covid_sobj$days_bw_sampling_ventstart <- patient_data_wOverlapStatus[covid_sobj$timepoint, ]$days_bw_sampling_ventstart
 
  }
 
  if(sum(is.na(covid_sobj@meta.data)) > 0) {
    stop("There are some NA values in metadata. Stopping here.")
  }
  
  # Find number of cells per sample and dex status
  pheatmap(table(covid_sobj$patient, covid_sobj$dex))
  
  ### QC the data
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  covid_sobj <- NormalizeData(covid_sobj) %>% 
    CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
  
  covid_sobj[["percent.mt"]] <- PercentageFeatureSet(covid_sobj, pattern = "^MT-")
  covid_sobj[["percent.ribo"]] <- PercentageFeatureSet(covid_sobj, pattern = "^RP[SL][[:digit:]]|^RPLP[[:digit:]]|^RPSA")
  p = VlnPlot(covid_sobj, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "S.Score", "G2M.Score"), pt.size = 0)
  pdf(paste0(out_dir, "/vlnplots_b4QC_", today, ".pdf"), width = 15, height = 6)
  print(p)
  dev.off()
  
  # Split the object by samples
  covid_sobj.list <- SplitObject(covid_sobj,split.by = "orig.ident")
  saveRDS(covid_sobj.list, file = paste0("covid_sobj.list_",today,".Rds"))
  rm(covid_sobj)
  
  # Make QC plots
  
  # Setting default cutoffs
  params = c("percent.mt.upper"=15, "percent.mt.lower"=0, "percent.ribo.upper"=60, "percent.ribo.lower"=0, "nFeature_RNA.upper"=5000, "nFeature_RNA.lower"=0, "nCount_RNA.upper"=20000, "nCount_RNA.lower"=0)
  
  # Make plots
  pdf(paste0(out_dir, "/qcplots_b4QC_", today, ".pdf"), width = 20, height = 10)
  for( i in 1:length(covid_sobj.list) ) {
    sobj = covid_sobj.list[[i]]
    plot_list = list()
    plot_list[[length(plot_list)+1]] = scatterhist("percent.ribo","percent.mt",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","percent.mt",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nFeature_RNA","percent.mt",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","percent.ribo",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nFeature_RNA","percent.ribo",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","nFeature_RNA",sobj,params)
    p = ggarrange(plotlist=plot_list, ncol=3,nrow=2)
    print(annotate_figure(p, top = text_grob(names(covid_sobj.list)[i], 
                                          color = "red", face = "bold", size = 14)))
    #print(p)
  }
  dev.off()
  stop("Done generating QC plots. Check the QC plots and set appropriate cutoffs for filtering low quality cells.")
}

if( length(processed_Rdata) == 0 ) {
  covid_sobj.list = readRDS(sobj_list_Rds[1])
  ### Filter the low quality cells
  # Print sample names for reference:
  samplenames = names(covid_sobj.list)
  samplenames
  
  # INPUT::: Determine sample-specific cutoffs as needed for the following parameters: c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo").
  cutoffs = data.frame(samplename = samplenames)
  cutoffs$nCount_RNA.lower = c( rep(1000, length(samplenames)) )
  cutoffs$nCount_RNA.upper = c( rep(500000000000, length(samplenames)) )
  #cutoffs$nFeature_RNA.lower = c( rep(150,8), 200, 200, rep(150,4), 200, rep(150,2), 200, rep(150,3), 200 )  # With HS450
  cutoffs$nFeature_RNA.lower = c( rep(200, length(samplenames)) )  # After droping HS450
  cutoffs$nFeature_RNA.upper = c( rep(6000, length(samplenames)) )
  #cutoffs$percent.mt.upper = c( rep(15,8), 25, rep(15,13) ) # With HS450
  cutoffs$percent.mt.upper = c( rep(10, length(samplenames)) ) # After droping HS450
  cutoffs$percent.ribo.upper = c( rep(1000, length(samplenames)) )
  cutoffs$percent.mt.lower = c( rep(0, length(samplenames)) )
  cutoffs$percent.ribo.lower = c( rep(0, length(samplenames)) )
  #cutoffs$S.Score = c()
  #cutoffs$G2M.Score = c()
  covid_sobj.list.copy = covid_sobj.list
  
  # Make plots
  pdf(paste0(out_dir, "/qcplots_b4QC_wSelectCutoff_", today, ".pdf"), width = 20, height = 10)
  for( i in 1:length(covid_sobj.list) ) {
    sobj = covid_sobj.list[[i]]
    plot_list = list()
    params = cutoffs[i,]
    plot_list[[length(plot_list)+1]] = scatterhist("percent.ribo","percent.mt",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","percent.mt",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nFeature_RNA","percent.mt",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","percent.ribo",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nFeature_RNA","percent.ribo",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","nFeature_RNA",sobj,params)
    p = ggarrange(plotlist=plot_list, ncol=3,nrow=2)
    print(annotate_figure(p, top = text_grob(names(covid_sobj.list)[i], 
                                             color = "red", face = "bold", size = 14)))
    #print(p)
  }
  dev.off()
  
  
  pdf(paste0(out_dir, "/qcplots_afterQC_", today, ".pdf"), width = 20, height = 10)
  for( i in 1:length(covid_sobj.list) ) {
    covid_sobj.list[[i]] = subset(covid_sobj.list[[i]], 
                                  nCount_RNA > cutoffs$nCount_RNA.lower[i] & 
                                    nCount_RNA < cutoffs$nCount_RNA.upper[i] & 
                                    nFeature_RNA > cutoffs$nFeature_RNA.lower[i] &
                                    nFeature_RNA < cutoffs$nFeature_RNA.upper[i] &
                                    percent.mt < cutoffs$percent.mt.upper[i] &
                                    percent.ribo < cutoffs$percent.ribo.upper[i]
    )
    p = VlnPlot(covid_sobj.list[[i]], features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "S.Score", "G2M.Score"), pt.size = 0)
    #print(p)
    sobj = covid_sobj.list[[i]]
    plot_list = list()
    plot_list[[length(plot_list)+1]] = scatterhist("percent.ribo","percent.mt",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","percent.mt",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nFeature_RNA","percent.mt",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","percent.ribo",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nFeature_RNA","percent.ribo",sobj,params)
    plot_list[[length(plot_list)+1]] = scatterhist("nCount_RNA","nFeature_RNA",sobj,params)
    p = ggarrange(plotlist=plot_list, ncol=3,nrow=2)
    print(annotate_figure(p, top = text_grob(names(covid_sobj.list)[i], 
                                             color = "red", face = "bold", size = 14)))
  }
  dev.off()
  
  rm(covid_sobj.list.copy)
  gc()
  
  
  ### Perform batch-correction
  # Normalize and identify variable features for each dataset independently
  covid_sobj.list <- lapply(X = covid_sobj.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  })
  
  # Select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = covid_sobj.list)
  
  # Integration
  covid.anchors <- FindIntegrationAnchors(object.list = covid_sobj.list, anchor.features = features)
  saveRDS(covid.anchors, file = paste0("COVIDintegrationanchors_", today , ".Rds"))
  
  #covid.anchors = readRDS("COVIDintegrationanchors_2022-03-20.Rds")
  
  covid.combined <- IntegrateData(anchorset = covid.anchors)
  saveRDS(covid.combined, file = paste0("covid.combined_",today,".Rds"))
  rm(covid.anchors)
  rm(covid_sobj.list)
  gc()
  #covid.combined = readRDS("covid.combined_2022-03-20.Rds")
  
  
  
  ### Analysis of the integrated data
  covid.combined <- ScaleData(covid.combined, vars.to.regress = c("orig.ident","nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "S.Score", "G2M.Score"))
  covid.combined <- RunPCA(covid.combined, npcs = 30)
  covid.combined <- RunUMAP(covid.combined, dims = 1:30)
  covid.combined <- FindNeighbors(covid.combined, reduction = "pca", dims = 1:30)
  covid.combined <- FindClusters(covid.combined, resolution = 1)
  DimPlot(covid.combined, pt.size = 0.1, group.by = "seurat_clusters", label = T)
  DimPlot(covid.combined, pt.size = 0.1, group.by = "orig.ident", label = T) + scale_color_manual(values=kelly(22))
  FeaturePlot(covid.combined, features=c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "S.Score", "G2M.Score"))
  
  saveRDS(covid.combined, file = paste0("covid.combined_",today,".Rds"))
  
  
  ### Automated cell-type annotation using SingleR
  MulticoreParam(6)
  
  # Download/Load the reference dataset
  #hpca.se <- celldex::HumanPrimaryCellAtlasData()
  blueprint <- celldex::BlueprintEncodeData()
  ref <- blueprint
  
  # Run SingleR
  DefaultAssay(covid.combined) <- "RNA"
  somecells <- GetAssayData(object = covid.combined, slot = "data")
  
  cell.annot <- SingleR(test = somecells, ref = ref,assay.type.test = 1,
                        labels = ref$label.fine,BPPARAM = MulticoreParam(6))
  
  # Assign annotation labels to the Seurat object; Remove the cells that are called NA in pruned.labels
  covid.combined$singleR.label <- cell.annot$pruned.labels
  covid.combined$major.label <- cell.annot$pruned.labels %>% str_remove(":.+")
  covid.combined <- subset(covid.combined,  cells=Cells(covid.combined)[!is.na(covid.combined$major.label)] )

  # Plot SingleR scores for cell types with >500 cells.
  #plotScoreHeatmap(cell.annot, cells.use=Cells(covid.combined)[ covid.combined$major.label %in% names(which(table(covid.combined$major.label) > 500))] )
  saveRDS(cell.annot, file = paste0("cell.annot_",today,".rds"))
  rm(somecells)
  rm(cell.annot)
  rm(ref)
  
  
  if(FALSE) { # Assessing SingleR scores.
    ct = "Monocyte"
    top_n = 12
    ct_index = which(covid.combined$major.label == ct)
    ct_index = which(covid.combined$singleR.label == ct)
    ct_score = apply(cell.annot$scores[ct_index,], 1, function(x) head(sort(x, decreasing=T), top_n) ) %>% apply(2,scale) %>% t() %>%
      data.frame() %>% 
      setNames(1:top_n) %>% 
      mutate(idx=1:nrow(.)) %>% 
      pivot_longer(names_to="rank", cols=-idx) %>% 
      mutate(rank = factor(rank, levels=as.character(sort(unique(as.numeric(rank))))))
    ct_names = apply(cell.annot$scores[ct_index,], 1, function(x) names(head(sort(x, decreasing=T), top_n)) ) %>% t() %>% 
      data.frame() %>% 
      setNames(1:top_n) %>% 
      mutate(idx=1:nrow(.)) %>% 
      pivot_longer(names_to="rank", cols=-idx) %>%
      group_by(rank, value) %>%
      summarise(freq = n()) %>% 
      ungroup() %>%
      mutate(rank = factor(rank, levels=as.character(sort(unique(as.numeric(rank))))))
    
    
    ggplot(ct_score, aes(rank, value, group=idx)) + geom_line() + theme(legend.position = "none")
    ggplot(ct_names, aes(rank,freq,fill=value)) + geom_bar(stat="identity", position="fill") + scale_fill_manual(values=rep(pals::kelly(),4))
  }
  
  
  save(covid.combined, file = paste0("covid.combined_",today,".Rdata"))
  # Some plots
  pdf("umaps_by_various_things_v1.pdf", width=15, height=12)
  print(DimPlot(covid.combined, pt.size = 0.1, group.by = "singleR.label", label = T) + NoLegend())
  print(DimPlot(covid.combined, pt.size = 0.1, group.by = "major.label", label = T))
  print(DimPlot(covid.combined, pt.size = 0.1, group.by = "dex", label = T))
  print(DimPlot(covid.combined, pt.size = 0.1, label = F, group.by = "major.label", split.by = "major.label", ncol=4))
  dev.off()

  Idents(covid.combined) = covid.combined$major.label
  marker_rna = FindAllMarkers( covid.combined, assay = "RNA", test.use = "MAST", only.pos = T,  max.cells.per.ident = 1000, random.seed = 1234, latent.vars = "patient")
  write.csv(marker_rna, file = paste0("marker_rna_",today,".csv"))

  top_5_genes = marker_rna %>% arrange(p_val_adj, -avg_log2FC) %>% group_by(cluster) %>% dplyr::filter(row_number() %in% 1:10) %>% arrange(cluster) %>% pull(gene) %>% unique()

  pdf("dotplot.pdf", width=30, height=10)
  p = DotPlot(covid.combined, features = unique(c(top_5_genes)) , cols = c("blue","red"), assay = "RNA") + theme(axis.text.x = element_text(angle=90)) + scale_color_distiller(palette = "RdBu")
  print(p)
  dev.off()
  
}


if(length(processed_Rdata) == 1)
  load(processed_Rdata, verbose = TRUE)


# Clean up the cell types.
celltype_change = read.csv("../../../cross_tissue/celltype_label_conversion.csv")
covid.combined$major.label.uncleaned = covid.combined$major.label
covid.combined$major.label = ifelse( covid.combined$major.label %in% names(which(table(covid.combined$major.label) > 100)), covid.combined$major.label, "other")
covid.combined$major.label = celltype_change$major.label2[ match(covid.combined$major.label, celltype_change$major.label) ]
if(sum(is.na(covid.combined$major.label)) > 0) {
  stop("NAs were introduced in the cleaned up major.labels. Fix the issue and rerun.")
}
save(covid.combined, file = paste0("covid.combined_",today,".Rdata"))

