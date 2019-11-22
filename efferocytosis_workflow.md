## g0_v_g1_ctrl_n_ifn_n_ifn_tun

library(Seurat)
library(ggplot2)
library(cowplot)
library(stringr)
library(xlsx)

###############################################################################
## Setup Seruat g0_ctrl, g1_ctrl, g0_ifn & g1_ifn
merkd.ctrl <- Read10X(data.dir = "../MerKD_Ctrl")
merkd_ctrl <- CreateSeuratObject(counts = merkd.ctrl, project = "MerKD_Ctrl")
merkd_ctrl$stim <- "MerKD_Ctrl"
merkd_ctrl <- SCTransform(merkd_ctrl)

merkd.2hr <- Read10X(data.dir = "../MerKD_2_hr")
merkd_2hr <- CreateSeuratObject(counts = merkd.2hr, project = "MerKD_2hr")
merkd_2hr$stim <- "MerKD_2hr"
merkd_2hr <- SCTransform(merkd_2hr)

merkd.6hr <- Read10X(data.dir = "../MerKD_6_hr")
merkd_6hr <- CreateSeuratObject(counts = merkd.6hr, project = "MerKD_6hr")
merkd_6hr$stim <- "MerKD_6hr"
merkd_6hr <- SCTransform(merkd_6hr)

wt.ctrl <- Read10X(data.dir = "../WT_Ctrl")
wt_ctrl <- CreateSeuratObject(counts = wt.ctrl, project = "WT_Ctrl")
wt_ctrl$stim <- "WT_Ctrl"
wt_ctrl <- SCTransform(wt_ctrl)

wt.2hr <- Read10X(data.dir = "../WT_2_hr")
wt_2hr <- CreateSeuratObject(counts = wt.2hr, project = "WT_2hr")
wt_2hr$stim <- "WT_2hr"
wt_2hr <- SCTransform(wt_2hr)

wt.6hr <- Read10X(data.dir = "../WT_6_hr")
wt_6hr <- CreateSeuratObject(counts = wt.6hr, project = "WT_6hr")
wt_6hr$stim <- "WT_6hr"
wt_6hr <- SCTransform(wt_6hr)

## Integrate Conditions
all.features <- SelectIntegrationFeatures(object.list = list(merkd_ctrl, merkd_2hr, merkd_6hr, wt_ctrl, wt_2hr, wt_6hr),
                                          nfeatures = 3000)
all.check <- PrepSCTIntegration(object.list = list(merkd_ctrl, merkd_2hr, merkd_6hr, wt_ctrl, wt_2hr, wt_6hr),
                                anchor.features = all.features)
all.anchors <- FindIntegrationAnchors(object.list = all.check, normalization.method = "SCT", anchor.features = all.features)
all.combined <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT")

## Run PCA and UMAP
all.combined <- RunPCA(object = all.combined)
all.combined <- RunUMAP(object = all.combined, dims = 1:30)

## Visualization
plots <- DimPlot(all.combined, group.by = c(), combine = FALSE)
plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 3,
                                                                                                              byrow = TRUE,
                                                                                                              override.aes = list(size = 3))))
CombinePlots(plots)

## Cell Cluster Specific Genes
all.combined <- FindNeighbors(all.combined, dims = 1:10)
all.combined <- FindClusters(all.combined, resolution = 0.4)
DimPlot(all.combined, group.by="seurat_clusters", label=TRUE)
DimPlot(all.combined, group.by='singler.Immgen', label=FALSE)

## Create FindConservedMarkers Function
DefaultAssay(object = all.combined) <- "RNA"
clus = 0
while (clus < length(levels(all.combined@meta.data$integrated_snn_res.0.4)))
{
    print(paste("Cluster", clus, sep=" "))
    markers <- FindConservedMarkers(object = all.combined, ident.1=clus, grouping.var="stim", verbose=FALSE)
    if (clus == 0) {
        write.xlsx(markers, file="conserved_markers.xlsx", sheetName = paste("Cluster", clus, sep="_"), append=FALSE)
    }
    else {
        write.xlsx(markers, file="conserved_markers.xlsx", sheetName = paste("Cluster", clus, sep="_"), append=TRUE)
    }
    clus = clus + 1
}

DimPlot(combined, group.by="singler", label=FALSE, cols=c('Macrophages (MF.II-480HI)'='#990000', 'B cells (B1a)'='#000099',
                                  'Macrophages (MF.103-11B+24-)'='#FF3333', 'B cells (B1A)'='#33FFFF',
                                  'T cells (T.CD4TESTCJ)'='#999900'))


#############################################################################
## SingleR Analysis
#############################################################################
library(Seurat)
library(SingleR)
library(scRNAseq)
library(scater)

Idents(all.combined) <- all.combined@meta.data$seurat_clusters

immgen.ref <- ImmGenData()
test <- as.SingleCellExperiment(all.combined)

common <- intersect(rownames(test), rownames(immgen.ref))
immgen.ref <- immgen.ref[common,]
test <- test[common,]
test <- logNormCounts(test)

clusters <- all.combined@active.ident

pred <- SingleR(test=test, ref=immgen.ref, labels=immgen.ref$label.fine, method="cluster", clusters=clusters)

## Add Seurat New Cluster IDs
new.cluster.ids <- pred$pruned.labels

names(new.cluster.ids) <- levels(all.combined)

all.combined <- RenameIdents(all.combined, new.cluster.ids)
DimPlot(all.combined, reduction="umap", label=TRUE)

combined <- all.combined
combined@meta.data$singler <- Idents(all.combined)

#############################################################################
## Monocle 3 Analysis
#############################################################################
library(Seurat)
library(monocle3)
library(reticulate)
library(hdf5r)
library(dplyr)

## Set Identities to cluster.id
Idents(combined) <- combined@meta.data$singler

## Extract counts, phenotype data, and feature data from SeuratObject
data <- as(as.matrix(combined@assays$integrated@data), 'sparseMatrix')
pd <- combined@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

## Construct Monocle CDS
cds <- new_cell_data_set(data,
                         cell_metadata = pd,
                         gene_metadata = fData)

cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds)

## Use cell.cluster Labels w/ Monocle3 UMAP
plot_cells(cds, color_cells_by="singler")

plot_cells(cds, color_cells_by="seurat_clusters", group_label_size=4)

cds <- cluster_cells(cds)

cds <- learn_graph(cds)

plot <- plot_cells(cds, color_cells_by = "singler", group_label_size=4,
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE, label_branch_points=FALSE) +
    scale_color_manual(values = c('Macrophages (MF.II-480HI)'='#990000', 'B cells (B1a)'='#000099',
                                  'Macrophages (MF.103-11B+24-)'='#FF3333', 'B cells (B1A)'='#33FFFF',
                                  'T cells (T.CD4TESTCJ)'='#999900'))
ggsave("monocle3_plot_singler_labeled.png", plot, width=40, height=40, units="cm")


cds <- order_cells(cds)

plot_cells(cds,
           color_cells_by = "pseudotime",
           graph_label_size=1.5)

plot_cells(cds,
           color_cells_by = "stim",
           graph_label_size=1.5,
           label_groups_by_cluster=FALSE,
           label_cell_groups=FALSE)

g1ifntun <- cds[,colData(cds)$stim == "G1_IFN_TUN"]
plot_cells(g1ifntun, color_cells_by="cluster.id",
           label_groups_by_cluster=FALSE, cell_size=1,
           group_label_size=4)

#############################################################################
## SWNE Plot w/ Seurat v3.0
library(Seurat)
library(swne)

## Extract integrated expression matrix
aligned.counts <- as.matrix(GetAssayData(combined, assay = "integrated"))

dim(aligned.counts)

hist(as.matrix(aligned.counts))

## Reconstruct matrix to remove negative values
aligned.counts <- t(apply(aligned.counts, 1, function(x) (x - min(x))/(max(x) - min(x))))

### Run NMF
k <- 20
n.cores <- 14
nmf.res <- RunNMF(aligned.counts, k = k, alpha = 0, init = "ica", n.cores = n.cores)

pc.emb <- Embeddings(combined, "pca")
snn <- CalcSNN(t(pc.emb), k = 20, prune.SNN = 0)

swne.embedding <- EmbedSWNE(nmf.res$H, snn, alpha.exp = 1.5, snn.exp = 1, n_pull = 3)

norm.counts <- ExtractNormCounts(combined, rescale = T, rescale.method = "log", batch = NULL)

nmf.res$W <- ProjectFeatures(norm.counts, nmf.res$H, n.cores = n.cores)

gene.factor.df <- SummarizeAssocFeatures(nmf.res$W, features.return = 1)
genes.embed <- unique(gene.factor.df$feature)

swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, genes.embed, n_pull = 3)
swne.embedding$H.coords$name <- ""

clusters <- combined@active.ident

swneplot <- PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = clusters, do.label = F,
                     label.size = 3.5, pt.size = 1.25, show.legend = F, seed = 42,
                     colors.use = c('Macrophages (MF.II-480HI)'='#990000', 'B cells (B1a)'='#000099',
                                    'Macrophages (MF.103-11B+24-)'='#FF3333', 'B cells (B1A)'='#33FFFF',
                                    'T cells (T.CD4TESTCJ)'='#999900'))
ggsave("swne_plot_colored.png", swneplot, width=40, height=40, units="cm")


### Cell Specific Markers
swne.embedding$feature.coords$name <- ""
genes.use <- "Cd226"
gene.expr <- norm.counts[genes.use,]
FeaturePlotSWNE(swne.embedding, gene.expr, genes.use, alpha.plot=0.4, label.size=3.5, pt.size=2, color.palette="Reds")

######################## Dot Plots ###############################################
Idents(combined) <- "stim"

DotPlot(combined, features = c("Cd3e", "Cd4"), dot.scale=30)  ## T Cell Markers

DotPlot(combined, features = c("Cd19", "Cd79a", "Cd80", "Ms4a1"), dot.scale=30)  ## B Cell Markers

DotPlot(combined, features = c("Itgam", "Itgax", "Adgre1", "H2-Ab1", "Mertk", "Icam2", "Il1b", "Ccr2", "Cd226"),
        dot.scale=30)  ## Macrophages Markers

DotPlot(combined, features = c("Cd3e", "Cd4", "Cd8", "Cd19", "Cd79a", "Cd80", "Ms4a1", "Itgam", "Itgax", "Adgre1",
                               "H2-Ab1", "Mertk", "Icam2", "Il1b", "Ccr2", "Cd226"), dot.scale=20)


########################################################################################################
## VlnPlots
Idents(combined) <- "singler"

VlnPlot(object = combined, features = c("Cd3e", "Cd4"), pt.size=0, assay="SCT", 
cols = c('Macrophages (MF.II-480HI)'='#990000', 'B cells (B1a)'='#000099',
'Macrophages (MF.103-11B+24-)'='#FF3333', 'B cells (B1A)'='#33FFFF',
'T cells (T.CD4TESTCJ)'='#999900'))

VlnPlot(object = combined, features = c("Cd19", "Cd79a"), pt.size=0, assay="SCT",
cols = c('Macrophages (MF.II-480HI)'='#990000', 'B cells (B1a)'='#000099',
'Macrophages (MF.103-11B+24-)'='#FF3333', 'B cells (B1A)'='#33FFFF',
'T cells (T.CD4TESTCJ)'='#999900'))

VlnPlot(object = combined, features = c("Cd80", "Ms4a1"), pt.size=0, assay="SCT",
cols = c('Macrophages (MF.II-480HI)'='#990000', 'B cells (B1a)'='#000099',
'Macrophages (MF.103-11B+24-)'='#FF3333', 'B cells (B1A)'='#33FFFF',
'T cells (T.CD4TESTCJ)'='#999900'))

VlnPlot(object = combined, features = c("Itgam", "Itgax"), pt.size=0, assay="SCT",
cols = c('Macrophages (MF.II-480HI)'='#990000', 'B cells (B1a)'='#000099',
'Macrophages (MF.103-11B+24-)'='#FF3333', 'B cells (B1A)'='#33FFFF',
'T cells (T.CD4TESTCJ)'='#999900'))

VlnPlot(object = combined, features = c("Adgre1", "H2-Ab1"), pt.size=0, assay="SCT",
cols = c('Macrophages (MF.II-480HI)'='#990000', 'B cells (B1a)'='#000099',
'Macrophages (MF.103-11B+24-)'='#FF3333', 'B cells (B1A)'='#33FFFF',
'T cells (T.CD4TESTCJ)'='#999900'))

VlnPlot(object = combined, features = c("Mertk", "Icam2"), pt.size=0, assay="SCT", 
cols = c('Macrophages (MF.II-480HI)'='#990000', 'B cells (B1a)'='#000099',
'Macrophages (MF.103-11B+24-)'='#FF3333', 'B cells (B1A)'='#33FFFF',
'T cells (T.CD4TESTCJ)'='#999900'))

VlnPlot(object = combined, features = c("Il1b", "Ccr2", "Cd226"), pt.size=0, assay="SCT",
cols = c('Macrophages (MF.II-480HI)'='#990000', 'B cells (B1a)'='#000099',
'Macrophages (MF.103-11B+24-)'='#FF3333', 'B cells (B1A)'='#33FFFF',
'T cells (T.CD4TESTCJ)'='#999900'))

######################## Dot Plots CAD Loci #########################################
musGenes <- read.table(file = 'cad.csv', header=FALSE)
musGenes <- as.character(musGenes$V1)
musGenes <- unique(musGenes)

Idents(combined) <- "stim"
new.cluster.ids <- c("MerKD Ctrl", "MerKD 2hr", "MerKD 6hr", "WT Ctrl", "WT 2hr", "WT 6hr")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
combined@meta.data$stim_dotplot <- combined@active.ident

Idents(combined) <- "singler"
DotPlot(combined, features = musGenes, cols = c("purple", "pink", "red", "royal blue", "cyan", "lightseagreen"), split.by="stim_dotplot", dot.scale=20, assay="SCT")

NCAD <- c("Pltp", "Rab5c", "Rac1", "Vegfa", "Cdkn1a", "Fndc3b", "Uba5")
CAD <- c("Fndc3b", "Osm", "Zeb2", "Sod2", "Fn1", "Ler3", "Lpl", "Furin", "Lipa", "Lrp1", "Camk2g", "Adam15")

musGenes <- c("Cd3e")
DotPlot(combined, features = musGenes, cols = c("purple", "pink", "red", "royal blue", "cyan", "lightseagreen"), split.by="stim_dotplot", dot.scale=20, assay="SCT")

######################## Dot Plots ###############################################
DotPlot(combined, features = factor(c("Cd3e", "Cd4", "Cd8"), levels=c("Cd8", "Cd4", "Cd3e")), dot.scale=30)  ## T Cell Markers

DotPlot(combined, features = factor(c("Cd19", "Cd79a", "Cd80", "Ms4a1"), levels=c("Cd19", "Cd79a", "Cd80", "Ms4a1")), dot.scale=30)  ## B Cell Markers

DotPlot(combined, features = c("Itgam", "Itgax", "Adgre1", "H2-Ab1", "Mertk", "Icam2", "Il1b", "Ccr2", "Cd226"),
        group.by=, dot.scale=30)  ## Macrophages Markers

tcell_markers <- c("Cd3e", "Cd4")
bcell_markers <- c("Cd19", "Cd79a", "Cd80", "Ms4a1")
macrophage_markers <- c("Adgre1", "Ccr2", "H2-Ab1", "Icam2", "Il1b", "Itgam", "Itgax")

DotPlot(combined, features = c(tcell_markers, bcell_markers, macrophage_markers), 
split.by= c(rep("T Cell Markers", length(tcell_markers)), 
rep("B Cell Markers", length(bcell_markers)), rep("Macrophage Markers", length(macrophage_markers))), dot.scale=20)

DotPlot(combined, features = c(tcell_markers, bcell_markers, macrophage_markers), group.by="singler", dot.scale=20,
cols=c("red", "blue", "green", "orange", "pink"))

DotPlot(combined, features = c(tcell_markers, bcell_markers, macrophage_markers), split.by="singler", dot.scale=20,
cols=c("red", "blue", "green", "orange", "pink"))

##################################################################################
## Diffexp Markers
Idents(combined) <- combined$integrated_snn_res.0.2
combined$celltype.stim <- paste(Idents(object = combined), combined$stim, sep = "_")
combined$celltype <- Idents(object = combined)
Idents(object = combined) <- "celltype.stim"

## Create Marker List for Each Cluster (1, 2, 3, 4)
for (clus in c('B cells (B1a)', 'Macrophages (MF.II-480HI)', 'T cells (T.CD4TESTCJ)', 'Macrophages (MF.103-11B+24-)', 'B cells (B1A)'))
{
    for (sample in c('1', '2', '3', '4', '5'))
    {
        if (sample == '1') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = combined, ident.1 = paste(clus, 'WT_Ctrl', sep='_'), ident.2 = paste(clus, 'MerKD_Ctrl', sep='_'),
                                verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'res0.2_markers.xlsx', sep='_'), sheetName = "wtctrl_v_merkdctrl", append=FALSE)
        }
        else if (sample == '2') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = combined, ident.1 = paste(clus, 'WT_2hr', sep='_'), ident.2 = paste(clus, 'MerKD_2hr', sep='_'),
                                verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'res0.2_markers.xlsx', sep='_'), sheetName = "wt2hr_v_merkd2hr", append=TRUE)
        }
        else if (sample == '3') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = combined, ident.1 = paste(clus, 'WT_6hr', sep='_'), ident.2 = paste(clus, 'MerKD_6hr', sep='_'),
                                verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'res0.2_markers.xlsx', sep='_'), sheetName = "wt6hr_v_merkd6hr", append=TRUE)
        }
        else if (sample == '4') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = combined, ident.1 = paste(clus, 'WT_Ctrl', sep='_'), ident.2 = paste(clus, 'MerKD_6hr', sep='_'),
                                verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'res0.2_markers.xlsx', sep='_'), sheetName = "wtctrl_v_merkd6hr", append=TRUE)
        }
        else if (sample == '5') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = combined, ident.1 = paste(clus, MerKD_Ctrl', sep='_'), ident.2 = paste(clus, 'MerKD_6hr', sep='_'),
                                verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'res0.2_markers.xlsx', sep='_'), sheetName = "merkdctrl_v_merkd6hr", append=TRUE)
        }
    }
}

## FindAllMarkers
markers <- FindAllMarkers(combined)

### Create XSLS File
for (clus in c('B cells (B1a)', 'Macrophages (MF.II-480HI)', 'T cells (T.CD4TESTCJ)', 'Macrophages (MF.103-11B+24-)', 'B cells (B1A)'))
{
    markers %>% dplyr::filter(cluster==clus & p_val_adj <= 0.05) -> val
    if (clus == 0) {
        rownames(val) <- val$gene
        write.xlsx(val, file="combined_markers_res0.2.xlsx", sheetName = clus, append=FALSE)
    }
    else {
        rownames(val) <- val$gene
        write.xlsx(val, file="combined_markers_res0.2.xlsx", sheetName = clus, append=TRUE)
    }
}



## Diffexp Markers
Idents(combined) <- combined@meta.data$singler
combined$celltype.stim <- paste(Idents(object = combined), combined$stim, sep = "_")
combined$celltype <- Idents(object = combined)
Idents(object = combined) <- "celltype.stim"

## Create Marker List for Each Cluster (1, 2, 3, 4)
for (clus in c('B cells (B1a)', 'Macrophages (MF.II-480HI)', 'T cells (T.CD4TESTCJ)', 'Macrophages (MF.103-11B+24-)', 'B cells (B1A)'))
{
    for (sample in c('1', '2', '3', '4', '5'))
    {
        if (sample == '1') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = combined, ident.1 = paste(clus, 'WT_Ctrl', sep='_'), ident.2 = paste(clus, 'MerKD_Ctrl', sep='_'),
                                verbose = FALSE, logfc.threshold=0.1)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'singler_res0.2_markers.xlsx', sep='_'), sheetName = "wtctrl_v_merkdctrl", append=FALSE)
        }
        else if (sample == '2') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = combined, ident.1 = paste(clus, 'WT_2hr', sep='_'), ident.2 = paste(clus, 'MerKD_2hr', sep='_'),
                                verbose = FALSE, logfc.threshold=0.1)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'singler_res0.2_markers.xlsx', sep='_'), sheetName = "wt2hr_v_merkd2hr", append=TRUE)
        }
        else if (sample == '3') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = combined, ident.1 = paste(clus, 'WT_6hr', sep='_'), ident.2 = paste(clus, 'MerKD_6hr', sep='_'),
                                verbose = FALSE, logfc.threshold=0.1)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'singler_res0.2_markers.xlsx', sep='_'), sheetName = "wt6hr_v_merkd6hr", append=TRUE)
        }
        else if (sample == '4') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = combined, ident.1 = paste(clus, 'WT_Ctrl', sep='_'), ident.2 = paste(clus, 'MerKD_6hr', sep='_'),
                                verbose = FALSE, logfc.threshold=0.1)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'singler_res0.2_markers.xlsx', sep='_'), sheetName = "wtctrl_v_merkd6hr", append=TRUE)
        }
        else if (sample == '5') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = combined, ident.1 = paste(clus, 'MerKD_Ctrl', sep='_'), ident.2 = paste(clus, 'MerKD_6hr', sep='_'),
                                verbose = FALSE, logfc.threshold=0.1)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'singler_res0.2_markers.xlsx', sep='_'), sheetName = "merkdctrl_v_merkd6hr", append=TRUE)
        }
    }
}

##################################################################################
## DiffExp VlnPlots for WT_Ctrl v. MerKD_Ctrl, MerKD_2hr, MerKD_6hr

## Set Idents Prior to Create Diff Exp Plots
Idents(combined) <- combined@meta.data$stim

## Diff Exp Violin Plots for Il1b
plots <- VlnPlot(object = combined, features = c("Il1b"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_il1b_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Il6
plots <- VlnPlot(object = combined, features = c("Il6"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"), split.by="stim",
                 cols=c('purple', 'grey'), pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_il6_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Tnf
plots <- VlnPlot(object = combined, features = c("Tnf"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"), split.by="stim",
                 cols=c('purple', 'grey'), pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_tnf_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Tgfb1
plots <- VlnPlot(object = combined, features = c("Tgfb1"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 split.by="stim", cols=c('purple', 'grey'), pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_tgfb1_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Il10
plots <- VlnPlot(object = combined, features = c("Il10"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"), split.by="stim",
                 cols=c('purple', 'grey'), pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_il10_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Gas2l3
plots <- VlnPlot(object = combined, features = c("Gas2l3"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_gas2l3_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Cdc42ep2
plots <- VlnPlot(object = combined, features = c("Cdc42ep2"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_cdc42ep2_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Ssh1
plots <- VlnPlot(object = combined, features = c("Ssh1"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"), split.by="stim",
                 cols=c('purple', 'grey'), pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_ssh1_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Pcdh7
plots <- VlnPlot(object = combined, features = c("Pcdh7"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 split.by="stim", cols=c('purple', 'grey'), pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_pcdh7_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Actr2
plots <- VlnPlot(object = combined, features = c("Actr2"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_actr2_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Nrg2
plots <- VlnPlot(object = combined, features = c("Nrg2"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"), split.by="stim",
                 cols=c('purple', 'grey'), pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_nrg2_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Elmo3
plots <- VlnPlot(object = combined, features = c("Elmo3"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_elmo3_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Pak2
plots <- VlnPlot(object = combined, features = c("Pak2"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"), split.by="stim",
                 cols=c('purple', 'grey'), pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_pak2_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Syne3
plots <- VlnPlot(object = combined, features = c("Syne3"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_syne3_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Tubb4b
plots <- VlnPlot(object = combined, features = c("Tubb4b"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_tubb4b_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Cd44
plots <- VlnPlot(object = combined, features = c("Cd44"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"), split.by="stim",
                 cols=c('purple', 'grey'), pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_cd44_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Vim
plots <- VlnPlot(object = combined, features = c("Vim"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"), split.by="stim",
                 cols=c('purple', 'grey'), pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_vim_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Swap70
plots <- VlnPlot(object = combined, features = c("Swap70"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_swap70_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Col16a1
plots <- VlnPlot(object = combined, features = c("Col16a1"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_col16a1_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Itag3
plots <- VlnPlot(object = combined, features = c("Itag3"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_itag3_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Ssx2ip
plots <- VlnPlot(object = combined, features = c("Ssx2ip"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_ssx2ip_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Dock4
plots <- VlnPlot(object = combined, features = c("Dock4"), , group.by="seurat_clusters", idents = c("WT_2hr", "MerKD_2hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_dock4_wt2hr_v_merkd2hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Il6
plots <- VlnPlot(object = combined, features = c("Il6"), , group.by="singler", idents = c("MerKD_Ctrl", "MerKD_6hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_Il6_merkdctrl_v_merkd6hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Serpine1
plots <- VlnPlot(object = combined, features = c("Serpine1"), , group.by="singler", idents = c("MerKD_Ctrl", "MerKD_6hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_Serpine1_merkdctrl_v_merkd6hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Tnf
plots <- VlnPlot(object = combined, features = c("Tnf"), , group.by="singler", idents = c("MerKD_Ctrl", "MerKD_6hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_Tnf_merkdctrl_v_merkd6hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Ccl5
plots <- VlnPlot(object = combined, features = c("Ccl5"), , group.by="singler", idents = c("MerKD_Ctrl", "MerKD_6hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_Ccl5_merkdctrl_v_merkd6hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Marco
plots <- VlnPlot(object = combined, features = c("Marco"), , group.by="singler", idents = c("MerKD_Ctrl", "MerKD_6hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_Marco_merkdctrl_v_merkd6hr_res0.2.png", width=40, height=40, units="cm")

## Diff Exp Violin Plots for Ly6a
plots <- VlnPlot(object = combined, features = c("Ccl4"), , group.by="singler", idents = c("MerKD_Ctrl", "MerKD_6hr"),
                 cols=c('purple', 'grey'), split.by="stim", pt.size=0, combine=FALSE, assay="SCT")
ggsave("diffexp_vlnplot_Ccl4_merkdctrl_v_merkd6hr_res0.2.png", width=40, height=40, units="cm")


########################################################################################################
## FeaturePlot
featureplot <- FeaturePlot(object = combined, features = c("Dcstamp"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_dcstamp_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Cd3e"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_cd3e_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Cd4"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_cd4_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Cd8"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_cd8_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Cd19"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_cd19_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Cd79a"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_cd79a_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Cd80"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_cd80_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Ms4a1"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_ms4a1_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Itgam"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_itgam_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Itgax"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_itgax_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Adgre1"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_adgre1_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("H2-Ab1"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_h2-ab1_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Mertk"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_mertk_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Icam2"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_icam2_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Il1b"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_il1b_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Ccr2"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_ccr2_res0.2.png", featureplot, width=40, height=40, units="cm")

featureplot <- FeaturePlot(object = combined, features = c("Cd226"), cols = c("grey", "red"), pt.size=0.2)
ggsave("featureplot_cd226_res0.2.png", featureplot, width=40, height=40, units="cm")

#########################################################################################################
library(Seurat)
library(pheatmap)
library(RColorBrewer)
library(tools)

## Create Heatmap w/ Average Expression: Combined
Idents(combined) <- "singler"
cluster.averages <- AverageExpression(combined, assays='integrated', return.seurat=TRUE)
avgs <- cluster.averages@assays$integrated@scale.data
avgs.df <- as.data.frame(avgs)

cadloci <- read.csv(file="/media/data/projects/LinThorp/vst_normalize_wout_regress/scanpy/cad_loci.csv", header=FALSE, sep=',')$V1
cad <- as.character(cadloci)
cad <- tolower(cad)
cad <- toTitleCase(cad)
write.table(cad, file="cad.csv", row.names=FALSE, col.names=FALSE)

novelcadloci <- read.csv(file="/media/data/projects/LinThorp/vst_normalize_wout_regress/scanpy/novel_cad_loci.csv", header=FALSE, sep=',')$V1
ncad <- as.character(novelcadloci)
ncad <- tolower(ncad)
ncad <- toTitleCase(ncad)
write.table(ncad, file="ncad.csv", row.names=FALSE, col.names=FALSE)

test <- avgs.df[cad,
                c('B cells (B1a)', 'B cells (B1A)', 'T cells (T.CD4TESTCJ)', 'Macrophages (MF.II-480HI)', 'Macrophages (MF.103-11B+24-)')]

test <- as.matrix(test)
test <- test[rowSums(is.na(test)) != ncol(test), ]
pheatmap(test, color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100), cluster_cols=FALSE, cluster_rows=FALSE)

## Create Heatmap w/ Average Expression: WT_CTRL, WT_2hr, WT_6hr, MerKD_Ctrl, MerKD_2hr, MerKD_6hr
Idents(combined) <- 'stim'
test <- SubsetData(combined, ident.use="MerKD_6hr")

Idents(test) <- 'singler'
cluster.averages <- AverageExpression(test, assays='integrated', return.seurat=TRUE)
avgs <- cluster.averages@assays$integrated@scale.data
avgs.df <- as.data.frame(avgs)

test <- avgs.df[cad,
                c('B cells (B1a)', 'B cells (B1A)', 'T cells (T.CD4TESTCJ)', 'Macrophages (MF.II-480HI)', 'Macrophages (MF.103-11B+24-)')]

test <- as.matrix(test)
test <- test[rowSums(is.na(test)) != ncol(test), ]
pheatmap(test, color = colorRampPalette(brewer.pal(n = 9, name ="YlOrRd"))(100), cluster_cols=FALSE, cluster_rows=FALSE)

#########################################################################################################
## Subset Macrophage Clusters
Idents(combined) <- "singler"
macrophage <- SubsetData(combined, ident.use=c('Macrophages (MF.II-480HI)', 'Macrophages (MF.103-11B+24-)'))
macrophage <- NormalizeData(object = macrophage, verbose = FALSE)
macrophage <- FindVariableFeatures(object = macrophage, selection.method = "vst", nfeatures = 2000)
macrophage <- ScaleData(macrophage)
macrophage <- RunPCA(macrophage, npcs=30)
macrophage <- FindNeighbors(macrophage, reduction="pca", dims=1:20)
macrophage <- FindClusters(macrophage, resolution=0.2)

Idents(macrophage) <- macrophage@meta.data$seurat_clusters

immgen.ref <- ImmGenData()
test <- as.SingleCellExperiment(macrophage)

common <- intersect(rownames(test), rownames(immgen.ref))
immgen.ref <- immgen.ref[common,]
test <- test[common,]
test <- logNormCounts(test)

clusters <- macrophage@active.ident

pred <- SingleR(test=test, ref=immgen.ref, labels=immgen.ref$label.fine, method="cluster", clusters=clusters)

## Add Seurat New Cluster IDs
new.cluster.ids <- pred$pruned.labels

names(new.cluster.ids) <- levels(macrophage)
macrophage <- RenameIdents(macrophage, new.cluster.ids)
DimPlot(macrophage, reduction="umap", label=TRUE)
macrophage@meta.data$singler <- Idents(macrophage)

Idents(macrophage) <- "singler"
DimPlot(macrophage, reduction="umap", label=FALSE)

##############################################
## Monocle 3: Set Identities to Seurat Clusters
Idents(macrophage) <- "seurat_clusters"

## Extract counts, phenotype data, and feature data from SeuratObject
data <- as(as.matrix(macrophage@assays$integrated@data), 'sparseMatrix')
pd <- macrophage@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

## Construct Monocle CDS
cds <- new_cell_data_set(data,
                         cell_metadata = pd,
                         gene_metadata = fData)

cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds)

## Use cell.cluster Labels w/ Monocle3 UMAP
plot_cells(cds, color_cells_by="singler")

plot_cells(cds, color_cells_by="seurat_clusters", group_label_size=4)

cds <- cluster_cells(cds)

cds <- learn_graph(cds)

plot <- plot_cells(cds, color_cells_by = "stim", group_label_size=4,
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE, label_branch_points=FALSE, cell_size=1, label_cell_groups=FALSE) +
    scale_color_manual(values = c('MerKD_2hr'='#990000', 'MerKD_6hr'='#D3D3D3', 'MerKD_Ctrl'='#D3D3D3',
                                  'WT_2hr'='#000099', 'WT_6hr'='#D3D3D3', 'WT_Ctrl'='#D3D3D3'))

plot <- plot_cells(cds, color_cells_by = "singler", group_label_size=4,
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE, label_branch_points=FALSE, cell_size=1)
ggsave("monocle3_plot_singler_labeled.png", plot, width=40, height=40, units="cm")

plot <- plot_cells(cds, color_cells_by = "stim", group_label_size=4,
                   label_groups_by_cluster=TRUE, label_cell_groups=FALSE,
                   label_leaves=FALSE, label_branch_points=FALSE, cell_size=1)
ggsave("monocle3_plot_singler_labeled.png", plot, width=40, height=40, units="cm")

## Create FindConservedMarkers Function
DefaultAssay(object = macrophage) <- "RNA"
clus = 0
while (clus < length(levels(macrophage@meta.data$seurat_clusters)))
{
    print(paste("Cluster", clus, sep=" "))
    markers <- FindConservedMarkers(object = macrophage, ident.1=clus, grouping.var="stim", verbose=FALSE)
    if (clus == 0) {
        write.xlsx(markers, file="macrophage_conserved_markers.xlsx", sheetName = paste("Cluster", clus, sep="_"), append=FALSE)
    }
    else {
        write.xlsx(markers, file="macrophage_conserved_markers.xlsx", sheetName = paste("Cluster", clus, sep="_"), append=TRUE)
    }
    clus = clus + 1
}

#############################################################################
## SWNE Plot For Macrophage Subcluster
library(Seurat)
library(swne)

## Extract integrated expression matrix
aligned.counts <- as.matrix(GetAssayData(macrophage, assay = "integrated"))

dim(aligned.counts)

hist(as.matrix(aligned.counts))

## Reconstruct matrix to remove negative values
aligned.counts <- t(apply(aligned.counts, 1, function(x) (x - min(x))/(max(x) - min(x))))

### Run NMF
k <- 20
n.cores <- 14
nmf.res <- RunNMF(aligned.counts, k = k, alpha = 0, init = "ica", n.cores = n.cores)

pc.emb <- Embeddings(macrophage, "pca")
snn <- CalcSNN(t(pc.emb), k = 20, prune.SNN = 0)

swne.embedding <- EmbedSWNE(nmf.res$H, snn, alpha.exp = 1.5, snn.exp = 1, n_pull = 3)

norm.counts <- ExtractNormCounts(macrophage, rescale = T, rescale.method = "log", batch = NULL)

nmf.res$W <- ProjectFeatures(norm.counts, nmf.res$H, n.cores = n.cores)

gene.factor.df <- SummarizeAssocFeatures(nmf.res$W, features.return = 1)
genes.embed <- unique(gene.factor.df$feature)

swne.embedding <- EmbedFeatures(swne.embedding, nmf.res$W, genes.embed, n_pull = 3)
swne.embedding$H.coords$name <- ""
swne.embedding$feature.coords$name <- ""

Idents(macrophage) <- "singler"
clusters <- macrophage@active.ident

swneplot = PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = clusters, do.label = F,
                label.size = 3.5, pt.size = 1.25, show.legend = T, seed = 42,
                colors.use = c('Macrophages (MF.II-480HI)'='#990000', 'DC (DC.103-11B+24+)'='#000099', 
                'Macrophages (MF.103-11B+24-)'='#FF3333'))
ggsave("swne_plot_colored.png", swneplot, width=40, height=40, units="cm")

################ Seperate SWNE Per Library ######################################
Idents(macrophage) <- "stim"
clusters <- macrophage@active.ident

merkd6hr <- clusters[clusters=="MerKD_6hr"]
Idents(macrophage) <- "seurat_clusters"
clusters <- macrophage@active.ident

test <- clusters[names(merkd6hr)]

Idents(macrophage) <- "singler"
clusters <- macrophage@active.ident

test <- clusters[names(merkd6hr)]
swneplot = PlotSWNE(swne.embedding, alpha.plot = 0.4, sample.groups = test, do.label = F,
                    label.size = 3.5, pt.size = 2, show.legend = T, seed = 42)

################################################################################
## Diffexp Markers
Idents(macrophage) <- "seurat_clusters"
macrophage$celltype.stim <- paste(Idents(object = macrophage), macrophage$stim, sep = "_")
macrophage$celltype <- Idents(object = macrophage)
Idents(object = macrophage) <- "celltype.stim"

## Create Marker List for Each Cluster (1, 2, 3, 4)
for (clus in c('0', '1', '2', '3'))
{
    for (sample in c('1', '2', '3'))
    {
        if (sample == '1') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = macrophage, ident.1 = paste(clus, 'WT_Ctrl', sep='_'), ident.2 = paste(clus, 'MerKD_Ctrl', sep='_'),
                                verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'macrophage_res0.2_markers.xlsx', sep='_'), sheetName = "wtctrl_v_merkdctrl", append=FALSE)
        }
        else if (sample == '2') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = macrophage, ident.1 = paste(clus, 'WT_Ctrl', sep='_'), ident.2 = paste(clus, 'MerKD_2hr', sep='_'),
                                verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'macrophage_res0.2_markers.xlsx', sep='_'), sheetName = "wtctrl_v_merkd2hr", append=TRUE)
        }
        else if (sample == '3') {
            print(paste('Cluster', clus, sep=' '))
            cluster <- paste('cluster', clus, sep='')
            mark <- FindMarkers(object = macrophage, ident.1 = paste(clus, 'WT_Ctrl', sep='_'), ident.2 = paste(clus, 'MerKD_6hr', sep='_'),
                                verbose = FALSE)
            mark %>% dplyr::mutate(genes = rownames(mark)) %>% dplyr::arrange(desc(abs(avg_logFC))) %>% dplyr::filter(p_val_adj <= 0.05) %>%
                dplyr::arrange(desc(avg_logFC)) -> c
            rownames(c) <- c$genes
            write.xlsx(c, file=paste(cluster, 'macrophage_res0.2_markers.xlsx', sep='_'), sheetName = "wtctrl_v_merkd6hr", append=TRUE)
        }
    }
}

##############################################################################################
## Subset Macrophage Subset into WT_Ctrl, WT_2hr, WT_6hr & MerKD_Ctrl, MerKD_2hr, MerKD_6hr
################################## Monocle3 ##################################################
library(Seurat)
library(monocle3)
library(reticulate)
library(hdf5r)
library(dplyr)

Idents(macrophage) <- "orig.ident"

wtcells <- subset(macrophage, idents=c("wt_ctrl", "wt_2hr", "wt_6hr"))

## Set Identities to cluster.id
Idents(wtcells) <- wtcells@meta.data$seurat_clusters

## Extract counts, phenotype data, and feature data from SeuratObject
data <- as(as.matrix(wtcells@assays$integrated@data), 'sparseMatrix')
pd <- wtcells@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

## Construct Monocle CDS
cds <- new_cell_data_set(data,
                         cell_metadata = pd,
                         gene_metadata = fData)

cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds)

## Use cell.cluster Labels w/ Monocle3 UMAP
plot_cells(cds, color_cells_by="singler")

plot_cells(cds, color_cells_by="seurat_clusters", group_label_size=4)

cds <- cluster_cells(cds)

cds <- learn_graph(cds)

plot <- plot_cells(cds, color_cells_by = "orig.ident", group_label_size=4,
                   label_groups_by_cluster=TRUE, label_cell_groups=FALSE,
                   label_leaves=FALSE, label_branch_points=FALSE, cell_size=1.5)
ggsave("monocle3_plot_singler_labeled.png", plot, width=40, height=40, units="cm")

###### MerKD Cells ######################################
merkdcells <- subset(macrophage, idents=c("MerKD_Ctrl", "merkd_2hr", "merkd_6hr"))

## Set Identities to cluster.id
Idents(merkdcells) <- merkdcells@meta.data$seurat_clusters

## Extract counts, phenotype data, and feature data from SeuratObject
data <- as(as.matrix(merkdcells@assays$integrated@data), 'sparseMatrix')
pd <- merkdcells@meta.data
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))

## Construct Monocle CDS
cds <- new_cell_data_set(data,
                         cell_metadata = pd,
                         gene_metadata = fData)

cds <- preprocess_cds(cds, num_dim = 100)

cds <- reduce_dimension(cds)

## Use cell.cluster Labels w/ Monocle3 UMAP
plot_cells(cds, color_cells_by="singler")

plot_cells(cds, color_cells_by="seurat_clusters", group_label_size=4)

cds <- cluster_cells(cds)

cds <- learn_graph(cds)

plot <- plot_cells(cds, color_cells_by = "seurat_clusters", group_label_size=4,
                   label_groups_by_cluster=TRUE, label_cell_groups=FALSE,
                   label_leaves=FALSE, label_branch_points=FALSE, cell_size=1.5)
ggsave("monocle3_plot_singler_labeled.png", plot, width=40, height=40, units="cm")

