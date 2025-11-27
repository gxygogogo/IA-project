#%% 导入R包
library(limma)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(WGCNA)
library(forcats)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(tidyr)
library(GOplot)
library(Seurat)
library(harmony)
library(RColorBrewer)
library(paletteer)
library(ggrastr)
library(CellChat)

#%% 导入数据
GSM5813881_count = Read10X_h5('/public3/Xinyu/IA/scRNA-seq/GSM5813881_sham_filtered_feature_bc_matrix.h5')
GSM5813883_count = Read10X_h5('/public3/Xinyu/IA/scRNA-seq/GSM5813883_formed_filtered_feature_bc_matrix.h5')
GSM5813885_count = Read10X_h5('/public3/Xinyu/IA/scRNA-seq/GSM5813885_ruptured_filtered_feature_bc_matrix.h5')

GSM5813881_obj = CreateSeuratObject(counts = GSM5813881_count, project = 'GSE193533_obj', min.cells = 3, min.features = 200)
GSM5813883_obj = CreateSeuratObject(counts = GSM5813883_count, project = 'GSE193533_obj', min.cells = 3, min.features = 200)
GSM5813885_obj = CreateSeuratObject(counts = GSM5813885_count, project = 'GSE193533_obj', min.cells = 3, min.features = 200)

GSM5813881_obj = PercentageFeatureSet(GSM5813881_obj, pattern = "^MT-", col.name = "percent.mt")
GSM5813881_obj = PercentageFeatureSet(GSM5813881_obj, pattern = "^HB[AB]-", col.name = "percent.hb")
GSM5813883_obj = PercentageFeatureSet(GSM5813883_obj, pattern = "^MT-", col.name = "percent.mt")
GSM5813883_obj = PercentageFeatureSet(GSM5813883_obj, pattern = "^HB[AB]-", col.name = "percent.hb")
GSM5813885_obj = PercentageFeatureSet(GSM5813885_obj, pattern = "^MT-", col.name = "percent.mt")
GSM5813885_obj = PercentageFeatureSet(GSM5813885_obj, pattern = "^HB[AB]-", col.name = "percent.hb")

GSM5813881_obj$log10GenesPerUMI = log10(GSM5813881_obj$nFeature_RNA) / log10(GSM5813881_obj$nCount_RNA)
GSM5813881_obj = subset(GSM5813881_obj, subset = nFeature_RNA > 200 & percent.mt < 10 & percent.hb < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.7)
GSM5813883_obj$log10GenesPerUMI = log10(GSM5813883_obj$nFeature_RNA) / log10(GSM5813883_obj$nCount_RNA)
GSM5813883_obj = subset(GSM5813883_obj, subset = nFeature_RNA > 200 & percent.mt < 10 & percent.hb < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.7)
GSM5813885_obj$log10GenesPerUMI = log10(GSM5813885_obj$nFeature_RNA) / log10(GSM5813885_obj$nCount_RNA)
GSM5813885_obj = subset(GSM5813885_obj, subset = nFeature_RNA > 200 & percent.mt < 10 & percent.hb < 5 & nCount_RNA > 1000 & log10GenesPerUMI > 0.7)

# GSM5813881
genes_use = rownames(GSM5813881_obj[['RNA']]@counts)
genes_use = genes_use[which(!grepl('^MT-', genes_use))]
genes_use = genes_use[which(!grepl('^RPS', genes_use))]
genes_use = genes_use[which(!grepl('^RPL', genes_use))]
GSM5813881_obj = subset(GSM5813881_obj, features=genes_use)
GSM5813881_obj = NormalizeData(GSM5813881_obj)
# GSM5813883
genes_use = rownames(GSM5813883_obj[['RNA']]@counts)
genes_use = genes_use[which(!grepl('^MT-', genes_use))]
genes_use = genes_use[which(!grepl('^RPS', genes_use))]
genes_use = genes_use[which(!grepl('^RPL', genes_use))]
GSM5813883_obj = subset(GSM5813883_obj, features=genes_use)
GSM5813883_obj = NormalizeData(GSM5813883_obj)
# GSM5813885
genes_use = rownames(GSM5813885_obj[['RNA']]@counts)
genes_use = genes_use[which(!grepl('^MT-', genes_use))]
genes_use = genes_use[which(!grepl('^RPS', genes_use))]
genes_use = genes_use[which(!grepl('^RPL', genes_use))]
GSM5813885_obj = subset(GSM5813885_obj, features=genes_use)
GSM5813885_obj = NormalizeData(GSM5813885_obj)

GSM5813881_obj = RenameCells(GSM5813881_obj, new.names = paste0("GSM5813881_", colnames(GSM5813881_obj)))
GSM5813883_obj = RenameCells(GSM5813883_obj, new.names = paste0("GSM5813883_", colnames(GSM5813883_obj)))
GSM5813885_obj = RenameCells(GSM5813885_obj, new.names = paste0("GSM5813885_", colnames(GSM5813885_obj)))

GSM5813881_obj$type = "sham"
GSM5813883_obj$type = "formed"
GSM5813885_obj$type = "ruptured"

#%% 合并数据
scRNAlist = list(GSM5813881_obj, GSM5813883_obj, GSM5813885_obj)
names(scRNAlist) = c('GSM5813881', 'GSM5813883', 'GSM5813885')

GSE193533_obj = merge(scRNAlist[[1]], y=scRNAlist[2:3])

#%% PCA
## PCA
GSE193533_obj = FindVariableFeatures(GSE193533_obj, selection.method = "vst", nfeatures = 2000)
GSE193533_obj = ScaleData(GSE193533_obj,features = row.names(GSE193533_obj))
GSE193533_obj = RunPCA(object = GSE193533_obj, assay = "RNA", npcs = 50, verbose = TRUE)

pdf("/public3/Xinyu/IA/scRNA-seq/PCA.pdf", width = 6, height = 5)
DimPlot(GSE193533_obj, reduction = "pca", group.by = 'type')
dev.off()

#%% 去除批次效应
GSE193533_obj_harmony = RunHarmony(GSE193533_obj, group.by.vars = "type")

pdf("/public3/Xinyu/IA/scRNA-seq/Harmony_pca.pdf", width = 6, height = 5)
DimPlot(GSE193533_obj_harmony, reduction = "harmony", group.by = "type")
dev.off()

#%% 聚类
GSE193533_obj_harmony = FindNeighbors(GSE193533_obj_harmony, dims = 1:10, reduction = "harmony")
GSE193533_obj_harmony = FindClusters(GSE193533_obj_harmony, resolution = 0.2)
GSE193533_obj_harmony = RunUMAP(GSE193533_obj_harmony, reduction = "harmony", dims = 1:10)

saveRDS(GSE193533_obj_harmony, file = '/public3/Xinyu/IA/scRNA-seq/GSE193533_obj_harmony.Rds')

umap = GSE193533_obj_harmony@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(seurat_clusters = GSE193533_obj_harmony@meta.data$seurat_clusters) # 注释后的label信息 ，改为cell_type
allcolour = distinctColorPalette(13)

p = ggplot(umap,aes(x= UMAP_1 , y = UMAP_2 ,color = seurat_clusters)) +  
        geom_point(size = 1 , alpha =1 )  +  
        scale_color_manual(values = allcolour)+
        theme(panel.grid.major = element_blank(), #主网格线
              panel.grid.minor = element_blank(), #次网格线
              panel.border = element_blank(), #边框
              axis.title = element_blank(),  #轴标题
              axis.text = element_blank(), # 文本
              axis.ticks = element_blank(),
              panel.background = element_rect(fill = 'white'), #背景色
              plot.background=element_rect(fill="white"),
              legend.title = element_blank(), #去掉legend.title 
              legend.key=element_rect(fill='white'), #
              legend.text = element_text(size=20), #设置legend标签的大小
              legend.key.size=unit(1,'cm'))+
        guides(color = guide_legend(override.aes = list(size=5)))+
        geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
        geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3), colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
        annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1", color="black",size = 3, fontface="bold" ) + 
        annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2", color="black",size = 3, fontface="bold" ,angle=90) 
ggsave(p, file = '/public3/Xinyu/IA/scRNA-seq/Cluster_umap.pdf', width = 6, height = 5)

#%% 细胞类型注释
markers = c('Myh11', 'Tagln', 'Acta2', 'Cd68', 'C1qa', 'C1qb', 'S100a9', 'S100a8', 'Lcn2', 'Dcn', 'Col1a1', 'Lum', 'Cdh5', 'Pecam1', 'Esam', 'Cd3d', 'Cd3g', 'Cd28', 'Rgs5',
           'Rgs16', 'Kcnj8', 'Cd79a', 'Cd79b', 'Ly6d', 'Cd209a', 'Ccr7', 'Ifitm1', 'Mcpt4', 'Cma1', 'Tpsb2', 'Igkc', 'Trbc1', 'Trbc2', 'Coro1a', 'Alas2', 'Fech', 'Epb41', 'Cox8a', 'Sirt2', 'Mdh1')

p_list = list()

for(i in 1:length(markers)){
    p = plot_density(GSE193533_obj_harmony, features = markers[i], pal = 'magma', raster = F, size = 0.8) +
            theme_blank() +
            theme(legend.frame = element_rect(colour = 'black'),
                  legend.ticks = element_line(colour = 'black', linewidth = 0),
                  legend.key.width = unit(0.3, 'cm'),
                  legend.key.height = unit(0.6, 'cm'),
                  legend.title = element_text(color = 'black', face = 'bold', size = 8),
                  plot.title = element_text(hjust = 0.5))
    p_list[[i]] = p
}

pdf('/public3/Xinyu/IA/scRNA-seq/FeaturesUMAP.pdf', width = 8, height = 12)
do.call(plot_grid, c(p_list, ncol = 3))
dev.off()

#%% 差异基因
cluster.markers = FindAllMarkers(GSE193533_obj_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0)
# 查看top marker
top10 = cluster.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 = cluster.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

#%% 细胞类型注释
# VSMC(1, 2, 5): Myh11, Tagln, Acta2
# Macrophage(0,11,10): Cd68, C1qa, C1qb
# Neutrophil(3): S100a8, S100a9, Lcn2
# Froblast(4, 5): Dcn, Col1a1, Lum
# Endothelium(7): Cdh5, Pecam1, Esam
# T lymphocyte(8): Cd3d, Cd3g, Cd28
# B lymphocyte(9): Cd79a, Cd79b, Ly6d
# Mast cells(12): Mcpt4, Cma1, Tpsb2
# Erythroid lineage(6): Alas2

markers = c('Myh11', 'Tagln', 'Acta2', 'Cd68', 'C1qa', 'C1qb', 'S100a8', 'S100a9', 'Lcn2', 'Dcn', 'Col1a1', 'Lum', 'Cdh5', 'Pecam1', 'Esam',
            'Cd3d', 'Cd3g', 'Cd28', 'Cd79a', 'Cd79b', 'Ly6d', 'Mcpt4', 'Cma1', 'Tpsb2', 'Alas2')

# markers = c('Myh11', 'Cd68', 'S100a8', 'Dcn', 'Cdh5', 'Cd3d', 'Cd79a', 'Mcpt4', 'Alas2')

GSE193533_obj_harmony@meta.data$CellType = NA
GSE193533_obj_harmony@meta.data[which(GSE193533_obj_harmony@meta.data$seurat_clusters == 3), "CellType"] = "Neutrophil"
GSE193533_obj_harmony@meta.data[which(GSE193533_obj_harmony@meta.data$seurat_clusters == 7), "CellType"] = "Endothelium"
GSE193533_obj_harmony@meta.data[which(GSE193533_obj_harmony@meta.data$seurat_clusters == 8), "CellType"] = "T lymphocyte"
GSE193533_obj_harmony@meta.data[which(GSE193533_obj_harmony@meta.data$seurat_clusters == 9), "CellType"] = "B lymphocyte"
GSE193533_obj_harmony@meta.data[which(GSE193533_obj_harmony@meta.data$seurat_clusters == 6), "CellType"] = "Erythroid lineage"
GSE193533_obj_harmony@meta.data[which(GSE193533_obj_harmony@meta.data$seurat_clusters == 12), "CellType"] = "Mast cells"
GSE193533_obj_harmony@meta.data[which(GSE193533_obj_harmony@meta.data$seurat_clusters %in% c(1, 2, 5)), "CellType"] = "VSMC"
GSE193533_obj_harmony@meta.data[which(GSE193533_obj_harmony@meta.data$seurat_clusters %in% c(0, 10, 11)), "CellType"] = "Macrophage"
GSE193533_obj_harmony@meta.data[which(GSE193533_obj_harmony@meta.data$seurat_clusters %in% c(4)), "CellType"] = "Firoblast"

#%% marker热图绘制
markers = c('Myh11', 'Tagln', 'Acta2', 'Cd68', 'C1qa', 'C1qb', 'S100a8', 'S100a9', 'Lcn2', 'Dcn', 'Col1a1', 'Lum', 'Cdh5', 'Pecam1', 'Esam',
            'Cd3d', 'Cd3g', 'Cd28', 'Cd79a', 'Cd79b', 'Ly6d', 'Mcpt4', 'Cma1', 'Tpsb2', 'Alas2')

#计算平均表达量
gene_cell_exp = AverageExpression(GSE193533_obj_harmony,features = markers, group.by = 'CellType',slot = 'data') 
gene_cell_exp = as.data.frame(gene_cell_exp$RNA)

#顶部细胞类型注释
df = data.frame(colnames(gene_cell_exp))
colnames(df) = 'class'
df$class = factor(df$class, levels = c("Neutrophil", "Endothelium", "T lymphocyte", "B lymphocyte", "Mast cells", "VSMC", "Macrophage", "Firoblast", "Erythroid lineage"))
df = as.data.frame(df[order(df$class), ])
colnames(df) = 'class'

top_anno = HeatmapAnnotation(df = df,#细胞名/cluster
                             border = T,
                             show_annotation_name = F,
                             gp = gpar(col = 'black'),
                             col = list(class = c("Neutrophil" = '#386CB0', 
                                                  "Endothelium" = '#BEAED4', 
                                                  "T lymphocyte" = '#FC8D62', 
                                                  "B lymphocyte" = '#FFFF99', 
                                                  "Mast cells" = '#FBB4AE', 
                                                  "VSMC" = '#E31A1C',
                                                  "Macrophage" = "#BFEFFF",
                                                  "Firoblast" = "#CAE1FF", 
                                                  "Erythroid lineage" = '#E78AC3')))
#数据标准化缩放一下
marker_exp = t(scale(t(gene_cell_exp),scale = T,center = T))

cell = c("Neutrophil", "Endothelium", "T lymphocyte", "B lymphocyte", "Mast cells", "VSMC", "Macrophage", "Firoblast", "Erythroid lineage")

# 确保只选择存在的行和列
markers = c('S100a8', 'S100a9', 'Lcn2', 'Cdh5', 'Pecam1', 'Esam', 'Cd3d', 'Cd3g', 'Cd28', 'Cd79a', 'Cd79b', 'Ly6d', 'Mcpt4', 'Cma1', 'Tpsb2', 'Myh11', 'Tagln', 'Acta2', 'Cd68', 'C1qa', 'C1qb', 'Dcn', 'Col1a1', 'Lum', 'Alas2')
valid_genes = markers[markers %in% rownames(marker_exp)]
valid_cells = cell[cell %in% colnames(marker_exp)]

# 子集选择
marker_exp_subset = marker_exp[valid_genes, valid_cells]

col_cluster = setNames(c(rep("#386CB0", 3), rep("#BEAED4", 3),
                         rep("#FC8D62", 3), rep("#FFFF99", 3), rep('#FBB4AE', 3),
                         rep('#E31A1C', 3), rep("#BFEFFF", 3),
                         rep("#CAE1FF", 3),
                         rep("#E78AC3", 1)),
                        rownames(marker_exp_subset))#设置对应标签颜色

row_info = rowAnnotation(foo = anno_text(rownames(marker_exp_subset), 
                         location = 0, 
                         just = "left",
                         gp = gpar(fill = col_cluster, col = "black"),
                         width = max_text_width(rownames(marker_exp_subset))*1.2))

pdf('/public3/Xinyu/IA/scRNA-seq/MarkerHeatmap.pdf', height = 9, width = 7)
Heatmap(marker_exp_subset,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = T,
        column_title = NULL,
        heatmap_legend_param = list(title=' '),
        col = colorRampPalette(c("#0000EF","black","#FDFE00"))(100),
        border = 'black',
        rect_gp = gpar(col = "black", lwd = 1),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        top_annotation = top_anno)+row_info
dev.off()

#%% 带轮廓的细胞类型UMAP图
plotData = as.data.frame(GSE193533_obj_harmony[['umap']]@cell.embeddings)
plotData$cluster = GSE193533_obj_harmony@meta.data$CellType

color_pal <- c(
"Neutrophil" = '#386CB0', 
"Endothelium" = '#BEAED4', 
"T lymphocyte" = '#FC8D62', 
"B lymphocyte" = '#FFFF99', 
"Mast cells" = '#FBB4AE', 
"VSMC" = '#E31A1C',
"Macrophage" = "#BFEFFF",
"Firoblast" = "#CAE1FF", 
"Erythroid lineage" = '#E78AC3'
)

plot = ggplot(plotData, aes(x = UMAP_1, y = UMAP_2, 
                            fill = GSE193533_obj_harmony@meta.data$CellType, 
                            color = GSE193533_obj_harmony@meta.data$CellType))+
            stat_unchull(alpha = 0.25, size = 0.25, delta = 0.35, lty = 2)+
            geom_point(size = 0.1)+
            theme(aspect.ratio = 1, 
                  panel.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_blank(),
                  axis.title = element_text(hjust = 0.05, face = 'italic', size = 14),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank())+
            labs(fill = 'CellType', color = 'CellType')+
            scale_x_continuous(breaks = NULL) +
            scale_y_continuous(breaks = NULL) +
            scale_fill_manual(values = color_pal) +
            scale_color_manual(values = color_pal) +
            # 手动绘制短坐标轴
            annotate("segment", x = -15, xend = -9, y = -15, yend = -15, color = "black", size = 0.8, arrow = arrow(type = "closed", length = unit(0.3, "cm"))) +
            annotate("segment", x = -15, xend = -15, y = -15, yend = -9, color = "black", size = 0.8, arrow = arrow(type = "closed", length = unit(0.3, "cm")))

pdf("/Volumes/GaoxyData/IAs_object/scRNA-seq/CellType_umap.pdf", width = 8, height = 8)
print(plot)
dev.off()
#%% 分面绘制UMAP图
plotData = as.data.frame(GSE193533_obj_harmony[['umap']]@cell.embeddings)
plotData$cluster = GSE193533_obj_harmony@meta.data$CellType
plotData$group = GSE193533_obj_harmony@meta.data$type
plotData$group = factor(plotData$group, levels = c('sham', 'formed', 'ruptured'))

color_pal <- c(
"Neutrophil" = '#386CB0', 
"Endothelium" = '#BEAED4', 
"T lymphocyte" = '#FC8D62', 
"B lymphocyte" = '#FFFF99', 
"Mast cells" = '#FBB4AE', 
"VSMC" = '#E31A1C',
"Macrophage" = "#BFEFFF",
"Firoblast" = "#CAE1FF", 
"Erythroid lineage" = '#E78AC3'
)

plot = ggplot(plotData, aes(x = UMAP_1, y = UMAP_2, fill = cluster, color = cluster))+
            stat_unchull(alpha = 0.25, size = 0.25, delta = 0.35, lty = 2)+
            geom_point(size = 0.1)+
            theme(aspect.ratio = 1, 
                  panel.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_blank(),
                  axis.title = element_text(hjust = 0.05, face = 'italic', size = 14),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank())+
            labs(fill = 'CellType', color = 'CellType')+
            scale_x_continuous(breaks = NULL) +
            scale_y_continuous(breaks = NULL) +
            scale_fill_manual(values = color_pal) +
            scale_color_manual(values = color_pal) +
            facet_wrap(~group)+
            # 手动绘制短坐标轴
            annotate("segment", x = -15, xend = -9, y = -15, yend = -15, color = "black", size = 0.8, arrow = arrow(type = "closed", length = unit(0.3, "cm"))) +
            annotate("segment", x = -15, xend = -15, y = -15, yend = -9, color = "black", size = 0.8, arrow = arrow(type = "closed", length = unit(0.3, "cm")))

pdf("/Volumes/GaoxyData/IAs_object/scRNA-seq/CellType_umap_facet.pdf", width = 22, height = 8)
print(plot)
dev.off()


saveRDS(GSE193533_obj_harmony, file = '/public3/Xinyu/IA/scRNA-seq/GSE193533_obj_harmony_anno.Rds')

#%% 细胞类型统计
# 假设 CellType 和 orig.ident 都在 meta.data 中
cell_counts <- GSE193533_obj_harmony@meta.data %>%
  group_by(type, CellType) %>%
  summarise(Count = n()) %>%
  ungroup()

# 计算每个样本中细胞类型的比例
cell_props <- cell_counts %>%
  group_by(type) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()
cell_props$type = factor(cell_props$type, levels = c('sham', 'formed', 'ruptured'))

pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/CellType_stat.pdf', width = 7, height = 8)
ggplot(cell_props, aes(x = type, y = Proportion, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = color_pal) +  # 用你之前定义的颜色
  labs(x = "Sample", y = "Proportion", fill = "Cell Type",title = "Cell Type Composition per Sample") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))
dev.off()


#%% 将人类基因SYMBOL转换为老鼠的SYMBOL
purple_up_hub_gene = c('COL1A1', 'COL5A1', 'CHRF', 'FSTL1', 'KDELR3', 'FBN1', 'ADAMTS2', 'RCN3', 'FKBP10', 'CREB3L1')
pink_up_hub_gene = c('KRT8', 'KRT18', 'MFAP5', 'KRT7', 'INHBA', 'LTBP2', 'SPATA18', 'PRRX2', 'PODNL1', 'ITGA11')

purple_up_mouse <- human2mouse(purple_up_hub_gene)
pink_up_mouse   <- human2mouse(pink_up_hub_gene)

FeaturePlot(GSE193533_obj_harmony, features = all_hub_gene, ncol = 3)

#%% hub gene 分析
Fibroblast_obj = subset(GSE193533_obj_harmony, subset = CellType == 'Firoblast')
VSMC_obj = subset(GSE193533_obj_harmony, subset = CellType == 'VSMC')

all_hub_gene = c(purple_up_mouse$mouseGene, 'Mfap5', 'Inhba', 'Ltbp2', 'Prrx2', 'Itga11')

FeaturePlot(GSE193533_obj_harmony, features = all_hub_gene, ncol = 3)

# 提取表达数据
expr_data <- FetchData(Fibroblast_obj, vars = all_hub_gene)
expr_data$Cell <- rownames(expr_data)
expr_data$SampleType <- Fibroblast_obj@meta.data$type

# 长格式转换
expr_long <- expr_data %>%
  pivot_longer(cols = all_of(all_hub_gene), names_to = "Gene", values_to = "Expression")
expr_long$SampleType = factor(expr_long$SampleType, levels = c('sham', 'formed', 'ruptured'))

expr_long$Gene = factor(expr_long$Gene, levels = c('Col5a1', 'Creb3l1', 'Fbn1', 'Fkbp10', 'Inhba', 'Kdelr3', 'Ltbp2', 'Mfap5', 'Rcn3', 'Itga11', 'Adamts2', 'Col1a1', 'Fstl1', 'Prrx2'))

pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/Hub_gene_express_Fibroblast.pdf', width = 16, height = 7)
ggplot(expr_long, aes(x = SampleType, y = Expression, fill = SampleType)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = c('#62D1B7', '#FFA42E', '#FF4F4F'))+
  facet_wrap(~ Gene, scales = "free_y", ncol = 7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank()) +
  labs(x = "Sample Type", y = "Expression", title = "Hub Gene Expression in Fibroblasts")
dev.off()


#%% 在聚类
Fibroblast_obj <- FindVariableFeatures(Fibroblast_obj, selection.method = "vst", nfeatures = 2000)

# 3. 缩放
Fibroblast_obj <- ScaleData(Fibroblast_obj)

# 4. PCA
Fibroblast_obj <- RunPCA(Fibroblast_obj)

# 5. 构建邻接图
Fibroblast_obj <- FindNeighbors(Fibroblast_obj, dims = 1:10)

# 6. 聚类（关键步骤）👉 设置 `resolution = 0.01` 等低值，通常可以得到2类
Fibroblast_obj <- FindClusters(Fibroblast_obj, resolution = 0.05)

# 7. 可视化聚类
Fibroblast_obj <- RunUMAP(Fibroblast_obj, dims = 1:10)

DimPlot(Fibroblast_obj, reduction = "umap", label = TRUE)
Fibroblast_obj$subType = NA
Fibroblast_obj@meta.data[which(Fibroblast_obj@meta.data$seurat_clusters == 0), 'subType'] = 'Fibroblast cluster1'
Fibroblast_obj@meta.data[which(Fibroblast_obj@meta.data$seurat_clusters != 0), 'subType'] = 'Fibroblast cluster2'


F1_obj = subset(Fibroblast_obj, subset = subType == 'Fibroblast cluster1')
F2_obj = subset(Fibroblast_obj, subset = subType == 'Fibroblast cluster2')

#%% UMAP图
plotData = as.data.frame(Fibroblast_obj[['umap']]@cell.embeddings)
plotData$cluster = Fibroblast_obj@meta.data$subType

color_pal <- c(
"Fibroblast cluster1" = '#55CCFF', 
"Fibroblast cluster2" = '#F26101'
)

plot = ggplot(plotData, aes(x = UMAP_1, y = UMAP_2, 
                            fill = Fibroblast_obj@meta.data$subType, 
                            color = Fibroblast_obj@meta.data$subType))+
            geom_point(size = 1)+
            theme(aspect.ratio = 1, 
                  panel.background = element_blank(),
                  panel.grid = element_blank(),
                  axis.line = element_blank(),
                  axis.title = element_text(hjust = 0.05, face = 'italic', size = 14),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank())+
            labs(fill = 'CellType', color = 'CellType')+
            scale_x_continuous(breaks = NULL) +
            scale_y_continuous(breaks = NULL) +
            scale_fill_manual(values = color_pal) +
            scale_color_manual(values = color_pal) +
            # 手动绘制短坐标轴
            annotate("segment", x = -10, xend = -4, y = -10, yend = -10, color = "black", size = 0.8, arrow = arrow(type = "closed", length = unit(0.3, "cm"))) +
            annotate("segment", x = -10, xend = -10, y = -10, yend = -4, color = "black", size = 0.8, arrow = arrow(type = "closed", length = unit(0.3, "cm")))

pdf("/Volumes/GaoxyData/IAs_object/scRNA-seq/Fibro_CellType_umap.pdf", width = 8, height = 8)
print(plot)
dev.off()

# 提取表达数据
expr_data <- FetchData(F2_obj, vars = all_hub_gene)
expr_data$Cell <- rownames(expr_data)
expr_data$SampleType <- F2_obj@meta.data$type

# 长格式转换
expr_long <- expr_data %>%
  pivot_longer(cols = all_of(all_hub_gene), names_to = "Gene", values_to = "Expression")
expr_long$SampleType = factor(expr_long$SampleType, levels = c('sham', 'formed', 'ruptured'))

expr_long$Gene = factor(expr_long$Gene, levels = c('Col5a1', 'Creb3l1', 'Fbn1', 'Fkbp10', 'Inhba', 'Kdelr3', 'Ltbp2', 'Mfap5', 'Rcn3', 'Itga11', 'Adamts2', 'Col1a1', 'Fstl1', 'Prrx2'))

pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/Hub_gene_express_F2.pdf', width = 16, height = 7)
ggplot(expr_long, aes(x = SampleType, y = Expression, fill = SampleType)) +
  geom_boxplot(outlier.size = 0.3) +
  scale_fill_manual(values = c('#62D1B7', '#FFA42E', '#FF4F4F'))+
  facet_wrap(~ Gene, scales = "free_y", ncol = 7) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank()) +
  labs(x = "Sample Type", y = "Expression", title = "Hub Gene Expression in Fibroblasts Cluster2")
dev.off()


#%% 比例
# 自定义表达阈值，默认 0（表示非 0 就算表达）
threshold <- 0

# 创建一个用于存储比例的数据框
prop_df <- data.frame()

for (gene in all_hub_gene) {
  for (type_group in unique(Fibroblast_obj@meta.data$type)) {
    cells_in_group <- WhichCells(Fibroblast_obj, expression = type == type_group)
    gene_expr <- FetchData(Fibroblast_obj, vars = gene)[cells_in_group, gene]
    
    total_cells <- length(gene_expr)
    positive_cells <- sum(gene_expr > threshold)
    proportion <- positive_cells / total_cells
    
    prop_df <- rbind(prop_df, data.frame(
      Gene = gene,
      SampleType = type_group,
      Proportion = proportion
    ))
  }
}

# 设置分组顺序
prop_df$SampleType <- factor(prop_df$SampleType, levels = c("sham", "formed", "ruptured"))

# 绘图
ggplot(prop_df, aes(x = SampleType, y = Proportion, fill = SampleType)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ Gene, ncol = 3) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_bw() +
  labs(title = "Proportion of Cells Expressing Each Hub Gene (VSMC)",
       y = "Proportion of Expressing Cells",
       x = "Sample Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


#%% hub gene 集表达打分
# all_hub_gene = c(purple_up_mouse$mouseGene, pink_up_mouse$mouseGene)

GSE193533_obj_harmony = AddModuleScore(GSE193533_obj_harmony, features = list(all_hub_gene), name = "all_hub_gene")

FeaturePlot(GSE193533_obj_harmony, features = 'all_hub_gene1', combine = FALSE)

## ggplot2 绘图
umap = GSE193533_obj_harmony@reductions$umap@cell.embeddings %>%  #坐标信息
        as.data.frame() %>%
        cbind(cluster = GSE193533_obj_harmony@meta.data$seurat_clusters) # 注释后的label信息 ，改为cluster
umap$embedding = rownames(umap)
EnrichScore = as.data.frame(GSE193533_obj_harmony$all_hub_gene1)
EnrichScore$embedding = rownames(EnrichScore)
colnames(EnrichScore) = c("EnrichScore", "embedding")
EnrichScore[EnrichScore$EnrichScore>1, 'EnrichScore'] = 1
EnrichScore[EnrichScore$EnrichScore<0, 'EnrichScore'] = 0

EnrichScore_df = merge(umap, EnrichScore, by = 'embedding')

p = ggplot(EnrichScore_df,aes(x= UMAP_1 , y = UMAP_2 ,color = EnrichScore))+
        geom_point(size = 1 , alpha = 1 ) +
        scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")), limits = c(0, 1)) +
        theme_bw()+
        theme(panel.grid.major = element_blank(), #主网格线
              panel.grid.minor = element_blank(), #次网格线
              #panel.border = element_blank(), #边框
              #axis.title = element_blank(),  #轴标题
              axis.text = element_blank(), # 文本
              axis.ticks = element_blank(),
              panel.background = element_rect(fill = 'white'), #背景色
              plot.background = element_rect(fill = "white"),
              legend.title = element_blank(), #去掉legend.title
              legend.key = element_rect(fill = 'white'), #
              legend.text = element_text(size = 18), #设置legend标签的大小
              legend.key.size = unit(1, 'cm'))
ggsave(p, file = '/Volumes/GaoxyData/IAs_object/scRNA-seq/AddModuleScore.pdf', width = 7, height = 6)

markers = c('Col1a1', 'Fbn1', 'Ltbp2')

p_list = list()

for(i in 1:length(markers)){
    p = plot_density(GSE193533_obj_harmony, features = markers[i], pal = 'magma', raster = F, size = 0.8) +
            theme_blank() +
            theme(legend.frame = element_rect(colour = 'black'),
                  legend.ticks = element_line(colour = 'black', linewidth = 0),
                  legend.key.width = unit(0.3, 'cm'),
                  legend.key.height = unit(0.6, 'cm'),
                  legend.title = element_text(color = 'black', face = 'bold', size = 8),
                  plot.title = element_text(hjust = 0.5))
    p_list[[i]] = p
}

pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/HubGene-UMAP.pdf', width = 21, height = 7)
do.call(plot_grid, c(p_list, ncol = 3))
dev.off()


#%% 将亚型信息放回到原始数据中进行细胞通讯分析
# 复制 CellType 到新列
GSE193533_obj_harmony$CellChat_Group <- GSE193533_obj_harmony$CellType

# 替换 Fibroblast 为亚型
fibro_subtype_df <- data.frame(
  cell = colnames(Fibroblast_obj),
  subType = Fibroblast_obj$subType
)

# 替换原始 CellChat_Group 中对应细胞的值
GSE193533_obj_harmony$CellChat_Group[fibro_subtype_df$cell] <- fibro_subtype_df$subType

# 检查
table(GSE193533_obj_harmony$CellChat_Group)

saveRDS(GSE193533_obj_harmony, file = '/Volumes/GaoxyData/IAs_object/scRNA-seq/GSE193533_obj_harmony_recluster.rds')
#%% CellChat细胞通讯
library(CellChat)
library(patchwork)
library(Seurat)

GSE193533_obj_harmony_recluster = readRDS('/public3/Xinyu/IA/scRNA-seq/GSE193533_obj_harmony_recluster.rds')

# 创建 CellChat 对象
data.input <- GetAssayData(GSE193533_obj_harmony, assay = "RNA", slot = "data")  # 使用标准化数据
meta <- data.frame(labels = GSE193533_obj_harmony$CellChat_Group, row.names = colnames(GSE193533_obj_harmony))

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")

# 设置数据库（人类或小鼠）
CellChatDB <- CellChatDB.mouse  # 如果是小鼠
cellchat@DB <- CellChatDB

# 预处理
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

## 查看互作总数

## 查看两个亚型与其他细胞的相互作用强度

## 查看某个信号通路下的细胞通讯

