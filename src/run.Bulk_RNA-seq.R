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

#%% 自定义函数
## pheatmap使用的是grid图形体系，不是ggplot2体系
save_pheatmap_pdf = function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


#%% 读取数据
GSE122897_raw = read.table(
  "/Volumes/GaoxyData/IAs_object/bulkRNA-seq/GSE122897_readCounts_raw.txt.gz",
  header = TRUE,           # 第一行为列名
  row.names = 1,           # 第一列是基因名
  sep = "\t",              # 制表符分隔
  check.names = FALSE      # 防止样本名被自动修改
)

# ID转换
mapped = mapIds(org.Hs.eg.db,
                keys = rownames(GSE122897_raw),
                column = "SYMBOL",
                keytype = "ENSEMBL",
                multiVals = "first") %>% as.data.frame() %>% na.omit()
colnames(mapped) = 'Symbol'

# 去重
mapped_unique = mapped %>% distinct(Symbol, .keep_all = TRUE)

com_ensemble = intersect(rownames(GSE122897_raw), rownames(mapped_unique))

GSE122897_raw_renamed = GSE122897_raw[com_ensemble, ]
rownames(GSE122897_raw_renamed) = mapped_unique[com_ensemble, 'Symbol']

#%% DESeq2筛选差异基因
GSE122897_filter = GSE122897_raw_renamed[rowSums(GSE122897_raw_renamed) > 100, ]

# 假设前16个样本为对照组 CO，后44个样本为实验组 IA
sample_group = rep(c("CO", "IA"), times = c(16, 44))
coldata = data.frame(condition = factor(sample_group))
rownames(coldata) = colnames(GSE122897_filter)

# 构建 DESeqDataSet
dds = DESeqDataSetFromMatrix(countData = GSE122897_filter,
                              colData = coldata,
                              design = ~ condition)

# 筛选差异基因
dds = DESeq(dds)
res = results(dds, contrast = c("condition", "IA", "CO"))

# 结果筛选
# 去除NA
res = na.omit(res)

# 按 padj 排序
res_ordered = res[order(res$padj), ]

# 筛选显著差异基因
res_sig = subset(res_ordered, padj < 0.05 & abs(log2FoldChange) > 1)

# 添加标签
res_ordered$express = 'None'
res_ordered[which(res_ordered$padj < 0.05 & res_ordered$log2FoldChange > 1), 'express'] = 'Up'
res_ordered[which(res_ordered$padj < 0.05 & res_ordered$log2FoldChange < -1), 'express'] = 'Down'

#%% PCA分析
vsd = vst(dds, blind = FALSE) 

# 选取方差最大的前1000个基因进行 PCA
top_var_genes = head(order(apply(assay(vsd), 1, var), decreasing = TRUE), 1000)
pca_res = prcomp(t(assay(vsd)[top_var_genes, ]), scale. = TRUE)

# 提取 PCA 结果
pca_df = as.data.frame(pca_res$x)

# 添加样本分组信息
pca_df$Sample = rownames(pca_df)
pca_df$Group = colData(vsd)$condition  # 请根据你的实验设计调整，如 "group"、"condition" 等

# 计算解释变异的比例
percentVar = round(100 * summary(pca_res)$importance[2, 1:2], 1)

# 画图
p = ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
        geom_point(size = 5) +
        scale_color_manual(values = c('#FF008C', '#63D3FF'))+
        xlab(paste0("PC1 (", percentVar[1], "%)")) +
        ylab(paste0("PC2 (", percentVar[2], "%)")) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          text = element_text(size = 12),
          panel.grid = element_blank(),
          legend.position = "right",
          plot.title = element_text(hjust = 0.5)
        ) +
        ggtitle("PCA")

ggsave(p, file = '/Volumes/GaoxyData/IAs_object/bulkRNA-seq/PCA.pdf', width = 7, height = 6)


#%% 火山图
res_df = as.data.frame(res_ordered)

# 上调前5
top5_up_gene = res_df %>% 
  filter(log2FoldChange > 0, padj < 0.05) %>%
  arrange(desc(log2FoldChange)) %>%
  slice(1:5)

# 下调前5
top5_down_gene = res_df %>% 
  filter(log2FoldChange < 0, padj < 0.05) %>%
  arrange(log2FoldChange) %>%
  slice(1:5)

# 合并前10基因
top10_gene = bind_rows(top5_up_gene, top5_down_gene)
top10_gene$gene = rownames(top10_gene)
# 绘图
pdf("/Volumes/GaoxyData/IAs_object/bulkRNA-seq/Volcano.pdf", width = 5, height = 5)
ggplot(res_ordered, aes(x = log2FoldChange, y = -log10(padj))) + 
  geom_point(aes(colour = express), size = 1, alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), lty = 4, lwd = 0.5, col = "#A1A1A6") +
  geom_hline(yintercept = -log10(0.05), lty = 4, lwd = 0.5, col = "#A1A1A6") +
  scale_color_manual(values = c("#277FF2", "#D9D8D7", "#FF9C08")) +
  geom_text_repel(
    data = top10_gene,
    aes(x = log2FoldChange, y = -log10(padj), label = gene),
    size = 3,
    box.padding = 0.3,
    max.overlaps = Inf,
    segment.color = "grey50"
  ) +
  labs(x = "log2(FoldChange)", y = "-log10(pvalue.adjust)") +
  xlim(-6.5, 6.5) +
  ylim(0, 10) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.title = element_text(face = "bold"),
        axis.text = element_text(colour = "black"))
dev.off()

#%% 绘制热图
# DESeq2结果转为data.frame
res_df = as.data.frame(res_ordered)

# 添加基因名列
res_df$gene = rownames(res_df)

write.csv(res_df, file = '/Volumes/GaoxyData/IAs_object/bulkRNA-seq/DESeq2_DEG.csv')

# 选出差异表达基因
DEG = res_df %>% filter(padj < 0.05, abs(log2FoldChange) > 1)

# 按 log2FoldChange 升序和降序分别取下调和上调基因
top_up_genes = res_df %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(log2FoldChange)) %>%
  slice_head(n = 10) %>%
  pull(gene)

top_down_genes = res_df %>%
  filter(log2FoldChange < 0) %>%
  arrange(log2FoldChange) %>%
  slice_head(n = 10) %>%
  pull(gene)

# 合并为需要展示行名的20个基因
top20_genes = c(top_up_genes, top_down_genes)
# 使用rlog或vst处理后的表达数据

rld = rlog(dds, blind = TRUE)
expr_matrix = assay(rld)  # 或 assay(vsd)

# 提取DEG的表达数据
DEG_expr = expr_matrix[rownames(expr_matrix) %in% DEG$gene, ]

# 创建一个行名向量，仅保留这20个基因，其它设置为空字符串
labels_row = ifelse(rownames(DEG_expr) %in% top20_genes, rownames(DEG_expr), "")


# 可选：对样本进行聚类
p = pheatmap(DEG_expr,
             scale = "row", 
             show_rownames = FALSE,
             fontsize_row = 5,
             cluster_cols = FALSE,
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             clustering_method = "complete",
             color = colorRampPalette(c("navy", "white", "firebrick3"))(50))

save_pheatmap_pdf(p, filename = "/Volumes/GaoxyData/IAs_object/bulkRNA-seq/DEG_heatmap.pdf", width = 10, height = 8)

#%% GO 和 KEGG富集分析
# GO 富集分析
df_sig = res_df[which(res_df$express != 'None'), c('express', 'gene')]
ALL_ego = compareCluster(
            gene~express, 
            data = df_sig, 
            keyType = 'SYMBOL',
            fun = "enrichGO", 
            OrgDb = "org.Hs.eg.db",
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.05
            )
# 去除冗余
ALL_ego_sim = simplify(ALL_ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
dotplot(ALL_ego_sim, showCategory=10, font.size = 8)

# 提取结果为数据框
df = as.data.frame(ALL_ego_sim)

# 保留上调 & 下调中 p.adjust 最小的 top 5
top_up = df %>%
  filter(Cluster == "Up") %>%
  slice_min(order_by = p.adjust, n = 5)

top_down = df %>%
  filter(Cluster == "Down") %>%
  slice_min(order_by = p.adjust, n = 5)

# 合并上下调
top_terms = bind_rows(top_up, top_down)

top_terms$Description = sapply(top_terms$Description, function(x) paste(strwrap(x, width = 40), collapse = "\n"))

# 计算 log10(p.adjust)，上调为正，下调为负
top_terms = top_terms %>%
  mutate(
    logP = -log10(p.adjust),
    logP = ifelse(Cluster == "Down", -logP, logP),  # Down向左
    Description = fct_reorder(Description, logP)    # 排序
  )

# 绘图
p = ggplot(top_terms) +
      geom_segment(aes(x = Description, xend = Description, y = 0, yend = logP), color = "#d3d3d3") +
      geom_point(aes(x = Description, y = logP, size = abs(logP), color = Cluster)) +
      geom_hline(yintercept = 0, color = "gray30", linewidth = 0.5) +
      scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5")) +
      scale_size_continuous(range = c(3, 10), breaks = c(1, 3, 5, 10, 20),name = expression("-log"[10]*"(adj. P)")) +
      coord_flip() +
      labs(x = "", y = expression("-log"[10]*"(adj. P value)"),title = "Top 5 Enriched GO BP Terms") +
      theme_classic() +
      theme(axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title.position = "plot",
            plot.title = element_text(size = 20, hjust = 0.5),
            legend.text = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14)
      )
ggsave(p, file = '/Volumes/GaoxyData/IAs_object/bulkRNA-seq/DEG_GO_plot.pdf', height = 8, width = 10)

# KEGG 富集分析
# SYMBOL 转 ENTREZID
df_entrez = bitr(df_sig$gene,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

# 合并表达类型
df_sig = merge(df_sig, df_entrez, by.x = "gene", by.y = "SYMBOL")

# KEGG 富集分析（使用 ENTREZID，不指定 keyType）
ALL_ekegg = compareCluster(
  ENTREZID ~ express, 
  data = df_sig, 
  fun = "enrichKEGG", 
  organism = "hsa",        # 指定物种
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)

# 画图
dotplot(ALL_ekegg, showCategory = 10, font.size = 8) + 
  ggtitle("KEGG Pathway Enrichment")

df = as.data.frame(ALL_ekegg)

# 保留上调 & 下调中 p.adjust 最小的 top 5
top_up = df %>%
  filter(Cluster == "Up") %>%
  slice_min(order_by = p.adjust, n = 5)

top_down = df %>%
  filter(Cluster == "Down") %>%
  slice_min(order_by = p.adjust, n = 5)

# 合并上下调
top_kegg = bind_rows(top_up, top_down)

top_kegg$Description = sapply(top_kegg$Description, function(x) paste(strwrap(x, width = 40), collapse = "\n"))

# 计算 log10(p.adjust)，上调为正，下调为负
top_kegg = top_kegg %>%
  mutate(
    logP = -log10(p.adjust),
    logP = ifelse(Cluster == "Down", -logP, logP),  # Down向左
    Description = fct_reorder(Description, logP)    # 排序
  )

p = ggplot(top_kegg) +
      geom_segment(aes(x = Description, xend = Description, y = 0, yend = logP), color = "#d3d3d3") +
      geom_point(aes(x = Description, y = logP, size = abs(logP), color = Cluster)) +
      geom_hline(yintercept = 0, color = "gray30", linewidth = 0.5) +
      scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5")) +
      scale_size_continuous(range = c(3, 10), breaks = c(1, 3, 5, 10, 20),name = expression("-log"[10]*"(adj. P)")) +
      coord_flip() +
      labs(x = "", y = expression("-log"[10]*"(adj. P value)"),title = "Top 5 Enriched KEGG pathways") +
      theme_classic() +
      theme(axis.title.y = element_blank(),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title.position = "plot",
            plot.title = element_text(size = 20, hjust = 0.5),
            legend.text = element_text(size = 12),
            axis.text.y = element_text(size = 12),
            axis.text.x = element_text(size = 12),
            axis.title = element_text(size = 14),
            legend.title = element_text(size = 14)
      )
ggsave(p, file = '/Volumes/GaoxyData/IAs_object/bulkRNA-seq/DEG_KEGG_plot.pdf', height = 8, width = 10)

#%% 进行WGCNA分析
# 转置为行是样本，列是基因
expr_t = t(expr_matrix) %>% as.data.frame()

# 检查缺失值和离群值
gsg = goodSamplesGenes(expr_t, verbose = 3);
gsg$allOK

# 样本进行聚类，检查缺失值和离群值
sampleTree = hclust(dist(expr_t), method = "average")
sizeGrWindow(12,9) #视图
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# 删除离群样本
abline(h = 160, col = "red") #划定需要剪切的枝长
clust = cutreeStatic(sampleTree, cutHeight = 160, minSize = 10)

table(clust)   
keepSamples = (clust==1)  #保留非离群(clust==1)的样本
datExpr = expr_t[keepSamples, ]  #去除离群值后的数据
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# 选择软阈值
powers = c(1:20)
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# 可视化
pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/SoftThreshold.pdf', width = 8, height = 8)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R²",
     main = 'Scaleindependence',
     type="n")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=0.9, col="red")
abline(h=0.9, col="red")
dev.off()


# 平均连接度
pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/MeanConnectivity.pdf', width = 8, height = 8)
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
dev.off()

# 软阈值选择
sft$powerEstimate

# 构建共表达网络并识别模块
net = blockwiseModules(datExpr, power = 7,
                        TOMType = "signed", minModuleSize = 30,
                        reassignThreshold = 0,
                        mergeCutHeight = 0.25,
                        numericLabels = TRUE,
                        pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        saveTOMFileBase = "TOM",
                        verbose = 3)


# 查看划分的模块数和每个模块里面包含的基因个数
table(net$colors)

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)

pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/ClusterDendrogram.pdf', height = 6, width = 10)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


# 保存分配模块和模块包含的基因信息。
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree, file = "/Volumes/GaoxyData/IAs_object/bulkRNA-seq/WGCNA.1.RData")

#%% 模块表型关联分析
# 表型数据获取
library(GEOquery)
gse = getGEO("GSE122897", GSEMatrix = TRUE)[[1]]
pheno = pData(gse)
#head(pheno)

pheno_df = pheno[, c('title', 'characteristics_ch1', 'characteristics_ch1.1', 'source_name_ch1')]
pheno_df$group = c(rep('CO', 16), rep('IA', 44))
rownames(pheno_df) = pheno_df$title

## 数值化
traitData = pheno_df %>%
  mutate(
    Sex = ifelse(grepl("Male", characteristics_ch1), 1, 0),
    Age = as.numeric(gsub("age: ", "", characteristics_ch1.1)),
    Group = ifelse(group == "IA", 1, 0)
  ) %>%
  dplyr::select(Sex, Age, Group)

rownames(traitData) = rownames(pheno_df)
datTraits = traitData[rownames(datExpr), ]

sampleTree2 = hclust(dist(datExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度
# 先把 source_name_ch1 转成因子，然后映射颜色
source_factor = factor(pheno_df$source_name_ch1)

my_source_colors = c(
  "Intracranial cortical artery" = "#B5F587",
  "Ruptured intracranial aneurysm" = "#F09294",
  "Unknown intracranial aneurysm" = "#F2F080",
  "Unknown rupture status" = "#90CBF0",
  "Unruptured intracranial aneurysm" = "#D2A4DB"
)
source_colors = my_source_colors[as.character(source_factor)]

# traitColors 是之前对数值型表型的颜色编码，组合起来
combined_colors = cbind(traitColors, Source = source_colors)
color_legend = unique(data.frame(color=source_colors, label=levels(source_factor)))

pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/WGCNA-Module-Trait.pdf', width = 12, height = 5)
# 第一块画主图
plotDendroAndColors(sampleTree2, 
                    combined_colors,
                    groupLabels = c(names(datTraits), "Type"),
                    main = "Sample dendrogram and trait heatmap with source_name_ch1")
dev.off()

pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/WGCNA-Module-Trait-legend.pdf', width = 12, height = 6)
plot.new()
legend("center", legend=levels(source_factor),
       fill=unique(source_colors),
       border=NA, bty="n", cex=1)
dev.off()


# 计算模块与样本特征的相关性
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# 重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 通过相关值对每个关联进行颜色编码
sizeGrWindow(10,6)
# 展示模块与表型数据的相关系数和 P值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

# 用热图的形式展示相关系数
pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/WGCNA-Module-Trait-Correlation.pdf', width = 5, height = 10)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               xLabelsAngle = 0,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


#%% 可视化加权网络
#计算TOM矩阵
TOM = TOMsimilarityFromExpr(datExpr, power = 7);
dissTOM = 1 - TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
plotTOM = dissTOM^7;
diag(plotTOM) = NA;
sizeGrWindow(9,9)
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

# 选取400个基因进行展示
nSelect = 200
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
sizeGrWindow(9,9)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
# TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

#改变热图的深色背景为白色背景：
library(gplots) # 需要先安装这个包才能加载。
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')

pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/WGCNA-Module-Trait-heatmap.pdf', width = 10, height = 10)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes", col = myheatcol)
dev.off()

#%% 可视化
# 重新计算模块的eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# 提取体重的表型数据
Condition = as.data.frame(datTraits$Group);
names(Condition) = "Condition"
# 加入到相应的模块
MET = orderMEs(cbind(MEs, Condition))
#画图
sizeGrWindow(6,6);
pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/WGCNA-engenier-Condition.pdf', width = 10, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/WGCNA-engenier-Condition-heatmap.pdf', width = 7, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 0)
dev.off()

# 重新计算模块的eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
# 提取体重的表型数据
Age = as.data.frame(datTraits$Age);
names(Age) = "Age"
# 加入到相应的模块
MET = orderMEs(cbind(MEs, Age))
#画图
sizeGrWindow(6,6);
pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/WGCNA-engenier-Age.pdf', width = 10, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()

pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/WGCNA-engenier-Age-heatmap.pdf', width = 7, height = 6)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 0)
dev.off()

#%% 保存网络列表
# 把连接矩阵转为边列表（导出给 Cytoscape）
# 选择棕色和红色的模块
modules = c("sienna3", "purple", "orange", "skyblue", "saddlebrown", 'darkred', 'royalblue', 'pink', 'violet', 'darkolivegreen', 'midnightblue');
probes = names(datExpr)
inModule = is.finite(match(moduleColors, modules));
modGenes = probes[inModule];

# 选择相关的 TOM矩阵
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modGenes, modGenes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("/Volumes/GaoxyData/IAs_object/bulkRNA-seq/", "CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("/Volumes/GaoxyData/IAs_object/bulkRNA-seq/", "CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])

#%% 每个模块的GO 
# probes: 基因名（或ENSEMBL ID，或者probe ID）
# moduleColors: 每个基因对应的模块颜色（比如 "blue", "turquoise"...）
modules = c("sienna3", "purple", "orange", "skyblue", "saddlebrown", 'darkred', 'royalblue', 'pink', 'violet', 'darkolivegreen', 'midnightblue')
probes = colnames(datExpr)  # 或 rownames(datExpr)

# 创建一个列表，每个模块对应一组基因
module_genes_list = lapply(modules, function(mod) {
  probes[moduleColors == mod]
})
names(module_genes_list) = modules


ego_list = lapply(names(module_genes_list), function(mod) {
  genes = module_genes_list[[mod]]
  enrichGO(gene         = genes,
           OrgDb        = org.Hs.eg.db,
           keyType      = "SYMBOL",
           ont          = "BP",
           pAdjustMethod= "BH",
           qvalueCutoff = 0.05,
           readable     = TRUE)
})
names(ego_list) = names(module_genes_list)

# 提取每个模块 top5 GO 结果，加上模块名
top_terms_per_module = lapply(names(ego_list), function(mod) {
  ego = ego_list[[mod]]
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  
  df = as.data.frame(ego) %>%
    slice_min(p.adjust, n = 5) %>%
    mutate(Module = mod)
})

# 合并所有模块
top_enrich_df = do.call(rbind, top_terms_per_module)

top_enrich_df$GeneRatio = sapply(top_enrich_df$GeneRatio, function(x) {
  eval(parse(text = x))
})
top_enrich_df$Module = factor(top_enrich_df$Module, levels = modules)

pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/WGCNA-module-GO.pdf', height = 8, width = 8)
ggplot(top_enrich_df, aes(x = Module, 
                          y = Description, 
                          size = GeneRatio, 
                          color = -log10(p.adjust))) +
  geom_point() +
  scale_color_gradient(low = "#4876FF", high = "violetred1") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "GO BP Terms",
       x = "Module", y = "GO Term",
       color = "-log10(adj. P value)",
       size = "Gene Ratio")
dev.off()

#%% 每个模块的KEGG
# 1. 获取模块及基因列表
modules = c("sienna3", "purple", "orange", "skyblue", "saddlebrown", 'darkred', 'royalblue', 'pink', 'violet', 'darkolivegreen', 'midnightblue')
probes = colnames(datExpr)  # 或 rownames(datExpr)

module_genes_list = lapply(modules, function(mod) {
  probes[moduleColors == mod]
})
names(module_genes_list) = modules

# 2. SYMBOL ➜ ENTREZID 转换
gene_df = bitr(unique(unlist(module_genes_list)), 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db)

# 3. 构建 ENTREZID 列表 per module
module_entrez_list = lapply(module_genes_list, function(genes) {
  gene_df$ENTREZID[match(genes, gene_df$SYMBOL)] %>% na.omit() %>% unique()
})

# 4. KEGG 富集分析
ekegg_list = lapply(module_entrez_list, function(genes) {
  enrichKEGG(gene         = genes,
             organism     = 'hsa',
             pvalueCutoff = 0.05)
})
names(ekegg_list) = names(module_entrez_list)

# 5. 提取每个模块的 top5 KEGG 通路
top_terms_per_module = lapply(names(ekegg_list), function(mod) {
  ego = ekegg_list[[mod]]
  if (is.null(ego) || nrow(ego@result) == 0) return(NULL)
  
  df = as.data.frame(ego) %>%
    slice_min(p.adjust, n = 5) %>%
    mutate(Module = mod)
})
top_enrich_df = do.call(rbind, top_terms_per_module)

# 6. 转换 GeneRatio 为数值
top_enrich_df$GeneRatio = sapply(top_enrich_df$GeneRatio, function(x) eval(parse(text = x)))
top_enrich_df$Module = factor(top_enrich_df$Module, levels = modules)

# 7. 可视化 KEGG dotplot
pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/WGCNA-module-KEGG.pdf', height = 8, width = 8)
ggplot(top_enrich_df, aes(x = Module, y = Description, size = GeneRatio, color = -log10(p.adjust))) +
  geom_point() +
  scale_color_gradient(low = "#4876FF", high = "violetred1") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(title = "KEGG Pathways",
       x = "Module", y = "KEGG Pathway",
       color = "-log10(adj. P value)",
       size = "Gene Ratio")
dev.off()

#%% DEG和WGCNA交集
modules = c("sienna3", "purple", "orange", "skyblue", "saddlebrown", 'darkred', 'royalblue', 'pink', 'violet', 'darkolivegreen', 'midnightblue')
probes = colnames(datExpr)  # 或 rownames(datExpr)

module_genes_list = lapply(modules, function(mod) { probes[moduleColors == mod] })
names(module_genes_list) = modules

# 计算上调交集
up_gene_list = lapply(module_genes_list, function(genes) {
  intersect(genes, res_df[res_df$express == 'Up', 'gene'])
})

# 计算下调交集
down_gene_list = lapply(module_genes_list, function(genes) {
  intersect(genes, res_df[res_df$express == 'Down', 'gene'])
})

# 统计每个模块中上调/下调基因数
up_counts = sapply(up_gene_list, length)
down_counts = sapply(down_gene_list, length)

# 设置下调为负值，使其在图中出现在下方
df_counts_sym = data.frame(Module = names(up_counts), Up = up_counts, Down = -down_counts)

# 转换为长格式
df_counts_long = pivot_longer(df_counts_sym, cols = c("Up", "Down"), names_to = "Direction", values_to = "Count")

# 绘图
p = ggplot(df_counts_long, aes(x = Module, y = Count, fill = Direction)) +
      geom_bar(stat = "identity", width = 0.6) +
      scale_fill_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5")) +
      geom_hline(yintercept = 0, color = "black") +
      labs(title = "Barplot of Up/Down-regulated Genes per Module", x = "Module", y = "Gene Count") +
      coord_flip() +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            axis.text.y = element_text(size = 10)
      )
ggsave(p, file = '/Volumes/GaoxyData/IAs_object/bulkRNA-seq/Module-DEG-stat.pdf', height = 8, width = 8)

#%% net网络文件构建
net_node = read.table('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/CytoscapeInput-nodes-sienna3-purple-orange-skyblue-saddlebrown-darkred-royalblue-pink-violet-darkolivegreen-midnightblue.txt', sep = '\t', header = TRUE)
net_edge = read.table('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/CytoscapeInput-edges-sienna3-purple-orange-skyblue-saddlebrown-darkred-royalblue-pink-violet-darkolivegreen-midnightblue.txt', sep = '\t', header = TRUE)

purple_up_net_edge = net_edge %>% filter(fromNode %in% up_gene_list[['purple']] & toNode %in% up_gene_list[['purple']])
purple_up_net_edge_unique <- purple_up_net_edge %>%
                                mutate(edge_pair = ifelse(fromNode < toNode,
                                                          paste(fromNode, toNode, sep = "_"),
                                                          paste(toNode, fromNode, sep = "_"))) %>%
                                distinct(edge_pair, .keep_all = TRUE) %>%
                                dplyr::select(-edge_pair)

pink_up_net_edge = net_edge %>% filter(fromNode %in% up_gene_list[['pink']] & toNode %in% up_gene_list[['pink']])
pink_up_net_edge_unique <- pink_up_net_edge %>%
                                mutate(edge_pair = ifelse(fromNode < toNode,
                                                          paste(fromNode, toNode, sep = "_"),
                                                          paste(toNode, fromNode, sep = "_"))) %>%
                                distinct(edge_pair, .keep_all = TRUE) %>%
                                dplyr::select(-edge_pair)
write.table(purple_up_net_edge_unique, file = '/Volumes/GaoxyData/IAs_object/bulkRNA-seq/purple_up_net_edge.txt', sep = '\t', col.names = T, row.names = F, quote = FALSE)
write.table(pink_up_net_edge_unique, file = '/Volumes/GaoxyData/IAs_object/bulkRNA-seq/pink_up_net_edge.txt', sep = '\t', col.names = T, row.names = F, quote = FALSE)



purple_down_net_edge = net_edge %>% filter(fromNode %in% down_gene_list[['purple']] & toNode %in% down_gene_list[['purple']])
purple_down_net_edge_unique <- purple_down_net_edge %>%
                                mutate(edge_pair = ifelse(fromNode < toNode,
                                                          paste(fromNode, toNode, sep = "_"),
                                                          paste(toNode, fromNode, sep = "_"))) %>%
                                distinct(edge_pair, .keep_all = TRUE) %>%
                                dplyr::select(-edge_pair)

pink_down_net_edge = net_edge %>% filter(fromNode %in% down_gene_list[['pink']] & toNode %in% down_gene_list[['pink']])
pink_down_net_edge_unique <- pink_down_net_edge %>%
                                mutate(edge_pair = ifelse(fromNode < toNode,
                                                          paste(fromNode, toNode, sep = "_"),
                                                          paste(toNode, fromNode, sep = "_"))) %>%
                                distinct(edge_pair, .keep_all = TRUE) %>%
                                dplyr::select(-edge_pair)
write.table(purple_down_net_edge_unique, file = '/Volumes/GaoxyData/IAs_object/bulkRNA-seq/purple_down_net_edge.txt', sep = '\t', col.names = T, row.names = F, quote = FALSE)
write.table(pink_down_net_edge_unique, file = '/Volumes/GaoxyData/IAs_object/bulkRNA-seq/pink_down_net_edge.txt', sep = '\t', col.names = T, row.names = F, quote = FALSE)

#%% 富集分析
purple_up_hub_gene = c('COL1A1', 'COL5A1', 'CHRF', 'FSTL1', 'KDELR3', 'FBN1', 'ADAMTS2', 'RCN3', 'FKBP10', 'CREB3L1')
pink_up_hub_gene = c('KRT8', 'KRT18', 'MFAP5', 'KRT7', 'INHBA', 'LTBP2', 'SPATA18', 'PRRX2', 'PODNL1', 'ITGA11')
purple_down_hub_gene = down_gene_list[['purple']]
pink_down_hub_gene = down_gene_list[['pink']]

df_go = data.frame(gene = c(purple_up_hub_gene,pink_up_hub_gene,purple_down_hub_gene,pink_down_hub_gene),
                   group = c(rep('Purple hub up gene', 10),
                             rep('Pink hub up gene', 10),
                             rep('Purple down gene', length(down_gene_list[['purple']])),
                             rep('Pink down gene', length(down_gene_list[['pink']]))))

hub_ego = compareCluster(
            gene~group, 
            data = df_go, 
            keyType = 'SYMBOL',
            fun = "enrichGO", 
            OrgDb = "org.Hs.eg.db",
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.05
            )
pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/Hub_go.pdf', height = 8, width = 7)
dotplot(hub_ego, showCategory = 10, font.size = 10, title = "GO Enrichment of Hub Genes")+
theme(plot.title = element_text(hjust = 0.5),
axis.title.x = element_blank(),
axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))
dev.off()


df_entrez = bitr(df_go$gene,
                 fromType = "SYMBOL",
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)

# 合并表达类型
df_sig = merge(df_go, df_entrez, by.x = "gene", by.y = "SYMBOL")

# KEGG 富集分析（使用 ENTREZID，不指定 keyType）
ALL_ekegg = compareCluster(
  ENTREZID ~ group, 
  data = df_sig, 
  fun = "enrichKEGG", 
  organism = "hsa",        # 指定物种
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/Hub_kegg.pdf', height = 8, width = 7)
dotplot(ALL_ekegg, showCategory = 10, font.size = 10, title = "GO Enrichment of Hub Genes")+
theme(plot.title = element_text(hjust = 0.5),
axis.title.x = element_blank(),
axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.5))
dev.off()

## 弦图绘制
purple_hub_up_FC = res_df[which(res_df$gene %in% purple_up_hub_gene), c('log2FoldChange', 'gene')]
purple_hub_up_term = hub_ego@compareClusterResult[hub_ego@compareClusterResult$Cluster == 'Purple hub up gene', c('Description', 'geneID')][1:5, ]
purple_hub_up_Clean = purple_hub_up_term %>% 
  separate_rows(geneID, sep = "/") %>%
  left_join(purple_hub_up_FC,by = c("geneID" = "gene"))
colnames(purple_hub_up_Clean) = c('Term', 'Genes', 'logFC')

purple_hub_up_Plot = chord_dat(data.frame(purple_hub_up_Clean),  process = unique(purple_hub_up_Clean$Term))

pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/purple_hub_up_chord.pdf', width = 7, height = 8)
GOChord(purple_hub_up_Plot,
        title = "GOChord plot",                          # 标题名称
        space = 0.02,                                    # 基因方格之间的空隙
        gene.order = "logFC",                            # 基因的排序方式，可以按照"logFC", "alphabetical", "none", 
        gene.space = 0.25,                               # 基因标签离图案的距离
        gene.size = 5,                                   # 基因标签的大小
        lfc.col = c('firebrick3', 'white','royalblue3'), # logFC图例的颜色
        ribbon.col = brewer.pal(ncol(dfPlot)-1, "Set2"), # 条带颜色设置
        process.label = 8                                # 图例标签的大小
)+
theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
dev.off()



pink_hub_up_FC = res_df[which(res_df$gene %in% pink_up_hub_gene), c('log2FoldChange', 'gene')]
pink_hub_up_term = hub_ego@compareClusterResult[hub_ego@compareClusterResult$Cluster == 'Pink hub up gene', c('Description', 'geneID')][1:5, ]
pink_hub_up_Clean = pink_hub_up_term %>% 
  separate_rows(geneID, sep = "/") %>%
  left_join(pink_hub_up_FC,by = c("geneID" = "gene"))
colnames(pink_hub_up_Clean) = c('Term', 'Genes', 'logFC')

pink_hub_up_Plot = chord_dat(data.frame(pink_hub_up_Clean),  process = unique(pink_hub_up_Clean$Term))

pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/pink_hub_up_chord.pdf', width = 10, height = 11)
GOChord(pink_hub_up_Plot,
        title = "GOChord plot",                          # 标题名称
        space = 0.02,                                    # 基因方格之间的空隙
        gene.order = "logFC",                            # 基因的排序方式，可以按照"logFC", "alphabetical", "none", 
        gene.space = 0.25,                               # 基因标签离图案的距离
        gene.size = 5,                                   # 基因标签的大小
        lfc.col = c('firebrick3', 'white','royalblue3'), # logFC图例的颜色
        ribbon.col = brewer.pal(ncol(dfPlot)-1, "Set2"), # 条带颜色设置
        process.label = 8                                # 图例标签的大小
)+
theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
dev.off()

pink_hub_down_FC = res_df[which(res_df$gene %in% pink_down_hub_gene), c('log2FoldChange', 'gene')]
pink_hub_down_term = hub_ego@compareClusterResult[hub_ego@compareClusterResult$Cluster == 'Pink down gene', c('Description', 'geneID')][1:5, ]
pink_hub_down_Clean = pink_hub_down_term %>% 
  separate_rows(geneID, sep = "/") %>%
  left_join(pink_hub_down_FC,by = c("geneID" = "gene"))
colnames(pink_hub_down_Clean) = c('Term', 'Genes', 'logFC')

pink_hub_down_Plot = chord_dat(data.frame(pink_hub_down_Clean),  process = unique(pink_hub_down_Clean$Term))

pdf('/Volumes/GaoxyData/IAs_object/bulkRNA-seq/pink_hub_down_chord.pdf', width = 10, height = 11)
GOChord(pink_hub_down_Plot,
        title = "GOChord plot",                          # 标题名称
        space = 0.02,                                    # 基因方格之间的空隙
        gene.order = "logFC",                            # 基因的排序方式，可以按照"logFC", "alphabetical", "none", 
        gene.space = 0.25,                               # 基因标签离图案的距离
        gene.size = 5,                                   # 基因标签的大小
        lfc.col = c('firebrick3', 'white','royalblue3'), # logFC图例的颜色
        ribbon.col = brewer.pal(ncol(dfPlot)-1, "Set2"), # 条带颜色设置
        process.label = 8                                # 图例标签的大小
)+
theme(plot.title = element_text(hjust = 0.5, face = 'bold'))
dev.off()

#%% 保存环境
save.image(file = '/Volumes/GaoxyData/IAs_object/bulkRNA-seq/bulk_image.RData')
