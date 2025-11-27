library(ggplot2)
library(dplyr)
library(tidyverse)
library(CellChat)
library(Seurat)

#%% 导入数据并提出亚集
GSE193533_obj_harmony = readRDS('/Volumes/GaoxyData/IAs_object/scRNA-seq/GSE193533_obj_harmony_recluster.rds')

F1_obj = subset(GSE193533_obj_harmony, subset = CellChat_Group %in% c("Fibroblast cluster1", "B lymphocyte", "Endothelium", "Erythroid lineage", "Macrophage", "Mast cells", "Neutrophil", "T lymphocyte", "VSMC"))
F2_obj = subset(GSE193533_obj_harmony, subset = CellChat_Group %in% c("Fibroblast cluster2", "B lymphocyte", "Endothelium", "Erythroid lineage", "Macrophage", "Mast cells", "Neutrophil", "T lymphocyte", "VSMC"))

Idents(F1_obj) = F1_obj$CellType
Idents(F2_obj) = F2_obj$CellType

#%% 创建CellChat对象
## F1
cellchat_F1 = createCellChat(F1_obj, group.by = 'CellType')
CellChatDB = CellChatDB.mouse
CellChatDB.use = CellChatDB
cellchat_F1@DB = CellChatDB.use
cellchat_F1 = subsetData(cellchat_F1)
cellchat_F1 = identifyOverExpressedGenes(cellchat_F1)
cellchat_F1 = identifyOverExpressedInteractions(cellchat_F1)
cellchat_F1 = computeCommunProb(cellchat_F1)
cellchat_F1 = aggregateNet(cellchat_F1)


groupSize = as.numeric(table(cellchat_F1@idents))
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat_F1@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = 'Number of interactions')
netVisual_circle(cellchat_F1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = 'Interaction weights/strength')

## F2
cellchat_F2 = createCellChat(F2_obj, group.by = 'CellType')
CellChatDB = CellChatDB.mouse
CellChatDB.use = CellChatDB
cellchat_F2@DB = CellChatDB.use
cellchat_F2 = subsetData(cellchat_F2)
cellchat_F2 = identifyOverExpressedGenes(cellchat_F2)
cellchat_F2 = identifyOverExpressedInteractions(cellchat_F2)
cellchat_F2 = computeCommunProb(cellchat_F2)
cellchat_F2 = aggregateNet(cellchat_F2)


groupSize = as.numeric(table(cellchat_F2@idents))
par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat_F2@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = 'Number of interactions')
netVisual_circle(cellchat_F2@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = 'Interaction weights/strength')


#%% 差异分析
# 获取目标统一分组
group.new <- union(levels(cellchat_F1@idents), levels(cellchat_F2@idents))
cellchat_F1 <- liftCellChat(cellchat_F1, group.new)
cellchat_F2 <- liftCellChat(cellchat_F2, group.new)

# 对两个对象统一 CellGroup 顺序
cellchat_F1 <- liftCellChat(cellchat_F1, group.new)
cellchat_F2 <- liftCellChat(cellchat_F2, group.new)

cellchat_F1 <- computeCommunProb(cellchat_F1)
cellchat_F1 <- computeCommunProbPathway(cellchat_F1)
cellchat_F1 <- aggregateNet(cellchat_F1)

cellchat_F2 <- computeCommunProb(cellchat_F2)
cellchat_F2 <- computeCommunProbPathway(cellchat_F2)
cellchat_F2 <- aggregateNet(cellchat_F2)


object.list = list(F1 = cellchat_F1, F2 = cellchat_F2)
cellchat.merge = mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/Number_interaction.pdf', width = 8, height = 5)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

#不同细胞群之间的相互作用数量或强度的差异
par(mfrow = c(1,2), xpd=TRUE)
pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/Differential_interaction.pdf', width = 8, height = 5)
netVisual_diffInteraction(cellchat.merge, comparison = c(1,2), weight.scale = T, color.edge = c("#F22613","#16B4F2"))
netVisual_diffInteraction(cellchat.merge, comparison = c(1,2), weight.scale = T, measure = "weight", color.edge = c("#F22613","#16B4F2"))
dev.off()





#%% 差异信号通路
cellchat_F1@idents <- droplevels(cellchat_F1@idents)
cellchat_F2@idents <- droplevels(cellchat_F2@idents)

cellchat_F1 = computeCommunProb(cellchat_F1)
cellchat_F1 = computeCommunProbPathway(cellchat_F1)
cellchat_F1 = netAnalysis_computeCentrality(cellchat_F1, slot.name = "netP") 


cellchat_F2 = computeCommunProb(cellchat_F2)
cellchat_F2 = computeCommunProbPathway(cellchat_F2)
cellchat_F2 = netAnalysis_computeCentrality(cellchat_F2, slot.name = "netP") 

group.new <- union(levels(cellchat_F1@idents), levels(cellchat_F2@idents))

# 对两个对象统一 CellGroup 顺序
cellchat_F1 <- liftCellChat(cellchat_F1, group.new)
cellchat_F2 <- liftCellChat(cellchat_F2, group.new)

object.list = list(F1 = cellchat_F1, F2 = cellchat_F2)
cellchat.merge2 = mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)

#计算和可视化通路距离
options(repr.plot.height = 12, repr.plot.width = 12)
rankSimilarity(cellchat.merge, type = "functional")
#比较每个信号通路的整体信息流
options(repr.plot.height = 10, repr.plot.width = 18)
gg1 = rankNet(cellchat.merge, mode = "comparison", stacked = T, do.stat = TRUE)+
  scale_fill_manual(values = c("#16B4F2","#F22613"))
gg2 = rankNet(cellchat.merge, mode = "comparison", stacked = F, do.stat = TRUE)
pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/DE_Signal_pathway.pdf', height = 17, width = 10)
gg1
dev.off()

ptm = Sys.time()
gg1 <- compareInteractions(cellchat.merge, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat.merge, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat.merge, weight.scale = T)
netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight")

g1 <- netVisual_heatmap(cellchat.merge)
#> Do heatmap based on a merged object
g2 <- netVisual_heatmap(cellchat.merge, measure = "weight")
#> Do heatmap based on a merged object
pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/DE_heatmap.pdf', height = 10, width = 20)
g1 + g2
dev.off()

gg <- lapply(names(object.list), function(nm) {
  netAnalysis_signalingRole_scatter(object.list[[nm]], title = nm)
})
patchwork::wrap_plots(gg)


cellchat_F1 <- netAnalysis_computeCentrality(cellchat_F1)
cellchat_F2 <- netAnalysis_computeCentrality(cellchat_F2)

object.list <- list(F1 = cellchat_F1, F2 = cellchat_F2)
cellchat.merge <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)

pg1 <- netAnalysis_signalingChanges_scatter(cellchat.merge, idents.use = "VSMC", signaling.exclude = "MIF")
pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/VSMC.pdf', height = 10, width = 12)
pg1
dev.off()

#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
pg2 <- netAnalysis_signalingChanges_scatter(cellchat.merge, idents.use = "Macrophage", signaling.exclude = c("MIF"))
pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/Macrophage.pdf', height = 10, width = 12)
pg2
dev.off()

#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/VSMC-Macro.pdf', height = 10, width = 24)
patchwork::wrap_plots(plots = list(pg1,pg2))
dev.off()
pg1 <- netAnalysis_signalingChanges_scatter(cellchat.merge, idents.use = "Firoblast", signaling.exclude = "MIF")
pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/Firoblast.pdf', height = 10, width = 12)
pg1
dev.off()

#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0
#> The following `from` values were not present in `x`: 0, -1
pg2 <- netAnalysis_signalingChanges_scatter(cellchat.merge, idents.use = "Endothelium", signaling.exclude = c("MIF"))
pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/Endothelium.pdf', height = 10, width = 12)
pg2
dev.off()

#> Visualizing differential outgoing and incoming signaling changes from NL to LS
#> The following `from` values were not present in `x`: 0, 2
#> The following `from` values were not present in `x`: 0, -1
pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/F-E.pdf', height = 10, width = 24)
patchwork::wrap_plots(plots = list(pg1,pg2))
dev.off()

#%% 绘图
## 气泡图
cellchat.merge@idents <- factor(cellchat.merge@meta$CellChat_Group)
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat.merge, group.dataset = "datasets", pos.dataset = 'F1', features.name = 'Fibroblast cluster2', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = 'Fibroblast')
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "F1",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "F2",ligand.logFC = -0.1, receptor.logFC = -0.1)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat.merge2)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat.merge2)


pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat.merge2, pairLR.use = pairLR.use.up, sources.use = c(3,4), targets.use = c(7:9), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat.merge2, pairLR.use = pairLR.use.down, sources.use = c(3,4), targets.use = c(7:9), comparison = c(1, 2),  angle.x = 45, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[1]))
#> Comparing communications on a merged object
pdf('D:/DE_signal_bubble.pdf', width = 12, height = 9)
gg1 + gg2
dev.off()


## 单个通路的弦图
pathways.show <- c("SPP1") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/SPP1_signal.pdf', width = 8, height = 5)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()


pathways.show <- c("SEMA3") 

pdf('/Volumes/GaoxyData/IAs_object/scRNA-seq/F1_SEMA3_signal.pdf', width = 8, height = 5)
netVisual_aggregate(cellchat_F1, signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
dev.off()


