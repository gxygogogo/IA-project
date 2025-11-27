rm(list = ls())
library(Seurat)
library(ggplot2)
library(CellChat)
library(ggpubr)
library(ggsignif)
library(tidyverse)
library(cowplot)

memory.limit(1024000)

## 导入Seurat对象
GSE193533_obj_harmony = readRDS('/Volumes/GaoxyData/IAs_object/scRNA-seq/GSE193533_obj_harmony_recluster.rds')

## 提取亚组
# GSE193533_obj_harmony@meta.data$CellType2 = GSE193533_obj_harmony@meta.data$CellType
# GSE193533_obj_harmony@meta.data[rownames(GSE193533_obj_harmony@meta.data[which(GSE193533_obj_harmony@meta.data$group == 'High mRNAsi'), ]), 'CellType2'] = 'High mRNAsi'
# GSE193533_obj_harmony@meta.data[rownames(GSE193533_obj_harmony@meta.data[which(GSE193533_obj_harmony@meta.data$group == 'Low mRNAsi'), ]), 'CellType2'] = 'Low mRNAsi'

Idents(GSE193533_obj_harmony) = GSE193533_obj_harmony$CellChat_Group
#sub_object = subset(GSE193533_obj_harmony, idents = c('Macrophage', 'Fibroblast', 'Meyloid Cells', 'NK Cells', 'Endothelial Cells', 'T Cells',
#                                           'Monocytes', 'Epithelial Cells', 'Mast Cells', 'High mRNAsi'))
sub_object = subset(GSE193533_obj_harmony, idents = c('Fibroblast cluster1'))

## 创建CellChat对象
Idents(sub_object) = sub_object$CellType
cellchat = createCellChat(sub_object, group.by = 'CellType')
CellChatDB = CellChatDB.mouse
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)

# CellChatDB.use = subsetDB(CellChatDB, search = 'Secreted Signaling')
CellChatDB.use = CellChatDB
cellchat@DB = CellChatDB.use

cellchat = subsetData(cellchat)
cellchat = identifyOverExpressedGenes(cellchat)
cellchat = identifyOverExpressedInteractions(cellchat)
cellchat = computeCommunProb(cellchat)
cellchat = filterCommunication(cellchat, min.cells = 10)
cellchat = aggregateNet(cellchat)

saveRDS(cellchat, file = 'D:/LowmRNAsi_CellChat_obj_subset.RDS')

groupSize = as.numeric(table(cellchat@idents))

par(mfrow = c(1,2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = 'Number of interactions')
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F, title.name = 'Interaction weights/strength')


########################################################################################################################################
## 差异分析
cellchat.HighmRNAsi = readRDS('D:/HighmRNAsi_CellChat_obj_subset.RDS')
cellchat.LowmRNAsi = readRDS('D:/LowmRNAsi_CellChat_obj_subset.RDS')


group.new = levels(cellchat.HighmRNAsi@idents)
cellchat.LowmRNAsi = liftCellChat(cellchat.LowmRNAsi, group.new)


object.list = list(HighmRNAsi = cellchat.HighmRNAsi, LowmRNAsi = cellchat.LowmRNAsi)
cellchat.merge = mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)

weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
pdf('D:/Number_interaction.pdf', width = 8, height = 5)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

#不同细胞群之间的相互作用数量或强度的差异
par(mfrow = c(1,2), xpd=TRUE)
pdf('D:/Differential_interaction.pdf', width = 8, height = 5)
netVisual_diffInteraction(cellchat.merge, comparison = c(1,2), weight.scale = T, color.edge = c("#F22613","#16B4F2"))
netVisual_diffInteraction(cellchat.merge, comparison = c(1,2), weight.scale = T, measure = "weight", color.edge = c("#F22613","#16B4F2"))
dev.off()


## 差异信号通路
cellchat.HighmRNAsi = computeCommunProb(cellchat.HighmRNAsi)
cellchat.HighmRNAsi = computeCommunProbPathway(cellchat.HighmRNAsi)
cellchat.HighmRNAsi = netAnalysis_computeCentrality(cellchat.HighmRNAsi, slot.name = "netP") 


cellchat.LowmRNAsi = computeCommunProb(cellchat.LowmRNAsi)
cellchat.LowmRNAsi = computeCommunProbPathway(cellchat.LowmRNAsi)
cellchat.LowmRNAsi = netAnalysis_computeCentrality(cellchat.LowmRNAsi, slot.name = "netP") 

group.new = levels(cellchat.HighmRNAsi@idents)
cellchat.LowmRNAsi = liftCellChat(cellchat.LowmRNAsi, group.new)


object.list = list(HighmRNAsi = cellchat.HighmRNAsi, LowmRNAsi = cellchat.LowmRNAsi)
cellchat.merge2 = mergeCellChat(object.list, add.names = names(object.list), cell.prefix = T)

#计算和可视化通路距离
options(repr.plot.height = 12, repr.plot.width = 12)
rankSimilarity(cellchat.merge2, type = "functional")
#比较每个信号通路的整体信息流
options(repr.plot.height = 10, repr.plot.width = 18)
gg1 = rankNet(cellchat.merge2, mode = "comparison", stacked = T, do.stat = TRUE)+
        scale_fill_manual(values = c("#16B4F2","#F22613"))
gg2 = rankNet(cellchat.merge2, mode = "comparison", stacked = F, do.stat = TRUE)
pdf('D:/DE_Signal_pathway.pdf', height = 17, width = 10)
gg1
dev.off()





## 气泡图
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat.merge2, group.dataset = "datasets", pos.dataset = 'HighmRNAsi', features.name = 'HighmRNAsi', only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = 'HighmRNAsi')
# extract the ligand-receptor pairs with upregulated ligands in LS
net.up <- subsetCommunication(cellchat, net = net, datasets = "HighmRNAsi",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in NL, i.e.,downregulated in LS
net.down <- subsetCommunication(cellchat, net = net, datasets = "LowmRNAsi",ligand.logFC = -0.1, receptor.logFC = -0.1)

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
pathways.show <- c("EGF") 
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
pdf('D:/EGF_signal.pdf', width = 8, height = 5)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
dev.off()







