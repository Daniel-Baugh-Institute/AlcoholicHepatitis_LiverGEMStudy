#### Figure Replication ####

library(dataVisEasy)
library(limma)
library(ica)
library(matrixStats)
library(pheatmap)
library(pcaMethods)
library(tidyverse)
library(dplyr)
library(reshape2)
library(stats)

setwd('Set working directory to location of supplemental files')

## Figure 1B ##

geneExp <- as.matrix(read.table('DiseaseData_TPMnorm.txt'))
rownames(geneExp) <- geneExp[,1]
geneExp <- geneExp[,-1]
colnames(geneExp) <- geneExp[1,]
geneExp <- geneExp[-1,]
geneExp <- geneExp[,c(3,1,2,4,8,7,5,6)]

annots <- NULL
annots$disease_state <- factor(unique(colnames(geneExp)), levels=c("explant.AH", "severe.AH", "nonsevere.AH", "early.ASH",
                                                                   "healthy.control", "comp.cirrhosis", "HCV", "NASH"))
annots <- as.data.frame(annots)
rownames(annots) <- annots$disease_state
ann_colors <- list(disease_state = c(explant.AH = '#ea472a', severe.AH = '#ea8f59', nonsevere.AH = '#f2e32b', early.ASH = '#7aaf3e',
                                     healthy.control = '#f280ae', comp.cirrhosis = '#7580aa', HCV = '#5e8f9c', NASH = '#814284'))

initiate_params()
set_annotations(annots)
set_annot_cols(ann_colors)
set_annot_samps(c("disease_state"))
set_scale.range(c(-1,1))

class(geneExp) <- "numeric"
geneExp_disease <- geneExp
geneExp_scale <- (geneExp - rowMeans(geneExp)) / rowSds(geneExp)
geneExp_scale <- geneExp_scale[-c(which(is.na(rowSums(geneExp_scale)))),]
geneExp_scale_disease <- geneExp_scale


metgenes = read.table('livermetabolism-genes-Human1.txt')
metgenes_disease <- geneExp_disease[na.omit(match(metgenes$V1,rownames(geneExp_disease))),]
metgenes_scale_disease <- geneExp_scale_disease[na.omit(match(metgenes$V1,rownames(geneExp_scale_disease))),]


# show.rownames and show.colnames must be set to T if running ExtractMatrix function in line 54
p <- myHeatmapByAnnotation(metgenes_scale_disease, main = "Metabolic Gene Expression - Disease", groupings = "disease_state",
                           show.rownames = T, gaps.row = FALSE, gaps.col = TRUE, clust.rows = TRUE,
                           show.colnames = T, row.groups = 2)

res <- ExtractMatrix(metgenes_scale_disease, p) # extract matrix row-ordered as in heatmap
clusts <- extractClusters(res, heatmap = p, to.extract = "genes", nclusters = 2, GeneGroup_Name = "Gene_Clusters") #extracts clusters from row.groups in heatmap
clusts <- cbind(rownames(clusts),clusts)

## Figure 2A Right ##

structComp = read.table('structComp.txt', sep = ',')
modelID=c("healthy.control","severe.AH", "nonsevere.AH", "explant.AH", "early.ASH", "HCV", "NASH", "comp.cirrhosis")
colnames(structComp) = modelID
rownames(structComp) = modelID

initiate_params()
set_scale.range(c(0.95,1))

callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

pheatmap(structComp, color = colorRampPalette(c("blue", "black","yellow"))(100),
         clustering_callback = callback, legend_breaks = c(0.95,0.96,0.97,0.98,0.99,1),
         legend_labels = c("0.95","0.96","0.97","0.98","0.99","1"), breaks = seq(0.95,1,by= 5e-4),
         main = 'Model Hamming Similarity')

## Figure 2B ##

subsystem_mat = read.table('subsystemmatrix.txt', sep = ',')
subsystem_ID = as.character(read.table('subsystemIDs.txt')[[1]])

subsystemmat_annots <-subsystem_mat
rownames(subsystemmat_annots) <- subsystem_ID
colnames(subsystemmat_annots) <- c("healthy.control","severe.AH","nonsevere.AH","explant.AH",
                                   "early.ASH","HCV","NASH","comp.cirrhosis")


subCoverage = (subsystemmat_annots - rowMeans(subsystemmat_annots)) / rowMeans(subsystemmat_annots) * 100

# only show those having at least a 10% difference in one or more GEMs
incl = which(apply(abs(subCoverage) > 10, 1, any))
Names = subsystem_ID[incl]

plotcov = subCoverage[incl,]
annots <- as.data.frame(colnames(plotcov))
rownames(annots) <- annots[,1]
annots$disease_state <- factor(colnames(plotcov), levels=c("explant.AH", "severe.AH", "nonsevere.AH", "early.ASH",
                                                           "healthy.control", "comp.cirrhosis", "HCV", "NASH"))

ann_colors <- list(disease_state = c(explant.AH = '#ea472a', severe.AH = '#ea8f59', nonsevere.AH = '#f2e32b', early.ASH = '#7aaf3e',
                                     healthy.control = '#f280ae', comp.cirrhosis = '#7580aa', HCV = '#5e8f9c', NASH = '#814284'))


initiate_params()
set_annotations(annots)
set_annot_cols(ann_colors)
set_annot_samps(c("disease_state"))
set_scale.range(c(-2,2))


plotcov_order <- subsystemmat_annots[match(rownames(plotcov),rownames(subsystemmat_annots)),]
plotcov_order$sum <- rowSums(plotcov_order)
plotcov_order <- plotcov_order[order(plotcov_order$sum, decreasing =TRUE),-9]

# remove "Miscellaneous","Isolated","Pool reactions" subsystems
plotcov_order <- plotcov_order[-c(match(c("Miscellaneous","Isolated","Pool reactions"),rownames(plotcov_order))),]

# show.rownames and show.colnames must be set to T if running ExtractMatrix function in line 126
p <- myHeatmapByAnnotation(plotcov_order[1:10,], clust.rows = FALSE, clust.cols = FALSE, 
                           main = 'Subsystem Tissue Coverage with >10% Difference', groupings = "disease_state",
                           scale.rows = 'zscore', show.colnames = TRUE)

res <- ExtractMatrix(plotcov_order,p) # extract matrix row-ordered as in heatmap


## Figure 2C ##

sphing_genes <- read.table('sphingolipidgenes.txt')
sphing_genes_diseasescale <- metgenes_scale_disease[match(sphing_genes[,1],rownames(metgenes_scale_disease)),]
rownames(sphing_genes_diseasescale) <- sphing_genes[,2]
p <- myHeatmapByAnnotation(sphing_genes_diseasescale, main = "Sphingolipid Gene Expression", groupings = "disease_state",
                           show.rownames = T, gaps.row = FALSE, gaps.col = TRUE, clust.rows = F,
                           show.colnames = T, gaps.row.spec = c(11,13,16), fontsize.row = 14)


## Figure 3B ##

FBA_results = read.table('DiseaseState_Fluxes.txt', sep = '\t', header = 1)
FBA_mat = data.matrix(FBA_results[,4:11])
Human1_Rxn = FBA_results[,3]
Recon_Rxn = FBA_results[,2]
Subsystem = FBA_results[,1]
rownames(FBA_mat) = Human1_Rxn

colnames(FBA_mat) <- c("explant.AH", "severe.AH", "nonsevere.AH", "early.ASH",
                       "healthy.control", "comp.cirrhosis", "HCV", "NASH")

annots <- NULL
annots$disease_state <- factor(colnames(FBA_mat), levels=c("explant.AH", "severe.AH", "nonsevere.AH", "early.ASH",
                                                           "healthy.control", "comp.cirrhosis", "HCV", "NASH"))
annots <- as.data.frame(annots)
rownames(annots) <- annots$disease_state
ann_colors <- list(disease_state = c(explant.AH = '#ea472a', severe.AH = '#ea8f59', nonsevere.AH = '#f2e32b', early.ASH = '#7aaf3e',
                                     healthy.control = '#f280ae', comp.cirrhosis = '#7580aa', HCV = '#5e8f9c', NASH = '#814284'))

initiate_params()
set_annotations(annots)
set_annot_cols(ann_colors)
set_annot_samps(c("disease_state"))
set_scale.range(c(-1,1))

# show.rownames and show.colnames must be set to T if running ExtractMatrix function in line 170
p <- myHeatmapByAnnotation(FBA_mat, main = "All Fluxes", groupings = "disease_state", scale.rows = "zscore", 
                           row.groups = 2, show.rownames = T,show.colnames = T)

res <- ExtractMatrix(FBA_mat,p) # extract matrix row-ordered as in heatmap
clusts <- extractClusters(res, heatmap = p, to.extract = "genes", nclusters = 2, GeneGroup_Name = "Gene_Clusters") #extracts clusters from row.groups in heatmap


## Figure 3C ##

subsys_res <- list()
kw_pval <- c()
for(s in unique(Subsystem)){
  flux_test <- abs(FBA_mat[Subsystem==s,])
  flux_test_long <- melt(flux_test)
  if (nrow(flux_test_long)==8){
    flux_test_long$Var2 <- rownames(flux_test_long)
  }
  flux_test_long$Var2 <- factor(flux_test_long$Var2)
  flux_test_long$Var2 <- relevel(flux_test_long$Var2, ref="healthy.control")
  kw_test = kruskal.test(value ~ Var2, data = flux_test_long)
  if(kw_test$p.value < 0.05){
    res_greater<-pairwise.wilcox.test(flux_test_long$value, flux_test_long$Var2, 
                                      p.adjust.method = "none", paired=T, alternative = "greater", exact = F)
    res_less<-pairwise.wilcox.test(flux_test_long$value, flux_test_long$Var2, 
                                   p.adjust.method = "none", paired=T, alternative = "less", exact = F)
    res_greater$p.value[is.na(res_greater$p.value)] <- 0
    res_less$p.value[is.na(res_less$p.value)] <- 0
    bindlower = cbind(lower.tri(res_greater$p.value, diag = T)*res_greater$p.value, rep(0,7))
    bindlower = rbind(rep(0,8),bindlower)
    colnames(bindlower)[8] = "NASH"; rownames(bindlower)[1] = "healthy.control"
    bindupper = cbind(rep(0,7), upper.tri(res_less$p.value, diag = T)*t(res_less$p.value))
    bindupper = rbind(bindupper,rep(0,8))
    colnames(bindupper)[1] = "healthy.control"; rownames(bindupper)[8] = "NASH"
    subsys_res[[s]] <- bindlower + bindupper
    kw_pval = c(kw_pval,kw_test$p.value)
  }
  else{print(s)} # non-significant subsystems
}

pvals <-cbind(names(subsys_res),kw_pval)


# take the average flux for each subsystem
count = 0
sum = 0
avg_mat = matrix(NA,nrow = length(unique(Subsystem)), ncol = ncol(FBA_mat))
for (i in 1:length(unique(Subsystem))){
  for (z in 1:ncol(FBA_mat)){
    for (j in 1:nrow(FBA_mat)){
      if (unique(Subsystem)[i] == Subsystem[j]){
        sum = sum + abs(FBA_mat[j,z])
        count = count + 1
      }
    }
    avg_mat[i,z] = sum / count
    sum = 0
    count = 0
  }
}

rownames(avg_mat) = unique(Subsystem)
colnames(avg_mat) = colnames(FBA_mat)

avg_mat = avg_mat[match(names(subsys_res),rownames(avg_mat)),]

# remove "Miscellaneous","Pool reactions" subsystems
avg_mat <- avg_mat[-c(match(c("Miscellaneous","Pool reactions"),rownames(avg_mat))),]

p <- myHeatmapByAnnotation(avg_mat, main = "Subsystems with Differential Fluxes", groupings = "disease_state", scale.rows = "zscore", 
                           show.rownames = T,show.colnames = F, gaps.col = T)

# number of reactions in subsystems from the above heatmap
numrxn <- as.data.frame(table(Subsystem[which(Subsystem%in%rownames(avg_mat))]))
numrxn


## Figure 4A ##

geneExp_scale_disease_pca <- pcaMethods::pca(t(geneExp_scale_disease))
geneExp_scale_disease_pcascrs <- pcaMethods::scores(geneExp_scale_disease_pca)
geneExp_scale_disease_pcaplot <- data.frame(geneExp_scale_disease_pcascrs,
                                            Samples = colnames(geneExp_scale_disease))
geneExp_scale_disease_pcavars <- round(geneExp_scale_disease_pca@R2 * 100, digits = 2)

pcolor<- c("#ea472a","#ea8f59","#f2e32b", "#7aaf3e","#f280ae","#7580aa","#5e8f9c", "#814284")
names(pcolor) <- c("explant.AH","severe.AH","nonsevere.AH", "early.ASH","healthy.control","comp.cirrhosis","HCV","NASH")
geneExp_scale_disease_pcaplot$Samples <- factor(geneExp_scale_disease_pcaplot$Samples, levels = names(pcolor))
pcolor <- colorRampPalette(pcolor)

geneExp_scale_disease_pcaplot$shapes = c(0,0,0,0,1,2,2,2)
geneExp_scale_disease_pcaplot$shapes <- factor(geneExp_scale_disease_pcaplot$shapes)

# put healthy point on top of plot
geneExp_scale_disease_pcaplot <- geneExp_scale_disease_pcaplot %>% slice(1:4,6:8, 5)

ggplot(geneExp_scale_disease_pcaplot, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = shapes, color=Samples), size=16)+
  scale_shape_manual(values=c(16, 17, 15)) + 
  scale_color_manual(values = pcolor(8)) +
  xlab(paste0("PC1 (",geneExp_scale_disease_pcavars[1],"%)")) +
  ylab(paste0("PC2 (",geneExp_scale_disease_pcavars[2],"%)")) +
  ggtitle("PCA – All Genes (23873 Genes)") +
  theme_bw() + theme(panel.grid = element_blank(), 
                     plot.title = element_text(hjust = 0.5, size = 35), 
                     axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
                     legend.position = "right", axis.text.x=element_text(colour="black"),
                     axis.text.y=element_text(colour="black"),
                     panel.border = element_rect(color = "black",
                                                 fill = NA,
                                                 size = 2))

## Figure 4B ##

nometgenes_scale_disease <- geneExp_scale_disease[-c(which(rownames(metgenes_scale_disease)%in% rownames(geneExp_scale_disease))),]
nometgenes_scale_disease_pca <- pcaMethods::pca(t(nometgenes_scale_disease))
nometgenes_scale_disease_pcascrs <- pcaMethods::scores(nometgenes_scale_disease_pca)
nometgenes_scale_disease_pcaplot <- data.frame(nometgenes_scale_disease_pcascrs,
                                               Samples = colnames(nometgenes_scale_disease))
nometgenes_scale_disease_pcavars <- round(nometgenes_scale_disease_pca@R2 * 100, digits = 2)


pcolor<- c("#ea472a","#ea8f59","#f2e32b", "#7aaf3e","#f280ae","#7580aa","#5e8f9c", "#814284")
names(pcolor) <- c("explant.AH","severe.AH","nonsevere.AH", "early.ASH","healthy.control","comp.cirrhosis","HCV","NASH")
nometgenes_scale_disease_pcaplot$Samples <- factor(nometgenes_scale_disease_pcaplot$Samples, levels = names(pcolor))
pcolor <- colorRampPalette(pcolor)

nometgenes_scale_disease_pcaplot$shapes = c(0,0,0,0,1,2,2,2)
nometgenes_scale_disease_pcaplot$shapes <- factor(nometgenes_scale_disease_pcaplot$shapes)

# put healthy point on top of plot
nometgenes_scale_disease_pcaplot <- nometgenes_scale_disease_pcaplot %>% slice(1:4,6:8, 5)

ggplot(nometgenes_scale_disease_pcaplot, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = shapes, color=Samples), size=16)+
  scale_shape_manual(values=c(16, 17, 15)) + 
  scale_color_manual(values = pcolor(8)) +
  xlab(paste0("PC1 (",nometgenes_scale_disease_pcavars[1],"%)")) +
  ylab(paste0("PC2 (",nometgenes_scale_disease_pcavars[2],"%)")) +
  ggtitle("PCA – No Metabolic Genes (20315 Genes)") +
  theme_bw() + theme(panel.grid = element_blank(), 
                     plot.title = element_text(hjust = 0.5, size = 35), 
                     axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
                     legend.position = "right", axis.text.x=element_text(colour="black"),
                     axis.text.y=element_text(colour="black"),
                     panel.border = element_rect(color = "black",
                                                 fill = NA,
                                                 size = 2))


## Figure 4C ##

metgenes_scale_disease_pca <- pcaMethods::pca(t(metgenes_scale_disease))
metgenes_scale_disease_pcascrs <- pcaMethods::scores(metgenes_scale_disease_pca)
metgenes_scale_disease_pcaplot <- data.frame(metgenes_scale_disease_pcascrs,
                                             Samples = colnames(metgenes_scale_disease))
metgenes_scale_disease_pcavars <- round(metgenes_scale_disease_pca@R2 * 100, digits = 2)


pcolor<- c("#ea472a","#ea8f59","#f2e32b", "#7aaf3e","#f280ae","#7580aa","#5e8f9c", "#814284")
names(pcolor) <- c("explant.AH","severe.AH","nonsevere.AH", "early.ASH","healthy.control","comp.cirrhosis","HCV","NASH")
metgenes_scale_disease_pcaplot$Samples <- factor(metgenes_scale_disease_pcaplot$Samples, levels = names(pcolor))
pcolor <- colorRampPalette(pcolor)

metgenes_scale_disease_pcaplot$shapes = c(0,0,0,0,1,2,2,2)
metgenes_scale_disease_pcaplot$shapes <- factor(metgenes_scale_disease_pcaplot$shapes)

# put healthy point on top of plot
metgenes_scale_disease_pcaplot <- metgenes_scale_disease_pcaplot %>% slice(1:4,6:8, 5)

ggplot(metgenes_scale_disease_pcaplot, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = shapes, color=Samples), size=16)+
  scale_shape_manual(values=c(16, 17, 15)) + 
  scale_color_manual(values = pcolor(8)) +
  xlab(paste0("PC1 (",metgenes_scale_disease_pcavars[1],"%)")) +
  ylab(paste0("PC2 (",metgenes_scale_disease_pcavars[2],"%)")) +
  ggtitle("PCA – Metabolic Genes Only (3558 Genes)") +
  theme_bw() + theme(panel.grid = element_blank(), 
                     plot.title = element_text(hjust = 0.5, size = 35), 
                     axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
                     legend.position = "right", axis.text.x=element_text(colour="black"),
                     axis.text.y=element_text(colour="black"),
                     panel.border = element_rect(color = "black",
                                                 fill = NA,
                                                 size = 2))


## Figure 5D ##

FBA_mat_norm <- (FBA_mat - rowSums(FBA_mat)) / rowSds(FBA_mat)
FBA_mat_norm_pca <- pcaMethods::pca(t(FBA_mat_norm))
FBA_mat_norm_pcascrs <- pcaMethods::scores(FBA_mat_norm_pca)
FBA_mat_norm_pcaplot <- data.frame(FBA_mat_norm_pcascrs,
                                   Samples = colnames(FBA_mat_norm))
FBA_mat_norm_pcavars <- round(FBA_mat_norm_pca@R2 * 100, digits = 2)

pcolor<- c("#ea472a","#ea8f59","#f2e32b", "#7aaf3e","#f280ae","#7580aa","#5e8f9c", "#814284")
names(pcolor) <- c("explant.AH","severe.AH","nonsevere.AH", "early.ASH","healthy.control","comp.cirrhosis","HCV","NASH")
FBA_mat_norm_pcaplot$Samples <- factor(FBA_mat_norm_pcaplot$Samples, levels = names(pcolor))
pcolor <- colorRampPalette(pcolor)

FBA_mat_norm_pcaplot$shapes = c(0,0,0,0,1,2,2,2)
FBA_mat_norm_pcaplot$shapes <- factor(FBA_mat_norm_pcaplot$shapes)

# put healthy point on top of plot
FBA_mat_norm_pcaplot <- FBA_mat_norm_pcaplot %>% slice(1:4,6:8, 5)

ggplot(FBA_mat_norm_pcaplot, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = shapes, color=Samples), size=16)+
  scale_shape_manual(values=c(16, 17, 15)) + 
  scale_color_manual(values = pcolor(8)) +
  xlab(paste0("PC1 (",FBA_mat_norm_pcavars[1],"%)")) +
  ylab(paste0("PC2 (",FBA_mat_norm_pcavars[2],"%)")) +
  ggtitle("PCA – Metabolic Fluxes (886)") +
  theme_bw() + theme(panel.grid = element_blank(), 
                     plot.title = element_text(hjust = 0.5, size = 35), 
                     axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
                     legend.position = "right", axis.text.x=element_text(colour="black"),
                     axis.text.y=element_text(colour="black"),
                     panel.border = element_rect(color = "black",
                                                 fill = NA,
                                                 size = 2))

## Figure 5E ##
nsca <- function(A)
{
  Dr <- apply(A, 1, sum)
  Dc <- apply(A, 2, sum)
  
  eig.res <- eigen(diag(1 / sqrt(Dr)) %*% A %*% diag(1 / sqrt(Dc)))
  r <- diag(1 / Dr) %*% (eig.res$vectors)[, 1:5]
  rownames(r) <- rownames(A)
  r
}

ica.dat <- icafast(geneExp_disease,5)
ica.matrix <- ica.dat$M
rownames(ica.matrix) <- colnames(geneExp_disease)
ica.matrix <- as.matrix(dist(ica.matrix))
dat.mst <- ape::mst(ica.matrix)
plot(dat.mst, graph = 'nsca')

pdata<-as.data.frame(nsca(dat.mst))
pdata$group <- colnames(geneExp_disease)
#     custom color
pcolor<- c("#ea472a","#ea8f59","#f2e32b", "#7aaf3e","#f280ae","#7580aa","#5e8f9c", "#814284")
names(pcolor) <- c("explant.AH","severe.AH","nonsevere.AH", "early.ASH","healthy.control","comp.cirrhosis","HCV","NASH")
pdata$group <- factor(pdata$group, levels = names(pcolor))
pcolor <- colorRampPalette(pcolor)

pdata$shapes = c(0,0,0,0,1,2,2,2)
pdata$shapes <- factor(pdata$shapes)
pdata$ptsize = c(1,1,1,1,2,1,1,1)

# put healthy point on the bottom of the plot
pdata <- pdata %>% slice(5,1:4,6:8)

ggplot(as.data.frame(pdata), aes(V2, V3))+ 
  geom_point(aes(shape = shapes, color=group, size=ptsize)) +
  scale_size_continuous(range = c(16, 24)) +
  scale_shape_manual(values=c(16, 17, 15)) + 
  scale_color_manual(values = pcolor(8)) +
  xlim(-0.4,0.6) + ylim(-0.6,0.4) +
  xlab("nsca axis 1") + ylab("nsca axis 2") +
  ggtitle("Minimum Spanning Tree - All Gene Space") + 
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank())


## Figure 5F ##

nometgenes_disease <- geneExp_disease[-c(match(rownames(metgenes_disease),rownames(geneExp_disease))),]
ica.dat <- icafast(nometgenes_disease,5)
ica.matrix <- ica.dat$M
rownames(ica.matrix) <- colnames(nometgenes_disease)
ica.matrix <- as.matrix(dist(ica.matrix))
dat.mst <- ape::mst(ica.matrix)
plot(dat.mst, graph = 'nsca')

pdata<-as.data.frame(nsca(dat.mst))
pdata$group <- colnames(nometgenes_disease)
#     custom color
pcolor<- c("#ea472a","#ea8f59","#f2e32b", "#7aaf3e","#f280ae","#7580aa","#5e8f9c", "#814284")
names(pcolor) <- c("explant.AH","severe.AH","nonsevere.AH", "early.ASH","healthy.control","comp.cirrhosis","HCV","NASH")
pdata$group <- factor(pdata$group, levels = names(pcolor))
pcolor <- colorRampPalette(pcolor)

pdata$shapes = c(0,0,0,0,1,2,2,2)
pdata$shapes <- factor(pdata$shapes)
pdata$ptsize = c(1,1,1,1,2,1,1,1)

# put healthy point on the bottom of the plot
pdata <- pdata %>% slice(5,1:4,6:8)

ggplot(as.data.frame(pdata), aes(V2, V3))+ 
  geom_point(aes(shape = shapes,color=group, size=ptsize))+
  scale_size_continuous(range = c(16, 24)) +
  scale_shape_manual(values=c(16, 17, 15)) + 
  scale_color_manual(values = pcolor(8)) +
  xlab("nsca axis 1") + ylab("nsca axis 2") +
  ggtitle("Minimum Spanning Tree - No Metabolic Genes") + 
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank())


## Figure 5G ##

ica.dat <- icafast(metgenes_disease,5)
ica.matrix <- ica.dat$M
rownames(ica.matrix) <- colnames(metgenes_disease)
ica.matrix <- as.matrix(dist(ica.matrix))
dat.mst <- ape::mst(ica.matrix)
plot(dat.mst, graph = 'nsca')

pdata<-as.data.frame(nsca(dat.mst))
pdata$group <- colnames(metgenes_disease)
#     custom color
pcolor<- c("#ea472a","#ea8f59","#f2e32b", "#7aaf3e","#f280ae","#7580aa","#5e8f9c", "#814284")
names(pcolor) <- c("explant.AH","severe.AH","nonsevere.AH", "early.ASH","healthy.control","comp.cirrhosis","HCV","NASH")
pdata$group <- factor(pdata$group, levels = names(pcolor))
pcolor <- colorRampPalette(pcolor)

pdata$shapes = c(0,0,0,0,1,2,2,2)
pdata$shapes <- factor(pdata$shapes)
pdata$ptsize = c(1,1,1,1,2,1,1,1)

# put healthy point on the bottom of the plot
pdata <- pdata %>% slice(5,1:4,6:8)

ggplot(as.data.frame(pdata), aes(V2, V3))+ 
  geom_point(aes(shape = shapes,color=group, size=ptsize))+
  scale_size_continuous(range = c(16, 24)) +
  scale_shape_manual(values=c(16, 17, 15)) + 
  scale_color_manual(values = pcolor(8)) +
  xlab("nsca axis 1") + ylab("nsca axis 2") +
  ggtitle("Minimum Spanning Tree - Metabolic Gene Space") + 
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank())


## Figure 5H ##

flux <- FBA_mat
ica.dat <- icafast(flux,8)
ica.matrix <- ica.dat$M
rownames(ica.matrix) <- colnames(flux)
ica.matrix <- as.matrix(dist(ica.matrix))
dat.mst <- ape::mst(ica.matrix)

pdata<-as.data.frame(nsca(dat.mst))
pdata$group <- colnames(flux)
#     custom color
pcolor <- c("#ea472a","#ea8f59","#f2e32b", "#7aaf3e","#f280ae","#7580aa","#5e8f9c","#814284")
names(pcolor) <- c("explant.AH","severe.AH","nonsevere.AH", "early.ASH", "healthy.control",
                   "comp.cirrhosis","HCV","NASH")
pdata$group <- factor(pdata$group, levels = names(pcolor))
pcolor <- colorRampPalette(pcolor)

pdata$shapes = c(0,0,0,0,1,2,2,2)
pdata$shapes <- factor(pdata$shapes)
pdata$ptsize = c(1,1,1,1,2,1,1,1)

# put healthy point on the bottom of the plot
pdata <- pdata %>% slice(5,1:4,6:8)

ggplot(as.data.frame(pdata), aes(V2, V3))+ 
  geom_point(aes(shape = shapes,color=group, size=ptsize))+
  scale_size_continuous(range = c(16, 24)) +
  scale_shape_manual(values=c(16, 17, 15)) + 
  scale_color_manual(values = pcolor(8)) +
  xlab("nsca axis 1") + ylab("nsca axis 2") +
  ylim(-0.6,0.4) + 
  ggtitle("Minimum Spanning Tree - Metabolic Fluxes") + 
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank())

## Figure 5A ##
gr <- read.table('transportReaction_genes.txt')
SLC_scale_disease <- geneExp_scale_disease[na.omit(match(gr$V1,rownames(geneExp_scale_disease))),]
SLC_scale_disease_pca <- pcaMethods::pca(t(SLC_scale_disease))
SLC_scale_disease_pcascrs <- pcaMethods::scores(SLC_scale_disease_pca)



SLC_scale_disease_pcaplot <- data.frame(SLC_scale_disease_pcascrs,
                                        Samples = colnames(SLC_scale_disease))
SLC_scale_disease_pcavars <- round(SLC_scale_disease_pca@R2 * 100, digits = 2)

pcolor<- c("#ea472a","#ea8f59","#f2e32b", "#7aaf3e","#f280ae","#7580aa","#5e8f9c", "#814284")
names(pcolor) <- c("explant.AH","severe.AH","nonsevere.AH", "early.ASH","healthy.control","comp.cirrhosis","HCV","NASH")
SLC_scale_disease_pcaplot$Samples <- factor(SLC_scale_disease_pcaplot$Samples, levels = names(pcolor))
pcolor <- colorRampPalette(pcolor)

SLC_scale_disease_pcaplot$shapes = c(0,0,0,0,1,2,2,2)
SLC_scale_disease_pcaplot$shapes <- factor(SLC_scale_disease_pcaplot$shapes)

# put healthy point on top of plot
SLC_scale_disease_pcaplot <- SLC_scale_disease_pcaplot %>% slice(1:4,6:8, 5)

ggplot(SLC_scale_disease_pcaplot, aes(x = PC1, y = PC2)) + 
  geom_point(aes(shape = shapes, color=Samples), size=16)+
  scale_shape_manual(values=c(16, 17, 15)) + 
  scale_color_manual(values = pcolor(8)) +
  xlab(paste0("PC1 (",geneExp_scale_disease_pcavars[1],"%)")) +
  ylab(paste0("PC2 (",geneExp_scale_disease_pcavars[2],"%)")) +
  ggtitle("PCA – (SLC genes)") +
  theme_bw() + theme(panel.grid = element_blank(), 
                     plot.title = element_text(hjust = 0.5, size = 35), 
                     axis.text = element_text(size = 20), axis.title = element_text(size = 25), 
                     legend.position = "right", axis.text.x=element_text(colour="black"),
                     axis.text.y=element_text(colour="black"),
                     panel.border = element_rect(color = "black",
                                                 fill = NA,
                                                 size = 2))

## Figure 5B ##

ica.dat <- icafast(SLC_scale_disease,5)
ica.matrix <- ica.dat$M
rownames(ica.matrix) <- colnames(geneExp_disease)
ica.matrix <- as.matrix(dist(ica.matrix))
dat.mst <- ape::mst(ica.matrix)
plot(dat.mst, graph = 'nsca')

pdata<-as.data.frame(nsca(dat.mst))
pdata$group <- colnames(geneExp_disease)
#     custom color
pcolor<- c("#ea472a","#ea8f59","#f2e32b", "#7aaf3e","#f280ae","#7580aa","#5e8f9c", "#814284")
names(pcolor) <- c("explant.AH","severe.AH","nonsevere.AH", "early.ASH","healthy.control","comp.cirrhosis","HCV","NASH")
pdata$group <- factor(pdata$group, levels = names(pcolor))
pcolor <- colorRampPalette(pcolor)

pdata$shapes = c(0,0,0,0,1,2,2,2)
pdata$shapes <- factor(pdata$shapes)
pdata$ptsize = c(1,1,1,1,2,1,1,1)

# put healthy point on the bottom of the plot
pdata <- pdata %>% slice(5,1:4,6:8)

ggplot(as.data.frame(pdata), aes(V2, V3))+ 
  geom_point(aes(shape = shapes, color=group, size=ptsize)) +
  scale_size_continuous(range = c(16, 24)) +
  scale_shape_manual(values=c(16, 17, 15)) + 
  scale_color_manual(values = pcolor(8)) +
  xlab("nsca axis 1") + ylab("nsca axis 2") +
  ggtitle("Minimum Spanning Tree - SLC Genes") + 
  theme(plot.title = element_text(hjust = 0.5), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.background = element_blank())


## Figure 5C ##
transport <- read.table('transportfluxes.txt', sep = '\t')
rownames(transportfluxes) <- transportfluxes[,3]
transportfluxes <- transport[,4:11]
colnames(transportfluxes) <- colnames(FBA_mat)

# show.rownames and show.colnames must be set to T if running ExtractMatrix function in line 257
p <- myHeatmapByAnnotation(abs(transportfluxes), main = "Transport Reaction Fluxes", groupings = "disease_state", scale.rows = "zscore", 
                           show.rownames = T,show.colnames = T, gaps.col = T, row.groups = 4)

res <- ExtractMatrix(transportfluxes,p) # extract matrix row-ordered as in heatmap
clusts <- extractClusters(res, heatmap = p, to.extract = "genes", nclusters = 4, GeneGroup_Name = "Gene_Clusters") #extracts clusters from row.groups in heatmap
clusts$Gene_Rules <- transport[match(rownames(clusts),transport[,3]),12]

## frequency tables of most common gene rule for each cluster

#cluster A: 230 total fluxes; 89 fluxes with most frequent gene rule ENSG00000103257
View(as.data.frame(table(clusts[which(clusts$Gene_Clusters=='A'),2]))) 

#cluster B: 120 total fluxes; 30 fluxes with most frequent gene rule ENSG00000021488 and ENSG00000103064 and ENSG00000138079
View(as.data.frame(table(clusts[which(clusts$Gene_Clusters=='B'),2]))) 

#cluster C: 227 total fluxes; 70 fluxes with most frequent gene rule ENSG00000134538
View(as.data.frame(table(clusts[which(clusts$Gene_Clusters=='C'),2]))) 

#cluster C: 64 total fluxes; 19 fluxes with most frequent gene rule ENSG00000084453 and ENSG00000111700
View(as.data.frame(table(clusts[which(clusts$Gene_Clusters=='D'),2]))) 

# reorder res matrix so reactions with most frequent gene rules appear first in the clusters
res <- as.data.frame(res)
res <- res %>% slice(which(clusts$Gene_Clusters=='A' & clusts$Gene_Rules=='ENSG00000103257'),
                     which(clusts$Gene_Clusters=='A' & !(clusts$Gene_Rules=='ENSG00000103257')),
                     which(clusts$Gene_Clusters=='B' & clusts$Gene_Rules=='ENSG00000021488 and ENSG00000103064 and ENSG00000138079'),
                     which(clusts$Gene_Clusters=='B' & !(clusts$Gene_Rules=='ENSG00000021488 and ENSG00000103064 and ENSG00000138079')),
                     which(clusts$Gene_Clusters=='C' & clusts$Gene_Rules=='ENSG00000134538'),
                     which(clusts$Gene_Clusters=='C' & !(clusts$Gene_Rules=='ENSG00000134538')),
                     which(clusts$Gene_Clusters=='D' & clusts$Gene_Rules=='ENSG00000084453 and ENSG00000111700'),
                     which(clusts$Gene_Clusters=='D' & !(clusts$Gene_Rules=='ENSG00000084453 and ENSG00000111700')),)

gene.groups <- c(rep(1,length(which(clusts$Gene_Clusters=='A' & clusts$Gene_Rules=='ENSG00000103257'))),
                 rep(2,length(which(clusts$Gene_Clusters=='A' & !(clusts$Gene_Rules=='ENSG00000103257')))),
                 rep(3,length(which(clusts$Gene_Clusters=='B' & clusts$Gene_Rules=='ENSG00000021488 and ENSG00000103064 and ENSG00000138079'))),
                 rep(4,length(which(clusts$Gene_Clusters=='B' & !(clusts$Gene_Rules=='ENSG00000021488 and ENSG00000103064 and ENSG00000138079')))),
                 rep(5,length(which(clusts$Gene_Clusters=='C' & clusts$Gene_Rules=='ENSG00000134538'))),
                 rep(6,length(which(clusts$Gene_Clusters=='C' & !(clusts$Gene_Rules=='ENSG00000134538')))),
                 rep(7,length(which(clusts$Gene_Clusters=='D' & clusts$Gene_Rules=='ENSG00000084453 and ENSG00000111700'))),
                 rep(8,length(which(clusts$Gene_Clusters=='D' & !(clusts$Gene_Rules=='ENSG00000084453 and ENSG00000111700')))))
gene.groups <- as.data.frame(gene.groups)
rownames(gene.groups) <- rownames(res)
colnames(gene.groups) <- 'group1'
gene.groups$group2 <- c(rep(1,length(which(clusts$Gene_Clusters=='A'))), rep(2,length(which(clusts$Gene_Clusters=='B'))),
                        rep(3,length(which(clusts$Gene_Clusters=='C'))), rep(4,length(which(clusts$Gene_Clusters=='D'))))

gene.groups$group1 <- factor(gene.groups$group1)
gene.groups$group2 <- factor(gene.groups$group2)



set_annotations.genes(gene.groups)
set_annot_genes(c("group1", "group2"))

p <- myHeatmapByAnnotation(abs(res), main = "Transport Reaction Fluxes", groupings = "disease_state", scale.rows = "zscore", 
                           show.rownames = T,show.colnames = T, gaps.col = T, row.groups = 4,
                           groupings.genes = c("group1","group2") , groupings.genes.gaps = c(1,2))


