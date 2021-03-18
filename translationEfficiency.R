library(reshape2)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(dendsort)
library(dendextend)
library(DESeq2)
library(pheatmap)
library(magrittr)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#Load count values for organs
organ_mRNA <- read.delim("F:/Ribosome_profiling/Mammals/Mouse/Injections/Translation_efficiency/organ_mRNA.txt", sep = "\t", header=TRUE, check.names = FALSE, stringsAsFactors = FALSE)
organ_Ribo <- read.delim("F:/Ribosome_profiling/Mammals/Mouse/Injections/Translation_efficiency/organ_Ribo.txt", sep = "\t", header=TRUE, check.names = FALSE, stringsAsFactors = FALSE)

#remove gene records with low read count (<10 reads in any of the mRNA-seq samples and <10 reads in any of the Ribo-seq samples)
organ_mRNA_sub <- apply(organ_mRNA[,-1], 1, function(x) {all(x >= 10)})
organ_Ribo_sub <- apply(organ_Ribo[,-1], 1, function(x) {all(x >= 10 )})
organ_mRNA  <- organ_mRNA[organ_mRNA_sub,]
organ_Ribo  <- organ_Ribo[organ_Ribo_sub,]

#Calculate TE
common_genes  <- match(organ_Ribo$geneID, organ_mRNA$geneID)
common_organs <- match(colnames(organ_Ribo), colnames(organ_mRNA))
common_mRNA <- organ_mRNA[na.omit(common_genes),as.integer(na.omit(common_organs))]
common_Ribo <- organ_Ribo[which(!is.na(common_genes)),which(!is.na(common_organs))]
organ_TE <- as.data.frame(cbind("geneID"=common_mRNA$geneID, log2(common_Ribo[,-1] / common_mRNA[,-1])))

#log transform mRNA-seq and Ribo-seq counts (TE are already log transformed)
organ_mRNA <- data.frame(cbind("geneID"=organ_mRNA$geneID, log10(organ_mRNA[,-1])  ))
organ_Ribo <- data.frame(cbind("geneID"=organ_Ribo$geneID, log10(organ_Ribo[,-1])  ))
organ_mRNA$geneID <- as.character(organ_mRNA$geneID)
organ_Ribo$geneID <- as.character(organ_Ribo$geneID)
organ_TE$geneID <- as.character(organ_TE$geneID)

# Use correlation between variables as a distance for TE, mRNA-seq and Ribo-seq
cormat_TE   <- cor(x=organ_TE[,2:ncol(organ_TE)], y=organ_TE[,2:ncol(organ_TE)], method="pearson")
cormat_mRNA <- cor(x=organ_mRNA[,2:ncol(organ_mRNA)], y=organ_mRNA[,2:ncol(organ_mRNA)], method="pearson")
cormat_Ribo <- cor(x=organ_Ribo[,2:ncol(organ_Ribo)], y=organ_Ribo[,2:ncol(organ_Ribo)], method="pearson")
dist_TE <- as.dist((1-cormat_TE)/2)
dist_mRNA <- as.dist((1-cormat_mRNA)/2)
dist_Ribo <- as.dist((1-cormat_Ribo)/2)
hc_TE <- hclust(dist_TE)
hc_mRNA <- hclust(dist_mRNA)
hc_Ribo <- hclust(dist_Ribo)
cormat_TE <-cormat_TE[rev(hc_TE$order),hc_TE$order]
cormat_mRNA <-cormat_mRNA[rev(hc_mRNA$order),hc_mRNA$order]
cormat_Ribo <-cormat_Ribo[rev(hc_Ribo$order),hc_Ribo$order]
melted_cormat_TE <- melt(cormat_TE)
melted_cormat_mRNA <- melt(cormat_mRNA)
melted_cormat_Ribo <- melt(cormat_Ribo)

#Organ TE correlation heatmap
cairo_pdf(filename = "Organ_TE_heatmap.pdf", height=5, width=5)
ggplot(data = melted_cormat_TE, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "black", size=0.3) +
  scale_fill_gradient(low = "#ffffcc", high = "#b10026") +
  theme_classic() +
  coord_equal() +
  theme(axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y =element_blank())
dev.off()

#Organ mRNA-seq correlation heatmap
cairo_pdf(filename = "Organ_mRNA_heatmap.pdf", height=5, width=5)
ggplot(data = melted_cormat_mRNA, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "black", size=0.3) +
  scale_fill_gradient(low = "#f7fbff", high = "#084594") +
  theme_classic() +
  coord_equal() +
  theme(axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y =element_blank())
dev.off()

#Organ Ribo-seq correlation heatmap
cairo_pdf(filename = "Organ_Ribo_heatmap.pdf", height=5, width=5)
ggplot(data = melted_cormat_Ribo, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile(color = "black", size=0.3) +
  scale_fill_gradient(low = "#ffffcc", high = "#006837") +
  theme_classic() +
  coord_equal() +
  theme(axis.text.x = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), axis.title.y =element_blank())
dev.off()

#Build TE heatmap for individual genes
gene_TE <- scale(organ_TE[,-1]) %>% as.matrix() %>% set_rownames(organ_TE$geneID)
dist_genes   <- as.dist(1-cor(t(gene_TE), method = 'pearson'))
dist_samples <- t(gene_TE) %>% dist(gene_TE, method = 'euclidean')

# prepare list of ribosomal genes
rps <- grepl('Rp[ls]\\d+', row.names(gene_TE))
row.names(gene_TE)[rps]

annotation_row <- data.frame(Ribosome = c(rep('TRUE', 500), rep('FALSE', nrow(gene_TE) - 500)) ) %>% set_rownames(row.names(gene_TE))

p <- pheatmap(gene_TE, 
              cluster_rows = hclust(dist_genes, method = 'average'),
              cluster_cols = hclust(dist_samples, method = 'average'),
              scale = 'row',
              #color = colorRampPalette(c("blue","black","yellow"))(256), #different color scheme
              color = rev(colorRampPalette(brewer.pal(11, 'RdBu'))(256)),
              #annotation_row = annotation_row,
              border_color = NA, show_rownames = FALSE, fontsize = 14, fontsize_col = 16)


png(filename = "gene_TE_pheatmap4.png", height=12, width=10, units="in", bg="white", res=300)
p
dev.off()

#pancreas has unusually low TE of ribosomal genes comparet do other organs
rps_df <- gene_TE[rps,]
png(filename = "ribosomal_TE_boxplot.png", height=5, width=10, units="in", bg="white", res=300)
boxplot(rps_df, las = 2)
dev.off()

#amino acyl tRNA synthetases
organ_mRNA_scaled <- scale(organ_mRNA[,-1]) %>% as.matrix() %>% set_rownames(organ_mRNA$geneID)
organ_Ribo_scaled <- scale(organ_Ribo[,-1]) %>% as.matrix() %>% set_rownames(organ_Ribo$geneID)

aa_synth <- c('Aars', 'Aimp1', 'Aimp2', 'Cars', 'Dars', 'Eprs', 'Farsa', 'Farsb', 'Gars', 'Hars', 'Iars', 'Kars', 'Lars', 'Mars', 'Nars', 'Ppa1', 'Ppa2', 'Qars', 'Rars', 'Sars', 'Tars', 'Wars', 'Vars', 'Yars' ) 
elongation <- c('Eef1a1', 'Eef1d', 'Eef1e1', 'Eef1g', 'Eef2')

aa_synth_df   <- organ_Ribo_scaled[row.names(organ_Ribo_scaled) %in% aa_synth,]
elongation_df <- organ_Ribo_scaled[row.names(organ_Ribo_scaled) %in% elongation,]

boxplot(aa_synth_df)
boxplot(elongation_df)

#extract gene order information from the heatmap
gene_TE_clustered <- gene_TE[p$tree_row$order, p$tree_col$order]

#interactive heatmap for cluster boundaries estimation
gene_TE_numbered <- gene_TE %>% set_rownames(c(1:nrow(.)))
dist_genes   <- as.dist(1-cor(t(gene_TE_numbered), method = 'pearson'))
dist_samples <- t(gene_TE_numbered) %>% dist(gene_TE_numbered, method = 'euclidean')
library(d3heatmap)
d3heatmap(gene_TE_numbered,
      Rowv = rev(as.dendrogram(hclust(dist_genes, method = 'average'))),
      Colv = as.dendrogram(hclust(dist_samples, method = 'average')),
      colors = rev(colorRampPalette(brewer.pal(11, 'RdBu'))(256)),
      show_grid = FALSE,
      scale = 'row'
)

#GSEA enrichment plot
GSEA_functions <- read.delim("F:/Ribosome_profiling/Mammals/Mouse/Injections/Translation_efficiency/GSEA_functions.txt", sep = "\t", header=FALSE, check.names = FALSE, stringsAsFactors = FALSE)
names(GSEA_functions) <- c("category", "score")
GSEA_functions$direction <- sapply(GSEA_functions$score, function(x) {if(x>0) {return("up")} else {return("down")}     })
GSEA_functions <- GSEA_functions[rev(order(abs(GSEA_functions$score))),]

png(filename = "GSEA.png", height=5, width=6, units="in", bg="white", res=600)
ggplot(data = GSEA_functions) +
  geom_bar(aes(x = reorder(category, rev(1:nrow(GSEA_functions))), y = score, fill = direction ), color = "black", stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values=c("#4D99C6", "#DC6E58")) +
  coord_flip() +
  theme_bw() +
  ylim(c(-6,6.5)) +
  theme(panel.grid = element_blank(), panel.border = element_rect(linetype = "solid", fill = NA, size = 1), axis.ticks = element_line(size = 1), axis.ticks.length = unit(0.2, "cm"), axis.title.y=element_blank())
dev.off()


#Correlation of mRNA-seq (expression) vs TE
par(mfrow=c(3,6))
for(i in 2:ncol(common_mRNA)) {
  plot(x=log10(common_mRNA[,i]), y=organ_TE[,i], type="p", pch=16, cex=0.5, col=alpha("black",0.5), ylab=colnames(organ_TE)[i])
  
}





