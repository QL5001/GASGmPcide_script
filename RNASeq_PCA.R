#Principal component analysis (PCA) on the RNASeq data of Streptococcus pyogenes HSC5 strain for determining its transcriptomic changes under the treatment of antimicrobial compounds.
#Begin#
library(knitr)
library(mixOmics)
library(dplyr)
library(stats)

RS_GP1 = read.csv(file = file.choose(), sep = ",", header = TRUE, stringsAsFactors = FALSE)    #open file "Differential_Expression_PCA" consisting of the differential expression profiel of S. pyogenes HSC5 strain under the treatment of antimicrobial compounds

X = RS_GP1[, 1:1764]  
X_mean = apply(X, 2, mean)
X_center = scale(X, center = X_mean, scale = FALSE)
Y = RS_GP1$Sample

pca.RS = pca(X_center, ncomp = 6, center = FALSE, scale = FALSE)
plot(pca.RS, ylim = c(0, 0.8), main = "PCA on WGS")

ev = pca.RS$prop_expl_var   #explained_variance
PCA_ev = data.frame(ev)
write.csv(PCA_ev, "PCA_ev.csv")  #export explained variance for each principal component into a csv file

sco = pca.RS$variates   #prinicipal component values
sco_PC = data.frame(SampleNo = RS_GP1$Sample, sco)
write.csv(sco_PC, "PCA_PC.csv")  #export the principal component values into a csv file

loa = data.frame(pca.RS$loadings)   #loading values
gene_name = colnames(X)   #assign gene names
loa_PC1 = data.frame(gene = gene_name, PC1loa = loa[, 1])
write.csv(loa_PC1, "PCA_loa_PC1.csv")   #export loading values associated with PC1 into a csv file
loa_PC2 = data.frame(gene = gene_name, PC2loa = loa[, 2])
write.csv(loa_PC2, "PCA_loa_PC2.csv")   #export loading values associated with PC2 into a csv file

plotIndiv(pca.RS, comp = c(1,2), group = RS_GP1$Sample, pch = c(15, 16), cex = 3, col = rainbow(2), ellipse = FALSE, ind.names = FALSE, legend = TRUE, size.legend = 15, size.xlabel = 18, size.ylabel = 18,  size.axis = 15,  title = "RS_PCA")

#Exported csv files will be used for further data analyses and drawing figures in Excel, Prism, and R programming
#End#
