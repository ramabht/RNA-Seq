RNASeq
================
RP Bhatia
17 November 2022

# Data analysis for -- Bhatia, R.P., Kirit, H.A., Predeus, A.V. et al.

# Transcriptomic profiling of Escherichia coli K-12 in response to a compendium of

# stressors. Sci Rep 12, 8788 (2022). <https://doi.org/10.1038/s41598-022-12463-3>

# library load

``` r
library(DESeq2)
```

``` r
library(data.table)
```

``` r
library(ggplot2)
```

``` r
library(ggrepel)
library(knitr)
```

``` r
library(reshape2)
```

``` r
library(RColorBrewer)
library(WGCNA)
```

``` r
library(viridis)
```

``` r
library(hrbrthemes)
library(extrafont)
```

``` r
library(grid)
library(gridExtra)
```

``` r
library(scales)
```

``` r
library(forcats)
```

``` r
library(dplyr)
```

``` r
library(tidyr)
```

``` r
library(magrittr)
```

``` r
loadfonts()
```

# Differential expression

``` r
counts <- read.csv("MG1655_counts.csv", header=TRUE, sep=',', row.names = 1)
colData <- read.csv("meta_data.txt.csv", header=TRUE, sep = ',')
ann_table <- read.table("MG1655.annotated.counts.tsv", header=T, check.names=F)
ann <- ann_table[,c(1:3)]
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=colData,
                              design = ~ Condition)
```

``` r
dds
```

``` r
dds <- DESeq(dds)
```

``` r
DE  <- list()
pairs <- fread("Contrasts.txt")
N <- nrow(pairs)
for (i in 1:N) {
  column   <- pairs$cond_col[i]
  cond1 <- pairs$cond1[i]
  cond2 <- pairs$cond2[i]
  DE[[i]]  <- results(dds,contrast=c(column,cond2,cond1))
  res <- as.data.frame(DE[[i]])
  baseMean   <- round(res$baseMean,digits=0)
  log2FC     <- round(res$log2FoldChange,digits=2)
  padj       <- formatC(res$padj, format="e", digits=2)
  lfcSE      <- round(res$lfcSE, digits=2)
  pvalue     <- formatC(res$pvalue, format="e", digits=2)
  stat       <- round(res$stat, digits=2)
  df <- data.frame(ann,baseMean,log2FC,padj,pvalue,lfcSE,stat)
  filename   <- paste(cond1,"_vs_",cond2,".deseq2.tsv",sep="")
  write.table(df[order(-log2FC), ], filename,
              quote=F,sep="\t", row.names = F)
}
```

# Volcano plots

``` r
# Rich M9 with Chloramphenicol
data1 <- read.csv('Rich_M9_ATC5_2_vs_Rich_M9_ATC5_CAM1.2.deseq2.tsv.csv')
data1$padj <- as.numeric(data1$padj)
```

``` r
data1$Diffexpressed <- "NS"
data1$Diffexpressed[data1$log2FC > 2 & data1$padj < 0.001] <- "UP"
data1$Diffexpressed[data1$log2FC < -2 & data1$padj < 0.001] <- "DOWN"
mycolors <- c("salmon", "cyan4", "azure4")
names(mycolors) <- c("DOWN", "UP", "NS")

p1 <- ggplot(data=data1, aes(x=log2FC, y=-log10(padj), 
                       col=Diffexpressed)) + 
  geom_point(size = 1, shape = 16) +
  geom_hline(yintercept=-log10(0.001), lty = 'dotted', size = 0.8) +
  geom_vline(xintercept = c(-2, 2), lty = 'dotted', size = 0.8) +
  scale_colour_manual(values = mycolors) + 
  scale_x_continuous(limits = c(-15,15), 
                     breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
  geom_label_repel(aes(label=Label), show.legend = FALSE, size = 2, 
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.3, "lines"), fontface = "bold") + 
  theme_bw() +
  theme(text = element_text(family = "Helvetica", color = "black", size=9),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1.5),
        axis.title=element_text(face="bold")) +
  theme(legend.position = "none") +
  xlab(expression('Log'[2]*' fold change')) +
  ylab(expression('-log'[10]*' FDR')) +
  labs(title = "Rich M9 with chloramphenicol") +
  annotate("text", x = 10, y = 40,
           label = "upregulated", color = "cyan4" ,fontface = "italic", size=4) +
  annotate("text", x = -10, y = 40,
           label = "downregulated", color = "salmon" ,fontface = "italic", size=4)+
  annotate("text", x = 0, y = 40,
           label = "ns", color = "azure4" ,fontface = "italic", size=4)
ggsave("CAM.pdf", p1, dpi = 300, height = 5, width = 4)
```

``` r
# Rich M9 with low oxygen

data2 <- read.csv('Rich_M9_ATC5_2_vs_Rich_M9_ATC5_NOX.deseq2.tsv.csv')
data2$padj <- as.numeric(data2$padj)
```

``` r
data2$Diffexpressed <- "NS"
data2$Diffexpressed[data2$log2FC > 2 & data2$padj < 0.001] <- "UP"
data2$Diffexpressed[data2$log2FC < -2 & data2$padj < 0.001] <- "DOWN"
mycolors <- c("salmon", "cyan4", "azure4")
names(mycolors) <- c("DOWN", "UP", "NS")

p2 <- ggplot(data=data2, aes(x=log2FC, y=-log10(padj), 
                             col=Diffexpressed)) + 
  geom_point(size = 1, shape = 16) +
  geom_hline(yintercept=-log10(0.001), lty = 'dotted', size =0.8) +
  geom_vline(xintercept = c(-2, 2), lty = 'dotted', size = 0.8) +
  scale_colour_manual(values = mycolors) +
  scale_x_continuous(limits = c(-15,15), 
                     breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
  geom_label_repel(aes(label=Label), show.legend = FALSE, size = 2, 
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.3, "lines"), fontface = "bold") + 
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", color = "black", size=9),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size=1.5),
    axis.title=element_text(face="bold")) +
  xlab(expression('Log'[2]*' fold change')) +
  ylab(expression('-log'[10]*' FDR')) +
  theme(legend.position = "none") +
  labs(title = "Rich M9 with low oxygen") +
  annotate("text", x = 10, y = 52,
           label = "upregulated", color = "cyan4" ,fontface = "italic", size=4) +
  annotate("text", x = -10, y = 52,
           label = "downregulated", color = "salmon" ,fontface = "italic", size=4)+
  annotate("text", x = 0, y = 52,
           label = "ns", color = "azure4" ,fontface = "italic", size=4)
ggsave("LOX.pdf", p2, dpi = 300, height = 5, width = 4)
```

``` r
# Rich M9 pH5

data3 <- read.csv('Rich_M9_ATC5_2_vs_Rich_M9_ATC5_PH5.deseq2.tsv.csv')
data3$padj <- as.numeric(data3$padj)
```

``` r
data3$Diffexpressed <- "NS"
data3$Diffexpressed[data3$log2FC > 2 & data3$padj < 0.001] <- "UP"
data3$Diffexpressed[data3$log2FC < -2 & data3$padj < 0.001] <- "DOWN"
mycolors <- c("salmon", "cyan4", "azure4")
names(mycolors) <- c("DOWN", "UP", "NS")

p3 <- ggplot(data=data3, aes(x=log2FC, y=-log10(padj), 
                             col=Diffexpressed)) + 
  geom_point(size = 1, shape = 16) +
  geom_hline(yintercept=-log10(0.001), lty = 'dotted', size = 0.8) +
  geom_vline(xintercept = c(-2, 2), lty = 'dotted', size = 0.8) +
  scale_colour_manual(values = mycolors) + 
  scale_x_continuous(limits = c(-15,15), 
                     breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
  geom_label_repel(aes(label=Label), show.legend = FALSE, size = 2, 
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.3, "lines"), fontface = "bold") + 
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", color = "black", size=9),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size=1.5),
    axis.title=element_text(face="bold")) +
  xlab(expression('Log'[2]*' fold change')) +
  ylab(expression('-log'[10]*' FDR')) +
  theme(legend.position = "none") +
  labs(title = "Rich M9 pH5") +
  annotate("text", x = 10, y = 50,
           label = "upregulated", color = "cyan4" ,fontface = "italic", size=4) +
  annotate("text", x = -10, y = 50,
           label = "downregulated", color = "salmon" ,fontface = "italic", size=4)+
  annotate("text", x = 0, y = 50,
           label = "ns", color = "azure4" ,fontface = "italic", size=4)

ggsave("pH5.pdf", p3, dpi = 300, height = 5, width = 4)
```

``` r
# Rich M9 with Trimethoprim

data4 <- read.csv('Rich_M9_ATC5_2_vs_Rich_M9_ATC5_TMP.deseq2.tsv.csv')
data4$padj <- as.numeric(data4$padj)
```

``` r
data4$Diffexpressed <- "NS"
data4$Diffexpressed[data4$log2FC > 2 & data4$padj < 0.001] <- "UP"
data4$Diffexpressed[data4$log2FC < -2 & data4$padj < 0.001] <- "DOWN"
mycolors <- c("salmon", "cyan4", "azure4")
names(mycolors) <- c("DOWN", "UP", "NS")

p4 <- ggplot(data=data4, aes(x=log2FC, y=-log10(padj), 
                             col=Diffexpressed)) + 
  geom_point(size = 1, shape = 16) +
  geom_hline(yintercept=-log10(0.001), lty = 'dotted', size = 0.8) +
  geom_vline(xintercept = c(-2, 2), lty = 'dotted', size = 0.8) +
  scale_colour_manual(values = mycolors) + 
  scale_x_continuous(limits = c(-15,15), 
                     breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
  geom_label_repel(aes(label=Label), show.legend = FALSE, size = 2, 
                   box.padding = unit(0.35, "lines"),
                   point.padding = unit(0.3, "lines"), fontface = "bold") + 
  theme_bw() +
  theme(
    text = element_text(family = "Helvetica", color = "black", size=9),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size=1.5),
    axis.title=element_text(face="bold")) +
  xlab(expression('Log'[2]*' fold change')) +
  ylab(expression('-log'[10]*' FDR')) +
  theme(legend.position = "none") +
  labs(title = "Rich M9 with trimethoprim") +
  annotate("text", x = 10, y = 68,
           label = "upregulated", color = "cyan4" ,fontface = "italic",size=4) +
  annotate("text", x = -10, y = 68,
           label = "downregulated", color = "salmon" ,fontface = "italic",size=4)+
  annotate("text", x = 0, y = 68,
           label = "ns", color = "azure4" ,fontface = "italic",size=4)

ggsave("TMP.pdf", p4, dpi = 300, height = 5, width = 4)
```

# Weighted Gene Coexpression Analysis (WGCNA)

# Data Input and Transformation

``` r
counts <- read.csv("MG1655.TPM.csv", header=TRUE, sep=',', row.names = 1)
colData <- read.csv("meta_data.txt.csv", header=TRUE, sep = ',')
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=colData,
                              design = ~ Condition)
```

``` r
dds
```

    ## class: DESeqDataSet 
    ## dim: 4566 20 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(4566): thrL thrA ... yqiD ysdD
    ## rowData names(0):
    ## colnames(20): Poor_M9_ATC5.rep1 Poor_M9_ATC5.rep2 ...
    ##   Rich_M9_No_ATC_2.rep1 Rich_M9_No_ATC_2.rep2
    ## colData names(2): Sample Condition

``` r
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds2 <- vst(dds, blind = T)
datExpr <- assay(dds2)
dim(datExpr)
```

``` r
datExpr0 <- as.data.frame(t(datExpr))
fix(datExpr0)
```

# Check Data for missing values

``` r
gsg = goodSamplesGenes(datExpr0, verbose = 3)
```

``` r
gsg$allOK
```

``` r
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
dim(datExpr0)
```

# Cluster to look for outliers

``` r
sampleTree = hclust(dist(datExpr0), method = "average")
sizeGrWindow(20,10)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
fix(datExpr0)
```

# Pick soft threshold

``` r
enableWGCNAThreads(nThreads = 2)
```

``` r
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5,
                        networkType = "signed")
```

    ## pickSoftThreshold: will use block size 4269.
    ##  pickSoftThreshold: calculating connectivity for given powers...
    ##    ..working on genes 1 through 4269 of 4269
    ##    Power SFT.R.sq slope truncated.R.sq mean.k. median.k. max.k.
    ## 1      1    0.261  3.44          0.744  2140.0    2140.0   2270
    ## 2      2    0.340 -3.76          0.416  1240.0    1210.0   1500
    ## 3      3    0.159 -1.39          0.662   788.0     766.0   1110
    ## 4      4    0.209 -1.06          0.821   535.0     518.0    874
    ## 5      5    0.326 -1.01          0.894   383.0     366.0    717
    ## 6      6    0.431 -1.02          0.941   284.0     267.0    608
    ## 7      7    0.534 -1.08          0.948   218.0     200.0    526
    ## 8      8    0.611 -1.11          0.962   172.0     154.0    462
    ## 9      9    0.681 -1.17          0.970   138.0     120.0    411
    ## 10    10    0.747 -1.24          0.975   113.0      95.7    369
    ## 11    12    0.792 -1.33          0.964    79.0      62.7    305
    ## 12    14    0.831 -1.41          0.965    57.8      42.7    258
    ## 13    16    0.836 -1.47          0.951    43.9      30.3    221
    ## 14    18    0.860 -1.49          0.956    34.2      22.1    193
    ## 15    20    0.868 -1.52          0.958    27.3      16.5    169

``` r
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.80,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```

# Network construction

``` r
#STEPWISE NETWORK CONSTRUCTION
softPower = 14
adjacency = adjacency(datExpr0, power = softPower, type = "signed")
TOM = TOMsimilarity(adjacency)
```

    ## ..connectivity..
    ## ..matrix multiplication (system BLAS)..
    ## ..normalization..
    ## ..done.

``` r
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average")
sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
minModuleSize = 30

# Module identification using dynamic tree cut
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
```

    ##  ..cutHeight not given, setting it to 0.992  ===>  99% of the (truncated) height range in dendro.
    ##  ..done.

``` r
table(dynamicMods)
```

    ## dynamicMods
    ##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
    ## 434 378 313 241 231 210 191 173 151 149 129 128 125 124 114 109 100  95  95  86 
    ##  21  22  23  24  25  26  27  28  29  30  31  32 
    ##  84  84  78  74  61  57  51  44  41  41  41  37

``` r
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
```

    ## dynamicColors
    ##         black          blue         brown          cyan     darkgreen 
    ##           191           378           313           124            84 
    ##      darkgrey    darkorange       darkred darkturquoise         green 
    ##            74            57            84            78           231 
    ##   greenyellow        grey60     lightcyan    lightgreen   lightyellow 
    ##           129           100           109            95            95 
    ##       magenta  midnightblue        orange paleturquoise          pink 
    ##           151           114            61            41           173 
    ##        purple           red     royalblue   saddlebrown        salmon 
    ##           149           210            86            41           125 
    ##       skyblue     steelblue           tan     turquoise        violet 
    ##            44            41           128           434            37 
    ##         white        yellow 
    ##            51           241

``` r
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
fix(MEs)

# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs)

# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")

# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25 #0.25 in previous run

# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")

# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, 
                          cutHeight = MEDissThres, verbose = 3)
```

    ##  mergeCloseModules: Merging modules whose distance is less than 0.25
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 32 module eigengenes in given set.
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 21 module eigengenes in given set.
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 20 module eigengenes in given set.
    ##    Calculating new MEs...
    ##    multiSetMEs: Calculating module MEs.
    ##      Working on set 1 ...
    ##      moduleEigengenes: Calculating 20 module eigengenes in given set.

``` r
# The merged module colors
mergedColors = merge$colors

# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs
sizeGrWindow(12, 9)

#Plot the dendogram with merged and merged modules
dendro = plotDendroAndColors(geneTree, cbind(dynamicColors, 
          mergedColors), c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
table(moduleColors)
```

    ## moduleColors
    ##         black          blue          cyan     darkgreen      darkgrey 
    ##           191           691           202           209            74 
    ##    darkorange       darkred         green   greenyellow    lightgreen 
    ##           230           503           472           563           181 
    ##   lightyellow       magenta  midnightblue        orange paleturquoise 
    ##           146           151           155            61            41 
    ##        purple   saddlebrown       skyblue           tan        violet 
    ##           149            41            44           128            37

``` r
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1

MEs = mergedMEs
fix(MEs)
write.csv(MEs) # add file path
```


# Module-Trait Relationship

``` r
meta_data <- read.csv("meta_data.csv")
head(meta_data)
```

    ##                Sample Lennox.Broth Poor.M9 Rich.M9 CAM TMP LOX ph5
    ## 1   Poor_M9_ATC5.rep1            0       1       0   0   0   0   0
    ## 2   Poor_M9_ATC5.rep2            0       1       0   0   0   0   0
    ## 3  Rich_LB_ATC12.rep1            1       0       0   0   0   0   0
    ## 4  Rich_LB_ATC12.rep2            1       0       0   0   0   0   0
    ## 5 Rich_M9_ATC5_1.rep1            0       0       1   0   0   0   0
    ## 6 Rich_M9_ATC5_1.rep2            0       0       1   0   0   0   0

``` r
traits <- meta_data[,c(2:8)]
head(traits)
```

    ##   Lennox.Broth Poor.M9 Rich.M9 CAM TMP LOX ph5
    ## 1            0       1       0   0   0   0   0
    ## 2            0       1       0   0   0   0   0
    ## 3            1       0       0   0   0   0   0
    ## 4            1       0       0   0   0   0   0
    ## 5            0       0       1   0   0   0   0
    ## 6            0       0       1   0   0   0   0

``` r
rownames(traits) <- rownames(datExpr0)
rownames(traits)==rownames(datExpr0)
```

    ##  [1] TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE TRUE
    ## [16] TRUE TRUE TRUE TRUE TRUE

``` r
fix(traits)
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPadj = p.adjust(moduleTraitPvalue, method = "fdr")
fix(moduleTraitPadj)
sizeGrWindow(10,6)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPadj, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

#bottom, left, top, right margins
par(mar = c(6, 12, 3, 2))

# Display the correlation values within a heatmap plot
pal <- colorspace::diverge_hcl(3,l= c(50,90), c = 100, power = 1)
coul <- colorRampPalette(pal)(25)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = coul,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1, 1),
               main = paste("Module-trait relationships"))
```

# Functional enrichment analysis plots (code is running for darkred and green modules)

``` r
mydat <- read.csv("bubbleplot_data1.csv", header = T)
mydat$GO_term <- factor(mydat$GO_term, 
                        levels = unique(mydat$GO_term[order(mydat$Module)]))
#First we group by Type and GO_term, and assign a "yes" to the first row
#and "no" to every other row of the grouping
mydat %<>% 
  group_by(Module, GO_term) %>%
  mutate(typefill = if_else(row_number() == 1, "yes", "no")) %>%
  ungroup()
#Then in the whole data.frame, typefill = "yes" will be replaced by the Type value
#from that row, and typefill = "no" will be replaced with NA
mydat %<>% mutate(typefill = ifelse(typefill == "yes", 
                                    as.character(Module), NA))

ggplot(mydat, aes(y = factor(GO_term), x = Gene_ratio)) +
  geom_tile(aes(width = Inf, fill = typefill), size = 0.01,
            alpha = 0.5) +
  geom_point(aes(size = Count), alpha = 0.5) +
  scale_fill_manual(values = c("darkred", "green4")) +
  theme_ipsum(axis_title_size = 14,
        axis_text_size = 11,
        axis_title_face = "bold")+
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(colour="black"),
        legend.text=element_text(size=12)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_size(range = c(.1, 15), name="Gene count") +
  ylab('Enriched GO terms') + 
  xlab('Gene Ratio') +
  guides(fill=guide_legend(title="Co-expression Modules", order = 1),
         size = guide_legend(order = 2))
```

![](rna-seq_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
ggsave("GO_enriched_terms_darkred_green.pdf", dpi = 300, height = 25,
       width = 30, units = "cm")
```

# Bar plot for KEGG Terms

``` r
data <- read.csv('kegg_plot_all.csv', header=T)
data$Description <- factor(data$Description, levels = unique(data$Description[order(data$Module)]))

ggplot(data, aes(x=Count, y=Description, fill=Module))+
  geom_bar(stat="identity", width = 0.8, alpha = 0.5) +
  scale_fill_manual(values=c("blue", "darkorange1", "darkred", "green4", "purple")) +
  theme_ipsum(axis_title_size = 14,
              axis_text_size = 12, 
              axis_title_face = "bold") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  theme(axis.line = element_line(colour = "black"),
        axis.text=element_text(colour="black"),
        legend.text=element_text(size=12)) +
  ylab('Enriched KEGG terms') + 
  xlab('Gene Count')
```

![](rna-seq_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
ggsave("KEGG_barplot.pdf", device = "pdf", dpi = 300, height = 24, 
       width = 30, units = "cm")
```
