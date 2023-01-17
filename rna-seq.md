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

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: GenomicRanges

    ## Loading required package: GenomeInfoDb

    ## Loading required package: SummarizedExperiment

    ## Loading required package: Biobase

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: DelayedArray

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'DelayedArray'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colMaxs, colMins, colRanges, rowMaxs, rowMins, rowRanges

    ## The following objects are masked from 'package:base':
    ## 
    ##     aperm, apply, rowsum

``` r
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following object is masked from 'package:SummarizedExperiment':
    ## 
    ##     shift

    ## The following object is masked from 'package:GenomicRanges':
    ## 
    ##     shift

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     shift

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, second

``` r
library(ggplot2)
```

    ## Warning: package 'ggplot2' was built under R version 4.0.4

``` r
library(ggrepel)
library(knitr)
```

    ## Warning: package 'knitr' was built under R version 4.0.5

``` r
library(reshape2)
```

    ## 
    ## Attaching package: 'reshape2'

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     dcast, melt

``` r
library(RColorBrewer)
library(WGCNA)
```

    ## Loading required package: dynamicTreeCut

    ## Loading required package: fastcluster

    ## 
    ## Attaching package: 'fastcluster'

    ## The following object is masked from 'package:stats':
    ## 
    ##     hclust

    ## 

    ## 
    ## Attaching package: 'WGCNA'

    ## The following object is masked from 'package:IRanges':
    ## 
    ##     cor

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     cor

    ## The following object is masked from 'package:stats':
    ## 
    ##     cor

``` r
library(viridis)
```

    ## Loading required package: viridisLite

``` r
library(hrbrthemes)
library(extrafont)
```

    ## Registering fonts with R

``` r
library(grid)
library(gridExtra)
```

    ## 
    ## Attaching package: 'gridExtra'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     combine

``` r
library(scales)
```

    ## 
    ## Attaching package: 'scales'

    ## The following object is masked from 'package:viridis':
    ## 
    ##     viridis_pal

``` r
library(forcats)
```

    ## Warning: package 'forcats' was built under R version 4.0.3

``` r
library(dplyr)
```

    ## Warning: package 'dplyr' was built under R version 4.0.3

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:gridExtra':
    ## 
    ##     combine

    ## The following objects are masked from 'package:data.table':
    ## 
    ##     between, first, last

    ## The following object is masked from 'package:matrixStats':
    ## 
    ##     count

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:GenomicRanges':
    ## 
    ##     intersect, setdiff, union

    ## The following object is masked from 'package:GenomeInfoDb':
    ## 
    ##     intersect

    ## The following objects are masked from 'package:IRanges':
    ## 
    ##     collapse, desc, intersect, setdiff, slice, union

    ## The following objects are masked from 'package:S4Vectors':
    ## 
    ##     first, intersect, rename, setdiff, setequal, union

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
```

    ## 
    ## Attaching package: 'tidyr'

    ## The following object is masked from 'package:reshape2':
    ## 
    ##     smiths

    ## The following object is masked from 'package:S4Vectors':
    ## 
    ##     expand

``` r
library(magrittr)
```

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

``` r
loadfonts()
```

    ## Agency FB already registered with pdfFonts().

    ## Algerian already registered with pdfFonts().

    ## Arial Black already registered with pdfFonts().

    ## Arial already registered with pdfFonts().

    ## Arial Narrow already registered with pdfFonts().

    ## Arial Rounded MT Bold already registered with pdfFonts().

    ## Bahnschrift already registered with pdfFonts().

    ## Baskerville Old Face already registered with pdfFonts().

    ## Bauhaus 93 already registered with pdfFonts().

    ## Bell MT already registered with pdfFonts().

    ## Berlin Sans FB already registered with pdfFonts().

    ## No regular (non-bold, non-italic) version of Berlin Sans FB Demi. Skipping setup for this font.

    ## Bernard MT Condensed already registered with pdfFonts().

    ## Blackadder ITC already registered with pdfFonts().

    ## No regular (non-bold, non-italic) version of Bodoni Bd BT. Skipping setup for this font.

    ## Bodoni Bk BT already registered with pdfFonts().

    ## Bodoni MT already registered with pdfFonts().

    ## Bodoni MT Black already registered with pdfFonts().

    ## Bodoni MT Condensed already registered with pdfFonts().

    ## Bodoni MT Poster Compressed already registered with pdfFonts().

    ## Book Antiqua already registered with pdfFonts().

    ## Bookman Old Style already registered with pdfFonts().

    ## Bookshelf Symbol 7 already registered with pdfFonts().

    ## Bradley Hand ITC already registered with pdfFonts().

    ## Britannic Bold already registered with pdfFonts().

    ## Broadway already registered with pdfFonts().

    ## No regular (non-bold, non-italic) version of Brush Script MT. Skipping setup for this font.

    ## Calibri already registered with pdfFonts().

    ## Calibri Light already registered with pdfFonts().

    ## Californian FB already registered with pdfFonts().

    ## Calisto MT already registered with pdfFonts().

    ## No regular (non-bold, non-italic) version of Cambria. Skipping setup for this font.

    ## Candara already registered with pdfFonts().

    ## Candara Light already registered with pdfFonts().

    ## Castellar already registered with pdfFonts().

    ## Centaur already registered with pdfFonts().

    ## Century already registered with pdfFonts().

    ## Century725 Cn BT already registered with pdfFonts().

    ## Century751 No2 BT already registered with pdfFonts().

    ## Century751 BT already registered with pdfFonts().

    ## No regular (non-bold, non-italic) version of Century751 SeBd BT. Skipping setup for this font.

    ## Century Gothic already registered with pdfFonts().

    ## CentSchbkCyrill BT already registered with pdfFonts().

    ## Century Schoolbook already registered with pdfFonts().

    ## Chiller already registered with pdfFonts().

    ## Clarendon Blk BT already registered with pdfFonts().

    ## Clarendon BT already registered with pdfFonts().

    ## Clarendon Lt BT already registered with pdfFonts().

    ## Colonna MT already registered with pdfFonts().

    ## Comic Sans MS already registered with pdfFonts().

    ## Consolas already registered with pdfFonts().

    ## Constantia already registered with pdfFonts().

    ## Cooper Black already registered with pdfFonts().

    ## Copperplate Gothic Bold already registered with pdfFonts().

    ## Copperplate Gothic Light already registered with pdfFonts().

    ## Corbel already registered with pdfFonts().

    ## Corbel Light already registered with pdfFonts().

    ## Courier New already registered with pdfFonts().

    ## Curlz MT already registered with pdfFonts().

    ## DeVinne Txt BT already registered with pdfFonts().

    ## DFGothic-EB already registered with pdfFonts().

    ## DFKaiSho-SB already registered with pdfFonts().

    ## DFMincho-SU already registered with pdfFonts().

    ## DFMincho-UB already registered with pdfFonts().

    ## DFMincho-W5 already registered with pdfFonts().

    ## DFPOP1-W9 already registered with pdfFonts().

    ## Dubai already registered with pdfFonts().

    ## Dubai Light already registered with pdfFonts().

    ## Dubai Medium already registered with pdfFonts().

    ## Ebrima already registered with pdfFonts().

    ## Edwardian Script ITC already registered with pdfFonts().

    ## Elephant already registered with pdfFonts().

    ## Embassy BT already registered with pdfFonts().

    ## EngraversGothic BT already registered with pdfFonts().

    ## Engravers MT already registered with pdfFonts().

    ## Eras Bold ITC already registered with pdfFonts().

    ## Eras Demi ITC already registered with pdfFonts().

    ## Eras Light ITC already registered with pdfFonts().

    ## Eras Medium ITC already registered with pdfFonts().

    ## No regular (non-bold, non-italic) version of Exotc350 Bd BT. Skipping setup for this font.

    ## No regular (non-bold, non-italic) version of Exotc350 DmBd BT. Skipping setup for this font.

    ## Felix Titling already registered with pdfFonts().

    ## Footlight MT Light already registered with pdfFonts().

    ## Forte already registered with pdfFonts().

    ## Franklin Gothic Book already registered with pdfFonts().

    ## Franklin Gothic Demi already registered with pdfFonts().

    ## Franklin Gothic Demi Cond already registered with pdfFonts().

    ## Franklin Gothic Heavy already registered with pdfFonts().

    ## Franklin Gothic Medium already registered with pdfFonts().

    ## Franklin Gothic Medium Cond already registered with pdfFonts().

    ## Freehand521 BT already registered with pdfFonts().

    ## Freestyle Script already registered with pdfFonts().

    ## French Script MT already registered with pdfFonts().

    ## Futura Md BT already registered with pdfFonts().

    ## Futura Bk BT already registered with pdfFonts().

    ## Gabriola already registered with pdfFonts().

    ## Gadugi already registered with pdfFonts().

    ## Garamond already registered with pdfFonts().

    ## More than one version of regular/bold/italic found for Geometr212 BkCn BT. Skipping setup for this font.

    ## Geometr415 Blk BT already registered with pdfFonts().

    ## Geometr706 BlkCn BT already registered with pdfFonts().

    ## GeoSlab703 Md BT already registered with pdfFonts().

    ## GeoSlab703 MdCn BT already registered with pdfFonts().

    ## Georgia already registered with pdfFonts().

    ## Gigi already registered with pdfFonts().

    ## Gill Sans Ultra Bold already registered with pdfFonts().

    ## Gill Sans Ultra Bold Condensed already registered with pdfFonts().

    ## Gill Sans MT already registered with pdfFonts().

    ## Gill Sans MT Condensed already registered with pdfFonts().

    ## Gill Sans MT Ext Condensed Bold already registered with pdfFonts().

    ## Gloucester MT Extra Condensed already registered with pdfFonts().

    ## Goudy Old Style already registered with pdfFonts().

    ## Goudy Stout already registered with pdfFonts().

    ## Haettenschweiler already registered with pdfFonts().

    ## No regular (non-bold, non-italic) version of Harlow Solid Italic. Skipping setup for this font.

    ## Harrington already registered with pdfFonts().

    ## High Tower Text already registered with pdfFonts().

    ## HoloLens MDL2 Assets already registered with pdfFonts().

    ## Humanst521 BT already registered with pdfFonts().

    ## Humanst521 Lt BT already registered with pdfFonts().

    ## Humnst777 Blk BT already registered with pdfFonts().

    ## Humnst777 BlkCn BT already registered with pdfFonts().

    ## Humnst777 Cn BT already registered with pdfFonts().

    ## Humnst777 Lt BT already registered with pdfFonts().

    ## Humnst777 BT already registered with pdfFonts().

    ## Impact already registered with pdfFonts().

    ## Imprint MT Shadow already registered with pdfFonts().

    ## Informal Roman already registered with pdfFonts().

    ## Ink Free already registered with pdfFonts().

    ## Javanese Text already registered with pdfFonts().

    ## Jokerman already registered with pdfFonts().

    ## Juice ITC already registered with pdfFonts().

    ## Kaufmann BT already registered with pdfFonts().

    ## Kristen ITC already registered with pdfFonts().

    ## Kunstler Script already registered with pdfFonts().

    ## Wide Latin already registered with pdfFonts().

    ## Leelawadee already registered with pdfFonts().

    ## Leelawadee UI already registered with pdfFonts().

    ## Leelawadee UI Semilight already registered with pdfFonts().

    ## More than one version of regular/bold/italic found for Lucida Bright. Skipping setup for this font.

    ## No regular (non-bold, non-italic) version of Lucida Calligraphy. Skipping setup for this font.

    ## Lucida Console already registered with pdfFonts().

    ## More than one version of regular/bold/italic found for Lucida Fax. Skipping setup for this font.

    ## No regular (non-bold, non-italic) version of Lucida Handwriting. Skipping setup for this font.

    ## More than one version of regular/bold/italic found for Lucida Sans. Skipping setup for this font.

    ## Lucida Sans Typewriter already registered with pdfFonts().

    ## Lucida Sans Unicode already registered with pdfFonts().

    ## No regular (non-bold, non-italic) version of Magneto. Skipping setup for this font.

    ## Maiandra GD already registered with pdfFonts().

    ## Malgun Gothic already registered with pdfFonts().

    ## Malgun Gothic Semilight already registered with pdfFonts().

    ## Marlett already registered with pdfFonts().

    ## Matura MT Script Capitals already registered with pdfFonts().

    ## Microsoft Himalaya already registered with pdfFonts().

    ## Microsoft Yi Baiti already registered with pdfFonts().

    ## Microsoft New Tai Lue already registered with pdfFonts().

    ## Microsoft PhagsPa already registered with pdfFonts().

    ## Microsoft Sans Serif already registered with pdfFonts().

    ## Microsoft Tai Le already registered with pdfFonts().

    ## Microsoft Uighur already registered with pdfFonts().

    ## Mistral already registered with pdfFonts().

    ## Modern No. 20 already registered with pdfFonts().

    ## Mongolian Baiti already registered with pdfFonts().

    ## Monotype Corsiva already registered with pdfFonts().

    ## MS Outlook already registered with pdfFonts().

    ## MS Reference Sans Serif already registered with pdfFonts().

    ## MS Reference Specialty already registered with pdfFonts().

    ## MT Extra already registered with pdfFonts().

    ## MV Boli already registered with pdfFonts().

    ## Myanmar Text already registered with pdfFonts().

    ## No regular (non-bold, non-italic) version of News701 BT. Skipping setup for this font.

    ## No regular (non-bold, non-italic) version of News706 BT. Skipping setup for this font.

    ## NewsGoth BT already registered with pdfFonts().

    ## NewsGoth Lt BT already registered with pdfFonts().

    ## NewsGoth Cn BT already registered with pdfFonts().

    ## Niagara Engraved already registered with pdfFonts().

    ## Niagara Solid already registered with pdfFonts().

    ## Nirmala UI already registered with pdfFonts().

    ## Nirmala UI Semilight already registered with pdfFonts().

    ## OCR-A BT already registered with pdfFonts().

    ## OCR A Extended already registered with pdfFonts().

    ## OCR-B 10 BT already registered with pdfFonts().

    ## Old English Text MT already registered with pdfFonts().

    ## Onyx already registered with pdfFonts().

    ## Palace Script MT already registered with pdfFonts().

    ## Palatino Linotype already registered with pdfFonts().

    ## Papyrus already registered with pdfFonts().

    ## Parchment already registered with pdfFonts().

    ## Perpetua already registered with pdfFonts().

    ## Perpetua Titling MT already registered with pdfFonts().

    ## Playbill already registered with pdfFonts().

    ## Poor Richard already registered with pdfFonts().

    ## Pristina already registered with pdfFonts().

    ## No regular (non-bold, non-italic) version of Rage Italic. Skipping setup for this font.

    ## Ravie already registered with pdfFonts().

    ## Rockwell already registered with pdfFonts().

    ## Rockwell Condensed already registered with pdfFonts().

    ## Rockwell Extra Bold already registered with pdfFonts().

    ## Schadow BT already registered with pdfFonts().

    ## Script MT Bold already registered with pdfFonts().

    ## Segoe MDL2 Assets already registered with pdfFonts().

    ## Segoe Print already registered with pdfFonts().

    ## Segoe Script already registered with pdfFonts().

    ## Segoe UI already registered with pdfFonts().

    ## Segoe UI Light already registered with pdfFonts().

    ## Segoe UI Semibold already registered with pdfFonts().

    ## Segoe UI Semilight already registered with pdfFonts().

    ## Segoe UI Black already registered with pdfFonts().

    ## Segoe UI Emoji already registered with pdfFonts().

    ## Segoe UI Historic already registered with pdfFonts().

    ## Segoe UI Symbol already registered with pdfFonts().

    ## Showcard Gothic already registered with pdfFonts().

    ## SimSun-ExtB already registered with pdfFonts().

    ## Snap ITC already registered with pdfFonts().

    ## Square721 BT already registered with pdfFonts().

    ## Square721 Cn BT already registered with pdfFonts().

    ## Stencil already registered with pdfFonts().

    ## Swis721 Blk BT already registered with pdfFonts().

    ## Swis721 BlkCn BT already registered with pdfFonts().

    ## More than one version of regular/bold/italic found for Swis721 WGL4 BT. Skipping setup for this font.

    ## More than one version of regular/bold/italic found for Swis721 BT. Skipping setup for this font.

    ## Swis721 Cn BT already registered with pdfFonts().

    ## Swis721 Hv BT already registered with pdfFonts().

    ## Swis721 Lt BT already registered with pdfFonts().

    ## Swis721 LtEx BT already registered with pdfFonts().

    ## Sylfaen already registered with pdfFonts().

    ## Symbol already registered with pdfFonts().

    ## Tahoma already registered with pdfFonts().

    ## Tempus Sans ITC already registered with pdfFonts().

    ## Times New Roman already registered with pdfFonts().

    ## Trebuchet MS already registered with pdfFonts().

    ## Tw Cen MT already registered with pdfFonts().

    ## Tw Cen MT Condensed already registered with pdfFonts().

    ## Tw Cen MT Condensed Extra Bold already registered with pdfFonts().

    ## TypoUpright BT already registered with pdfFonts().

    ## Verdana already registered with pdfFonts().

    ## Viner Hand ITC already registered with pdfFonts().

    ## No regular (non-bold, non-italic) version of Vivaldi. Skipping setup for this font.

    ## Vladimir Script already registered with pdfFonts().

    ## Webdings already registered with pdfFonts().

    ## Wingdings already registered with pdfFonts().

    ## Wingdings 2 already registered with pdfFonts().

    ## Wingdings 3 already registered with pdfFonts().

    ## ZWAdobeF already registered with pdfFonts().

# Differential expression

``` r
counts <- read.csv("MG1655_counts.csv", header=TRUE, sep=',', row.names = 1)
colData <- read.csv("meta_data.txt.csv", header=TRUE, sep = ',')
ann_table <- read.table("MG1655.annotated.counts.tsv", header=T, check.names=F)
ann <- ann_table[,c(1:3)]
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=colData,
                              design = ~ Condition)
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

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
dds <- DESeq(dds)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

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

    ## Warning: NAs introduced by coercion

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

    ## Warning: Removed 196 rows containing missing values (geom_point).

    ## Warning: Removed 196 rows containing missing values (geom_label_repel).

``` r
# Rich M9 with low oxygen

data2 <- read.csv('Rich_M9_ATC5_2_vs_Rich_M9_ATC5_NOX.deseq2.tsv.csv')
data2$padj <- as.numeric(data2$padj)
```

    ## Warning: NAs introduced by coercion

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

    ## Warning: Removed 368 rows containing missing values (geom_point).

    ## Warning: Removed 368 rows containing missing values (geom_label_repel).

``` r
# Rich M9 pH5

data3 <- read.csv('Rich_M9_ATC5_2_vs_Rich_M9_ATC5_PH5.deseq2.tsv.csv')
data3$padj <- as.numeric(data3$padj)
```

    ## Warning: NAs introduced by coercion

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

    ## Warning: Removed 282 rows containing missing values (geom_point).

    ## Warning: Removed 282 rows containing missing values (geom_label_repel).

``` r
# Rich M9 with Trimethoprim

data4 <- read.csv('Rich_M9_ATC5_2_vs_Rich_M9_ATC5_TMP.deseq2.tsv.csv')
data4$padj <- as.numeric(data4$padj)
```

    ## Warning: NAs introduced by coercion

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

    ## Warning: Removed 109 rows containing missing values (geom_point).

    ## Warning: Removed 109 rows containing missing values (geom_label_repel).

# Weighted Gene Coexpression Analysis (WGCNA)

# Data Input and Transformation

``` r
counts <- read.csv("MG1655.TPM.csv", header=TRUE, sep=',', row.names = 1)
colData <- read.csv("meta_data.txt.csv", header=TRUE, sep = ',')
dds <- DESeqDataSetFromMatrix(countData=round(counts), colData=colData,
                              design = ~ Condition)
```

    ## converting counts to integer mode

    ## Warning in DESeqDataSet(se, design = design, ignoreRank): some variables in
    ## design formula are characters, converting to factors

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

    ## [1] 4269   20

``` r
datExpr0 <- as.data.frame(t(datExpr))
fix(datExpr0)
```

# Check Data for missing values

``` r
gsg = goodSamplesGenes(datExpr0, verbose = 3)
```

    ##  Flagging genes and samples with too many missing values...
    ##   ..step 1

``` r
gsg$allOK
```

    ## [1] TRUE

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

    ## [1]   20 4269

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

    ## Allowing parallel execution with up to 2 working processes.

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

    ## "","MEmagenta","MEpurple","MEviolet","MEdarkgrey","MEdarkred","MEskyblue","MEdarkgreen","MEgreen","MEblack","MElightyellow","MEsaddlebrown","MEdarkorange","MEtan","MEcyan","MEmidnightblue","MElightgreen","MEblue","MEgreenyellow","MEorange","MEpaleturquoise"
    ## "Poor_M9_ATC5.rep1",0.148023450299891,0.513763543886181,0.236168995791729,-0.437266510318311,-0.48045301197659,-0.122425218901012,-0.228282185958609,0.000914476722855968,-0.370971960629496,0.0461130177322646,-0.0811320489664395,-0.0498826204102785,-0.313344543206031,-0.376095038538396,-0.265689015620539,0.311873090932001,0.637356587125749,0.401583572542743,-0.0194475220002,0.0256434736364665
    ## "Poor_M9_ATC5.rep2",0.169493201118876,0.535861633141382,0.262270095047881,-0.467808633378297,-0.541451699899909,-0.104345538158782,-0.224191325261602,0.00242490765046779,-0.362995428128174,0.0139703919650327,-0.0521095854853681,-0.0486214123871205,-0.300596663631245,-0.396274286007914,-0.305507707685475,0.328597068828337,0.668394429233862,0.424305810650601,0.00312635548921592,0.0173771721891663
    ## "Rich_LB_ATC12.rep1",0.122591623802232,0.24240282347377,0.270981728263887,0.0605203032256452,0.25678635788645,0.201496191686585,0.398757506424917,0.516528082916718,0.160686584456651,-0.0897892636842958,-0.115813555717165,-0.50327630562187,-0.225348276534829,-0.300287821449602,-0.136160695686651,-0.319133841164133,-0.129632922756852,-0.389572752634236,-0.190380372384948,-0.194290875042254
    ## "Rich_LB_ATC12.rep2",0.193242367709555,0.246427418220796,0.283166429039667,0.102386440093965,0.324453682234328,0.262316033248041,0.430548520038785,0.515861663514874,0.100855648595244,-0.203432824234886,-0.0749342516503449,-0.523467826880861,-0.252813427703952,-0.288652676236914,-0.160629822785542,-0.330059602528283,-0.165457689905191,-0.495841722520448,-0.240363194691589,-0.150784580218306
    ## "Rich_M9_ATC5_1.rep1",0.13909722868439,-0.0732989304382445,-0.41540763128061,-0.338532706351172,0.0382347992661086,-0.112631654036493,-0.387895552488464,-0.0930920064685771,0.536798429494338,0.644132913850015,-0.193217235252406,-0.0051757155380644,-0.0654287647677226,0.299019755585199,-0.0512325063587671,0.0486746083430298,-0.0643988229469377,0.0655760668779891,-0.161729316626695,-0.52699386141369
    ## "Rich_M9_ATC5_1.rep2",-0.24583251452153,-0.22682182860371,-0.0403136385874573,0.282100459098036,0.174176476105069,-0.151520086059796,0.0594207997611345,-0.17396897739642,-0.0549213849003893,-0.0548400125362267,0.000412174744022095,0.112066899133751,0.2047896943159,0.265666049392568,0.384290358180616,-0.26376358511916,-0.138528267763316,0.00494826835506209,-0.131245880061458,0.19111512014274
    ## "Rich_M9_ATC5_2.rep1",-0.356438768838969,-0.219910979815513,-0.0800233746062153,0.289866315267419,0.134319772367001,-0.122569406728857,0.0695173363896034,-0.184700997031859,-0.140883372986085,-0.122968691067525,0.0576104291892831,0.151997043427963,0.316512024036798,0.193553377238027,0.390328621354428,-0.238830628864843,-0.107310317852597,0.0313065790027887,-0.119042811291744,0.258607850529942
    ## "Rich_M9_ATC5_2.rep2",-0.323583179857974,-0.195018416472255,-0.0684694058754851,0.214124897153948,0.047293570544495,-0.135227377317559,0.0255408801677487,-0.1865770773611,-0.180812405547395,-0.0979365819398014,0.0228509381133877,0.12381996855194,0.300498342160827,0.157096151033232,0.381354339520477,-0.259696258952964,-0.0853119752506599,0.0756228297319332,-0.110755401846816,0.196215014860879
    ## "Rich_M9_ATC5_CAM1.2.rep1",0.248247383264992,-0.0965162209579781,0.0696244156815559,0.0304761694573496,0.109341554585466,0.0542430657551801,0.0990774757231679,-0.0646638129814952,-0.040580853808708,-0.144817586751518,-0.0802254712289454,0.104074598392292,-0.182880996364123,0.0516730556181205,-0.0876069267181285,0.183650760843953,-0.0620560346447901,-0.132031748462104,0.168708368198258,0.139631020654767
    ## "Rich_M9_ATC5_CAM1.2.rep2",0.230654558635825,-0.0977014603775222,0.0457398990051597,0.0597407847156936,0.117552336339167,0.0603124044020759,0.0898983030636187,-0.0632910026837467,-0.0345749336066853,-0.126764883342108,-0.0773836818939145,0.111628227814705,-0.170565992884838,0.0593702889164918,-0.104588809808462,0.180069199086458,-0.0702498610484946,-0.127515254617303,0.190982602730283,0.131697469573614
    ## "Rich_M9_ATC5_NOX.rep1",-0.348134469318627,-0.133237091844213,0.0158595204199617,0.0847831005942691,-0.0414117719908144,-0.194355756241287,-0.0474639691542433,-0.120670732991466,-0.0107476366151399,-0.0832980096863974,0.650300223561675,0.140794117484287,0.0904416588211077,0.0808747198049944,0.170190467931642,0.0648758049061078,-0.0496282203302242,0.0463157013688738,0.123362949406381,0.0842528328914855
    ## "Rich_M9_ATC5_NOX.rep2",-0.378906760908176,-0.135930144838909,0.0405229298565349,0.0314752091048196,-0.0558074554588377,-0.220888576324285,-0.0383077360221135,-0.100514333688627,-0.0274601190926278,-0.0968217910573332,0.64621721116651,0.1248145445106,0.093171493185592,0.0720083908001592,0.153101321674046,0.0484201977649033,-0.0384859167199985,0.0507348056536244,0.160675911001989,0.0720150347976436
    ## "Rich_M9_ATC5_PH5.rep1",0.200479534141156,-0.0854486247298491,0.112909834977734,-0.110233116230937,-0.0532924191736733,-0.206156525876924,-0.133677888783897,-0.0917027997988925,0.006228794386817,0.0144013035638299,-0.0384216999896679,0.074907786696628,-0.210290987283898,0.0482976741797717,-0.0866378975722963,0.105869927224938,0.0160160328687548,0.104011714727614,-0.153791383460545,-0.114594401031536
    ## "Rich_M9_ATC5_PH5.rep2",0.197687561636876,0.254706381513786,0.240614607621414,0.155981880777785,0.297446256940531,0.311612811362021,0.412042764524221,0.517621278302878,0.122828842362881,-0.146767863192822,-0.0730552299032978,-0.500499754887412,-0.242206858760017,-0.27715012306956,-0.128687842973141,-0.302122023406132,-0.146978591636923,-0.405991717828903,-0.169597510116299,-0.110130909208937
    ## "Rich_M9_ATC5_TMP.rep1",0.178945401979597,-0.0825096253301265,0.0967481654800135,0.150574414999873,-0.194152590759726,-0.0222354918751077,0.0218050252724963,-0.0799114805850314,-0.117446919247673,-0.0583698892768351,-0.0976357957238723,0.170088684262542,0.0557840678382624,-0.139010519700196,-0.184893372207999,0.148591948546694,0.0151057807787496,0.117144254708935,0.557255286873427,0.229851019657892
    ## "Rich_M9_ATC5_TMP.rep2",0.166586309749613,-0.0680233137626465,0.0840971657499327,0.114492090903681,-0.206074046480866,-0.0473738720909033,0.0126333493014412,-0.0818871709410604,-0.137966545165005,-0.0673176757949885,-0.110950315518037,0.167801067182278,0.0660111320094841,-0.151026717632182,-0.180519986873022,0.153482704104291,0.0157842576740884,0.119958898448019,0.557044054065571,0.213503416523635
    ## "Rich_M9_No_ATC_1.rep1",0.0244016805532689,-0.0603758568284231,-0.409129555134506,-0.313471187232917,0.0211865937936163,-0.15786107534586,-0.389754975379843,-0.0953873847412387,0.543498441482191,0.648128016804381,-0.151640125038695,-0.0155089805731577,-0.0247866986871272,0.288756351803907,-0.00789903117812725,-0.0100756346712671,-0.0575014624421977,0.0748582409397373,-0.144442045699664,-0.569665907642934
    ## "Rich_M9_No_ATC_1.rep2",-0.218402907477569,-0.187465783602488,-0.00728482168528045,0.234236604948883,0.198309194815452,-0.258356543006706,0.0653746786969956,-0.161557638843878,-0.0347254066238234,-0.0266922052561276,-0.00454772685014345,0.115008577101722,0.171950290949507,0.268820379031386,0.405241544147649,-0.275425386428431,-0.126028231321822,-0.0203824012655673,-0.0790172479675186,0.119444603810996
    ## "Rich_M9_No_ATC_2.rep1",-0.135486852392272,-0.0570433808125654,-0.389239396627237,-0.112903054305338,-0.0812210107330977,0.450903893484182,-0.122966279048948,-0.0176594866793138,0.0508930833682732,-0.00503944413523487,-0.100795236488473,0.111415222811722,0.381247528228862,0.0520187758055778,-0.107853729032672,0.197529258722069,-0.0435643904315998,0.0405869236278252,-0.107943890484224,-0.0666904836498079
    ## "Rich_M9_No_ATC_2.rep2",-0.0126648482611545,-0.073860141821471,-0.348835963138679,-0.0305434625243957,-0.0652365884041723,0.515062722025487,-0.112076727266408,-0.0377655069150882,-0.00770285779519434,-0.0418889219594246,-0.125529017068107,0.138015878928334,0.307856978277444,0.0913422134253299,-0.0765993083080344,0.227472391832429,-0.0675243826296002,0.0143819306928144,-0.133398951133424,0.0537969889382366

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
