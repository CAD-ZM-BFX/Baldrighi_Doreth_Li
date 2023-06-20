#!/usr/local/bin/Rscript
# CAD_xl315_0001, Proteomic data analysis 
# BUT,NUT,PUT,BT,NT,PT each has 3 replicates
# PTvsPUT revised analysis, 
# 1) only PTvsPUT samples perform normalisation
# 2) add BT,BUT samples to PT,PUT comparison
#
# Analysis Performed by Xiaohui Zhao
# School of Medicine, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
# 23/12/2021
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

message("+-------------------------------------------------------------------------------+")
message("+------------ General set-up,libraries calling, metafile generating           --+")
message("+-------------------------------------------------------------------------------+")

cran.pkgs  <- c("DT", "tidyverse", "dplyr", "gridExtra","corrplot", "data.table", "ggplot2",
                "ggfortify","lazyeval", "lubridate", "pheatmap", "reshape2","readr", "rlang",
                "tibble", "tidyr", "wesanderson","WGCNA")
bioc.pkgs  <- c("ncdf4", "MSnbase", "pRoloc", "qPLEXanalyzer", "svach", "mzR",
                "GO.db", "impute", "preprocixOmics")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(bioc.pkgs)
install.packages(cran.pkgs)
devtools::install_github("https://github.com/EvaYiwenWang/PLSDAbatch")

suppressMessages({
  library(MSnbase)
  library(DT)
  library(pRoloc)
  library(tidyverse)
  library(qPLEXanalyzer)
  library(sva)
  library(dplyr)
  library(gridExtra)
  library(readr)
  library(knitr) 
  library(dplyr) 
  library(tibble) 
  library(ggplot2)
  library(xlsx)
  library(ggrepel)
  library(reshape2)
  library(circlize)
  library(cowplot)
})
data(mouse_anno) 


message("+------------        metafile generating           -------------+")

dat.dir     <- "/Users/xz289/Documents/Xuan_Li/xl305_0001"
setwd(dat.dir)

metadata    <- read.csv("P756_metadata.csv", header=T)
Protein.dat <- read.csv("P756_Mergerd_proteins_summary_rmLow_rmMulti_allsamples_28_Dec_2021.csv", header=T)

message("+---- Remove the missings and keep complete                                 ------------------+")
message("+---- The missing patterns did not have untreated/treated completely missing------------------+")

missmat         <- apply(Protein.dat[,c(23:26,31:37,40)], 2, function(x) ifelse(is.na(x),1,0))
missind         <- which(rowSums(missmat)!=0)
Protein.dat.new <- Protein.dat[-missind,-c(27:30,38:39)]  ## 3578 after remove missing

message("+------------        selmarkers gene symbol and protein accessions            -------------+")
message("+------------        Check in mouse_anno availability - reviewed/non-reviewed -------------+")

selgenes    <- c("Nlrp3", "Plk1", "Ddx3x", "Btk", "Pcm1", "Gbp5", "Cdk1", "Pcnt", "Cep192", "Cenpe", 
                 "Hdac1", "Hdac2", "Hdac5", "Hdac6")
selproteinsL <- list(c("Q8R4B8", "Q5SXI6"), c("A0A0U1RNM9", "A0A0U1RQ54", "Q07832"), c("Q62167"),
                     c("P35991", "A2BDW0"), c("Q9R0L6"), c("Q8CFB4", "Q8BMN7"), c("P11440", "Q14AX6", "D3Z2T9"),
                     c("A0A1W2P737", "P48725"), c("E9Q4Y4", "A0A286YDK4"), c("Q6RT24", "A0A0G2JG58", "E9QKK1"),
                     c("O09106"), c("P70288", "A0A0R4J008","Q8BQ10"), c("B7ZDF4", "Q9Z2V6"), 
                     c("Q9Z2V5", "Q3U4Q5", "B1AUA9", "A0A1B0GX25", "Q3UG37", "Q8CGC3"))
selProdat   <- mouse_anno[mouse_anno$GeneSymbol%in%selgenes,]
selproteins <- mouse_anno[mouse_anno$GeneSymbol%in%selgenes,'Accession']

message("+--------------------------------------------------------------------+")
message("+-----DES analysis use the protein list 7216, 09/11/2021         ----+")
message("+-----Use qPLEXanalyzer package                                 -----+")
message("+--------------------------------------------------------------------+")

metadata        <- metadata[-grep(pattern = 'NT|NUT', x =metadata$SampleGroup),]
Protein.uniM    <- Protein.dat.new[Protein.dat.new$Accession%in%mouse_anno$Accessions,]  ## 1249

message("+---- merge isoforms for abundance, use the mean values                      ------------------+")

Protein.dat.new$NewAccession <- gsub("-[1-9]", "", Protein.dat.new$Accession)
NewAccession.Names <- unique(Protein.dat.new$NewAccession[which(duplicated(Protein.dat.new$NewAccession)==T)])

Protein.dat1    <- Protein.dat.new[which(Protein.dat.new$NewAccession%in%NewAccession.Names==F),] ## 2741
Protein.dat2    <- NULL
for(i in 1:length(NewAccession.Names)){
  print(i)
  subT  <- Protein.dat.new[Protein.dat.new$NewAccession==NewAccession.Names[i],]
  merC  <- colMeans(subT[,23:34])
  nsub  <- subT[which(subT$Accession==NewAccession.Names[i]),]
  nsub[23:34] <- merC
  Protein.dat2 <- rbind(Protein.dat2, nsub)
  Protein.dat2
}

Protein.datM   <- rbind(Protein.dat1, Protein.dat2)  ## 3063

message("+---- Generate MSnSet object of the protein data                     ------------------+")
message("+---- After the above filtering, only High FDR confidence keep in the analysis   ------+")

MSnSet_data       <- convertToMSnset(Protein.datM,
                                     metadata=metadata,
                                     indExpData=c(23:34),
                                     Sequences=NULL,
                                     type="protein",
                                     Accessions=5, 
                                     rmMissing=T) ## option rmMissing default is true, 3053 features.
table(fData(MSnSet_data)[,23]) ## High 3063

MSnSet_norm       <-  normalizeScaling(MSnSet_data)

pData(MSnSet_norm)$Plex <- pData(MSnSet_norm)$RepName
MSnSet_norm_batch <- MSnSet_norm
batch      <- pData(MSnSet_norm_batch)$Plex
library(Biobase)
y2         <- ComBat(exprs(MSnSet_norm), batch)
exprs(MSnSet_norm_batch) <- y2 + abs(min(y2))+1

p1 <- pcaPlot(MSnSet_norm, labelColumn="Plex", pointsize=3, x.nudge=1,transform=TRUE,labelsize=3)
p2 <- pcaPlot(MSnSet_norm, labelColumn="SampleGroup", pointsize=3, x.nudge=1,transform=TRUE,labelsize=3)
grid.arrange(p1, p2, ncol = 2, nrow = 1)

p1 <- pcaPlot(MSnSet_norm_batch, labelColumn="Plex", pointsize=3, x.nudge=1,transform=TRUE,labelsize=3)
p2 <- pcaPlot(MSnSet_norm_batch, labelColumn="SampleGroup", pointsize=3, x.nudge=1,transform=TRUE,labelsize=3)
grid.arrange(p1, p2, ncol = 2, nrow = 1)

corrPlot(MSnSet_norm_batch,addValues=FALSE)
hierarchicalPlot(MSnSet_norm_batch)


contrasts <- c(
  BT_vs_BUT = "BT - BUT",
  PT_vs_PUT = "PT - PUT",
  PT_vs_BT = "PT - BT",
  PUT_vs_BUT = "PUT - BUT"
)

diffstats <- computeDiffStats(MSnSet_norm_batch, contrasts=contrasts)
diffexp <- list()
for (i in 1:length(contrasts))
  diffexp[[i]] <- getContrastResults(diffstats=diffstats, contrast=contrasts[i], writeFile= TRUE)
for (i in 1:length(contrasts))
{
  toplot <- diffexp[[i]]$Accessions[diffexp[[i]]$adj.P.Val < 0.05]
  if(length(toplot) > 10)
    toplot <- toplot[1:10]
  print(maVolPlot(diffstats, contrast = contrasts[i], plotType="MA", title= contrasts[i], selectedGenes =
                    toplot, fdrCutOff=0.05))
}

i = 2
PTvsPUT_ori <- diffexp[[2]]
PTvsPUT_ori <- PTvsPUT_ori[, c("Accessions", "GeneSymbol", "log2FC", "adj.P.Val", "t", "B",
                               "PT_1", "PT_2", "PT_3", "PUT_1", "PUT_2", "PUT_3", "Description")]
PTvsPUT_ori <- PTvsPUT_ori[order(PTvsPUT_ori[,4], -PTvsPUT_ori[,3] ),]

write.csv(PTvsPUT_ori, file = "Protein_190920_N7216_F3063_normOri_PTvsPUT_withBP_23_Dec_21.csv")

message("+----Full reviewed+nonreviewed list PTvsPUT add the GO information from ensEMBL2id.mer --+")

load("Ensembl_gene_id_symbol_entrezid_goid_name1006_GRCm39.RData")
colnames(ensEMBL2id.mer)[2] <- "GeneSymbol"
PTvsPUT_ori.merE <- merge(PTvsPUT_ori, ensEMBL2id.mer, by = "GeneSymbol")
PTvsPUT_ori.merE.selGO <- PTvsPUT_ori.merE[grep("microtubule", PTvsPUT_ori.merE$name_1006), ]
Qgo_anno <- read.csv("oriGO_mictubule.csv", header=T)
PTvsPUT_ori.microtubu  <- merge(PTvsPUT_ori.merE.selGO, Qgo_anno, by = "go_id")[,-c(16,19:20)]
PTvsPUT_ori.microtubu  <- unique(PTvsPUT_ori.microtubu)
## 526 itmes. 241 unique proteins, 55 GO ID involved with microtubule,108 genes
write.csv(PTvsPUT_ori.microtubu, file ="PTvsPUT_N3063_all_withBP_GOmicrotubu_Relating_Gene_Summary_23_Dec_2021.csv", row.names=F )

message("+----Volcano plot and highlight the sigDEs microtubule markers in the plots --+")

options(ggrepel.max.overlaps = Inf)
functionPlotDEVolcano_v1<- function(results, sig_cut, logfc_cut, title, selmarkers, xlabel, ylabel) {
  Sel.sub  <- results[results$GeneSymbol%in%selmarkers,]
  volc.plt <- ggplot(data=results, aes(x=log2FC, y=-log10(adj.P.Val), label=GeneSymbol)) +
    geom_vline(xintercept = logfc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
    geom_vline(xintercept = -(logfc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
    geom_hline(yintercept = -log10(sig_cut), colour="black", linetype = "dashed", alpha=0.5) +
    geom_point(data=subset(results, abs(log2FC) < logfc_cut | adj.P.Val > sig_cut), alpha=0.75, size=0.3, colour="grey") +
    geom_point(data=subset(results, adj.P.Val<=sig_cut & log2FC >= logfc_cut),      alpha=0.75, size=0.8, colour="red") +
    geom_point(data=subset(results, adj.P.Val<=sig_cut & log2FC <= -(logfc_cut)),   alpha=0.75, size=0.8, colour="blue") +
    geom_text_repel( data= subset(Sel.sub, adj.P.Val<=sig_cut & log2FC >= logfc_cut),
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    geom_text_repel( data= subset(Sel.sub, adj.P.Val<=sig_cut & log2FC <= -(logfc_cut)),
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    xlab(xlabel) + ylab(ylabel) +
    theme(aspect.ratio=1) +
    ggtitle(title) +
    theme_update(plot.title = element_text(size=16, face="bold", hjust=0.5),
                 axis.title.x = element_text(size=12, face= "bold"),
                 axis.text.x = element_text(size=12, face="bold"),
                 axis.title.y.left = element_text(size=12, face= "bold"),
                 axis.text.y = element_text(size=12, face="bold")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white", colour = NA)) 
  volc.plt
}
sig_cut       <- 0.05
logfc_cut     <- 0.4
selmarkersVol <- unique(c(PTvsPUT_ori.merE.selGO$GeneSymbol, selgenes))
ylabel        <- bquote("-log"[10]~"(adj.p.value)")
xlabel        <- "log2FC (Activated/Primed_BASU-PLK1)"
title         <- "PT vs PUT ori"
results.ori   <- PTvsPUT_ori[, c("GeneSymbol", "log2FC", "adj.P.Val")]
results.ori   <- results.ori[-which(duplicated(results.ori$GeneSymbol)==T),]
PTvsPUT.oriVolt <- functionPlotDEVolcano_v1(results.ori, sig_cut, logfc_cut, title, selmarkersVol, xlabel, ylabel) 

pdf("Protein_190920_N7216_F3063_normOri_PTvsPUT_withBP_VolcanoPlot_p05_l2fc04_23_Dec_21.pdf")
print(PTvsPUT.oriVolt)
dev.off()


message("+-------------------------------------------------------------------+")
message("+---------Use all reviewed proteins to perform the DEG analysis-----+")
message("+    and use normalise instead of group normalise to check          +")
message("+---------   Revised on 23/12/2021                         ---------+")

uniprotList            <- read.xlsx("uniprot-filtered-organism-Mus_03_2021.xlsx", sheetIndex=1)
colnames(uniprotList)  <- c("Accession", "UniProt")
Protein.datR  <- Protein.datM[Protein.datM$NewAccession%in%uniprotList$Accession,] ## 1252
## unique gene 1249, plus 2 missing and 1 Naca with (NACA, NACAM proteins)
Protein.datR.rmNA <- Protein.datR[complete.cases(Protein.datR[,23:34]), ] 

MSnSet_dataR       <- convertToMSnset(Protein.datR.rmNA,
                                      metadata=metadata,
                                      indExpData=c(23:34),
                                      Sequences=NULL,
                                      type="protein",
                                      Accessions=5, 
                                      rmMissing=T) ## option rmMissing default is true, 1252
table(fData(MSnSet_dataR)[,3]) ; 
## High 1239, Medium 13
table(fData(MSnSet_dataR)[,23])
## Group 1, 1252

message("+---- normalised all samples by median ------------------+")

MSnSet_normR       <-  normalizeScaling(MSnSet_dataR)

pData(MSnSet_normR)$Plex <- pData(MSnSet_normR)$RepName
MSnSet_normR_batch <- MSnSet_normR
batch      <- pData(MSnSet_normR_batch)$Plex
y2         <- ComBat(exprs(MSnSet_normR), batch)
exprs(MSnSet_normR_batch) <- y2 + abs(round(min(y2))) + 1

pdf("PTvsPUT_reviewed_withBP_N1252_rmBatch_PCA_23_Dec_2021.pdf")
p1 <- pcaPlot(MSnSet_normR, labelColumn="Plex", pointsize=3, x.nudge=1,transform=TRUE,labelsize=3)
p2 <- pcaPlot(MSnSet_normR, labelColumn="SampleGroup", pointsize=3, x.nudge=1,transform=TRUE,labelsize=3)
grid.arrange(p1, p2, ncol = 2, nrow = 1)

p1 <- pcaPlot(MSnSet_normR_batch, labelColumn="Plex", pointsize=3, x.nudge=1,transform=TRUE,labelsize=3)
p2 <- pcaPlot(MSnSet_normR_batch, labelColumn="SampleGroup", pointsize=3, x.nudge=1,transform=TRUE,labelsize=3)
grid.arrange(p1, p2, ncol = 2, nrow = 1)
dev.off()

message("+---- Produce the box plot for PT PUT each plex----------+")
## non-normalised boxplot
MSnSetObj <- MSnSet_dataR[,grep("PT|PUT", colnames(exprs(MSnSet_dataR)))]
plotDat <- exprs(MSnSetObj) %>% 
           as.data.frame() %>% pivot_longer(names_to = "SampleName", 
           values_to = "Intensity", everything()) %>% 
           mutate(logInt = log2(Intensity)) %>% 
           filter(is.finite(logInt)) %>% left_join(pData(MSnSetObj), 
           "SampleName") 
plotDat$SampleGroup <- ifelse(plotDat$SampleGroup=="PUT", "Primed", "Activated")
plotDat$SampleName  <- paste0(plotDat$SampleGroup, "_", plotDat$BioRep)
plotDat$SampleName  <- factor(plotDat$SampleName,
                       levels = c("Primed_1", "Primed_2", "Primed_3",
                                  "Activated_1","Activated_2","Activated_3"),
                       ordered = TRUE)

NonnormBox <- ggplot(plotDat, aes(x =SampleName, y = logInt, fill=SampleGroup)) + 
  geom_boxplot()  + 
  scale_fill_manual(values = c("red2", "tan2")) +
  labs(x = "", y = "log2(Intensity)", title = title) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, 
   hjust = 1), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = FALSE) +
  theme(text = element_text(size = 14), 
        plot.title = element_text(size = 14, hjust = 0.5), aspect.ratio = 1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
## normalised
MSnSetObj <- MSnSet_normR[,grep("PT|PUT", colnames(exprs(MSnSet_normR)))]
plotDat <- exprs(MSnSetObj) %>% 
  as.data.frame() %>% pivot_longer(names_to = "SampleName", 
                                   values_to = "Intensity", everything()) %>% 
  mutate(logInt = log2(Intensity)) %>% 
  filter(is.finite(logInt)) %>% left_join(pData(MSnSetObj), 
                                          "SampleName") 
plotDat$SampleGroup <- ifelse(plotDat$SampleGroup=="PUT", "Primed", "Activated")
plotDat$SampleName  <- paste0(plotDat$SampleGroup, "_", plotDat$BioRep)
plotDat$SampleName  <- factor(plotDat$SampleName,
                              levels = c("Primed_1", "Primed_2", "Primed_3",
                                         "Activated_1","Activated_2","Activated_3"),
                              ordered = TRUE)
normBox <- ggplot(plotDat, aes(x =SampleName, y = logInt, fill=SampleGroup)) + 
  geom_boxplot()  + 
  scale_fill_manual(values = c("red2", "tan2")) +
  labs(x = "", y = "log2(Intensity)", title = title) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, 
                                                hjust = 1), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = FALSE) +
  theme(text = element_text(size = 14), 
        plot.title = element_text(size = 14, hjust = 0.5), aspect.ratio = 1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
## ## normalised+batch
MSnSetObj <- MSnSet_normR_batch[,grep("PT|PUT", colnames(exprs(MSnSet_normR_batch)))]
plotDat <- exprs(MSnSetObj) %>% 
  as.data.frame() %>% pivot_longer(names_to = "SampleName", 
                                   values_to = "Intensity", everything()) %>% 
  mutate(logInt = log2(Intensity)) %>% 
  filter(is.finite(logInt)) %>% left_join(pData(MSnSetObj), 
                                          "SampleName") 
plotDat$SampleGroup <- ifelse(plotDat$SampleGroup=="PUT", "Primed", "Activated")
plotDat$SampleName  <- paste0(plotDat$SampleGroup, "_", plotDat$BioRep)
plotDat$SampleName  <- factor(plotDat$SampleName,
                              levels = c("Primed_1", "Primed_2", "Primed_3",
                                         "Activated_1","Activated_2","Activated_3"),
                              ordered = TRUE)
normBatchBox <- ggplot(plotDat, aes(x =SampleName, y = logInt, fill=SampleGroup)) + 
  geom_boxplot()  + 
  scale_fill_manual(values = c("red2", "tan2")) +
  labs(x = "", y = "log2(Intensity)", title = title) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45, 
                                                hjust = 1), plot.title = element_text(hjust = 0.5)) + 
  guides(fill = FALSE) +
  theme(text = element_text(size = 14), 
        plot.title = element_text(size = 14, hjust = 0.5), aspect.ratio = 1) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA))
###------------------------SFig5c------------------------------------#####
pdf("PTvsPUT-NonNorm_Norm_NormBatch_Logintens_boxplot_26_jan_2022.pdf")
NonnormBox
normBox
normBatchBox
dev.off()
###------------------------SFig5c------------------------------------#####



message("+customised a good pca plot for paper. revised 26/01/2022 by changing the groups names,to activated and primed.+")

customised_PCA <- function (MSnSetObj,transform = TRUE, colourBy = "SampleGroup", title = "", labelColumn = "BioRep", 
          labelsize = 4, pointsize = 3, x.nudge = 4, x.PC = 1) 
{
 intensities <- exprs(MSnSetObj) %>% as.data.frame() %>% na.omit()
  if (!nrow(intensities)) {
    return(NULL)
  }
  if (transform) {
    intensities <- log2(intensities+1)
  }
  pca <- intensities %>% t() %>% prcomp()
  pcaVariance <- round((pca$sdev^2/sum(pca$sdev^2)) * 100)
  
  plotDat <- as.data.frame(pca$x) %>% rownames_to_column("SampleName") %>% 
    left_join(pData(MSnSetObj), "SampleName") 
  
  xLab <- str_c("PC", x.PC, ", ", pcaVariance[x.PC], "% variance")
  yLab <- str_c("PC", x.PC + 1, ", ", pcaVariance[x.PC + 1], 
                "% variance")
  xPC <- str_c("PC", x.PC) %>% sym()
  yPC <- str_c("PC", x.PC + 1) %>% sym()
  plotDat$SampleGroup <- factor(plotDat$SampleGroup, 
                                   levels= c("Activated_BASU", "Primed_BASU",
                                             "Activated_BASU-PLK1", "Primed_BASU-PLK1"))
  pcPlot <- ggplot(plotDat, aes(x = PC1, y = PC2, colour=SampleGroup,  shape=RepName)) +
    geom_point(size = pointsize) + 
    scale_color_manual(name="Group", values = c("blue3", "turquoise2", "red2", "tan2")) + 
    labs(x = xLab, y = yLab, fill = NULL, title = title) + 
    scale_shape_manual(name="Batch",values=c(0,1,2)) + 
    geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.2) +
    geom_vline(xintercept=0, linetype="dashed", color = "black", size=0.2) +
    theme(axis.title.x = element_text(face= "bold",size=16),
          axis.text.x = element_text(face= "bold",size=20),
          axis.title.y.left = element_text(face= "bold",size=16),
          axis.text.y = element_text(face= "bold",size=20),
          legend.title = element_text(face= "bold",size=16),
          legend.text = element_text(face= "bold",size=16)) +
    #geom_text_repel(aes(label=SampleName), show.legend = FALSE, size=2) +
    theme_bw() + theme(text = element_text(face= "bold",size=16), plot.title = element_text(size = 16, hjust = 0.5), aspect.ratio = 1) +
      theme( panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  pcPlot
}


#######----------------------------SFig5b-----------------------------------------------------------#########
pdf("PTvsPUT_reviewed_withBP_N1252_rmBatch_customised_PCA_23_Dec_2021.pdf")
pData(MSnSet_normR_batch)$SampleGroup <- c("Primed_BASU", "Activated_BASU",
                                           "Primed_BASU-PLK1", "Activated_BASU-PLK1",
                                           "Activated_BASU-PLK1", "Primed_BASU-PLK1",
                                           "Activated_BASU", "Primed_BASU", 
                                           "Activated_BASU", "Primed_BASU-PLK1",
                                           "Activated_BASU-PLK1", "Primed_BASU")

PCA_cust <- customised_PCA(MSnSet_normR_batch,transform = TRUE, colourBy = "SampleGroup", title = "", labelColumn = "BioRep", 
                            labelsize = 4, pointsize = 4, x.nudge = 4, x.PC = 1) 
print(PCA_cust)
dev.off()
#######--------------------------------------------------------------------------------------------#########


diffstats <- computeDiffStats(MSnSet_normR_batch, contrasts=contrasts)
diffexp <- list()
for (i in 1:length(contrasts))
  diffexp[[i]] <- getContrastResults(diffstats=diffstats, contrast=contrasts[i], writeFile= TRUE)
for (i in 1:length(contrasts))
{
  toplot <- diffexp[[i]]$Accessions[diffexp[[i]]$adj.P.Val < 0.05]
  if(length(toplot) > 10)
    toplot <- toplot[1:10]
  print(maVolPlot(diffstats, contrast = contrasts[i], plotType="MA", title= contrasts[i], selectedGenes =
                    toplot, fdrCutOff=0.05))
}

i = 1
BTvsBUT_oriR <- diffexp[[i]]
BTvsBUT_oriR <- BTvsBUT_oriR[, c("Accessions", "GeneSymbol", "log2FC", "adj.P.Val", "t", "B",
                                 "BT_1", "BT_2", "BT_3", "BUT_1", "BUT_2", "BUT_3", "Description")]
i = 2
PTvsPUT_oriR <- diffexp[[i]]
PTvsPUT_oriR <- PTvsPUT_oriR[, c("Accessions", "GeneSymbol", "log2FC", "adj.P.Val", "t", "B",
                                 "PT_1", "PT_2", "PT_3", "PUT_1", "PUT_2", "PUT_3", "Description")]
i = 3
PTvsBT_oriR <- diffexp[[i]]
PTvsBT_oriR <- PTvsBT_oriR[, c("Accessions", "GeneSymbol", "log2FC", "adj.P.Val", "t", "B",
                               "PT_1", "PT_2", "PT_3", "BT_1", "BT_2", "BT_3", "Description")]
i = 4
PUTvsBUT_oriR <- diffexp[[i]]
PUTvsBUT_oriR <- PUTvsBUT_oriR[, c("Accessions", "GeneSymbol", "log2FC", "adj.P.Val", "t", "B",
                                   "PUT_1", "PUT_2", "PUT_3", "BUT_1", "BUT_2", "BUT_3", "Description")]


write.csv(PTvsPUT_oriR, file = "Protein_190920_N7216_F1252_reviewed_withBP_PTvsPUT_23_Dec_21.csv", row.names=F)
write.csv(PTvsBT_oriR, file = "Protein_190920_N7216_F1252_reviewed_withBP_PTvsBT_23_Dec_21.csv", row.names=F)
write.csv(PUTvsBUT_oriR, file = "Protein_190920_N7216_F1252_reviewed_withBP_PUTvsBUT_23_Dec_21.csv", row.names=F)
write.csv(BTvsBUT_oriR, file = "Protein_190920_N7216_F1252_reviewed_withBP_BTvsBUT_23_Dec_21.csv", row.names=F)


message("+---- Volcano Plot for reviewed proteins N1245 -------------+")

logfc_cut     <- 0.40
title         <- "PT vs PUT"

functionPlotDEVolcano_v2 <- function(results, sig_cut, logfc_cut, title, selmarkers, xlabel, ylabel) {
  Sel.sub  <- results[results$GeneSymbol%in%selmarkers,]
  Sel.sub$GeneSymbol <- c("DDX3X", "BTK", "BTK", "SYK", "SYK", "PKR", "NLRP3", "GBP5", "GBP5", "NLRP3")
  Sel.sub  <- Sel.sub[-c(3,5,8,10),]
  volc.plt <- ggplot(data=results, aes(x=log2FC, y=-log10(adj.P.Val), label=GeneSymbol)) +
    geom_vline(xintercept = logfc_cut,     colour="black", linetype = "dashed", alpha=0.5, size=1.2) +
    geom_vline(xintercept = -(logfc_cut),  colour="black", linetype = "dashed", alpha=0.5, size=1.2) +
    geom_hline(yintercept = -log10(sig_cut), colour="black", linetype = "dashed", alpha=0.5, size=1.2) +
    geom_point(data=subset(results,  abs(log2FC) < logfc_cut | adj.P.Val > sig_cut), alpha=0.5, size=0.8, colour="grey") +
    geom_point(data=subset(results, adj.P.Val<=sig_cut & log2FC >= logfc_cut),      alpha=0.3, size=0.8, colour="red") +
    geom_point(data=subset(results, adj.P.Val<=sig_cut & log2FC <= -(logfc_cut)),   alpha=0.5, size=0.8, colour="blue") +
    geom_point(data=Sel.sub,  size=3, colour="deeppink4", shape=17) +
    geom_text_repel(data= Sel.sub,show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=6) +
    ##geom_text_repel( data= subset(Sel.sub, adj.P.Val<=sig_cut & log2FC <= -(logfc_cut)),
    ##                 show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    xlab(xlabel) + ylab(ylabel) +
    theme(aspect.ratio=1) +
    theme(axis.title.x = element_text(face= "bold",size=18),
          axis.title.y = element_text(face= "bold",size=18),
          axis.text.x = element_text(face= "bold",size=18),
          axis.text.y = element_text(face= "bold",size=18)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white", colour = NA)) 
  volc.plt
}

selnewMarkers <- c("Ddx3x", "Btk", "Syk", "Eif2ak2", "Mark2", "Gbp5", "Nlrp3", "Plk1",
                   "Alb", "Pc", "Pcca", "Pccb", "Banf1", "Acaca", "Hist1h1b", "Mccc1",
                   "Krt16")
PTvsPUT_oriRVolt <- functionPlotDEVolcano_v2(PTvsPUT_oriR, sig_cut, logfc_cut=0.4, title, selnewMarkers, xlabel, ylabel) 

## Check the PTvsBT significant genes and highlight with dark green color
PTvsBT_oriR.sig   <- subset(PTvsBT_oriR, adj.P.Val < 0.05 & abs(log2FC)  >=0.4) ## 
PUTvsBUT_oriR.sig <- subset(PUTvsBUT_oriR, adj.P.Val < 0.05 & abs(log2FC) >=0.4)
PTvsPUT_oriR.sig  <- subset(PTvsPUT_oriR, adj.P.Val < 0.05 & abs(log2FC) >=0.4) ## 

sigmer     <- merge(PTvsBT_oriR.sig[,c("GeneSymbol", "log2FC", "adj.P.Val")],
                  PTvsPUT_oriR.sig[,c("GeneSymbol", "log2FC", "adj.P.Val")], by = "GeneSymbol")
sigmer.sub <- subset(sigmer, sign(log2FC.x)!=sign(log2FC.y))
enrichBT     <- sigmer.sub[sigmer.sub$log2FC.x<0&sigmer.sub$log2FC.y>0, "GeneSymbol"]
newdata_oriR <- PTvsPUT_oriR.sig[PTvsPUT_oriR.sig$GeneSymbol%in%enrichBT,]
## there is no PT enrich ones
enrichPT     <- sigmer.sub[sigmer.sub$log2FC.x>0&sigmer.sub$log2FC.y>0, "GeneSymbol"] ## 0
newdata_oriRP <- PTvsPUT_oriR.sig[PTvsPUT_oriR.sig$GeneSymbol%in%enrichPT,]


## PUT enrich compare to BUT and PT
sigmerU     <- merge(PUTvsBUT_oriR.sig[,c("GeneSymbol", "log2FC", "adj.P.Val")],
                    PTvsPUT_oriR.sig[,c("GeneSymbol", "log2FC", "adj.P.Val")], by = "GeneSymbol")
sigmer.subU <- subset(sigmerU, sign(log2FC.x)!=sign(log2FC.y))
enrichBUT     <- sigmer.subU[sigmer.subU$log2FC.x<0&sigmer.subU$log2FC.y<0, "GeneSymbol"]
newdata_oriRU <- PTvsPUT_oriR.sig[PTvsPUT_oriR.sig$GeneSymbol%in%enrichBUT,] ## 0
enrichPUT     <- sigmer.subU[sigmer.subU$log2FC.x>0&sigmer.subU$log2FC.y<0, "GeneSymbol"]
newdata_oriRPU <- PTvsPUT_oriR.sig[PTvsPUT_oriR.sig$GeneSymbol%in%enrichPUT,] 
write.xlsx(newdata_oriR, file="BT_PUT_enrich_volcano_N40_N37_lists_Jan_2022.xlsx", sheetName = "BTenrich", 
           col.names = TRUE, row.names = TRUE, append = FALSE)
write.xlsx(newdata_oriRPU, file="BT_PUT_enrich_volcano_N40_N37_lists_Jan_2022.xlsx", sheetName = "PUTenrich", 
           col.names = TRUE, row.names = TRUE, append = TRUE)

####----------------------------------------------Fig3e---------------------------------------------####
pdf("Protein_190920_N7216_F1252_normOri_withBP_PTvsPUT_VolcanoPlot_p05_l2fc040_FC1.32_24_Jan_21.pdf")
selnewMarkers1 <- c("Ddx3x", "Btk", "Syk", "Eif2ak2", "Gbp5", "Nlrp3")
PTvsPUT_oriRVolt1 <- functionPlotDEVolcano_v2(PTvsPUT_oriR, sig_cut, logfc_cut=0.4, title="NULL", selnewMarkers1, xlabel, ylabel) 
print(PTvsPUT_oriRVolt1)
dev.off()
####-----------------------------------------------------------------------------------------------####

pdf("Protein_190920_N7216_F1252_normOri_withBP_PTvsPUT_VolcanoPlot_p05_l2fc040_FC1.32_23_Dec_21.pdf")

newVolt <- PTvsPUT_oriRVolt + geom_point(data=newdata_oriR, aes(x=log2FC, y=-log10(adj.P.Val)),
                                           alpha=0.75, size=1.2, colour="purple")+
           geom_text(x=-1, y=10.5, label="N=172", color="blue") +
           geom_text(x=1.5, y=10.5, label="N=247", color="red") +
           geom_text(x=2, y=9, label="NB=40", color="purple") +
           geom_point(data=newdata_oriRPU, aes(x=log2FC, y=-log10(adj.P.Val)),
                       alpha=0.75, size=1.2, colour="turquoise2")+
           geom_text(x=-2, y=9, label="NP=37", color="turquoise2") 

           
print(newVolt)

dev.off()

enrichBT_0.4     <- sigmer.sub[sigmer.sub$log2FC.x<0&sigmer.sub$log2FC.y>0, ]
colnames(enrichBT_0.4) <- c("GeneSymbol", "log2FC(PTvsBT)", "adj.P.Val(PTvsBT)",
                            "log2FC(PTvsPUT)", "adj.P.Val(PTvsPUT)")
##
PTvsBT_oriR.sig  <- subset(PTvsBT_oriR, adj.P.Val < 0.05 & abs(log2FC)  >=0.6) 
PTvsPUT_oriR.sig <- subset(PTvsPUT_oriR, adj.P.Val < 0.05 & abs(log2FC) >=0.6) 
sigmer     <- merge(PTvsBT_oriR.sig[,c("GeneSymbol", "log2FC", "adj.P.Val")],
                      PTvsPUT_oriR.sig[,c("GeneSymbol", "log2FC", "adj.P.Val")], by = "GeneSymbol")
sigmer.sub <- subset(sigmer, sign(log2FC.x)!=sign(log2FC.y))
enrichBT_0.6     <- sigmer.sub[sigmerBP.sub$log2FC.x<0&sigmer.sub$log2FC.y>0, ]
colnames(enrichBT_0.6) <- c("GeneSymbol", "log2FC(PTvsBT)", "adj.P.Val(PTvsBT)",
                            "log2FC(PTvsPUT)", "adj.P.Val(PTvsPUT)")

write.csv(enrichBT_0.4, file = "Protein_190920_N7216_F1252_normOri_withBP_PTvsPUT_reviewed_enrichBT_enrichPT_l2fc0.4_23_Dec_2021.csv", row.names=F)
write.csv(enrichBT_0.6, file = "Protein_190920_N7216_F1252_normOri_withBP_PTvsPUT_reviewed_enrichBT_enrichPT_l2fc0.6_23_Dec_2021.csv", row.names=F)


message("+--- Highlight microtubule relating genes volcano plot-------------------------+")

microRelat <- read.csv("PTvsPUT_N1252_withBP_GOmicrotubu_Relating_Gene_Summary_23_Dec_2021.csv", header=T)
selmicroMarkers <- unique(microRelat$GeneSymbol[microRelat$adj.P.Val<=0.05&abs(microRelat$log2FC)>=0.4])
selmicroMarkers <- selmicroMarkers[-16]
pdf("Protein_190920_N7216_F1252_normOri_withBP_PTvsPUT_microtubule_marked_VolcanoPlot_p05_l2fc040_FC1.32_23_Dec_21.pdf")
microVolt <- functionPlotDEVolcano_v2(PTvsPUT_oriR, sig_cut, logfc_cut=0.4, title, selmicroMarkers, xlabel, ylabel) 
print(microVolt)
dev.off()


message("+-------- GO pathway analysis using clusterProfiler---------------------+")

library("clusterProfiler")
library("DOSE")
library("GSEABase")
library("AnnotationHub")
library("org.Mm.eg.db")
library("gage")
library("gageData")
library("enrichplot")
library("ggraph")
library("ggforce")
library("enrichR")

colnames(ensEMBL2id.mer)   <- c("ENSEMBL", "SYMBOL", "ENTREZID", "GO_ID", "GO_Desc", "Desc", "GO_func")

PTvsPUT_oriR.sig   <- subset(PTvsPUT_oriR, adj.P.Val <= 0.05 & abs(log2FC) >= 0.4) ## 419
PTvsPUT_oriR.sig    <- PTvsPUT_oriR.sig[, c("GeneSymbol", "log2FC")]
colnames(PTvsPUT_oriR.sig)  <- c("SYMBOL", "L2FC")

PTvsPUT_oriR.sigM  <- merge(PTvsPUT_oriR.sig, ensEMBL2id.mer, by = "SYMBOL")
PTvsPUT_oriR.sigM  <- unique(PTvsPUT_oriR.sigM[,c(4,3,1,2)])
PTvsPUT_oriR.sigM  <- PTvsPUT_oriR.sigM[order(-PTvsPUT_oriR.sigM$L2FC), ] ## 400

Presdf  <- PTvsPUT_oriR.sigM
PgeneL  <- PTvsPUT_oriR.sigM$L2FC
names(PgeneL) <- PTvsPUT_oriR.sigM$ENTREZID
PgeneLS <- PTvsPUT_oriR.sigM$L2FC
names(PgeneLS) <- PTvsPUT_oriR.sigM$SYMBOL

group.GOBP <- groupGO(gene = as.character(PTvsPUT_oriR.sigM$ENTREZID), OrgDb = org.Mm.eg.db, ont = "BP", level = 4, readable = TRUE)
group.GOCC <- groupGO(gene = as.character(PTvsPUT_oriR.sigM$ENTREZID), OrgDb = org.Mm.eg.db, ont = "CC", level = 4, readable = TRUE)
group.GOMF <- groupGO(gene = as.character(PTvsPUT_oriR.sigM$ENTREZID), OrgDb = org.Mm.eg.db, ont = "MF", level = 4, readable = TRUE)

GO.all  <- enrichGO(gene = Presdf$SYMBOL, OrgDb = org.Mm.eg.db,
                    keyType = "SYMBOL", ont = "ALL",
                    pAdjustMethod = "BH", pvalueCutoff  = 0.05 )

table(GO.all$ONTOLOGY)
## BP  CC  MF 
## 552 120  107 
GO.KEGG.all <-  enrichKEGG(Presdf$ENTREZID, organism = "mmu")  
GO.KEGG.all <- setReadable(GO.KEGG.all, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

PgeneL.new  <- PgeneL[!is.na(names(PgeneL))]
gseGO.all   <- gseGO(geneList = PgeneL.new, OrgDb  = org.Mm.eg.db, ont = "All", nPermSimple = 10000,
      minGSSize = 5, maxGSSize = 350, pvalueCutoff = 0.05, verbose = TRUE)
gseGO.all   <- setReadable(gseGO.all, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
gseKEGG.all <- gseKEGG(geneList = PgeneL.new, organism = "mmu", keyType = "kegg",nPermSimple = 10000,
                       minGSSize =5, maxGSSize = 350, pvalueCutoff = 0.05, pAdjustMethod = "BH")
gseKEGG.all <- setReadable(gseKEGG.all, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

save(group.GOBP,group.GOCC,group.GOMF,GO.all,GO.KEGG.all, Presdf,gseGO.all,gseKEGG.all,
     file ="PTvsPUT_withBASU_GeneOntology_data_23_Dec_21.RData")

write.csv(GO.all, file = "PTvsPUT_withBASU_p05_l2fc04_pathway_BP552_MF120_CC107_summary_23_Dec_2021.csv")
write.csv(GO.KEGG.all, file = "PTvsPUT_withBASU_p05_l2fc04_KEGG55_summary_23_Dec_2021.csv")

message("+---- Generate a barplot for microtubule relating pathways plot with identified sig DEs---+")

pathway.dat <- GO.all[grep("microtubule", GO.all$Description),] ## three BP and 2 CC
pgenes      <- lapply(pathway.dat$geneID, function(x) strsplit(x, split="/")[[1]])
BP_list_up  <- list();
BP_list_tot <- list()
for (i in 1:nrow(pathway.dat)){
  df_tmp  <- Presdf[Presdf$SYMBOL%in%pgenes[[i]],]
  tmp_up  <- length(subset(df_tmp$L2FC, df_tmp$L2FC > 0))
  tmp_tot <- dim(df_tmp)[1]
  BP_list_up[[i]]  <- tmp_up
  BP_list_tot[[i]] <- tmp_tot
}

pathway.dat$UP    <- as.numeric(unlist(BP_list_up))
pathway.dat$DOWN  <- as.numeric(unlist(BP_list_tot)) - pathway.dat$UP
pathway.dat$DOWN  <- -pathway.dat$DOWN
pathway.dat$Count <- as.numeric(unlist(BP_list_tot))
pathway.dat_molten <- melt(pathway.dat[,c(3,7,10:12)],
                           id.vars=c("Description","Count", "p.adjust") )
pathway.dat_molten$p.adjust <- as.numeric(pathway.dat_molten$p.adjust)

brks <- seq(-5, 15, 5)
lbls = as.character(c(seq(5, 0, -5), seq(5, 15, 5)))


pdf("Protein_190920_N7216_F1252_reviewed_withBP_GO_BP3_CC2_microtubule_pathways_barplot_23_Dec_2021.pdf", width=8, height=5)
barplotSum <- ggplot(pathway.dat_molten, aes(x = reorder(Description, -p.adjust), y = value, fill = variable)) +   
  geom_bar(stat = "identity", width = 0.8) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls,
                     limits=c(-20,20)) + # Labels
  coord_flip() +  
  xlab("") +
  ylab("") +
  guides(fill=FALSE) +
  theme(plot.title = element_text(hjust = .2),
        axis.ticks = element_blank()
  ) +   
  scale_fill_manual(values = c("red", "blue"))+  # Color palette
  theme_bw() +
  theme_update(axis.title.x = element_text(size=12, face= "bold"),
               axis.text.x = element_text(size=12, face="bold"),
               axis.title.y = element_text(size=12, face= "bold"),
               axis.text.y.left = element_text(size=12, face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white", colour = NA)) 

barplotSum

dev.off()
message("+---- Revised by Xuan 10/01/2022 and try to combine BP pathways within groups selected ---+")
## 1) similarity first semantic after removing some pathways from go.all
## 2) check the overlap with Xuan's group defination
rmPaths <- paste0("GO:", c("0019079", "0016032", "0045069", "0019058", "0009615", "0051607", "1903900",
                         "0050792", "0007051", "0007059", "0000910", "0075522", "0050688", "0045070",
                         "0039694", "0055001", "0060236", "1901989", "0098586", "0019080", "1902850",
                         "0030098", "0048525", "0061640", "0050691", "0010971", "0039528", "0051146",
                         "0035088", "0061245", "0002706", "0002285", "0002708", "0030217", "0002726",
                         "0051546", "0000082", "0002369", "0086003", "0060048", "0021915", "0055013",
                         "0045197", "0050870", "0099092", "0099091", "0098871", "0005876", "0042383"))
GO.alldat    <- as.data.frame(GO.all)
newPaths.all <- GO.alldat[-which(GO.alldat$ID%in%rmPaths==T),]
table(newPaths.all$ONTOLOGY)
## BP  CC  MF 
## 508 115 107 
library("pathview")
library("rrvgo")
library("GOSemSim")

GOSeSimdata.BP   <- godata(OrgDb = "org.Mm.eg.db", keytype = "ENTREZID", ont="BP", computeIC = TRUE)
GOSeSimdata.CC   <- godata(OrgDb = "org.Mm.eg.db", keytype = "ENTREZID", ont="CC", computeIC = TRUE)
GOSeSimdata.MF   <- godata(OrgDb = "org.Mm.eg.db", keytype = "ENTREZID", ont="MF", computeIC = TRUE)
GOSemSim_Data_function <- function(goterm, gofun, model, goName, threshcut, Project, max.overlap, GOSeSimdata){
  
  ego.dat       <- goterm
  ego.dat       <- subset(ego.dat, ONTOLOGY==gofun) 
  ego.dat.dim   <- dim(ego.dat)[1]
  
  simMatrix     <- calculateSimMatrix(ego.dat$ID, orgdb="org.Mm.eg.db", ont=gofun, semdata=GOSeSimdata, method="Rel")
  scores        <- setNames(-log10(ego.dat$qvalue), ego.dat$ID)
  reducedTerms  <- reduceSimMatrix(simMatrix, scores, threshold=threshcut, orgdb="org.Mm.eg.db")
  maxclust      <- max(reducedTerms$cluster)
  message("+---------               Supplementary data                                   --------------+")
  
  write.csv(reducedTerms, file = paste0(Project, "-", model, "_GOSeSim_", goName, "_", gofun, "_N", ego.dat.dim, "_thresh_",
                                       threshcut, "_cluster", maxclust, "_summaryTable.csv"), row.names=F)
  save(reducedTerms, simMatrix, file = paste0(Project,"-", model,"_GOSeSim_", goName, "_", gofun, ".RData"))
  
  options(ggrepel.max.overlaps = max.overlap)
  
  message("+---------               Supplementary Fig                                   --------------+")
  
  pdf(paste0(Project,"-", model, "_GOSeSim_", goName, "_", gofun, "_", ego.dat.dim, "_thresh_", threshcut, "_cluster", maxclust, "_scatterPlot.pdf"),
      width=8, height=6)
  scatterPlot(simMatrix, reducedTerms, size="score", labelSize = 3)
  dev.off()
}
model <- "PTvsPUT"
allBP.sem  <- GOSemSim_Data_function(newPaths.all,"BP", model, "allBP", 0.9, Project, 508, GOSeSimdata.BP)
load("CAD_xl315_0001-PTvsPUT_GOSeSim_allBP_BP.RData")
G1.BP <- reducedTerms[grep("RNA", reducedTerms$term),]  ## 45 items C1-33, C3-4, C5-1, C6-1,C14-6
## (4 items-mRNA metabolic C1-regulation of RNA)
G2.BP <- reducedTerms[grep("metabolic process", reducedTerms$term),] ## 59 items-mRNA metabolic C1
G3.BP <- reducedTerms[grep("cell shape|cell morphogenesis|cell polarity", reducedTerms$term),] ## 4 terms--cell shape&polarity C4-2, C10-2
G4.BP <- reducedTerms[grep("circadian", reducedTerms$term),]  ## 3 items--circadian rhythm C8
G5.BP <- reducedTerms[grep("phosphorylation", reducedTerms$term),] ## 7items-mRNA metabolic C1
G6.BP <- reducedTerms[grep("actin|intermediate filament|microtubule", reducedTerms$term),] ## 53 items. C6-10, C7-14, C9-3, C13-2, C18-2
G7.BP <- reducedTerms[grep("microtubule|dynactin", reducedTerms$term),] ## 2 microtubule C18
G8.BP <- reducedTerms[grep("protein localization", reducedTerms$term),] ## 27 --protein localization to nucleus: C3
G9.BP <- reducedTerms[grep("immune|immunity|inflammasome|inflammatory|NF-kappaB|interferon|Toll-like receptor|cytokine|interleukin", reducedTerms$term),] 
## 54, C1-22,C5-15,C16-14 are the main clusters. C6-2, C12-1
G10.BP <- reducedTerms[grep("inflammasome|pattern recognition receptor", reducedTerms$term),] ## C5-3&C6-2 5 items
G11.BP <- reducedTerms[grep("transport", reducedTerms$term),] ## 23 C3, protein localization to nucleus
G12.BP <- reducedTerms[grep("dissambly|degradation", reducedTerms$term),] ## no terms, possible in CC or MF
BP.sel <- list(G1.BP, G2.BP, G3.BP, G4.BP, G5.BP, G6.BP, G7.BP, G8.BP, G9.BP, G10.BP, G11.BP)
BP.selN <- paste0("G", c(1:12)," ", c("RNA", "Metabolic", "Cell", "Circadian", "Phosphorylation", "Cytoskeleton", 
 "Microtubule", "Localization", "Immune", "Inflammation", "Transport", "Dissamb"))
write.xlsx(BP.sel[[1]], file="Selected_GOBP_Group_10_Jan_2022.xlsx", sheetName = BP.selN[1], 
           col.names = TRUE, row.names = TRUE, append = FALSE)
for(i in 2:11){
  print(i)
  write.xlsx(BP.sel[[i]], file="Selected_GOBP_Group_10_Jan_2022.xlsx", sheetName = BP.selN[i], 
             col.names = TRUE, row.names = TRUE, append = TRUE)
}
allCC.sem  <- GOSemSim_Data_function(newPaths.all,"CC", model, "allCC", 0.9, Project, 115, GOSeSimdata.CC)
allMF.sem  <- GOSemSim_Data_function(newPaths.all,"MF", model, "allMF", 0.9, Project, 107, GOSeSimdata.MF)

## selected for CC
load("CAD_xl315_0001-PTvsPUT_GOSeSim_allCC_CC.RData")
G6.CC <- reducedTerms[grep("actin|intermediate filament|microtubule", reducedTerms$term),] 
## 10items, C2-8, C6-1, C8-1
G7.CC <- reducedTerms[grep("microtubule|dynactin", reducedTerms$term),]  ## C2-2 
load("CAD_xl315_0001-PTvsPUT_GOSeSim_allMF_MF.RData")
G1.MF <- reducedTerms[grep("RNA", reducedTerms$term),]  ## 15 items, C1-10, C2-2, C5-2, C7-1
G6.MF <- reducedTerms[grep("actin|intermediate filament|microtubule", reducedTerms$term),] ## 7 items, C2-5, C7-2
G11.MF <- reducedTerms[grep("transport", reducedTerms$term),] ## C15-3
PCM.sel  <- list(G6.CC, G7.CC, G1.MF, G6.MF, G11.MF)
PCM.selN <- paste0("G", c("6C","7C","1M","6M","11M")," ", c("Cytoskeleton",  "Microtubule", "RNA", "Cytoskeleton", "Transport"))
write.xlsx(PCM.sel[[1]], file="Selected_GOCCMF_Group_10_Jan_2022.xlsx", sheetName = PCM.selN[1], 
           col.names = TRUE, row.names = TRUE, append = FALSE)
for(i in 2:5){
  print(i)
  write.xlsx(PCM.sel[[i]], file="Selected_GOCCMF_Group_10_Jan_2022.xlsx", sheetName = PCM.selN[i], 
             col.names = TRUE, row.names = TRUE, append = TRUE)
}


## Discuss and selected GOs

G1.BP <- G1.BP[G1.BP$go!="GO:0050658"&G1.BP$go!="GO:0071359", ]
G5.BP <- G5.BP[G5.BP$go!="GO:0006165"&G1.BP$go!="GO:0046939", ]
G6.BP <- G5.BP[G5.BP$go!="GO:0098974", ]
G9.BP <- G5.BP[G5.BP$go!="GO:0000281"&G1.BP$go!="GO:0002312", ]

message("+---- Merged related grouped pathways as above defination and generate barplot -----+")

G1.subdat <- newPaths.all[newPaths.all$ID%in%c(G1.BP$go, G1.MF$go),]
G2.subdat <- newPaths.all[newPaths.all$ID%in%G2.BP$go,]
G3.subdat <- newPaths.all[newPaths.all$ID%in%G3.BP$go,]
G4.subdat <- newPaths.all[newPaths.all$ID%in%G4.BP$go,]
G5.subdat <- newPaths.all[newPaths.all$ID%in%G5.BP$go,]
G6.subdat <- newPaths.all[newPaths.all$ID%in%c(G6.BP$go,G6.MF$go,G6.CC$go),]
G7.subdat <- newPaths.all[newPaths.all$ID%in%c(G7.BP$go,G7.CC$go),]
G8.subdat <- newPaths.all[newPaths.all$ID%in%G8.BP$go,]
G9.subdat <- newPaths.all[newPaths.all$ID%in%G9.BP$go,]
G10.subdat <- newPaths.all[newPaths.all$ID%in%G10.BP$go,]
G11.subdat <- newPaths.all[newPaths.all$ID%in%G11.BP$go,]
Term  <- c("RNA", "Metabolic ", "Cell ",
           "Circadian", "Phosphorylation", "Cytoskeleton", 
           "Microtubule", "Localization", "Immune",
           "Inflammasome", "Transport")
## Reviewed on 24/Jan/2022, remove the G4.subdat
GOsel.dat   <- list(newPaths.all, G1.subdat,G2.subdat,G3.subdat,G4.subdat,G5.subdat,G6.subdat,
                    G7.subdat,G8.subdat,G9.subdat,G10.subdat,G11.subdat)
GOsel.datN  <- c("Enrich GO Pathways", paste0("G", c(1:11), " ", Term)) 
write.xlsx(GOsel.dat[[1]], file="Supplementary_GeneOntology_Group_10_Jan_2022.xlsx", sheetName = GOsel.datN[1], 
           col.names = TRUE, row.names = TRUE, append = FALSE)
for(i in 2:12){
  print(i)
  write.xlsx(GOsel.dat[[i]], file="Supplementary_GeneOntology_Group_10_Jan_2022.xlsx", sheetName = GOsel.datN[i], 
             col.names = TRUE, row.names = TRUE, append = TRUE)
}

BP_list_up  <- list();
BP_list_tot <- list()
newGOsel.dat <- list(G1.subdat,G2.subdat,G3.subdat,G5.subdat,G6.subdat,
                     G7.subdat,G8.subdat,G9.subdat,G10.subdat,G11.subdat)
for (i in 1:length(newGOsel.dat)){
  pgenes  <- unique(unlist(lapply(newGOsel.dat[[i]]$geneID, function(x) strsplit(x, split="/")[[1]])))
  df_tmp  <- Presdf[Presdf$SYMBOL%in%pgenes,]
  tmp_up  <- length(subset(df_tmp$L2FC, df_tmp$L2FC > 0))
  tmp_tot <- dim(df_tmp)[1]
  BP_list_up[[i]]  <- tmp_up
  BP_list_tot[[i]] <- tmp_tot
}

UP    <- as.numeric(unlist(BP_list_up))
DOWN  <- as.numeric(unlist(BP_list_tot)) - UP
DOWN  <- -DOWN
Count <- as.numeric(unlist(BP_list_tot))
NTerm  <- c("RNA processing,stablization \nand localization", "Metabolic process", "Cell shape and cell polarity",
           "Protein phosphorylation", "Cytoskeleton organization", 
           "Microtubule organization and \nmicrotubule-based motor activity", "Protein localization", "Immune response",
           "Regulation of PRRs", "Intracellular transport")
MerPath.dat <- data.frame(Term=NTerm, Count=Count, UP=UP, DOWN=DOWN)

MerPath.dat_molten <- melt(MerPath.dat, id.vars=c("Term","Count") )
MerPath.dat_molten <- MerPath.dat_molten[order(-MerPath.dat_molten$Count),]
brks <- seq(-80, 100, 20)
lbls = as.character(c(seq(80, 0, -20), seq(20, 100, 20)))

#####------------------------------------------------------------Fig 3f---------------------------------------##
pdf("Protein_190920_N7216_F1252_reviewed_withBP_Grouped_pathways_barplot_17_Jan_2022.pdf", width=8, height=5)
barplotSum <- ggplot(MerPath.dat_molten, aes(x = reorder(Term, Count), y = value, fill = variable)) +   
  geom_bar(stat = "identity", width = 0.8) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls,
                     limits=c(-80,100)) + # Labels
  coord_flip() +  
  xlab("") +
  ylab("Number of genes") +
  guides(fill=FALSE) +
  theme(plot.title = element_text(hjust = .2),
        axis.ticks = element_blank()
  ) +   
  scale_fill_manual(values = c("red", "blue"))+  # Color palette
  theme_bw() +
  theme_update(axis.title.x = element_text(size=14, face= "bold"),
               axis.text.x = element_text(size=14, face="bold"),
               axis.title.y = element_text(size=14, face= "bold"),
               axis.text.y.left = element_text(size=14, face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white", colour = NA)) 

barplotSum

dev.off()
#####------------------------------------------------------------------------------------------------------##


message("+-----Heatmap with the sigDEs relating to microtubule pathway   ------------------------- +")

library("ComplexHeatmap")
library("seriation")
breaksList.deg           <- seq(6, 15, by = 0.001)
ScaleCols.deg            <- colorRampPalette(colors = c("gold1", "brown2"))(length(breaksList.deg))
ScaleCols.l2             <- colorRamp2(c(-2.5, 0, 1), c("purple4", "white", "darkgreen"))

phgenes        <- unique(unlist(pgenes))
hlt.matSig     <- PTvsPUT_oriR[PTvsPUT_oriR$GeneSymbol%in%phgenes, 
                               c("GeneSymbol", paste0(rep(c("PT", "PUT"), each=3),"_", rep(c(1:3), length=6)))]
heatmat.Sig    <- hlt.matSig
rownames(heatmat.Sig)  <- hlt.matSig[,1]
heatmat.Sig    <- hlt.matSig[,-1]


basemeans      <- Presdf[Presdf$SYMBOL%in%rownames(heatmat.Sig),c("SYMBOL", "L2FC")]
L2FC           <- basemeans[order(-basemeans$L2FC), ]
L2FC           <- as.matrix(L2FC[,"L2FC"], ncol=1)
colnames(L2FC) <- "L2FC"
rownames(L2FC) <- basemeans$SYMBOL

## match the two matrix rownames order.
genomic_idx    <- match(rownames(L2FC), rownames(heatmat.Sig))
heatmat.Sig    <- heatmat.Sig[genomic_idx,]

condition   <- rep(c("PT", "PUT"), c(3,3))
ha  = HeatmapAnnotation(Condition = condition,  col = list(Condition = c("PT" = "red", "PUT"= "blue")), 
                        annotation_name_side = "left")

###
o2mat     <- seriate(dist(t(heatmat.Sig)), method = "GW")

pht_heatSig = Heatmap(as.matrix(heatmat.Sig), col = heat.colors(1000),  
                      name = "Heatmap", show_row_names=T,
                      row_names_side="left",
                      show_column_names = T, width = unit(9, "cm"),
                      heatmap_legend_param =list(title = "Normalised Abundance",
                                                 title_position = "leftcenter-rot",
                                                 title_gp = gpar(fontsize = 10, fontface="bold"),
                                                 labels_gp = gpar(fontsize = 10),
                                                 legend_width = unit(0.8,"cm"),
                                                 legend_height= unit(4, "cm")),
                      top_annotation = ha,
                      cluster_rows = F, 
                      show_column_dend = T,
                      column_title=NULL,
                      row_title_rot = 0,
                      show_row_dend=F,
                      column_names_rot = 45,
                      cluster_columns = reorder(as.dendrogram(o2mat[[1]]), c(4:6,1:3)), 
                      column_names_gp = gpar(fontsize = 10, col = rep(c("red","blue"), c(3,3))))+
  Heatmap(L2FC, name = "l2FC", col=ScaleCols.l2, width = unit(5, "mm"),
          heatmap_legend_param = list(title = "Log2FoldChange",title_position = "leftcenter-rot",
                                      legend_width = unit(0.8,"cm"),
                                      legend_height= unit(4, "cm")), 
          show_row_names=F,row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 9,  fontface="bold"),
          cluster_rows=T, show_row_dend=F) 


pdf("PTvsPUT_microtubule_selMarkers_withBP_Heatmap_23_Dec_2021.pdf", width=7, height= 10)
draw(pht_heatSig, heatmap_legend_side = "right", merge_legend=T)
dev.off()

## plot the mean of PT and PUT replicates
heatmat.SigM <- cbind(rowMeans(heatmat.Sig[,1:3]), rowMeans(heatmat.Sig[,4:6]))
colnames(heatmat.SigM)  <- c("PT", "PUT")
haM  = HeatmapAnnotation(Condition = c("PT", "PUT"),  col = list(Condition = c("PT" = "red", "PUT"= "blue")), 
                         annotation_name_side = "left")
o2matM     <- seriate(dist(t(heatmat.SigM)), method = "GW")
pht_heatSigM = Heatmap(as.matrix(heatmat.SigM), col = heat.colors(1000),  
                       name = "Heatmap", show_row_names=T,
                       row_names_side="left",
                       show_column_names = T, width = unit(4, "cm"),
                       heatmap_legend_param =list(title = "Normalised exprs.",
                                                  title_position = "leftcenter-rot",
                                                  title_gp = gpar(fontsize = 10, fontface="bold"),
                                                  labels_gp = gpar(fontsize = 10),
                                                  legend_width = unit(0.8,"cm"),
                                                  legend_height= unit(3, "cm")),
                       top_annotation = haM,
                       cluster_rows = F, 
                       show_column_dend = F,
                       column_title=NULL,
                       row_title_rot = 0,
                       show_row_dend=F,
                       column_names_rot = 45,
                       cluster_columns = reorder(as.dendrogram(o2matM[[1]]), c(2,1)), 
                       column_names_gp = gpar(fontsize = 10, col = rep(c("red","blue"), c(1,1))))+
  Heatmap(L2FC, name = "l2FC", col=ScaleCols.l2, width = unit(5, "mm"),
          heatmap_legend_param = list(title = "Log2FoldChange",title_position = "leftcenter-rot",
                                      legend_width = unit(0.8,"cm"),
                                      legend_height= unit(4, "cm")), 
          show_row_names=F,row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 9,  fontface="bold"),
          cluster_rows=T, show_row_dend=F) 
pdf("PTvsPUT_microtubule_selMarkers_withBP_Mean_Heatmap_23_Dec_2021.pdf", width=6, height= 10)
draw(pht_heatSigM, heatmap_legend_side = "right", merge_legend=T)
dev.off()


message("+-----Heatmap with all genes relating to microtubule pathway, sep sig/nonsig, revised 25Jan   -----------+")

G1.mKey <- c("microtubule", "microtubule cytoskeleton", "cytoplasmic microtubule", "microtubule-based process")
G3.mKey <- c("microtubule cytoskeleton organization", "regulation of microtubule cytoskeleton organization",
             "cytoplasmic microtubule organization")
G2.mKey <- c("microtubule organizing center", "microtubule anchoring at centrosome")
G5.mKey <- c("positive regulation of microtubule nucleation")
#G5.mKey <- c("establishment or maintenance of microtubule cytoskeleton polarity")
#G6.mKey <- c("positive regulation of microtubule polymerization")
G4.mKey <- c("microtubule motor activity", "ATP-dependent microtubule motor activity",
             "mitochondrion transport along microtubule", "organelle transport along microtubule",
             "regulation of organelle transport along microtubule")
#G8.mKey <- c("minus-end-drirected vesicle transport along microtubule",
#             "ATP-dependent microtubule motor activity, minus-end-directed")

load("Ensembl_gene_id_symbol_entrezid_goid_name1006_GRCm39.RData")
colnames(ensEMBL2id.mer)[2] <- c("GeneSymbol")
PTvsPUT_oriR.merGO        <- merge(PTvsPUT_oriR, ensEMBL2id.mer, by = "GeneSymbol")
PTvsPUT_oriR.microtubu    <- PTvsPUT_oriR.merGO[grep("microtubule", PTvsPUT_oriR.merGO$name_1006),]
PTvsPUT_oriR.microtubu    <- merge(PTvsPUT_oriR.microtubu, Qgo_anno, by = "go_id")[,-c(16,19:20)]
PTvsPUT_oriR.microtubu    <- unique(PTvsPUT_oriR.microtubu)
## 219 items, involve 99 genes, 99 proteins, 52 GO relating to microtubules.
write.csv(PTvsPUT_oriR.microtubu, file = "PTvsPUT_N1252_withBP_GOmicrotubu_Relating_Gene_Summary_23_Dec_2021.csv", row.names=F)

unique(PTvsPUT_oriR.microtubu$GeneSymbol[which(PTvsPUT_oriR.microtubu$adj.P.Val<0.05)]) 
## 54 genes have padj < 0.05

message("+---- Choose the subgroup of microtubule relating as above 5 groups----------+")

PTvsPUT_oriR.mSig <- subset(PTvsPUT_oriR.microtubu,adj.P.Val<0.05) ## 125
G1.micro <- PTvsPUT_oriR.mSig[PTvsPUT_oriR.mSig$name_1006%in%G1.mKey,]
G2.micro <- PTvsPUT_oriR.mSig[PTvsPUT_oriR.mSig$name_1006%in%G2.mKey,]
G3.micro <- PTvsPUT_oriR.mSig[PTvsPUT_oriR.mSig$name_1006%in%G3.mKey,]
G4.micro <- PTvsPUT_oriR.mSig[PTvsPUT_oriR.mSig$name_1006%in%G4.mKey,]
G5.micro <- PTvsPUT_oriR.mSig[PTvsPUT_oriR.mSig$name_1006%in%G5.mKey,]


message("+----generate heatmap data and then plot it --------+")

microGN <- unique(c(G1.micro$GeneSymbol,G2.micro$GeneSymbol, G3.micro$GeneSymbol,
                    G4.micro$GeneSymbol,G5.micro$GeneSymbol))

microDL <- list(G1.micro[,c("GeneSymbol", "log2FC")],
                G2.micro[,c("GeneSymbol", "log2FC")],
                G3.micro[,c("GeneSymbol", "log2FC")],
                G4.micro[,c("GeneSymbol", "log2FC")],
                G5.micro[,c("GeneSymbol", "log2FC")])
multi_mer <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "GeneSymbol", all = TRUE),
                    microDL)
multi_mer <- unique(multi_mer)
rownames(multi_mer) <- multi_mer$GeneSymbol
colnames(multi_mer) <- c("GeneSymbol", "G1", "G2", "G3", "G4", "G5")
microHeat.Mat <- multi_mer[,-1]
microHeat.Mat[is.na(microHeat.Mat)] <- 0
roat.Mat <- t(as.matrix(microHeat.Mat))
rownames(roat.Mat) <- c("MT and MT-based process",
                        "MT organizing center",
                        "MT cytoskeleton organization",
                        "MT-based movement \nand motor activity",
                        "Positive regulation of \nMT nucleation")
colnames(roat.Mat) <- multi_mer$GeneSymbol
roat.Mat <- roat.Mat[,-c(26,36)]
ProtNameF <- read.xlsx("Gene_Protein_Name_Mar_2022.xlsx", sheetIndex=1)
colnames(roat.Mat)%in%ProtNameF[,1]
colnames(roat.Mat) <- ProtNameF[,2]

##----------------------------SFig5e-------------------------------------###
pdf("PTvsPUT_selMicrotubule_G1-5_Heatamap_jan_2022.pdf", width= 10, height= 3)
microHeat <- Heatmap(roat.Mat, name = "log2FC", show_row_names=T,
                       row_names_side="left", cluster_rows=F, show_row_dend = F, 
                       cluster_columns = T, show_column_dend = F, 
                       row_names_gp = gpar(fontsize = 12, fontface = "bold"),
                       column_names_gp = gpar(fontsize = 9, fontface = "bold"))
                      #,
                      # width = ncol(microHeat.Mat)*unit(30, "mm"), 
                      # height = nrow(microHeat.Mat)*unit(0.8, "mm"))
print(microHeat)
dev.off()
##-----------------------------------------------------------------------###


message("+-----localisation analysis of existing 1251 proteins   -----------+")

## convert mouse gene to human gene 
require("biomaRt")
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ss    <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", mart = mouse, 
                values=PTvsPUT_oriR$GeneSymbol,
                attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
colnames(ss) <- c("GeneSymbol", "Gene.name")

subcellular <- read.table("subcellular_location.tsv", sep="\t", header=T)
subcellular.over <- merge(ss, subcellular, by = "Gene.name")
subcellular.sel  <- subcellular.over[,c("Gene.name", "GeneSymbol", "Reliability", "Main.location")]
subcellular.sel  <- unique(subcellular.sel)

## merge the subcellular.over with the counts from PT and PUT samples
colnames(PTvsPUT_oriR)[2] <- "GeneSymbol"
PTvsPUT_locs <- merge(subcellular.sel, PTvsPUT_oriR, by = "GeneSymbol")
## 1019
PTvsPUT_locs <- PTvsPUT_locs[-which(duplicated(PTvsPUT_locs$GeneSymbol)==T),]
## 999, "Btk"   "Cdk1"  "Cenpe" "Ddx3x" "Hdac1" "Pcm1"  "Pcnt" included.
write.csv(PTvsPUT_locs, file = "PTvsPUT_localisation_N999_23_Dec_2021.csv", row.names=F)

message("+----- Density plot to show the difference between PT and PUT -----------+")

Actinfilam <- PTvsPUT_locs[grep("Actin filaments", PTvsPUT_locs$Main.location),]
Celljunc   <- PTvsPUT_locs[grep("Cell Junctions", PTvsPUT_locs$Main.location),]
Cytokbridge<- PTvsPUT_locs[grep("Cytokinetic bridge", PTvsPUT_locs$Main.location),]
Cytoplasbod<- PTvsPUT_locs[grep("Cytoplasmic bodies", PTvsPUT_locs$Main.location),]
FAS        <- PTvsPUT_locs[grep("Focal adhesion sites", PTvsPUT_locs$Main.location),]
Lipiddrop  <- PTvsPUT_locs[grep("Lipid droplets", PTvsPUT_locs$Main.location),]
Peroxisomes<- PTvsPUT_locs[grep("Peroxisomes", PTvsPUT_locs$Main.location),]
NuclearM    <- PTvsPUT_locs[grep("Nuclear membrane", PTvsPUT_locs$Main.location),] 
NuclearS    <- PTvsPUT_locs[grep("Nuclear speckles", PTvsPUT_locs$Main.location),] 
Lysosomes   <- PTvsPUT_locs[grep("Lysosomes", PTvsPUT_locs$Main.location),] 
Mitoticchr <- PTvsPUT_locs[grep("Mitotic chromosome", PTvsPUT_locs$Main.location),]

SelList1   <- list(Actinfilam, Celljunc, Cytokbridge, Cytoplasbod, FAS, Lipiddrop, Peroxisomes,
                NuclearM, NuclearS, Lysosomes, Mitoticchr)
SelName1   <- c("Actin filaments", "Cell Junctions", "Cytokinetic bridge", "Cytoplasmic bodies",
                "Focal adhesion sites", "Lipid droplets", "Peroxisomes", "Nuclear membrane",
                "Nuclear speckles", "Lysosomes", "Mitotic chromosome")
SelLen1    <- unlist(lapply(SelList1, function(x) dim(x)[1]))
## [1]  7  7  4  5  8  3  4 11 44  5 6

Centro     <- PTvsPUT_locs[grep("Centro", PTvsPUT_locs$Main.location),]
##Centriolar satellite,Centrosome
Golgi      <- PTvsPUT_locs[grep("Golgi", PTvsPUT_locs$Main.location),]
Mitochondria    <- PTvsPUT_locs[grep("Mitochondria", PTvsPUT_locs$Main.location),]
Cytosol    <- PTvsPUT_locs[grep("Cytosol", PTvsPUT_locs$Main.location),]
ER         <- PTvsPUT_locs[grep("Endoplasmic reticulum", PTvsPUT_locs$Main.location),]
Endosomes  <- PTvsPUT_locs[grep("Endosomes", PTvsPUT_locs$Main.location),]
Interfilam <- PTvsPUT_locs[grep("Intermediate filaments", PTvsPUT_locs$Main.location),]
Micrtotubu <- PTvsPUT_locs[grep("Microtubules", PTvsPUT_locs$Main.location),]
Nucleoli   <- PTvsPUT_locs[grep("Nucleoli", PTvsPUT_locs$Main.location),] 
## incl. Nucleoli and Nucleoli rim &Nucleoli fibrillar center
Nucleoplasm <- PTvsPUT_locs[grep("Nucleoplasm", PTvsPUT_locs$Main.location),]
Nuclearbod <- PTvsPUT_locs[grep("Nuclear bodies", PTvsPUT_locs$Main.location),]
Plasmamemb <- PTvsPUT_locs[grep("Plasma membrane", PTvsPUT_locs$Main.location),]
Vesicles   <- PTvsPUT_locs[grep("Vesicles", PTvsPUT_locs$Main.location),]

SelList2   <- list(Centro, Golgi, Mitochondria, Cytosol, ER, Endosomes, Interfilam, Micrtotubu,
                Nucleoli, Nucleoplasm, Mitoticchr, Nuclearbod, Plasmamemb, Vesicles)
SelName2   <- c("Centro", "Golgi", "Mitochondria", "Cytosol", "Endoplasmic reticulum", 
                "Endosomes", "Intermediate filaments", "Microtubules", "Nucleoli",
                "Nucleoplasm",  "Nuclear bodies","Plasma membrane",
                "Vesicles")
SelLen2    <- unlist(lapply(SelList2, function(x) dim(x)[1]))
## [1]  11  28  91 362  70   7  13  15  81 330 15  77  76

Micrtotubu$PT_mean <- rowMeans(Micrtotubu[,grep("PT", colnames(Micrtotubu))])
Micrtotubu$PUT_mean <- rowMeans(Micrtotubu[,grep("PUT", colnames(Micrtotubu))])
MicPLK1_ratio <- mean(Micrtotubu$PT_mean)/mean(Micrtotubu$PUT_mean)

Centro$PT_mean <- rowMeans(Centro[,grep("PT", colnames(Centro))])
Centro$PUT_mean <- rowMeans(Centro[,grep("PUT", colnames(Centro))])
CenPLK1_ratio <- mean(Centro$PT_mean)/mean(Centro$PUT_mean)
message("+--- Generate the data for density plot, -----------------------+")
data_fun <- function(input, inputLen, inputName){
  PTmeans  <- rowMeans(input[,grep("PT_", colnames(input))])
  PUTmeans <- rowMeans(input[,grep("PUT_", colnames(input))])
  teststat <- t.test(PTmeans, PUTmeans, paired = TRUE, alternative = "greater")$p.value
  plab     <- paste0("p=",round(teststat, digits=4))
  Group  <- rep(c("Activated", "Primed"), each=inputLen)
  Intensity <- c(PTmeans, PUTmeans)
  testdens  <- data.frame(Group=Group, Intensity=Intensity)
  mu        <- plyr::ddply(testdens, "Group", summarise, grp.mean=mean(Intensity))
  dens.plt  <- ggplot(testdens, aes(x=Intensity, color=Group)) +
               geom_density()+
               scale_color_manual(values = c("red2", "tan2"))+
               geom_vline(data=mu, aes(xintercept=grp.mean, color=Group),
                           linetype="dashed") +
               ggtitle(paste0(inputName, ":", plab)) +
               theme_bw() +
               theme(axis.title.x = element_text(size=14, face= "bold"),
                     axis.text.x = element_text(size=14, face="bold"),
                     axis.title.y = element_text(size=14, face= "bold"),
                     axis.text.y.left = element_text(size=14, face="bold"))+
               theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     axis.line = element_line(colour = "black"),
                     panel.background = element_rect(fill = "white", colour = NA), 
                     plot.title = element_text(size=12, face="bold")) 
               
  dens.plt
}

dens.pltList2 <- list()
inputNames2   <- paste0(SelName2, ":N=", SelLen2)
for(i in 1:length(SelName2)){
  dens.pltList2[[i]] <- data_fun(SelList2[[i]], SelLen2[i], inputNames2[i])
  dens.pltList2
}

dens.pltList1 <- list()
inputNames1   <- paste0(SelName1, ":N=", SelLen1)
for(i in 1:length(SelName1)){
  dens.pltList1[[i]] <- data_fun(SelList1[[i]], SelLen1[i], inputNames1[i])
  dens.pltList1
}
###--------------------------SFig5D--------------------------------------------------------------------#####
pdf("CAD_xl315_0001-Localization_groupDensity_PTvsPUT_plot_Selected_N4_25_Jan_2022.pdf", width = 8, height = 6)
prow <- plot_grid(dens.pltList1[[1]]+ theme(legend.position="none"), 
                dens.pltList2[[7]]+ theme(legend.position="none"), 
                dens.pltList2[[8]]+ theme(legend.position="none"), 
                dens.pltList2[[1]]+ theme(legend.position="none"), nrow=2, ncol=2, byrow=T)
legend_b <- get_legend(
  dens.pltList1[[1]] + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))

dev.off()
###---------------------------------------------------------------------------------------------------#####



pdf("CAD_xl315_0001-Localization_groupDensity_PTvsPUT_plot_Selected_N13_14_Jan_2022.pdf", width = 10, height = 6)
plot_grid(dens.pltList2[[1]], dens.pltList2[[2]], dens.pltList2[[3]], dens.pltList2[[4]], nrow=2)
plot_grid(dens.pltList2[[5]], dens.pltList2[[6]], dens.pltList2[[7]], dens.pltList2[[8]], nrow=2)
plot_grid(dens.pltList2[[9]], dens.pltList2[[10]], dens.pltList2[[11]], dens.pltList2[[12]], nrow=2)
plot_grid(dens.pltList2[[13]], nrow=2, ncol=2)
dev.off()
pdf("CAD_xl315_0001-Localization_groupDensity_PTvsPUT_plot_NonSelected_N11_14_Jan_2022.pdf", width = 10, height = 6)
plot_grid(dens.pltList1[[1]], dens.pltList1[[2]], dens.pltList1[[3]], dens.pltList1[[4]], nrow=2)
plot_grid(dens.pltList1[[5]], dens.pltList1[[6]], dens.pltList1[[7]], dens.pltList1[[8]], nrow=2)
plot_grid(dens.pltList1[[9]], dens.pltList1[[10]], dens.pltList1[[11]], nrow=2, ncol=2)
dev.off()



message("+-----------Plk1 normalisation compared to get the selected genes L2FC and padj ----+")

MSnSet_normR_batch_nplk1  <- MSnSet_normR_batch
Plk1.val <- exprs(MSnSet_normR_batch_nplk1)[grep("Q07832", rownames(exprs(MSnSet_normR_batch_nplk1))),grep("PT|PUT", colnames(exprs(MSnSet_normR_batch_nplk1)))]
exprs(MSnSet_normR_batch_nplk1)[,grep("PT|PUT", colnames(exprs(MSnSet_normR_batch_nplk1)))] <- 
  sweep(exprs(MSnSet_normR_batch_nplk1)[,grep("PT|PUT", colnames(exprs(MSnSet_normR_batch_nplk1)))], 2, as.numeric(Plk1.val), "/")

MSnSetObj.plk1 <- MSnSet_normR_batch_nplk1[,grep("PT|PUT", colnames(exprs(MSnSet_normR_batch_nplk1)))]
plotDat <- exprs(MSnSetObj.plk1) %>% 
  as.data.frame() %>% pivot_longer(names_to = "SampleName", 
                                   values_to = "Intensity", everything()) %>% 
  mutate(logInt = log2(Intensity)) %>% 
  filter(is.finite(logInt)) %>% left_join(pData(MSnSetObj.plk1), 
                                          "SampleName") 
plotDat$SampleName  <- paste0(plotDat$SampleGroup, "_", plotDat$BioRep)
plotDat$SampleName  <- factor(plotDat$SampleName,
                              levels = c("Primed_BASU-PLK1_1", "Primed_BASU-PLK1_2", "Primed_BASU-PLK1_3",
                                         "Activated_BASU-PLK1_1","Activated_BASU-PLK1_2","Activated_BASU-PLK1_3"),
                              ordered = TRUE)
normBatchBox <- ggplot(plotDat, aes(x =SampleName, y = logInt, fill=SampleGroup)) + 
  geom_boxplot()  + 
  scale_fill_manual(values = c("red2", "tan2")) +
  labs(x = "", y = "log2(Intensity)", title = title) + 
  theme_bw() + theme(axis.text.x = element_text(angle = 45,  hjust = 1, size=8, face="bold"),
                     axis.title.y = element_text(size=12, face= "bold"),
                     axis.text.y.left = element_text(size=12, face="bold")) + 
  guides(fill = FALSE) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_blank()) 
  

pdf("PTvsPUT-NonNorm_Norm_NormBatch_Logintens_normPlk1_boxplot_26_jan_2022.pdf", width=3, height=4.5)
normBatchBox
dev.off()



##---------------FINISH 23Dec,2021-----------------------------------------------------------##















