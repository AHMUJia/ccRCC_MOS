#before run the code, you should install all the packages, espically "MOVICS"
#install "MOVICS" through this page: https://github.com/xlucpu/MOVICS

tumor.path <- "set your own path";setwd(tumor.path) #create dir
data.path   <- file.path(tumor.path, "InputData")
fig1.path    <- file.path(tumor.path, "Figure1")
fig2.path    <- file.path(tumor.path, "Figure2")
fig5.path    <- file.path(tumor.path, "Figure5")
fig6.path    <- file.path(tumor.path, "Figure6")
figS.path    <- file.path(tumor.path, "FigureS")
tabl.path    <- file.path(tumor.path, "Tables")
Scripts    <- file.path(tumor.path, "Scripts")

if (!file.exists(tumor.path)) { dir.create(tumor.path) }
if (!file.exists(data.path)) { dir.create(data.path) }
if (!file.exists(fig1.path)) { dir.create(fig1.path) }
if (!file.exists(fig2.path)) { dir.create(fig2.path) }
if (!file.exists(fig5.path)) { dir.create(fig5.path) }
if (!file.exists(fig6.path)) { dir.create(fig6.path) }
if (!file.exists(figS.path)) { dir.create(figS.path) }
if (!file.exists(tabl.path)) { dir.create(tabl.path) }
if (!file.exists(Scripts)) { dir.create(Scripts) }

library(MOVICS)
library(aplot)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(cowplot)
library("wesanderson")
#
#--------------------------#
#---------Figure S1--------#
#--------------------------#
load("./InputData/TCGA_KIRC.fpkm.rda")
load("./InputData/mo.data.rda")
optk.KIRC <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,F,T), # note: the 4th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.path    = figS.path,
                         fig.name    = "Figure S1A.CLUSTER NUMBER OF TCGA-KIRC")

# load R data
load("./InputData/moic.res.list.rda")

cmoic.KIRC <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               fig.path    = figS.path,
                               distance      = "euclidean",
                               linkage       = "average")

getSilhouette(sil      = cmoic.KIRC$sil, # a sil object returned by getConsensusMOIC()
              fig.path = figS.path,
              fig.name = "Figure S1B.SILHOUETTE",
              height   = 5.5,
              width    = 5)

#--------------------------#
#---------Figure 1--------#
#--------------------------#

indata <- mo.data
indata$meth <- log2(indata$meth / (1 - indata$meth))

# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,T,F)) # no scale for mutation

feat1  <- feat[which(feat$dataset == "mRNA"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "cna"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "meth"),][1:10,"feature"]
feat5  <- feat[which(feat$dataset == "mut"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4,feat5)

#set the colors for mult-omics
mRNA.col   <- c("#31B29E", "black","#EB6C5A")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
cna.col  <-  c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
meth.col   <- c("#0099CC", "white","#CC0033")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, cna.col,meth.col, mut.col)

surv.info<-read.table(file.path(data.path,"KIRC.surv.txt"),header=T,sep="\t",row.names=1,check.names=F)
# comprehensive heatmap (may take a while)

annCol    <- surv.info[,c("Age","Grade","Stage","New_tumor"), drop = FALSE]
# generate corresponding colors for sample annotation
annColors <- list(Age    = circlize::colorRamp2(breaks = c(min(annCol$Age),
                                                           max(annCol$Age)), 
                                                colors = c("#7FC97F",  "#E00115")),
                  New_tumor  = c("TUMOR FREE" = "#27B1EB",
                                 "WITH TUMOR"   = "#F26622",
                                 "unknow" = "#C2BEBC"),
                  Grade  = c("G1" = "#AAD1DF",
                             "G2"   = "#55A3C0",
                             "G3"   = "#004E6B",
                             "G4"   = "#002735",
                             "unknow" = "#C2BEBC"),
                  Stage = c("Stage I"    = "#C3AEA5",
                            "Stage II"    = "#A6877A",
                            "Stage III"    = "#7B4B38",
                            "Stage IV"    = "#513225", 
                            "unknow"    = "#C2BEBC"))

getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA","lncRNA","CNA","Methylation","Mutation"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("mRNA.expr","lncRNA.expr","CNA","M value","Mutated"),
             clust.res     = cmoic.KIRC$clust.res, # consensusMOIC results
             clust.dend    = NULL, # show no dendrogram for samples
             show.rownames = c(F,F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             show.row.dend = c(F,F,F,F,F), # show no dendrogram for features
             annRow        = annRow, # no selected features
             color         = col.list,
             annCol        = annCol, # annotation for samples
             annColors     = annColors, # annotation color
             width         = 10, # width of each subheatmap
             height        = 3, # height of each subheatmap
             fig.name      = "Figure 1A.COMPREHENSIVE HEATMAP OF CONSENSUSMOIC",
             fig.path      = fig1.path)



surv.info=read.table(file.path(data.path,"KIRC.surv.txt"),header=T,sep="\t",row.names=1,check.names=F)

surv.info$futime<-surv.info$PFI.time
surv.info$fustat<-surv.info$PFI

surv.KIRC.PFI <- compSurv(moic.res         = cmoic.KIRC,
                      surv.info        = surv.info,
                      surv.cut         = 150,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h",
                      fig.path         = fig1.path,
                      fig.name         = "Figure 1C. PFI KAPLAN-MEIER CURVE")
surv.info$futime<-surv.info$OS.time
surv.info$fustat<-surv.info$OS

surv.KIRC.OS <- compSurv(moic.res         = cmoic.KIRC,
                      surv.info        = surv.info,
                      surv.cut         = 150,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h",
                      fig.path         = fig1.path,
                      fig.name         = "Figure 1B. OS KAPLAN-MEIER CURVE")

#compare the different clinical features among subtypes
clin.KIRC <- compClinvar(moic.res      = cmoic.KIRC,
                         var2comp      = surv.info, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("OS","PFI"), # features that are considered categorical variables
                         #nonnormalVars = "futime", # feature(s) that are considered using nonparametric test
                         #exactVars     = "pstage", # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         tab.name      = "Table 1.SUMMARIZATION OF CLINICAL FEATURES",
                         res.path      = tabl.path)

#--------------------------#
#---------Figure S2--------#
#--------------------------#

#------------------DEGs  different pathways
#load the mRNA expression of TCGA-KIRC

MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)

gsea.up <- runGSEA(moic.res     = cmoic.KIRC,
                   dea.method   = "limma", # name of DEA method
                   prefix       = "TCGA-KIRC", # MUST be the same of argument in runDEA()
                   dat.path     = data.path, # path of DEA files
                   res.path     = figS.path, # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = fpkm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.15, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "Figure S2A. Go UPREGULATED PATHWAY HEATMAP")


library("clusterProfiler")
options(connectionObserver = NULL)
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

print(gsea.up$gsea.list$CS1[1:6,])

Go1<- pairwise_termsim(gsea.up$gsea.list$CS1)

pdf(file=file.path(figS.path,"Go-net-CS1.pdf"),width = 10,height = 8)

cnetplot(Go1,
         showCategory = 5, 
         #foldChange = ID, 
         node_label="all", 
         colorEdge = TRUE,
         categorySize="pvalue")


dev.off()

Go2<- pairwise_termsim(gsea.up$gsea.list$CS2)

pdf(file=file.path(figS.path,"Go-net-CS2.pdf"),width = 10,height = 8)

cnetplot(Go2,
         showCategory = 5, 
         #foldChange = ID, 
         node_label="all", 
         colorEdge = TRUE,
         categorySize="pvalue")

dev.off()

Go3<- pairwise_termsim(gsea.up$gsea.list$CS3)

pdf(file=file.path(figS.path,"Go-net-CS3.pdf"),width = 14,height = 8)

cnetplot(Go3,
         showCategory = 5, 
         #foldChange = ID, 
         node_label="all", 
         colorEdge = TRUE,
         categorySize="pvalue")

dev.off()

#--------------------------#
#----------Figure 2--------#
#--------------------------#

####Figure 2B######
library(MOVICS)

GSET.FILE <- "./InputData/immunesuppress.gmt"

gsva.res <- 
  runGSVA(moic.res      = cmoic.KIRC,
          norm.expr     = fpkm,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          annCol        = annCol,
          annColors     = annColors,
          color         =c("#0099df", "white","#CC0033"),
          fig.path      = fig2.path,
          fig.name      = "Figure 2B. immune GENE SETS OF INTEREST HEATMAP",
          height        = 8,
          width         = 10)

####Figure 2C-D######
#--------immune checkpoints 3 group
library(ggpubr)
finalsam<-cmoic.KIRC$clust.res$samID

ICB<-tpm[c("PDCD1","CTLA4"),finalsam]
ICB<-log(ICB+1)

output<-tpm[,finalsam]
group<-as.data.frame(cmoic.KIRC$clust.res)

#write.table(output,"expr225.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#write.table(group,"group225.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

ICB<-cbind(t(ICB),group)

p <- ggboxplot(ICB, x = "clust", y = "CTLA4",
               color = "clust", palette = c("#2EC4B6","#E71D36","#FF9F1C"),
               add = "jitter")
# Change method
p + stat_compare_means(method = "anova")

ggsave(file=file.path(fig2.path,"Figure 2C. CTLA4.pdf"), height = 4, width = 4)

p <- ggboxplot(ICB, x = "clust", y = "PDCD1",
               color = "clust", palette = c("#2EC4B6","#E71D36","#FF9F1C"),
               add = "jitter")
# Change method
p + stat_compare_means(method = "anova")

ggsave(file=file.path(fig2.path,"Figure 2D. PDCD1.pdf"), height = 4, width = 4)

####Figure 2E######
#------------------------------------------------------------------------------------------#
#For Figure 2E, we predicted the potential response to immune therapy via SubMap Analysis,
#trough the online GenePattern, the subsequent code just show the process of how to draw a heatmap
#The results of SubMap analysis saved in the file named "SubMap_SubMapResult.txt"
#-----------------------------------------------------------------------------------------#

# Read content from the file
lines <- readLines(file.path(data.path,"SubMap_SubMapResult.txt"), warn = FALSE)

# Concatenate the character vector into a single string
txt_content <- paste(lines, collapse = "\n")

# Extract content from $SA.matrix and $nominal.p.matrix.Fisher, then get all the numbers
get_numbers_from_section <- function(section_name) {
  section_txt <- sub(paste0(".*\\$", section_name, "\\n"), paste0("$", section_name, "\n"), txt_content)
  section_txt <- sub("\\n\\$.*", "", section_txt)
  regmatches(section_txt, gregexpr("\\b\\d+\\.?\\d*\\b", section_txt))[[1]]
}

sa_numbers <- as.numeric(get_numbers_from_section("SA\\.matrix"))
fisher_numbers <- as.numeric(get_numbers_from_section("nominal\\.p\\.matrix\\.Fisher"))

# Create matrices
sa_matrix <- matrix(sa_numbers, ncol=4, byrow=TRUE)
fisher_matrix <- matrix(fisher_numbers, ncol=4, byrow=TRUE)

# Set row and column names
colnames(sa_matrix) <- colnames(fisher_matrix) <- c("CTAL4-noR", "CTLA4-R", "PD1-noR", "PD1-R")
rownames(sa_matrix) <- c("MoS1-b","MoS2-b","MoS3-b")
rownames(fisher_matrix) <- c("MoS1","Mos2","MoS3")

# Combine both matrices
final_matrix <- rbind(fisher_matrix, sa_matrix)

# Assign the result to the tmp variable
tmp <- final_matrix

library(pheatmap)
heatmap.YlGnPe <- c("#0D2735","#21526C","#25769A","#539CB5","#9BC6D4")
cherry    <- "#383D49"
  lightgrey <- "#dcddde"
    
pheatmap(tmp, cellwidth = 30, cellheight = 30,
           cluster_rows = F,cluster_cols = F,
           color = heatmap.YlGnPe[5:1],
           gaps_row = 3,
           annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Nominal p value","Bonferroni corrected","Bonferroni corrected","Bonferroni corrected"),row.names = rownames(tmp)),
           annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
           filename = file.path(fig2.path, "Figure 2E.heatmap_submap.pdf"))

#####Figure 2F#######
#----------------CheckMate-prepare validate cohort
load("./InputData/TCGA.marker.rda")

CheckMate.expr<-read.table(file.path(data.path,"CheckMate.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)
CheckMate.clin<-read.table(file.path(data.path,"CheckMate.clin.txt"),header=T,sep="\t",row.names=1,check.names=F)
CheckMate.expr2<-CheckMate.expr[,rownames(CheckMate.clin)]

CheckMate.ntp.pred <- runNTP(expr      = CheckMate.expr2,
                             templates = marker.up$templates, # the template has been already prepared in runMarker()
                             scaleFlag = TRUE,
                             centerFlag= TRUE,
                             nPerm     = 1000,
                             distance  = "cosine",
                             seed      = 123456,
                             verbose   = TRUE,
                             doPlot    = TRUE, # to generate heatmap
                             fig.name  = "Figure 2F.NTP HEATMAP FOR CheckMate",
                             fig.path  = fig2.path) 
# Intersect and subset
merge<-intersect(rownames(CheckMate.clin),CheckMate.ntp.pred$clust.res$samID)
CheckMate.clin<-CheckMate.clin[merge,]
IMclust <- CheckMate.ntp.pred$clust.res[merge, ]
IMout<-cbind(CheckMate.clin,IMclust)

# Chi-square test
test.data<-print(table(IMout$clust,IMout$Benefit))
ax<-chisq.test(test.data)$p.value

# Ordering and subsetting
test.data2 <- data.frame(test.data)
test.data2 <- test.data2[order(test.data2[,1]),] 

subsets <- lapply(1:3, function(i) {
  subset <- test.data2[(2*i-1):(2*i), ]
  subset$pct <- subset$Freq / sum(subset$Freq)
  return(subset)
})
# Combine subsets
test.data3<-rbind.data.frame(subset1,subset2,subset3)
test.data3

# Generate the plot
p<- ggplot(test.data3,aes(x=Var1, y=pct, fill=Var2)) + 
  geom_bar(stat="identity",position="stack", colour="black")+
  guides(fill=guide_legend(reverse=TRUE)) +
  geom_text(aes(label =scales::percent (pct)), position = position_stack(vjust = .5), color="black", size=5)+
  labs(x="", y="Percentage", fill="",size=15) + 
  theme(plot.title = element_text(size=25, margin=margin(t=20, b=30)))+
  theme(axis.text.x = element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(legend.text=element_text(size=15,colour='black'),legend.position = 'right')+
  scale_fill_manual(values = wes_palette("BottleRocket2", n = 2))

pdf(file.path(fig2.path,"Figure 2F. CheckMate Immunotherapy response.pdf"), width=5,height=4,onefile = FALSE)
p
dev.off()

#####Figure 2G#######
#----------------Miao-prepare validate cohort
library("wesanderson")
Miao.expr<-read.table(file.path(data.path,"Miao.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)
Miao.expr<-log2(Miao.expr+1)
Miao.clin<-read.table(file.path(data.path,"Miao.clin.txt"),header=T,sep="\t",row.names=1,check.names=F)

Miao.expr<-Miao.expr[,rownames(Miao.clin)]

Miao.ntp.pred <- runNTP(expr      = Miao.expr,
                        templates = marker.up$templates, # the template has been already prepared in runMarker()
                        scale     = TRUE, # scale input data (by default)
                        center    = TRUE, # center input data (by default)
                        doPlot    = TRUE, # to generate heatmap
                        fig.name  = "Figure 2F. NTP HEATMAP FOR Miao",
                        fig.path  = fig2.path) 


# Intersect and subset
merge<-intersect(rownames(Miao.clin),Miao.ntp.pred$clust.res$samID)
Miao.clin<-Miao.clin[merge,]
Miaoclust <- Miao.ntp.pred$clust.res[merge, ]
Miaoclust<-cbind(Miao.clin,Miaoclust)

# Chi-square test
test.data<-print(table(Miaoout$clust,Miaoout$response_category))
ax<-chisq.test(test.data)$p.value

# Ordering and subsetting
test.data2 <- data.frame(test.data)
test.data2 <- test.data2[order(test.data2[,1]),] 

subsets <- lapply(1:3, function(i) {
  subset <- test.data2[(2*i-1):(2*i), ]
  subset$pct <- subset$Freq / sum(subset$Freq)
  return(subset)
})
# Combine subsets
test.data3 <- do.call(rbind, subsets)

# Generate the plot
p<- ggplot(test.data3,aes(x=Var1, y=pct, fill=Var2)) + 
  geom_bar(stat="identity",position="stack", colour="black")+
  guides(fill=guide_legend(reverse=TRUE)) +
  geom_text(aes(label =scales::percent (pct)), position = position_stack(vjust = .5), color="black", size=5)+
  labs(x="", y="Percentage", fill="",size=15) + 
  theme(plot.title = element_text(size=25, margin=margin(t=20, b=30)))+
  theme(axis.text.x = element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(legend.text=element_text(size=15,colour='black'),legend.position = 'right')+
  scale_fill_manual(values = wes_palette("BottleRocket2", n = 2))

pdf(file.path(fig2.path,"Figure 2G. Miao Immunotherapy response.pdf"), width=5,height=4,onefile = FALSE)
p
dev.off()

#--------------------------#
#----------Figure 3--------#
#--------------------------#



#--------------------------#
#----------Figure 5--------#
#--------------------------#
####Figure 5A####
drug.KIRC <- compDrugsen(moic.res    = cmoic.KIRC,
                         norm.expr   = log2(tpm+1)[,cmoic.KIRC$clust.res$samID], # double guarantee sample order
                         drugs       = c("Sunitinib"), # a vector of names of drug in GDSC
                         tissueType  = "all", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "Figure 5A.BOXVIOLIN OF ESTIMATED IC50",
                         fig.path    = fig5.path) 

####Figure 5B####

EMTAB3267.expr<-read.table(file.path(data.path,"EMTAB3267.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)
EMTAB3267.expr<-t(EMTAB3267.expr)
EMTAB3267.clin<-read.table(file.path(data.path,"EMTAB3267.clin.txt"),header=T,sep="\t",row.names=1,check.names=F)
EMTAB3267.expr<-EMTAB3267.expr[,rownames(EMTAB3267.clin)]
EMTAB3267.clin$`sunitinib-response`<-ifelse(EMTAB3267.clin$`sunitinib-response`=="PR","DPR",EMTAB3267.clin$`sunitinib-response`)


EMTAB3267.ntp.pred <- runNTP(expr      = EMTAB3267.expr,
                             templates = marker.up$templates, # the template has been already prepared in runMarker()
                             scale     = TRUE, # scale input data (by default)
                             center    = TRUE, # center input data (by default)
                             doPlot    = TRUE, # to generate heatmap
                             fig.name  = "Figure 5B.NTP HEATMAP FOR EMTAB3267",
                             fig.path  = fig5.path) 

surv.EMTAB3267 <- compSurv(moic.res         = EMTAB3267.ntp.pred,
                           surv.info        = EMTAB3267.clin,
                           convt.time       = "m", # switch to year
                           surv.median.line = "hv", # switch to both
                           fig.name         = "Figure 5B. KAPLAN-MEIER CURVE OF NTP FOR EMTAB3267",
                           fig.path         = fig5.path) 


###Figure 5C#####
library(ggplot2)
df<-read.delim(file.path(data.path,'MoS drug prediction.txt'), sep = '\t', stringsAsFactors = FALSE,row.names = 1)
green <- "#2EC4B6";cyan <- "#E71D36";blue <- "#FF9F1C"
  
p.val <- kruskal.test(Axitinib ~ MoS,
                      data = df)
p.lab<-paste0("P = ",formatC(p.val$p.value, format = "e", digits = 3),sep = "")      
p_top <- ggplot(df, aes(x = Axitinib, color = MoS, fill = MoS)) +
      geom_density() +
      scale_color_manual(values = c(alpha(green,0.7),alpha(cyan,0.7),alpha(blue,0.7))) + # 设置透明色
      scale_fill_manual(values = c(alpha(green,0.7),alpha(cyan,0.7),alpha(blue,0.7))) +
      theme_classic() + 
      xlab(paste0("Estimated IC50 of ", unique(df$Drug))) + ylab(NULL) + 
      theme(legend.position = "none", 
            legend.title = element_blank(),
            axis.text.x = element_text(size = 12,color = "black"),
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) +
      geom_rug()
    p_top
    
    p_bot <- ggplot(df, aes(MoS, Axitinib, fill = MoS)) + 
      geom_boxplot(aes(col = MoS)) + 
      scale_fill_manual(values = c(green, cyan, blue)) + 
      scale_color_manual(values = c(green, cyan, blue)) + 
      xlab(NULL) + ylab("Estimated IC50") + 
      theme_void() +
      theme(legend.position = "right",
            legend.title = element_blank(),
            axis.text.x = element_blank(), # 
            axis.text.y = element_text(size = 11,color = "black"),
            panel.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()) + 
      annotate(geom="text",
               x = 1.5,
               hjust = 1,
               y = 5,
               size = 4, angle = 270, fontface = "bold",
               label = p.lab) +
      coord_flip() 
    
    dat <- ggplot_build(p_bot)$data[[1]]
    p_bot <- p_bot + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
    p_bot
    
 p <- p_top %>% insert_bottom(p_bot, height = 0.4)
pdf(file = file.path(fig5.path,"Figure 5C. Axitinib boxdensity.pdf"), width = 8,height = 6)
p
invisible(dev.off())
#---------------------------

p.val <- kruskal.test(GDC0941 ~ MoS,
                      data = df)
p.lab<-paste0("P = ",formatC(p.val$p.value, format = "e", digits = 3),sep = "")     

p_top <- ggplot(df, aes(x = GDC0941, color = MoS, fill = MoS)) +
  geom_density() +
  scale_color_manual(values = c(alpha(green,0.7),alpha(cyan,0.7),alpha(blue,0.7))) + # 设置透明色
  scale_fill_manual(values = c(alpha(green,0.7),alpha(cyan,0.7),alpha(blue,0.7))) +
  theme_classic() + 
  xlab(paste0("Estimated IC50 of ", unique(df$Drug))) + ylab(NULL) + 
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12,color = "black"),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_rug()
p_top

p_bot <- ggplot(df, aes(MoS, GDC0941, fill = MoS)) + 
  geom_boxplot(aes(col = MoS)) + 
  scale_fill_manual(values = c(green, cyan, blue)) + 
  scale_color_manual(values = c(green, cyan, blue)) + 
  xlab(NULL) + ylab("Estimated IC50") + 
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text.x = element_blank(), # 
        axis.text.y = element_text(size = 11,color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  annotate(geom="text",
           x = 1.5,
           hjust = 1,
           y = 7,
           size = 4, angle = 270, fontface = "bold",
           label = p.lab) +
  coord_flip() 

dat <- ggplot_build(p_bot)$data[[1]]
p_bot <- p_bot + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p_bot

p <- p_top %>% insert_bottom(p_bot, height = 0.4)
pdf(file = file.path(fig5.path,"Figure 5C. GDC0941 boxdensity.pdf"), width = 8,height = 6)
p
invisible(dev.off())
#----------------
p.val <- kruskal.test(Dimethyloxalylglycine ~ MoS,
                      data = df)
p.lab<-paste0("P = ",formatC(p.val$p.value, format = "e", digits = 3),sep = "")     

p_top <- ggplot(df, aes(x = Dimethyloxalylglycine, color = MoS, fill = MoS)) +
  geom_density() +
  scale_color_manual(values = c(alpha(green,0.7),alpha(cyan,0.7),alpha(blue,0.7))) + # 设置透明色
  scale_fill_manual(values = c(alpha(green,0.7),alpha(cyan,0.7),alpha(blue,0.7))) +
  theme_classic() + 
  xlab(paste0("Estimated IC50 of ", unique(df$Drug))) + ylab(NULL) + 
  theme(legend.position = "none", 
        legend.title = element_blank(),
        axis.text.x = element_text(size = 12,color = "black"),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_rug()
p_top

p_bot <- ggplot(df, aes(MoS, Dimethyloxalylglycine, fill = MoS)) + 
  geom_boxplot(aes(col = MoS)) + 
  scale_fill_manual(values = c(green, cyan, blue)) + 
  scale_color_manual(values = c(green, cyan, blue)) + 
  xlab(NULL) + ylab("Estimated IC50") + 
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_blank(),
        axis.text.x = element_blank(), # 
        axis.text.y = element_text(size = 11,color = "black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  annotate(geom="text",
           x = 1.5,
           hjust = 1,
           y = 15,
           size = 4, angle = 270, fontface = "bold",
           label = p.lab) +
  coord_flip() 

dat <- ggplot_build(p_bot)$data[[1]]
p_bot <- p_bot + geom_segment(data=dat, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", inherit.aes = F)
p_bot

p <- p_top %>% insert_bottom(p_bot, height = 0.4)
pdf(file = file.path(fig5.path,"Figure 5C. Dimethyloxalylglycine boxdensity.pdf"), width = 8,height = 6)
p
invisible(dev.off())

###Figure 5D#####    
drugs <- read.delim(file.path(data.path,'SETD2 export.txt'), sep = '\t', stringsAsFactors = FALSE)
head(drugs)

drugs[which(drugs$P.value < 0.05 & drugs$Effect.size <= 0),'sig'] <- 'Sensitive'
#drugs[which(drugs$P.value < 0.05 & drugs$Effect.size >= 0),'sig'] <- 'Up'
drugs[which(drugs$P.value >= 0.05 | abs(drugs$Effect.size) < 0.3),'sig'] <- 'None'

highlighted_drugs <- c("AZD8186", "AZD5363", "Alpelisib")

p <- ggplot(drugs, aes(x = Effect.size, y = -log10(P.value), color = sig)) +
  geom_point(alpha = 1.2, size = 3) +
  scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c('Sensitive', 'Down', 'None')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
  theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
  geom_vline(xintercept = c(0), color = 'gray', size = 0.3) +
  geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.3) +
  geom_label(data = subset(drugs, Drug %in% highlighted_drugs), aes(label = Drug), size = 3, nudge_y = 0.2, show.legend = FALSE) + # 新添加的代码来标记三个药物
  xlim(-1, 1) + ylim(0, 3.1) +
  labs(x = '\nIC50 Effect', y = 'log10(p-value)\n', color = '', title = 'Drugs for SETD2 mut\n')

ggsave(file.path(fig5.path,'Figure 5D. Drugs for SETD2.pdf'), p, width = 4, height = 4)

####Figure 5E#####
load("./InputData/PI3KAKT inh predi.rda")

#适合展示两种药物
p1 <- plot_grid(plotp[[1]],plotp[[2]],nrow = 1) # title可以AI下拉到合适位置，就如例文所示

ggsave(file.path(fig5.path,"Figure 5E. boxplot of predicted IC50.pdf"), width = 8, height = 4)


#--------------------------#
#----------Figure 6--------#
#--------------------------#

#####GSE22541
GSE22541.expr<-read.table(file.path(data.path,"GSE22541.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)
GSE22541.clin<-read.table(file.path(data.path,"GSE22541.clin.txt"),header=T,sep="\t",row.names=1,check.names=F)
GSE22541.expr<-GSE22541.expr[,rownames(GSE22541.clin)]
GSE22541.expr<-log2(GSE22541.expr+1)

GSE22541.ntp.pred <- runNTP(expr      = GSE22541.expr,
                            templates = marker.up$templates, # the template has been already prepared in runMarker()
                            scale     = TRUE, # scale input data (by default)
                            center    = TRUE, # center input data (by default)
                            doPlot    = TRUE, # to generate heatmap
                            fig.name  = "Figure 6A.NTP HEATMAP FOR GSE22541",
                            fig.path    = fig6.path) 

surv.GSE22541 <- compSurv(moic.res         = GSE22541.ntp.pred,
                          surv.info        = GSE22541.clin,
                          convt.time       = "m", # switch to year
                          surv.median.line = "hv", # switch to both
                          fig.name         = "Figure 6A.KAPLAN-MEIER CURVE OF NTP FOR GSE22541",
                          fig.path    = fig6.path) 

####GSE40435
GSE40435.expr<-read.table(file.path(data.path,"GSE40435.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)
GSE40435.clin<-read.table(file.path(data.path,"GSE40435.clin.txt"),header=T,sep="\t",row.names=1,check.names=F)
GSE40435.expr2<-GSE40435.expr[,rownames(GSE40435.clin)]

GSE40435.ntp.pred <- runNTP(expr      = GSE40435.expr2,
                            templates = marker.up$templates, # the template has been already prepared in runMarker()
                            scaleFlag = TRUE,
                            centerFlag = TRUE,
                            nPerm = 1000,
                            distance = "cosine",
                            seed = 123456,
                            verbose = TRUE,
                            doPlot    = TRUE, # to generate heatmap
                            fig.name  = "Figure 6B.NTP HEATMAP FOR GSE40435",
                            fig.path    = fig6.path) 

merge<-intersect(rownames(GSE40435.clin),GSE40435.ntp.pred$clust.res$samID)

GSE40435.clin<-GSE40435.clin[merge,]
IMclust<-GSE40435.ntp.pred$clust.res
IMclust<-IMclust[merge,]

IMout<-cbind(GSE40435.clin,IMclust)

test.data<-print(table(IMout$clust,IMout$`1:GRADE`))

test.data2<-data.frame(test.data)
test.data2<-test.data2[order(test.data2[,1]),] 
test.data2

subset1 = test.data2[1:4,]
subset1$pct = subset1$Freq/sum(subset1$Freq)
subset1

subset2 = test.data2[5:8,]
subset2$pct = subset2$Freq/sum(subset2$Freq)
subset2

subset3 = test.data2[9:12,]
subset3$pct = subset3$Freq/sum(subset3$Freq)
subset3

test.data3<-rbind.data.frame(subset1,subset2,subset3)
test.data3

p<- ggplot(test.data3,aes(x=Var1, y=pct, fill=Var2)) + 
  geom_bar(stat="identity",position="stack", colour="black")+
  guides(fill=guide_legend(reverse=TRUE)) +
  geom_text(aes(label =scales::percent (pct)), position = position_stack(vjust = .5), color="black", size=5)+
  labs(x="", y="Percentage", fill="",size=15) + 
  theme(plot.title = element_text(size=25, margin=margin(t=20, b=30)))+
  theme(axis.text.x = element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(legend.text=element_text(size=15,colour='black'),legend.position = 'right')+
  scale_fill_brewer(palette="Blues",direction = 1)

pdf(file.path(fig6.path,"Figure 6B.GSE40435 Group and stage.pdf"), width=5,height=4,onefile = FALSE)
p
dev.off()

#####GSE53757
#----------------GSE53757-prepare validate cohort

GSE53757.expr<-read.table(file.path(data.path,"GSE53757.expr.txt"),header=T,sep="\t",row.names=1,check.names=F)
GSE53757.clin<-read.table(file.path(data.path,"GSE53757.clin.txt"),header=T,sep="\t",row.names=1,check.names=F)
GSE53757.expr2<-GSE53757.expr[,rownames(GSE53757.clin)]

GSE53757.ntp.pred <- runNTP(expr      = GSE53757.expr2,
                            templates = marker.up$templates, # the template has been already prepared in runMarker()
                            scaleFlag = TRUE,
                            centerFlag = TRUE,
                            nPerm = 1000,
                            distance = "cosine",
                            seed = 123456,
                            verbose = TRUE,
                            doPlot    = TRUE, # to generate heatmap
                            fig.name  = "Figure 6C.NTP HEATMAP FOR GSE53757",
                            fig.path    = fig6.path) 


merge<-intersect(rownames(GSE53757.clin),GSE53757.ntp.pred$clust.res$samID)

GSE53757.clin<-GSE53757.clin[merge,]
IMclust<-GSE53757.ntp.pred$clust.res
IMclust<-IMclust[merge,]
IMout<-cbind(GSE53757.clin,IMclust)

test.data<-print(table(IMout$clust,IMout$STAGE))
test.data2<-data.frame(test.data)
test.data2<-test.data2[order(test.data2[,1]),] 
test.data2

subset1 = test.data2[1:4,]
subset1$pct = subset1$Freq/sum(subset1$Freq)
subset1

subset2 = test.data2[5:8,]
subset2$pct = subset2$Freq/sum(subset2$Freq)
subset2

subset3 = test.data2[9:12,]
subset3$pct = subset3$Freq/sum(subset3$Freq)
subset3

test.data3<-rbind.data.frame(subset1,subset2,subset3)
test.data3

p<- ggplot(test.data3,aes(x=Var1, y=pct, fill=Var2)) + 
  geom_bar(stat="identity",position="stack", colour="black")+
  guides(fill=guide_legend(reverse=TRUE)) +
  geom_text(aes(label =scales::percent (pct)), position = position_stack(vjust = .5), color="black", size=5)+
  labs(x="", y="Percentage", fill="",size=15) + 
  theme(plot.title = element_text(size=25, margin=margin(t=20, b=30)))+
  theme(axis.text.x = element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(axis.text.y = element_text(size = 12, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 15, color = "black", face = "plain", vjust = 0.5, hjust = 0.5))+
  theme(legend.text=element_text(size=15,colour='black'),legend.position = 'right')+
  scale_fill_brewer(palette="Blues",direction = 1)

pdf(file.path(fig6.path,"Figure 6C. GSE53757 Group and stage.pdf"), width=5,height=4,onefile = FALSE)
p
dev.off()

#--------------------------#
#----------Figure S3--------#
#--------------------------#
#------------------GSE22541
# run DEA with limma
runDEA(dea.method = "limma",
       expr       = GSE22541.expr, # normalized expression data
       moic.res   = GSE22541.ntp.pred,
       prefix     = "GSE22541",
       res.path   = figS.path )

MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)
gsea.up <- runGSEA(moic.res     = GSE22541.ntp.pred,
                   dea.method   = "limma", # name of DEA method
                   prefix       = "GSE22541", # MUST be the same of argument in runDEA()
                   dat.path      = figS.path, # path of DEA files
                   res.path      = figS.path, # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = GSE22541.expr, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.15, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "Figure S3. GSE22541 Go UPREGULATED PATHWAY HEATMAP",
                   fig.path      = figS.path)

GSET.FILE <- 
  system.file("extdata", "up pannal signals.gmt", package = "MOVICS", mustWork = TRUE)

gsva.res <- 
  runGSVA(moic.res      = GSE22541.ntp.pred,
          norm.expr     = GSE22541.expr,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          #annCol        = annCol,
          #annColors     = annColors,
          color         =c("#3F89C9","white", "#D31D24"),
          fig.path      = figS.path,
          fig.name      = "Figure S3.GSE22541 immune GENE SETS OF INTEREST HEATMAP",
          height        = 8,
          width         = 10)

##GSE40435
runDEA(dea.method = "limma",
       expr       = GSE40435.expr, # normalized expression data
       moic.res   = GSE40435.ntp.pred,
       prefix     = "GSE40435",
       res.path   = figS.path )

MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)
gsea.up <- runGSEA(moic.res     = GSE40435.ntp.pred,
                   dea.method   = "limma", # name of DEA method
                   prefix       = "GSE40435", # MUST be the same of argument in runDEA()
                   dat.path      = figS.path, # path of DEA files
                   res.path      = figS.path, # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = GSE40435.expr, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.15, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "Figure S3. GSE40435 Go UPREGULATED PATHWAY HEATMAP",
                   fig.path      = figS.path)

GSET.FILE <- 
  system.file("extdata", "up pannal signals.gmt", package = "MOVICS", mustWork = TRUE)

gsva.res <- 
  runGSVA(moic.res      = GSE40435.ntp.pred,
          norm.expr     = GSE40435.expr,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          #annCol        = annCol,
          #annColors     = annColors,
          color         =c("#3F89C9","white", "#D31D24"),
          fig.path      = figS.path,
          fig.name      = "Figure S3.GSE40435 immune GENE SETS OF INTEREST HEATMAP",
          height        = 8,
          width         = 10)
##GSE53757

runDEA(dea.method = "limma",
       expr       = GSE53757.expr, # normalized expression data
       moic.res   = GSE53757.ntp.pred,
       prefix     = "GSE53757",
       res.path   = figS.path )

MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)
gsea.up <- runGSEA(moic.res     = GSE53757.ntp.pred,
                   dea.method   = "limma", # name of DEA method
                   prefix       = "GSE53757", # MUST be the same of argument in runDEA()
                   dat.path      = figS.path, # path of DEA files
                   res.path      = figS.path, # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = GSE53757.expr, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.15, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "Figure S3. GSE53757 Go UPREGULATED PATHWAY HEATMAP",
                   fig.path      = figS.path)

GSET.FILE <- 
  system.file("extdata", "up pannal signals.gmt", package = "MOVICS", mustWork = TRUE)

gsva.res <- 
  runGSVA(moic.res      = GSE53757.ntp.pred,
          norm.expr     = GSE53757.expr,
          gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
          gsva.method   = "gsva", # method to calculate single sample enrichment score
          #annCol        = annCol,
          #annColors     = annColors,
          color         =c("#3F89C9","white", "#D31D24"),
          fig.path      = figS.path,
          fig.name      = "Figure S3.GSE53757 immune GENE SETS OF INTEREST HEATMAP",
          height        = 8,
          width         = 10)
