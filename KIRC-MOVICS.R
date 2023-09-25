setwd("E:\\2021ÂÛÎÄ\\KIRCmulti")

load("mo.data.RData")

library("MOVICS")

optk.KIRC <- getClustNum(data        = mo.data,
                         is.binary   = c(F,F,F,F,T), # note: the 4th data is somatic mutation which is a binary matrix
                         try.N.clust = 2:8, # try cluster number from 2 to 8
                         fig.name    = "CLUSTER NUMBER OF TCGA-KIRC")




iClusterBayes.res <- getiClusterBayes(data        = mo.data,
                                      N.clust     = 3,
                                      type        = c("gaussian","gaussian","gaussian","gaussian","binomial"),
                                      n.burnin    = 1800,
                                      n.draw      = 1200,
                                      prior.gamma = c(0.5, 0.5, 0.5, 0.5,0.3),
                                      sdev        = 0.05,
                                      thin        = 3)

moic.res.list <- getMOIC(data        = mo.data,
                         methodslist = list("SNF", "PINSPlus", "NEMO", "KIRC", "LRAcluster", "ConsensusClustering", "IntNMF", "CIMLR", "MoCluster"),
                         N.clust     = 3,
                         type        = c("gaussian", "gaussian", "gaussian", "gaussian","binomial"))



moic.res.list <- append(moic.res.list, 
                        list("iClusterBayes" = iClusterBayes.res))

# save moic.res.list to local path
save(moic.res.list, file = "moic.res.list.rda")


cmoic.KIRC <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")

getSilhouette(sil      = cmoic.KIRC$sil, # a sil object returned by getConsensusMOIC()
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)

indata <- mo.data
indata$meth <- log2(indata$meth / (1 - indata$meth))

# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,T,F)) # no scale for mutation

feat   <- iClusterBayes.res$feat.res
feat1  <- feat[which(feat$dataset == "mRNA"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "lncRNA"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "cna"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "meth"),][1:10,"feature"]
feat5  <- feat[which(feat$dataset == "mut"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4,feat5)

test<-intersect(colnames(lncrna.tpm),rownames(surv.info))

mRNA.col   <- c("#31B29E", "black","#EB6C5A")
lncRNA.col <- c("#6699CC", "white"  , "#FF3C38")
cna.col  <-  c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
meth.col   <- c("#0099CC", "white","#CC0033")
mut.col    <- c("grey90" , "black")
col.list   <- list(mRNA.col, lncRNA.col, cna.col,meth.col, mut.col)

# comprehensive heatmap (may take a while)
surv.info=read.table("KIRCsurv.txt",header=T,sep="\t",row.names=1,check.names=F)
surv.info <- rt222[,c("futime","fustat","Age","Grade","Stage","New_tumor")]

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
             fig.name      = "2COMPREHENSIVE HEATMAP OF CONSENSUSMOIC")





t1<-moic.res.list$LRAcluster$clust.res$clust
t2<-moic.res.list$SNF$clust.res$clust
t3<-moic.res.list$PINSPlus$clust.res$clust
t4<-moic.res.list$NEMO$clust.res$clust
t5<-moic.res.list$COCA$clust.res$clust
t6<-moic.res.list$ConsensusClustering$clust.res$clust
t7<-moic.res.list$IntNMF$clust.res$clust
t8<-moic.res.list$CIMLR$clust.res$clust
t9<-moic.res.list$MoCluster$clust.res$clust
t10<-moic.res.list$iClusterBayes$clust.res$clust
ttotal<-rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10)
rownames(ttotal)<-c("LRA","SNF","PINSPlus","NEMO","COCA","ConsensusClustering","IntNMF","CIMLR","MoCluster","iClusterBayes")
colnames(ttotal)<-moic.res.list$LRAcluster$clust.res$samID

write.table(ttotal,file="10methodcluster.txt",sep="\t",row.names=T,quote=F)

library(dplyr)
library(pheatmap)  

rt3<-lncrna.tpm[1:10,]

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
        outdata=t(scale(t(rt3), center=centerFlag, scale=scaleFlag))
        if (!is.null(halfwidth)) {
                outdata[outdata>halfwidth]=halfwidth
                outdata[outdata<(-halfwidth)]= -halfwidth
        }
        return(outdata)
}

rt3 <- standarize.fun(rt3,halfwidth =0.8)
rt3[1:10,1:5]
rt3<-rt3[,row.names(rt4)]


rt5<-rbind(ttotal,rt3[,colnames(ttotal)])

rt5<-t(rt5)
rt5<-as.data.frame(rt5)

rt5<-rt5[order(rt5$LRA,decreasing = F),] 

rt6<-t(rt5[,11:20])

rt4<-rt5[,1:10]


        pdf("ClustOutPut.pdf",height=3,width=10)
p_heatmap<-pheatmap(rt6, annotation=rt4, 
                    color = colorRampPalette(c("#0099CC", "white","#CC0033"))(50),
                    cluster_cols =F,
                    fontsize=8,
                    fontsize_row=8,
                    scale="row",
                    show_colnames=F,
                    fontsize_col=3)
dev.off()

p_heatmap<-pheatmap(rt3, 
                    #annotation=rt4, 
                    #color = colorRampPalette(c("#0099CC", "white","#CC0033"))(50),
                    cluster_cols =F,
                    fontsize=8,
                    fontsize_row=8,
                    scale="row",
                    show_colnames=F,
                    fontsize_col=3)

p_heatmap



surv.info=read.table("KIRCsurv.txt",header=T,sep="\t",row.names=1,check.names=F)

surv.info$futime<-surv.info$PFI.time
surv.info$fustat<-surv.info$PFI


surv.KIRC <- compSurv(moic.res         = cmoic.KIRC,
                      surv.info        = surv.info,
                      surv.cut         = 150,
                      convt.time       = "m", # convert day unit to month
                      surv.median.line = "h",
                      fig.name         = "PFI KAPLAN-MEIER CURVE OF CONSENSUSMOIC")

print(surv.KIRC)


clin.KIRC <- compClinvar(moic.res      = cmoic.KIRC,
                         var2comp      = surv.info, # data.frame needs to summarize (must has row names of samples)
                         strata        = "Subtype", # stratifying variable (e.g., Subtype in this example)
                         factorVars    = c("PAM50","pstage","fustat"), # features that are considered categorical variables
                         nonnormalVars = "futime", # feature(s) that are considered using nonparametric test
                         exactVars     = "pstage", # feature(s) that are considered using exact test
                         doWord        = TRUE, # generate .docx file in local path
                         tab.name      = "SUMMARIZATION OF CLINICAL FEATURES")

print(clin.KIRC$compTab)



mut.KIRC <- compMut(moic.res     = cmoic.KIRC,
                    mut.matrix   = mut2, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
                    p.cutoff = 0.2,
                    #p.adj.cutoff = 0.2, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                    innerclust   = TRUE, # perform clustering within each subtype
                    annCol       = annCol, # same annotation for heatmap
                    annColors    = annColors, # same annotation color for heatmap
                    width        = 11, 
                    height       = 4,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")

print(mut.KIRC)

head(maf)


tmb.KIRC <- compTMB(moic.res     = cmoic.KIRC,
                    maf          = maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")

head(tmb.KIRC$TMB)


segment$value<-segment$seg.mean
segment$chrom<-segment$chromosome
head(segment)

fga.KIRC <- compFGA(moic.res     = cmoic.KIRC,
                    segment      = segment,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "BARPLOT OF FGA",
                    width = 7,
                    height = 4)
head(fga.KIRC$summary)

drug.KIRC <- compDrugsen(moic.res    = cmoic.KIRC,
                         norm.expr   = fpkm[,cmoic.KIRC$clust.res$samID], # double guarantee sample order
                         drugs       = c("AZ628","AKT inhibitor VIII"), # a vector of names of drug in GDSC
                         tissueType  = "all", # choose specific tissue type to construct ridge regression model
                         test.method = "nonparametric", # statistical testing method
                         prefix      = "BOXVIOLIN OF ESTIMATED IC50") 
head(drug.KIRC$Sorafenib)


surv.info$pstage <- factor(surv.info$pstage, levels = c("TX","T1","T2","T3","T4"))

# agreement compMiaoon (support up to 6 classifications include current subtype)
agree.KIRC <- compAgree(moic.res  = cmoic.KIRC,
                        subt2comp = surv.info[,"Grade",drop=FALSE],
                        doPlot    = TRUE,
                        box.width = 0.2,
                        fig.name  = "AGREEMENT OF CONSENSUSMOIC WITH Grade")
#> --all samples matched.


print(agree.KIRC)


#------------------DEGs  different pathways
# run DEA with limma
runDEA(dea.method = "limma",
       expr       = fpkm, # normalized expression data
       moic.res   = cmoic.KIRC,
       prefix     = "TCGA-KIRC")
#> --all samples matched.
#> --you choose limma and please make sure a microarray profile or a normalized expression data [FPKM or TPM without log2 transformation is recommended] was provided.
#> --log2 transformation done for expression data.
#> limma of CS1_vs_Others done...
#> limma of CS2_vs_Others done...
#> limma of CS3_vs_Others done...
#> limma of CS4_vs_Others done...
#> limma of CS5_vs_Others done...

marker.up <- runMarker(moic.res      = cmoic.KIRC,
                       dea.method    = "limma", # name of DEA method
                       prefix        = "TCGA-KIRC", # MUST be the same of argument in runDEA()
                       dat.path      = getwd(), # path of DEA files
                       res.path      = getwd(), # path to save marker files
                       p.cutoff      = 0.05, # p cutoff to identify significant DEGs
                       p.adj.cutoff  = 0.05, # padj cutoff to identify significant DEGs
                       dirct         = "up", # direction of dysregulation in expression
                       n.marker      = 100, # number of biomarkers for each subtype
                       doplot        = TRUE, # generate diagonal heatmap
                       norm.expr     = fpkm, # use normalized expression as heatmap input
                       annCol        = annCol, # sample annotation in heatmap
                       annColors     = annColors, # colors for sample annotation
                       show_rownames = FALSE, # show no rownames (biomarker name)
                       color = c("#004E6B", "white", "#FF0707"),
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")

tcga.ntp.pred <- runNTP(expr      = mrna.tpm[,colnames(mut2)],
                        templates = marker.up$templates, # the template has been already prepared in runMarker()
                        scale     = TRUE, # scale input data (by default)
                        center    = TRUE, # center input data (by default)
                        doPlot    = TRUE, # to generate heatmap
                        fig.name  = "NTP HEATMAP FOR KIRC") 
runKappa(subt1     = cmoic.KIRC$clust.res$clust,
         subt2     = tcga.ntp.pred$clust.res$clust,
         subt1.lab = "CMOIC",
         subt2.lab = "NTP",
         fig.name  = "CONSISTENCY HEATMAP FOR TCGA between CMOIC and NTP")


MSIGDB.FILE <- system.file("extdata", "reactome.gmt", package = "MOVICS", mustWork = TRUE)
gsea.up <- runGSEA(moic.res     = cmoic.KIRC,
                   dea.method   = "limma", # name of DEA method
                   prefix       = "TCGA-KIRC", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = fpkm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.15, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "reactome UPREGULATED PATHWAY HEATMAP")
head(round(gsea.up$grouped.es,3))



MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)
gsea.up <- runGSEA(moic.res     = cmoic.KIRC,
                   dea.method   = "limma", # name of DEA method
                   prefix       = "TCGA-KIRC", # MUST be the same of argument in runDEA()
                   dat.path     = getwd(), # path of DEA files
                   res.path     = getwd(), # path to save GSEA files
                   msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                   norm.expr    = fpkm, # use normalized expression to calculate enrichment score
                   dirct        = "up", # direction of dysregulation in pathway
                   p.cutoff     = 0.05, # p cutoff to identify significant pathways
                   p.adj.cutoff = 0.15, # padj cutoff to identify significant pathways
                   gsva.method  = "gsva", # method to calculate single sample enrichment score
                   norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                   fig.name     = "Go UPREGULATED PATHWAY HEATMAP")
head(round(gsea.up$grouped.es,3))


library("clusterProfiler")
options(connectionObserver = NULL)
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

print(gsea.up$gsea.list$CS1[1:6,])

Go1<- pairwise_termsim(gsea.up$gsea.list$CS1)

pdf(file="Go-net-CS1.pdf",width = 10,height = 8)

cnetplot(Go1,
         showCategory = 5, 
         #foldChange = ID, 
         node_label="all", 
         colorEdge = TRUE,
         categorySize="pvalue")


dev.off()

Go2<- pairwise_termsim(gsea.up$gsea.list$CS2)

pdf(file="Go-net-CS2.pdf",width = 10,height = 8)

cnetplot(Go2,
         showCategory = 5, 
         #foldChange = ID, 
         node_label="all", 
         colorEdge = TRUE,
         categorySize="pvalue")

dev.off()

Go3<- pairwise_termsim(gsea.up$gsea.list$CS3)

pdf(file="Go-net-CS3.pdf",width = 14,height = 8)

cnetplot(Go3,
         showCategory = 5, 
         #foldChange = ID, 
         node_label="all", 
         colorEdge = TRUE,
         categorySize="pvalue")

dev.off()


GSET.FILE <- 
        system.file("extdata", "28immunemarker.gmt", package = "MOVICS", mustWork = TRUE)

gsva.res <- 
        runGSVA(moic.res      = cmoic.KIRC,
                norm.expr     = fpkm,
                gset.gmt.path = GSET.FILE, # ABSOLUTE path of gene set file
                gsva.method   = "gsva", # method to calculate single sample enrichment score
                annCol        = annCol,
                annColors     = annColors,
                color         =c("#30A9A0","white", "#E42053"),
                fig.path      = getwd(),
                fig.name      = "28 immune GENE SETS OF INTEREST HEATMAP",
                height        = 8,
                width         = 10)


#--------immune checkpoints 3 group
library(ggpubr)
finalsam<-cmoic.KIRC$clust.res$samID

ICB<-tpm[c("CD274","PDCD1","PDCD1LG2","CTLA4"),finalsam]
ICB<-log(ICB+1)

output<-tpm[,finalsam]
group<-as.data.frame(cmoic.KIRC$clust.res)

write.table(output,"expr225.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(group,"group225.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

ICB<-cbind(t(ICB),group)

p <- ggboxplot(ICB, x = "clust", y = "CTLA4",
               color = "clust", palette = c("#2EC4B6","#E71D36","#FF9F1C"),
               add = "jitter")
# Change method
p + stat_compare_means(method = "anova")

ggsave( "CTLA4  L.pdf", height = 4, width = 4)

p <- ggboxplot(ICB, x = "clust", y = "PDCD1",
               color = "clust", palette = c("#2EC4B6","#E71D36","#FF9F1C"),
               add = "jitter")
# Change method
p + stat_compare_means(method = "anova")

ggsave( "PDCD1  L.pdf", height = 4, width = 4)


#-----------metabolism selected

# load R package
library(GSVA)
library(ComplexHeatmap)
library(gplots)
library(clusterProfiler)
library(tidyr)
library(ggplot2)
library(estimate)
library(ggpubr)
source("twoclasslimma.R")
library("cogena")


# load meta signature
meta.sig <- gmt2list("metabolism signatures.gmt")
meta.sig <- sapply(meta.sig, function(x) setdiff(x,""))
meta.class <- NULL
for (i in names(meta.sig)) {
        tmp <- meta.sig[i]
        for (j in tmp) {
                meta.class <- rbind.data.frame(meta.class,
                                               data.frame(gene = j,
                                                          path = i,
                                                          stringsAsFactors = F),
                                               stringsAsFactors = F)
        }
}
expr = as.matrix(log2(mrna.tpm[,colnames(mut2)] + 1))

# calculate GSVA enrichment score
meta.score <- gsva(expr = expr,
                   gset.idx.list = meta.sig,
                   method = "gsva")
write.table(meta.score,"gsva enrichment score of metabolism signature in tcga cohort.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
meta.score<-t(meta.score)
# differential analysis
## CS1 vs others

tmp <- as.data.frame(cbind(cmoic.KIRC$clust.res$samID,cmoic.KIRC$clust.res$clust))

colnames(tmp)<-c("ID","Multi_omic_subtype")
rownames(tmp)<-tmp$ID

tmp$Multi_omic_subtype <- ifelse(tmp$Multi_omic_subtype == "1","CS1","Others")
subt <- data.frame(condition = tmp$Multi_omic_subtype,
                   row.names = rownames(tmp))


twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
              featmat  = meta.score[,rownames(subt)], # expression file (fill detect data scale automatically)
              treatVar = "CS1", # name of treatment group
              ctrlVar  = "Others", # name of control group
              prefix   = "tcga_Multi_omic_subtype", # prefix of the DE file
              overwt   = TRUE, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE) # if showing verbose result) # path for result


## CS2 vs Others
tmp<-as.data.frame(cbind(cmoic.KIRC$clust.res$samID,cmoic.KIRC$clust.res$clust))
colnames(tmp)<-c("ID","Multi_omic_subtype");rownames(tmp)<-tmp$ID
tmp$Multi_omic_subtype <- ifelse(tmp$Multi_omic_subtype == "2","CS2","Others")
subt <- data.frame(condition = tmp$Multi_omic_subtype,
                   row.names = rownames(tmp))
twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
              featmat  = meta.score[,rownames(subt)], # expression file (fill detect data scale automatically)
              treatVar = "CS2", # name of treatment group
              ctrlVar  = "Others", # name of control group
              prefix   = "tcga_Multi_omic_subtype", # prefix of the DE file
              overwt   = TRUE, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE, # if showing verbose result
              res.path =) # path for result

## CS3 vs Others
tmp<-as.data.frame(cbind(cmoic.KIRC$clust.res$samID,cmoic.KIRC$clust.res$clust))
colnames(tmp)<-c("ID","Multi_omic_subtype");rownames(tmp)<-tmp$ID
tmp$Multi_omic_subtype <- ifelse(tmp$Multi_omic_subtype == "3","CS3","Others")
subt <- data.frame(condition = tmp$Multi_omic_subtype,
                   row.names = rownames(tmp))
twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
              featmat  = meta.score[,rownames(subt)], # expression file (fill detect data scale automatically)
              treatVar = "CS3", # name of treatment group
              ctrlVar  = "Others", # name of control group
              prefix   = "tcga_Multi_omic_subtype", # prefix of the DE file
              overwt   = TRUE, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE, # if showing verbose result
              res.path = ) # path for result

# extract group specific pathways
tmp1 <- read.table("tcga_Multi_omic_subtype_limma_test_result.CS1_vs_Others.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tmp2 <- read.table("tcga_Multi_omic_subtype_limma_test_result.CS2_vs_Others.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tmp3 <- read.table("tcga_Multi_omic_subtype_limma_test_result.CS3_vs_Others.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

dep1 <- tmp1[which(tmp1$log2fc > 0.20 & tmp1$padj < 0.01),]
dep2 <- tmp2[which(tmp2$log2fc > 0.2 & tmp2$padj < 0.01),]
dep3 <- tmp3[which(tmp3$log2fc > 0.1 & tmp3$pvalue < 0.05),]

dep <- c(rownames(dep1),rownames(dep2),rownames(dep3))
dep <- setdiff(dep,dep[duplicated(dep)])

dep1 <- intersect(dep,rownames(dep1))
dep2 <- intersect(dep,rownames(dep2))
dep3 <- intersect(dep,rownames(dep3))

# generate metabolism pathway heatmap

Multi_omic_subtype <- as.data.frame(cbind(cmoic.KIRC$clust.res$samID,cmoic.KIRC$clust.res$clust))
colnames(Multi_omic_subtype)<-c("ID","Multi_omic_subtype")

annCol <- Multi_omic_subtype
annCol <- Multi_omic_subtype[order(Multi_omic_subtype$Multi_omic_subtype),]
annRow <- data.frame(path = rep(c("CS1-specific","CS2-specific","CS3-specific"),c(length(dep1),length(dep2),length(dep3))),
                     row.names = c(dep1,dep2,dep3))

annColors <- list()
annColors[["Multi_omic_subtype"]] <- c("1" = "#2EC4B6", "2" = "#E71D36", "3" = "#FF9F1C")
annColors[["path"]] <- c("CS1-specific" = "#2EC4B6", "CS2-specific" = "#E71D36", "CS3-specific" = "#FF9F1C")

indata <- meta.score[c(dep1,dep2,dep3),annCol$ID]

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
        outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
        if (!is.null(halfwidth)) {
                outdata[outdata>halfwidth]=halfwidth
                outdata[outdata<(-halfwidth)]= -halfwidth
        }
        return(outdata)
}

indata <- standarize.fun(indata, halfwidth = 1)
heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
pdf("metabolism genesets.pdf",width = 10,height = 4)
pheatmap(indata,
         border_color = NA,
         color = colorRampPalette(c("#31B29E", "black","#EB6C5A"))(50),
         annotation_col = annCol[,"Multi_omic_subtype",drop = F],
         #annotation_row = annRow,
         annotation_colors = annColors,
         show_rownames = T,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         gaps_row = cumsum(table(annRow$path))[1:2])

dev.off()

#-----------GO terms selected

# load R package
library(GSVA)
library(ComplexHeatmap)
library(gplots)
library(clusterProfiler)
library(tidyr)
library(ggplot2)
library(estimate)
library(ggpubr)
source("twoclasslimma.R")
library("cogena")


# load meta signature
meta.sig <- gmt2list("immune.gmt")
meta.sig <- sapply(meta.sig, function(x) setdiff(x,""))
meta.class <- NULL
for (i in names(meta.sig)) {
        tmp <- meta.sig[i]
        for (j in tmp) {
                meta.class <- rbind.data.frame(meta.class,
                                               data.frame(gene = j,
                                                          path = i,
                                                          stringsAsFactors = F),
                                               stringsAsFactors = F)
        }
}
expr = as.matrix(log2(mrna.tpm[,colnames(mut2)] + 1))

# calculate GSVA enrichment score
meta.score <- gsva(expr = expr,
                   gset.idx.list = meta.sig,
                   method = "gsva")
write.table(meta.score,"immune gsva enrichment score of metabolism signature in tcga cohort.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
meta.score<-t(meta.score)
# differential analysis
## CS1 vs others

tmp <- as.data.frame(cbind(cmoic.KIRC$clust.res$samID,cmoic.KIRC$clust.res$clust))

colnames(tmp)<-c("ID","Multi_omic_subtype")
rownames(tmp)<-tmp$ID

tmp$Multi_omic_subtype <- ifelse(tmp$Multi_omic_subtype == "1","CS1","Others")
subt <- data.frame(condition = tmp$Multi_omic_subtype,
                   row.names = rownames(tmp))


twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
              featmat  = meta.score[,rownames(subt)], # expression file (fill detect data scale automatically)
              treatVar = "CS1", # name of treatment group
              ctrlVar  = "Others", # name of control group
              prefix   = "immune_Multi_omic_subtype", # prefix of the DE file
              overwt   = TRUE, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE) # if showing verbose result) # path for result


## CS2 vs Others
tmp<-as.data.frame(cbind(cmoic.KIRC$clust.res$samID,cmoic.KIRC$clust.res$clust))
colnames(tmp)<-c("ID","Multi_omic_subtype");rownames(tmp)<-tmp$ID
tmp$Multi_omic_subtype <- ifelse(tmp$Multi_omic_subtype == "2","CS2","Others")
subt <- data.frame(condition = tmp$Multi_omic_subtype,
                   row.names = rownames(tmp))
twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
              featmat  = meta.score[,rownames(subt)], # expression file (fill detect data scale automatically)
              treatVar = "CS2", # name of treatment group
              ctrlVar  = "Others", # name of control group
              prefix   = "immune_Multi_omic_subtype", # prefix of the DE file
              overwt   = TRUE, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE, # if showing verbose result
              res.path =) # path for result

## CS3 vs Others
tmp<-as.data.frame(cbind(cmoic.KIRC$clust.res$samID,cmoic.KIRC$clust.res$clust))
colnames(tmp)<-c("ID","Multi_omic_subtype");rownames(tmp)<-tmp$ID
tmp$Multi_omic_subtype <- ifelse(tmp$Multi_omic_subtype == "3","CS3","Others")
subt <- data.frame(condition = tmp$Multi_omic_subtype,
                   row.names = rownames(tmp))
twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
              featmat  = meta.score[,rownames(subt)], # expression file (fill detect data scale automatically)
              treatVar = "CS3", # name of treatment group
              ctrlVar  = "Others", # name of control group
              prefix   = "immune_Multi_omic_subtype", # prefix of the DE file
              overwt   = TRUE, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE, # if showing verbose result
              res.path = ) # path for result

# extract group specific pathways
tmp1 <- read.table("immune_Multi_omic_subtype_limma_test_result.CS1_vs_Others.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tmp2 <- read.table("immune_Multi_omic_subtype_limma_test_result.CS2_vs_Others.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tmp3 <- read.table("immune_Multi_omic_subtype_limma_test_result.CS3_vs_Others.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

dep1 <- tmp1[which(tmp1$log2fc > 0.2 & tmp1$pvalue < 0.05),]
dep2 <- tmp2[which(tmp2$log2fc > 0.1 & tmp2$pvalue < 0.05),]
dep3 <- tmp3[which(tmp3$log2fc > 0.1 & tmp3$pvalue < 0.05),]

dep <- c(rownames(dep1),rownames(dep2),rownames(dep3))
dep <- setdiff(dep,dep[duplicated(dep)])

dep1 <- intersect(dep,rownames(dep1))
dep2 <- intersect(dep,rownames(dep2))
dep3 <- intersect(dep,rownames(dep3))

# generate metabolism pathway heatmap

Multi_omic_subtype <- as.data.frame(cbind(cmoic.KIRC$clust.res$samID,cmoic.KIRC$clust.res$clust))
colnames(Multi_omic_subtype)<-c("ID","Multi_omic_subtype")

annCol <- Multi_omic_subtype
annCol <- Multi_omic_subtype[order(Multi_omic_subtype$Multi_omic_subtype),]
annRow <- data.frame(path = rep(c("CS1-specific","CS2-specific","CS3-specific"),c(length(dep1),length(dep2),length(dep3))),
                     row.names = c(dep1,dep2,dep3))

annColors <- list()
annColors[["Multi_omic_subtype"]] <- c("1" = "#2EC4B6", "2" = "#E71D36", "3" = "#FF9F1C")
annColors[["path"]] <- c("CS1-specific" = "#2EC4B6", "CS2-specific" = "#E71D36", "CS3-specific" = "#FF9F1C")

indata <- meta.score[c(dep1,dep2,dep3),annCol$ID]

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
        outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
        if (!is.null(halfwidth)) {
                outdata[outdata>halfwidth]=halfwidth
                outdata[outdata<(-halfwidth)]= -halfwidth
        }
        return(outdata)
}

indata <- standarize.fun(indata, halfwidth = 1)
heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
pdf("immune28 genesets.pdf",width = 10,height = 4)
pheatmap(indata,
         border_color = NA,
         color = colorRampPalette(c("#2295FF", "#FFFFFF","#E00115"))(50),
         annotation_col = annCol[,"Multi_omic_subtype",drop = F],
         #annotation_row = annRow,
         annotation_colors = annColors,
         show_rownames = T,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         gaps_row = cumsum(table(annRow$path))[1:2])

dev.off()



#-----------------prepare validate cohort

MSKCC.ntp.pred <- runNTP(expr      = MSKCC.expr,
                         templates = marker.up$templates, # the template has been already prepared in runMarker()
                         scale     = TRUE, # scale input data (by default)
                         center    = TRUE, # center input data (by default)
                         doPlot    = TRUE, # to generate heatmap
                         fig.name  = "NTP HEATMAP FOR MSKCC") 




MSKCC.expr<-read.table("MSKCC_PCa_mRNA_log2_data.txt",header=T,sep="\t",row.names=1,check.names=F)
MSKCC.clin<-read.table("MSKCC patients clinical.txt",header=T,sep="\t",row.names=1,check.names=F)
MSKCC.expr<-MSKCC.expr[,rownames(MSKCC.clin)]


MSKCC.ntp.pred <- runNTP(expr      = MSKCC.expr,
                         templates = marker.up$templates, # the template has been already prepared in runMarker()
                         scale     = TRUE, # scale input data (by default)
                         center    = TRUE, # center input data (by default)
                         doPlot    = TRUE, # to generate heatmap
                         fig.name  = "NTP HEATMAP FOR MSKCC") 

surv.MSKCC <- compSurv(moic.res         = MSKCC.ntp.pred,
                       surv.info        = MSKCC.clin,
                       convt.time       = "m", # switch to year
                       surv.median.line = "hv", # switch to both
                       fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR MSKCC") 


runDEA(dea.method = "limma",
       expr       = MSKCC.expr, # normalized expression data
       moic.res   = MSKCC.ntp.pred,
       prefix     = "MSKCC")




MSKCC.gsea.up <- runGSEA(moic.res     = MSKCC.ntp.pred,
                         dea.method   = "limma", # name of DEA method
                         prefix       = "MSKCC", # MUST be the same of argument in runDEA()
                         dat.path     = getwd(), # path of DEA files
                         res.path     = getwd(), # path to save GSEA files
                         msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                         norm.expr    = MSKCC.expr, # use normalized expression to calculate enrichment score
                         dirct        = "up", # direction of dysregulation in pathway
                         p.cutoff     = 0.05, # p cutoff to identify significant pathways
                         p.adj.cutoff = 0.15, # padj cutoff to identify significant pathways
                         gsva.method  = "gsva", # method to calculate single sample enrichment score
                         norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                         fig.name     = "MSKCC UPREGULATED PATHWAY HEATMAP")

MSKCC.gsea.dn <- runGSEA(moic.res     = MSKCC.ntp.pred,
                         dea.method   = "limma",
                         prefix       = "MSKCC",
                         msigdb.path  = MSIGDB.FILE,
                         norm.expr    = MSKCC.expr,
                         dirct        = "down",
                         p.cutoff     = 0.05,
                         p.adj.cutoff = 0.15,
                         gsva.method  = "ssgsea", # switch to ssgsea
                         norm.method  = "mean", # switch to median
                         fig.name     = "MSKCC DOWNREGULATED PATHWAY HEATMAP") 


#---------multi-survival
library("tidyverse")
library(survival)
library(survminer)

CSG<-MSKCC.ntp.pred$clust.res$clust
CGST<-cbind(MSKCC.clin,CSG)
CGST[CGST==""]<-NA
CGST<-CGST %>% 
        drop_na()

CGST$CSG<-paste0("CS",CGST$CSG)
CGST$CSG<-ifelse(CGST$CSG=="CS3","ACS3",CGST$CSG)
CGST$Gleason<-paste0("X",CGST$Gleason)
CGST$PSA<-ifelse(CGST$PSA>4,">4","<=4")

#æž„å»ºæ¨¡åž‹
model <- coxph( Surv(futime, fustat) ~., data =CGST,na.action=na.exclude )
multiCoxSum = summary(model)
multiCoxSum

outTab=data.frame()
outTab=cbind(
        HR=multiCoxSum$conf.int[,"exp(coef)"],
        HR.95L=multiCoxSum$conf.int[,"lower .95"],
        HR.95H=multiCoxSum$conf.int[,"upper .95"],
        pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)

write.table(outTab,file="MSKCCmultiCox.xls",sep="\t",row.names=F,quote=F)


#-----------------prepare validate cohort

GSE29609.expr<-read.table("GSE29609.expr.txt",header=T,sep="\t",row.names=1,check.names=F)
GSE29609.clin<-read.table("GSE29609.clin.txt",header=T,sep="\t",row.names=1,check.names=F)

GSE29609.ntp.pred <- runNTP(expr      = GSE29609.expr,
                             templates = marker.up$templates, # the template has been already prepared in runMarker()
                             scale     = TRUE, # scale input data (by default)
                             center    = TRUE, # center input data (by default)
                             doPlot    = TRUE, # to generate heatmap
                             fig.name  = "NTP HEATMAP FOR GSE29609") 

surv.GSE29609 <- compSurv(moic.res         = GSE29609.ntp.pred,
                           surv.info        = GSE29609.clin,
                           convt.time       = "m", # switch to year
                           surv.median.line = "hv", # switch to both
                           fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR GSE29609") 




#---------multi-survival
library("tidyverse")
library(survival)
library(survminer)

CSG<-GSE29609.ntp.pred$clust.res$clust
CGST<-cbind(GSE29609.clin,CSG)
CGST[CGST==""]<-NA
CGST[CGST=="Unknown"]<-NA
CGST<-CGST %>% 
        drop_na()

CGST<-CGST[,c(10,11,14,15,16,18,20)]
colnames(CGST)<-c("Gleason","PSA","fustat","Tstage","futime","Age","CSG")
CGST<-CGST[,c(1,3:7)]

CGST$CSG<-paste0("CS",CGST$CSG)
CGST$CSG<-ifelse(CGST$CSG=="CS3","ACS3",CGST$CSG)
CGST$Gleason<-paste0("X",CGST$Gleason)
CGST$PSA<-ifelse(CGST$PSA>4,">4","<=4")

#æž„å»ºæ¨¡åž‹
model <- coxph( Surv(futime, fustat) ~., data =CGST,na.action=na.exclude )
multiCoxSum = summary(model)
multiCoxSum

outTab=data.frame()
outTab=cbind(
        HR=multiCoxSum$conf.int[,"exp(coef)"],
        HR.95L=multiCoxSum$conf.int[,"lower .95"],
        HR.95H=multiCoxSum$conf.int[,"upper .95"],
        pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)

write.table(outTab,file="GSE29609multiCox.xls",sep="\t",row.names=F,quote=F)


#-----------------prepare validate cohort

GSE22541.expr<-read.table("GSE22541.expr.txt",header=T,sep="\t",row.names=1,check.names=F)
GSE22541.clin<-read.table("GSE22541.clin.txt",header=T,sep="\t",row.names=1,check.names=F)
GSE22541.expr<-GSE22541.expr[,rownames(GSE22541.clin)]

GSE22541.ntp.pred <- runNTP(expr      = GSE22541.expr,
                            templates = marker.up$templates, # the template has been already prepared in runMarker()
                            scale     = TRUE, # scale input data (by default)
                            center    = TRUE, # center input data (by default)
                            doPlot    = TRUE, # to generate heatmap
                            fig.name  = "NTP HEATMAP FOR GSE22541") 

surv.GSE22541 <- compSurv(moic.res         = GSE22541.ntp.pred,
                          surv.info        = GSE22541.clin,
                          convt.time       = "m", # switch to year
                          surv.median.line = "hv", # switch to both
                          fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR GSE22541") 




#---------multi-survival
library("tidyverse")
library(survival)
library(survminer)
GSE22541.clin<-read.table("GSE22541_clin.txt",header=T,sep="\t",row.names=1,check.names=F)


CSG<-GSE22541.ntp.pred$clust.res$clust
CGST<-cbind(GSE22541.clin,CSG)
CGST[CGST==""]<-NA
CGST[CGST=="unknown"]<-NA
CGST[CGST=="UNKNOWN"]<-NA
CGST<-CGST[,c(1:2,4:8)]
CGST<-CGST %>% 
        drop_na()

colnames(CGST)<-c("fustat","futime","Gleason","PSA","Tstage","Margin","CSG")

CGST$CSG<-paste0("CS",CGST$CSG)
CGST$CSG<-ifelse(CGST$CSG=="CS3","ACS3",CGST$CSG)
CGST$Gleason<-paste0("X",CGST$Gleason)

CGST$Gleason<-ifelse(CGST$Gleason=="X5","X6",CGST$Gleason)
CGST$PSA<-ifelse(CGST$PSA>4,">4","<=4")

#æž„å»ºæ¨¡åž‹
model <- coxph( Surv(futime, fustat) ~., data =CGST,na.action=na.exclude )
multiCoxSum = summary(model)
multiCoxSum

outTab=data.frame()
outTab=cbind(
        HR=multiCoxSum$conf.int[,"exp(coef)"],
        HR.95L=multiCoxSum$conf.int[,"lower .95"],
        HR.95H=multiCoxSum$conf.int[,"upper .95"],
        pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)

write.table(outTab,file="GSE22541multiCox.xls",sep="\t",row.names=F,quote=F)


#-----------------prepare validate cohort

EMTAB3267.expr<-read.table("EMTAB3267.expr.txt",header=T,sep="\t",row.names=1,check.names=F)
EMTAB3267.expr<-t(EMTAB3267.expr)
EMTAB3267.clin<-read.table("EMTAB3267.clin.txt",header=T,sep="\t",row.names=1,check.names=F)
EMTAB3267.expr<-EMTAB3267.expr[,rownames(EMTAB3267.clin)]
EMTAB3267.clin$`sunitinib-response`<-ifelse(EMTAB3267.clin$`sunitinib-response`=="PR","DPR",EMTAB3267.clin$`sunitinib-response`)


EMTAB3267.ntp.pred <- runNTP(expr      = EMTAB3267.expr,
                            templates = marker.up$templates, # the template has been already prepared in runMarker()
                            scale     = TRUE, # scale input data (by default)
                            center    = TRUE, # center input data (by default)
                            doPlot    = TRUE, # to generate heatmap
                            fig.name  = "NTP HEATMAP FOR EMTAB3267") 

surv.EMTAB3267 <- compSurv(moic.res         = EMTAB3267.ntp.pred,
                          surv.info        = EMTAB3267.clin,
                          convt.time       = "m", # switch to year
                          surv.median.line = "hv", # switch to both
                          fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR EMTAB3267") 


test.data<-table(EMTAB3267.ntp.pred$clust.res$clust,EMTAB3267.clin$`sunitinib-response`)
test.data
fisher.test(test.data)
fisher.test(test.data)$expected
ax<-fisher.test(test.data)$p.value

pvalue= paste(pvalue= ifelse (ax <0.001, " P < 0.001", paste ("P = ",round(ax,3),sep = "")))

test.data2<-data.frame(test.data)
test.data2

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
        scale_fill_brewer(palette="Dark2")

p

p2<-p+geom_text(aes(x=2,y=0.8,label= pvalue), vjust=0.5, hjust=0.5, size=6, fontface="italic",color="black")
p2


pdf("EMTAB3267.sunitinib response.pdf",width=6,height=4,onefile = FALSE)
p2
dev.off()


runDEA(dea.method = "limma",
       expr       = EMTAB3267.expr, # normalized expression data
       moic.res   = EMTAB3267.ntp.pred,
       prefix     = "EMTAB3267")




EMTAB3267.gsea.up <- runGSEA(moic.res     = EMTAB3267.ntp.pred,
                            dea.method   = "limma", # name of DEA method
                            prefix       = "EMTAB3267", # MUST be the same of argument in runDEA()
                            dat.path     = getwd(), # path of DEA files
                            res.path     = getwd(), # path to save GSEA files
                            msigdb.path  = MSIGDB.FILE, # MUST be the ABSOLUTE path of msigdb file
                            norm.expr    = EMTAB3267.expr, # use normalized expression to calculate enrichment score
                            dirct        = "up", # direction of dysregulation in pathway
                            p.cutoff     = 0.05, # p cutoff to identify significant pathways
                            p.adj.cutoff = 0.15, # padj cutoff to identify significant pathways
                            gsva.method  = "gsva", # method to calculate single sample enrichment score
                            norm.method  = "mean", # normalization method to calculate subtype-specific enrichment score
                            fig.name     = "EMTAB3267 UPREGULATED PATHWAY HEATMAP")

EMTAB3267.gsea.dn <- runGSEA(moic.res     = EMTAB3267.ntp.pred,
                            dea.method   = "limma",
                            prefix       = "EMTAB3267",
                            msigdb.path  = MSIGDB.FILE,
                            norm.expr    = EMTAB3267.expr,
                            dirct        = "down",
                            p.cutoff     = 0.05,
                            p.adj.cutoff = 0.15,
                            gsva.method  = "ssgsea", # switch to ssgsea
                            norm.method  = "mean", # switch to median
                            fig.name     = "EMTAB3267 DOWNREGULATED PATHWAY HEATMAP") 

#---------multi-survival
library("tidyverse")
library(survival)
library(survminer)

EMTAB3267.clin<-read.table("EMTAB3267.clin.txt",header=T,sep="\t",row.names=1,check.names=F)

CSG<-EMTAB3267.ntp.pred$clust.res$clust
CGST<-cbind(EMTAB3267.clin,CSG)

CGST[CGST==""]<-NA
CGST[CGST=="unknown"]<-NA
CGST<-CGST %>% 
        drop_na()

colnames(CGST)<-c("fustat","futime","Gleason","Age","PSA","Tstage","Margin","CSG")

CGST$CSG<-paste0("CS",CGST$CSG)
CGST$CSG<-ifelse(CGST$CSG=="CS3","ACS3",CGST$CSG)
CGST$Gleason<-paste0("X",CGST$Gleason)
CGST$PSA<-ifelse(CGST$PSA>4,">4","<=4")

#æž„å»ºæ¨¡åž‹
model <- coxph( Surv(futime, fustat) ~., data =CGST,na.action=na.exclude )
multiCoxSum = summary(model)
multiCoxSum

outTab=data.frame()
outTab=cbind(
        HR=multiCoxSum$conf.int[,"exp(coef)"],
        HR.95L=multiCoxSum$conf.int[,"lower .95"],
        HR.95H=multiCoxSum$conf.int[,"upper .95"],
        pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)

write.table(outTab,file="EMTAB3267multiCox.xls",sep="\t",row.names=F,quote=F)

#----------------Miao-prepare validate cohort

Miao.expr<-read.table("Miao.expr.txt",header=T,sep="\t",row.names=1,check.names=F)
Miao.expr<-log2(Miao.expr+1)
Miao.clin<-read.table("Miao.clin.txt",header=T,sep="\t",row.names=1,check.names=F)

Miao.expr<-Miao.expr[,rownames(Miao.clin)]

Miao.ntp.pred <- runNTP(expr      = Miao.expr,
                        templates = marker.up$templates, # the template has been already prepared in runMarker()
                        scale     = TRUE, # scale input data (by default)
                        center    = TRUE, # center input data (by default)
                        
                        doPlot    = TRUE, # to generate heatmap
                        fig.name  = "NTP HEATMAP FOR Miao") 


Miao.group<-cbind(Miao.clin,Miao.ntp.pred$clust.res[rownames(Miao.clin),])

table(Miao.group$response_category,Miao.group$clust)

surv.Miao <- compSurv(moic.res         = Miao.ntp.pred,
                      surv.info        = Miao.clin,
                      convt.time       = "m", # switch to year
                      surv.median.line = "hv", # switch to both
                      fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR Miao") 

merge<-intersect(rownames(Miao.clin),Miao.ntp.pred$clust.res$samID)

Miao.clin<-Miao.clin[merge,]
Miaoclust<-Miao.ntp.pred$clust.res
Miaoclust<-Miaoclust[merge,]

Miaoout<-cbind(Miao.clin,Miaoclust)

test.data<-print(table(Miaoout$clust,Miaoout$response_category))

ax<-fisher.test(test.data)$p.value

test.data2<-data.frame(test.data)
test.data2<-test.data2[order(test.data2[,1]),] 
test.data2

subset1 = test.data2[1:2,]
subset1$pct = subset1$Freq/sum(subset1$Freq)
subset1

subset2 = test.data2[3:4,]
subset2$pct = subset2$Freq/sum(subset2$Freq)
subset2

subset3 = test.data2[5:6,]
subset3$pct = subset3$Freq/sum(subset3$Freq)
subset3


test.data3<-rbind.data.frame(subset1,subset2,subset3)
test.data3

library(wesanderson)
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

p


pvalue= paste(pvalue= ifelse (ax <0.001, " P < 0.001", paste ("P = ",round(ax,3),sep = "")))
p2<-p+geom_text(aes(x=2,y=0.1,label= pvalue), vjust=0.5, hjust=0.5, size=6, fontface="italic",color="black")
p2


pdf("Miao response.pdf", width=5,height=4,onefile = FALSE)
p
dev.off()

#----------------CheckMate-prepare validate cohort

CheckMate.expr<-read.table("CheckMate.expr.txt",header=T,sep="\t",row.names=1,check.names=F)


CheckMate.clin<-read.table("CheckMate.clin.txt",header=T,sep="\t",row.names=1,check.names=F)

CheckMate.expr2<-CheckMate.expr[,rownames(CheckMate.clin)]

CheckMate.ntp.pred <- runNTP(expr      = CheckMate.expr2,
                              templates = marker.up$templates, # the template has been already prepared in runMarker()
                             scaleFlag = TRUE,
                             centerFlag = TRUE,
                             nPerm = 1000,
                             distance = "cosine",
                             seed = 123456,
                             verbose = TRUE,
                            
                              doPlot    = TRUE, # to generate heatmap
                              fig.name  = "NTP HEATMAP FOR CheckMate-2") 
surv.CheckMate <- compSurv(moic.res         = CheckMate.ntp.pred,
                      surv.info        = CheckMate.clin,
                      convt.time       = "m", # switch to year
                      surv.median.line = "hv", # switch to both
                      fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR CheckMate") 


merge<-intersect(rownames(CheckMate.clin),CheckMate.ntp.pred$clust.res$samID)

CheckMate.clin<-CheckMate.clin[merge,]
IMclust<-CheckMate.ntp.pred$clust.res
IMclust<-IMclust[merge,]

IMout<-cbind(CheckMate.clin,IMclust)




test.data<-print(table(IMout$clust,IMout$Benefit))

ax<-chisq.test(test.data)$p.value

test.data2<-data.frame(test.data)
test.data2<-test.data2[order(test.data2[,1]),] 
test.data2

subset1 = test.data2[1:2,]
subset1$pct = subset1$Freq/sum(subset1$Freq)
subset1

subset2 = test.data2[3:4,]
subset2$pct = subset2$Freq/sum(subset2$Freq)
subset2

subset3 = test.data2[5:6,]
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
        scale_fill_manual(values = wes_palette("BottleRocket2", n = 2))

p


pvalue= paste(pvalue= ifelse (ax <0.001, " P < 0.001", paste ("P = ",round(ax,3),sep = "")))
p2<-p+geom_text(aes(x=2,y=0.1,label= pvalue), vjust=0.5, hjust=0.5, size=6, fontface="italic",color="black")
p2


pdf("CheckMate-1 Immunotherapy response.pdf", width=5,height=4,onefile = FALSE)
p
dev.off()



#----------------GSE40435-prepare validate cohort

GSE40435.expr<-read.table("GSE40435.expr.txt",header=T,sep="\t",row.names=1,check.names=F)


GSE40435.clin<-read.table("GSE40435.clin.txt",header=T,sep="\t",row.names=1,check.names=F)

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
                        fig.name  = "NTP HEATMAP FOR GSE40435-2") 
surv.GSE40435 <- compSurv(moic.res         = GSE40435.ntp.pred,
                      surv.info        = GSE40435.clin,
                      convt.time       = "m", # switch to year
                      surv.median.line = "hv", # switch to both
                      fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR GSE40435") 


merge<-intersect(rownames(GSE40435.clin),GSE40435.ntp.pred$clust.res$samID)

GSE40435.clin<-GSE40435.clin[merge,]
IMclust<-GSE40435.ntp.pred$clust.res
IMclust<-IMclust[merge,]

IMout<-cbind(GSE40435.clin,IMclust)




test.data<-print(table(IMout$clust,IMout$`1:GRADE`))

ax<-chisq.test(test.data)$p.value

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

p


pvalue= paste(pvalue= ifelse (ax <0.001, " P < 0.001", paste ("P = ",round(ax,3),sep = "")))
p2<-p+geom_text(aes(x=2,y=0.1,label= pvalue), vjust=0.5, hjust=0.5, size=6, fontface="italic",color="black")
p2


pdf("GSE40435 Group and stage.pdf", width=5,height=4,onefile = FALSE)
p
dev.off()

#----------------GSE53757-prepare validate cohort

GSE53757.expr<-read.table("GSE53757.expr.txt",header=T,sep="\t",row.names=1,check.names=F)


GSE53757.clin<-read.table("GSE53757.clin.txt",header=T,sep="\t",row.names=1,check.names=F)

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
                            fig.name  = "NTP HEATMAP FOR GSE53757-2") 


merge<-intersect(rownames(GSE53757.clin),GSE53757.ntp.pred$clust.res$samID)

GSE53757.clin<-GSE53757.clin[merge,]
IMclust<-GSE53757.ntp.pred$clust.res
IMclust<-IMclust[merge,]

IMout<-cbind(GSE53757.clin,IMclust)




test.data<-print(table(IMout$clust,IMout$STAGE))

ax<-fisher.test(test.data)$p.value

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

p


pvalue= paste(pvalue= ifelse (ax <0.001, " P < 0.001", paste ("P = ",round(ax,3),sep = "")))
p2<-p+geom_text(aes(x=2,y=0.1,label= pvalue), vjust=0.5, hjust=0.5, size=6, fontface="italic",color="black")
p2


pdf("GSE53757 Group and stage.pdf", width=5,height=4,onefile = FALSE)
p
dev.off()


#----------Cibersort--------------
write.table(log2(output+1),file="exprMat.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

source("CIBERSORT/CIBERSORT.R") 
# Ö÷ÒªÊÇ read.table µÄ²ÎÊý¿ÉÒÔËæÒâÐÞ¸Ä£¬¸ù¾ÝÄã×Ô¼ºµÄ±í´ï¾ØÕótxtÎÄ¼þÊÊÓ¦ÐÔµ÷Õû¼´¿É
## sig_matrix <- read.table("CIBERSORT/LM22.txt",header=T,sep="\t",row.names=1,check.names=F)
sig_matrix <- "CIBERSORT/LM22.txt"
mixture_file = 'exprMat.txt'
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=10, QN=TRUE)


immunocyte<-as.data.frame(cbind(cmoic.KIRC$clust.res,res_cibersort[cmoic.KIRC$clust.res$samID,]))

rt<-immunocyte[order(immunocyte$clust),]
table(rt$clust)
rt<-rt[,c(3:6,8:24)]

write.table(rt,file="Cibersortgroup.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
library(ggplot2)
library(ggpubr)
conNum=58                                                     
treatNum=86
treatNum2=81                                                 

type=c(rep("CS1",conNum),rep("CS2",treatNum),rep("CS3",treatNum2))

data=data.frame()
for(i in colnames(rt)){
        data=rbind(data,cbind(expression=rt[,i],gene=i,type))
}
write.table(data,file="data.txt",sep="\t",row.names=F,quote=F)

data=read.table("data.txt",sep="\t",header=T,check.names=F)      
p=ggboxplot(data, x="gene", y="expression", notch = F,
            lwd=0.5, alpha = 0.6,color = "type", fill = 'type', 
            ylab="Immunocyte Infiltration",
            xlab="",
            palette = c("#2EC4B6","#E71D36","#FFD121") )
p=p+rotate_x_text(60)
pdf(file="boxplot.pdf",width=10,height=5)                     
p+stat_compare_means(aes(group=type),method = "kruskal.test", symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),label = "p.signif")
dev.off()




#-------------ssGSEAÍ¨Â·Í¼----------

cellMarker <- data.table::fread("mitoch.csv",data.table = F)
table(cellMarker$celltype)
length(table(cellMarker$celltype))
colnames(cellMarker)[2] <- "celltype"

cellMarker <- lapply(split(cellMarker,cellMarker$celltype), function(x){
        dd = x$Metagene
        unique(dd)
})


### 2.×¼??????ï¿½ï¿½????
### 
expr <- mrna.tpm
#rownames(expr) <- expr[,1]
expr <- expr[,-1]
expr <- as.matrix(expr)

### 3.Ê¹??ssGSEAï¿½ï¿½?????ß½???
library(GSVA)
gsva_data <- gsva(expr,cellMarker, method = "ssgsea")
gsva_data <- t(gsva_data)

write.table(gsva_data,file="MOVICS mitoch.xls",sep = "\t",row.names = T,quote = F)

clust<-cmoic.KIRC$clust.res
rownames<-intersect(rownames(gsva_data), rownames(clust))
gsva_data<-gsva_data[rownames,]
clust<-clust[rownames,]

clust$Group<-as.factor(clust$clust)
gsva_data<-cbind(clust,gsva_data)


outTab <- NULL

pl<-list()
for(gene in colnames(gsva_data)[4:ncol(gsva_data)]){
        a<- gsva_data[,c(gene,"Group")]
        colnames(a)[1] <- "gene"
        Ktest <- kruskal.test(gene~Group, data=a)
        outTab=rbind.data.frame(outTab, data.frame(gene = gene,
                                                   p = Ktest$p.value,
                                                   stringsAsFactors = F),
                                stringsAsFactors = F)
        
}
outTab
write.table(outTab,file="MOVICS 28immunocytes K-M 3group.xls",sep = "\t",row.names = T,quote = F)



library(dplyr)
library(pheatmap)  

test<-gsva_data[order(gsva_data$Group),]
test<-test[,3:ncol(test)]

rt3<-test[,2:ncol(test)]
rt3<-t(rt3)

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
        outdata=t(scale(t(rt3), center=centerFlag, scale=scaleFlag))
        if (!is.null(halfwidth)) {
                outdata[outdata>halfwidth]=halfwidth
                outdata[outdata<(-halfwidth)]= -halfwidth
        }
        return(outdata)
}

rt3 <- standarize.fun(rt3,halfwidth =0.8)

rt4 <- test %>% 
        dplyr::select(c(Group))

cluster<-rt4
pdf("ClustOutPut.pdf",height=3,width=10)
p_heatmap<-pheatmap(rt3, annotation_col=cluster, 
                    color = colorRampPalette(c("#0099CC", "white","#CC0033"))(50),
                    cluster_cols =T,
                    fontsize=8,
                    fontsize_row=8,
                    scale="row",
                    show_colnames=F,
                    fontsize_col=3)
dev.off()

p_heatmap


#------------------heatmap---------------------------


library(dplyr)
library(pheatmap)                   #å¼•ç”¨åŒ?


rt=read.table("Comparesig\\PAM50-KIRC.txt",header=T,sep="\t",row.names=1,check.names=F)
rt1=read.table("Comparesig\\PCa molecular subtype.txt",header=T,sep="\t",row.names=1,check.names=F)
PAM50merge<-intersect(rownames(rt1),colnames(rt))

rt<-log2(rt+1)
rt2<-rt[,PAM50merge]


rt2 <- scale(rt2) #scaleæ ‡åŒ–

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
        outdata=t(scale(t(rt2), center=centerFlag, scale=scaleFlag))
        if (!is.null(halfwidth)) {
                outdata[outdata>halfwidth]=halfwidth
                outdata[outdata<(-halfwidth)]= -halfwidth
        }
        return(outdata)
}

rt3 <- standarize.fun(rt2,halfwidth =1)

head(rt3[1:3,1:3])

#è¯»å–åˆ†ç»„æ–‡ä»¶


rt4 <- rt1[PAM50merge,]

#rt4$riskscore<- ifelse(rt4$riskscore>=median(rt4$riskscore), "High-risk", "Low-risk")

#rt4$fustat<-rt4[order(rt4$fustat,decreasing = T),] 

cluster=rt4
library("pheatmap")


annColors <- list(Multi_omic_subtype = c("CS1"="#2EC4B6",
                                         "CS2"="#E71D36",
                                         "CS3"="#FF9F1C"),
                  PAM50 = c("LumA"="#9BCS2E6",
                            "LumB"="#7030A0",
                            "Basal"="#FFD966"),
                  Immunophenotype = c("Non-Immune"="#D8D8D8",
                                      "Suppressed"="#70AD47",
                                      "Activated"="#FF0000"),
                  ImmuneSubtype = c("CS1"="#E71D36",
                                    "CS2"="yellow",
                                    "CS3"="#00B050",
                                    "C4"="#00B0F0",
                                    "Unknown"="#D8D8D8")
)

#ç»˜åˆ¶çƒ­å›¾
pdf("Comparesig\\multi-KIRC_ClustOutPut.pdf",height=4,width=9)
p_heatmap<-pheatmap(rt3, annotation=cluster, 
                    color = colorRampPalette(c("#31B29E", "black","#EB6C5A"))(50),
                    cluster_cols =F,
                    fontsize=8,
                    fontsize_row=8,
                    scale="row",
                    annotation_colors=annColors,
                    show_colnames=F,
                    fontsize_col=3)
dev.off()



#---------multi-survival


st=read.table("Comparesig\\PCa molecular subtype.txt",header=T,sep="\t",row.names=1,check.names=F)
st<-st[rownames(surv.info),]
surv.info<-surv.info[rownames(st),]

newsurv<-cbind(surv.info,st$Multi_omic_subtype)

library(survival)
library(survminer)

newsurv$`st$Multi_omic_subtype`<-ifelse(newsurv$`st$Multi_omic_subtype`=="CS3","ACS3",newsurv$`st$Multi_omic_subtype`)
newsurv$Gleason<-ifelse(newsurv$Gleason=="X9","X9+10",newsurv$Gleason)
newsurv$Gleason<-ifelse(newsurv$Gleason=="X10","X9+10",newsurv$Gleason)
#æž„å»ºæ¨¡åž‹
model <- coxph( Surv(futime, fustat) ~., data =newsurv,na.action=na.exclude )
multiCoxSum = summary(model)
multiCoxSum

outTab=data.frame()
outTab=cbind(
        HR=multiCoxSum$conf.int[,"exp(coef)"],
        HR.95L=multiCoxSum$conf.int[,"lower .95"],
        HR.95H=multiCoxSum$conf.int[,"upper .95"],
        pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)

write.table(outTab,file="KIRCmultiCox.xls",sep="\t",row.names=F,quote=F)

pdf(file="Comparesig\\TOTAL clinical forest.pdf",onefile = FALSE,
    width = 6,             #å›¾ç‰‡çš„å®½åº?
    height = 6,            #å›¾ç‰‡çš„é«˜åº?
)

ggforest(model,
         data=newsurv,
         main = "Entire Cohort",
         cpositions = c(0.01,0.14,0.36), 
         fontsize = 1, 
         refLabel = "reference", 
         noDigits = 3)
dev.off()

#------sankey
df<-read.table("PCa molecular subtype.txt",header=T,sep="\t",row.names=1,check.names=F)

subdf <- df
library(ggalluvial)
library(dplyr)
mycol <- rep(c("#2EC4B6","#E71D36","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#029149","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767"),2)
subdf <- subdf %>% 
        group_by(Multi_omic_subtype, PAM50, Immunophenotype, ImmuneSubtype) %>% 
        tally(name = "Freq") %>% 
        as.data.frame()

ggplot(as.data.frame(subdf),
       aes(y = Freq,
           axis1 = Multi_omic_subtype,
           axis2 = PAM50,
           axis3 = Immunophenotype, 
           axis4 = ImmuneSubtype))+
        scale_fill_manual(values = mycol) + 
        ggalluvial::geom_flow(stat = "alluvium",width = 1/8,aes(fill = Multi_omic_subtype)) +
        scale_y_reverse()+
        
        coord_flip() + 
        
        geom_stratum(width = 1/8, reverse = T) +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)),
                  reverse = T) +
        scale_x_continuous(breaks = 1:4, labels = c("Multi_omic_subtype","PAM50","Immunophenotype", "ImmuneSubtype")) +
        theme(legend.position = "bottom", 
              legend.title = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.text.x = element_text(size = 12, face = "bold", color = "black")) +
        ggtitle("")
ggsave("sankey update.pdf",height   = 10, width    = 4)


# load R package
library(GSVA)
library(ComplexHeatmap)
library(gplots)
library(clusterProfiler)
library(tidyr)
library(ggplot2)
library(estimate)
library(ggpubr)
source("twoclasslimma.R")
library("cogena")


# load meta signature
meta.sig <- gmt2list("metabolism signatures.gmt")
meta.sig <- sapply(meta.sig, function(x) setdiff(x,""))
meta.class <- NULL
for (i in names(meta.sig)) {
        tmp <- meta.sig[i]
        for (j in tmp) {
                meta.class <- rbind.data.frame(meta.class,
                                               data.frame(gene = j,
                                                          path = i,
                                                          stringsAsFactors = F),
                                               stringsAsFactors = F)
        }
}
expr = as.matrix(log2(mrna.tpm + 1))
write.table(expr,"KIRC 486.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# calculate GSVA enrichment score
meta.score <- gsva(expr = as.matrix(log2(mrna.tpm + 1)),
                   gset.idx.list = meta.sig,
                   method = "gsva")
write.table(meta.score,"gsva enrichment score of metabolism signature in tcga cohort.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
meta.score<-t(meta.score)
# differential analysis
## CS1 vs others
cluster<-read.table("PCa molecular subtype.txt",header=T,sep="\t",row.names=1,check.names=F)

tmp <- test1

colnames(tmp)<-c("ID","Multi_omic_subtype")

tmp$Multi_omic_subtype <- ifelse(tmp$Multi_omic_subtype == "1","CS1","Others")
subt <- data.frame(condition = tmp$Multi_omic_subtype,
                   row.names = rownames(tmp))
twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
              featmat  = meta.score[,rownames(subt)], # expression file (fill detect data scale automatically)
              treatVar = "CS1", # name of treatment group
              ctrlVar  = "Others", # name of control group
              prefix   = "tcga_metacluster", # prefix of the DE file
              overwt   = TRUE, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE, # if showing verbose result
              res.path = ) # path for result


maf <- read_tsv("data_mutations_extended.txt", comment = "#")

label <- c("Tumor_Sample_Barcode",
           "Hugo_Symbol",
           "Chromosome",
           "Start_Position",
           "End_Position",
           "Variant_Classification",
           "Variant_Type",
           "Reference_Allele",
           "Tumor_Seq_Allele1",
           "Tumor_Seq_Allele2")
maf <- maf[,label]

## CS2 vs Others
tmp <- test1;colnames(tmp)<-c("ID","Multi_omic_subtype")
tmp$Multi_omic_subtype <- ifelse(tmp$Multi_omic_subtype == "2","CS2","Others")
subt <- data.frame(condition = tmp$Multi_omic_subtype,
                   row.names = rownames(tmp))
twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
              featmat  = meta.score[,rownames(subt)], # expression file (fill detect data scale automatically)
              treatVar = "CS2", # name of treatment group
              ctrlVar  = "Others", # name of control group
              prefix   = "tcga_Multi_omic_subtype", # prefix of the DE file
              overwt   = TRUE, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE, # if showing verbose result
              res.path =) # path for result

## CS3 vs Others
tmp <- test1;colnames(tmp)<-c("ID","Multi_omic_subtype")
tmp$Multi_omic_subtype <- ifelse(tmp$Multi_omic_subtype == "3","CS3","Others")
subt <- data.frame(condition = tmp$Multi_omic_subtype,
                   row.names = rownames(tmp))
twoclasslimma(subtype  = subt, # subtype information (must contain a column named 'condition')
              featmat  = meta.score[,rownames(subt)], # expression file (fill detect data scale automatically)
              treatVar = "CS3", # name of treatment group
              ctrlVar  = "Others", # name of control group
              prefix   = "tcga_Multi_omic_subtype", # prefix of the DE file
              overwt   = TRUE, # whether over write files that already existed
              sort.p   = TRUE, # if sorting the results by adjusted p value
              verbose  = TRUE, # if showing verbose result
              res.path = ) # path for result

# extract group specific pathways
tmp1 <- read.table("tcga_Multi_omic_subtype_limma_test_result.CS1_vs_Others.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tmp2 <- read.table("tcga_Multi_omic_subtype_limma_test_result.CS2_vs_Others.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
tmp3 <- read.table("tcga_Multi_omic_subtype_limma_test_result.CS3_vs_Others.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

dep1 <- tmp1[which(tmp1$log2fc > 0.19 & tmp1$padj < 0.01),]
dep2 <- tmp2[which(tmp2$log2fc > 0.15 & tmp2$padj < 0.01),]
dep3 <- tmp3[which(tmp3$log2fc > 0.24 & tmp3$padj < 0.01),]

dep <- c(rownames(dep1),rownames(dep2),rownames(dep3))
dep <- setdiff(dep,dep[duplicated(dep)])

dep1 <- intersect(dep,rownames(dep1))
dep2 <- intersect(dep,rownames(dep2))
dep3 <- intersect(dep,rownames(dep3))

# generate metabolism pathway heatmap

Multi_omic_subtype <- test1;colnames(Multi_omic_subtype)<-c("ID","Multi_omic_subtype")


annCol <- Multi_omic_subtype
annCol <- Multi_omic_subtype[order(Multi_omic_subtype$Multi_omic_subtype),]
annRow <- data.frame(path = rep(c("CS1-specific","CS2-specific","CS3-specific"),c(length(dep1),length(dep2),length(dep3))),
                     row.names = c(dep1,dep2,dep3))

annColors <- list()
annColors[["Multi_omic_subtype"]] <- c("1" = "#2EC4B6", "2" = "#E71D36", "3" = "#FF9F1C")
annColors[["path"]] <- c("CS1-specific" = "#2EC4B6", "CS2-specific" = "#E71D36", "CS3-specific" = "#FF9F1C")

indata <- meta.score[c(dep1,dep2,dep3),rownames(annCol)]

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
        outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
        if (!is.null(halfwidth)) {
                outdata[outdata>halfwidth]=halfwidth
                outdata[outdata<(-halfwidth)]= -halfwidth
        }
        return(outdata)
}

indata <- standarize.fun(indata, halfwidth = 1)
heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
pdf("metabolism genesets.pdf",width = 10,height = 4)
pheatmap(indata,
         border_color = NA,
         color = colorRampPalette(c("#31B29E", "black","#EB6C5A"))(50),
         annotation_col = annCol[,"Multi_omic_subtype",drop = F],
         #annotation_row = annRow,
         annotation_colors = annColors,
         show_rownames = T,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         gaps_row = cumsum(table(annRow$path))[1:2])
invisible(dev.off())

dev.off()




cluster<-read.table("PCa molecular subtype.txt",header=T,sep="\t",row.names=1,check.names=F)


library(survival)
library("survminer")

rt=cluster
rt$futime=rt$futime/30.5                                       #å¦‚æžœä»¥æœˆä¸ºå•ä½ï¼Œé™¤ä»¥30ï¼›ä»¥å¹´ä¸ºå•ä½ï¼Œé™¤ä»?365
outTab=data.frame()


fit <- survfit(Surv(futime, fustat) ~ PAM50, data = rt)

pdf(file="PAM50.pdf", width=4.1,height=5,onefile = FALSE)

ggsurvplot(fit,
           pval = TRUE,  pval.coord = c(0, 0.15), conf.int = F,
           conf.int.style = "step",#ç½®ä¿¡åŒºé—´çš„ç±»åž‹ï¼Œè¿˜å¯æ”¹ä¸ºribbon,step
           censor = T, #ä¸æ˜¾ç¤ºè§‚å¯Ÿå€¼æ‰€åœ¨çš„ä½ç½®
           risk.table = T, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "solid", # Change line type by groups
           ggtheme = theme_classic(), # Change ggplot2 theme
           palette = c("#FFD966", "#9BC2E6","#7030A0"), #
           ylim = c(0,1),
           font.legend = 12,
           font.main = c(14, "bold", "darkblue"),
           font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"),
           font.tickslab = c(12, "plain", "black"))

dev.off()

LumB<- rownames(cluster[which(cluster$PAM50 %in% "LumB"),])

LumB.surv <- cluster[which(cluster$PAM50 %in% "LumB"),]

expr<-t(expr)

testx<-intersect(LumB,rownames(expr))

LumB.expr <- t(expr[testx,])

LumB.group<-LumB.surv[testx,]

MMRs<-c("MLH1","MSH2","MSH6","PMS2","EPCAM")

MMRs.expr<-LumB.expr[MMRs,]


write.table(LumB.expr,file="LumB.expr.txt",sep = "\t",row.names = T,quote = F)

write.table(LumB.group,file="LumB.group.txt",sep = "\t",row.names = T,quote = F)


rt2=as.data.frame(LumB.surv)
rt2$futime=rt2$futime/30.5                                   #å¦‚æžœä»¥æœˆä¸ºå•ä½ï¼Œé™¤ä»¥30ï¼›ä»¥å¹´ä¸ºå•ä½ï¼Œé™¤ä»?365
outTab=data.frame()

rt2$Multi_omic_subtype<-ifelse(rt2$Multi_omic_subtype=="CS2","CS2","CS1+CS3")

fit <- survfit(Surv(futime, fustat) ~ Multi_omic_subtype, data = rt2)

pdf(file="LumB.pdf", width=4.8,height=5,onefile = FALSE)

ggsurvplot(fit,
           pval = TRUE,  pval.coord = c(0, 0.5), conf.int = F,
           conf.int.style = "step",#ç½®ä¿¡åŒºé—´çš„ç±»åž‹ï¼Œè¿˜å¯æ”¹ä¸ºribbon,step
           censor = T, #ä¸æ˜¾ç¤ºè§‚å¯Ÿå€¼æ‰€åœ¨çš„ä½ç½®
           risk.table = T, # Add risk table
           risk.table.col = "strata", # Change risk table color by groups
           linetype = "solid", # Change line type by groups
           ggtheme = theme_classic(), # Change ggplot2 theme
           palette = c("#5DAC54", "#E71D36","#FF9F1C"), #
           ylim = c(0,1),
           font.legend = 12,
           font.main = c(14, "bold", "darkblue"),
           font.x = c(14, "bold", "black"),
           font.y = c(14, "bold", "black"),
           font.tickslab = c(12, "plain", "black"))

dev.off()


#-----------MMRs heatmap


MMRs<-c("MLH1","MSH2","MSH6","PMS2","EPCAM")

MMRs.expr<-LumB.expr[MMRs,]

MMRsmerge<-intersect(rownames(LumB.group),colnames(MMRs.expr))

LumB.group<-LumB.group[MMRsmerge,]
MMRs.expr<-t(MMRs.expr[,MMRsmerge])

MMRstotal<-cbind(LumB.group,MMRs.expr)
MMRstotal$Multi_omic_subtype<-ifelse(MMRstotal$Multi_omic_subtype=="CS2","CS2","CS1+CS3")
MMRstotal<-MMRstotal[order(MMRstotal$Multi_omic_subtype,decreasing = F),] 

rt<-log2(MMRstotal[,5:9]+1)

rt2 <- scale(t(rt)) #scaleæ ‡åŒ–

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
        outdata=t(scale(t(rt2), center=centerFlag, scale=scaleFlag))
        if (!is.null(halfwidth)) {
                outdata[outdata>halfwidth]=halfwidth
                outdata[outdata<(-halfwidth)]= -halfwidth
        }
        return(outdata)
}

rt3 <- standarize.fun(rt2,halfwidth =1)

head(rt3[1:3,1:3])

#è¯»å–åˆ†ç»„æ–‡ä»¶
rt4 <- MMRstotal[,1:4]

cluster=rt4
library("pheatmap")


annColors <- list(Multi_omic_subtype = c("CS1+CS3"="#5DAC54",
                                         "CS2"="#E71D36"),
                  PAM50 = c("LumB"="#7030A0"),
                  Immunophenotype = c("Non-Immune"="#D8D8D8",
                                      "Suppressed"="#70AD47",
                                      "Activated"="#FF0000"),
                  ImmuneSubtype = c("C1"="#E71D36",
                                    "C2"="yellow",
                                    "C3"="#00B050",
                                    "C4"="#00B0F0",
                                    "Unknown"="#D8D8D8")
)

#ç»˜åˆ¶çƒ­å›¾
pdf("LumB-KIRC_ClustOutPut.pdf",height=4,width=9)
p_heatmap<-pheatmap(rt3, annotation=cluster, 
                    color = colorRampPalette(c("#31B29E", "black","#EB6C5A"))(50),
                    cluster_cols =F,
                    fontsize=8,
                    fontsize_row=8,
                    scale="row",
                    annotation_colors=annColors,
                    show_colnames=F,
                    fontsize_col=3)
dev.off()

table(cluster$Multi_omic_subtype)


library(ggplot2)
library(ggpubr)
conNum=114                                                      #normalç»„æ ·å“æ•°ç›?
treatNum=106                                         #tumorç»„æ ·å“æ•°ç›?
rt=rt

type=c(rep("CS1+CS3",conNum),rep("CS2",treatNum))

#å‡†å¤‡ç®±çº¿å›¾çš„è¾“å…¥æ–‡ä»¶
data=data.frame()
for(i in colnames(rt)){
        data=rbind(data,cbind(expression=rt[,i],gene=i,type))
}
write.table(data,file="data.txt",sep="\t",row.names=F,quote=F)

#ç»˜åˆ¶ç®±åž‹å›?
data=read.table("data.txt",sep="\t",header=T,check.names=F)  
#ç»˜åˆ¶ç®±åž‹å›?
p=ggboxplot(data, x="gene", y="expression", notch = F,
            lwd=0.5, alpha = 0.6,color = "type", fill = 'type', 
            ylab="Immunocyte Infiltration",
            xlab="",
            palette = c( "#4EBC97","#E7B800") )
p=p+rotate_x_text(60)
pdf(file="boxplot TCGA-KIRC PCIPI.pdf",width=10,height=6)                          #è¾“å‡ºå›¾ç‰‡æ–‡ä»¶  #t.test
p+stat_compare_means(aes(group=type),method = "wilcox.test", symnum.args=list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")),label = "p.signif")
dev.off()

#---------------PVIT 
PVIT<-read.table("PVITCNA.txt",header=T,sep="\t",row.names=1,check.names=F)
data<-PVIT

data$number <- 1
library(ggplot2)
library(plyr)

test.data<-table(data$cluster,data$CAN)
test.data
chisq.test(test.data)
chisq.test(test.data)$expected
ax<-chisq.test(test.data)$p.value
ax
pvalue= paste(pvalue= ifelse (ax <0.001, " P < 0.001", paste ("P = ",round(ax,3),sep = "")))



test.data2<-data.frame(test.data)
test.data2

test.data2<-test.data2[order(test.data2[,2]),] 

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
        scale_fill_brewer(palette="PRGn")

p

p2<-p+geom_text(aes(x=2,y=0.8,label= pvalue), vjust=0.5, hjust=0.5, size=6, fontface="italic",color="black")
p2

pdf("PVIT CNA.pdf",width=5.8,height=4,onefile = FALSE)
p2
dev.off()

PVIT3<-PVIT;PVIT3$PVTI<-log2(PVIT3$PVTI+1)
library(ggpubr)
library("ggplot2")
p<-ggplot(PVIT3, aes(x=CAN, y=PVTI, fill=CAN)) +
        geom_boxplot()+scale_x_discrete(limits=c("4Deletion", "3Diploid","2Gain","1Amplification"))+
        stat_compare_means()
p


ggboxplot(PVIT3, x = "CAN", y = "PVTI",
          fill = "CAN", palette = "jco")+
        stat_compare_means()
ggsave("PVIT CNA expr.pdf",width=3.5,height=4.5,onefile = FALSE)

PVIT4<-PVIT3
PVIT4$CAN<-ifelse(PVIT4$CAN=="4Deletion","3Diploid",PVIT4$CAN)
#PVIT3$CAN<-ifelse(PVIT3$CAN=="2Gain","1Amplification",PVIT3$CAN)

library("survminer")
require("survival")
fit <- survfit(Surv(futime, fustat) ~ CAN, data =PVIT4)
ggsurvplot(fit, data = PVIT4,
           pval = TRUE,
           palette =c("#7D3397", "#C2A5CF","#A6DBA0"))

ggsave("PVIT CNA surv2.pdf",width=3.5,height=4.5,onefile = FALSE)



#---------------10åˆ†ç»„ ç»Ÿä¸€å±•ç¤º

t1<-moic.res.list$LRAcluster$clust.res$clust
t2<-moic.res.list$SNF$clust.res$clust
t3<-moic.res.list$PINSPlus$clust.res$clust
t4<-moic.res.list$NEMO$clust.res$clust
t5<-moic.res.list$COCA$clust.res$clust
t6<-moic.res.list$ConsensusClustering$clust.res$clust
t7<-moic.res.list$IntNMF$clust.res$clust
t8<-moic.res.list$CIMLR$clust.res$clust
t9<-moic.res.list$MoCluster$clust.res$clust
t10<-moic.res.list$iClusterBayes$clust.res$clust
ttotal<-rbind(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,cmoic.KIRC$clust.res$clust)
rownames(ttotal)<-c("LRA","SNF","PINSPlus","NEMO","COCA","ConsensusClustering","IntNMF","CIMLR","MoCluster","iClusterBayes","Cmoic")
colnames(ttotal)<-moic.res.list$LRAcluster$clust.res$samID

write.table(ttotal,file="10methodcluster.txt",sep="\t",row.names=T,quote=F)

library(dplyr)
library(pheatmap)  

rt3<-lncrna.tpm[1:10,moic.res.list$LRAcluster$clust.res$samID]

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
        outdata=t(scale(t(rt3), center=centerFlag, scale=scaleFlag))
        if (!is.null(halfwidth)) {
                outdata[outdata>halfwidth]=halfwidth
                outdata[outdata<(-halfwidth)]= -halfwidth
        }
        return(outdata)
}

rt3 <- standarize.fun(rt3,halfwidth =0.8)
rt3[1:10,1:5]

rt4<-ifelse(ttotal=="1","1",ifelse(ttotal=="2","2","3"))
rt4<-as.data.frame(t(ttotal))

rtx<-

        pdf("ClustOutPut.pdf",height=3,width=10)
p_heatmap<-pheatmap(rt3, annotation=rt4, 
                    color = colorRampPalette(c("#0099CC", "white","#CC0033"))(50),
                    cluster_cols =F,
                    fontsize=8,
                    fontsize_row=8,
                    scale="row",
                    show_colnames=F,
                    fontsize_col=3)
dev.off()

p_heatmap

pheatmap(rt3,annotation=rt4)


annColors <- list(Multi_omic_subtype = c("CS1"="#2EC4B6",
                                         "CS2"="#E71D36",
                                         "CS3"="#FF9F1C"),
                  PAM50 = c("LumA"="#9BCS2E6",
                            "LumB"="#7030A0",
                            "Basal"="#FFD966"),
                  Immunophenotype = c("Non-Immune"="#D8D8D8",
                                      "Suppressed"="#70AD47",
                                      "Activated"="#FF0000"),
                  ImmuneSubtype = c("CS1"="#E71D36",
                                    "CS2"="yellow",
                                    "CS3"="#00B050",
                                    "C4"="#00B0F0",
                                    "Unknown"="#D8D8D8")
)

#ç»˜åˆ¶çƒ­å›¾
pdf("Comparesig\\multi-PRAD_ClustOutPut.pdf",height=4,width=9)
p_heatmap<-pheatmap(rt3, annotation=rt4, 
                    color = colorRampPalette(c("#31B29E", "black","#EB6C5A"))(50),
                    cluster_cols =F,
                    fontsize=8,
                    fontsize_row=8,
                    scale="row",
                    #annotation_colors=annColors,
                    show_colnames=F,
                    fontsize_col=3)
dev.off()



#-------------------------library(ggplot2)

##Ò»¸öÕ¹Ê¾»ùÒò²îÒì±í´ïË®Æ½µÄ»ðÉ½Í¼Ê¾Àý
gene <- read.delim('SETD2 export.txt', sep = '\t', stringsAsFactors = FALSE)
head(gene)
#ÀýÈçÕâÀï×Ô¶¨Òå¸ù¾Ý |log2FC| >= 1 ºÍ adj.P.Val < 0.01 ±ê¼Ç²îÒìÀàÐÍ
gene[which(gene$P.value < 0.05 & gene$Effect.size <= 0),'sig'] <- 'Sensitive'
#gene[which(gene$P.value < 0.05 & gene$Effect.size >= 0),'sig'] <- 'Up'
gene[which(gene$P.value >= 0.05 | abs(gene$Effect.size) < 0.3),'sig'] <- 'None'

#ºáÖá log2FC£¬×ÝÖá -log10(adj.P.Val)£¬ÑÕÉ«±íÊ¾²îÒì
p <- ggplot(gene, aes(x = Effect.size, y = -log10(P.value), color = sig)) +
        geom_point(alpha = 1.2, size = 3) +
        scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c('Sensitive', 'Down', 'None')) +
        theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
        theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
        geom_vline(xintercept = c(0), color = 'gray', size = 0.3) +
        geom_hline(yintercept = -log(0.05, 10), color = 'gray', size = 0.3) +
        xlim(-1, 1) + ylim(0, 3) +
        labs(x = '\nIC50 Effect', y = 'log10(p-value)\n', color = '', title = 'Drugs for SETD2 mut\n')

p

ggsave('Drugs for SETD2.pdf', p, width = 4, height = 4)

