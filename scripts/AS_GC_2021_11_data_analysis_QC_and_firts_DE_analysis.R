library(tidyverse)
library(DESeq2)
library(stringr)
library("RColorBrewer")
library("pheatmap")

getwd()
primeDirectory = paste0("/Users/anita/Box Sync/BulkRNAseq/GC006/")
setwd(primeDirectory)
getwd()



GC.sample.info <- read_csv("sampleinfo_GC.csv")
GC.sample.info

# reading summary and creating QC graphs
pre.OC.GC <- read.table("GC.counts.txt.summary", header = TRUE)
QC.table.GC <- as.data.frame(t(pre.OC.GC[,-1]))
colnames(QC.table.GC) <- pre.OC.GC$Status
rm(pre.OC.GC)
QC.table.GC$sample <- str_remove(row.names(QC.table.GC),".star.Aligned.sortedByCoord.out.bam")
QC.table.GC$sample <- str_remove(QC.table.GC$sample,"X.project.Owens_Rivanna.00.Raw.and.Aligned.Files.bulk.RNA.seq.2021.11.in.vitro.GC.03.align.")
QC.table.GC <- select(QC.table.GC, sample, Assigned:Unassigned_NoFeatures)
QC.table.GC <- inner_join(GC.sample.info, QC.table.GC, by="sample")
QC.table.GC <- mutate(QC.table.GC, counts = Assigned+Unassigned_MultiMapping+Unassigned_NoFeatures)
QC.table.GC <- mutate(QC.table.GC, percent_assigned = (Assigned/counts)*100)
QC.table.GC <- mutate(QC.table.GC, percent_unassigned = 100-percent_assigned)

#Plotting graphs and saving them
ggplot(QC.table.GC, aes(condition,percent_assigned)) + geom_bar(stat="identity", aes(fill=sample), position="dodge", show.legend = FALSE) + labs(y="Percent of assigned reads", x="") + theme(axis.text.x = element_text(angle=45, hjust = 1))
tiff(filename = "2021.11.QC.Percent.Assigned.by.Condition.tiff", width = 1920, height = 1080)
ggplot(QC.table.GC, aes(condition,percent_assigned)) + geom_bar(stat="identity", aes(fill=sample), position="dodge", show.legend = FALSE) + 
  labs(y="Percent of assigned reads", x="Condition") + 
  theme(axis.text.x = element_text(angle=45, hjust = 1, size = "24")) + 
  theme(axis.text.y = element_text(size = 24)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30))
dev.off()

ggplot(QC.table.GC, aes(condition,percent_unassigned)) + geom_bar(stat="identity", aes(fill=sample), position="dodge", show.legend = FALSE) + labs(y="Percent of unassigned reads", x="") + theme(axis.text.x = element_text(angle=45, hjust = 1))
tiff(filename = "2021.11.QC.Percent.Unassigned.by.Condition.tiff", width = 1920, height = 1080)
ggplot(QC.table.GC, aes(condition,percent_unassigned)) + geom_bar(stat="identity", aes(fill=sample), position="dodge", show.legend = FALSE) +
  labs(y="Percent of unassigned reads", x="Condition") +
  theme(axis.text.x = element_text(angle=45, hjust = 1, size = "24")) + 
  theme(axis.text.y = element_text(size = 24)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30))
dev.off()


# Import & pre-process ----------------------------------------------------

## Import data from featureCounts
## featureCounts -a genes.gtf -o counts.txt -T 12 -t exon -g gene_name *sam
countGCdata <- read.table("GC.counts.txt", header=TRUE, row.names=1, sep="\t")
colnames(countGCdata)

## Remove first five columns (chr, start, end, strand, length)
countGCdata <- countGCdata[ ,6:ncol(countGCdata)]

## Remove cruft from filenames
colnames(countGCdata) <- colnames(countGCdata) %>% str_remove(".star.Aligned.sortedByCoord.out.bam")
colnames(countGCdata) <- colnames(countGCdata) %>% str_remove("X.project.Owens_Rivanna.00.Raw.and.Aligned.Files.bulk.RNA.seq.2021.11.in.vitro.GC.03.align.")

colnames(countGCdata)

## Convert to matrix (mot sure if this is still necessary in the newer version, but it does not hurt)
countGCdata <- as.matrix(countGCdata)
head(countGCdata)

## Get colASE081data
colGCdata <- data.frame(row.names=GC.sample.info$sample, GC.sample.info %>% select(-sample))
# Convert character columns to factor
(colGCdata$condition <- factor(colGCdata$condition))
colGCdata


# Check to see if the data is set properly. It is absolutely necessary to have the columns in the count data to be in the same order as the row names
# in col data. First we check if all of them are present and then if they are in the same order
all(rownames(colGCdata) %in% colnames(countGCdata))

all(rownames(colGCdata) == colnames(countGCdata))

# Now, we create the DESeq2 data frame using the info that we previously arranged
ddsGC <- DESeqDataSetFromMatrix(countData=countGCdata, 
                                       colData=colGCdata, 
                                       design=~condition) # if we have more info, you can add here, such as different genotypes, time points, etc
ddsGC

rldGC <- rlogTransformation(ddsGC)
plotPCA(rldGC)

  tiff(filename = "2021.11.PCA.GC1.tiff", width = 1920, height = 1080)
  plotPCA(rldGC)+
    geom_point(size = 10, aes(shape=group)) +
    scale_shape_manual(values=1:nlevels(rldGC$condition)) +
    theme(axis.text.x = element_text(size = 24)) +
    theme(axis.text.y = element_text(size = 24)) +
    theme(axis.title.x = element_text(size = 30)) +
    theme(axis.title.y = element_text(size = 30)) +
    theme(legend.title = element_text(size=30)) +
    theme(legend.text = element_text(size=30))
  dev.off()
#######
sampleDists <- dist(t(assay(rldGC)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rldGC$condition
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "BuGn")) )(255)

tiff(filename = "2021.11.Sample.Distance.Heatmap.tiff", width = 1920, height = 1080)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, fontsize = 16)
dev.off()

# Start the Differential analysis
ddsGC <- DESeq(ddsGC)


genename <-"Spp1"
#singlegenes analysis
t <- plotCounts(ddsASE081, gene=genename, intgroup="conditions", returnData = TRUE)
#t %>% View() 
tiff(filename = paste0(genename), width = 400, height = 290) 
ggplot(t, aes(x = conditions, y = count, color = conditions)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0), size=5) +
  #geom_text_repel(aes(label = rownames(t))) + 
  theme_bw()+
  theme(axis.text.x =element_blank())+
  ggtitle(genename) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position= "none")+
  #theme(axis.text.x = element_text(angle=70, hjust = 1))+
  labs(y="normalized count", x="Condition")
dev.off()
?ggplot
rm(t)

sum(ddsASE081)










rldASE0811 <- rlogTransformation(ddsASE081)
plotPCA(rldASE0811)

tiff(filename = "2020.01.PCA.All.Samples1.tiff", width = 1920, height = 1080)
plotPCA(rldASE0811) +
  geom_point(size = 10, aes(shape=group)) +
  scale_shape_manual(values=1:nlevels(rldASE081$condition)) +
  theme(axis.text.x = element_text(size = 24)) +
  theme(axis.text.y = element_text(size = 24)) +
  theme(axis.title.x = element_text(size = 30)) +
  theme(axis.title.y = element_text(size = 30)) +
  theme(legend.title = element_text(size=30)) +
  theme(legend.text = element_text(size=30))
dev.off()

sampleDists1 <- dist(t(assay(rldASE081)))

sampleDistMatrix1 <- as.matrix(sampleDists1)
rownames(sampleDistMatrix1) <- rldASE0811$condition
colnames(sampleDistMatrix1) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "BuGn")) )(255)

tiff(filename = "2020.01.Sample.Distance.Heatmap1.tiff", width = 1920, height = 1080)
pheatmap(sampleDistMatrix1,
         clustering_distance_rows=sampleDists1,
         clustering_distance_cols=sampleDists1,
         col=colors, fontsize = 16)
dev.off()

rownames(colASE081data)

# Different contrasts
results.chol.alone.vs.vehicle <- results(ddsASE081, contrast = c("condition", "Cholesterol","Vehicle"))
sum(results.chol.alone.vs.vehicle$padj < 0.05, na.rm=TRUE)
results.chol.alone.vs.vehicle
results.chol.alone.vs.vehicle <- as.data.frame(results.chol.alone.vs.vehicle)
results.chol.alone.vs.vehicle$gene <- rownames(results.chol.alone.vs.vehicle)
results.chol.alone.vs.vehicle <- select(results.chol.alone.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(results.chol.alone.vs.vehicle, "2020.01.results.chol.alone.vs.vehicle.csv")
sig_results.chol.alone.vs.vehicle <- subset(results.chol.alone.vs.vehicle, padj <0.05)
sig_results.chol.alone.vs.vehicle <- as.data.frame(sig_results.chol.alone.vs.vehicle)
sig_results.chol.alone.vs.vehicle$gene <- rownames(sig_results.chol.alone.vs.vehicle)
sig_results.chol.alone.vs.vehicle <- select(sig_results.chol.alone.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.chol.alone.vs.vehicle, "2020.01.significant-altered-genes_results.chol.alone.vs.vehicle.csv")

results.chol.PDGF.TGFB1.vs.vehicle <- results(ddsASE081, contrast = c("condition", "Cholesterol+PDGF+TGFB1","Vehicle"))
sum(results.chol.PDGF.TGFB1.vs.vehicle$padj < 0.05, na.rm=TRUE)
results.chol.PDGF.TGFB1.vs.vehicle
results.chol.PDGF.TGFB1.vs.vehicle <- as.data.frame(results.chol.PDGF.TGFB1.vs.vehicle)
results.chol.PDGF.TGFB1.vs.vehicle$gene <- rownames(results.chol.PDGF.TGFB1.vs.vehicle)
results.chol.PDGF.TGFB1.vs.vehicle <- select(results.chol.PDGF.TGFB1.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(results.chol.PDGF.TGFB1.vs.vehicle, "2020.01.results.chol.PDGF.TGFB1.vs.vehicle.csv")
sig_results.chol.PDGF.TGFB1.vs.vehicle <- subset(results.chol.PDGF.TGFB1.vs.vehicle, padj <0.05)
sig_results.chol.PDGF.TGFB1.vs.vehicle <- as.data.frame(sig_results.chol.PDGF.TGFB1.vs.vehicle)
sig_results.chol.PDGF.TGFB1.vs.vehicle$gene <- rownames(sig_results.chol.PDGF.TGFB1.vs.vehicle)
sig_results.chol.PDGF.TGFB1.vs.vehicle <- select(sig_results.chol.PDGF.TGFB1.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.chol.PDGF.TGFB1.vs.vehicle, "2020.01.significant-altered-genes_results.chol.PDGF.TGFB1.vs.vehicle.csv")

results.CPI613.alone.vs.vehicle <- results(ddsASE081, contrast = c("condition", "CPI613","Vehicle"))
sum(results.CPI613.alone.vs.vehicle$padj < 0.05, na.rm=TRUE)
results.CPI613.alone.vs.vehicle
results.CPI613.alone.vs.vehicle <- as.data.frame(results.CPI613.alone.vs.vehicle)
results.CPI613.alone.vs.vehicle$gene <- rownames(results.CPI613.alone.vs.vehicle)
results.CPI613.alone.vs.vehicle <- select(results.CPI613.alone.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(results.CPI613.alone.vs.vehicle, "2020.01.results.CPI613.alone.vs.vehicle.csv")
sig_results.CPI613.alone.vs.vehicle <- subset(results.CPI613.alone.vs.vehicle, padj <0.05)
sig_results.CPI613.alone.vs.vehicle <- as.data.frame(sig_results.CPI613.alone.vs.vehicle)
sig_results.CPI613.alone.vs.vehicle$gene <- rownames(sig_results.CPI613.alone.vs.vehicle)
sig_results.CPI613.alone.vs.vehicle <- select(sig_results.CPI613.alone.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.CPI613.alone.vs.vehicle, "2020.01.significant-altered-genes_results.CPI613.alone.vs.vehicle.csv")

results.Galloflavin.alone.vs.vehicle <- results(ddsASE081, contrast = c("condition", "Galloflavin","Vehicle"))
sum(results.Galloflavin.alone.vs.vehicle$padj < 0.05, na.rm=TRUE)
results.Galloflavin.alone.vs.vehicle
results.Galloflavin.alone.vs.vehicle <- as.data.frame(results.Galloflavin.alone.vs.vehicle)
results.Galloflavin.alone.vs.vehicle$gene <- rownames(results.Galloflavin.alone.vs.vehicle)
results.Galloflavin.alone.vs.vehicle <- select(results.Galloflavin.alone.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(results.Galloflavin.alone.vs.vehicle, "2020.01.results.Galloflavin.alone.vs.vehicle.csv")
sig_results.Galloflavin.alone.vs.vehicle <- subset(results.Galloflavin.alone.vs.vehicle, padj <0.05)
sig_results.Galloflavin.alone.vs.vehicle <- as.data.frame(sig_results.Galloflavin.alone.vs.vehicle)
sig_results.Galloflavin.alone.vs.vehicle$gene <- rownames(sig_results.Galloflavin.alone.vs.vehicle)
sig_results.Galloflavin.alone.vs.vehicle <- select(sig_results.Galloflavin.alone.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.Galloflavin.alone.vs.vehicle, "2020.01.significant-altered-genes_results.Galloflavin.alone.vs.vehicle.csv")

results.PDGF.alone.vs.vehicle <- results(ddsASE081, contrast = c("condition", "PDGF-BB","Vehicle"))
sum(results.PDGF.alone.vs.vehicle$padj < 0.05, na.rm=TRUE)
results.PDGF.alone.vs.vehicle
results.PDGF.alone.vs.vehicle <- as.data.frame(results.PDGF.alone.vs.vehicle)
results.PDGF.alone.vs.vehicle$gene <- rownames(results.PDGF.alone.vs.vehicle)
results.PDGF.alone.vs.vehicle <- select(results.PDGF.alone.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(results.PDGF.alone.vs.vehicle, "2020.01.results.PDGF.alone.vs.vehicle.csv")
sig_results.PDGF.alone.vs.vehicle <- subset(results.PDGF.alone.vs.vehicle, padj <0.05)
sig_results.PDGF.alone.vs.vehicle <- as.data.frame(sig_results.PDGF.alone.vs.vehicle)
sig_results.PDGF.alone.vs.vehicle$gene <- rownames(sig_results.PDGF.alone.vs.vehicle)
sig_results.PDGF.alone.vs.vehicle <- select(sig_results.PDGF.alone.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.PDGF.alone.vs.vehicle, "2020.01.significant-altered-genes_results.PDGF.alone.vs.vehicle.csv")

results.PDGF.TGFB1.vs.vehicle <- results(ddsASE081, contrast = c("condition", "PDGF+TGFB1","Vehicle"))
sum(results.PDGF.TGFB1.vs.vehicle$padj < 0.05, na.rm=TRUE)
results.PDGF.TGFB1.vs.vehicle
results.PDGF.TGFB1.vs.vehicle <- as.data.frame(results.PDGF.TGFB1.vs.vehicle)
results.PDGF.TGFB1.vs.vehicle$gene <- rownames(results.PDGF.TGFB1.vs.vehicle)
results.PDGF.TGFB1.vs.vehicle <- select(results.PDGF.TGFB1.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(results.PDGF.TGFB1.vs.vehicle, "2020.01.results.PDGF.TGFB1.vs.vehicle.csv")
sig_results.PDGF.TGFB1.vs.vehicle <- subset(results.PDGF.TGFB1.vs.vehicle, padj <0.05)
sig_results.PDGF.TGFB1.vs.vehicle <- as.data.frame(sig_results.PDGF.TGFB1.vs.vehicle)
sig_results.PDGF.TGFB1.vs.vehicle$gene <- rownames(sig_results.PDGF.TGFB1.vs.vehicle)
sig_results.PDGF.TGFB1.vs.vehicle <- select(sig_results.PDGF.TGFB1.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.PDGF.TGFB1.vs.vehicle, "2020.01.significant-altered-genes_results.PDGF.TGFB1.vs.vehicle.csv")

results.PDGF.TGFB1.CPI613.vs.vehicle <- results(ddsASE081, contrast = c("condition", "PDGF+TGFB1+CPI613","Vehicle"))
sum(results.PDGF.TGFB1.CPI613.vs.vehicle$padj < 0.05, na.rm=TRUE)
results.PDGF.TGFB1.CPI613.vs.vehicle
results.PDGF.TGFB1.CPI613.vs.vehicle <- as.data.frame(results.PDGF.TGFB1.CPI613.vs.vehicle)
results.PDGF.TGFB1.CPI613.vs.vehicle$gene <- rownames(results.PDGF.TGFB1.CPI613.vs.vehicle)
results.PDGF.TGFB1.CPI613.vs.vehicle <- select(results.PDGF.TGFB1.CPI613.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(results.PDGF.TGFB1.CPI613.vs.vehicle, "2020.01.results.PDGF.TGFB1.CPI613.vs.vehicle.csv")
sig_results.PDGF.TGFB1.CPI613.vs.vehicle <- subset(results.PDGF.TGFB1.CPI613.vs.vehicle, padj <0.05)
sig_results.PDGF.TGFB1.CPI613.vs.vehicle <- as.data.frame(sig_results.PDGF.TGFB1.CPI613.vs.vehicle)
sig_results.PDGF.TGFB1.CPI613.vs.vehicle$gene <- rownames(sig_results.PDGF.TGFB1.CPI613.vs.vehicle)
sig_results.PDGF.TGFB1.CPI613.vs.vehicle <- select(sig_results.PDGF.TGFB1.CPI613.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.PDGF.TGFB1.CPI613.vs.vehicle, "2020.01.significant-altered-genes_results.PDGF.TGFB1.CPI613.vs.vehicle.csv")

results.PDGF.TGFB1.Galloflavin.vs.vehicle <- results(ddsASE081, contrast = c("condition", "PDGF+TGFB1+Galloflavin","Vehicle"))
sum(results.PDGF.TGFB1.Galloflavin.vs.vehicle$padj < 0.05, na.rm=TRUE)
results.PDGF.TGFB1.Galloflavin.vs.vehicle
results.PDGF.TGFB1.Galloflavin.vs.vehicle <- as.data.frame(results.PDGF.TGFB1.Galloflavin.vs.vehicle)
results.PDGF.TGFB1.Galloflavin.vs.vehicle$gene <- rownames(results.PDGF.TGFB1.Galloflavin.vs.vehicle)
results.PDGF.TGFB1.Galloflavin.vs.vehicle <- select(results.PDGF.TGFB1.Galloflavin.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(results.PDGF.TGFB1.Galloflavin.vs.vehicle, "2020.01.results.PDGF.TGFB1.Galloflavin.vs.vehicle.csv")
sig_results.PDGF.TGFB1.Galloflavin.vs.vehicle <- subset(results.PDGF.TGFB1.Galloflavin.vs.vehicle, padj <0.05)
sig_results.PDGF.TGFB1.Galloflavin.vs.vehicle <- as.data.frame(sig_results.PDGF.TGFB1.Galloflavin.vs.vehicle)
sig_results.PDGF.TGFB1.Galloflavin.vs.vehicle$gene <- rownames(sig_results.PDGF.TGFB1.Galloflavin.vs.vehicle)
sig_results.PDGF.TGFB1.Galloflavin.vs.vehicle <- select(sig_results.PDGF.TGFB1.Galloflavin.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.PDGF.TGFB1.Galloflavin.vs.vehicle, "2020.01.significant-altered-genes_results.PDGF.TGFB1.Galloflavin.vs.vehicle.csv")

results.TGFB1.alone.vs.vehicle <- results(ddsASE081, contrast = c("condition", "TGFB1","Vehicle"))
sum(results.TGFB1.alone.vs.vehicle$padj < 0.05, na.rm=TRUE)
results.TGFB1.alone.vs.vehicle
results.TGFB1.alone.vs.vehicle <- as.data.frame(results.TGFB1.alone.vs.vehicle)
results.TGFB1.alone.vs.vehicle$gene <- rownames(results.TGFB1.alone.vs.vehicle)
results.TGFB1.alone.vs.vehicle <- select(results.TGFB1.alone.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(results.TGFB1.alone.vs.vehicle, "2020.01.results.TGFB1.alone.vs.vehicle.csv")
sig_results.TGFB1.alone.vs.vehicle <- subset(results.TGFB1.alone.vs.vehicle, padj <0.05)
sig_results.TGFB1.alone.vs.vehicle <- as.data.frame(sig_results.TGFB1.alone.vs.vehicle)
sig_results.TGFB1.alone.vs.vehicle$gene <- rownames(sig_results.TGFB1.alone.vs.vehicle)
sig_results.TGFB1.alone.vs.vehicle <- select(sig_results.TGFB1.alone.vs.vehicle, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.TGFB1.alone.vs.vehicle, "2020.01.significant-altered-genes_results.TGFB1.alone.vs.vehicle.csv")

### Now comparing drug and/combinations
results.chol.PDGF.TGFV1.vs.chol.alon <- results(ddsASE081, contrast = c("condition", "Cholesterol+PDGF+TGFB1","Cholesterol"))
sum(results.chol.PDGF.TGFV1.vs.chol.alon$padj < 0.05, na.rm=TRUE)
results.chol.PDGF.TGFV1.vs.chol.alon
results.chol.PDGF.TGFV1.vs.chol.alon <- as.data.frame(results.chol.PDGF.TGFV1.vs.chol.alon)
results.chol.PDGF.TGFV1.vs.chol.alon$gene <- rownames(results.chol.PDGF.TGFV1.vs.chol.alon)
results.chol.PDGF.TGFV1.vs.chol.alon <- select(results.chol.PDGF.TGFV1.vs.chol.alon, "gene", "log2FoldChange", "padj", everything())
write_csv(results.chol.PDGF.TGFV1.vs.chol.alon, "2020.01.results.chol.PDGF.TGFV1.vs.chol.alon.csv")
sig_results.chol.PDGF.TGFV1.vs.chol.alon <- subset(results.chol.PDGF.TGFV1.vs.chol.alon, padj <0.05)
sig_results.chol.PDGF.TGFV1.vs.chol.alon <- as.data.frame(sig_results.chol.PDGF.TGFV1.vs.chol.alon)
sig_results.chol.PDGF.TGFV1.vs.chol.alon$gene <- rownames(sig_results.chol.PDGF.TGFV1.vs.chol.alon)
sig_results.chol.PDGF.TGFV1.vs.chol.alon <- select(sig_results.chol.PDGF.TGFV1.vs.chol.alon, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.chol.PDGF.TGFV1.vs.chol.alon, "2020.01.significant-altered-genes_results.chol.PDGF.TGFV1.vs.chol.alon.csv")

results.PDGF.TGFB1.vs.PDGF <- results(ddsASE081, contrast = c("condition", "PDGF+TGFB1","PDGF-BB"))
sum(results.PDGF.TGFB1.vs.PDGF$padj < 0.05, na.rm=TRUE)
results.PDGF.TGFB1.vs.PDGF
results.PDGF.TGFB1.vs.PDGF <- as.data.frame(results.PDGF.TGFB1.vs.PDGF)
results.PDGF.TGFB1.vs.PDGF$gene <- rownames(results.PDGF.TGFB1.vs.PDGF)
results.PDGF.TGFB1.vs.PDGF <- select(results.PDGF.TGFB1.vs.PDGF, "gene", "log2FoldChange", "padj", everything())
write_csv(results.PDGF.TGFB1.vs.PDGF, "2020.01.results.PDGF.TGFB1.vs.PDGF.csv")
sig_results.PDGF.TGFB1.vs.PDGF <- subset(results.PDGF.TGFB1.vs.PDGF, padj <0.05)
sig_results.PDGF.TGFB1.vs.PDGF <- as.data.frame(sig_results.PDGF.TGFB1.vs.PDGF)
sig_results.PDGF.TGFB1.vs.PDGF$gene <- rownames(sig_results.PDGF.TGFB1.vs.PDGF)
sig_results.PDGF.TGFB1.vs.PDGF <- select(sig_results.PDGF.TGFB1.vs.PDGF, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.PDGF.TGFB1.vs.PDGF, "2020.01.significant-altered-genes_results.PDGF.TGFB1.vs.PDGF.csv")

results.PDGF.TGFB1.vs.TGFB1 <- results(ddsASE081, contrast = c("condition", "PDGF+TGFB1","TGFB1"))
sum(results.PDGF.TGFB1.vs.TGFB1$padj < 0.05, na.rm=TRUE)
results.PDGF.TGFB1.vs.TGFB1
results.PDGF.TGFB1.vs.TGFB1 <- as.data.frame(results.PDGF.TGFB1.vs.TGFB1)
results.PDGF.TGFB1.vs.TGFB1$gene <- rownames(results.PDGF.TGFB1.vs.TGFB1)
results.PDGF.TGFB1.vs.TGFB1 <- select(results.PDGF.TGFB1.vs.TGFB1, "gene", "log2FoldChange", "padj", everything())
write_csv(results.PDGF.TGFB1.vs.TGFB1, "2020.01.results.PDGF.TGFB1.vs.TGFB1.csv")
sig_results.PDGF.TGFB1.vs.TGFB1 <- subset(results.PDGF.TGFB1.vs.TGFB1, padj <0.05)
sig_results.PDGF.TGFB1.vs.TGFB1 <- as.data.frame(sig_results.PDGF.TGFB1.vs.TGFB1)
sig_results.PDGF.TGFB1.vs.TGFB1$gene <- rownames(sig_results.PDGF.TGFB1.vs.TGFB1)
sig_results.PDGF.TGFB1.vs.TGFB1 <- select(sig_results.PDGF.TGFB1.vs.TGFB1, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.PDGF.TGFB1.vs.TGFB1, "2020.01.significant-altered-genes_results.PDGF.TGFB1.vs.TGFB1.csv")

results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1 <- results(ddsASE081, contrast = c("condition", "PDGF+TGFB1+CPI613","PDGF+TGFB1"))
sum(results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1$padj < 0.05, na.rm=TRUE)
results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1
results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1 <- as.data.frame(results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1)
results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1$gene <- rownames(results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1)
results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1 <- select(results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1, "gene", "log2FoldChange", "padj", everything())
write_csv(results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1, "2020.01.results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1.csv")
sig_results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1 <- subset(results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1, padj <0.05)
sig_results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1 <- as.data.frame(sig_results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1)
sig_results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1$gene <- rownames(sig_results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1)
sig_results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1 <- select(sig_results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1, "2020.01.significant-altered-genes_results.PDGF.TGFB1.CPI613.vs.PDGF.TGFB1.csv")

results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1 <- results(ddsASE081, contrast = c("condition", "PDGF+TGFB1+Galloflavin","PDGF+TGFB1"))
sum(results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1$padj < 0.05, na.rm=TRUE)
results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1
results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1 <- as.data.frame(results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1)
results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1$gene <- rownames(results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1)
results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1 <- select(results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1, "gene", "log2FoldChange", "padj", everything())
write_csv(results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1, "2020.01.results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1.csv")
sig_results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1 <- subset(results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1, padj <0.05)
sig_results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1 <- as.data.frame(sig_results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1)
sig_results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1$gene <- rownames(sig_results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1)
sig_results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1 <- select(sig_results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1, "gene", "log2FoldChange", "padj", everything())
write_csv(sig_results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1, "2020.01.significant-altered-genes_results.PDGF.TGFB1.Calloflavin.vs.PDGF.TGFB1.csv")

write.csv(counts(ddsGC, normalized=TRUE), file="2021.11.GC.experiment.normalized-counts.csv")

# Raw counts for each genotype
metabolic.experiment.raw.counts.clean.names <- as.data.frame(countASE081data)
metabolic.experiment.raw.counts.clean.names$gene <- rownames(countASE081data)

metabolic.experiment.raw.counts.clean.names.vehicle <- select(metabolic.experiment.raw.counts.clean.names, "gene", starts_with("A"))
write.csv(metabolic.experiment.raw.counts.clean.names.vehicle, file = "2020.01.metabolic.experiment.raw.counts.clean.names.vehicle.csv")

metabolic.experiment.raw.counts.clean.names.Galloflavin <- select(metabolic.experiment.raw.counts.clean.names, "gene", starts_with("B"))
write.csv(metabolic.experiment.raw.counts.clean.names.Galloflavin, file = "2020.01.metabolic.experiment.raw.counts.clean.names.Galloflavin.csv")

metabolic.experiment.raw.counts.clean.names.CPI613 <- select(metabolic.experiment.raw.counts.clean.names, "gene", starts_with("C"))
write.csv(metabolic.experiment.raw.counts.clean.names.CPI613, file = "2020.01.metabolic.experiment.raw.counts.clean.names.CPI613.csv")

metabolic.experiment.raw.counts.clean.names.PDGFBB <- select(metabolic.experiment.raw.counts.clean.names, "gene", starts_with("D"))
write.csv(metabolic.experiment.raw.counts.clean.names.PDGFBB, file = "2020.01.metabolic.experiment.raw.counts.clean.names.PDGFBB.csv")

metabolic.experiment.raw.counts.clean.names.TGFB1 <- select(metabolic.experiment.raw.counts.clean.names, "gene", starts_with("E"))
write.csv(metabolic.experiment.raw.counts.clean.names.TGFB1, file = "2020.01.metabolic.experiment.raw.counts.clean.names.TGFB1.csv")

metabolic.experiment.raw.counts.clean.names.PDGF.TGFB1 <- select(metabolic.experiment.raw.counts.clean.names, "gene", starts_with("F"))
write.csv(metabolic.experiment.raw.counts.clean.names.PDGF.TGFB1, file = "2020.01.metabolic.experiment.raw.counts.clean.names.PDGF.TGFB1.csv")

metabolic.experiment.raw.counts.clean.names.PDGF.TGFB1.Galloflavin <- select(metabolic.experiment.raw.counts.clean.names, "gene", starts_with("G"))
write.csv(metabolic.experiment.raw.counts.clean.names.PDGF.TGFB1.Galloflavin, file = "2020.01.metabolic.experiment.raw.counts.clean.names.PDGF.TGFB1.Galloflavin.csv")

metabolic.experiment.raw.counts.clean.names.PDGF.TGFB1.CPI613 <- select(metabolic.experiment.raw.counts.clean.names, "gene", starts_with("H"))
write.csv(metabolic.experiment.raw.counts.clean.names.PDGF.TGFB1.CPI613, file = "2020.01.metabolic.experiment.raw.counts.clean.names.PDGF.TGFB1.CPI613.csv")

metabolic.experiment.raw.counts.clean.names.Cholesterol <- select(metabolic.experiment.raw.counts.clean.names, "gene", starts_with("I"))
write.csv(metabolic.experiment.raw.counts.clean.names.Cholesterol, file = "2020.01.metabolic.experiment.raw.counts.clean.names.Cholesterol.csv")

metabolic.experiment.raw.counts.clean.names.Cholesterol.PDGF.TGFB1 <- select(metabolic.experiment.raw.counts.clean.names, "gene", starts_with("J"))
write.csv(metabolic.experiment.raw.counts.clean.names.Cholesterol.PDGF.TGFB1, file = "2020.01.metabolic.experiment.raw.counts.clean.names.Cholesterol.PDGF.TGFB1.csv")
