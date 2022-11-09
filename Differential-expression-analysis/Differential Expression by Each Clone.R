library("tximport")
library("readr")
library("jsonlite")
library("DESeq2")

clone_list <- c("AroMoose", "LK16", "PA42", "POVI4", "SW4")
stage1 <- "LA"
stage2 <- "EA"

for (clone in clone_list) {

dir <- "E:/School Stuff/Lab Documents/RNA_Seq_Analysis/Refined Analysis/Salmon_Quant_Full_Real_Name"
sample_names <- list.files(dir)

samples_daphnia <- data.frame("Sample_Dir" = c(sample_names))
samples_daphnia$Clone <- vapply(strsplit(as.character(samples_daphnia$Sample_Dir), '-'), '[', 1, FUN.VALUE=character(1))
samples_daphnia$Stage <- vapply(strsplit(as.character(samples_daphnia$Sample_Dir), '-'), '[', 2, FUN.VALUE=character(1))
samples_daphnia$Replica <- vapply(strsplit(as.character(samples_daphnia$Sample_Dir), '-'), '[', 3, FUN.VALUE=character(1))
rownames(samples_daphnia) <- samples_daphnia$Sample_Dir

files <- file.path(dir, samples_daphnia$Sample_Dir, "quant.sf")

names(files) <- samples_daphnia$Sample_Dir

tx2gene <- read_csv("E:/School Stuff/Lab Documents/RNA_Seq_Analysis/Refined Analysis/Refined_Analysis/tx_file.cvs")
txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples_daphnia,
                                   design = ~ Stage)
ddsTxi$group <- factor(paste0(ddsTxi$Clone, ddsTxi$Stage))
design(ddsTxi) <- ~ group

dds <- DESeq(ddsTxi)
res <- results(dds, contrast=c("group", paste0(clone,stage1), paste0(clone,stage2)))

write.table(res[ which(res$padj < 0.05 & res$log2FoldChange > 0.5849625007211562) ,],paste0(stage1,"v",stage2,"_", clone,"_up_0.05padj_1.5FC.txt", sep=""), quote=FALSE)
write.table(res[ which(res$padj < 0.05 & res$log2FoldChange < 0.5849625007211562) ,],paste0(stage1,"v",stage2,"_", clone,"_down_0.05padj_1.5FC.txt", sep=""), quote=FALSE)

}
