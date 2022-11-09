library("tximport")
library("readr")
library("jsonlite")
library("DESeq2")

compare1 = "LM"
compare2 = "LA"
output_dir = "E:/School Stuff/Lab Documents/RNA_Seq_Analysis/Refined Analysis/Refined_Analysis/Clone_Stage_Model/EMvEA"

dir <- "E:/School Stuff/Lab Documents/RNA_Seq_Analysis/Refined Analysis/Salmon_Quant_Full_Real_Name"
sample_names <- list.files(dir, pattern = paste("-", compare1, "-|-", compare2, "-", sep = ""))


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
                                   design = ~ Clone + Stage)

dds <- DESeq(ddsTxi)
res <- results(dds, contrast=c("Stage", compare1, compare2))

write.table(res,
            file.path(output_dir,paste(compare1,"v",compare2,"_","clone_stage_all.txt", sep=""))
            ,sep="\t", quote=FALSE)

write.table(res[ which(res$padj < 0.05 & abs(res$log2FoldChange) > 0.5849625007211562),],
            file.path(output_dir,paste(compare1,"v",compare2,"_","clone_stage","_0.05padj_1.5FC.txt", sep=""))
            ,sep="\t", quote=FALSE)

write.table(res,
            file.path(output_dir,paste(compare1,"v",compare2,"_","clone_stage","_all.txt", sep=""))
            ,sep="\t", quote=FALSE)

tbl <- res[ which(res$padj < 0.05 & abs(res$log2FoldChange) > 0.5849625007211562),]
sigDE <- rownames(tbl) 
write.table(sigDE, file.path(output_dir, "AroMoose&LK16&PA42&POVI4&SW4.txt"), sep="\t", quote=FALSE, 
            row.names = FALSE, col.names = FALSE)


library("pcaExplorer")
pcaExplorer(dds = dds)
