library('topGO')

# Specify the conditions to compare (EM/EA/LM/LA)
compare1 = "EM"
compare2 = "EA"
# Specify analysis of up-regulated (Up) or down-regulated (Down)
direction <- "Down"

geneID2GO <- readMappings(file = choose.files(caption = "Select the GO terms table"))
geneNames <- names(geneID2GO)

comparison <- paste(compare1, 'v', compare2, '_', direction, sep="")

myInterestingGenes <- scan(paste0(choose.dir(caption = "Select query DE gene folder"),"\\",
                                  compare1, "v", compare2,"\\", direction, "regulated\\",
                                  "genes_of_interest.txt"), "character")
myInterestingGenes
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))

names(geneList) <- geneNames

GOdata <- new("topGOdata", description = comparison, ontology = "BP",
              allGenes = geneList, nodeSize = 5, annot = annFUN.gene2GO,  gene2GO = geneID2GO)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultWeightFisher <- runTest(GOdata, algorithm='weight01', statistic='fisher')

allGO <- usedGO(GOdata)
allRes <- GenTable(GOdata, weightFisher = resultWeightFisher,
                   orderBy = "weightFisher", topNodes=length(allGO))

p.adj <- p.adjust(allRes$weightFisher,method="holm")
allRes_final <- cbind(allRes,p.adj)
allRes_final <- allRes_final[order(allRes_final$p.adj),]
allRes_final$weightFisher <- as.numeric(allRes_final$weightFisher)

results.table.p = allRes_final[which(allRes_final$weightFisher<=0.05),]
results.table.bh =allRes_final[which(allRes_final$p.adj<=0.10),]

write.table(results.table.p, paste(comparison, "_GO_enrichment.tsv", sep = ""), sep="\t", quote=FALSE, row.names = FALSE) # use for graphing top 10 GO terms
write.table(results.table.bh, paste(comparison, "_GO_enrichment_holm.tsv", sep = ""), sep="\t", quote=FALSE, row.names = FALSE)
