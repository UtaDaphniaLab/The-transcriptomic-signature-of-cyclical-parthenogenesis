library(dplyr)
library(tidyr)

# change to match number of genes in reference used (18440 for PA42-3.0)
gene.ref.num = 18440
filter.pway.min.genes = 5

#load gene pathway table
gene.pway.table = read.table(choose.files(caption = "Select Kegg pathway table"), sep="\t", header=FALSE) #choose gene-pathway-table
gene.pway.table = gene.pway.table[gene.pway.table$V3 > 5,]
#fix(gene.pway.table)
num.genes=nrow(gene.pway.table %>% count(V1)) #num.genes is 6282. 6282 genes mapped to pathway
num.genes
mapped.genes = length(gene.pway.table$V1) #30319 genes (including genes repeated many times) in pathways
mapped.genes

pway.table = gene.pway.table[,-1]
#fix(pway.table)
pway.table = distinct(pway.table, V2, .keep_all = TRUE)
npathway = length(pway.table$V2) #385 pathways identified in Daphnia pulex genome
npathway


#load gene list
df1=read.table(choose.files(caption = "Select query gene list"), sep="\t", header=F) #choose gene-list.txt or your own gene list
#replaces all "mRNA" prefix with "gene" prefix in query gene table to match gene.pway.table
df1 <- data.frame(lapply(df1, function(x) {gsub("mRNA","gene",x)}))

df6 = merge (df1, gene.pway.table, by.x="V1", by.y="V1", all.x=FALSE, all.y=FALSE) #find pathway for each gene. remove genes not mapped to pathway
#fix(df6)

write.table(df6, file="df6.txt", sep="\t", col.names=TRUE, row.names = FALSE, quote=FALSE)
df7 = df6 %>% count(V2)
#fix(df7) 
dim(df7) # genes are in 60 pathways
test.table = merge(df7, pway.table, by.x= "V2", by.y="V2")



#fix(test.table)
colnames(test.table) = c("pathway", "wht.drawn", "wht.in.urn")
# no. of genes: 19
# black in urn: number of genes not mapped to pathway of interest
test.table$blk.in.urn = num.genes - test.table$wht.in.urn

# yields the number of unique DEGs/genes of interest that are mapped to a pathway
test.table$total.draw = nrow(df6 %>% count(V1))
test.table$total.draw
#fix(test.table)
#perform hypergeometric test
result = apply(test.table, 1, function(x) phyper(as.numeric(x[2])-1, as.numeric(x[3]), as.numeric(x[4]), as.numeric(x[5]),lower.tail= FALSE))
#fix(result)
test.table = data.frame(test.table, result)
#fix(test.table)


test.table$p.adjust = round(p.adjust(test.table$result, "BH"),4)

#only sorts by asc p value
test.table = test.table[order(test.table$result), ]
# line below will filter out all genes not significant by the adjusted p value
# test.table = test.table[which(test.table$p.adjust < 0.05),]
write.table(test.table, file=file.choose(), sep = "\t", col.names=TRUE, row.names = FALSE, quote=FALSE)
