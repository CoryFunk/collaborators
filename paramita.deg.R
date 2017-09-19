# DEGs for Paramita
library(edgeR)
library(dplyr)

#read in the data
paramita <- read.csv("~/Downloads/paramita/paramita.rnaseq.csv",header=TRUE,sep=",",stringsAsFactors=FALSE,check.names=FALSE)

# subsone_four the data to gone_four just the counts
rownames(paramita) <- paramita[,1]
paramita2 <- paramita[,3:ncol(paramita)]

# assign groups based on Paramita's xlsx sheone_four
group <- c("1","1","1","2","2","3","3","4","4","4","5","5","5")
group <- as.matrix(group, ncol = 1)
rownames(group) <- colnames(paramita2)
colnames(group) <- "group"

# convert to cpm and drop the lowely expressed genes
cpm_x <- cpm(paramita2)
keep <- (rowSums(cpm_x>=1) >= 6) 
paramita2 <- paramita2[keep,]

# calculate the dispersions
y0 <- DGEList(counts=paramita2, group=group)
y0 <- estimateCommonDisp(y0)
y0 <- estimateTagwiseDisp(y0)

# do a pairwise comparison
one_four <- exactTest(y0, pair = c(1, 4))
four_five <- exactTest(y0, pair = c(4, 5))
one_three <- exactTest(y0, pair = c(1, 3))

# pull out the pvalue and do an FDR
one_four$table$FDR <- p.adjust(one_four$table$PValue, "fdr")
four_five$table$FDR <- p.adjust(four_five$table$PValue, "fdr")
one_three$table$FDR <- p.adjust(one_three$table$PValue, "fdr")

# add back the gene names
# one_four
one_four$table$ensmbl_id <- rownames(one_four$table)
names(paramita)[1] <- "ensmbl_id" 
one_four$table <- left_join(one_four$table, paramita, by = "ensmbl_id") %>% select(1:6)
one_four$table <- one_four$table[,c(5,6,1,2,3,4)]

# four_five
four_five$table$ensmbl_id <- rownames(four_five$table)
names(paramita)[1] <- "ensmbl_id" 
four_five$table <- left_join(four_five$table, paramita, by = "ensmbl_id") %>% select(1:6)
four_five$table <- four_five$table[,c(5,6,1,2,3,4)]

# one_three
one_three$table$ensmbl_id <- rownames(one_three$table)
names(paramita)[1] <- "ensmbl_id" 
one_three$table <- left_join(one_three$table, paramita, by = "ensmbl_id") %>% select(1:6)
one_three$table <- one_three$table[,c(5,6,1,2,3,4)]

# identify the genes that are above a fold change and FDR
one_four$table$threshold = as.factor(abs(one_four$table$logFC) > 2 & one_four$table$FDR < 0.05)
four_five$table$threshold = as.factor(abs(four_five$table$logFC) > 2 & four_five$table$FDR < 0.01)
one_three$table$threshold = as.factor(abs(one_three$table$logFC) > 2 & one_three$table$FDR < 0.05)

DEG_1vs4 <- head(arrange(one_four$table, desc(threshold), FDR), 200)
write.table(DEG_1vs4, file = "~/Downloads/DEG_1vs4.txt", sep = '\t', quote=FALSE, row.names=FALSE)

DEG_4vs5 <- head(arrange(four_five$table, desc(threshold), FDR), 1000)
write.table(DEG_4vs5, file = "~/Downloads/DEG_4vs5.txt", sep = '\t', quote=FALSE, row.names=FALSE)

DEG_1vs3 <- head(arrange(one_three$table, desc(threshold), FDR), 200)
write.table(DEG_1vs3, file = "~/Downloads/DEG_1vs3.txt", sep = '\t', quote=FALSE, row.names=FALSE)


#----------------------------------------------------------------------------------------------------
#plots
# one_four
g = ggplot(data=one_four$table, aes(x=logFC, y=-log10(PValue), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-10, 10)) + ylim(c(0, 30)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  scale_colour_discrete(name = "Significant") +
  geom_text(aes(label=ifelse(abs(logFC)>2.0 & threshold=="TRUE", as.character(gene_id),'')),hjust=1.2,vjust=-0.2)
g

# four_five
g = ggplot(data=four_five$table, aes(x=logFC, y=-log10(PValue), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-10, 10)) + ylim(c(0, 100)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  scale_colour_discrete(name = "Significant") +
  geom_text(aes(label=ifelse(abs(logFC)>4.0 & threshold=="TRUE", as.character(gene_id),'')),hjust=1.2,vjust=-0.2)
g

# one_three
g = ggplot(data=one_three$table, aes(x=logFC, y=-log10(PValue), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-10, 10)) + ylim(c(0, 70)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  scale_colour_discrete(name = "Significant") +
  geom_text(aes(label=ifelse(abs(logFC)>4.0 & threshold=="TRUE", as.character(gene_id),'')),hjust=1.2,vjust=-0.2)
