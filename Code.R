# read gene expression data from WD

csvfile <- "COREAD_CMS13_gex.csv"

x1 <- read.csv(csvfile, row.names=1)


# Fix the gene names in the data frame


rownames(x1) <- sapply(strsplit(rownames(x1), "\\|"), function(x) x[1])

# Output a table

library(knitr)

knitr::kable(head(t(head(x1))), caption="Example gene expression data (head)")

write.csv(x1, "Example gene expression data.csv", row.names = TRUE)


#read mutation data from WD


csvfile <- "COREAD_CMS13_muts.csv"

x2 <- read.csv(csvfile, row.names=1)

x2[x2>0]=1

knitr::kable(head(t(head(x2))), caption="Example mutation data (head)")

write.csv(x1, "Example mutation data.csv", row.names = TRUE)

#read copy number data from WD

csvfile <- "COREAD_CMS13_cnv.csv"

x3 <- read.csv(csvfile, row.names=1)

knitr::kable(head(t(head(x3))), caption="Example copy number data for CRC samples")

write.csv(x1, "Example copy number data for CRC samples.csv", row.names = TRUE)

#read subtypes data from WD

csvfile <- "COREAD_CMS13_subtypes.csv"

covariates <- read.csv(csvfile, row.names=1)

# Fix the TCGA identifiers so they match up with the omics data

rownames(covariates) <- gsub(pattern = '-', replacement = '\\.',
                             rownames(covariates))


covariates <- covariates[colnames(x1),]

# create a dataframe which will be used to annotate later graphs

anno_col <- data.frame(cms=as.factor(covariates$cms_label))
rownames(anno_col) <- rownames(covariates)


#plot a heatmap and cluster the tumors, while displaying a color annotation on top of the heatmap, indicating which subtype each tumor belongs to.


library(pheatmap)

pheatmap::pheatmap(x1,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Gene expression data")

# heatmap for mutation data

pheatmap::pheatmap(x2,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Mutation data")




# heatmap for CNV data



pheatmap::pheatmap(x3,
                   annotation_col = anno_col,
                   show_colnames = FALSE,
                   show_rownames = FALSE,
                   main="Copy number data")





#https://compgenomr.github.io/book/use-case-multi-omics-data-from-colorectal-cancer.html


























