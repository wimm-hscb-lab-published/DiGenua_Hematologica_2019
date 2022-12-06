# Load packages
library(org.Mm.eg.db) # 3.12.0
library(GO.db) # 3.12.1
library(GOstats) # 2.56.0

# Read R object
path <- "/Users/seanwen/Documents/Cristina/Github/DiGenua_Hematologica_2019/data/"
file <- "object.rdata"
load(paste(path, file, sep=""))

# Retrieve Entrez Gene ID for all genes
    # Retrieve DE results table
    diff <- object$de.table

    # Retrieve
    ID_all <- select(org.Mm.eg.db, keys=diff$Geneid, columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL")

    # Remove non-matches
    ID_all <- ID_all[which(!is.na(ID_all$ENTREZID)), "ENTREZID"]

# Subset significant genes
up <- diff[which(diff$FDR < 0.05 & diff$logFC > 0), "Geneid"]

# Retrieve Entrez Gene ID for significant genes
    # Retrieve
    ID <- select(org.Mm.eg.db, keys=up, columns=c("ENTREZID", "SYMBOL"), keytype="SYMBOL")

    # Remove non-matches
    ID <- ID[which(!is.na(ID$ENTREZID)), "ENTREZID"]

# Test of over-representation
    # Set up parameters for analysis
    params <- new('GOHyperGParams', geneIds=ID, universeGeneIds=ID_all, annotation='org.Mm.eg.db', ontology="BP", pvalueCutoff=0.05, conditional=TRUE, testDirection="over")

    # Analyse
    go <- hyperGTest(params)

    # Generate result table
    go.table <- summary(go)

    # Filter significant terms after adjustment
    go.table$bonferroni <- p.adjust(go.table$Pvalue, method="bonferroni")
    
    # Indicate no. of genes up-regulated
    go.table$Total_genes_upregulated <- length(up)
    
    # Calculate precentage of hits
    go.table$Count_pct <- round((go.table$Count/go.table$Total_genes_upregulated)*100, 2)
    
    # Include comparison label
    go.table$Comparison <- "AKMvsAM"
    
    # Reorder columns
    go.table <- go.table[, c(1, 3:5, 9:10, 6:7, 2, 8, 11)]

# Keep adjusted p-values < 0.05
#go.table <- go.table[which(go.table$bonferroni < 0.10), ]

# Write file
path <- "/Users/seanwen/Documents/Cristina/Github/DiGenua_Hematologica_2019/tables/"
file <- "supp_table_4.pdf"
write.table(go.table, paste(path, file, sep=""), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
