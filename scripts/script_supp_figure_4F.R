# Load packages
library(genefilter)
library(pvclust)

# Read R object
path <- "/Users/seanwen/Documents/Cristina/Github/DiGenua_Hematologica_2019/data/"
file <- "object.rdata"
load(paste(path, file, sep=""))

# Retrieve transformed RPKM table
df <- object$rpkm.log2.matrix

# Set factor levels
df$Group <- factor(df$Group, levels=c("CON", "AM", "KM", "AKM"))

# Identify significant gene sets
  # Perform F-tests
  ftest <- colFtests(as.matrix(df[,-1]), fac=df$Group, var.equal=FALSE)

  # Correct for multiple testing
  ftest$fdr <- p.adjust(ftest$p.value, method="fdr")

  # Subset for significant genes
  ftest <- ftest[which(ftest$fdr < 0.05),]
  ftest$Gene <- row.names(ftest)

# Retrieve significant genes
df <- df[, names(df) %in% c("Group", ftest$Gene)]

# Transpose data frame
n <- df$Group
df <- as.data.frame(t(df[,-1]))
colnames(df) <- n

# Give unique names to genotypes
df <- df[, order(names(df))]

# Reorder rows for aesthetic purpose
conam <- df[, names(df) %in% c("CON", "CON.1", "CON.2", "CON.3", "CON.4", "CON.5", "AM", "AM.1", "AM.2", "AM.3", "AM.4", "AM.5")]
df$conam_mean <- rowMeans(conam, na.rm=TRUE)
df <- df[order(df$conam_mean), ]
df$conam_mean <- NULL

# Run clustering
result <- pvclust(df[ ,c(6:11, 12:17, 1:5, 18:23)], method.dist="euclidean", method.hclust="complete", nboot=1000)

# Plot dendrogram
path <- "/Users/seanwen/Documents/Cristina/Github/DiGenua_Hematologica_2019/figures/"
file <- "supp_figure_4F.pdf"
pdf(paste(path, file, sep=""))

plot(result)
#pvrect(result, border="none")

dev.off()
