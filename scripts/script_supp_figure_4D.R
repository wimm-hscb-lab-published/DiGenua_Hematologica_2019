# Load packages
library(genefilter)
library(ggplot2)

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

# PCA
    # Reduce dimensions
    pca <- prcomp(df[,-1], center=TRUE, scale.=TRUE)

    # Retrieve proportion of variance explained
    pca_var <- summary(pca)
    
    # Scatterplot
        # Definitions
        data <- data.frame(pca$x[,1], pca$x[,1])
        x <- pca$x[,1]
        y <- pca$x[,2]
        group <- df$Group
        xtitle <- paste("PC1 ", "(", signif(pca_var$importance[2,1]*100, 2), "%)", sep="")
        ytitle <- paste("PC2 ", "(", signif(pca_var$importance[2,2]*100, 2), "%)", sep="")
        legend.title <- "Group"
        color_group <- c("CON", "AM", "KM", "AKM")
        color <- c("blue", "red", "green", "purple")

        # Plot
        plot <- ggplot() +
           geom_point(data, mapping=aes(x=x, y=y, color=group), size=2) +
           labs(x=xtitle, y=ytitle) +
           scale_color_manual(breaks=color_group, values=color, name=legend.title) +
           theme_classic() +
           theme(panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank(),
                 plot.title = element_text(size=12, hjust=0.5),
                 axis.line=element_line(colour = "black"),
                 axis.title=element_text(size=12),
                 axis.text.x=element_text(size=10, colour="black"),
                 axis.text.y=element_text(size=10, colour="black"),
                 legend.title=element_text(size=8),
                 legend.text=element_text(size=8)
                 )

        # Save plot
        path <- "/Users/seanwen/Documents/Cristina/Github/DiGenua_Hematologica_2019/figures/"
        file <- "supp_figure_4D.pdf"
        ggsave(paste(path, file, sep=""), plot, width=3.5, height=3)
