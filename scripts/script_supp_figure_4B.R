# Load packages
library(ggplot2)

# Read R object
path <- "/Users/seanwen/Documents/Cristina/Github/DiGenua_Hematologica_2019/data/"
file <- "object.rdata"
load(paste(path, file, sep=""))

# Retrieve transformed RPKM table
df <- object$rpkm.log2.matrix

# Transpose data frame
n <- df$Group
df <- as.data.frame(t(df[,-1]))
colnames(df) <- n

# Identify highly variable genes
group <- c("CON", "AM", "KM", "AKM")

genes.var_list <- list()

for(i in 1:length(group)) {

    # Subset group
    sub <- df[ , which(names(df)==group[i])]

    # Compute mean
    means <- rowMeans(sub, na.rm=TRUE)

    # Compute variance
    vars <- NULL
    for(m in 1:nrow(sub)) {
        vars[m] <- var(as.numeric(sub[m, ]), na.rm=TRUE)
    }

    # Compute coefficient of variance
    CV <- vars/means^2

    # Compile mean, var, CV and perform filtering
    summary <- data.frame(means, vars, CV)
    summary <- summary[which(!is.nan(summary$CV)), ]
    summary <- summary[which(summary$means > 1), ]

    # Retrieve top 10% most variable genes
    summary <- summary[order(summary$CV, decreasing=TRUE), ]
    genes.var_list[[i]] <- head(row.names(summary), 0.1*nrow(summary))

}

genes.var <- unlist(genes.var_list)

# Remove duplicate genes
genes.var <- unique(genes.var)
length(genes.var)

# Subset original data frame
df <- df
df$Gene <- row.names(df)
df <- df[which(df$Gene %in% genes.var), ]
df <- df[, c(ncol(df), 1:(ncol(df)-1)) ]

# Transpose data frame
n <- df$Gene
df <- as.data.frame(t(df[,-1]))
colnames(df) <- n
df$Group <- row.names(df)
df <- df[, c(ncol(df), 1:(ncol(df)-1)) ]
df$Group <- gsub("\\.[0-9]", "", df$Group)

# Set factor levels
df$Group <- factor(df$Group, levels=c("CON", "AM", "KM", "AKM"))

# PCA
    # Reduce dimensions
    pca <- prcomp(df[,-1], center=TRUE, scale.=TRUE)

    # Retrieve proportion of variance dflained
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
        file <- "supp_figure_4B.pdf"
        ggsave(paste(path, file, sep=""), plot, width=3.5, height=3)
