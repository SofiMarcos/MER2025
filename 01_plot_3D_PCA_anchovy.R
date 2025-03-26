########################################################################
######                                                            ######
######                      Plot a 3D PCA                         ######
######                                                            ######
########################################################################


library(ggplot2)
library(plot3D)

# set working directory
setwd("C:/Users/smarcos007/Documents/04_github/MER2025/anchovy_PCA")

# read bamlist used with ANGSD
bams <- read.table("bamlist")[,1]

# extract sample names from bamlist
samples1 <- sub("/.*/", "", bams)
samples <- sub("_removed_duplicates.bam", "", samples1)

# read covariance matrix generated with PCAngsd
salsal_cov <- as.matrix(read.table("snps.ld_pruned.hwe_filter.cov"))

# append sample names as row and column names to covariance matrix
dimnames(salsal_cov) <- list(samples, samples)

# perform PCA
pca <- prcomp(salsal_cov, scale = TRUE)

# scree plot
eigenval <- pca$sdev^2
explained_var <- 100*(eigenval/sum(eigenval))
df1 <- data.frame(prin_comp = c(seq(1, length(eigenval))),
                  explained_var)

ggplot(df1, aes(prin_comp, explained_var)) +
  geom_col(fill = c(rep("steelblue", 3),
                    rep("grey40", length(eigenval)-3))) +
  xlab("Principal components") +
  ylab("Explained variance") +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
# save scree plot in PDF
ggsave("scree_plot.pdf", width = 170, height = 92, units = "mm")

# create data frame
df2 <- as.data.frame(pca$x)

# Add column for populations
pop <- read.table("anchovy_localities.csv", sep = ";", header = TRUE)

# Check if sample order in the metadata file is the same as in the pca dataframe
pop$Sample == rownames(df2)

# Add two columns to the pca dataframe
df2$subpop <- pop$Subpopulation
df2$pop <- pop$Population

# change the order of groups
df2$pop <- factor(df2$pop,
levels = c("Russia", "North America", "Iceland", "Baltic Sea","Europe"))

df2$subpop <- factor(df2$subpop,
                    levels = c("Russia", 
                               "North America",
                               "Iceland",
                               "Baltic Sea", 
                               "Scotland", "Ireland", "Norway", "Sweden"))

# set color blind friendly palette. Information about colours in ggplot: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
cbPalette <- c("Russia" = "#F72585", #56B4E9
               "North America" = "#B5179E", #D55E00
               "Iceland" = "#7209B7",  #CC79A7
               "Baltic Sea" = "#03045E",  #009E73
               "Scotland" = "#0077B6",  #006046
               "Ireland" = "#4361EE", #0096C7 #F0E442
               "Norway" = "#4895EF", #48CAE4 #0072B2
               "Sweden" = "#4CC9F0")  #ADE8F4 #E69F00

# set point shapes
shapes <- c(18, 15, 17, 8, 16) # Here you can find information about ggplot point shapes: http://www.sthda.com/english/wiki/ggplot2-point-shapes


# simple PCA with 2 components
pdf("anchovy_simple_pca_plot.pdf", width = 8, height = 6)

ggplot(df2, aes(PC1, PC2, color = subpop, shape = pop)) +
  geom_point(size = 5) +
  scale_color_manual(values = cbPalette) +
  scale_shape_manual(values = c(18, 15, 17, 8, 16))

dev.off()


# 3D PCA plot
pdf("anchovy_3dpca_plot.pdf", width = 8, height = 6)

layout(matrix(c(1, 1, 1, 0,
                1, 1, 1, 0,
                1, 1, 1, 2,
                1, 1, 1, 2),
              nrow = 4, ncol = 4, byrow = TRUE))
par(mar = c(0, 0, 0, 0))
scatter3D(df2$PC1, df2$PC2, df2$PC3, bty = "g", theta = 30, phi = 45,
          colvar = NULL, colkey = FALSE, col = cbPalette[as.factor(df2$subpop)],
          pch = shapes[as.factor(df2$pop)], cex = 2.5, cex.lab = 2,
          xlab = paste("PC1 (", round(explained_var[1], 2), "%)", sep = ""),
          ylab = paste("PC2 (", round(explained_var[2], 2), "%)", sep = ""),
          zlab = paste("PC3 (", round(explained_var[3], 2), "%)", sep = ""))
plot.new()
legend("left",
       c("Russia", 
         "North America",
         "Iceland",
         "Baltic Sea", 
         "Scotland", "Ireland", "Norway", "Sweden"),
       col = cbPalette,
       # pch = c(18, 15, 17, 8, 16, 16, 16, 16), # Should be the same order than in line 88
       pch = c(18, 15, 17, 8, rep(16, 4)), # Instead writing '16' four times for the four European subpopulations we write 'rep(16,4)'
       pt.cex = 2.5, cex = 1.5, y.intersp = 1.5, bty = "n")

invisible(dev.off())

