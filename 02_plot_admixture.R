###########################################################################
######                                                               ######
######      Plot ADMIXTURE                                           ######
######                                                               ######
###########################################################################
# Multiple plot function
#
# Source: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# load libraries
library(ggplot2)
library(reshape2)

# set working directory
setwd("C:/Users/smarcos007/Documents/04_github/MER2025/admixture_data")

# read bamlist used with ANGSD
bams <- read.table("bamlist")[,1]

# extract sample names from bamlist
samples_1 <- sub("/.*/", "", bams)
samples <- sub("_removed_duplicates.bam", "", samples_1)

# create a rank to keep the sample order of the bamlist
rank <- c(seq(1, length(samples)))

# get input file names
file.names <- dir("./", pattern = ".qopt")

# generate an admixture barplot for each K
for(i in 1:length(file.names)){
  # read input Q-matrix
  q.matrix <- read.table(paste0("./", file.names[i]))
  # create a data frame with rank, sample names, and Q-matrix
  data <- cbind(rank, samples, q.matrix)
  # reformat the data frame for plotting
  data.m <- melt(data, id = c("rank", "samples"),
                 value.name = "proportion",
                 variable.name = "ancestry")
  
  # assign plot to a variable
  assign(paste0("p", i),
         # read data frame and map aesthetics for plotting
         ggplot(data.m, aes(x = rank, y = proportion, fill = ancestry)) +
           # add barplot
           geom_bar(stat = "identity") +
           # add solid vertical lines between species
           geom_vline(xintercept = c(14.5, 18.5, 23.5),
                      color = "black",
                      lwd = 1.25) +
           # add dashed vertical lines between species
           geom_vline(xintercept = c(27.5, 35.5, 45.5, 59.5),
                      color = "black",
                      lty = 2) +
           # add y label
           ylab(paste("K =", i)) +
           # add y label
           xlab(paste("samples")) +
           # change axis elements, remove legend and background panel
           theme(axis.title.x = element_blank(),
                 axis.title.y = element_text(size = 16),
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 legend.position = "none",
                 panel.background = element_blank()))
}

# render multiple plots in a single page (colors set by trial and error)
multiplot(p2 + scale_fill_manual(values = c("#F72585", "#4361EE")),
          p3 + scale_fill_manual(values = c("#F72585", "#4361EE", "#B5179E")),
          p4 + scale_fill_manual(values = c("#B5179E", "#4361EE", "#F72585", "#03045E")),
          p5 + scale_fill_manual(values = c("#4361EE", "#03045E", "#F72585", "#7209B7", "#B5179E")),
          p6 + scale_fill_manual(values = c("#F72585", "#B5179E", "#7209B7", "#009E73", "#4361EE", "#03045E")),
          p7 + scale_fill_manual(values = c("#7209B7", "#B5179E", "#0077B6", "#F72585", "#03045E", "#009E73", "#4361EE")),
          p8 + scale_fill_manual(values = c("#03045E", "#009E73", "#F72585", "#4CC9F0", "#4361EE", "#0077B6", "#7209B7", "#B5179E")),
          p9 + scale_fill_manual(values = c("#F72585", "#4361EE", "#009E73", "#7209B7", "#0077B6", "#4CC9F0", "#B5179E", "#4895EF", "#03045E")),
          p1 + scale_fill_manual(values = c("#7209B7", "#4895EF", "#009E73", "#03045E", "#F72585", "#CC79A7", "#4361EE", "#4CC9F0", "#0077B6", "#B5179E")) +
          # fix y label
          ylab("K = 10")
          )

#F0E442 - yellow
#009E73 - green
#56B4E9 - blue
#CC79A7 - pink
#E69F00 - orange
#006046 - dark green      
#D55E00 - Dark orange         
#0072B2 - Dark blue
#666666 - Grey
#873b00 - Brown

