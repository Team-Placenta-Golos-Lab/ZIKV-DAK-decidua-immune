
library(Spectre)
Spectre::package.check()    # Check that all required packages are installed
Spectre::package.load()  

setwd("~/Desktop/Cytokines_decidua")

data <- read.csv("cytokines_extrapolated_decidua_only.csv")
sum.dat <- read.csv("cytokines_extrapolated_decidua_only.csv")
sum.dat
as.matrix(names(sum.dat))
sample.col <- "Sample"
group.col <- "Group"
treatment.col <- "Treatment"
#gdgroup.col <- "gd.group"
gd.col <- "gd"

as.matrix(names(sum.dat))

annot.cols <- c(group.col, treatment.col)

plotperct.cols <- names(sum.dat)[c(6:17)]
plotperct.cols
#sum.dat <- do.reorder(sum.dat, group.col, grp.order)
#sum.dat[,c(1:3)]
group.col
grp.order <- c("7-dpi-Control", "7-dpi-ZIKV", "14-dpi-Control","14-dpi-ZIKV")
grp.order
sum.dat <- do.reorder(sum.dat, group.col, grp.order)
sum.dat[,c(1:3)]

make.pheatmap(sum.dat, 
              sample.col = sample.col, 
              plot.cols = plotperct.cols, 
              plot.title = 'Decidual Chemokines',
              annot.cols = annot.cols,
              dendrograms = "both",
              normalise = TRUE,
              standard.colours = "Blues",
              file.name = "Blues_no_gd.png")
make.pheatmap(sum.dat, 
              sample.col = sample.col, 
              plot.cols = plotperct.cols, 
              plot.title = 'Decidual Chemokines',
              annot.cols = annot.cols,
              dendrograms = "both",
              dendrograms.sort = TRUE,
              normalise = TRUE,
              standard.colours = "Blues",
              file.name = "Blues_more_annots_sort.png")
#file.name DEFAULT = paste0("Pheatmap by ", sample.col, ".png").
#@param cell.size DEFAULT = 15. Numeric.
#@param standard.colours DEFAULT = "BuPu". Character. Can also be "RdYlBu", "YlGnBu", "viridis", "magma", "inferno", "spectral", "Blues", "Reds", "Greys", or "rev(RdBu)"