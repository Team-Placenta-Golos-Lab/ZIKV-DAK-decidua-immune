##########################################################################################################
#### 1. Load packages, and set working directory
##########################################################################################################

### Load libraries

library(Spectre)
Spectre::package.check()    # Check that all required packages are installed
Spectre::package.load()     # Load required packages

### Set PrimaryDirectory

setwd("~/Desktop/all_Decidua_8Dec2022_used")
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

### Set 'input' directory

setwd(PrimaryDirectory)
setwd("~/Desktop/all_Decidua_8Dec2022_used/data/")
InputDirectory <- getwd()
setwd(PrimaryDirectory)

### Set 'metadata' directory

setwd(PrimaryDirectory)
setwd("~/Desktop/all_Decidua_8Dec2022_used/metadata/")
MetaDirectory <- getwd()
setwd(PrimaryDirectory)

### Create output directory

dir.create("Output_Spectre", showWarnings = FALSE)
setwd("Output_Spectre")
OutputDirectory <- getwd()
setwd(PrimaryDirectory)

##########################################################################################################
#### 2. Import and prep data
##########################################################################################################

### Import data

setwd(InputDirectory)
list.files(InputDirectory, ".csv")
##################Note: the files that I imported were unmixed in Specter and pre-gated in Flowjo. Flowjo gating selected live CD45+ cells and excluded debris, dead cells and obvious clumps of cells, aiming to only get single cells.I celled this population all immune cells 
################### Samples that had more than 50,000 'all immune cells' were downsampled to 50,000 cells. Using the flowjo plugin DownSampleV3. Samples were exported with 'compensated' channels as .CSV files, scaled values. It is important to export as scaled vlaues if you are doing the transformation in r as I do in this script. 

### Check the data

check <- do.list.summary(data.list)

check$name.table # Review column names and their subsequent values
check$ncol.check # Review number of columns (features, markers) in each sample
check$nrow.check # Review number of rows (cells) in each sample

data.list[[1]]

### Merge data  ###go back to here to try differnt co-factor for transformation 

cell.dat <- Spectre::do.merge.files(dat = data.list)
cell.dat

### Read in metadata  

setwd(MetaDirectory)

meta.dat <- fread("sample.details.csv")
meta.dat

##########################################################################################################
#### 3. Data transformation
##########################################################################################################
setwd(OutputDirectory)
dir.create("Output 1 - transformed plots")
setwd("Output 1 - transformed plots")


### Arcsinh transformation
##First, check the column names of the dataset.
as.matrix(names(cell.dat))
##[2,] "SSC-A"   
##[3,] "DC-SIGN" 15,000 
##[4,] "HLA-DR" 8,000
##[5,] "CD56"    5000 1
##[6,] "CD45"    
##[7,] "CCR6"    15000 
##[8,] "CD45RA"  5000 2
##[9,] "CD4"     5000 3
##[10,] "CD69"    5000 4
##[11,] "CD8a"    8000
##[12,] "FoxP3"   5000 5
##[13,] "CD16"    8000
##[14,] "CD14"    5000 6
##[15,] "CD49a"   5000 7 
##[16,] "CD163"   5000 8
##[17,] "CD3"     5000 9
##[18,] "CCR7"    15000 
##[19,] "RORyt"   15000
##[20,] "T-bet"   5000 10
##[21,] "NKp46"   15000
##[22,] "CD11c"   5000 11
##[23,] "CD127"   5000 12
##[24,] "CD86"    5000 13
##[25,] "CD20"    5000 14
##[26,] "Eomes"   5000 15
##[27,] "Time"    
##[28,] "FileName"
##[29,] "FileNo"
## specify the columns to apply arcsinh transformation to

as.matrix(names(cell.dat))
#Define the cofactor 
###set #1
to.asinh <- names(cell.dat)[c(5,8:10,12,14:17,20,22:26)] 
  
to.asinh
#In this worked example we will use a cofactor of 500, the recommended vlaue to start with for flow data
cofactor <- 5000
##You can also choose a column to use for plotting the transformed result – ideally something that is expressed on a variety of cell types in your dataset.

plot.against <- "CD45"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 20000), i, plot.against)
}

plot.against <- "CD3_asinh"
for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 20000), i, plot.against)
}
###set #2
to.asinh <- names(cell.dat)[c(3,7,18,19,21)] 
to.asinh
cofactor <- 15000

plot.against <- "CD45"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 20000), i, plot.against)
}
plot.against <- "CD3_asinh"
for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 20000), i, plot.against)
}

###set #3
to.asinh <- names(cell.dat)[c(4,13,11)] 
to.asinh
cofactor <- 8000
plot.against <- "CD45"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 20000), i, plot.against)
}
plot.against <- "CD3_asinh"
for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 20000), i, plot.against)
}

### Add metadata to data.table

meta.dat
sample.info <- meta.dat[,c(1:4)]
sample.info

meta.dat
counts <- meta.dat[,c(2,5)]
counts

cell.dat <- do.add.cols(cell.dat, "FileName", sample.info, "Filename", rmv.ext = TRUE)
cell.dat

### Columns

as.matrix(names(cell.dat))
##[,1]           
##[1,] "FSC-A"        
#... 
##[30,] "CD56_asinh"   
##[31,] "CD45RA_asinh" 
##[32,] "CD4_asinh"    
##[33,] "CD69_asinh"   
##[34,] "FoxP3_asinh"  
##[35,] "CD14_asinh"   
##[36,] "CD49a_asinh"  
##[37,] "CD163_asinh"  
##[38,] "CD3_asinh"    
##[39,] "T-bet_asinh"  
##[40,] "CD11c_asinh"  
##[41,] "CD127_asinh"  
##[42,] "CD86_asinh"   
##[43,] "CD20_asinh"   
##[44,] "Eomes_asinh"  
##[45,] "DC-SIGN_asinh"
##[46,] "CCR6_asinh"   
##[47,] "CCR7_asinh"   
##[48,] "RORyt_asinh"  
##[49,] "NKp46_asinh"  
##[50,] "HLA-DR_asinh" 
##[51,] "CD16_asinh"   
##[52,] "CD8a_asinh"   
##[53,] "Sample"       
##[54,] "Group"        
##[55,] "Batch"    
cellular.cols <- names(cell.dat)[c(30:52)]
as.matrix(cellular.cols)
     
cluster.cols <- names(cell.dat)[c(30:32,35:41,43:47,49:52)] 
as.matrix(cluster.cols)

exp.name <- "77-MRK_Decidua_8Dec2022"
sample.col <- "Sample"
group.col <- "Group"
batch.col <- "Batch"

### Subsample targets per group

data.frame(table(cell.dat[[group.col]])) # Check number of cells per sample. ###its per group
##my results

##Var1   Freq
##1 14 dpi DAK 195839
##2 14 dpi PBS 101058
##3  7 dpi DAK 169957
##4  7 dpi PBS 166297

unique(cell.dat[[group.col]])

sub.targets <- c(50000,50000,50000,50000) 

##########################################################################################################
#### 5. Clustering and dimensionality reduction
##########################################################################################################

setwd(OutputDirectory)
dir.create("Output 2 - clustering")
setwd("Output 2 - clustering")

### Clustering

cell.dat <- run.flowsom(cell.dat, cluster.cols, meta.k = 25)
cell.dat

exp <- do.aggregate(cell.dat, cellular.cols, by = "FlowSOM_metacluster")
make.pheatmap(exp, "FlowSOM_metacluster", cluster.cols, dendrograms = "row",dendrograms.sort = FALSE,file.name = "Cluster_cols.png",plot.title = "Cluster_Cols")
make.pheatmap(exp, "FlowSOM_metacluster", cellular.cols, dendrograms = "row",file.name = "Cellular_cols.png",plot.title = "Cellular_Cols") 


### Dimensionality reduction

cell.sub <- do.subsample(cell.dat, sub.targets, group.col)
cell.sub


cell.sub <- run.tsne(cell.sub,
                     cluster.cols,
                     tsne.x.name = "tSNE_X",
                     tsne.y.name = "tSNE_Y",
                     perplexity = 30)

make.colour.plot(cell.sub, "tSNE_X", "tSNE_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE)
make.multi.plot(cell.sub, "tSNE_X", "tSNE_Y", cellular.cols)
make.multi.plot(cell.sub, "tSNE_X", "tSNE_Y", "FlowSOM_metacluster", group.col, col.type = 'factor')



##########################################################################################################
#### 6. Annotate clusters
##########################################################################################################

setwd(OutputDirectory)
dir.create("Output 3 - annotation")
setwd("Output 3 - annotation")

### Make specific annoations

annots <- list("B cells" = c(20),
               "CD8 T cells" = c(9,13,8,11),
               "CD4 T cells" = c(1),
               "Macrophages" = c(16,17,4,7,19,12,5,3),
               "DCs" = c(18,21,6,15),
               "ILCs" = c(22,14,24,23,2,10,25)
)

annots <- do.list.switch(annots)
names(annots) <- c("Values", "Population")
setorderv(annots, 'Values')
annots


### Add annotations

cell.dat <- do.add.cols(cell.dat, "FlowSOM_metacluster", annots, "Values")
cell.dat

cell.sub <- do.add.cols(cell.sub, "FlowSOM_metacluster", annots, "Values")
cell.sub


make.colour.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", col.type = 'factor', add.label = TRUE)
make.colour.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", col.type = 'factor', add.label = FALSE, title = "Decidua")
make.multi.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", group.col, col.type = 'factor')
make.multi.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", batch.col, col.type = 'factor')


### Expression heatmap

rm(exp)
exp <- do.aggregate(cell.dat, cellular.cols, by = "Population",)
make.pheatmap(exp, "Population", cellular.cols, dendrograms = "row", file.name = "Cellular_cols.png",plot.title = "Cellular_Cols" )

exp <- do.aggregate(cell.dat, cluster.cols, by = "Population",)
make.pheatmap(exp, "Population", cluster.cols, dendrograms = "row", file.name = "Cluster_cols.png",plot.title = "Cluster_Cols" )

#make.pheatmap(exp, "Population", cluster.cols, dendrograms = "row", cutree_rows = 3)


### Write FCS files

setwd(OutputDirectory)
setwd("Output 3 - annotation")

fwrite(cell.dat, "Annotated.data.csv")
fwrite(cell.sub, "Annotated.data.DR.csv")

dir.create('FCS files')
setwd('FCS files')

write.files(cell.dat,
            file.prefix = exp.name,
            divide.by = sample.col,
            write.csv = FALSE,
            write.fcs = TRUE)
write.files(cell.dat,
            file.prefix = exp.name,
            divide.by = sample.col,
            write.csv = TRUE,
            write.fcs = FALSE)
##########################################################################################################
#### 7. Summary data and statistical analysis
##########################################################################################################

setwd(OutputDirectory)
dir.create("Output 4 - summary data")
setwd("Output 4 - summary data")

### Setup

variance.test <- 'anova'
pairwise.test <- "t.test"

#Var1  Freq
#1 14 dpi DAK 71812
#2 14 dpi PBS 44743
#3  7 dpi DAK 69672
#4   7dpi PBS 72147  
comparisons <- list(c("14 dpi DAK", "14 dpi PBS"),
                    c("7 dpi DAK", "7 dpi PBS"),
                    c("14 dpi DAK", "7 dpi DAK"),
                    c("14 dpi PBS", "7 dpi PBS")
)
comparisons

grp.order <- c("7 dpi PBS", "7 dpi DAK", "14 dpi PBS","14 dpi DAK")
grp.order

### Select columns to measure MFI

as.matrix(cellular.cols)
MFI.cols <- cellular.cols[c(1:23)]
MFI.cols

### Create summary tables

sum.dat <- create.sumtable(dat = cell.dat, 
                           sample.col = sample.col,
                           pop.col = "Population",
                           use.cols = cellular.cols, 
                           annot.cols = c(group.col, batch.col), 
                           counts = counts)

### Review summary data

sum.dat
as.matrix(names(sum.dat))

annot.cols <- c(group.col, batch.col)

plotperct.cols <- names(sum.dat)[c(4:9)]
plotperct.cols

### Reorder summary data and SAVE

sum.dat <- do.reorder(sum.dat, group.col, grp.order)
sum.dat[,c(1:3)]

fwrite(sum.dat, 'sum.dat.csv')
#######heat map of summary data percentages
make.pheatmap(sum.dat, 
              sample.col = sample.col, 
              plot.cols = plotperct.cols, 
              plot.title = 'Frequency of Decidual Leukocytes',
              annot.cols = annot.cols,
              dendrograms = "both",
              normalise = FALSE)
make.pheatmap(sum.dat, 
              sample.col = sample.col, 
              plot.cols = plotperct.cols, 
              plot.title = 'Frequency of Decidual Leukocytes',
              annot.cols = annot.cols,
              dendrograms = "row",
              normalise = FALSE,
              file.name = 'other-colors_inferno_only_rows.png',
              standard.colours = "inferno",
              fold.colours = 'Spectre'
              )

rm(exp)
exp <- do.aggregate(cell.dat, cellular.cols, by = "Population",)
make.pheatmap(exp, "Population", cluster.cols, dendrograms = "row", file.name = "Cluster_cols_spectral.png",plot.title = "Cluster_Cols", standard.colours = 'Reds' )
#standard.colours = "BuPu",
#fold.colours = 'Spectre',
#file.name = paste0("Pheatmap by ", sample.col, ".png"),
#@param standard.colours DEFAULT = "BuPu". Character. Can also be "RdYlBu", "YlGnBu", "viridis", "magma", "inferno", "spectral", "Blues", "Reds", "Greys", or "rev(RdBu)".
#@param fold.colours DEFAULT = "Spectre". Character, a custom blue:black:red fold-change colour scheme. Can also be "RdYlBu", "YlGnBu", "viridis", "magma", "inferno", "spectral", "Blues", "Reds", "Greys", or "rev(RdBu)".

for(i in plotperct.cols){
  
  measure <- gsub("\\ --.*", "", i)
  measure
  
  pop <- gsub("^[^--]*.-- ", "", i)
  pop
  
  make.autograph(sum.dat,
                 x.axis = group.col,
                 y.axis = i,
                 y.axis.label = measure,
                 
                 grp.order = grp.order,
                 my_comparisons = comparisons,
                 
                 Variance_test = variance.test,
                 Pairwise_test = pairwise.test,
                 
                 title = pop,
                 subtitle = measure,
                 filename = paste0(i, '.pdf'))
  
}

setwd(OutputDirectory)
dir.create("Output - info", showWarnings = FALSE)
setwd("Output - info")

sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
session_info()
sink()

write(cellular.cols, "cellular.cols.txt")
write(cluster.cols, "cluster.cols.txt")
