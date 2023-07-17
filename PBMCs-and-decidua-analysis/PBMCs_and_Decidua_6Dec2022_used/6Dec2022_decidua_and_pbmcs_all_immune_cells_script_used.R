##########################################################################################################
#### 1. Load packages, and set working directory
##########################################################################################################

### Load libraries

library(Spectre)
Spectre::package.check()    # Check that all required packages are installed
Spectre::package.load()     # Load required packages

### Set PrimaryDirectory

setwd("~/Desktop/PBMCs_and_Decidua_6Dec2022_used")
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

### Set 'input' directory

setwd(PrimaryDirectory)
setwd("~/Desktop/PBMCs_and_Decidua_6Dec2022_used/data/")
InputDirectory <- getwd()
setwd(PrimaryDirectory)

### Set 'metadata' directory

setwd(PrimaryDirectory)
setwd("~/Desktop/PBMCs_and_Decidua_6Dec2022_used/metadata/")
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

data.list <- Spectre::read.files(file.loc = InputDirectory,
                                 file.type = ".csv",
                                 do.embed.file.names = TRUE)

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
##[3,] "DC-SIGN" 
##[4,] "HLA-DR"
##[5,] "CD56"    
##[6,] "CD45"    
##[7,] "CCR6"    
##[8,] "CD45RA"  
##[9,] "CD4"     
##[10,] "CD69"    
##[11,] "CD8a"    
##[12,] "FoxP3"   
##[13,] "CD16"    
##[14,] "CD14"    
##[15,] "CD49a"   
##[16,] "CD163"   
##[17,] "CD3"     
##[18,] "CCR7"    
##[19,] "RORyt"
##[20,] "T-bet"   
##[21,] "NKp46"   
##[22,] "CD11c"   
##[23,] "CD127"   
##[24,] "CD86"    
##[25,] "CD20"    
##[26,] "Eomes"   
##[27,] "Time"    
##[28,] "FileName"
##[29,] "FileNo"
## specify the columns to apply arcsinh transformation to
##[3,] "DC-SIGN" ##cf 12000
##[5,] "CD56"  ##cf 8000
##[19,] "RORyt"##cf 12000
as.matrix(names(cell.dat))
#Define the cofactor 
###set #1
to.asinh <- names(cell.dat)[c(4,7:18,20:26)] 
  
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

###set #2
to.asinh <- names(cell.dat)[c(3,19)] 
to.asinh
cofactor <- 12000

plot.against <- "CD45"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 20000), i, plot.against)
}
###set #3
to.asinh <- names(cell.dat)[c(5)] 
to.asinh
cofactor <- 8000
plot.against <- "CD45"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

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

cellular.cols <- names(cell.dat)[c(30:52)]
as.matrix(cellular.cols)
##[,1]           
###....     
##[29,] "FileNo"       
##[30,] "HLA-DR_asinh" 
##[31,] "CCR6_asinh"   
##[32,] "CD45RA_asinh" 
##[33,] "CD4_asinh"    
##[34,] "CD69_asinh"   
##[35,] "CD8a_asinh"   
##[36,] "FoxP3_asinh"  
##[37,] "CD16_asinh"   
##[38,] "CD14_asinh"   
##[39,] "CD49a_asinh"  
##[40,] "CD163_asinh"  
##[41,] "CD3_asinh"    
##[42,] "CCR7_asinh"   
##[43,] "T-bet_asinh"  
##[44,] "NKp46_asinh"  
##[45,] "CD11c_asinh"  
##[46,] "CD127_asinh"  
##[47,] "CD86_asinh"   
##[48,] "CD20_asinh"   
##[49,] "Eomes_asinh"  
##[50,] "DC-SIGN_asinh"
##[51,] "RORyt_asinh"  
##[52,] "CD56_asinh"   
##[53,] "Sample"       
##[54,] "Group"        
##[55,] "Batch"       
cluster.cols <- names(cell.dat)[c(30:33,35,37:46,48:50, 52)] 
as.matrix(cluster.cols)

exp.name <- "77-MRK_Decidua_and_PBMCs_6Dec2022"
sample.col <- "Sample"
group.col <- "Group"
batch.col <- "Batch"

### Subsample targets per group

data.frame(table(cell.dat[[group.col]])) # Check number of cells per sample. ###its per group
##my results
##Var1   Freq
##1 Decidua 633151
##2   PBMCs 746652

unique(cell.dat[[group.col]])

sub.targets <- c(50000,50000) 

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

annots <- list("B cells" = c(21,19),
               "CD8 T cells" = c(3,2),
               "CD4 T cells" = c(1,4),
               "DP T cells" = c(5),
               "DN T cells" = c(6),
               "Macrophages/Monocytes" = c(24,13,22,25,23,14),
               "ILCs" = c(10,11,9,12,16,8,7),
               "DCs" = c(20,15,17,18)
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
make.multi.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", group.col, col.type = 'factor')
make.colour.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", col.type = 'factor', add.label = FALSE, )


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

data.frame(table(cell.dat[[group.col]])) 

comparisons <- list(c("Decidua", "PBMCs")
)
comparisons

grp.order <- c("PBMCs", "Decidua")
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

plotperct.cols <- names(sum.dat)[c(4:11)]
plotperct.cols

### Reorder summary data and SAVE

sum.dat <- do.reorder(sum.dat, group.col, grp.order)
sum.dat[,c(1:3)]

fwrite(sum.dat, 'sum.dat.csv')
#######heat map of summary data percentages
make.pheatmap(sum.dat, 
              sample.col = sample.col, 
              plot.cols = plotperct.cols, 
              plot.title = 'Frequency of Major Immune Cells Populations in Decidua Samples',
              annot.cols = annot.cols,
              dendrograms = "both",
              normalise = FALSE)

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
