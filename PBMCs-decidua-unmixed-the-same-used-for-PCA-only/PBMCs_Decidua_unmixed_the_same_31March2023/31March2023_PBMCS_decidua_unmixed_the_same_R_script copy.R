##########################################################################################################
#### 1. Load packages, and set working directory
##########################################################################################################

### Load libraries

library(Spectre)
Spectre::package.check()    # Check that all required packages are installed
Spectre::package.load()     # Load required packages

### Set PrimaryDirectory

setwd("~/Desktop/PBMCs_Decidua_unmixed_the_same_31March2023")
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

### Set 'input' directory

setwd(PrimaryDirectory)
setwd("~/Desktop/PBMCs_Decidua_unmixed_the_same_31March2023/data/")
InputDirectory <- getwd()
setwd(PrimaryDirectory)

### Set 'metadata' directory

setwd(PrimaryDirectory)
setwd("~/Desktop/PBMCs_Decidua_unmixed_the_same_31March2023/metadata/")
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
##my results
##[1,] "FSC-A"   
##[2,] "SSC-A"   
##[3,] "DC-SIGN" 
##[4,] "HLA-DR"  
##[5,] "CD56"    
##[6,] "CCR6"    
##[7,] "CD45RA"  
##[8,] "CD4"     
##[9,] "CD69"    
##[10,] "CD8a"    
##[11,] "FoxP3"   
##[12,] "CD16"    
##[13,] "CD14"    
##[14,] "CD49a"   
##[15,] "CD163"   
##[16,] "CD3"     
##[17,] "CCR7"    
##[18,] "RORyt"   
##[19,] "T-bet"   
##[20,] "NKp46"   
##[21,] "CD11c"   
##[22,] "CD127"   
##[23,] "CD86"    
##[24,] "CD20"    
##[25,] "Eomes"   
##[26,] "FileName"
##[27,] "FileNo"  
as.matrix(names(cell.dat))
#Define the cofactor 
###set #1
to.asinh <- names(cell.dat)[c(4,6:17,19:25)] 
#hla-dr, ccr6, ra, 4, 69, 8, fox, 16, 49a, 163, 3, ccr7 t-bet, nkp46, 11c, 127, 86, eomes
to.asinh
#In this worked example we will use a cofactor of 500, the recommended vlaue to start with for flow data
cofactor <- 5000
##You can also choose a column to use for plotting the transformed result – ideally something that is expressed on a variety of cell types in your dataset.

plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 20000), i, plot.against)
}

###set #2
to.asinh <- names(cell.dat)[c(3,18)] 
#dc-sign, roryt
to.asinh
cofactor <- 12000

plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 20000), i, plot.against)
}
###set #3
to.asinh <- names(cell.dat)[c(5)] 
to.asinh
#CD56
cofactor <- 8000
plot.against <- "FSC-A"

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
#[,1] 
#......
#[26,] "FileName"     
#[27,] "FileNo"       
#[28,] "HLA-DR_asinh" 
#[29,] "CCR6_asinh"   
#[30,] "CD45RA_asinh" 
#[31,] "CD4_asinh"    
#[32,] "CD69_asinh"   
#[33,] "CD8a_asinh"   
#[34,] "FoxP3_asinh"  
#[35,] "CD16_asinh"   
#[36,] "CD14_asinh"   
#[37,] "CD49a_asinh"  
#[38,] "CD163_asinh"  
#[39,] "CD3_asinh"    
#[40,] "CCR7_asinh"   
#[41,] "T-bet_asinh"  
#[42,] "NKp46_asinh"  
#[43,] "CD11c_asinh"  
#[44,] "CD127_asinh"  
#[45,] "CD86_asinh"   
#[46,] "CD20_asinh"   
#[47,] "Eomes_asinh"  
#[48,] "DC-SIGN_asinh"
#[49,] "RORyt_asinh"  
#[50,] "CD56_asinh"   
#[51,] "Sample"       
#[52,] "Group"        
#[53,] "Batch"      
cellular.cols <- names(cell.dat)[c(28:50)]
as.matrix(cellular.cols)

cluster.cols <- names(cell.dat)[c(28:31,33,35:44,46:48,50)] 
as.matrix(cluster.cols)

exp.name <- "77-MRK_Decidua_and_PBMCs_both_unmixed_like_decidua_31March2023"
sample.col <- "Sample"
group.col <- "Group"
batch.col <- "Batch"

### Subsample targets per group

data.frame(table(cell.dat[[group.col]])) # Check number of cells per sample. ###its per group
#Var1   Freq
#1 DAK Decidua 364763
#2    DAK PBMC 379271
#3 PBS Decidua 263378
#4   PBS PBMCs 342055
#data.frame(table(cell.dat[[batch.col]]))



unique(cell.dat[[group.col]])

sub.targets <- c(30000,30000,30000,30000) 

##########################################################################################################
#### 5. Clustering and dimensionality reduction
##########################################################################################################

setwd(OutputDirectory)
dir.create("Output 2 - clustering")
setwd("Output 2 - clustering")

### Clustering

cell.dat <- run.flowsom(cell.dat, cluster.cols, meta.k = 30)
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

annots <- list("CD4 T cells" = c(24),
               "CD8 T cells" = c(19,21,29),
               "DP T cells" = c(23),
               "DN T cells" = c(22,30),
               "B cells" = c(28,27,25),
               "ILCs" = c(1,13,14,15,16,26),
               "Macrophages-Monocytes" = c(3,6,7,8,9,10,17,18,20),
               "DCs" = c(2,4,5,11,12)
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
make.colour.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", col.type = 'factor', add.label = FALSE, title = "PBMCs and Decidua")


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

grp.order <- c("PBS PBMCs", "DAK PBMC", "PBS Decidua", "DAK Decidua")
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
            plot.title = 'Frequency of Major Immune Cells Populations',
           annot.cols = annot.cols,
            dendrograms = "both",
            normalise = FALSE)

#for(i in plotperct.cols){

#  measure <- gsub("\\ --.*", "", i)
#  measure

#  pop <- gsub("^[^--]*.-- ", "", i)
#  pop

#  make.autograph(sum.dat,
#               x.axis = group.col,
#               y.axis = i,
#               y.axis.label = measure,
#
#              grp.order = grp.order,
#               my_comparisons = comparisons,
#
#             Variance_test = variance.test,
#            Pairwise_test = pairwise.test,
#              
#            title = pop,
#            subtitle = measure,
#            filename = paste0(i, '.pdf'))
#
#}

setwd(OutputDirectory)
dir.create("Output - info", showWarnings = FALSE)
setwd("Output - info")

sink(file = "session_info.txt", append=TRUE, split=FALSE, type = c("output", "message"))
session_info()
sink()

write(cellular.cols, "cellular.cols.txt")
write(cluster.cols, "cluster.cols.txt")
