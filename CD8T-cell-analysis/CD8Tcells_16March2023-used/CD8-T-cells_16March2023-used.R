##########################################################################################################
#### 1. Load packages, and set working directory
##########################################################################################################

### Load libraries

library(Spectre)
Spectre::package.check()    # Check that all required packages are installed
Spectre::package.load()     # Load required packages

### Set PrimaryDirectory

setwd("~/Desktop/CD8Tcells_16March2023-used")
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

### Set 'input' directory

setwd(PrimaryDirectory)
setwd("~/Desktop/CD8Tcells_16March2023-used/data")
InputDirectory <- getwd()
setwd(PrimaryDirectory)

### Set 'metadata' directory

setwd(PrimaryDirectory)
setwd("~/Desktop/CD8Tcells_16March2023-used/metadata/")
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
#### 3. Data transformation- applying selected cofactors to parameters of choice
##########################################################################################################
setwd(OutputDirectory)
##### name your output for the co-factor you will try
dir.create("Output 1 - transformed plots")
setwd("Output 1 - transformed plots")


### Arcsinh transformation
##First, check the column names of the dataset.
as.matrix(names(cell.dat))
###My results and notes on what co-factor is going to be used for each 
##my results
##[,1]      
##[1,] "FSC-A"   not cellular
##[2,] "SSC-A"   not cellular
##[3,] "CCR6"   not cellular 
##[4,] "CCR7"    CF 8,000
##[5,] "CD127"   CF 3,000
##[6,] "CD4"     CF 5,000
##[7,] "CD45RA"  CF 2,000
##[8,] "CD49a"   CF 4,000
##[9,] "CD56"    CF 4,000
##[10,] "CD69"    CF 5,000
##[11,] "CD86"    CF 3,000
##[12,] "CD8a"    CF 4,000
##[13,] "Eomes"   CF 8,000
##[14,] "FoxP3"   not cellular
##[15,] "HLA-DR"  CF 5,000
##[16,] "T-bet"   CF 2,000
##[17,] "FileName"
##[18,] "FileNo" 


to.asinh <- names(cell.dat)[c(7,16)] 

to.asinh
cofactor <- 2000
plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 60000), i, plot.against, title = cofactor)
}

###Set 2 CF 3000
to.asinh <- names(cell.dat)[c(11,5)] 

to.asinh
cofactor <- 3000
plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 60000), i, plot.against, title = cofactor)
}
###Set 3 CF 4000
to.asinh <- names(cell.dat)[c(8,9,12)] 

to.asinh
cofactor <- 4000
plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 60000), i, plot.against, title = cofactor)
}

###Set 4 CF 5000
to.asinh <- names(cell.dat)[c(6,10,15)] 

to.asinh
cofactor <- 5000
plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 60000), i, plot.against, title = cofactor)
}
###Set 5 CF 8000
to.asinh <- names(cell.dat)[c(4,13)] 

to.asinh
cofactor <- 8000
plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")

for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 60000), i, plot.against, title = cofactor)
}
############################################################
### Add metadata to data.table ####This next set of code does not change until you get to the part where you are selecting your cellular cols

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
##My results
##[,1]         
##[1,] "FSC-A"       
##[2,] "SSC-A"       
##[3,] "CCR6"        
##[4,] "CCR7"        
##[5,] "CD127"       
##[6,] "CD4"         
##[7,] "CD45RA"      
##[8,] "CD49a"       
##[9,] "CD56"        
##[10,] "CD69"        
##[11,] "CD86"        
##[12,] "CD8a"        
##[13,] "Eomes"       
##[14,] "FoxP3"       
##[15,] "HLA-DR"      
##[16,] "T-bet"       
##[17,] "FileName"    
##[18,] "FileNo"      
##[19,] "CD45RA_asinh"
##[20,] "T-bet_asinh" 
##[21,] "CD86_asinh"  
##[22,] "CD127_asinh" 
##[23,] "CD49a_asinh" 
##[24,] "CD56_asinh"  
##[25,] "CD8a_asinh"  
##[26,] "CD4_asinh"   
##[27,] "CD69_asinh"  
##[28,] "HLA-DR_asinh"
##[29,] "CCR7_asinh"  
##[30,] "Eomes_asinh" 
##[31,] "Sample"      
##[32,] "Group"       
##[33,] "Batch"       
################################################################

################################################################

################################################################

###identify your cellular.cols, I choose the parameters I will want to visualize and or cluster with. Use your transformed cols for any parameter you transformed. 
cellular.cols <- names(cell.dat)[c(19:30)]
as.matrix(cellular.cols)
 ####### select the parameters you want to use for clustering. I exclude activation markers from this list (HLA-DR, CD69, CD86) I will also try clustering with or without CD45RA    
cluster.cols <- names(cell.dat)[c(20,22:26,29,30)] 
as.matrix(cluster.cols)

exp.name <- "77-MRK_CD8T-cells-global_16March2023"
sample.col <- "Sample"
group.col <- "Group"
batch.col <- "Batch"

### Subsample targets per group

data.frame(table(cell.dat[[group.col]])) # Check number of cells per sample. ###its per group
##my results ###with notes on 50% and 25% couts
##Var1 Freq
##1 14 dpi DAK 61062
##2 14 dpi PBS 29776
##3  7 dpi DAK 61537
##4  7 dpi PBS 60055


unique(cell.dat[[group.col]])
###The order of these results is the order that you should put in the next line of code
###my results: [1] "7 dpi DAK"  "7 dpi PBS"  "14 dpi DAK" "14 dpi PBS"
####You are selecting what you want to down sample to for tsne visualization per group. You can do this different ways. 
##### Including more than 100,000 total doesn't really do anything. Including less than 50,000 total may make it so you have to change
##### the perplexitiy to make it look nicer. increase perplex when you have small cell numbers
##### YOu can take a percentage of each group, however, I have found that it makes it difficult to visulize differences in the cluters between groups
##### If you take differnt amount from each group, Therefore I like taking the same amount for each group. There are good arguments for taking more from the groups that have more cells too. 
##### I will sometimes start with taking a percentage from each group when figureing out the parameters for the tsne map and what markers to use for clutering then use equal numbers once I have figured that out

#### I have tiny amounts of cells so I will be taking all but one
#### all cell counts of groups in order:  61537, 60055, 61062, 29776
sub.targets <- c(25000,25000,25000,25000) 

##########################################################################################################
#### 5. Clustering and dimensionality reduction
##########################################################################################################

setwd(OutputDirectory)
dir.create("Output 2 - clustering")
setwd("Output 2 - clustering")

### Clustering
#### The general advise is to over cluster and then reduce your cluters in the anotation step. I find that in gneral I have re-do the clustering processes a few times to select the markers, clusters, and other setting that work best
#### I find that looking at the flow data to guess how many cluters there are and then inbcluding 2-5 extra clusters is a good place to start. The tsne will help you too. 
cell.dat <- run.flowsom(cell.dat, cluster.cols, meta.k = 18)
cell.dat
### Dimensionality reduction

cell.sub <- do.subsample(cell.dat, sub.targets, group.col)
cell.sub


cell.sub <- run.tsne(cell.sub,
                     cluster.cols,
                     tsne.x.name = "tSNE_X",
                     tsne.y.name = "tSNE_Y",
                     perplexity = 50,
                     max_iter = 1000)

make.colour.plot(cell.sub, "tSNE_X", "tSNE_Y", "FlowSOM_metacluster", col.type = 'factor', add.label = TRUE)
make.multi.plot(cell.sub, "tSNE_X", "tSNE_Y", cellular.cols)
make.multi.plot(cell.sub, "tSNE_X", "tSNE_Y", "FlowSOM_metacluster", group.col, col.type = 'factor')
exp <- do.aggregate(cell.dat, cellular.cols, by = "FlowSOM_metacluster")
make.pheatmap(exp, "FlowSOM_metacluster", cluster.cols, dendrograms = "row",dendrograms.sort = FALSE,file.name = "Cluster_cols.png",plot.title = "Cluster_Cols")
make.pheatmap(exp, "FlowSOM_metacluster", cellular.cols, dendrograms = "row",file.name = "Cellular_cols.png",plot.title = "Cellular_Cols") 






##########################################################################################################
#### 6. Annotate clusters
##########################################################################################################

setwd(OutputDirectory)
dir.create("Output 3 - annotation")
setwd("Output 3 - annotation")

### Make specific annoations
####Come back to the step below if you need to fix your annotations, if you have not 'added your annotations' yet

#annots <- list("CD8_like_CD4+CD8-_T_cells_" = c(18),
#               "T-bet++_Eomes-low_CD8_T_cells" = c(10),
 #              "T-bet-low_Eomes+_CD8_T_cells" = c(4,5),
 #              "CD56+_CD8_T_cells" = c(6,7,12,13,9),
 #              "CCR7+_CD8_T_cells" = c(8),
 #              "CD127+_CD56+_CD8_T_cells" = c(11),
 #              "CD8_like_Double_Negative_T_cells" = c(14,15,16),
 #              "CD8_like_Double_Postitive_T_cells" = c(17),
#               "CD127+_CD56-_CD8_T_cells" = c(1,2),
#               "CD49+CD56-_CD8_T_cells" = c(3)
#)
annots <- list("CD8-+CD4+_T_cells" = c(18),
               "Double-negative-T-cells" = c(14,15,16),
               "Double-postitive-T-cells" = c(17),
               "CD8T-C1" = c(10),
               "CD8T-C2" = c(4,5),
               "CD8T-C3" = c(6,7,12,13,9),
               "CD8T-C4" = c(8),
               "CD8T-C5" = c(11),
               "CD8T-C6" = c(1,2),
               "CD8T-C7" = c(3)
)
annots <- do.list.switch(annots)
names(annots) <- c("Values", "Population")
setorderv(annots, 'Values')
####### Check that your annotations are correct and that you did not assign any cluster to more than one annotated population or forget to assign a cluster to a population
annots
#####
###Example of results that are done correctly
##The numbers in the first two columns should match up


######If there are mistakes go back to annotations and fix them before proceeding. You do not need to start over again. 
### Add annotations

cell.dat <- do.add.cols(cell.dat, "FlowSOM_metacluster", annots, "Values")
cell.dat

cell.sub <- do.add.cols(cell.sub, "FlowSOM_metacluster", annots, "Values")
cell.sub


make.colour.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", col.type = 'factor', add.label = TRUE)
make.multi.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", group.col, col.type = 'factor')
make.multi.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", batch.col, col.type = 'factor')
make.colour.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", col.type = 'factor', add.label = FALSE, title = "CD8 T cells")

make.multi.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", group.col, col.type = 'factor')

### Expression heatmap

rm(exp)
exp <- do.aggregate(cell.dat, cellular.cols, by = "Population",)
make.pheatmap(exp, "Population", cellular.cols, dendrograms = "row", file.name = "Cellular_cols.png",plot.title = "Cellular_Cols" )

exp <- do.aggregate(cell.dat, cluster.cols, by = "Population",)
make.pheatmap(exp, "Population", cluster.cols, dendrograms = "row", file.name = "Cluster_cols.png",plot.title = "Cluster_Cols" )




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

#Setup

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
MFI.cols <- cellular.cols[c(1:12)]
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

plotperct.cols <- names(sum.dat)[c(4:13)]
plotperct.cols

### Reorder summary data and SAVE

sum.dat <- do.reorder(sum.dat, group.col, grp.order)
sum.dat[,c(1:3)]

fwrite(sum.dat, 'sum.dat.csv')
#######heat map of summary data percentages
make.pheatmap(sum.dat, 
              sample.col = sample.col, 
              plot.cols = plotperct.cols, 
              plot.title = 'Frequency of CD8 T cell Pops.,Decidua',
              annot.cols = annot.cols,
              dendrograms = "both",
              file.name = "Frequency of CD8 T cell Pops normalized to column.png")
make.pheatmap(sum.dat, 
              sample.col = sample.col, 
              plot.cols = plotperct.cols, 
              plot.title = 'Frequency of CD8 T cell Pops, Decidua',
              annot.cols = annot.cols,
              dendrograms = "both",
              normalise = FALSE,
              file.name = "Frequency of CD8 T cell Pops pops.png")
#make.pheatmap(sum.dat, 
 #            sample.col = sample.col, 
 #            plot.cols = plotperct.cols, 
#             plot.title = 'Frequency of DC Pops, Decidua',
#             annot.cols = annot.cols,
 #            dendrograms = "both",
 #            normalise = FALSE,
  #           cutree_rows = 2,
 #            file.name = "Frequency of DC Pops cutree.png")


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

