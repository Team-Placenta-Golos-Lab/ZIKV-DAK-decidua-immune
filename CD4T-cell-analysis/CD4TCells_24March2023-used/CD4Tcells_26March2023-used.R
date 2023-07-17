#########################################################################################################
#### 1. Load packages, and set working directory
##########################################################################################################

### Load libraries

library(Spectre)
Spectre::package.check()    # Check that all required packages are installed
Spectre::package.load()     # Load required packages

### Set PrimaryDirectory

setwd("~/Desktop/CD4TCells_24March2023-used")
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

### Set 'input' directory

setwd(PrimaryDirectory)
setwd("~/Desktop/CD4TCells_24March2023-used/data/")
InputDirectory <- getwd()
setwd(PrimaryDirectory)

### Set 'metadata' directory

setwd(PrimaryDirectory)
setwd("~/Desktop/CD4TCells_24March2023-used/metadata/")
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
#### 3. Data transformation- testing different cofactors
##########################################################################################################
setwd(OutputDirectory)
##### name your output for the co-factor you will try
dir.create("Output 1 - transformed plots")
setwd("Output 1 - transformed plots")

### Arcsinh transformation
##First, check the column names of the dataset.
as.matrix(names(cell.dat))
##my results
##[,1]      
##[1,] "FSC-A"   
##[2,] "SSC-A"   
#[3,] "CCR6"    
#[4,] "CCR7"    
#[5,] "CD11c"   
#[6,] "CD127"   
#[7,] "CD4"     
#[8,] "CD45RA"  
#[9,] "CD49a"   
#[10,] "CD56"    
#[11,] "CD69"    
#[12,] "CD86"    
#[13,] "CD8a"    
#[14,] "Eomes"   
#[15,] "FoxP3"   
#[16,] "HLA-DR"  
#[17,] "RORyt"   
#[18,] "T-bet"   
#[19,] "FileName"
#[20,] "FileNo" 
###select the parameters you want to transform, I select all the cellular columns but FSC & SSC for testing

###set 1
to.asinh <- names(cell.dat)[c(8,15)] 

to.asinh

cofactor <- 1500
##You can also choose a column to use for plotting the transformed result – ideally something that is expressed on a variety of cell types in your dataset.

plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")
for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 30000), i, plot.against, title = cofactor)
}
###set 2
to.asinh <- names(cell.dat)[c(7,10,16)] 

to.asinh

cofactor <- 2000
##You can also choose a column to use for plotting the transformed result – ideally something that is expressed on a variety of cell types in your dataset.

plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")
for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 30000), i, plot.against, title = cofactor)
}
###set 3
to.asinh <- names(cell.dat)[c(18)] 

to.asinh

cofactor <- 2500

plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")
for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 30000), i, plot.against, title = cofactor)
}
###set 4
to.asinh <- names(cell.dat)[c(12)] 

to.asinh

cofactor <- 3000

plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")
for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 30000), i, plot.against, title = cofactor)
}
###set 5
to.asinh <- names(cell.dat)[c(3,5,6,9,11,13)] 

to.asinh

cofactor <- 4000

plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")
for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 30000), i, plot.against, title = cofactor)
}
###set 6
to.asinh <- names(cell.dat)[c(4)] 

to.asinh

cofactor <- 4500

plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")
for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 30000), i, plot.against, title = cofactor)
}
###set 7
to.asinh <- names(cell.dat)[c(14,17)] 

to.asinh

cofactor <- 6000

plot.against <- "FSC-A"

cell.dat <- do.asinh(cell.dat, to.asinh, cofactor = cofactor)
transformed.cols <- paste0(to.asinh, "_asinh")
for(i in transformed.cols){
  make.colour.plot(do.subsample(cell.dat, 30000), i, plot.against, title = cofactor)
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
##[5,] "CD11c"       
##[6,] "CD127"       
##[7,] "CD4"         
##[8,] "CD45RA"      
##[9,] "CD49a"       
##[10,] "CD56"        
##[11,] "CD69"        
##[12,] "CD86"        
##[13,] "CD8a"        
##[14,] "Eomes"       
##[15,] "FoxP3"       
##[16,] "HLA-DR"      
##[17,] "RORyt"       
##[18,] "T-bet"       
##[19,] "FileName"    
##[20,] "FileNo"      
##[21,] "CD45RA_asinh"
##[22,] "FoxP3_asinh" 
##[23,] "CD4_asinh"   
##[24,] "CD56_asinh"  
##[25,] "HLA-DR_asinh"
##[26,] "T-bet_asinh" 
##[27,] "CD86_asinh"  
##[28,] "CCR6_asinh"  
##[29,] "CD11c_asinh" 
##[30,] "CD127_asinh" 
##[31,] "CD49a_asinh" 
##[32,] "CD69_asinh"  
##[33,] "CD8a_asinh"  
##[34,] "CCR7_asinh"  
##[35,] "Eomes_asinh" 
##[36,] "RORyt_asinh" 
##[37,] "Sample"      
##[38,] "Group"       
##[39,] "Batch"      
################################################################
###identify your cellular.cols, I choose the parameters I will want to visualize and or cluster with. Use your transformed cols for any parameter you transformed. 
cellular.cols <- names(cell.dat)[c(21:36)]
as.matrix(cellular.cols)
####### select the parameters you want to use for clustering. I exclude activation markers from this list (HLA-DR, CD69, CD86) I will also try clustering with or without CD45RA    
cluster.cols <- names(cell.dat)[c(22:24,26,28:31,33,34,35,36)] 
as.matrix(cluster.cols)

exp.name <- "77-MRK_CD4Tcells_26March2023"
sample.col <- "Sample"
group.col <- "Group"
batch.col <- "Batch"

### Subsample targets per group

data.frame(table(cell.dat[[group.col]])) # Check number of cells per sample. ###its per group
##my results ###with notes on 50% and 25% couts
#Var1  Freq
#1 14 dpi DAK 12058
#2 14 dpi PBS 15927
#3  7 dpi DAK  9806
#4  7 dpi PBS 13861


unique(cell.dat[[group.col]])
###The order of these results is the order that you should put in the next line of code
###my results: [1] "7 dpi DAK"  "7 dpi PBS"  "14 dpi DAK" "14 dpi PBS"
####You are selecting what you want to down sample to for tsne visualization per group. You can do this different ways. 
##### Including more than 100,000 total doesn't really do anything. Including less than 50,000 total may make it so you have to change
##### the perplexitiy to make it look nicer. increase perplex when you have small cell numbers
##### YOu can take a percentage of each group, however, I have found that it makes it difficult to visulize differences in the cluters between groups
##### If you take differnt amount from each group, Therefore I like taking the same amount for each group. There are good arguments for taking more from the groups that have more cells too. 
##### I will sometimes start with taking a percentage from each group when figureing out the parameters for the tsne map and what markers to use for clutering then use equal numbers once I have figured that out

#### If I have tiny amounts of cells so I will be taking all but one
#### all cell counts of groups in order:  9806, 13861, 12058, 15927
sub.targets <- c(9805,13860,12057,15926) 

##########################################################################################################
#### 5. Clustering and dimensionality reduction
##########################################################################################################

setwd(OutputDirectory)
dir.create("Output 2 - clustering")
setwd("Output 2 - clustering")

### Clustering
#### The general advise is to over cluster and then reduce your cluters in the anotation step. I find that in gneral I have re-do the clustering processes a few times to select the markers, clusters, and other setting that work best
#### I find that looking at the flow data to guess how many cluters there are and then inbcluding 2-5 extra clusters is a good place to start. The tsne will help you too. 
cell.dat <- run.flowsom(cell.dat, cluster.cols, meta.k = 16)
cell.dat
### Dimensionality reduction

cell.sub <- do.subsample(cell.dat, sub.targets, group.col)
cell.sub


cell.sub <- run.tsne(cell.sub,
                     cluster.cols,
                     tsne.x.name = "tSNE_X",
                     tsne.y.name = "tSNE_Y",
                     perplexity = 100,
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

#annots <- list("CD4T-1" = c(9), #Tbet+CD49a-
#               "CD4T-2" = c(5,11,12), #Tbet- CD127+
#               "CD4T-3" = c(4), #Tbet+ CD49a+
#               "CD4T-4" = c(14), #Tregs
#               "CD4T-5" = c(1,3,8), #CD56+ 
#               "CD4T-6" = c(13), #CD11c+
#              "CD4T_like_Double-Negative_Tcells" = c(16,15),
#              "CD4T_like_Double-Positive_Tcells" = c(2,6,7,10)
#)
#annots <- list("CD4T-1" = c(9), #Tbet+CD49a-
#               "CD4T-2" = c(5,11,12), #Tbet- CD127+
#               "CD4T-3" = c(4), #Tbet+ CD49a+
#               "CD4T-4" = c(14), #Tregs
#               "CD4T-5" = c(1,3,8), #CD56+ 
#               "CD4T-6" = c(13), #CD11c+
#              "CD4T_like_Double-Negative_Tcells" = c(16,15),
#              "CD4T_like_Double-Positive_Tcells" = c(2,6,7,10)
#)
annots <- list("Th1-C1" = c(9), #Tbet+CD49a-
               "Naive-C2" = c(5,11,12), #Tbet- CD127+
               "Th1-C2" = c(4), #Tbet+ CD49a+
               "Tregs" = c(14), #Tregs
               "Th1-C3" = c(1,3,8), #CD56+ 
               "Naive-C1" = c(13), #CD11c+
               "Double-negative_T-cells" = c(16,15),
               "Double-positive-T-cells" = c(2,6,7,10)
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
make.colour.plot(cell.sub, "tSNE_X", "tSNE_Y", "Population", col.type = 'factor', add.label = FALSE, title = "CD4 T cells")

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
MFI.cols <- cellular.cols[c(1:16)]
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
             plot.title = 'Frequency of CD4T Pops.,Decidua',
             annot.cols = annot.cols,
              dendrograms = "both",
             file.name = "Frequency of CD4T Pops normalized to column.png")
make.pheatmap(sum.dat, 
             sample.col = sample.col, 
             plot.cols = plotperct.cols, 
             plot.title = 'Frequency of CD4T cell Pops, Decidua',
             annot.cols = annot.cols,
             dendrograms = "both",
             normalise = FALSE,
             file.name = "Frequency of CD4T cell Pops pops.png")
#make.pheatmap(sum.dat, 
#            sample.col = sample.col, 
#            plot.cols = plotperct.cols, 
#           plot.title = 'Frequency of DC Pops, Decidua',
#            annot.cols = annot.cols,
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
