#-----------------------------#
####     Setup options     ####
#-----------------------------#
###### Restart R #####
## This helps with certain memory related issues
## You only need to run this once before starting
source("https://utexas.box.com/shared/static/reml2t600ugc7ub3ogo69qhg8zrp32px.r")
cleanEnv()

###### Setup options you will likely change ######
## Create the options object
opt<-list()

## Setup files from filtering via GDC website
# My files are hosted in a way to let you easily use them as examples
# You will likely want to comment out the URL based download steps,
#    and just modify the File names / paths to match your data
opt$project    <- "TCGA-STAD"
opt$sample.type<- "Primary Tumor"
opt$group1Desc <- "MUC16 mutations with AJCC Stage < iii"
opt$group2Desc <- "MUC16 mutations with AJCC Stage >= iii"
opt$grp1File<- "less_than_iii.tsv"
opt$grp2File<- "greater_than_iii.tsv"



###### Setup options you will likely change if you are not on the Edupod ######
## Setup package and GDC directory options
opt$wd <- path.expand(file.path("~","bio321g","DiffExp"))
opt$pckDir  <- c(
  file.path(opt$wd,"R_gdc"),           #This is where any subsequent packages you install will end up
  "/stor/scratch/Bio321G_NB_Fall2023/R" #This is my package library made for this assignment with (largely) compatible software versions.
)
opt$gdcPath <- "/stor/scratch/Bio321G_NB_Fall2023/GDCdata" #We can all download GDC files to this file. Please don't put other files here.



###### Setup options you will NOT likely change #####
## Options for adjusting the package and gdc paths when needed
opt$createPckDir <- T # This toggles the if statement that adjusts the package directory
opt$createGdcDir <- T # This toggles the if statement that adjusts the GDC data directory

## Package options
# These are the minimum packages you need to run this script file
# These will not include my exact script for creating volcano and box plots
opt$biocPackages <- c(
  "ggplot2","ggrepel","TCGAbiolinks","DESeq2"
)



##### Other options that may or may not change ####
## Set query filtering options
opt$minPerGroup <- 12           # Used to control warning about small sample size
opt$maxPerGroup <- 50           # Used to control warning about large sample size
opt$sampleToMaxLogic <- T       # Used to control whether query is sampled to a smaller size when large
opt$sampleToMax.SeedLogic <- T  # Used to control whether seed is set prior to sampling
opt$sampleToMax.Seed <- 10112023# Used to specify the seed for random number generation

## Set cutoffs for determining significance
opt$fdr.max <- 0.01 ## Maximum false discovery rate adjusted p-value cutoff to be a DEG
opt$afc.min <- 1.50 ## Minimum absolute(Log2 fold change) cutoff to be considered DEG

## For visualization
opt$visScriptFile <- "deseq2visPlotFuns.R" # You won't have this file as it is part of your homework
opt$nPerSide      <- 15                   # How many significant gene names to print per x-axis side. Used in the visScriptFile.



#-----------------------------#
####   Setup environment   ####
#-----------------------------#
##### Setup location to store packages #####
## Setup the working directory
if(!dir.exists(opt$wd)){
  dir.create(opt$wd,recursive = T,showWarnings = T)
}
setwd(opt$wd)



##### Install packages #####
## Adjust package directory if needed
if(!dir.exists(opt$pckDir[1])){
  message("Creating home package directory.")
  dir.create(file.path(opt$wd,"R_gdc"),recursive = T,showWarnings = T)
}
if(!all(dir.exists(opt$pckDir))){
  message("Changing package directory to:")
  opt$pckDir <- file.path(opt$wd,"R_gdc")
  message(opt$pckDir)
}

## Tell R where to find packages
.libPaths(opt$pckDir)
message(paste0(.libPaths(),sep = "\n"))

## Install the packages
# Depending on your system this may take more than an hour & ~0.55 Gb of storage...
# However, I have already done this on the EduPod.
# That said, stuff often goes wrong here... Common things:
## You likely want to say no to restarting R if asked.
## You likely want to say no to "compile from source" if asked.
# Many other steps may go wrong; you may need to do additional troubleshooting / seek help
# Check if it succeeded by loading the libraries one at a time.
if(!all(opt$biocPackages%in%installed.packages())){
  ## Install an installer tool
  # Installation goes into custom library path due to file permission issues with some packages
  if(!"BiocManager"%in%installed.packages()){
    install.packages("BiocManager", lib = opt$pckDir)
  }
  
  ## Update the installed packages
  update.packages(instlib = opt$pckDir,ask = F)
  
  ## Install on modern version of R
  # Windows users will need to install RTools program first (not a package)
  BiocManager::install(
    opt$biocPackages, lib = opt$pckDir,
    ask = F,update = T
  )
}



##### Load packages #####
# Don't forget to change your package dir using .libPaths()
# Just check the following for warnings like the following:
#    Ex: Error in library(gibberish) : there is no package called ‘gibberish’
# Most other warning messages mean the package was installed, but it gave feedback.
if(!all(opt$pckDir%in%.libPaths())){
  print("Setting lib paths...")
  .libPaths(opt$pckDir)
}
library(ggplot2)
library(ggrepel)
library(TCGAbiolinks) #Takes a bit to load
library(DESeq2)       #Produces lots of "warnings"





#---------------------------#
#### Download TCGA Files ####
#---------------------------#
###### Search for TCGA Files ######
#Read more about this using browseVignettes("TCGAbiolinks")
query1 <- GDCquery(
  project = opt$project,
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  data.type = "Gene Expression Quantification"
)

###### Pull out targetted IDs from initial search ######
# Pull sample specific data out of complicated query object
samDf <- query1$results[[1]]
# Did this work and do you have the expected number of sample records?
#   For me, the following indicated 166 Primary Tumor sample records
table(samDf$sample_type)


###### Load the IDs identified on the GDC ######
## This requires you to have downloaded sample sheets from GDC
## I have used links to my examples using TCGA-READ samples with mutation data
##   collected into groups of with and without APC mutations.
## You will likely skip this step entirely for your data!
# download.file(url = opt$grp1Url,destfile = opt$grp1File)
# download.file(url = opt$grp2Url,destfile = opt$grp2File)

## Once these files (or your own!) are in the working directory, you can load them into R
group1 <- read.table(opt$grp1File,fill=T,header = 1)
group2 <- read.table(opt$grp2File,fill=T,header = 1)

## Check group size and if the exclusion analysis worked
if(nrow(group1)<opt$minPerGroup){warning("Group1: You should have a larger number of samples for this kind of analysis")}
if(nrow(group2)<opt$minPerGroup){warning("Group2: You should have a larger number of samples for this kind of analysis")}

if(nrow(group1)>opt$maxPerGroup){warning("Group1: You should have a smaller number of samples if working on the edupod")}
if(nrow(group2)>opt$maxPerGroup){warning("Group2: You should have a smaller number of samples if working on the edupod")}

if(!all(!group1$ID.1%in%group2$ID.1)&all(!group2$ID.1%in%group1$ID.1)){
  stop("The group IDs were not exclusive...")
}



###### Compare the browser identified IDs to the & GDC IDs ######
## Students commonly mess up the previous section, so I added the following check:
# You will need to adjust your files / code if this does not make sense!
# Your IDs should print off and represent sample ids present in samDf$sample.submitter_id
if(!"ID.1"%in%colnames(group1)|!"ID.1"%in%colnames(group2)){
  stop("Your column names in group1 and/or group2 are not as expected!")
}
message("The head and tail of your group1 and group2 IDs are as follows:")
head(group1$ID.1)
tail(group1$ID.1)
head(group2$ID.1)
tail(group2$ID.1)
message("These should represent the sample submitter ids present in samDf:")
head(samDf$sample.submitter_id)

## Subset samDf to just those rows with "sample.submitter_id" values present in groupIds
groupIds<- c(group1$ID.1,group2$ID.1)
samDf   <- samDf[samDf$sample.submitter_id%in%groupIds,]

## Quick check to head off common issues here.
if(nrow(samDf)==0){stop("Something went wrong when comparing IDs!")}



###### Sub-sample inputs if query is too large for EduPod ######
if(opt$sampleToMaxLogic){
  samDf.1 <- samDf[samDf$sample.submitter_id%in%group1$ID.1,]
  samDf.2 <- samDf[samDf$sample.submitter_id%in%group2$ID.1,]
  if(nrow(samDf.1)>opt$maxPerGroup){
    if(opt$sampleToMax.SeedLogic){set.seed(opt$sampleToMax.Seed)}
    samDf.1 <- samDf.1[sample(1:nrow(samDf.1),size = opt$maxPerGroup),]
  }
  if(nrow(samDf.2)>opt$maxPerGroup){
    if(opt$sampleToMax.SeedLogic){set.seed(opt$sampleToMax.Seed*2+1)} # Just avoiding using numbers I am likely to change opt$sampleToMax.Seed to
    samDf.2 <- samDf.2[sample(1:nrow(samDf.2),size = opt$maxPerGroup),]
  }
  samDf<-rbind(samDf.1,samDf.2)
}



###### Isolate barcodes corresponding to GDC data ######
desiredBarcodes <- samDf$cases
if(length(desiredBarcodes)==0){stop("Something went wrong getting barcodes!")}



###### Search for just these barcodes ######
query2 <- GDCquery(
  project = opt$project,
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  data.type = "Gene Expression Quantification",
  barcode = desiredBarcodes,
  sample.type = opt$sample.type
)



##### Download and prepare the data ######
## Adjust the storage directory if needed
if(!file.exists(opt$gdcPath)){
  if(opt$createGdcDir){
    warning("Creating directory to store GDC data in the wd.")
    opt$gdcPath <- "./GDCdata/"
    dir.create(opt$gdcPath,recursive = T,showWarnings = T)
  }else{
    stop("Oh, no! The file you are trying to store GDCdata in doesn't exist.")
  }
}

## Download data
# This commonly fails due to server
# This may require retrying at a later time
# For my data, this downloads about 322 Mb
# However, I have shared a directory with you as well.
GDCdownload(query2,method = "api",
            files.per.chunk = 10,
            directory = opt$gdcPath)



###### Prepare the downloaded data ######
# This step downloads and structures data into a "SummarizedExperiment"
# This is a class of data that organizes complicated multi-dimensional data
# See the notes in the next section for more details
dds <- GDCprepare(query = query2,directory = opt$gdcPath)

## Check if it worked
# && only checks the right side if the left side was TRUE
ddsLogic <-
  exists("dds")&&
  class(dds)=="RangedSummarizedExperiment"&&
  nrow(as.data.frame(colData(dds)))>0
if(!ddsLogic){stop("Something went wrong in GDCprepare()!")}






#---------------------------------------------#
#### Filter the samples and loci analyzed #####
#---------------------------------------------#
###### Notes on exploring the summarized experiment object
## Visualize **some** of the clinical data
View(as.data.frame(colData(dds)))


## Visualize the genomic data
# Organized into rows of "genomic features".
# Look at the column "gene_type" to see why I don't just say genes:
# gene_type defintions: https://www.gencodegenes.org/pages/biotypes.html
View(as.data.frame(rowRanges(dds)))

rowRanges_df <- as.data.frame(rowRanges(dds))

# Write to csv file 
write.csv(rowRanges_df, file = "rowRanges.csv", row.names = FALSE)

## See what data has been added to the observation section
#### These are the various downloaded datasets
names(assays(dds))
#### For my data:
##[1] "unstranded"      #The raw count data for paired-end reads, the one DESeq2 uses
##[2] "stranded_first"  #The raw count data for first read
##[3] "stranded_second" #The raw count data for second read
##[4] "tpm_unstrand"    #The transcripts per million. "for every 1,000,000 RNA molecules in the RNA-seq sample, x came from this gene/transcript." https://wiki.arrayserver.com/wiki/index.php?title=TPM
##[5] "fpkm_unstrand"   #The fragments per kilobase of transcript per million fragments mapped
##[6] "fpkm_uq_unstrand"#The upper-quartile modified FPKM calculation "in which the total protein-coding read count is replaced by the 75th percentile read count value for the sample." #https://www.biostars.org/p/213439/


## See what the non-normalized counts look like.
#    These are what we will use in DESeq2
View(as.data.frame(assays(dds)[["unstranded"]]))


## See what the tpm-transformed values look like
#    We won't use these because DESeq2 does its own tranformations
View(as.data.frame(assays(dds)[["tpm_unstrand"]]))


###### Check for multiple observations / undesired tissues, etc #####
if(sum(duplicated(dds$sample_submitter_id))>0){
  warning("Some IDs present more than once. Consider subsetting if this was unintentional.")
}
if(length(unique(dds$sample_type))>1){
  warning("More than one type of sample tissue type present. Consider subsetting if this was unintentional.")
}


###### Identify the groupings of samples ######
## Create new "columns" in the sample data based on grouping
dds$group1 <- dds$sample_submitter_id%in%group1$ID.1
dds$group2 <- dds$sample_submitter_id%in%group2$ID.1


## These should be the opposite of each other if barcode based subsetting was successful
if(!all(dds$group1==!dds$group2)){
  stop("Your groupings are not mutually exclusive")
}


# The way you communicate categorical groupings to DESeq is via a
#    column name corresponding to a factor vector.
# Deseq2 uses the order of the factor's levels to determine order
#    of levels in the calculations. Then there are only two levels:
#    The first level is set as the "control".
#    The second level is set as the "manipulated".
# Our experiments likely won't have true "controls", so we have to
#    pay attention to this solely to determine the order of
#    subsequent calculations such as log2 fold change
# The order for log2 fold change will be: level 2 - level 1
dds$comp <- factor(dds$group2,levels = c("FALSE","TRUE"))

## Fixing a rare issue:
# DESeq doesn't like spaces in levels, so you can remove them at the "level" level.
# This isn't needed if not using levels with spaces (e.g. FALSE and TRUE)
# Very few students will encounter this issues, but it is otherwise hard
#    to debug based on the error message.
levels(dds$comp)<-gsub(" ","_",levels(dds$comp))



###### Use DESeq2 to filter the data more! ######
## Convert to DESeq object
dds <- DESeqDataSet(dds, design = ~ comp)

## Normalize the read depth values
dds <- estimateSizeFactors(dds)

## Filter loci with extremely low RD
## Isolate the raw counts
rawCounts  <- as.data.frame(counts(dds, normalized = F))

## Filter based on % of each grouping with 0 RD
grp1PctWith0 <- rowMeans((rawCounts[,dds$group1]==0))
grp2PctWith0 <- rowMeans((rawCounts[,dds$group2]==0))
maxPct0Cutoff <- 0.9

## Visualize low rd loci
hist(c(grp1PctWith0,grp2PctWith0),1000)
abline(v = maxPct0Cutoff,col="red")

##Do the subset
pctWith0Logic <- grp1PctWith0<maxPct0Cutoff&grp2PctWith0<maxPct0Cutoff
dds <- dds[which(pctWith0Logic),]


##### Filter samples based on unnormalized RD #####
## Omit samples with extremely different RD
# The barcodes we will use to remove samples are in the column names which colMeans returns.
hist(colMeans(rawCounts),1000,border="blue")
colMeans(rawCounts)[which.max(colMeans(rawCounts))]
badSamples <- c(
  "TCGA-AH-6903-01A-11R-1928-07"# From: which.max(colMeans(rawCounts))# based on histogram
  #Customize this / add multiple samples based on your data's distribution
)

##### Remove redundant count data ####
## This is just to save memory
rm("rawCounts")


##### Filter samples based on normalized RD ####
## Isolate the normalized read depths
normCounts <- as.data.frame(counts(dds, normalized = T))

## Determining the most variable rows
perLocusVar   <- apply(normCounts,1,var)

## Calculating the 25% quantile of this data
pcaLocusLogic <- perLocusVar>quantile(perLocusVar,probs=0.25)

# Create PNG graphics device 
png(filename = "pca_plot.png", width = 480, height = 480)

## Use a principal component analysis to visualize the distribution 
pca <- prcomp(t(normCounts[which(pcaLocusLogic),]),scale. = T)

## Plot PCA
plot(pca$x[,1],pca$x[,2],asp = 1)  
abline(v = pca1Cutoffs,col="red")

# Save plot and close device
dev.off()

#Use logic to isolate and name "bad" samples
pcaLogic <- pca$x[,1]>pca1Cutoffs[1]&pca$x[,1]<pca1Cutoffs[2]
names(which(!pcaLogic))
badSamples <- c(
  badSamples, # Add to the previous badSamples object
  "TCGA-CD-A48A-01A-12R-A36D-31"# From the above code based on the plot
)

##### Utilize the bad sample vector to subset dds ####
dds <- dds[,which(!dds$barcode%in%badSamples)]

##### Remove redundant count data ####
rm("normCounts")






#-------------------------------------------#
#### Do differential expression analysis ####
#-------------------------------------------#
###### Convert to DESeq object ######
## This is redundant with a previous step,
## However, adding helps to make sure the object is intact after filtering.
dds <- DESeqDataSet(dds, design = ~ comp)

## Because we have significantly changed the loci used, recalculate normalization
dds <- estimateSizeFactors(dds)

## Run the Deseq2 package analysis
dds <- DESeq(dds)
res <- results(dds) # Organizes the results


##### Accumulate data into a data.frame #####
## Add in gene data from the rowRanges section of the SummarizedExperiment object
## Adding columns with the same name causes problems with ggplot2, so
## It only adds other columns
colsAdded <- !colnames(as.data.frame(rowRanges(dds)))%in%colnames(as.data.frame(res))
resOutput <- cbind(
  as.data.frame(rowRanges(dds))[,colsAdded],
  as.data.frame(res)
)



#--------------------------#
#### Visualize analysis ####
#--------------------------#
##### Note ####
#None of the following will work!!!!!!!!!!!!
#  until you have created a script file that codes
#  volcPlot() and geneViolin()

##### Make a volcano plot ####
# The following loads a script file that you will make as part of a HW question
# This will not work for your system
ro<-resOutput
volcPlot<-function(ro){
  library(ggplot2)
  nonMissingCount<-sum(
    !is.na(ro$log2FoldChange)&
      !is.na(ro$padj)&
      !is.na(ro$baseMean)
  )
  ggplot(ro,aes(log2FoldChange,-log10(padj),color=log10(baseMean)))+
    geom_point()+
    scale_color_viridis_c()+
    labs(title=paste0("Volcano plot of ",nonMissingCount," total gene features"))
}
volcPlot(resOutput)

volcPlot<-function(ro){
  library(ggplot2)
  nonMissingCount<-sum(
    !is.na(ro$log2FoldChange)&
      !is.na(ro$padj)&
      !is.na(ro$baseMean)
  )
  p1 <- ggplot(data = ro,
               mapping = aes(
                 x = log2FoldChange,
                 y = -log10(padj),
                 color=log10(baseMean)
               ))+
    geom_point()+
    scale_color_viridis_c()+
    labs(title=paste0("Volcano plot of ",nonMissingCount," total gene features"))
  return(p1)
}


volc1 <- volcPlot(
  dds_fun = dds,
  resPlotted = resOutput,
  afc_min = opt$afc.min,
  fdr_max = opt$fdr.max,
  n_per_side = opt$nPerSide,
  project = opt$project,
  groupAText = opt$group1Desc,
  groupBText = opt$group2Desc
)
volc1

volcPlot <- function(dds_fun, resPlotted, afc_min, fdr_max, n_per_side,  
                     project, groupAText, groupBText) {
  
  resPlotted$sig <- resPlotted$padj < fdr_max & 
    abs(resPlotted$log2FoldChange) > afc_min
  
  top_up   <- head(rownames(resPlotted[order(resPlotted$log2FoldChange, decreasing = T),]), n_per_side) 
  top_down <- head(rownames(resPlotted[order(resPlotted$log2FoldChange, decreasing = F),]), n_per_side)
  
  ggplot(resPlotted, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
    geom_point() +
    scale_color_manual(values = c("grey50","firebrick4")) + 
    geom_text_repel(data = subset(resPlotted, sig), 
                    aes(label = gene_name), size = 2) +   
    geom_vline(xintercept = c(-afc_min, afc_min), linetype = "dashed") +
    geom_hline(yintercept = -log10(fdr_max), linetype = "dashed") +
    labs(title = paste0(project, " MUC16 Gene Differential Expression"),
         x = "Log2 Fold Change",
         y = "-Log10 Adjusted P-value") +
    theme(legend.position = "none")
}


### Violin

# Filter the results for the most significant DE gene
significant_gene <- res[which.min(res$padj), ]

# Extract raw counts for the selected gene
gene_id <- rownames(significant_gene)
expression_data <- counts(dds)[gene_id, ]

# Log2-transform the expression values
expression_data <- log2(expression_data + 1)  # Adding 1 to avoid log(0)

# Create a data frame for plotting
plot_data <- data.frame(
  Group = factor(dds$comp),
  Expression = expression_data
)

# Load the required library
library(ggplot2)

# Create the violin plot
viol1 <- ggplot(plot_data, aes(x = Group, y = Expression)) +
  geom_violin(trim = FALSE, scale = "width", fill = "lightblue") +
  geom_jitter(position = position_jitter(.2), color = "black", size = 2) +
  labs(title = paste("Violin Plot for Gene:", gene_id),
       x = "Group",
       y = "Count") +
  theme_minimal()


#---------------------------------------------#
#### Save R environment for later analysis ####
#---------------------------------------------#
##### Save the plot ####
ggsave(filename = "muc16Plot_volc.png",volc1,
       width = 9,height = 9,dpi = 600)

ggsave(filename = "apcPlot_viol1.png",viol1,
       width = 9,height = 9,dpi = 600)


#You should be able to download my version of the image from the website
## I am first going to remove certain redundant data to decrease file size

save.image(paste0(
  "afterDeseq_",
  Sys.Date(),
  "_",
  as.numeric(Sys.time()),
  ".rimage")
)