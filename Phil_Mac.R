library(CINSignatureQuantification)

#Load example Data
data("TCGA_478_Samples_SNP6_GOLD")
head(TCGA_478_Samples_SNP6_GOLD)

## The first step is to load the data into the 'CNQuant' object. This is a S4 class object which standardises the input and additionl meta data fro downstream analysis
# Here the 'data' argument can be a loaded R data.frame object (as shown in this vignette), a file path to a delimited file, or a 'QDNAseqCopyNumbers' object from QDNASeq which has been segmented. A name for the experiment and the genome build to use are also specified. Viewing this object will print details of the contained data

cnobj <- createCNQuant(data = TCGA_478_Samples_SNP6_GOLD, experimentName = "VignetteExample", build = "hg19")
cnobj

### CNQuant object
# Additional meta data can be retrieved with getExperiment()

getExperiment(cnobj)

#The 'CNQuant' object conatins slots for all the required prerequisite data for copy number signature analysis. Of these,
# getSegments() returns the original segment table
segmentTable <- getSegments(cnobj)
head(segmentTable)

#getSamples() returns all sample IDs
getSamples(cnobj)[1:10]

# getSamplefeatures() returns sample features computed on CNQuant initialisation
Samplefeatures <- getSamplefeatures(cnobj)
head(Samplefeatures)

#If you have additional sample features or clinical data you can add those to the object using the `addsampleFeatures()` function
data(test.sample.features)
head(test.sample.features)
#id.col specifies the sample ID in the new sample features object
cnobj <- addsampleFeatures(object = cnobj, sample.data = test.sample.features, id.col = "sample")
NewSampleFeatures <- getSamplefeatures(cnobj)
head(NewSampleFeatures)

#Sub setting og `CNQuant` objects is implemented using the native R bracket `[` notation. This can be performed using either numerical indexing or specifying available sample identifiers. This sub setting functions at any point in the analysis pipeline

cnobj[1:10]

cnobj[c(1,24,5,77,100)]

cnobj[getSamples(cnobj)[1:10]]

#For the sake of computational efficiency, we will subset the full data set to 20 samples.

cnobj <- cnobj[1:20]
cnobj

############CALCULATE FEATURES##################
#Feature distributions can be calculated using the `calculateFeatures()` function.

cnobj <- calculateFeatures(object = cnobj,method = "mac", cores = 1)
cnobj

#These features can be retrieved as follows
#getfeats() returns list of data.frames containing the computed copy number features

feats <- getFeatures(cnobj)
head(feats[[1]])

########### Calculate samply-by-component

cnobj <- calculateSampleByComponentMatrix((object = cnobj))
cnobj

#The sample by component matrix can be retrieved as follows
# getSxC() returns the computed sample-by-component matrix

SxC <- getSampleByComponent(cnobj)
head(SxC)

########## Calculate signature activities

cnobj <- calculateActivity(object = cnobj)
cnobj

#The sample-by-signature activity matrix can be retrived as follows
# getActivity() returns the computed sample-by-component matrix

sigAct <- getActivities(cnobj, type = "threshold")
head(sigAct)

######### One-click function
#A wrapper function is provided to allow the entire pipeline to run as a single command. Here, the output object is a `SigQuant` object with all slots and calculations performed.

# Subset samples from the total dataset
subsample <- as.character(sample(TCGA_478_Samples_SNP6_GOLD$sample,size = 20))
TCGA_478_Samples_SNP6_GOLD_subset <- TCGA_478_Samples_SNP6_GOLD[TCGA_478_Samples_SNP6_GOLD$sample %in% subsample,]

cnobj <- quantifyCNSignatures(object = TCGA_478_Samples_SNP6_GOLD_subset,
                              experimentName = "VignetteExample",
                              method = "mac",
                              cores = 1,
                              build = "hg19")
cnobj

#Plotting Functions
#Sample can be specified as a numerical index or sample name
plotSegments(object = cnobj, sample = 1, cn.max = 8)

#Additional arguments can be provided for the heatmap() function
plotSampleByComponent(object = cnobj)

#Custom colours for signatures can be provided with the cols argument, the length of which must match the number of signatures
plotActivities(object = cnobj,type = "threshold")
