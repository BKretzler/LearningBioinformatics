# ---------------------------------------------------------
# First attempt with statgenGWAS package
# 10 Aug 2021
# Bailey M Kretzler
# Data + EU DROPS project
# ---------------------------------------------------------


library(Rcpp) #if this acts up reinstall
library(statgenGWAS)
library(dplyr)

#---------------------------------------------------------#
#                 Getting data into gData form                      
#---------------------------------------------------------#

#---------------reading/inspecting data-------------------#


# reads the data in 
data("dropsMarkers") #assuming this contains gene markers
data("dropsMap") # assuming this is a chrom map?
data("dropsPheno") # ?????

rownames(dropsMarkers)
colnames(dropsMarkers)
dropsMarkers[(1:10),(1:5)]
  #lists number of recessive (?) alleles for each SNP in the columns
  #row are genotypes
  #basically gives you the layout of what is/isn't there

rownames(dropsMap)
colnames(dropsMap) 
head(dropsMap)
  # contains name of SNPs, chromosome#, position, + alleles
  #rows are basically SNPs
  # then the detainls of where they are located
  #confused about the allele stuff


rownames(dropsPheno)
colnames(dropsPheno)
head(dropsPheno)
  #simply phenotypic data
  
#------------------preparing the data---------------------#


### dropsMarkers ###

#make row names genotype    
rownames(dropsMarkers) <- dropsMarkers[["Ind"]]

#drop ind col
dropsMarkers<- dropsMarkers[colnames(dropsMarkers)!= "Ind"]

dropsMarkers[(1:10),(1:5)] #check -- looks good

### dropsMap ###

#rownames <- snp names
rownames(dropsMap) <- dropsMap$SNP.names

#rename chromosome and position columns
dropsMap = dropsMap %>% 
  rename(chr = Chromosome, pos = Position)

head(dropsMap)

### dropsPheno ###

#rename varietyID col
dropsPheno <- dropsPheno %>%
  rename(genotype = Variety_ID)

head(dropsPheno)

#convert relevant data to a list
dropsPhenoList <- split(x = dropsPheno[c("genotype", "grain.yield", 
                                         "grain.number", "seed.size",
                                         "anthesis", "silking", 
                                         "plant.height", "tassel.height",
                                         "ear.height")], 
                        f = dropsPheno[["Experiment"]])


#------------------converting to gData-----------------------#

gDataDrops <- createGData(geno = dropsMarkers, 
                          map = dropsMap, 
                          pheno = dropsPhenoList)
    #geno is the markes data (sheet with SNPs), GT must be the row names
    #map is obvs chrom map (where snps go), row names must be snp names
    #pheno is phenotypic data and must be as a list
    # i think..

View(gDataDrops)

#---------------------------------------------------------#
#                      Remove duplicate SNPs                      
#---------------------------------------------------------#


gDataDropsCL <- codeMarkers(gDataDrops, impute = FALSE, verbose = TRUE)


#---------------------------------------------------------#
#                      running the GWAS                      
#---------------------------------------------------------#


GWASdrops <- runSingleTraitGwas(gData = gDataDropsCL,
                                trials = "Mur13W",
                                traits = c("grain.yield", "anthesis"))

plot(GWASdrops, plotType = "qq", trait = "grain.yield")

plot(GWASdrops, plotType = "manhattan", trait = "grain.yield")

plot(GWASdrops, plotType = "qtl", yThr = 4, normalize = TRUE)


