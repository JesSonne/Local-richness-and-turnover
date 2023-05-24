# visualyzing sampling algorithm
require(terra) # for handling the geographical data and computations

setwd("...Downloads/Local-richness-and-turnover-main")

#loading sampling functions
source("R functions.R")

#loading reference raster
reference_raster=rast("Data files/reference_raster.tif")

#global mountain regions downloaded from Rahbek et al. 2019: https://macroecology.ku.dk/resources/mountain_regions
mount=vect("...mountainregions2/CMEC_Mountains_Enh2018.shp")
mount=terra::project(mount,reference_raster)

i=86 # example with the Andes
mount$Name[i]
mountain_raster=mask(reference_raster,mount[i,])
mountain_raster=trim(mountain_raster)

# running the spatial sampling algorithm
set.seed(100)
sampling_field=spred_dye(size=50)

#plotting the results 
zero=reference_raster
zero[]=NA
zero[sampling_field]=1
zero=trim(zero)
plot(mount[i,],col="lightgray",border="lightgray",axes=F)
plot(zero,col="steelblue",add=T,legend=F,axes=F)

################################################
i=41 # example with the Eastern Central Andes (here named "AndreanYungas")
mount$Name[i]
mountain_raster=mask(reference_raster,mount[i,])
mountain_raster=trim(mountain_raster)

# running the spatial sampling algorithm
set.seed(53)
sampling_field=spred_dye(size=50)

#plotting the results 
zero=reference_raster
zero[]=NA
zero[sampling_field]=1
zero=trim(zero)
plot(mount[i,],col="lightgray",border="lightgray",axes=F)
plot(zero,col="steelblue",add=T,legend=F,axes=F)



