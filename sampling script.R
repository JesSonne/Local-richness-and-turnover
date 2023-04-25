# visualyzing sampling algorithm
require(terra) # for handling the geographical data and computations

setwd("/Users/jespersonne/Documents/GitHub/Local-richness-and-turnover")
reference_raster=rast("Data files/reference_raster.tif")

#global mountain regions downloaded from Rahbek et al. 2019: https://macroecology.ku.dk/resources/mountain_regions
setwd("/Users/jespersonne/Downloads")
mount=vect("mountainregions2/CMEC_Mountains_Enh2018.shp")
mount=project(mount,reference_raster)

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
plot(mount[po1,],col="lightgray",border="lightgray",axes=F)
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



