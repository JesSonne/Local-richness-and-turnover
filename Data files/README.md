## Description of data files
**data per mountain region.csv:** data frame listing the species diversity metrics and environmental factors for each mountain region following Rahbek et al. 2019.

Columns in the spreadsheet:

1. 'mount.Name' name each mountain region.
2. 'mount.Area_km2' Total area of each mountain region in square. kilometres.
3. 'n_GridsCells' Number of 0.25 degree grid cells covering each mountain region.
4. 'total_richness' Total number of species found in each mountain region.
5. 'local_richness' Average species richness in 0.25 degree grid cells.
6. 'Species turnover' Whittaker’s effective species turnover = local_rcihness/total_richness-1.
7. 'MAT' Mean annual temperature (degrees celcius*10).
8. 'MAP' Annual precipitation (mm).
9. 'Climate volume' Area of the minimum convex polygon surrounding all conditions of MAT and MAP found within a focal mountain region.
10. 'TopographicComplexity' standard deviation in elevations within a moutrain region.
11. 'NPP' Average Net Primary productivity in units of kg carbon per m2 per year at 500 metre pixel resolution.
12. 'average_ADF_vol' Assemblage dispersion field volume averaged per 0.25 degree grid cell (0.001 scale factor).

**hab_list.RData:** Frequency distribution of habitat types samples within each mountain region, stored as a list of numeric vectors. See Jung et al. (2020) for a reference key.     

**std_local richness_50gridcells.csv:** 1000 replicates of a biogeographical null model, standardising the pattern of local species richness by sampling species from a coherent field of 50 grid cells.

**std_species turnover_50gridcells.csv:** 1000 replicates of a biogeographical null model, standardising the pattern of species turnover by sampling species from a coherent field of 50 grid cells.

**std_Sorensen dissimilarity_50gridcells.csv:** 1000 replicates of a biogeographical null model, standardising the pattern of Sørensen's multiple-site dissimilarity by sampling species from a coherent field of 50 grid cells.

**std_regional species richness_50gridcells.csv:** 1000 replicates of a biogeographical null model, standardising the pattern of regional species richness by sampling species from a coherent field of 50 grid cells.

**reference_raster.tif:** Spatial raster containing reference id-values for each 0.25 degree grid cell in the World. Used for running the biogeographical null model

### References
Rahbek, C., M. K. Borregaard, B. Hermansen, D. Nogues-Bravo, and J. Fjeldså. 2019. Definition and description of the Montane Regions of the World. https://macroecology.ku.dk/resources/mountain_regions.




