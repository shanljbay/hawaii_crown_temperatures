#copied from early exploration script called "lidR_attempt.R"
#trying to add canopy traits to first paper 

library(raster)
#install.packages("lidR")
library(lidR)
#install.packages("rgdal", configure.args = c("--with-proj-lib=/usr/local/lib/", "--with-proj-include=/usr/local/include/"))
library(rgdal)
library(sf)
library(sp)
library(rgeos)
library(maptools)
library(ggplot2)
library(exactextractr)

#library(rlas)

las = readLAS("~/Hawaii_Data_From_The_Sky/March2022_Laupahoehoe/Lidar/24_Final_Products/GatorEye_SubjectTerms_20220307-201126_final_CHM_merged_dem.las", filter="-drop_z_below 0")
#select = "xyzic"
chm <- raster("~/Hawaii_Data_From_The_Sky/March2022_Laupahoehoe/Lidar/24_Final_Products/GatorEye_SubjectTerms_20220307-201126_final_CHM_merged_dem.tif")
#WGS84 UTM5N
st_crs(las) <- 32605
crs(chm) <- 32605

print(las)
summary(las)
str(las)
plot(las)

crowns <- readOGR("Crowns.shp")

#upload crowns first below
lascrown <- clip_roi(las, crowns)

summary(lascrown)
str(lascrown)
#make summary object first then call from that... 
summarylascrown <- summary(lascrown)
str(summarylascrown)
str(lascrown)

plot(lascrown[[2]])
lc1 <- lascrown[1]
summary(lc1)

col <- pastel.colors(200)

crowns.df <- crowns@data
View(crowns.df)
crowns.df <- crowns.df[,-c(14:15)] #remove empty columns at end
crowns.df <- cbind(crowns.df, c(1:770))
colnames(crowns.df)[14] <- "crown.reference"
crowns.df <- cbind(crowns.df, area_density[,c(6,7)])

chm.c <- crop(chm, crowns)
plot(chm.c, col=height.colors(50), asp=NA)
plot(crowns, add=TRUE)

chm.c@data

chm.stats <- exact_extract(chm.c, crowns, c("mean", "stdev", "coefficient_of_variation", "variance", "min", "max", "median", "mode", "count"))

# add an ID column for 1:770 crowns that matches those from Lau-4-byIndividualCrowns...
chm.stats <- cbind(chm.stats, c(1:770))
colnames(chm.stats)[10] <- "crown.reference"
crowns$crown.reference <- c(1:770)

#rumple index:
?rumple_index
?polygon_metrics

#rumple index is 1 for a perfectly flat surface
rumple_index(chm.c) #6.338

#but want separate by species... by individual crown
library(tidyverse)
library(raster)
library(sf)
library(sp)

metpol.crowns <- subset(crowns, crowns$species=="METPOL")
metpol.chm <- mask(chm.c, metpol.crowns)
rumple_index(metpol.chm) #4.072239

acakoa.crowns <- subset(crowns, crowns$species=="ACAKOA")
acakoa.chm <- mask(chm.c, acakoa.crowns)
rumple_index(acakoa.chm) #5.284669

chetri.crowns <- subset(crowns, crowns$species=="CHETRI")
chetri.chm <- mask(chm.c, chetri.crowns)
rumple_index(chetri.chm) #4.546374

coprhy.crowns <- subset(crowns, crowns$species=="COPRHY")
coprhy.chm <- mask(chm.c, coprhy.crowns)
rumple_index(coprhy.chm) #4.358118


#### make each crown a separate layer in rasterstack
rp <- do.call(stack,lapply(1:nrow(crowns), function(x) 
  raster::mask(crop(chm.c, extent(crowns)), crowns[x,])))

plot(rp[[1]])
plot(rp[[2]])
#YAY!
plot(rp)

rumple_index(rp[[1]])
rumple_index(rp[[2]])

# Define a function to calculate the rumple index for a single layer
rumple_index_single <- function(layer) {
  return(rumple_index(layer))
}

# Create an empty list to store results
rumple_results <- list()

# Loop through each layer in the RasterStack
for (i in 1:nlayers(rp)) {
  layer <- rp[[i]]
  rumple_result <- rumple_index_single(layer)
  rumple_results[[i]] <- rumple_result
}

#add to chm.stats
chm.stats$rumple_results <- unlist(rumple_results)
head(chm.stats)

#Calculate Height Variability:
#A proxy for leaf clumping can be the variation in canopy height within each polygon (crown). A higher variation in heights would imply more structural complexity (possibly leaf clumping), whereas a lower variation implies a more even, uniform canopy.

# Extract height values within each polygon and compute standard deviation
leaf_clumping <- exact_extract(chm, crowns, fun = function(values, cov_frac) {
  stddev_height <- sd(values, na.rm = TRUE)
  return(stddev_height)
})

chm.stats$leaf_clumping <- unlist(leaf_clumping)
head(chm.stats)

#Alternative: Coefficient of Variation (CV): If you want a normalized measure of variability (i.e., to compare between crowns of different sizes or mean heights), you can use the coefficient of variation (CV), which is the ratio of the standard deviation to the mean height.
leaf_clumping_cv <- exact_extract(chm, crowns, fun = function(values, cov_frac) {
  mean_height <- mean(values, na.rm = TRUE)
  stddev_height <- sd(values, na.rm = TRUE)
  cv_height <- stddev_height / mean_height
  return(cv_height)
})
chm.stats$leaf_clumping_cv <- unlist(leaf_clumping_cv)
head(chm.stats)
#Higher clumping index (standard deviation or CV) implies a more uneven or "clumped" distribution of canopy heights (greater structural variability).
#Lower clumping index implies a more uniform canopy height, suggesting a more evenly distributed leaf structure.
par(mfrow=c(1,2))
hist(leaf_clumping)
hist(leaf_clumping_cv)

#voxel-based
voxel_data <- voxel_metrics(las, length(Z))  # did not specify "res" 
head(voxel_data)

# Convert voxel_data to a spatial data frame (sf object) for intersection
voxel_sf <- st_as_sf(voxel_data, coords = c("X", "Y", "Z"), crs = crs(las))

crowns_sf <- st_as_sf(crowns)

# Perform spatial intersection to get voxels within crowns
voxel_within_crowns <- st_intersection(voxel_sf, crowns_sf)

library(dplyr)
# Check the result
head(voxel_within_crowns)
# Summarize voxel data by crown (e.g., mean density per crown)
crown_metrics_voxel <- voxel_within_crowns %>%
  group_by(crown.reference) %>%  # each crown
  summarize(mean_density = mean(V1), na.rm = TRUE)  # V1 contains point counts (density)

# Check the result
print(crown_metrics_voxel)
head(crown_metrics_voxel)
crown_metrics_voxel$crown.reference #CHECKIGN WHY ONLY 769 instead of 770 -- crown.reference #541 is missing
hist(crown_metrics_voxel$mean_density)


chm.stats <- merge(crowns.df, chm.stats, "crown.reference")

chm.stats <- merge(chm.stats, crown_metrics_voxel, "crown.reference")

chm.stats <- chm.stats[,-c(30:31)]
head(chm.stats)

names(chm.stats)[names(chm.stats) == "mean_density"] <- "mean_density_voxels"
head(chm.stats)

chm.stats <- data.frame(chm.stats)

### order of plots and filter bad crowns
str(chm.stats)
chm.stats$species <- as.factor(chm.stats$species)
chm.stats$confidence <- as.factor(chm.stats$confidence)
levels(chm.stats$species)
chm.stats$species <- factor(chm.stats$species, levels=c('METPOL','ACAKOA','CHETRI','COPRHY','TROPASH', 'dead'))
levels(chm.stats$species)

# Remove NA species and "dead" crowns in a single step
chm.stats <- chm.stats[!is.na(chm.stats$species) & chm.stats$species != 'dead', ]

# Create a filtered dataset while retaining the original chm.stats
chm.stats <- subset(chm.stats, 
                      !confidence %in% c('low', 'lowmed', 'med-diffcolor', 'dead') & 
                        !species %in% c('TROPASH', 'dead'))


metpol.color <- "#009E73"#"#8ede61"
acakoa.color <- "#CC79A7"#"#b061de" 
chetri.color <- "#E69F00"#"#de7261"
coprhy.color <- "#56B4E9"#"#61cdde"

library(ggplot2)
library(ggpubr)

rumple.plot <- ggplot(chm.stats, aes(species, rumple_results, fill=species))+
  geom_point(aes(colour = species), size=0.6,
             position = position_jitterdodge())+
  geom_boxplot(notch=T, alpha=0.4, outlier.shape=NA)+
  coord_cartesian(ylim=c(1, 10))+
  ylab("Rumple Index")+
  xlab("")+
  scale_y_continuous(breaks=c(0, 2, 4, 6, 8, 10))+
  scale_fill_manual(values=c(metpol.color, acakoa.color, chetri.color, coprhy.color))+
  scale_color_manual(values=c(metpol.color, acakoa.color, chetri.color, coprhy.color))+
  theme_classic()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.position = "none")
rumple.plot

height.plot <- ggplot(chm.stats, aes(species, Z_Mean, fill=species))+
  geom_point(aes(colour = species), size=0.6,
             position = position_jitterdodge())+
  geom_boxplot(notch=T, alpha=0.4, outlier.shape=NA)+
  ylab("Mean Crown Height (m)")+
  xlab("")+
  scale_y_continuous(breaks=c(0, 5, 10, 15, 20, 25))+
  scale_fill_manual(values=c(metpol.color, acakoa.color, chetri.color, coprhy.color))+
  scale_color_manual(values=c(metpol.color, acakoa.color, chetri.color, coprhy.color))+
  theme_classic()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.position = "none")
height.plot


height.plot.max <- ggplot(chm.stats, aes(species, Z_Max, fill=species))+
  geom_point(aes(colour = species), size=0.6,
             position = position_jitterdodge())+
  geom_boxplot(notch=T, alpha=0.4, outlier.shape=NA)+
  ylab("Maximum Crown Height (m)")+
  xlab("")+
  scale_y_continuous(breaks=c(0, 5, 10, 15, 20, 25, 30))+
  scale_fill_manual(values=c(metpol.color, acakoa.color, chetri.color, coprhy.color))+
  scale_color_manual(values=c(metpol.color, acakoa.color, chetri.color, coprhy.color))+
  theme_classic()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.position = "none")
height.plot.max

clumping.plot <- ggplot(chm.stats, aes(species, leaf_clumping_cv, fill=species))+
  geom_point(aes(colour = species), size=0.6,
             position = position_jitterdodge())+
  geom_boxplot(notch=T, alpha=0.4, outlier.shape=NA)+
  ylab("Leaf Clumping (CV)")+
  xlab("")+
  coord_cartesian(ylim=c(0, 0.4))+
  #scale_y_continuous(breaks=c(0, 5, 10, 15, 20, 25))+
  scale_fill_manual(values=c(metpol.color, acakoa.color, chetri.color, coprhy.color))+
  scale_color_manual(values=c(metpol.color, acakoa.color, chetri.color, coprhy.color))+
  theme_classic()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.position = "none")
clumping.plot

density.plot.voxels <- ggplot(chm.stats, aes(species, mean_density_voxels, fill=species))+
  geom_point(aes(colour = species), size=0.6,
             position = position_jitterdodge())+
  geom_boxplot(notch=T, alpha=0.4, outlier.shape=NA)+
  ylab(expression("Mean Density (pts/m"^3*")"))+
  xlab("")+
  #scale_y_continuous(breaks=c(0, 5, 10, 15, 20, 25))+
  scale_fill_manual(values=c(metpol.color, acakoa.color, chetri.color, coprhy.color))+
  scale_color_manual(values=c(metpol.color, acakoa.color, chetri.color, coprhy.color))+
  theme_classic()+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        legend.position = "none")
density.plot.voxels




###trait correlation plots for each species

library(corrplot)
library(RColorBrewer)

# Calculate correlation matrix
str(cdf.f.metpol)
print(colnames(cdf.f.metpol))
corr_matrix.metpol <- cor(cdf.f.metpol[, c("Z_Mean.x", "Z_Max.x", "area_m2","density_ptsm2", "rumple_results", "leaf_clumping", "leaf_clumping_cv","mean_density_voxels")])
# Create correlation plot with customization
par(mfrow=c(1,1))
corrplot(corr_matrix.metpol, method = "color",
         col = brewer.pal(n = 8, name = "RdYlBu"),
         addCoef.col = "black",  # Add correlation coefficients
         tl.col = "black",       # Text label color
         tl.srt = 45)            # Text label rotation

corr_matrix.acakoa <- cor(cdf.f.acakoa[, c("Z_Mean.x", "Z_Max.x", "area_m2","density_ptsm2", "rumple_results", "leaf_clumping", "leaf_clumping_cv", "mean_density_voxels")])
# Create correlation plot with customization
par(mfrow=c(1,1))
corrplot(corr_matrix.acakoa, method = "color",
         col = brewer.pal(n = 8, name = "RdYlBu"),
         addCoef.col = "black",  # Add correlation coefficients
         tl.col = "black",       # Text label color
         tl.srt = 45)            # Text label rotation

corr_matrix.chetri <- cor(cdf.f.chetri[, c("Z_Mean.x", "Z_Max.x", "area_m2","density_ptsm2", "rumple_results", "leaf_clumping", "leaf_clumping_cv", "mean_density_voxels")])
# Create correlation plot with customization
par(mfrow=c(1,1))
corrplot(corr_matrix.chetri, method = "color",
         col = brewer.pal(n = 8, name = "RdYlBu"),
         addCoef.col = "black",  # Add correlation coefficients
         tl.col = "black",       # Text label color
         tl.srt = 45)            # Text label rotation

corr_matrix.coprhy <- cor(cdf.f.coprhy[, c("Z_Mean.x", "Z_Max.x", "area_m2", "density_ptsm2", "rumple_results",  "leaf_clumping", "leaf_clumping_cv",  "mean_density_voxels")])
# Create correlation plot with customization
par(mfrow=c(1,1))
corrplot(corr_matrix.coprhy, method = "color",
         col = brewer.pal(n = 8, name = "RdYlBu"),
         addCoef.col = "black",  # Add correlation coefficients
         tl.col = "black",       # Text label color
         tl.srt = 45)            # Text label rotation

corr_matrix.metpol4 <- cor(cdf.f.metpol[, c("Z_Mean.x", "rumple_results", "leaf_clumping_cv","mean_density_voxels")])
par(mfrow=c(1,1))
mp4.plot <- corrplot(corr_matrix.metpol4, method = "color",
                     col = brewer.pal(n = 8, name = "RdYlBu"),
                     addCoef.col = "black",
                     tl.col = "black",       
                     tl.srt = 45)            
mp4.plot

corr_matrix.acakoa4 <- cor(cdf.f.acakoa[, c("Z_Mean.x", "rumple_results", "leaf_clumping_cv","mean_density_voxels")])
ak4.plot <- corrplot(corr_matrix.acakoa4, method = "color",
                     col = brewer.pal(n = 8, name = "RdYlBu"),
                     addCoef.col = "black",
                     tl.col = "black",       
                     tl.srt = 45)            
ak4.plot

corr_matrix.chetri4 <- cor(cdf.f.chetri[, c("Z_Mean.x", "rumple_results", "leaf_clumping_cv","mean_density_voxels")])
ct4.plot <- corrplot(corr_matrix.chetri4, method = "color",
                     col = brewer.pal(n = 8, name = "RdYlBu"),
                     addCoef.col = "black",
                     tl.col = "black",       
                     tl.srt = 45)            
ct4.plot

corr_matrix.coprhy4 <- cor(cdf.f.coprhy[, c("Z_Mean.x", "rumple_results", "leaf_clumping_cv","mean_density_voxels")])
cr4.plot <- corrplot(corr_matrix.coprhy4, method = "color",
                     col = brewer.pal(n = 8, name = "RdYlBu"),
                     addCoef.col = "black",
                     tl.col = "black",       
                     tl.srt = 45)            
cr4.plot

mp4.plot
ak4.plot
ct4.plot
cr4.plot



