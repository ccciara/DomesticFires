install.packages("tidyverse")
install.packages("units")
install.packages("dplyr")
install.packages("sf")
install.packages("tmap")
install.packages("sp")
install.packages("spdep")
install.packages("stringr")
library("stringr")
library("tidyverse")
library("units")
library("sf")
library("tmap")
library("dplyr")
library("sp")
library("spdep")
library("ggplot2")

setwd("SET YOUR WORKING DIRECTORY HERE")

#IMD score________________________________________________________________

# get IMD and other independent factors, aggregate to ward level
imdscore <- read.csv("imdscore.csv")

# remove extra data from IMD (mapped to 2011 boundaries)
columns_to_keep <- c("LSOA.code..2011.",
                     "Index.of.Multiple.Deprivation..IMD..Score",
                     "Living.Environment.Score",
                     "")
imdscore <- dplyr::select(imdscore, one_of(columns_to_keep))
colnames(imdscore)[colnames(imdscore) == "LSOA.code..2011."] <- "lsoa_code"
colnames(imdscore)[colnames(imdscore) == "Index.of.Multiple.Deprivation..IMD..Score"] <- "imd_score"
colnames(imdscore)[colnames(imdscore) == "Living.Environment.Score"] <- "living_env"

# download and clean LSOA shapefile, merge with LSOA level factors
lsoa <- st_read("LSOA_2011_London_gen_MHW.shp")
columns_to_keep <- c("LSOA11CD")
lsoa <- dplyr::select(lsoa, one_of(columns_to_keep))
colnames(lsoa)[colnames(lsoa) == "LSOA11CD"] <- "lsoa_code"
lsoa <- merge(lsoa, imdscore, by = "lsoa_code")

# map to check everything
tm_shape(lsoa) + tm_fill("living_env", style = "cont") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 1.5, legend.text.size = 1)

# download 2022 wards, set LSOA crs to ward version
wards <- st_read("england_wd_2022_bgc.shp")
lsoa <- st_transform(lsoa, st_crs(wards))

# clean wards
columns_to_keep <- c("wd22cd")
wards <- dplyr::select(wards, one_of(columns_to_keep))
colnames(wards)[colnames(wards) == "wd22cd"] <- "ward_code"

# Calculate intersections
intersections <- st_intersection(lsoa, wards)

# prevent scientific notation for future steps
options(scipen = 999)

# Calculate area of overlap
intersections$overlap_area <- st_area(intersections)

# Calculate area of LSOAs
lsoa$area <- st_area(lsoa)

# Join the area back to intersections
intersections <- intersections %>%
  mutate(LSOA_code = st_drop_geometry(.)$lsoa_code) %>%
  left_join(st_drop_geometry(lsoa)[, c("lsoa_code", "area")], by = "lsoa_code") %>%
  mutate(overlap_percentage = overlap_area / area)

#aggregate data, weighted by LSOA area contribution to ward
ward_data <- intersections %>%
  group_by(ward_code) %>%
  summarise(
    imd_score = sum(imd_score * overlap_area) / sum(overlap_area),
    living_env = sum(living_env * overlap_area) / sum(overlap_area),
    .groups = 'drop'
  )

# add proper ward boundaries back on
ward_factors <- st_drop_geometry(ward_data)
ward_factors <- merge(wards, ward_factors, by = "ward_code")

# add household # and overcrowding data___________________________________
#get number of households
households <- read.csv("households.csv")
colnames(households)[colnames(households) == "LSOA.code"] <- "lsoa_code"
colnames(households)[colnames(households) == "All.Households"] <- "sum_hh"
columns_to_keep <- c("lsoa_code", "sum_hh")
households <- dplyr::select(households, one_of(columns_to_keep))

overcrowd <- read.csv("overcrowding21.csv")
colnames(overcrowd)[colnames(overcrowd) == "Table.1c..Occupancy.rating..bedrooms...lower.layer.super.output.areas..England.and.Wales..Census.2021"] <- "lsoa_code"
colnames(overcrowd)[colnames(overcrowd) == "X.1"] <- "overcrowded"

# drop first two rows
overcrowd <- overcrowd[-c(1, 2), ]

# drop last three columns
overcrowd <- overcrowd[, -c((ncol(overcrowd)-2):ncol(overcrowd))]

# merge with households
overcrowd <- merge(overcrowd, households, by = "lsoa_code")

# make columns numeric
last_cols <- tail(names(overcrowd), 2)
overcrowd[last_cols] <- lapply(overcrowd[last_cols], as.numeric)

# drop NAs
overcrowd <- na.omit(overcrowd)

# get proportion of overcrowded
overcrowd$overcrowd <- overcrowd$overcrowded / overcrowd$sum_hh

# drop extra columns
columns_to_keep <- c("lsoa_code", "sum_hh", "overcrowd")
overcrowd <- dplyr::select(overcrowd, one_of(columns_to_keep))

# aggregate to ward level__________________
#get 2021 lsoa boundaries
lsoa2021 <- st_read("england_lsoa_2021_bgc.shp")
columns_to_keep <- c("lsoa21cd")
lsoa2021 <- dplyr::select(lsoa2021, one_of(columns_to_keep))
colnames(lsoa2021)[colnames(lsoa2021) == "lsoa21cd"] <- "lsoa_code"
lsoa2021 <- merge(lsoa2021, overcrowd, by = "lsoa_code")

# map to check everything
tm_shape(lsoa2021) + tm_fill("overcrowd", style = "cont") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 1.5, legend.text.size = 1)

# download 2022 wards, set LSOA crs to ward version
wards <- st_read("england_wd_2022_bgc.shp")
lsoa2021 <- st_transform(lsoa2021, st_crs(wards))

# clean wards
columns_to_keep <- c("wd22cd")
wards <- dplyr::select(wards, one_of(columns_to_keep))
colnames(wards)[colnames(wards) == "wd22cd"] <- "ward_code"

# Calculate intersections
intersections <- st_intersection(lsoa2021, wards)

# prevent scientific notation for future steps
options(scipen = 999)

# Calculate area of overlap
intersections$overlap_area <- st_area(intersections)

# Calculate area of LSOAs
lsoa2021$area <- st_area(lsoa2021)

# Join the area back to intersections
intersections <- intersections %>%
  mutate(lsoa_code = st_drop_geometry(.)$lsoa_code) %>%
  left_join(st_drop_geometry(lsoa2021)[, c("lsoa_code", "area")], by = "lsoa_code") %>%
  mutate(overlap_percentage = overlap_area / area)

#aggregate data, weighted by LSOA area contribution to ward
overcrowd_data <- intersections %>%
  group_by(ward_code) %>%
  summarise(
    sum_hh = sum(sum_hh * overlap_area) / sum(overlap_area),
    overcrowd = sum(overcrowd * overlap_area) / sum(overlap_area),
    .groups = 'drop'
  )

# add proper ward boundaries back on
overcrowd_data <- st_drop_geometry(overcrowd_data)

# merge with imd
ward_factors <- merge(ward_factors, overcrowd_data, by = "ward_code")

tm_shape(ward_factors) + tm_fill("living_env", style = "cont") +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("right", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 1.5, legend.text.size = 1)


# add fire data____________________________________________________________
#fire data
fires <- read.csv("fire_incidents.csv")

#remove non-2021 years, non dwellings and non fires
fires2021 <- fires[grep("2021", fires$CalYear), ]
fires2021 <- fires2021[grep("Fire", fires2021$IncidentGroup), ]
fires2021 <- fires2021[grep("Dwelling", fires2021$PropertyCategory), ]

#drop all but locational data
columns_to_keep <- c("IncGeo_WardCode")
fires2021 <- dplyr::select(fires2021, one_of(columns_to_keep))

#change ward code name
colnames(fires2021)[colnames(fires2021) == "IncGeo_WardCode"] <- "ward_code"

#add count to fires
fires2021$fire_count <- 1

#aggregate to ward level
fires2021_aggregated <- fires2021 %>%
  group_by(ward_code) %>%
  summarise(across(everything(), sum, na.rm = TRUE),
            .groups = 'drop')

# merge with other ward level data
ward_factors <- merge(ward_factors, fires2021_aggregated, by = "ward_code")

# omit nas
ward_factors <- na.omit(ward_factors)

# plot factors ________________________________________________
boroughs <- st_read("London_Borough_Excluding_MHW.shp")

tm_shape(ward_factors) + tm_fill("living_env", style = "cont", title="Living Environment") +
  tm_shape(boroughs) + tm_polygons(alpha = 0, border.alpha = 0.9, border.col = "black") +
  tm_text("NAME", size = 0.5) +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("left", "bottom")) +
  tm_layout(frame = FALSE, legend.title.size = 1.3, legend.text.size = 1)

# save files for analysis____________________________________________________

#separate and save london ward shapefile
london_wards <- ward_factors %>%
  select(ward_code)
st_write(london_wards, "gl_wards.shp", driver = "ESRI Shapefile", append = FALSE)

# drop geometry
ward_factors <- st_drop_geometry(ward_factors)

# save file
write.csv(ward_factors, "london_fire_data.csv", row.names=FALSE)


