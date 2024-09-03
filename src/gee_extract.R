# Earth Engine data extraction
# Author: Dario Schulz
# Date: 2021-09-29
# Purpose: This script extracts MODIS NDVI, population density and travel time 
# data from Earth Engine for the camera locations in the study area.

# Note: This script requires the rgee package to be installed and authenticated.
# for guidance on how to do this, see https://github.com/r-spatial/rgee
# I manually reprojected the camera GPS locations to EPSG:4326 in QGIS. 
# The cleaned data is stored in the file "repr_cam_locs.geojson". 


#### Setup ####

setwd("C:/Users/DaSchulz/OneDrive - European Forest Institute/Dokumente/research/ag2bio_sb")

# libraries
library(sf)
library(rgee) # earth engine
library(tidyverse)


# load reprojected camera locations

cm_locs <- st_read("./data/repr_cam_locs.geojson") %>%
  dplyr::select(all_of("Unique_Station_Name"))


ee_Initialize()


# prepare spatial objects for earth engine
hulls_ee <- sf_as_ee(cm_locs)
aoi_ee <- sf_as_ee(st_as_sfc(st_bbox(cm_locs)))

##### Extract data from Earth Engine #####

#### MODIS ####

# also get MODIS average and maximum annual NDVI
modis <- ee$ImageCollection("MODIS/061/MOD13A2")$
  select("EVI")$
  map(function(image) {
    image$multiply(0.0001)$copyProperties(image, image$propertyNames())
  })


# make annual means
years <- ee$List$sequence(2000, 2023)

# Function to calculate annual mean NDVI
calculate_annual_mean <- function(year) {
  year <- ee$Number(year)
  start_date <- ee$Date$fromYMD(year, 1, 1)
  end_date <- start_date$advance(1, "year")
  
  # Filter MODIS collection for the specific year and calculate the mean NDVI
  yearly_ndvi <- modis$
    filterDate(start_date, end_date)$
    mean()$
    set("year", year)
  
  return(yearly_ndvi)
}

# Apply the function over the list of years
annual_mean_ndvi <- ee$ImageCollection$fromImages(
  years$map(ee_utils_pyfunc(calculate_annual_mean))
)

# combine bands into one image
annual_mean_ndvi <- annual_mean_ndvi$toBands()

modis_ndvi_mean <- ee_extract(
  x = annual_mean_ndvi,
  y = hulls_ee,
  scale = 250,
  fun = ee$Reducer$mean()
)

# to long format
modis_ndvi_mean_long <- modis_ndvi_mean %>%
  pivot_longer(cols = ends_with("EVI"),
               names_to = "year",
               values_to = "ndvi") %>%
  mutate(year = as.numeric(gsub("X|_EVI", "", year)) + 2000)


# facet plot by cluster
modis_ndvi_mean_long %>%
  group_by(certification, site, year, cluster) %>%
  summarize(mean_ndvi = mean(ndvi),
            min_ndvi = min(ndvi),
            max_ndvi = max(ndvi)) %>%
  ggplot(aes(x = year, y = mean_ndvi, group = site)) +
  geom_line(aes(linetype = factor(certification)), lwd = 1.5) +
  geom_ribbon(aes(ymin = min_ndvi, ymax = max_ndvi, fill = site), alpha = 0.2) +
  theme_minimal() +
  labs(title = "MODIS annual mean NDVI in FSC and non-FSC sites",
       x = "Year",
       y = "NDVI") +
  theme(legend.position = "bottom") +
  facet_wrap(~cluster, ncol = 2)



#### Population density ####
pop_ee <- ee$ImageCollection("CIESIN/GPWv411/GPW_Population_Density")$
  select("population_density")

pop_mean <- ee_extract(
  x = pop_ee,
  y = hulls_ee,
  scale = 250,
  fun = ee$Reducer$mean()
)

pop_mean <- pop_mean %>%
  rename(pop2000 = "gpw_v4_population_density_rev11_2000_30_sec_population_density",
         pop2005 = "gpw_v4_population_density_rev11_2005_30_sec_population_density",
         pop2010 = "gpw_v4_population_density_rev11_2010_30_sec_population_density",
         pop2015 = "gpw_v4_population_density_rev11_2015_30_sec_population_density",
         pop2020 = "gpw_v4_population_density_rev11_2020_30_sec_population_density")


#### travel time ####
travel_ee <- ee$Image("Oxford/MAP/accessibility_to_cities_2015_v1_0")$
  select("accessibility")

access_mean <- ee_extract(
  x = travel_ee,
  y = hulls_ee,
  scale = 250,
  fun = ee$Reducer$mean()
)


#### Combine and save data ####

geedat <- modis_ndvi_mean %>%
  left_join(., pop_mean, by = "Unique_Station_Name") %>%
  left_join(., access_mean, by = "Unique_Station_Name")

save(geedat, file = "geedat.RData")
