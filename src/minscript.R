library(sf)
library(tidyverse)
library(lwgeom)

load("~/Dropbox/Apps/Texpad/Assessment of current 3Ps Cod Stock Assessment/crap.RData")

ageGrowth <- campelen$raw.data$age.growth %>%
    mutate(survey.year = as.factor(survey.year))

setLatLon <- campelen$raw.data$set.details %>%
    dplyr::select(shiptrip,survey.year,season,day,month,lat.start,long.start,strat,set) %>%
    mutate(survey.year = as.factor(survey.year),
           long.start = -long.start)

##UTM Zone for Newfoundland, very very important
st_crs(32621)
##This exact proj4string is what we want
##Because units
utm_proj <- "+proj=utm +zone=21 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"


ageGrowth <- ageGrowth %>% inner_join(setLatLon) %>%
    st_as_sf(coords=c("long.start","lat.start"),crs=4326) %>%
    st_transform_proj(crs=utm_proj) %>%
    filter(sex !="unsexed") %>%
     mutate(utmE.start = st_coordinates(.)[,2],
            utmN.start = st_coordinates(.)[,1],
            sex = as.factor(sex)) 

cag2012 <- filter(ageGrowth,survey.year == 2012,!is.na(age))


library(barrierALK)
data(barrierMesh)


mesh = barrierMesh$mesh
tri = barrierMesh$barrier_tri

print("Getting Ready...")

yo = system.time(largerRSA <- barrierALK(age ~ length + gf(cag2012$utmN.start,cag2012$utmE.start,
                                                      mesh=mesh,model='barrier',barrier.triangles=tri,
                                                      rangeKey=c(1,1,2,2,2,3,3,3,3),
                                                      sdUkey=c(1,1,2,2,2,3,3,3,3)),equalSlopes=FALSE,
                                       minAge=1,maxAge=10,silent=FALSE,data=cag2012,runSA=TRUE,
                                    control=list(eval.max=5000,iter.max=5000)))

print(yo)
