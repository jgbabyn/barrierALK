library(sf)
library(tidyverse)
library(lwgeom)


#compile("spatialBarrier.cpp","-O0 -g")
#dyn.load(dynlib('spatialBarrier'))


load("~/Dropbox/Apps/Texpad/Assessment of current 3Ps Cod Stock Assessment/crap.RData")

ageGrowth <- campelen$raw.data$age.growth %>%
    mutate(survey.year = as.factor(survey.year))

setLatLon <- campelen$raw.data$set.details %>%
    dplyr::select(shiptrip,survey.year,season,day,month,lat.start,long.start,strat,set) %>%
    mutate(survey.year = as.factor(survey.year),
           long.start = -long.start)

longLF <- campelen$raw.data$set.details %>%
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

setLatLon <- setLatLon %>%
    st_as_sf(coords=c("long.start","lat.start"),crs=4326) %>%
    st_transform_proj(crs=utm_proj) %>%
     mutate(utmE.start = st_coordinates(.)[,2],
            utmN.start = st_coordinates(.)[,1])

longLF <- longLF %>%
    st_as_sf(coords=c("long.start","lat.start"),crs=4326) %>%
    st_transform_proj(crs=utm_proj) %>%
    mutate(utmE.start = st_coordinates(.)[,2],
           utmN.start = st_coordinates(.)[,1]) %>%
    select(-c(af1:afNA))

longLF <- gather(longLF,"length","length_freq",lf4:lf130) %>%
    filter(length_freq != 0) %>%
    mutate(length = as.numeric(gsub("lf","",length)))

setLatLon2012 <- setLatLon %>%
    filter(survey.year == "2012")


cag2012 <- filter(ageGrowth,survey.year == 2012,!is.na(age))
cag1997 <- filter(ageGrowth,survey.year == 1997,!is.na(age))
cag <- filter(ageGrowth,!is.na(age))

library(barrierALK)
data(barrierMesh)

tmbD <- barrierALK:::interpret.salk(age ~ length:survey.year,FALSE,data=as.data.frame(cag),1,12)

library(TMB)

obj <- MakeADFun(tmbD$data,tmbD$parm)


library(optimParallel)

library(optimx)

cl <- makeCluster(spec=detectCores(),type="FORK",outfile="")
setDefaultCluster(cl=cl)
setDefaultCluster(cl=NULL)
stopCluster(cl)

library(rbenchmark)

benchmark("scale1"={
    opt1 <- nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=5000,iter.max=5000))},
    "scale10"={
        opt2 <- nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=5000,iter.max=5000),scale=2)
        },replications=5)
    

system.time(opt1 <- nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=5000,iter.max=5000)))
system.time(opt2 <- nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=5000,iter.max=5000),scale=10))



system.time(optPar <- optimParallel:::optimParallel(obj$par,obj$fn,obj$gr,control=list(maxit=5000,maxfeval=5000)))


opt <- nlminb(obj$par,obj$fn,obj$gr,control=list(eval.max=5000,iter.max=5000))

opt <- optimr(obj$par,obj$fn,obj$gr,method="nlminb",control=list(maxit=5000,maxfeval=5000))

system.time(opt2 <- optimr(obj$par,obj$fn,obj$gr,method="nlminb",control=list(maxit=5000,maxfeval=5000)))



mesh = barrierMesh$mesh
tri = barrierMesh$barrier_tri


lfCount <- dplyr::select(campelen$raw.data$set.details,lf4:lf130)
lfCount97 <- dplyr::filter(campelen$raw.data$set.details,survey.year=="1997") %>%
    select(lf4:lf130)
lb <- as.numeric(sub("lf","",names(lfCount)))





system.time(basic <- barrierALK(age ~ I(length/10):survey.year,1,10,data=as.data.frame(cag),equalSlopes=FALSE,control=list(eval.max=5000,iter.max=5000)))

library(rbenchmark)

benchmark("divide"={
    barrierALK(age ~ I(length/10):survey.year,1,10,data=as.data.frame(cag),equalSlopes=FALSE,control=list(eval.max=5000,iter.max=5000))},
    "log"={
        barrierALK(age ~ log(length):survey.year,1,10,data=as.data.frame(cag),equalSlopes = FALSE,control=list(eval.max=5000,iter.max=5000))},
    "regular"={
        barrierALK(age ~ length:survey.year,1,10,data=as.data.frame(cag),equalSlopes=FALSE,control=list(eval.max=5000,iter.max=5000))},replications=10)

system.time(predict(basic))

system.time(goop <-generateNatA(basic,lb,setLatLon,as.matrix(lfCount)))

system.time(basicNA <- getNatAge(basic,lb,setLatLon,lfCount))

system.time(basicALKs <- createALKs(basic,lb,setLatLon))


plot(basicALKs[[1]])
rage <- cag2012$age
rage[rage > 12] = 12

rage97 <- cag1997$age
rage97[rage97 > 10] = 10

rageF <- cag$age
rageF[rageF > 12 ] = 12

basY <- predict(basic,type="response")

bAcc <- sum(basY == rage97)/length(basY)


cag2 = cag

cag2$utmN.start = cag2$utmN.start/10
cag2$utmE.start = cag2$utmE.start/100

mesh2 = mesh
mesh2$loc[,1] = mesh2$loc[,1]/10
mesh2$loc[,2] = mesh2$loc[,2]/100


system.time(largerRSA <- barrierALK(age ~ length + gf(cag$utmN.start,cag$utmE.start,
                                                      mesh=mesh,model='barrier',barrier.triangles=tri,
                                                      rangeKey=c(1,1,2,2,2,3,3,3,3),
                                                      sdUkey=c(1,1,2,2,2,3,3,3,3)),equalSlopes=FALSE,
                                       minAge=1,maxAge=10,silent=FALSE,data=as.data.frame(cag),runSA=TRUE,
                                    control=list(eval.max=5000,iter.max=5000)))



system.time(largerRSA1 <- barrierALK(age ~ I(length*10) + gf(cag$utmN.start,cag$utmE.start,
                                                      mesh=mesh,model='barrier',barrier.triangles=tri,
                                                      rangeKey=c(1,1,2,2,2,3,3,3,3),
                                                      sdUkey=c(1,1,2,2,2,3,3,3,3)),equalSlopes=FALSE,
                                       minAge=1,maxAge=10,silent=FALSE,data=as.data.frame(cag),runSA=TRUE,
                                    control=list(eval.max=5000,iter.max=5000)))


system.time(largerRSA2 <- barrierALK(age ~ I(length/10) + gf(cag2$utmN.start,cag2$utmE.start,
                                                      mesh=mesh2,model='barrier',barrier.triangles=tri,
                                                      rangeKey=c(1,1,2,2,2,3,3,3,3),
                                                      sdUkey=c(1,1,2,2,2,3,3,3,3)),equalSlopes=FALSE,
                                       minAge=1,maxAge=10,silent=FALSE,data=as.data.frame(cag2),runSA=TRUE,
                                    control=list(eval.max=5000,iter.max=5000)))


probs = predict(largerRSA,newdata=as.data.frame(longLF),type='class')

testo = ageFish(largerRSA,longLF)

setLatLonY = split(setLatLon,setLatLon$survey.year)

system.time(lagerALKs <- lapply(setLatLonY,function(x) createALKs(largerRSA,lb,x)))

system.time(largePred <- createALKs(largerRSA,lb,setLatLonY[[1]]))

testo <- rbind(setLatLonY[[1]],setLatLonY[[2]],setLatLonY[[3]],setLatLonY[[4]])

load("needed.Rdata")

frT <- list()
fT <- list()

for(i in 1:100){
fT[[i]] <- system.time(fart <- barrierALK:::quickPredict(needed$X,needed$beta,needed$ranX,needed$A,needed$cohort,needed$ranXind))[1]
}
fT <- unlist(fT,recursive = FALSE)

for(i in 1:100){
frT[[i]] <- system.time(fartR <-barrierALK:::Rpred(needed$X,needed$beta,needed$ranX,needed$A,needed$cohort,needed$ranXind))[1]
}
frT <- unlist(frT,recursive=FALSE)

system.time(largePred <- createALKs(largerRSA,lb,testo))


system.time(largerALKs <- createALKs(largerRSA,lb,setLatLon))


larY <- predict(largerRSA,type="class")

system.time(larNA <- getNatAge(largerRSA,lb,setLatLon,lfCount))

longLF <- cbind(setLatLon,lfCount)

##how long will it take to predict this stuff?

lAcc <- sum(larY == rage97)/length(larY)

##Run for meeting
library(sf)
commArea <- st_read("~/Downloads/statReport")
threeLoc = grep("3PS",commArea$AREAID)
wheree = commArea$AREAID[threeLoc]
commArea <- commArea[threeLoc,]
