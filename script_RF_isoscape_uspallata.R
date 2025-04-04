################################################################################
############ RANDOM FOREST AND ENSEMBLE MCHINE LEARNING ANALYSES  ##############
################################################################################

# The script is adapated from the original script of Bataille CP, Crowley BE, Wooller MJ and GJ Bowen. 2020.
# Advances in global bioavailable strontium isoscapes. Palaeogeography, Palaeoclimatology, Palaeoecology, 555: 109849.
# DOI: https://doi.org/10.1016/j.palaeo.2020.109849


setwd("D:/.../RF_isoscape")

##########################SET LIBRARIES##################################################
  library(parallel)
  library(doParallel)
  library(raster)
  library(terra)
  library(sf)
  library(randomForest)
  library(readxl)
  library(exactextractr)
  library(caret)
  library(VSURF)
  library(ranger)
  library(dplyr)

  # Choose what is relevant for your OS
  os="unix"  # this includes osx on a Mac
  os="win"

  # Helper function to store (print) files in sub folders in different OS's
  os_path <- function(folder_path,filename){
    if (os=="win") return(paste(folder_path,"/",filename,sep="")) else
      return(paste(folder_path,"/",filename,sep="")) }

  #!! Create or use sub folders to keep things tidy!
   Out_Path <- "Output"
   Ras_Path <- "Projected_rasters"

    mydir  = getwd()
    dir.create(file.path(mydir,Ras_Path, fsep = .Platform$file.sep), showWarnings = FALSE)
    dir.create(file.path(mydir,Out_Path, fsep = .Platform$file.sep), showWarnings = FALSE)



################################################################################
############################## DATA PROCESSING #################################
################################################################################
                                     bioavailable_sr_dataset

# load bioavailable Sr database
sr_orig <- readxl::read_excel("bioavailable_sr_dataset.xlsx",col_names=TRUE, na="NA", sheet="data_sr")
sr_orig$X87Sr_86Sr<-as.numeric(sr_orig$X87Sr86Sr)
sr_orig$Latitude<-as.numeric(sr_orig$Latitude)
sr_orig$Longitude<-as.numeric(sr_orig$Longitude)

#remove duplicated data based on latitude, longitude and Sr values
sr_orig <- sr_orig %>%
  distinct(Latitude, Longitude, X87Sr86Sr, .keep_all = TRUE)
sr_orig_0<-sr_orig[!is.na(sr_orig$Latitude),]                #remove rows with missing XY coordinates
sr_orig_0<-sr_orig_0[!is.na(sr_orig0$Longitude),]
sr_orig_1<-sr_orig_0[!is.na(sr_orig_0$X87Sr86Sr),]         #remove rows with missing observed Sr data

# Project Sr data
sr_proj<-st_as_sf(sr_orig_1, coords=c("Longitude","Latitude"), crs=st_crs(4326))%>% 
                 st_transform(crs="+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")


#################################INPUT GLOBAL COVARIATE RASTERS#####################################################
###See Table 1 for raster list and references

###The projected covariates are available at: 
###https://drive.google.com/drive/folders/1g9rCGo3Kd3hz2o5JKkSbgNsGJclvsuQm?usp=sharing
###or on request to Mael Le Corre (mael.lecorre@abdn.ac.uk)
###Use this link and download the Projected_rasters folder as a .zip file into your working directory (35GO)
###Unzip the file using the 7zip program https://www.7-zip.org/ which can handle large .zip file.
###Verify that the unzip folder has a name of "Projected_rasters" as below and if necessary rename the folder to the correct name
###Once this is done the following raster should load in your script

  r.m1           =terra::rast("Projected_rasters\\srsrmed_qcupdate.tif")
  r.srsrq1       =terra::rast("Projected_rasters\\srsrq1_qcupdate.tif")
  r.srsrq3       =terra::rast("Projected_rasters\\srsrq3_qcupdate.tif")
  r.meanage_geol =terra::rast("Projected_rasters\\agemean_qcupdate.tif")
  r.minage_geol  =terra::rast("Projected_rasters\\agemin_qcupdate.tif")
  r.maxage_geol  =terra::rast("Projected_rasters\\agemax_qcupdate.tif")
  r.age          =terra::rast("Projected_rasters\\basement_age_reproj.tif")
  r.bouger       =terra::rast("Projected_rasters\\bouger_reproj.tif")
  r.elevation    =terra::rast("Projected_rasters\\elevation_reproj.tif")

  r.mat          =terra::rast("Projected_rasters\\mat_reproj.tif")
  r.map          =terra::rast("Projected_rasters\\map_reproj.tif")
  r.dust         =terra::rast("Projected_rasters\\dust_dep_reproj.tif")
  r.salt         =terra::rast("Projected_rasters\\seasalt_dep_reproj.tif")
  r.ai           =terra::rast("Projected_rasters\\ai_reproj.tif")
  r.pet          =terra::rast("Projected_rasters\\pet_reproj.tif")
  r.distance     =terra::rast("Projected_rasters\\distance.tif")
  r.dist<-resample(r.distance,r.srsrq1,method="near")                        #problem with r.dist extent -> need to be set to r.srsrq1 extent
  r.volc         =terra::rast("Projected_rasters\\volc_dep_reproj.tif")
  r.fire         =terra::rast("Projected_rasters\\fire_dep_reproj.tif")
  r.foss         =terra::rast("Projected_rasters\\fossilfuels_dep_reproj.tif")

  r.clay         =terra::rast("Projected_rasters\\r.clay_reproj.tif")
  r.ph           =terra::rast("Projected_rasters\\r.ph_reproj.tif")
  r.cec          =terra::rast("Projected_rasters\\r.cec_reproj.tif")
  r.bulk         =terra::rast("Projected_rasters\\r.bulk_reproj.tif")
  r.phkcl        =terra::rast("Projected_rasters\\phkcl.tif")
  r.ocs          =terra::rast("Projected_rasters\\r.ocs_reproj.tif")
  r.GUM          =terra::rast("Projected_rasters\\gum_mask3.tif")


# Function to extract covariables at the sampling sites or at the closest of the sampling site
# if the site fall within an empty cell

extract_nearest_non_na <- function(coordinates, raster) {
  # Convert the coordinates data.frame to an sf points object
  points <- st_as_sf(coordinates, coords = c("X", "Y"), crs = st_crs(raster))
  
  # Extract raster values under the points
  values <- extract(raster, points,method="simple")
  values<-values[,2]
  
  # Identify NA values
  na_indices <- is.na(values)
  
  # Find the indices of NA values
  na_indices <- which(na_indices)
  
  # Loop over each NA value
  for (i in na_indices) {
    # Get the current point
    point <- points[i, ]
    st_crs(point) <- st_crs(raster)
    # Create a buffer around the point with a 50km radius
    buffer <- st_buffer(point, dist = 50000)
    
    # Convert the buffer to the same CRS as the raster
    st_crs(buffer) <- st_crs(raster)
    
    # Crop the raster to the buffered extent
    cropped_raster <- crop(raster, buffer)
    
    # Extract the non-NA values within the buffer
    non_na_values <- cropped_raster[]
    non_na_values <- non_na_values[!is.na(non_na_values)]
    
    # Find the nearest non-NA value
    nearest_value <- non_na_values[which.min(st_distance(point, buffer))]
    
    # Assign the nearest non-NA value to replace the NA value
    values[i] <- nearest_value
  }
  
  # Return the extracted values
  return(data.frame(values))
}

# data extraction 
    m1xy<-extract_nearest_non_na(sr_proj,r.m1)
    srsrq1xy<-extract_nearest_non_na(sr_proj,r.srsrq1)
    srsrq3xy<-extract_nearest_non_na(sr_proj,r.srsrq3)
    meanage_geolxy<-extract_nearest_non_na(sr_proj,r.meanage_geol)
    minage_geolxy<-extract_nearest_non_na(sr_proj,r.minage_geol)
    maxage_geolxy<-extract_nearest_non_na(sr_proj,r.maxage_geol)
    agexy<-extract_nearest_non_na(sr_proj,r.age)
    bougerxy<-extract_nearest_non_na(sr_proj,r.bouger)
    elevationxy<-extract_nearest_non_na(sr_proj,r.elevation)

    mapxy<-extract_nearest_non_na(sr_proj,r.map)
    matxy<-extract_nearest_non_na(sr_proj,r.mat)
    dustxy<-extract_nearest_non_na(sr_proj,r.dust)
    saltxy<-extract_nearest_non_na(sr_proj,r.salt)
    aixy<-extract_nearest_non_na(sr_proj,r.ai)
    petxy<-extract_nearest_non_na(sr_proj,r.pet)
    distxy<-extract_nearest_non_na(sr_proj,r.distance)
    volcxy<-extract_nearest_non_na(sr_proj,r.volc)
    firexy<-extract_nearest_non_na(sr_proj,r.fire)
    fossxy<-extract_nearest_non_na(sr_proj,r.foss)

    clayxy<-extract_nearest_non_na(sr_proj,r.clay)
    phxy<-extract_nearest_non_na(sr_proj,r.ph)
    cecxy<-extract_nearest_non_na(sr_proj,r.cec)
    bulkxy<-extract_nearest_non_na(sr_proj,r.bulk)
    GUMxy<-extract_nearest_non_na(sr_proj,r.GUM)
    phkclxy<-extract_nearest_non_na(sr_proj,r.phkcl)
    ocsxy<-extract_nearest_non_na(sr_proj,r.ocs)

    sr_xy<- data.frame(st_coordinates(sr_proj))

### Append all extracted data and change names of column
  sr_proj_xy <- data.frame(sr_orig_1$ID,sr_orig_1$Latitude,sr_orig_1$Longitude, sr_orig_1$X87Sr86Sr,
                         sr_xy,m1xy,srsrq1xy,srsrq3xy,meanage_geolxy,minage_geolxy,maxage_geolxy,
                         agexy,bougerxy,elevationxy,
                         matxy,mapxy,dustxy,saltxy,aixy,petxy,distxy,volcxy,firexy,fossxy,
                         clayxy,phxy,cecxy,bulkxy,GUMxy,phkclxy,ocsxy)

  colnames(sr_proj_xy)<-c("ID","Latitude","Longitude","X87Sr86Sr",
                          "X","Y","r.m1","r.srsrq1","r.srsrq3","r.meanage_geol","r.minage_geol","r.maxage_geol",
                          "r.age","r.bouger","r.elevation",
                          "r.mat","r.map","r.dust","r.salt","r.ai","r.pet","r.dist","r.volc","r.fire","r.foss",
                          "r.clay","r.ph","r.cec","r.bulk","r.GUM","r.phkcl","r.ocs")

 write.table(sr_proj_xy,"sr_proj_xy_alldata.txt",sep="\t",col.names=T,row.names=F,quote=F)

# compute mean Sr value for each location
 sr_proj_xy2<-sr_proj_xy[!is.na(sr_proj_xy$r.m1),]                
 sr_agg1<-aggregate(sr_proj_xy2,by=list(sr_proj_xy2$Latitude,sr_proj_xy2$Longitude), FUN=median,na.rm=TRUE)
 sr_agg1<-cbind(sr_agg1[,c(7,8,6)],sr_agg1[,9:ncol(sr_agg1)])
# check correlation between covariates
 round(cor(sr_agg1[,-(1:3)]),2)
 sr_agg<-sr_agg1[,-c(4,7)]                                  #remove r.m1 and meanage_geol, because R>0.9



 write.table(sr_agg,"sr_agg_alldata.txt",sep="\t",col.names=T,row.names=F,quote=F)

################################################################################
############################# RANDOM FOREST ####################################
################################################################################


  sr_agg<-read.table("sr_agg_alldata.txt",sep="\t",h=T)

###Variable filtering using parallelized VSURF algorithm
  names(sr_agg)
  set.seed(1)

  sr_agg0<-sr_agg
  sr_agg1.vsurf<-VSURF(sr_agg0[,4:ncol(sr_agg0)],sr_agg0$X87Sr86Sr, RFimplem = "ranger", parallel = TRUE, ncores = detectCores() - 1, clusterType = "PSOCK")   
  sr_agg1.vsurf$varselect.pred
  sr_agg1.sub<-sr_agg0[,4:ncol(sr_agg0)]
  sr_agg1_VSURF <- sr_agg1.sub[c(sr_agg1.vsurf$varselect.pred)] 
  sr_agg1_VSURF<-cbind(sr_agg[,1:3],sr_agg1_VSURF)

  training_rf <- sr_agg1_VSURF[,-c(1,2)]

###Parallelize random forest modeling

  cluster <-parallel::makeCluster(detectCores() - 1) # convention to leave 1 core for OS
  registerDoParallel(cluster)

# Splitting the data for repeated cross validation
  fitControl <- trainControl(## 10-fold Crossvalidation
    method = "repeatedcv",
    number = 10,
    ## repeated ten times
    repeats = 5,
    verboseIter=FALSE ,
    returnResamp="final",
    savePredictions="all",
    # With parallel backend
    allowParallel=TRUE
    )

  bestmtry <- tuneRF(training_rf, training_rf$X87Sr86Sr, stepFactor=1, improve=1e-7, ntree=1000)
  mtry <- bestmtry[1]                                                           
  tunegrid <- expand.grid(.mtry=mtry)
  metric<-"Accuracy"

#variables selected: "r.srsrq1","r.maxage_geol","r.age","r.mat","r.dust","r.volc","r.foss","r.fire"
set.seed(1)
  RF_all <- train(X87Sr86Sr ~ ., data = training_rf, ntree=3000,method = "rf", importance=TRUE, tuneGrid=tunegrid,trControl= fitControl)

save(RF_all,file="RF_all")


# RMSE         Rsquared   MAE        
# 0.003403908  0.6803016  0.001771844


###variable importance and partial dependance plots
ly<-matrix(c(1,1,1,0,2,2,3,3,4,4,
             1,1,1,0,2,2,3,3,4,4,
             1,1,1,0,5,5,6,6,7,7,
             1,1,1,0,5,5,6,6,7,7,
             1,1,1,0,8,8,9,9,10,10,
             1,1,1,0,8,8,9,9,10,10),
             6,10,byrow=T)

layout(ly)
par(mar=c(5,4,2,0))
  varImpPlot(RF_all$finalModel,type=2,main="Variable Importance")
mtext("a)",cex=1.5,side=3,at=-0.0055)
mtext("b)",cex=1.5,side=3,at=0.04)
mtext(expression(paste(""^87,"Sr/"^86,"Sr predicted")),side=4,cex=1.5,line=8)
par(mar=c(4,2,1,1))

  partialPlot(RF_all$finalModel, training_rf, x.var = "r.srsrq1",main=NA,xlab="r.srsrq1",cex.axis=0.9,ylim=c(0.710,0.718))
  partialPlot(RF_all$finalModel, training_rf, x.var = "r.maxage_geol",main=NA,xlab="r.minage_geol",cex.axis=0.9,ylim=c(0.710,0.718))
  partialPlot(RF_all$finalModel, training_rf, x.var = "r.age",main=NA,xlab="r.age",cex.axis=0.9,ylim=c(0.710,0.718))
  partialPlot(RF_all$finalModel, training_rf, x.var = "r.mat",main=NA,xlab="r.mat",cex.axis=0.9,ylim=c(0.710,0.718))
  partialPlot(RF_all$finalModel, training_rf, x.var = "r.dust",main=NA,xlab="r.dust",cex.axis=0.9,ylim=c(0.710,0.718))
  partialPlot(RF_all$finalModel, training_rf, x.var = "r.volc",main=NA,xlab="r.volc",cex.axis=0.9,ylim=c(0.710,0.718))
  partialPlot(RF_all$finalModel, training_rf, x.var = "r.foss",main=NA,xlab="r.foss",cex.axis=0.9,ylim=c(0.710,0.718))
  partialPlot(RF_all$finalModel, training_rf, x.var = "r.fire",main=NA,xlab="r.fire",cex.axis=0.9,ylim=c(0.710,0.718))
  savePlot("Output\\importance_pplot.tif",type="tif")


###prediction map
#more filter
#"r.srsrq1","r.foss","r.volc","r.age","r.dust","r.fire","r.mat","r.maxage_geol"
stack_all<-c(r.srsrq1,r.foss,r.volc,r.age,r.dust,r.fire,r.mat,r.maxage_geol)
names(stack_all)<-c("r.srsrq1","r.foss","r.volc","r.age","r.dust","r.fire","r.mat","r.maxage_geol")


ex<-ext(-6500000,-5500000,-4750000,-3250000)
stack_all<-crop(stack_all, ex, snap='near')

rast.grid<-stack_all[[1]]/stack_all[[1]]

rf1 <- predict(stack_all, RF_all, ext=rast.grid, na.rm=TRUE, overwrite=TRUE)
writeRaster(rf1, filename="Output\\rf_all_pred.tif", overwrite=TRUE)

###spatial uncertainty using Quantile regression forest
training_qrf<-read.table("sr_agg_alldata.txt",sep="\t",h=T)
qrf_mod <- ranger(X87Sr86Sr ~ r.srsrq1+r.foss+r.volc+r.dust+r.maxage_geol+r.age+r.fire+r.mat,
                data = training_qrf, quantreg=TRUE, num.trees=3000,mtry=3)         #mtry same as for RF

# dataset need to be divided in several part to avoid memory problems
stack_all<-raster::stack(stack_all)
all_pxl<-as(stack_all,"SpatialPixelsDataFrame")
newdat0<-as.data.frame(all_pxl)
cc<-which(complete.cases(all_pxl@data))
lcc<-length(cc)
cc_div<-c(seq(0,lcc,100000),lcc)                                                #set size of the sub-dataset, need to be changed according available memory 

tab_pred<-as.data.frame(matrix(data=NA,ncol=3,nrow=0))
for (i in 1:(length(cc_div)-1))
    {
    print(paste(i,(length(cc_div)-1),Sys.time(),sep=" / "))
    cc_sub<-cc[(cc_div[i]+1):cc_div[i+1]]
    newdat1<-newdat0[cc_sub,]
    sr.rfd_low <- predict(qrf_mod, data=newdat1, type="quantiles", quantiles=0.159, fun = function(model, ...) predict(model, ...)$predictions)
    gc()
    sr.rfd_high <- predict(qrf_mod, data=newdat1, type="quantiles", quantiles=0.841, fun = function(model, ...) predict(model, ...)$predictions)
    gc()
    pred<-cbind(cc_sub,sr.rfd_low$predictions,sr.rfd_high$predictions)
    tab_pred<-rbind(tab_pred,pred)
    rm(sr.rfd_low)
    rm(sr.rfd_high)
    rm(pred)
    }
colnames(tab_pred)<-c("pxlnum","pred_low","pred_high")
tab_pred$error<-(tab_pred[,3]-tab_pred[,2])/2
all_pxl@data$se.sr<-NA
all_pxl@data$se.sr[cc]<-tab_pred$error
rast.se<-raster(all_pxl[,"se.sr"])
writeRaster(rast.se, filename="Output\\rf_all_sd.tif", format="GTiff", overwrite=TRUE)

######comparison predicted vs observed
RF_all<-load("RF_all")
isosr<-terra::rast("Output\\rf_all_pred.tif")

sampdb<-read.table("sr_proj_xy_alldata.txt",sep="\t",h=T)
samp_study<-sampdb[(nrow(sampdb)-53):nrow(sampdb),]
samp_proj<-st_as_sf(samp_study, coords=c("Longitude","Latitude"), crs=st_crs(4326))%>% 
                    st_transform(crs="+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs")
samp_xy<- data.frame(st_coordinates(samp_proj))

samp_pred<-extract_nearest_non_na(samp_xy,isosr)
tabcomp<-cbind(samp_xy,samp_study[,4],samp_pred)
colnames(tabcomp)[c(3,4)]<-c("sr_obs","sr_pred")

(rmse_all<-sqrt(mean((tabcomp$sr_obs-tabcomp$sr_pred)^2)))               



par(mfrow=c(1,2),mar=c(5,5,1,1))
plot(RF_all$pred$pred,RF_all$pred$obs,pch=15, cex=0.4, xlim=c(0.703,0.780),ylim=c(0.703,0.780),xlab=expression(paste(""^87,"Sr/"^86,"Sr predicted")),ylab=expression(paste(""^87,"Sr/"^86,"Sr samples")),
    cex.lab=1.5, cex.axis=1, las=1)
lm1<-lm(RF_all$pred$obs~RF_all$pred$pred)
abline(lm1)
text(x=0.7,y=0.78,"a)",cex=2,pos=4)

plot(tabcomp$sr_pred,tabcomp$sr_obs,pch=19,col="grey", xlim=c(0.703,0.711),ylim=c(0.703,0.711),xlab=expression(paste(""^87,"Sr/"^86,"Sr isoscape")),ylab=expression(paste(""^87,"Sr/"^86,"Sr samples")),
    cex.lab=1.5, cex.axis=1,las=1)
lm2<-lm(tabcomp$sr_obs~tabcomp$sr_pred)
abline(lm2,col="grey")

savePlot("Output\\correlation_plot.tif",type="tif")
 
#Call:
#lm(formula = tabcomp$sr_obs ~ tabcomp$sr_pred)
#
#Residuals:
#       Min         1Q     Median         3Q        Max 
#-2.053e-03 -4.714e-05  9.002e-05  3.295e-04  9.447e-04 
#
#Coefficients:
#                Estimate Std. Error t value Pr(>|t|)    
#(Intercept)     -0.04655    0.05841  -0.797    0.432    
#tabcomp$sr_pred  1.06568    0.08262  12.898 2.66e-13 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 0.0007103 on 28 degrees of freedom
#Multiple R-squared:  0.8559,    Adjusted R-squared:  0.8508 
#F-statistic: 166.4 on 1 and 28 DF,  p-value: 2.664e-13


#


  