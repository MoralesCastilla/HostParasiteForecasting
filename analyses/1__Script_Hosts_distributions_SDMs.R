#'#############################################################################################################
#' Script and functions to run:
#' * load and model current and future distribution for ungulate hosts 
#'
#'  Ignacio Morales-Castilla, et al.
#'  started October 2018
#'#############################################################################################################



## to start
rm(list=ls())
options(stringsAsFactors=FALSE)


## load packages
#install.packages("remotes")
#remotes::install_version("SDMTools", "1.1-221")


packs.to.extract<-list('raster','ncdf4','maptools','sp','foreach','dismo',
                       'doParallel','abind','stringr','rgdal','foreign','SDMTools','RColorBrewer'
                       ,'dismo',"letsR","rgeos","rworldmap","randomForest")

lapply(packs.to.extract,require, character.only=T)

#### STEP 1. Load data, transform and format data ####


## load source data for ungulate distributions
load("data/ungulates.lets.RData")


## explore data
plot(ungulates.lets)


## cropping data for Northamerica

## get worldmap and crop pres-abs object to Northamerica
worldmap <- getMap(resolution="high")
Northamerica <- worldmap[worldmap@data$NAME %in% c("United States","Canada","Mexico"),]
e <- c(-178.19450,  0.00000,  14.54539,   83.11612)
Northamerica <- crop(Northamerica,e)

ungulates.Nam <- lets.pamcrop(ungulates.lets,Northamerica)
ungulates.Nam$Richness_Raster <- crop(ungulates.lets$Richness_Raster,Northamerica)
plot(ungulates.Nam)
ungulate.species <- ungulates.Nam$Species_name


#### Add climate data from worldclim ####

## download current data
wclim<-getData('worldclim', var='bio', res=10)

## adding current climate variables
ungulates.clim<-lets.addvar(ungulates.Nam,wclim)
dim(ungulates.clim)


## getting future data (for CMIP5, RCP8.5)
wclim.fut.HD.85<-getData('CMIP5', var='bio', res=10, rcp=85, model='HD', year=70)

## adding future climate variables
projection(wclim.fut.HD.85) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
ungulates.clim.fut<-lets.addvar(ungulates.Nam,wclim.fut.HD.85)
ungulates.clim.fut<-cbind(ungulates.clim,ungulates.clim.fut[,16:34])
rm(ungulates.clim,wclim,wclim.fut.HD.85,Northamerica) ## free space


## remove collinear env.variables (in case we don't want to input all climatic variables)
cor.clim<-cor(ungulates.clim.fut[,16:34],use="pairwise.complete.obs")
drop.temp<-which(cor.clim[,"bio1_mean"]<1 & cor.clim[,"bio1_mean"]>0.8)
drop.prec<-which(cor.clim[,"bio12_mean"]<1 & cor.clim[,"bio12_mean"]>0.82)  


## defining environmental data
env.current<-ungulates.clim.fut[,16:34]
env.future<-ungulates.clim.fut[,35:53]
colnames(env.future)<-colnames(env.current)
env.current<-env.current[,-c(drop.temp,drop.prec)]
env.future<-env.future[,-c(drop.temp,drop.prec)]

## remove cells with NAs in both distributions and climate
cells.to.rem<-which(is.na(env.future[,1]))
ungulates.clim.fut<-ungulates.clim.fut[-cells.to.rem,]
env.current<-env.current[-cells.to.rem,]
env.future<-env.future[-cells.to.rem,]




#### model current and future distributions using randomForests ####


PA<-ungulates.clim.fut[,3:15]
current.future.dists<-array(NA,dim=c(nrow(PA),ncol(PA),8))
model.evaluations<-array(NA,dim=c(12,ncol(PA),2)) #R2, importance, thresh, AUC, Kappa
plotting=F

for (i in 1:ncol(PA)){#i=1
  print(i)
  target = colnames(PA)[i]
  
  sdmdata <- data.frame(pres.abs=PA[,target],env.current)
  sdmdata.fut <- data.frame(env.future)
  prevalence <- sum(sdmdata$pres.abs) 
  absvals <- sdmdata[sdmdata$pres.abs==0,] 
  presvals <- sdmdata[sdmdata$pres.abs==1,]
  
  ## run several iterations for each model
  list.probs.current <- list()
  list.probs.future <- list()
  list.tr.current <- list()
  list.tr.future <- list()
  list.evals <- list()
  nrun <- 10
  model.evals <- array(NA, dim= c(12,nrun))
  
  for (j in 1:nrun){#j=1
  print(paste(target,j))
  
  group <- kfold(presvals, 4)
  pres_train <- presvals[group != 1, ] 
  pres_test <- presvals[group == 1, ] 
  
  group.abs <- kfold(absvals, 4)
  abs_train <- absvals[group.abs != 1, ]
  abs_test <- absvals[group.abs == 1, ]
  
  train <- rbind(pres_train, abs_train)
  test <- rbind(pres_test, abs_test)
  
  
  ## calibrate models on train data
  m.curr.train<-randomForest(x=train[,2:9],y=train[,1] 
                               ,na.action=na.omit,ntree=500,importance=T) # env only
  
  eval <- evaluate(pres_test, abs_test, m.curr.train)
  pred.rf.curr <- predict(m.curr.train,sdmdata[2:9])
  pred.rf.fut <- predict(m.curr.train,sdmdata.fut)
  
  thres <- threshold(eval, 'spec_sens')
  distrib.binaria.curr <- pred.rf.curr > thres
  distrib.binaria.fut <- pred.rf.fut > thres
  
  
  ## store results
  ## saving modelled distributions w/o probabilities
  list.probs.current[[j]] <- pred.rf.curr
  list.probs.future[[j]] <- pred.rf.fut
  
  pred.rf.curr[pred.rf.curr<thres]<-0
  pred.rf.curr[pred.rf.curr>=thres]<-1
  pred.rf.fut[pred.rf.fut<thres]<-0
  pred.rf.fut[pred.rf.fut>=thres]<-1
  
  list.tr.current[[j]] <- pred.rf.curr
  list.tr.future[[j]] <- pred.rf.fut
  
  
  ## saving model evaluation
  model.evals[1,j] <- max(m.curr.train$rsq)
  model.evals[2:9,j] <- m.curr.train$importance[,1]
  model.evals[10,j] <- thres
  model.evals[11:12,j] <- c(eval@auc, max(eval@kappa))
  
  }
  
  ## summarize model results
  model.evaluations[,i,1] <- apply(model.evals,1,mean,na.rm=T) 
  model.evaluations[,i,2] <- apply(model.evals,1,sd,na.rm=T) 
  
  
  current.future.dists[,i,1] <- apply(do.call(cbind,list.probs.current),1,mean,na.rm=T)
  current.future.dists[,i,2] <- apply(do.call(cbind,list.probs.future),1,mean,na.rm=T)
  current.future.dists[,i,3] <- ifelse(apply(do.call(cbind,list.tr.current),1,sum,na.rm=T)>nrun*0.7,1,0)
  current.future.dists[,i,4] <- ifelse(apply(do.call(cbind,list.tr.future),1,sum,na.rm=T)>nrun*0.7,1,0)
  
  current.future.dists[,i,5] <- apply(do.call(cbind,list.probs.current),1,sd,na.rm=T)
  current.future.dists[,i,6] <- apply(do.call(cbind,list.probs.future),1,sd,na.rm=T)
  
  
 
}




#### retrieve and save results ####

## model evaluations
eval.means <- model.evaluations[,,1]
eval.sds <- model.evaluations[,,2]

rownames(eval.means) <- rownames(eval.sds) <- c("R²","BIO1","BIO2","BIO4","BIO7",
                                                "BIO8","BIO12","BIO15","BIO18",
                                                "tr",'AUC',"Kappa")
colnames(eval.means) <- colnames(eval.sds) <- colnames(PA)

t(eval.means)
t(eval.sds)



## saving host species current and future distributions as predicted by SDMs####
save(current.future.dists,file = "data/RFhost_fut_curr_distrbs.RData")

# plot some species 
i <- 3
ff <- current.future.dists[,i,1]
fff <- current.future.dists[,i,2]
gg <- current.future.dists[,i,3]
ggg <- current.future.dists[,i,4]


par(mfrow=c(2,2),mar=c(3,3,1,1))
  cols1<-hcl.colors(100,"RdYlBu",rev = T)[as.numeric(cut(ff,breaks = 100))]
  colsfut<-hcl.colors(100,"RdYlBu",rev = T)[as.numeric(cut(fff,breaks = 100))]
  cols1bin<-hcl.colors(100,"RdYlBu",rev = T)[as.numeric(cut(gg,breaks = 100))]
  colsfutbin<-hcl.colors(100,"RdYlBu",rev = T)[as.numeric(cut(ggg,breaks = 100))]
  
  plot(ungulates.clim.fut[,1],ungulates.clim.fut[,2],col=cols1)
  plot(ungulates.clim.fut[,1],ungulates.clim.fut[,2],col=colsfut)
  plot(ungulates.clim.fut[,1],ungulates.clim.fut[,2],col=cols1bin)
  plot(ungulates.clim.fut[,1],ungulates.clim.fut[,2],col=colsfutbin)


# plot some uncertainties 
  ff <- rowSums(current.future.dists[,,3],na.rm=T)
  fff <- rowSums(current.future.dists[,,4],na.rm=T)
  gg <- rowSums(current.future.dists[,,7],na.rm=T)
  ggg <- rowMeans(current.future.dists[,,8],na.rm=T)

  par(mfrow=c(2,2),mar=c(3,3,1,1))
  cols1<-hcl.colors(100,"RdYlBu",rev = T)[as.numeric(cut(ff,breaks = 100))]
  colsfut<-hcl.colors(100,"RdYlBu",rev = T)[as.numeric(cut(fff,breaks = 100))]
  cols1bin<-hcl.colors(100,"RdYlBu",rev = T)[as.numeric(cut(gg,breaks = 100))]
  colsfutbin<-hcl.colors(100,"RdYlBu",rev = T)[as.numeric(cut(ggg,breaks = 100))]
  
  plot(ungulates.clim.fut[,1],ungulates.clim.fut[,2],col=cols1)
  plot(ungulates.clim.fut[,1],ungulates.clim.fut[,2],col=colsfut)
  plot(ungulates.clim.fut[,1],ungulates.clim.fut[,2],col=cols1bin)
  plot(ungulates.clim.fut[,1],ungulates.clim.fut[,2],col=colsfutbin)
  
  
## save results as a raster
  fut.tr <- apply(current.future.dists[,,4],2,function(x){ifelse(x>0.69,1,0)})
  fut.uncert <- rowMeans(current.future.dists[,,8],na.rm=T)
  currentfut.richness <- cbind(ungulates.clim.fut[,1:2],
                            rowSums(ungulates.clim.fut[,3:15]),
                            apply(fut.tr,1,sum,na.rm=T),
                            fut.uncert,
                            fut.tr)
  
  curr.fut.rich <- rasterFromXYZ(currentfut.richness,
                             res=c(1,1),
                             crs=crs(ungulates.Nam$Richness_Raster))
  plot(curr.fut.rich)
  save(curr.fut.rich, file = "/curr.fut.rich.RData")
  
  
  #future.host.rich
  curr.fut.rich[[2]]
  
  #future.host.uncert
  curr.fut.rich[[3]]
  
  
## retrieve and plot parasite richness 
 
  dim(current.future.dists)[,,3]
  
  # retrieve thresholded distributions 
  ff <- rowSums(current.future.dists[,,3],na.rm=T)
  fff <- rowSums(current.future.dists[,,4],na.rm=T)
  gg <- rowSums(current.future.dists[,,7],na.rm=T)
  ggg <- rowMeans(current.future.dists[,,8],na.rm=T)
  
  
  
  
  
  
  
## end
  