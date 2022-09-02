#'#############################################################################################################
#' Script and functions to run:
#' * compute probabilities of sharing parasites for co-occurring hosts 
#'
#'  Ignacio Morales-Castilla, et al.
#'  started October 2018
#'#############################################################################################################



## to start
rm(list=ls())
options(stringsAsFactors=FALSE)


## load packages
packs.to.extract<-list('raster','ncdf4','maptools','sp','foreach','dismo',
                       'doParallel','abind','stringr','rgdal','foreign','SDMTools','RColorBrewer'
                       ,'dismo',"letsR","rgeos","rworldmap","randomForest")
lapply(packs.to.extract,require, character.only=T)





#### STEP 1. Load all data  ####

## load modelled ungulate distributions
load("data/RFhost_fut_curr_distrbs.RData")
load("data/ungulates.clim.fut.RData")


## data on host distributions
#h_distributions<-current.future.dists

## load co-occurrence meta-network matrices
load("data/coocurr_networks.RData")


## load host-parasite link data 


## load data on interactions - already curated from GMPD - all ungulates
hpint <- read.table("data/parasite_host_interactions.txt", sep="\t", header = TRUE)

## get names hosts
hostsnames<-unique(hpint[,3])

## get names parasites
parasnames<-unique(hpint[,2])


## hosts similarity data - all ungulates
hh_similarity <- read.table("data/host_host_overlap_distance.txt", sep="\t", header = TRUE)
head(hh_similarity)


## hosts in NAm (with distribution data)
hostsinNAM<-colnames(ungulates.clim.fut)[3:15]


#### STEP 2. Generate interaction matrices  ####

## for all ungulates
interaction.matrix<-array(NA, dim=c(length(hostsnames),length(parasnames)))
colnames(interaction.matrix)<-parasnames
rownames(interaction.matrix)<-hostsnames
calculate=T
if(calculate){
for(j in 1:length(parasnames)){
  print(j)  
  parasj<-parasnames[j]
  parasjints<-hpint[which(hpint[,2]==parasj),]
  for(i in 1:length(hostsnames)){
    hosti<-hostsnames[i]
    if(hosti%in%parasjints[,3]){
      interaction.matrix[i,j]<-1
    } else {
      interaction.matrix[i,j]<-0
    }  
  }
}
}

#heatmap(interaction.matrix)
image(interaction.matrix,col=c("white","black"))
dim(interaction.matrix)

## for NAm ungulates
## prunning and plotting interaction matrix for NAM
interaction.matrix.nam<-interaction.matrix[which(rownames(interaction.matrix)%in%hostsinNAM),]
interaction.matrix.nam<-interaction.matrix.nam[,which(colSums(interaction.matrix.nam)>0)]
dim(interaction.matrix.nam)
image(interaction.matrix.nam,col=c("white","black"))



#### STEP 3. Model parasite-sharing probability based on phylo.dist  ####

## function to rescale vector (0-1)
normalit<-function(m){return 
  (m - min(m))/(max(m)-min(m))}

## for parasite jth and modelled across hosts - FUNCTION
parasite.sharing <- function(hh_sim, int.mat, ngroups){

#hh_sim <- hh_similarity
#int.mat <- interaction.matrix
#hh_sim <- hh_similarity.nam
#int.mat <- interaction.matrix.nam
#ngroups <- 3  
resps <- array(NA, dim=c(nrow(hh_sim),2)) #npairs x npars x 2 (resp, atleast1)  
int.mat <- int.mat[,-which(colSums(int.mat[,])<ngroups)]
shrdint.pair.host.probs <- array(NA,dim=c(nrow(hh_sim),ncol(int.mat),4))
eval.ints <- array(NA, dim=c(7,ncol(int.mat),2))
row.names(eval.ints) <- c("Geo","Phy","Trt","R²","AUC","Kappa","tr")
colnames(eval.ints) <- colnames(shrdint.pair.host.probs) <- colnames(int.mat)

## loop across each parasite jth
for(j in 1:ncol(int.mat)){#j=46
print(paste(j,colnames(int.mat)[j]))  
  
    for(rowh in 1:nrow(hh_sim)){
      h1<-hh_sim[rowh,1]
      h2<-hh_sim[rowh,2]
      #print(paste(j,h1,h2))
      resps[rowh,1] <- ifelse(sum(int.mat[c(h1,h2),j])>1,1,0)
      resps[rowh,2] <- ifelse(sum(int.mat[c(h1,h2),j])>0,1,0)
      }
  
  ## modelling shared interactions
  dats <- data.frame(Lij=resps[,1],
                     Geo=scale(hh_sim[,6]),
                     Phy=scale(hh_sim[,8]),
                     Trt=scale(hh_sim[,10]))#,
                     #In1=resps[,2])
  
  absvals <- dats[dats$Lij==0,] 
  presvals <- dats[dats$Lij==1,]
  
  ## run several iterations for each model
  list.probs <- list()
  list.tr <- list()
  nrun <- 50
  model.evals <- array(NA, dim= c(7,nrun))
  
  for (reps in 1:nrun){#reps=1
    print(reps)
    group <- kfold(presvals, ngroups)
    pres_train <- presvals[group != 1, ] 
    pres_test <- presvals[group == 1, ] 
    
    if(nrow(absvals)>nrow(presvals)*2){
    absvalstr <- absvals[sample(1:nrow(absvals),round(nrow(presvals)*1.5,0)),]
    } else {
      absvalstr <- absvals
    }
    group.abs <- kfold(absvalstr, ngroups)
    abs_train <- absvals[group.abs != 1, ]
    abs_test <- absvals[group.abs == 1, ]
    
    train <- rbind(pres_train, abs_train)
    test <- rbind(pres_test, abs_test)
  
  mod1 <- glm(Lij~.,data=train,family=binomial)
  #mod1 <- randomForest(Lij~.,data=train,na.action=na.omit,
   #                  mtry=2,ntree=300,importance=T)
  
  
  ## predict on test data
  preds.test = predict(mod1,newdata=dats[2:4],type="response")
  
  
  ## evaluate and threshold
  threses = optim.thresh(obs = dats[,1], pred = preds.test)
  thres <- threses$`max.sensitivity+specificity`[1]
  eval <- evaluate(pres_test, abs_test, mod1)
  
  thres.pred <- preds.test
  thres.pred[preds.test<thres] <- 0

  list.probs[[reps]] <- preds.test
  list.tr[[reps]] <- normalit(thres.pred)
  
  
    model.evals[1:3,reps] <- coefficients(mod1)[2:4]
    model.evals[4,reps] <- 1-mod1$deviance/mod1$null.deviance
    model.evals[5,reps] <- eval@auc
    model.evals[6,reps] <- max(eval@kappa)
    model.evals[7,reps] <- thres
    
  
  }

  ## summarize model results
  eval.ints[,j,1] <- apply(model.evals,1,mean,na.rm=T) 
  eval.ints[,j,2] <- apply(model.evals,1,sd,na.rm=T) 
  
  
  shrdint.pair.host.probs[,j,1] <- apply(do.call(cbind,list.probs),1,mean,na.rm=T)
  shrdint.pair.host.probs[,j,2] <- apply(do.call(cbind,list.tr),1,mean,na.rm=T)
  shrdint.pair.host.probs[,j,3] <- apply(do.call(cbind,list.probs),1,sd,na.rm=T)
  shrdint.pair.host.probs[,j,4] <- apply(do.call(cbind,list.tr),1,sd,na.rm=T)
  
  
  }

return(list(evaluations=eval.ints,probabilities=shrdint.pair.host.probs))

}


## apply function first for all ungulates

mod.probs.allung <- parasite.sharing(hh_similarity, interaction.matrix, 3)

#mod.probs.allung$evaluations <- mod.probs.allung$evaluations[,1:272,]

mod.probs.allung$evaluations[4:7,,1]


#save(mod.probs.allung,file = "~/modelsHP.all.ungulates2.RData")
int.mat <- interaction.matrix
ngroups <- 3  
int.mat <- int.mat[,-which(colSums(int.mat[,])<ngroups)]


## load modelled probabilities ####
load("/modelsHP.all.ungulates2.RData")

evals.allungmean <- mod.probs.allung$evaluations[,1:320,1]
evals.allungsd <- mod.probs.allung$evaluations[,1:320,2]
colnames(evals.allungmean) <- colnames(evals.allungsd) <- colnames(int.mat)

t(evals.allungmean)
t(evals.allungsd)
 
## apply function for NAm ungulates only
int.mat <- interaction.matrix.nam
ngroups <- 3  
int.mat <- int.mat[,-which(colSums(int.mat[,])<ngroups)]

hostsnam <- row.names(interaction.matrix.nam)
hh_similarity.nam <- subset(hh_similarity,Host.1%in%hostsnam & Host.2%in%hostsnam)
parasites.nam <- colnames(interaction.matrix.nam)

mod.probs.allung.nam <- parasite.sharing(hh_similarity.nam, 
                                         interaction.matrix.nam,
                                         ngroups=3)

#save(modelsHP.NAm.ungulates.RData")
load("modelsHP.NAm.ungulates.RData")


evals.namungmean <- mod.probs.allung.nam$evaluations[,1:58,1]
evals.namungsd <- mod.probs.allung.nam$evaluations[,1:58,2]
colnames(evals.namungmean) <- colnames(evals.namungsd) <- colnames(int.mat)

#modelsHP.NAmung.evalsmeans
t(evals.namungmean)


#modelsHP.NAmung.evalssds
t(evals.namungsd)


#### STEP 4. Plotting summary of evaluations and coefficients  ####

# representing evals for northam vs. all ungulates
mat.evals.nam <- t(mod.probs.allung.nam$evaluations[c(4,6,5),1:58,1])
mat.evals <- t(mod.probs.allung$evaluations[c(4,6,5),1:320,1])

mat.evals.all <- cbind(mat.evals[,1],mat.evals.nam[,1],rep(NA,nrow(mat.evals)),rep(NA,nrow(mat.evals)),
                       mat.evals[,2],mat.evals.nam[,2],rep(NA,nrow(mat.evals)),rep(NA,nrow(mat.evals)),
                       mat.evals[,3],mat.evals.nam[,3])
dev.off()
par(mfrow=c(1,1))
boxplot(mat.evals.all,ylim=c(0,1),ylab="Parasite-sharing accuracy",
        col=c("cyan4","goldenrod1"),cex.lab=1.3,
        pars = list(boxwex = 0.5, staplewex = 0.0, outwex = 0.0),
        yaxt="n",xaxt="n",outline=F,lty=1,
        staplewex=0,whisklty=1)
axis(1,c(1.5,5.5,9.5),c("D²","Kappa","AUC"))        
axis(2,seq(0,1,.2),seq(0,1,.2),las=2)        
legend(7,0.2,c("World","N. America"),c("cyan4","goldenrod1"))


# representing predictor coeffs for northam vs. all ungulates
mat.evals.nam <- t(mod.probs.allung.nam$evaluations[1:3,1:58,1])
mat.evals <- t(mod.probs.allung$evaluations[1:3,1:320,1])

mat.evals.all <- cbind(mat.evals[,1],mat.evals.nam[,1],rep(NA,nrow(mat.evals)),rep(NA,nrow(mat.evals)),
                       mat.evals[,2],mat.evals.nam[,2],rep(NA,nrow(mat.evals)),rep(NA,nrow(mat.evals)),
                       mat.evals[,3],mat.evals.nam[,3])

par(mfrow=c(1,2))
boxplot(mat.evals,ylab="Parasite-sharing coefficients",
        col=c("cyan4"),cex.lab=1.3,
        pars = list(boxwex = 0.5, staplewex = 0.0, outwex = 0.0),
        #yaxt="n",
        #xaxt="n",
        outline=F,lty=1,
        staplewex=0,whisklty=1)
abline(h=0,lty=2,col="grey")

boxplot(mat.evals.nam,ylab="Parasite-sharing coefficients",
        col=c("goldenrod1"),cex.lab=1.3,
        pars = list(boxwex = 0.5, staplewex = 0.0, outwex = 0.0),
        #yaxt="n",
        #xaxt="n",
        outline=F,lty=1,
        staplewex=0,whisklty=1)
abline(h=0,lty=2,col="grey")


#### STEP 5. Plotting interaction matrices with host probabilities ####

hostnames <- unique(c(hh_similarity.nam$Host.1,hh_similarity.nam$Host.2))
nsps <- nrow(interaction.matrix.nam)
npar <- dim(mod.probs.allung.nam$probabilities[,,2])[2]
sharing.prbs <- array(NA,dim=c(nsps,nsps,npar))
sharing.prbs.sd <- array(NA,dim=c(nsps,nsps,npar))

colnames(sharing.prbs) <- colnames(sharing.prbs) <- hostnames
colnames(sharing.prbs.sd) <- colnames(sharing.prbs.sd) <- hostnames

for(parasite in 1:npar){#parasite=5
for(i in 1:nsps){#i=1
  for(j in 1:nsps){#j=1
    h1 <- hostnames[i]
    h2 <- hostnames[j]
    print(paste(h1,h2,parasite))
    pos.i.j <- which(hh_similarity.nam$Host.1==h1 & hh_similarity.nam$Host.2==h2)
    
    
    if(i==j){
      sharing.prbs[i,j,parasite] <- 1
      sharing.prbs.sd[i,j,parasite] <- 0
    } else {
      if(length(pos.i.j)>0){
        sharing.prbs[i,j,parasite] <- mod.probs.allung.nam$probabilities[pos.i.j,parasite,2]
        sharing.prbs.sd[i,j,parasite] <- mod.probs.allung.nam$probabilities[pos.i.j,parasite,4]
        
        }
    }
    
  }
  
}
}


## plot an example for blue tongue disease

for(k in 1:58){print(paste(k,sum(sharing.prbs[,,k],na.rm=T)))}

m <- sharing.prbs[,,which(colnames(int.mat)=="Orbivirus_Bluetongue_virus")]
s <- sharing.prbs.sd[,,which(colnames(int.mat)=="Orbivirus_Bluetongue_virus")]
d <- as.matrix(cooccur[,,3])
dfut <- as.matrix(cooccur[,,4])

d[lower.tri(d)]<-NA
dfut[lower.tri(dfut)]<-NA

dev.off()
#hcl.pals()
par(mfrow=c(2,2),mar=c(1,1,1,1))
image(m*d, 
      axes=FALSE, useRaster=TRUE, col=hcl.colors(30))
image(m*dfut, 
      axes=FALSE, useRaster=TRUE, col=hcl.colors(30))
image(s*d, 
      axes=FALSE, useRaster=TRUE, col=hcl.colors(30,"Plasma"))
image(s*dfut, 
      axes=FALSE, useRaster=TRUE, col=hcl.colors(30,"Plasma"))

## save resulting matrices

save(sharing.prbs,file = "/parasite_sharing_prob_networks.RData")
save(sharing.prbs.sd,file = "/parasite_sharing_probsds_networks.RData")






#### STEP 6. Getting info about the networks and how they change ####

parasnames.nam <- gsub("_"," ",colnames(interaction.matrix.nam))

gmpd.taxa <- read.csv("data/GMPD taxonomy.csv")
gmpd.taxa <- subset(gmpd.taxa, ParasiteCorrectedName %in% parasnames.nam)
unique(gmpd.taxa$ParasiteCorrectedName)

countingpars <- data.frame(par.name=parasnames.nam)
countingpars$n <- colSums(interaction.matrix.nam)
countingpars$par.type <- ifelse(countingpars$par.name %in% gmpd.taxa$ParasiteCorrectedName,
                                gmpd.taxa$ParType,NA)

table(countingpars$par.type)





#### STEP 7. Obtaining network metrics ####
library(igraph)

## generate graph for network of coocurring mammals
HP_network_nam <- graph_from_incidence_matrix(int.mat)


## generate graph for network of sharing parasites
## it first need to reconstruct the matrices of pairwise parasite sharing
dim(hh_similarity.nam)

HP_network_nam <- graph_from_incidence_matrix(int.mat)

parshare.intmat <- array (NA,dim=c(nrow(hh_similarity.nam),ncol(int.mat),2))
colnames(parshare.intmat) <- colnames(int.mat)
row.names(parshare.intmat)<-paste(hh_similarity.nam[,1],hh_similarity.nam[,2],sep=".")
d <- as.matrix(cooccur[,,1])
dfut <- as.matrix(cooccur[,,2])
d[is.na(d)] <- 0 
dfut[is.na(dfut)] <- 0 

for(i in 1:ncol(int.mat)){#i=1
  mat.i <- sharing.prbs[,,i]
  row.names(mat.i)<-colnames(mat.i)
  
  
  PShcurr <- mat.i*d
  PShfut <- mat.i*dfut
  
  vec.i <- c(PShcurr[upper.tri(PShcurr)])  
  vec.i.fut <- c(PShfut[upper.tri(PShfut)])  
  
  parshare.intmat[,i,1] <- vec.i 
  parshare.intmat[,i,2] <- vec.i.fut 
  
    
  } 


# get networks
HPshare_nam <- graph_from_incidence_matrix(parshare.intmat[,,1],
                                                   weighted = T)
HPshare_nam_fut <- graph_from_incidence_matrix(parshare.intmat[,,2],
                                                       weighted = T)


dev.off()
#plot(HPshare_nam)
#plot(HPshare_nam_fut)


# get network centrality measurements (how connected each species is)
HPsh_centrality <- data.frame(
  degree=degree(HPshare_nam)[1:78],
  betweenness=betweenness(HPshare_nam)[1:78],
  eigen_centrality=eigen_centrality(HPshare_nam)$vector[1:78],
  closeness=closeness(HPshare_nam)[1:78]
)

HPsh_centrality_fut <- data.frame(
  degree=degree(HPshare_nam_fut)[1:78],
  betweenness=betweenness(HPshare_nam_fut)[1:78],
  eigen_centrality=eigen_centrality(HPshare_nam_fut)$vector[1:78],
  closeness=closeness(HPshare_nam_fut)[1:78]
)

colMeans(HPsh_centrality_fut)/
colMeans(HPsh_centrality)



#### STEP 8. Identify novel interactions ####

# check individual host species
#interaction.matrix.nam['Odocoileus_hemionus',]


shiftsinprobs <- parshare.intmat[,,1]-parshare.intmat[,,2]

shiftsinprobs <- shiftsinprobs[which(-0.01>apply(shiftsinprobs,1,min)),
                               which(-0.01>apply(shiftsinprobs,2,min))]


dim(shiftsinprobs)
shiftsinprobs<-as.data.frame(shiftsinprobs)

shiftsinprobs$host1 <- unlist(lapply(strsplit(row.names(shiftsinprobs),
                                              ".",fixed = T),function(x)x[1]))
shiftsinprobs$host2 <- unlist(lapply(strsplit(row.names(shiftsinprobs),
                                              ".",fixed = T),function(x)x[2]))

#novel_interactions

shiftsinprobs



## get truly new interactions

#novelints <- array(NA, dim=dim(interaction.matrix.nam))


getnovelints <- function(shiftsinprobs,interaction.matrix.nam,Pr){
novelints <- list()

for(i in 1:(ncol(shiftsinprobs)-2)){
  #i=13
  print(i)
  par.i <- colnames(shiftsinprobs)[i]
  ints.par.i <- interaction.matrix.nam[,par.i]
  
  novelints.j<-list()
  for(j in 1:8){#j=1
  print(j)
  prob.ij <-  shiftsinprobs[j,c(i,26:27)]
    
  host1.i <- shiftsinprobs[j,26]
  host2.i <- shiftsinprobs[j,27]
  
  if(Pr > prob.ij[1,1] & (ints.par.i[host1.i]==0 | ints.par.i[host2.i]==0)){
    
    if(ints.par.i[host1.i]==0 & ints.par.i[host2.i]==0){
    novelints.j[[j]] <- c(par.i,unlist(c(prob.ij[1,2:3])))}
    if(ints.par.i[host1.i]==0 & ints.par.i[host2.i]!=0){
      novelints.j[[j]] <- c(par.i,unlist(c(prob.ij[1,2])),NA)}
    if(ints.par.i[host1.i]!=0 & ints.par.i[host2.i]==0){
      novelints.j[[j]] <- c(par.i,unlist(c(prob.ij[1,3])),NA)}
    
  } else {
    novelints.j[[j]] <- rep(NA,3)
  }
  }
  
  sum.nov.ints <- do.call(rbind,novelints.j)
  to.keep <- apply(sum.nov.ints,1,function(x)ifelse(sum(is.na(x))==3,F,T))
  sum.nov.ints <- sum.nov.ints[to.keep,]
  
  if(is.null(dim(sum.nov.ints))){
  
  unihost <- sum.nov.ints[2]
    
  } else {
  
  unihost <- unique(c(sum.nov.ints[,2:3]))
  unihost <- unihost[!is.na(unihost)]
  
  }
  
  
  novelints[[i]] <- data.frame(host=unihost,parasite=rep(par.i,length(unihost)))
  
}
  

novelints <- do.call(rbind,novelints)
novelints <- novelints[order(novelints$host),]

return(novelints)
}

novel.interactions <- getnovelints(shiftsinprobs,interaction.matrix.nam,-0.1)


#Novel_ints
novel.interactions




## end
  