##############################################################################################################
#' Script and functions to run:
#' * load modelled distributions for ungulate hosts and get co-occurrences among hosts
#'
#'  Ignacio Morales-Castilla, et al.
#'  started October 2018
#'#############################################################################################################



## to start
rm(list=ls())
options(stringsAsFactors=FALSE)


## load packages
NA

#### STEP 1. Load data for hosts modelled distributions ####

## load modelled ungulate distributions
load("data/RFhost_fut_curr_distrbs.RData")
load("data/ungulates.clim.fut.RData")

PA <- ungulates.clim.fut[,3:15]
PAfut <- current.future.dists[,,4]
PAfut[PAfut>0.6]<-1
PAfut[PAfut!=1]<-0

hostsnames <- colnames(PA)
nsps <- length(hostsnames)


#### STEP 2. Obtain co-occurrence matrices current and future ####

## threshold binary model
cooccur.curr <- array(NA,dim=c(nsps,nsps))
cooccur.fut <- array(NA,dim=c(nsps,nsps))
colnames(cooccur.curr) <- colnames(cooccur.fut) <- hostsnames
rownames(cooccur.curr) <- rownames(cooccur.fut) <- hostsnames

for(i in 1:nsps){#i=1
  for(j in 1:nsps){#j=2
    print(paste(i,j))
    cooccur.curr[i,j] <- ifelse(max(rowSums(cbind(PA[,i],PA[,j])))>1,1,0)
    cooccur.fut[i,j] <- ifelse(max(rowSums(cbind(PAfut[,i],PAfut[,j])))>1,1,0)
    
    }
  
}


m <- as.matrix(cooccur.curr)
m[lower.tri(m)]<-NA
m.fut <- as.matrix(cooccur.fut)
m.fut[lower.tri(m.fut)]<-NA

par(mfrow=c(1,2))
image(m, 
      axes=FALSE, useRaster=TRUE, col=c("white","cyan4"))
image(m.fut, 
      axes=FALSE, useRaster=TRUE, col=c("white","cyan4"))


## probabilistic model computing summed P(Hi)*P(Hj)

## function to rescale predictions (0-1)
PA <- current.future.dists[,,3]
PAfut <- current.future.dists[,,4]

cooccur.curr.prb <- array(NA,dim=c(nsps,nsps))
cooccur.fut.prb <- array(NA,dim=c(nsps,nsps))
colnames(cooccur.curr.prb) <- colnames(cooccur.fut.prb) <- hostsnames
rownames(cooccur.curr.prb) <- rownames(cooccur.fut.prb) <- hostsnames

for(i in 1:nsps){#i=3
  for(j in 1:nsps){#j=3
    print(paste(i,j))
    mult.curr <- PA[,i]*PA[,j]
    mult.fut <- PAfut[,i]*PAfut[,j]
    
    table(PA[,i])
    mult.curr <- mult.curr[mult.curr > 0.01]
    mult.fut <- mult.fut[mult.fut > 0.01]
    
    if(i==j){
    cooccur.curr.prb[i,j] <- quantile(mult.curr)[5]
    cooccur.fut.prb[i,j] <- quantile(mult.fut)[5]
    } else {
    cooccur.curr.prb[i,j] <- quantile(mult.curr)[4]
    cooccur.fut.prb[i,j] <- quantile(mult.fut)[4]
      
    }
    
  }
  
}

m <- as.matrix(cooccur.curr.prb)
m[lower.tri(m)]<-NA
m.fut <- as.matrix(cooccur.fut.prb)
m.fut[lower.tri(m.fut)]<-NA

dev.off()
par(mfrow=c(1,2))
image(m, 
      axes=FALSE, useRaster=TRUE, col=hcl.colors(30))
image(m.fut, 
      axes=FALSE, useRaster=TRUE, col=hcl.colors(30))


## save resulting matrices
cooccur <- abind::abind(cooccur.curr,cooccur.fut,
                        cooccur.curr.prb,cooccur.fut.prb,
                        along=3)

save(cooccur,file = "data/coocurr_networks.RData")




## end
  