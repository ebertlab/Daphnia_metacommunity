#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# 
# ## local arguments ##
# args = rep(0,3)
# args[1] = "M"
# args[2] =  2
# args[2] = 250

# species = as.character(args[1])
model = 'model_15b'
# nch = 3
nbin = as.numeric(as.character(args[1]))
niter = as.numeric(as.character(args[2]))
#####################

library(R2jags)

#  path_data <- "../3._Data/5b._Data_clean/" 
# path_data <- "./3._data/"
path_data <- "./data/"

#  path_model <- "../4._Analyses/1._Models/" 
# path_model <- "./1._models/"
path_model = "./model/"



get_occupancies<-function(){
  occup = list()
  occup[[1]] = readRDS(paste(path_data,"occupancies_longispina_82-17_core.RDS",sep=''))
  occup[[2]] = readRDS(paste(path_data,"occupancies_magna_82-17_core.RDS", sep = ''))
  occup[[3]] = readRDS(paste(path_data,"occupancies_pulex_82-17_core.RDS", sep =''))
  occup
}


get_revisits <- function(){
  revisit = list()
  revisit[[1]] = readRDS(paste(path_data,"revisit.longispina.RDS", sep =''))
  revisit[[2]] = readRDS(paste(path_data, "revisit.magna.RDS", sep =''))
  revisit[[3]] = readRDS(paste(path_data, "revisit.pulex.RDS", sep = ''))
  revisit
}


############ Run model 15 #######################

wet_dry_matrices_tmp = readRDS(paste(path_data,"dry_wet_82-17_core.RDS", sep =''))
distance.euclidean_tmp = readRDS(paste(path_data, "distances_matrix.RDS", sep =''))

# Keep only two first axes
cov_tmp = readRDS(paste(path_data, "environment_PCA.RDS", sep =''))[,1:2]

groups = readRDS(paste(path_data, "islandsGroups.RDS", sep = ''))

nsites <- unlist(lapply(groups, length))
nyears <- 36
nislands = length(groups)
  
#Subsetting matrices (by islands)
wet_dry_matrices = structure(rep(NA, 2 * length(groups) * max(nsites) * nyears), .Dim = c(length(groups), 2, max(nsites), nyears))
distance.euclidean = structure( rep(NA, length(groups) * max(nsites) * max(nsites)), .Dim = c(length(groups) , max(nsites) , max(nsites)))
cov = structure(rep(NA, length(groups) * dim(cov_tmp)[2] * max(nsites)), .Dim = c(length(groups), max(nsites), dim(cov_tmp)[2]))

occupancy <- get_occupancies()
revisit_tmp <- get_revisits()
for(i in 1:3) { revisit_tmp[[i]][[2]]$`2010` <- as.integer(revisit_tmp[[i]][[2]]$`2010`) }

occ <- structure(rep(NA, length(groups)*3*2*max(nsites)*nyears),.Dim = c(length(groups), 3, 2, max(nsites), nyears))

revisit <- structure(rep(NA, length(groups) * 3 * 2 * max(nsites)  * 4),.Dim = c(length(groups), 3, 2, max(nsites), 4))

for(i in 1:length(groups)){
  wet_dry_matrices[i, 1, 1:nsites[i], ] = wet_dry_matrices_tmp[[1]][ groups[[i]] ,]
  wet_dry_matrices[i, 2, 1:nsites[i], ] = wet_dry_matrices_tmp[[2]][ groups[[i]] ,]
  
  # Considering infinite distance between a pond and itself, such that own occupancy is not take into acount for colonization  
  distance.euclidean[i,1:nsites[i],1:nsites[i]] = distance.euclidean_tmp[groups[[i]], groups[[i]]] + (diag(nsites[i])*1000000)
  
  cov[i, 1:nsites[i], ] = cov_tmp[groups[[i]], ]  
  
  for(sp in 1:3){
    occ[i, sp, 1, 1:nsites[i], ] = occupancy[[sp]][[1]][groups[[i]], ]
    occ[i, sp, 2, 1:nsites[i], ] = occupancy[[sp]][[2]][groups[[i]], ]
  
    revisit[i, sp, 1, 1:nsites[i], ] = as.matrix(revisit_tmp[[sp]][[1]][groups[[i]], ])
    revisit[i, sp, 2, 1:nsites[i], ] = as.matrix(revisit_tmp[[sp]][[2]][groups[[i]], ])
  }
}


mask <- (diag(nrow = 3)-1)^2

parameters <- c("alpha","beta","d", "delta", "gamma_e", "gamma_c")

datax <- list(obs = occ,
              revisit = revisit,
              X=wet_dry_matrices,
              distM=distance.euclidean,
              mask = mask,
              cov = cov, 
              nsites = nsites,
              nislands = nislands,
              nyears = nyears)

for(g in 1:length(groups)){
  for(sp in 1:3){
    for(p in 1:2){
      for(i in 1:nsites[g]){
        for(j in 1:4){
          if(!is.na(revisit[g,sp,p,i,j]) & revisit[g,sp,p,i,j]==1 ){
            occ[g,sp,p,i,27+j] <- 1
          }
        }
      }
    }
  }
}

occ.inits1 = occ

for(i in 1:nislands){
  for(s in 1:3){
    for(j in 1:2){
      for(k in 1:nsites[i]){
        for(l in 1:nyears){
          if(is.na(occ.inits1[i,s,j,k,l])){
            occ.inits1[i,s,j,k,l] = sample(c(0,1),1)
          }
        }
      }
    }
  }
}


occ.inits1[which(occ==0)] <- round(runif(which(occ==0)))


inits <- list(list(alpha = runif(3,0.005,0.3),
                   beta = structure(rnorm(2*3*2), .Dim = c(2,3,2)),
                   delta = structure(rnorm(2*3*4), .Dim = c(2,3,4)),
                   gamma_e = structure(rnorm(2*3*3), .Dim = c(2,3,3)),
                   gamma_c = structure(rnorm(2*3*3), .Dim = c(2,3,3)),
                   occ = occ.inits1))

modelPATH = paste(path_model,model,sep='')

prt = proc.time()
jmodel <- jags(model.file = modelPATH,
     data = datax,
     n.chains = 1,
     inits = inits,
     n.iter = niter,
     n.burnin = nbin,
     n.thin = 2,
     parameters.to.save = parameters)
tim = proc.time()-prt
tim

idRAND = round(runif(1,100000,100000000))

saveRDS(jmodel,paste("results/jsample.full.allsp.model_15_",idRAND,".RDS", sep =''))

