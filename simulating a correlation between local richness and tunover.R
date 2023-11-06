# Simmulations showing how a positive correlation beteween local species richness and turnover can emerge
#if there are big differences in regional species pool size among the regions to be compared (Tuomisto and Ruokolainen 2012)

require(terra)
require(parallel)

setwd("...Downloads/Local-richness-and-turnover-main")
source("R functions.R")


#prepare the bacground raster 
rr=rast(matrix(0,nrow=52,ncol=52))
ras_id=rr
ras_id[]=1:length(ras_id[])
coords=crds(ras_id)
sub_ras=ras_id
sub_ras[1,]=NA
sub_ras[,1]=NA
sub_ras[,52]=NA
sub_ras[52,]=NA


set.seed(123)
n_rep=100 # state the number of simmulations

#specify the frequency distribution of total species richness (i.e., following a log-normal distribution)
rich=round(rlnorm(n = n_rep,meanlog = 5,sdlog = 1))
hist(rich)

#specify the frequency distribution of species' range sizes (i.e., following a log-normal distribution)
rs=round(rlnorm(n = max(rich),meanlog = 5,sdlog = 0.5))
hist(rs)


#running simmulations using paralell computation. 
#might take a long time to run...
dat=data.frame(gamma=rep(NA,n_rep),alpha=NA,beta=NA)
for(i in 1:n_rep){
  #sample species richness and and their range sizes
  sp_rich=rich[i]
  ra=sample(rs,sp_rich,prob = rs,replace = F)
  
  #simmulate distributions (optional: specify the number of cores allocated)
  res=mclapply(X=ra,FUN = spred_dye,dis_mat = NULL,sf1 = sub_ras,mc.cores = 8)
  
    
  #calculate local species richness and species turnover
  res=unlist(res)    
  sub_vol=ras_id
  sub_vol[]=0
  id=as.numeric(names(table(res)))
  alpha=as.numeric(table(res))
  sub_vol[id]=alpha
  sub_vol[which(is.na(ras_id[]))]=NA
  #plot(sub_vol,axes=F)
  
  dat$gamma[i]=sp_rich
  dat$alpha[i]=mean(sub_vol[])
  dat$beta[i]=sp_rich/mean(sub_vol[])
  
  
  plot(dat$alpha,dat$beta-1,xlab="local species richness",ylab="Species turnover")
  print(i)
}





