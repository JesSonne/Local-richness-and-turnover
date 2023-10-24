############################################## calculating overlap probability
setwd("C:/Users/hzc973/Downloads/Local-richness-and-turnover-main/Local-richness-and-turnover-main")

#Functions used for calculating similarities between two groups with respect to binary, categorical and continuous data types
#All functions were retrived from Geange, Shane W., et al. "A unified analysis of niche overlap incorporating data of different types." Methods in Ecology and Evolution 2.2 (2011): 175-184.
source("R functions.R") 

#loading summarized data for each mountain region
loc=read.csv("Data files/data per mountain region.csv")

# Define measurements of local species richness and turnover (i.e. alpha and beta)

## empirical patterns 
gamma=loc$total_richness
pos_rm=which(loc$total_richness<100) # excluding mountain regions with less than 100 species
beta=loc$Species.turnover
alpha=loc$local_richness
gamma[pos_rm]=NA
alpha[pos_rm]=NA
beta[pos_rm]=NA
pattern="Raw"

###OR

#results from spatial rarefaction analyses
gamma=rowMeans(read.csv("Data files/std_regional species richness_50gridcells.csv")[,-1])
beta=rowMeans(read.csv("Data files/std_species turnover_50gridcells.csv")[,-1])
alpha=rowMeans(read.csv("Data files/std_local richness_50gridcells.csv")[,-1])
pattern="standardised"

###OR

#results from spatial rarefaction analyses using alternative metrics for beta diversity
gamma=rowMeans(read.csv("Data files/std_regional species richness_50gridcells.csv")[,-1])
beta=rowMeans(read.csv("Data files/std_Sorensen dissimilarity_50gridcells.csv")[,-1])
alpha=rowMeans(read.csv("Data files/std_local richness_50gridcells.csv")[,-1])
pattern="Sørensen & standardised"


# define proportion used to shortlist mountain regions with the highest local species richness and turnover (10% used in the study)
prop=0.1
nn=round(length(which(!is.na(alpha)))*prop)

nr1=order(beta,decreasing = T)[1:nn]
nr2=order(alpha,decreasing = T)[1:nn]
nr3=c(nr1,nr2)[which(duplicated(c(nr1,nr2)))]
nr2=nr2[which(!nr2 %in% nr3)]
nr1=nr1[which(!nr1 %in% nr3)]

######################Testing for geographical overlap between the 10% of mountain regions with highest local species richness and turnover
dat=data.frame(geo=rep(0,nrow(loc)),lab="none")
dat$lab[c(nr1)]="beta"
dat$lab[c(nr2)]="alpha"
dat$geo[c(nr1)]=1
dat=rbind(subset(dat,dat$lab=="alpha"),subset(dat,dat$lab=="beta"))
dat=rbind(dat,data.frame(geo=rep(1,length(nr3)),lab="alpha"))
dat=rbind(dat,data.frame(geo=rep(0,length(nr3)),lab="beta"))


#Quantifying overlap for the empirical data using the equation (eqn.) 1 applicable for binary data (Geange et al. 2011) 
obs_overlap=no.bin.fn(sp = dat$lab,yy = dat$geo)[1,2]

#randomly swapping mountain regions and recalculating the overlap
set.seed(123) 
null=rep(NA,nrow = 1000)
for(i in 1:1000){
  sam=sample(dat$lab)
  null[i]=no.bin.fn(sp = sam,yy = dat$geo)[1,2]
}

#calculating the null model permutation p-value (P)
P_geo=1-length(which(null>obs_overlap))/1000
P_geo


#######################Testing for similarities in climate conditions between the 10% of mountain regions with highest local species richness and turnover
dat=data.frame(MAT=loc$MAT,MAP=loc$MAP,lab="none") # data on mean annual temperature (MAT) and annual precipitation (MAP)
dat$lab[nr1]="beta"
dat$lab[nr2]="alpha"
dat=subset(dat,!dat$lab=="none")


#######################calcularing the similarity in Temperature and precipitation separately following eqn. 3 in Geange et al 2011
NO_emp1=no.cts.fn(sp = dat$lab,yy = dat$MAT)[1,2]
NO_emp2=no.cts.fn(sp = dat$lab,yy = dat$MAP)[1,2]

#Then taking the mean following eqn. 5 applicable for continuous data (Geange et al. 2011)
obs_clim=mean(NO_emp1,NO_emp2)

#randomly swapping mountain regions and recalculating the similarities in climate
set.seed(123)
null=matrix(NA,nrow = 1000,ncol=2)
for(i in 1:1000){
 sam=sample(dat$lab)
 null[i,1]=no.cts.fn(sp = sam,yy = dat$MAT)[1,2]
 null[i,2]=no.cts.fn(sp = sam,yy = dat$MAT)[1,2]
}

#Then taking the mean following eqn. 5 in Geange et al 2011
null_mean=apply(null,1,mean)

#calculating the null model permutation p-value (P)
P_clim=1-length(which(null_mean>obs_clim))/1000
P_clim


######################Habitat analyses
#loading list with habitat types for each mountain regions, retrieved from Jung et al. (2020)
load("Data files/hab_list.RData")
clims=c(400,600,800,100,200,300,500)

#assembling frequency distribution of habitat types for the 10% of mountain regions with highest local species richness and turnover
ex_nr1=unlist(hab_list[nr1])
ex_nr2=unlist(hab_list[nr2])
ex_nr1=ex_nr1[which(ex_nr1<1400)]
ex_nr2=ex_nr2[which(ex_nr2<1400)]
ex_nr1=round(ex_nr1,-2)
ex_nr2=round(ex_nr2,-2)
ex_un1=c(table(ex_nr1))
ex_un2=c(table(ex_nr2))
hab_dat=data.frame(freq_beta=c(ex_un1),freq_alpha=ex_un2)
hab_dat=apply(hab_dat,2,function(x){x/sum(x)})
hab_dat=t(hab_dat)

#calculating the similarity in habitat compositions following eqn. 2 applicable for categorical data (Geange et al. 2011)
obs_hab=no.cat.fn(hab_dat)[1,2] 

#randomly swapping mountain regions and recalculating the similarities in habitat composition
set.seed(123)
comb=c(nr1,nr2)
null=rep(NA,1000)
for(y in 1:1000){
  
  nr1_null=sample(comb,length(nr1))
  nr2_null=comb[which(!comb %in% nr1_null)]
  ex_nr1=unlist(hab_list[nr1_null])
  ex_nr2=unlist(hab_list[nr2_null])
  ex_nr1=ex_nr1[which(ex_nr1<1400)]
  ex_nr2=ex_nr2[which(ex_nr2<1400)]
  ex_nr1=round(ex_nr1,-2)
  ex_nr2=round(ex_nr2,-2)
  ex_un1=c(table(ex_nr1))
  ex_un2=c(table(ex_nr2))
  
  hab_dat=data.frame(freq_beta=c(ex_un1),freq_alpha=ex_un2)
  hab_dat=apply(hab_dat,2,function(x){x/sum(x)})
  hab_dat=t(hab_dat)


  null[y]=no.cat.fn(hab_dat)[1,2]
  print(y)
  
}

#calculating the null model permutation p-value (P)
P_hab=1-length(which(null>obs_hab))/1000
P_hab



###Investigating Explanatory factors for the highest regional levels of avian species richness using pairwise two-tailed Mann-Whitney U tests 

#definding groups among the 10% of mountain regiosn with highest regional species richness
group=rep("Others",nrow(loc))
nr4=order(gamma,decreasing = T)[1:nn]
group[nr4]="top 10% total richness"
group[nr4[which(nr4 %in% c(nr1,nr3))]]="top 10% total richness and turnover"
loc1=loc
loc1[pos_rm,]=NA

pairwise.wilcox.test(loc1$TopographicComplexity, group, p.adjust.method = "holm")
pairwise.wilcox.test(loc1$Clim_vol_convex, group, p.adjust.method = "holm")
pairwise.wilcox.test(loc1$NPP, group, p.adjust.method = "holm")
pairwise.wilcox.test(loc1$average_ADF_vol, group, p.adjust.method = "holm")


