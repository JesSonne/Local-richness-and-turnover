# R functions 

#functions to compute pairwise distances between grid cells in a mountaion region
distance_matrix=function(mountain_raster_trimmed=mountain_raster){
  cor=crds(mountain_raster_trimmed)
  cor=cbind(cor,extract(mountain_raster_trimmed,cor))
  dis_mat=matrix(0,nrow=nrow(cor),ncol=nrow(cor))
  for(i in 1:nrow(cor)){
    for(j in 1:nrow(cor)){
      a=(abs(cor[j,1]-cor[i,1]))^2
      b=(abs(cor[j,2]-cor[i,2]))^2
      dis_mat[i,j]=sqrt(a+b)
    }
  }
  
  result=list()
  result[[1]]=cor
  result[[2]]=dis_mat
  
  return(result)
}

# helper search function
find_id=function(coor){
  return(mountain_raster[coor])
}



#function to construct the geographical sampling fields used in the spatial rarefaction analyses
spred_dye=function(size, # number of grid cells to sample
                   dis_mat=distance_matrix(mountain_raster) #geographical distance matrix occationally used for redirecting the search algorithm
                   ){
  
  #set.seed(135)
  sf1=mountain_raster
  empt=sf1
  empt[which(!is.na(sf1[]))]=0
  start=sample(which(!is.na(sf1[])),1)
  empt[start]=1
  nas=which(is.na(sf1[]))
  
  path=vector()
  
  look_up=which(empt[]==1)
  result=c(find_id(look_up))[[1]]
  
  if(size >1){
    for(x in 1:(size-1)){
      empt1=empt
      empt1[]=0
      nr=which(empt[]==1)
      empt1[nr]=1
      
      
      neigh_mat=terra::adjacent(empt,cells=nr)
      neigh=c(neigh_mat)
      neigh=neigh[complete.cases(neigh)]
      
      
      rm=which(!neigh %in% nas)
      neigh=neigh[rm]
      rm=which(!neigh %in% nr)
      neigh=neigh[rm]
      
      
      #if grid cells are surrounded, they automatically gets sampled
      sur=table(neigh)
      if(length(neigh)>0 & 4 %in% sur){
        neigh=as.numeric(names(sur))[which(sur==4)[1]]
        #print("hej")
      }
      
      
      if(length(neigh)==0){
        cor=dis_mat[[1]]
        dis_mat1=dis_mat[[2]]
        ee1=mountain_raster[nr]
        ee1=ee1$top_q_1
        pos=which(cor[,3] %in% ee1)
        dd1=dis_mat1[pos,-pos]
        dd2=apply(dd1,2,min)
        val=cor[-pos,][which(dd2==min(dd2)),][,3][1]
        neigh=which(mountain_raster[]==val)
      }
      
      get=sample(c(neigh,neigh),1)
      empt[get]=1
      
      path=append(path,get)
      #plot(empt)
      
      #print(x) 
    }
    #return(path)
    look_up=which(empt[]==1)
    result=c(find_id(look_up))[[1]]
  } #end major if
  
  
  return(result)
}

