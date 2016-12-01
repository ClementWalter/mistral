## -----------------------------------------------------------------------------
## Fonction clusterize
## -----------------------------------------------------------------------------
##    Copyright (C) 2013
##    Developpement : C. WALTER
##    CEA
## -----------------------------------------------------------------------------

clusterize = function(data,Ncluster,inMargin=TRUE) {
  
  data <- as.matrix(data)
  if(Ncluster>dim(data)[2]){
    message('More cluster that data points requested, all samples are kept')
    kmean_tmp = data
  }
  else {
    kmean_tmp = kmeans(t(data), centers=Ncluster, iter.max=30)
    if(inMargin==TRUE) {
      cluster = kmean_tmp$cluster
      cluster.sort = sort(cluster,index.return=TRUE)
      cluster.size = kmean_tmp$size
      index = c(1:Ncluster)
      for(iclust in 1:Ncluster) {
        index[iclust] = sample(1:cluster.size[iclust],1)+sum(cluster.size[0:(iclust-1)])
      }
      index = cluster.sort$ix[index]
      kmean_tmp = data[,index]
    }
    else {kmean_tmp = as.matrix(t(kmean_tmp$centers))}
  }
  
  return(kmean_tmp)
}