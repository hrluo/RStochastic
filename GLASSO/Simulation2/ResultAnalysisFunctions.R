##############################

#function that finds connected genes and reduces pmat dim
red.pmat.fn<-function(pmat.par){
  pmat.copy<-pmat.par
  for (i in 1:nrow(pmat.copy)){
    pmat.copy[i,i]<-0
  }
  n0offdiag<-which(abs(pmat.copy)>10^(-15),arr.ind = T)     #use abs and make cutoff value small!!!
  con.ind.num<-unique(n0offdiag[,1])
  con.index<-union.index[con.ind.num]
  red.pmat<-pmat.par[con.ind.num,con.ind.num]
  red.pmat[abs(red.pmat)<10^(-15)]<-0      #don't forget abs()
  row.names(red.pmat)<-con.index
  colnames(red.pmat)<-con.index
  result<-list("reduced.precision.matrix"=red.pmat,"non0offdiag.index"=n0offdiag)
  return (result)
}

#graph function
XFgraph.fn<-function(G,geneInPath,path.num){
  library(fastclime)
  if (nrow(G)>0){
    diag(G) = 0
    pathwayIndicator<-vector("list",path.num)
    Matrix(G, sparse = TRUE)
    g = graph.adjacency(as.matrix(G != 0), mode = "undirected", 
                        diag = F,add.colnames='label')
    V(g)$color<-"white"
    V(g)$size<-1
    for(i in path.num:1){
      pathwayIndicator[[i]]<-colnames(G) %in% geneInPath[[i]]
      V(g)[pathwayIndicator[[i]]]$color<-rainbow(path.num)[i]
      V(g)[pathwayIndicator[[i]]]$size<-path.num+2-i
    }
    
    
    # V(g)[pathwayIndicator]$color<-"gold"
    # V(g)[pathwayIndicator]$size<-3
    # V(g)[!pathwayIndicator]$color<-"blue"
    # V(g)[!pathwayIndicator]$size<-1
    layout.grid = layout.circle(g)
    par(mfrow = c(1, 1))
    par(mar=c(0.5,1,3.5,1))
    
    g2<-graph(edges=E(g)[1:2],n=length(V(g)),directed=F)
    
    if(dim(G)[1]>30){
      plot(g, layout = layout.grid, edge.color = "gray50", 
            vertex.label=V(g)$name,vertex.label.cex=0.6)
    } else {
      plot(g, layout = layout.grid, edge.color = "gray50", 
            vertex.label=V(g)$name,vertex.label.cex=0.8)
    }
    
    rm(g)
  } else {
    plot.new()
    mtext("Zero Dim Squared Matrix",side=1,line=-2)
  }
}

# function that finds the shared connected genes btw individuals
find.samepair.fn<-function(pmat.n0offdiag.par1,pmat.n0offdiag.par2){
  bind.n0offdiag<-rbind(pmat.n0offdiag.par1,pmat.n0offdiag.par2)
  duplicate.index<-which(duplicated(bind.n0offdiag)==T)
  if (length(duplicate.index)>0){
    duplicate.value<-bind.n0offdiag[duplicate.index,]
    unique.row<-unique(duplicate.value[,1])
    sh.mat<-matrix(0,nrow=length(unique.row),ncol=length(unique.row))
    row.names(sh.mat)<-unique.row
    colnames(sh.mat)<-unique.row
    for (i in 1:nrow(duplicate.value)){
      sh.mat[row.names(sh.mat)==duplicate.value[i,1],colnames(sh.mat)==duplicate.value[i,2]]<-1
    }
    unique.row.name<-union.index[unique.row]
    row.names(sh.mat)<-unique.row.name
    colnames(sh.mat)<-unique.row.name
  } else {
    sh.mat<-matrix(0,nrow=0,ncol=0)
  }
  return (sh.mat)
}





#function that finds connected genes and reduces pmat dim
easy.pmat.fn<-function(pmat.par,pmat.names,base.genes){
  pmat.copy<-pmat.par
  colnames(pmat.copy)<-pmat.names
  row.names(pmat.copy)<-pmat.names
  diag(pmat.copy)<-0
  
  base.index<-which(pmat.names %in% base.genes)
  connect.base.index<-unique(which(abs(pmat.copy[base.index,])>10^(-15),arr.ind = T)[,2])
  
  all.index<-union(base.index,connect.base.index)
  
  pmat.new<-pmat.copy[all.index,all.index]
  
  return (pmat.new)

}








XFgraph.new.fn<-function(G,geneInPath,path.num,oldG,pathwayName,basePathNum){
  library(fastclime)
  if (nrow(G)>0){
    diag(G) = 0
    pathwayIndicator<-vector("list",path.num)
    Matrix(G, sparse = TRUE)
    g = graph.adjacency(as.matrix(G != 0), mode = "undirected", 
                        diag = F,add.colnames='label')
    V(g)$color<-"gold"
    V(g)$size<-1
    for(i in path.num:1){
      pathwayIndicator[[i]]<-colnames(G) %in% geneInPath[[i]]
      V(g)[pathwayIndicator[[i]]]$color<-rainbow(path.num)[i]
    }
    
    V(g)[pathwayIndicator[[basePathNum]]]$color<-rainbow(path.num)[basePathNum]
    
    connect.freq<-table(which(abs(oldG)>10^(-15),arr.ind = T)[,1])
    big.genes<-connect.freq
    big.index<-as.numeric(names(big.genes))
    names(big.genes)<-colnames(oldG)[big.index]
    
    
    for(i in 1:length(big.genes)){
      tmp.index<-V(g)$label %in% names(big.genes)[i]
      V(g)[tmp.index]$size<-as.numeric(big.genes[i])*0.31+1
    }
    
    
    
    
    
    mylayout = layout_components(g)
    par(mfrow = c(1, 1))
    par(mar=c(0.5,2,3.5,2))
    
    V(g)$frame.color<-NA
    
    if(length(E(g))>6){
      plot(g, layout = mylayout, edge.color = "gray50", vertex.label.color="black",
           vertex.label=V(g)$name,vertex.label.cex=0.9)
    } else {
      V(g)$size<-V(g)$size*2
      plot(g, layout = mylayout, edge.color = "gray50", vertex.label.color="black",
           vertex.label=V(g)$name,vertex.label.cex=0.9,xlim=c(-2,2))
    }
    legend("bottomright",pch=16,legend=c(pathwayName[1:path.num],"Genes Not in the Enriched Pathways"),
           col=c(rainbow(path.num)[1:path.num],"gold"),cex = 0.6)
    
    rm(g)
  } else {
    plot.new()
    mtext("Zero Dim Squared Matrix",side=1,line=-2)
  }
}






