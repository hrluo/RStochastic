XFgraph.noPath.fn<-function(G){
  library(fastclime)
  if (nrow(G)>0){
    diag(G) = 0
    Matrix(G, sparse = TRUE)
    g = graph.adjacency(as.matrix(G != 0), mode = "undirected", 
                        diag = F,add.colnames='label')
    V(g)$color<-"gold"
    V(g)$size<-1

    
    connect.freq<-table(which(abs(G)>10^(-15),arr.ind = T)[,1])
    big.genes<-connect.freq
    big.index<-as.numeric(names(big.genes))
    names(big.genes)<-colnames(G)[big.index]
    
    
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
    
    rm(g)
  } else {
    plot.new()
    mtext("Zero Dim Squared Matrix",side=1,line=-2)
  }
}






