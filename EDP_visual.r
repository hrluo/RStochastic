require(ggplot2)
require(gridExtra)
matrix2list<-function(mat){
  nRow=dim(mat)[1];
  nCol=dim(mat)[2];
  if(nRow!=nCol){
    return(FALSE)
  }else{
    n=nRow;
    lis<-matrix(NA,nrow=n*n,ncol=3);
    colnames(lis)<-c("row","col","value");
    counter=1;
    for(i in 1:n){
      for(j in 1:n){
        lis[counter,]<-c(i,j,mat[i,j]);
        counter=counter+1;
      }
    }
    return(lis)
  }
}
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
analyze_plot<-function(rho,x.data,y.data){
  partitions<-paste(t(rho)[,1],t(rho)[,2])
  nums<-seq(1,dim(rho)[2],1)
  visdf<-as.data.frame(cbind(nums,x.data,y.data,partitions))
  require(ggplot2)
  x.data.astext<-round(x.data[,1],digits = 1)
  for(col_sel in 2:(dim(x.data)[2]) ){
    x.data.astext<-paste(x.data.astext,round(x.data[,col_sel],digits=1),sep=",")
  }
  x.data.astext<-paste("(",x.data.astext,")",sep="")
  for(i in 1:(dim(visdf)[2]-1)){
    visdf[,i]<-as.numeric(
      as.character(visdf[,i])
      )
  }
  g3<-ggplot(data=visdf)+
      geom_text(aes(x=nums,y=y.data),
                label=as.character(x.data.astext),size=5,alpha=.5, angle=45)+
      geom_point(aes(x=nums,y=y.data,colour=partitions),size=2)+
      xlab("# Observations")+
      scale_x_continuous()+
      theme_minimal()
  #require(multiplot)
  print(g3)
  
  plots <- list()
  for(xcol_sel in 1:p){
    g4<-ggplot(data=visdf)+
      geom_point(aes(x=as.numeric(as.character(visdf[,xcol_sel]))
                     ,y=y.data,colour=partitions),size=2)+
          ylab("Y") + xlab(paste("X_",xcol_sel))
    plots[[xcol_sel]] <- g4
    print(g4)
  }
  #multiplot(plotlist=plots, cols = 1, main = "Main title")
}


