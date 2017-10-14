#This shows how we could update the Ising model using Gibbs sampler.
size=6;
dfg<-as.data.frame(
  matrix(NA,nrow=size^2,ncol=4)
)
colnames(dfg)<-c('vertex','x','y','obs');

vertex_prob   <-matrix(runif(size^2),ncol=1,nrow=size^2);
edge_interact <-matrix(runif(size^4),ncol=size^2,nrow=size^2);
dfg$vertex<-seq(1:(size^2));
dfg$obs<-sample(c(0,1),size=36,replace=T);
#initial observations
ct<-1;
for(i in 1:size){
  for(j in 1:size){
    dfg[ct,]$x<-i;
    dfg[ct,]$y<-j;
    ct<-ct+1;
  }
}
#Ising plot
Isingplot<-function(df=dfg,s=size,draw=FALSE){
  if(draw==TRUE){plot(0,0,xlim=c(0,size),ylim=c(0,size));}
  mat<-matrix(NA,ncol=size,nrow=size);
  
  for(s in 1:dim(df)[1]){
    if(draw==TRUE){points(x=df[s,]$x,y=df[s,]$y,pch=df[s,]$obs*15);}
    mat[df[s,]$x,df[s,]$y]<-df[s,]$obs;
  }
  
  print(mat)
  
}
#Gibbs sampler below
M<-100;
for(k in 1:M){
  choose_vertex<-sample(dfg$vertex,size=1,replace=F);
  u<-runif(1);
  s<-choose_vertex;
  Isingplot(dfg,size)
    ver_x<-dfg[s,]$x;
    ver_y<-dfg[s,]$y;
    nbd1<-dfg[dfg$x==ver_x-1 & dfg$y==ver_y,]
    nbd2<-dfg[dfg$x==ver_x+1 & dfg$y==ver_y,]
    nbd3<-dfg[dfg$x==ver_x & dfg$y==ver_y+1,]
    nbd4<-dfg[dfg$x==ver_x & dfg$y==ver_y-1,]
    summ=0;
    if(dim(nbd1)[1]>0){summ<-summ+edge_interact[nbd1$vertex,dfg[s,]$vertex]*dfg[s,]$obs};
    if(dim(nbd2)[1]>0){summ<-summ+edge_interact[nbd2$vertex,dfg[s,]$vertex]*dfg[s,]$obs};
    if(dim(nbd3)[1]>0){summ<-summ+edge_interact[nbd3$vertex,dfg[s,]$vertex]*dfg[s,]$obs};
    if(dim(nbd4)[1]>0){summ<-summ+edge_interact[nbd4$vertex,dfg[s,]$vertex]*dfg[s,]$obs};
    
    if(u < 1/(1 + exp(-(vertex_prob[s]+summ) ) ) ){
      dfg[s,]$obs<-1;
    }else{
      dfg[s,]$obs<-0;
    }
  
}
