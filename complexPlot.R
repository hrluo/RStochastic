#This function is used to show how you can illustrate the complex transformation on a complex plane in R.
globalTop<-2*pi;
globalBottom<--2*pi;
globalLeft<--2*pi;
globalRight<-2*pi;
gridscale<-.1;
ftop<-function(x,y){
  if(y<2){
    return(1)
    #This is the top boundary curve in the complex surface.
  }else{
    return(0)
  }
}
fbottom<-function(x,y){
  if(y>1){
    return(1)
    #This is the bottom boundary curve in the complex surface.
  }else{
    return(0)
  }
}
fleft<-function(x,y){
  if(1){
    return(1)
    #This is the left boundary curve in the complex surface.
  }else{
    return(0)
  }
}
fright<-function(x,y){
  if(1){
    return(1)
    #This is the right boundary curve in the complex surface.
  }else{
    return(0)
  }
}
#Following loop form a grid that is decided by the four boundary functions.
xgrid<-seq(globalLeft,globalRight,gridscale)
ygrid<-seq(globalBottom,globalTop,gridscale)
testgrid<-expand.grid(xgrid,ygrid)
plot(0,0,xlim=c(globalLeft,globalRight),ylim=c(globalBottom,globalTop))
for(k in 1:dim(testgrid)[1]){
  testx<-testgrid[k,1];
  testy<-testgrid[k,2];
  if(ftop(testx,testy)*
     fbottom(testx,testy)*
     fleft(testx,testy)*
     fright(testx,testy)!=1){
    testgrid[k,]<-c(NA,NA)
  }else{
    points(testx,testy,xlim=c(globalLeft,globalRight),ylim=c(globalBottom,globalTop),cex=.5,pch=20)
  }
}
#Define the conformal map we are to investigate.
conformalMap1<-function(x,y){
  #For example we need to write the complex number in form of matrix and perform operations thereafter.
  matForm<-x*diag(2)+y*matrix(c(0,-1,1,0),nrow=2);
  #Here we construct the exponential mapping.
  library(expm)
  res<-diag(2)*0;
  for(k in 0:10){
    res<-res+matForm%^%k/factorial(k)
  }
  return( c(res[1,1],res[2,1]) )
}

conformalMap2<-function(x,y){
  #For example we need to write the complex number in form of matrix and perform operations thereafter.
  matForm<-x*diag(2)+y*matrix(c(0,-1,1,0),nrow=2);
  #Here we construct the exponential mapping.
  library(expm)
  matLinearTransform<-matrix(c(2,1,3,4));
  #This matrix represents the complex linear transform;
  res<-matLinearTransform%*%matForm;
  return( c(res[1,1],res[2,1]) )
}

#Following is a function that applies the complex mapping to the domain grid.
applyMap<-function(gridlist,bivariatefun){
  gridlistNew<-gridlist;
  for(k in 1:dim(gridlist)[1]){
    rec<-as.numeric(gridlist[k,]);
    applied<-bivariatefun(rec[1],rec[2]);
    gridlistNew[k,]<-applied;
  }
  return(gridlistNew)
}

transformedgrid<-applyMap(testgrid,conformalMap1);

#plot(0,0,xlim=c(-5,5),ylim=c(-5,5))
for(k in 1:dim(transformedgrid)[1]){
  newx<-as.numeric(transformedgrid[k,1]);
  newy<-as.numeric(transformedgrid[k,2]);
  points(newx,newy,xlim=c(globalLeft,globalRight),ylim=c(globalBottom,globalTop),cex=.5,pch=20,col="red");
}
