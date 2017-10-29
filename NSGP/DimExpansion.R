#The data set I used is the same as
#Bornn, Luke, Gavin Shaddick, and James V. Zidek. "Modeling nonstationary processes through dimension expansion." Journal of the American Statistical Association 107.497 (2012): 281-289.
#Section 3.1(Solar radiation data set)
#Hay, John E. "An assessment of the mesoscale variability of solar radiation at the earth's surface." Solar Energy 32.3 (1984): 425-434.
#According to a personal email correspondence with Mr.Bornn, he directed me to p.115 of
#Sampson, P. D., & Guttorp, P. (1992). Nonparametric estimation of nonstationary spatial covariance structure. Journal of the American Statistical Association, 87(417), 108-119.
#and said they never used the full solar data set from Hay's paper.
#I reproduce this solar radiation data set from [Sampson&Guttorp]p.115
#dispersion matrix for solar stations
A<-NA;
solar_dispersion<-as.matrix(rbind(
  c( A,34,36,50,41,50,63,38,42,47,49,49),
  c(34, A, 5,14, 9,15,26,16,18,31,31,29),
  c(36, 5, A,12, 6,12,25,14,17,31,31,29),
  c(50,14,12, A, 7, 4,12,22,21,41,36,32),
  c(41, 9, 6, 7, A, 6,17,14,15,33,30,27),
  c(50,15,12, 4, 6, A,10,19,17,37,31,29),
  c(63,26,25,12,17,10, A,26,24,40,34,31),
  c(38,16,14,22,14,19,26, A, 6,14,14,14),
  c(42,18,17,21,15,17,24, 6, A,15,11,10),
  c(47,31,31,41,33,37,40,14,15, A, 6, 9),
  c(49,31,31,36,30,31,34,14,11, 6, A, 3),
  c(49,29,29,32,27,29,31,14,10, 9, 3, A)
)
)
colnames(solar_dispersion)<-c(1,2,3,12,11,4,5,6,10,7,8,9)
rownames(solar_dispersion)<-c(1,2,3,12,11,4,5,6,10,7,8,9)
solar_dispersion<-(solar_dispersion/100)
#Since the term in table1 in [Sampson&Guttorp] is 100*(d_ij^2).

v_star<-function(i,j){
  return(solar_dispersion[as.character(i),as.character(j)])
}
#Now we consider the real dataset from [Hay]1984, Table 2(p.429)
locations2D<-as.data.frame(rbind(
  c(1,"UBC CLI. STN",49.16,123.15,93),
  c(2,"Langara",49.13,123.06,63),
  c(3,"Vancouver Airport",49.11,123.10,5),
  c(4,"B.C. Hydro",49.16,123.07,122),
  c(5,"Tsawwassen Ferry",49.00,123.08,3),
  c(6,"North Mount",49.19,123.04,114),
  c(7,"Grouse Mountain",49.23,123.05,1128),
  c(8,"Pitt Meadows Airport",49.13,122.42,5),
  c(9,"Mission - Habitat Apts",49.08,122.17,125),
  c(10,"Abbotsford (Library)",49.02,122.17,60),
  c(11,"Abbotsford Airport",49.01,122.22,61),
  c(12,"Langley",49.06,122.38,11) 
)
)
colnames(locations2D)<-c("num","station","lat","long","elev")


x_coord<-as.numeric(as.character(locations2D$lat))
y_coord<-as.numeric(as.character(locations2D$long))
elev_real<-as.numeric(as.character(locations2D$elev))
Z<-c()
library(rgdal)
#Use rgdal to convert the longitude-latitude into Universal Transverse Mercator coordinate system.
LatLong <- data.frame(X =x_coord, Y = -y_coord)
names(LatLong) <- c("Latitude","Longitude")
coordinates(LatLong) <- ~ Longitude+Latitude  # longitude first
# Add a coordinate reference system, here I used the default choice of WGS1984.
proj4string(LatLong) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
Utm <- spTransform(LatLong, CRS("+proj=utm +zone=10 ellps=WGS84"))
#The data set is collected in UTM zone 10.
UTM_X<-as.vector(Utm$Longitude/1000);
UTM_Y<-as.vector(Utm$Latitude/1000);
#Make the coordinate unit into km as [Sampson&Guttorp]

#gamma function
gamma<-function(d,phi1,phi2,phi3){
  res<-0;
  res<-phi1*( 1-exp(-d/phi2) )+phi3
  return(res);
}
#distance function
d_ij<-function(vec1,vec2){
  return( sqrt(sum((vec1-vec2)^2)) )
}
#objective function
OBJ<-function(phi1,phi2,phi3,Z1=rep(0,12),Z2=rep(0,12),lambda1){
  res<-0;
  #X is the known coordinates of locations in lower dimension from [Hay]1984.
  augmented_dist<-matrix(NA,nrow=12,ncol=12);
  for(a in 1:12){
    for(b in 1:12){
      LOC1<-c(UTM_X[a],UTM_Y[a],Z1[a],Z2[a]);
      LOC2<-c(UTM_X[b],UTM_Y[b],Z1[b],Z2[b]);
      augmented_dist[a,b]<-d_ij(LOC1,LOC2);
    }
  }
  #print(augmented_dist)
  for(j in 2:12){
    #there are 12 observations.
    for(i in 1:(j-1) ){
      #message(i,j)
      temp_dist_ij<-augmented_dist[i,j]
      res<-res+( v_star(i,j)-gamma(d=temp_dist_ij,phi1=phi1,phi2=phi2,phi3=phi3) )^2
    }
  }
  #Penalize by L1 norm of the expanded dimension.
  res<-res+lambda1*sum(abs(Z1))+lambda1*sum(abs(Z2))
  return(res)
}
OBJ05<-function(vec){
  #This is a interface function for optim procedure.
  #We must write such an interface because optim procedure only optimize one vector parameter at one run.
  #vec[1]~vec[3] correspond to phi1~phi3
  #vec[4]~vec[15] correspond to Z, which is a 12-dim vector corresponds to elevations of 12 locations.
  return( OBJ(phi1=vec[1],phi2=vec[2],phi3=vec[3],Z1=vec[4:15],Z2=vec[16:27],lambda1=0.5) )
  #lambda=0.5 is chosen from [Bornn et.al]2012
}

EXPAND<-function(){
  #Specify a location using observed data points in lower dimension.
  OPT<-optim(par=c(1,2,0,rep(0,12),rep(0,12)),fn=OBJ05,method="Nelder-Mead")
  #phi1 and phi2 initial value is chosen to be 1,2 because they are multipliers in covariance.
  #phi3 initial value is chosen to be 0 because it is additive constant in covariance.
  #Z initial value is chosen to be 0 because it is elevation.
  #According to [Bornn]2012, Nelder-Mead and BFGS works well for small amount of observations.
  #print(OPT)
  #optimiser phi
  phi_hat<-OPT$par[1:3]
  #optimiser Z
  Z1_hat<-OPT$par[4:15]
  Z2_hat<-OPT$par[16:27]
  #optimiser augmented [X,Z]
  result<-list(phi=phi_hat,Z1=Z1_hat,Z2=Z2_hat)
  return(result)
}
EXPAND()
Z1<-EXPAND()$Z1*1000;
Z2<-EXPAND()$Z2*1000;
library(fields)
lambda2<-10^-4;
elev_fit1<- Tps(cbind(UTM_X,UTM_Y), Z1,lambda=lambda2);
#The roughness penalty is built into the Tps function.
#surface(elev_fit)
out.p1<- predictSurface(elev_fit1,extrap=T)
surface(out.p1, type="C") 
points(UTM_X,UTM_Y)

elev_fit2<- Tps(cbind(UTM_X,UTM_Y), Z2,lambda=lambda2);
out.p2<- predictSurface(elev_fit2,extrap=T)
surface(out.p2, type="C") 
points(UTM_X,UTM_Y)
