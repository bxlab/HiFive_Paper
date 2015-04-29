#change matrix into vector
u<-as.matrix(counts)
features = as.matrix(features)
u_vec<-u[upper.tri(u,diag=F)]

#get cov matrix
len_m<-as.matrix(log(features[,3]%o%features[,3]))
gcc_m<-as.matrix(log(features[,4]%o%features[,4]))
map_m<-as.matrix(log(features[,5]%o%features[,5]))

#centralize cov matrix of enz, gcc
len_m<-(len_m-mean(c(len_m)))/sd(c(len_m))
gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m))

#change matrix into vector
len_vec<-len_m[upper.tri(len_m,diag=F)]
gcc_vec<-gcc_m[upper.tri(gcc_m,diag=F)]
map_vec<-map_m[upper.tri(map_m,diag=F)]

#fit Poisson regression: u~len+gcc+offset(map)
fit<-glm(u_vec~len_vec+gcc_vec+offset(map_vec),family="poisson")

#user can use the following two lines to fit negative binomial regression: u~len+gcc+offset(map).
#library("MASS")
#fit<-glm.nb(u_vec~len_vec+gcc_vec+offset(map_vec))

#summary(fit)
coeff<-round(fit$coeff,4)
res<- round(u/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)

#output normalized cis contact map, user can change the name of this output file
write.table(res, file=out_fname, row.names=F, col.names=F, sep="\t", quote=F)
