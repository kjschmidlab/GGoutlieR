# The functions in this script are the internal functions for `ggoutlier_geneticKNN`
# The functions are used to search K nearest neighbors, do KNN prediction and calculate Dgeo values
#---------------------------------
# search KNN
## a function to find KNN
find_gen_knn <- function(pgdM, k){
  indx <- 1:ncol(pgdM)
  res <- matrix(NA, nrow = nrow(pgdM), ncol = k)
  for(i in indx){
    res[i,] <- FastKNN::k.nearest.neighbors(i = i, distance_matrix = pgdM, k = k)
  }
  return(res)
} # find_gen_knn end

## a function to predict with KNN
pred_geo_coord_knn <- function(geo_coord, pgdM, knn.indx, w_power){
  res <- matrix(NA, nrow = nrow(geo_coord), ncol = ncol(geo_coord))
  for(j in 1:nrow(pgdM)){
    tmp.indx <- knn.indx[j,]
    tmp.geo_coord <- geo_coord[tmp.indx,]
    tmp.d <- (pgdM[tmp.indx,j]) ^ w_power
    w <- (1/tmp.d)/(sum(1/tmp.d))
    res[j,] <- apply(tmp.geo_coord, 2, function(x){weighted.mean(x, w)})
  }
  return(res)
} # pred_geo_coord_knn end

## a function to calculate Dg
cal_Dgeo <- function(pred.geo_coord, geo_coord, scalar){
  geodist <- (sapply(1:nrow(geo_coord),function(a){geosphere::distm(x = geo_coord[a,], y = pred.geo_coord[a,])/scalar}))
  return(geodist)
} # cal_Dgeo end
