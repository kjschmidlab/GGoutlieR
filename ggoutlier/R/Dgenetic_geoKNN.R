# The functions in this script are the internal functions for `ggoutlier_geoKNN`
# The functions are used to search K nearest neighbors, do KNN prediction and calculate Dg values
#-------------------------------------------
# search KNN
## a function to find KNN
find_geo_knn <- function(geo.dM, k, min_nn_dist){
  indx <- 1:ncol(geo.dM)
  res <- matrix(NA, nrow = nrow(geo.dM), ncol = k)
  for(i in indx){
    if(any(is.null(min_nn_dist), min_nn_dist == 0)){
      res[i,] <- FastKNN::k.nearest.neighbors(i = i, distance_matrix = geo.dM, k = k)
    }else{
      tmp.d <- cbind(indx, geo.dM[,i])
      tmp.d <- tmp.d[order(tmp.d[,2], decreasing = F),]
      # get the k nearest neighbors but not within a distance of `min_nn_dist`
      res[i,] <- tmp.d[tmp.d[,2] > min_nn_dist,][1:k,1]
    }
  }
  return(res)
} # find_geo_knn end

## a function to predict with KNN
pred_q_knn <- function(geo_coord, gen_coord, geo.dM, knn.indx, w_power){
  res <- matrix(NA, nrow = nrow(geo_coord), ncol = ncol(gen_coord))
  for(j in 1:nrow(gen_coord)){
    tmp.indx <- knn.indx[j,]
    tmp.q <- gen_coord[tmp.indx,]
    tmp.d <- geo.dM[tmp.indx,j] ^ w_power
    w <- (1/tmp.d)/(sum(1/tmp.d))
    res[j,] <- apply(tmp.q, 2, function(x){weighted.mean(x, w)})
  }
  return(res)
} # pred_q_knn end

## a function to calculate Dg
cal_Dg <- function(pred.q, gen_coord){
  apply(pred.q - gen_coord, 1, function(x){
    mean(x^2)
  })
} # cal_Dg end
