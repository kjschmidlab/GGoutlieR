#' obtain an adjacency matrix to make a network graph
#' @description `get_GGNet_adjacency_matrix` calculates p-values based on the KNNs and heuristic Gamma distribution obtained in the outlier identification processes. The matrices of p-values are then multiplied with the given genetic similarity matrix to form adjacency matrices.
#' @param GeoSP_knn_res an output from `detect_outlier_in_GeoSpace`
#' @param GenSP_knn_res an output from `detect_outlier_in_GeneticSpace`
#' @param geo_coord a two-column matrix or data.frame. the first column is longitude and the second one is latitude.
#' @param gen_coord a matrix of "coordinates in a genetic space". Users can provide ancestry coefficients or eigenvectors for calculation. If, for example, ancestry coefficients are given, each column corresponds to an ancestral population. Samples are ordered in rows as in `geo_coord`. Users have to provide `pgdM` if `gen_coord` is not given.
#' @param pgsM a matrix of pairwise genetic similarity with n x n dimensions. Weights of edges are extracted from this matrix according to the KNN results. In our demonstration, we used the output of `anc_coeff_to_GeneticSimilarityMatrix`. Users can consider other approaches to characterize genetic similarity.
#' @param mutual logic. If a multi-stage test is used in the outlier identification, some samples could not be a mutual neighbor with its K nearest neighbors. In this case, setting `mutual=TRUE` can force those samples to become mutual neighbors in the output adjacency matrix.
#' @return a list consisting of four matrices that can be used in building network graphs. The default is `TRUE`
#' `GeoSP_adjm` is a matrix describing the strength of edges as genetic similarity times p values from the empirical Gamma distribution identified by `detect_outlier_in_GeoSpace`
#' `GenSP_adjm` is a matrix describing the strength of edges as genetic similarity times p values from the empirical Gamma distribution identified by `detect_outlier_in_GeneticSpace`
#' `GeoSP_pvalue` is a matrix describing the strength of edges as p values from the empirical Gamma distribution identified by `detect_outlier_in_GeoSpace`
#' `GenSP_pvalue` is a matrix describing the strength of edges as p values from the empirical Gamma distribution identified by `detect_outlier_in_GeneticSpace`
#' @examples
#' ## load data
#' GeoSP_knn_res <- readRDS("data/sus_MLE_GeoSP_Gamma_samples.RDS")
#' GenSP_knn_res <- readRDS("data/sus_MLE_GenSP_Gamma_samples.RDS")
#' geo_coord <- read.table("./data/georef_1661ind_geo_coord_for_locator.txt", header = T, stringsAsFactors = F)
#' rownames(geo_coord) <- geo_coord[,1]
#' geo_coord <- geo_coord[,-1]
#' gen_coord <- t(as.matrix(read.csv("./data/alstructure_Q_hat_1661inds.csv", header = F, stringsAsFactors = F)))
#' ## get pairwise genetic similarity matrix
#' pgsM <- anc_coeff_to_GeneticSimilarityMatrix(anc_coef = gen_coord)
#' adjm <- get_GGNet_adjacency_matrix(GeoSP_knn_res = GeoSP_knn_res,
#'                                   GenSP_knn_res = GenSP_knn_res,
#'                                   gen_coord = gen_coord,
#'                                   geo_coord = geo_coord,
#'                                   pgsM = pgsM)
#' @export
# TO DO: make this function to take any ggoulier output
get_GGNet_adjacency_matrix <- function(ggoutlier_res,
    #GeoSP_knn_res,
    #                                   GenSP_knn_res,
                                       geo_coord,
                                       gen_coord,
                                       pgsM = NULL,
                                       mutual = FALSE
){
  if(attributes(ggoutlier_res)$model == "composite"){
    GeoSP_knn_res <- ggoutlier_res$geoKNN_result
    GenSP_knn_res <- ggoutlier_res$geneticKNN_result
    ggoutlier_method <- "composite"
    n = nrow(GeoSP_knn_res$statistics)
  }
  if(attributes(ggoutlier_res)$model == "ggoutlier_geoKNN"){
    GeoSP_knn_res <- ggoutlier_res
    ggoutlier_method <- "ggoutlier_geoKNN"
  }
  if(attributes(geneticknn_res)$model == "ggoutlier_geneticKNN"){
    GenSP_knn_res <- ggoutlier_res
    ggoutlier_method <- "ggoutlier_geneticKNN"
  }


  # check input data
  if(any(names(GeoSP_knn_res) == "combined_result")){
    GeoSP_knn_res <- GeoSP_knn_res$combined_result # extract the combined results (if setting `keep_all_stg_res=T` in previous analyses)
  }
  if(any(names(GenSP_knn_res) == "combined_result")){
    GenSP_knn_res <- GenSP_knn_res$combined_result # extract the combined results (if setting `keep_all_stg_res=T` in previous analyses)
  }
  if(nrow(GeoSP_knn_res$statistics) != nrow(GenSP_knn_res$statistics)){
    stop("the sample sizes of `GeoSP_knn_res` and `GenSP_knn_res` are not equal. please check if correct inputs are provided.")
  }
  if(is.null(pgsM)){
    pgsM <- anc_coeff_to_GeneticSimilarityMatrix(anc_coef = gen_coord)
  }
  pgsM <- as.matrix(unname(pgsM))
  message("obtaining p values of comparisons with k nearest neighbors\n")
  GeoSP.knn.p <- get_knn_pvalue(knn_res = GeoSP_knn_res, gen_coord = gen_coord)
  GenSP.knn.p <- get_knn_pvalue(knn_res = GenSP_knn_res, geo_coord = geo_coord)

  message("pairwise genetic similarity matrix `pgsM` is used as weights of edges\n")
  GeoSP.knn.pgsM <- matrix(NA, ncol = ncol(GeoSP_knn_res$knn_index), nrow = n)
  GenSP.knn.pgsM <- matrix(NA, ncol = ncol(GenSP_knn_res$knn_index), nrow = n)

  # get weights for KNNs
  for(i in 1:n){
    tmp.geosp.indx <- GeoSP_knn_res$knn_index[i,]
    tmp.gensp.indx <- GenSP_knn_res$knn_index[i,]
    GeoSP.knn.pgsM[i,] <- c(pgsM[i,tmp.geosp.indx])
    GenSP.knn.pgsM[i,] <- c(pgsM[i,tmp.gensp.indx])
  }
  # remove edges if shared ancestral allele proportion < 0.8
  #GeoSP.knn.pgsM[GeoSP.knn.pgsM < 0.8] <- 0
  #GenSP.knn.pgsM[GenSP.knn.pgsM < 0.9] <- 0
  # forming a weighted adjacency matrix according to geographical KNN
  GeoSP.adjm <- matrix(0, ncol = n, nrow = n)
  GenSP.adjm <- matrix(0, ncol = n, nrow = n)
  GeoSP.pm <- matrix(0, ncol = n, nrow = n)
  GenSP.pm <- matrix(0, ncol = n, nrow = n)
  if(!is.null(rownames(GeoSP_knn_res$statistics))){
    colnames(GeoSP.adjm) = rownames(GeoSP.adjm) =
      colnames(GenSP.adjm) = rownames(GenSP.adjm) =
      colnames(GeoSP.pm) = rownames(GeoSP.pm) =
      colnames(GenSP.pm) = rownames(GenSP.pm) =
      rownames(GeoSP_knn_res$statistics)
  }

  ## in a adjacency matrix, direction of arrow is from samples on rows to samples on columns
  ## so we fill in values by rows (although I decided to let GGNet as an undirected graph in the end to make plots neater)
  ## the weights are defined as "probablility of a sample not belong to the group of its KNNs" times "genetic similarity"
  ## the way of calculating genetic similarity can be customized by users
  for(i in 1:n){
    # likelihood multiplied with genetic similarity
    GeoSP.adjm[i,GeoSP_knn_res$knn_index[i,]] <- (GeoSP.knn.p[i,]) * GeoSP.knn.pgsM[i,]
    GenSP.adjm[i,GenSP_knn_res$knn_index[i,]] <- (GenSP.knn.p[i,]) * GenSP.knn.pgsM[i,]
    # likelihood only
    GeoSP.pm[i,GeoSP_knn_res$knn_index[i,]] <- GeoSP.knn.p[i,]
    GenSP.pm[i,GenSP_knn_res$knn_index[i,]] <- GenSP.knn.p[i,]
    #GenSP.adjm[i,GenSP_knn_res$knn_index[i,]] <- (GenSP.knn.p[i,])
  }

  ## convert the output to mutual neighborhood matrix
  if(mutual){
    make_mutual_Mat <- function(Mat){
      indx <- Mat != t(Mat)
      #all(apply(cbind(Mat[indx], t(Mat)[indx]), 1, function(x){sum(x != 0)}) == 1)
      Mat[indx] <- apply(cbind(Mat[indx], t(Mat)[indx]), 1, max)
      if(any(Mat != t(Mat))){stop("Fail to convert the given matrix to mutual neighborhood matrix")}
      return(Mat)
    }

    GeoSP.adjm <- make_mutual_Mat(GeoSP.adjm)
    GenSP.adjm <- make_mutual_Mat(GenSP.adjm)
    GeoSP.pm <- make_mutual_Mat(GeoSP.pm)
    GenSP.pm <- make_mutual_Mat(GenSP.pm)
  }

  return(list(GeoSP_adjm = GeoSP.adjm,
              GenSP_adjm = GenSP.adjm,
              GeoSP_pvalue = GeoSP.pm,
              GenSP_pvalue = GenSP.pm))


} # get_GGNet_adjacency_matrix end




#-----------------------------------------------------------------------------
# NOTE: the `get_knn_pvalue` here is the first version.
# `get_knn_pvalue` is used in `get_GGNet_adjacency_matrix`
# obtain p values by computing D statistics between a sample and its KNNs individually with the null Gamma distribution procured in the KNN outlier testing
# status: finished
# arguments:
# knn_res: the output from `detect_outlier_in_GeoSpace` or `detect_outlier_in_GeneticSpace`
# geo_coord: a two column matrix or data.frame. the first column is longitude and the second one is latitude.
# gen_coord: a matrix of "coordinates in a genetic space". Users can provide ancestry coefficients or eigenvectors for calculation. If, for example, ancestry coefficients are given, each column corresponds to an ancestral population. Samples are ordered in rows as in `geo_coord`. Users have provide `pgdM` if `gen_coord` is not given.
# test_type: type of test. It can be `NULL`, `GeoSP` or `GenSP`
get_knn_pvalue <- function(knn_res, geo_coord = NULL, gen_coord = NULL, test_type = NULL){
  library(geosphere)
  # check inputs
  if(any(names(knn_res) == "combined_result")){
    knn_res <- knn_res$combined_result # extract the combined results (if setting `keep_all_stg_res=T` in previous analyses)
  }

  ## identify which type of KNN test results is
  if(is.null(test_type)){
    if("Dg" %in% colnames(knn_res$statistics)){
      test_type = "GeoSP"
    }else{
      if("Dgeo" %in% colnames(knn_res$statistics)){
        test_type = "GenSP"
      }else{
        stop("cannot recognize the data format of `knn_res`. is it the output of `detect_outlier_in_GeoSpace` or `detect_outlier_in_GeneticSpace`?\n")
      }
    }
  }

  if(test_type == "GeoSP"){
    if(is.null(gen_coord)){stop("`get_knn_pvalue` requires `gen_coord` since you provide the output from `detect_outlier_in_GeoSpace`")}
    if(nrow(gen_coord) != nrow(knn_res$knn_index)){stop("the sample sizes of `knn_res` and `gen_coord` are unequal. please check if you provide correct inputs.")}
  }
  if(test_type == "GenSP"){
    if(is.null(geo_coord)){stop("`get_knn_pvalue` requires `geo_coord` since you provide the output from `detect_outlier_in_GeneticSpace`")}
    if(nrow(geo_coord) != nrow(knn_res$knn_index)){stop("the sample sizes of `knn_res` and `geo_coord` are unequal. please check if you provide correct inputs.")}
    if(ncol(geo_coord) != 2){stop("please ensure you provide correct format of geographical coordinates: a dataframe with longitude in the 1st column and latitude in the 2nd column.")}
  }

  # get the parameters of Gamma distribution
  a <- unname(knn_res$gamma_parameter[1])
  b <- unname(knn_res$gamma_parameter[2])

  # calculating p values
  n = nrow(knn_res$knn_index) # get sample size
  k = ncol(knn_res$knn_index) # get the number of neighbors
  knn.p <- matrix(NA, ncol = k, nrow = n)
  rownames(knn.p) <- rownames(knn_res$statistics)
  if(test_type == "GeoSP"){
    for(i in 1:n){
      tmp.gen_coord <- unname(gen_coord[i,])
      knn.gen_coord <- gen_coord[knn_res$knn_index[i,],]
      knn.Dg <- apply(knn.gen_coord, 1, function(x){
        return(sum((x - tmp.gen_coord)^2))
      })
      knn.p[i,] <- 1 - pgamma(unname(knn.Dg), shape = a, rate = b)
    }
  }else{
    s <- knn_res$scalar
    for(i in 1:n){
      tmp.geo_coord <- unname(geo_coord[i,])
      knn.geo_coord <- geo_coord[knn_res$knn_index[i,],]
      knn.Dgeo <- apply(knn.geo_coord, 1, function(a){
        return(distm(x = tmp.geo_coord, y = a)/s)
      })
      knn.p[i,] <- 1 - pgamma(unname(knn.Dgeo), shape = a, rate = b)
    }
  }

  return(knn.p)
  #plot(-log10(apply(knn.p,1,mean)), -log10(knn_res$statistics$p.value),
  #     xlim = c(0,10), ylim=c(0,10))
  #cor((apply(knn.p,1,mean)), (knn_res$statistics$p.value))
} # get_knn_pvalue end


#-----------------------------------------------------------------------------
# `get_knn_pvalue_v2`
# amplify p values by computing D statistics between a sample and its KNNs individually with the null Gamma distribution procured in the KNN outlier testing
# status: finished
# arguments:
# knn_res: the output from `detect_outlier_in_GeoSpace` or `detect_outlier_in_GeneticSpace`
# geo_coord: a two column matrix or data.frame. the first column is longitude and the second one is latitude.
# gen_coord: a matrix of "coordinates in a genetic space". Users can provide ancestry coefficients or eigenvectors for calculation. If, for example, ancestry coefficients are given, each column corresponds to an ancestral population. Samples are ordered in rows as in `geo_coord`. Users have provide `pgdM` if `gen_coord` is not given.
# test_type: type of test. It can be `NULL`, `GeoSP` or `GenSP`
get_knn_pvalue_v2 <- function(knn_res, geo_coord = NULL, gen_coord = NULL, test_type = NULL){
  library(geosphere)
  # check inputs
  if(any(names(knn_res) == "combined_result")){
    knn_res <- knn_res$combined_result # extract the combined results (if setting `keep_all_stg_res=T` in previous analyses)
  }

  ## identify which type of KNN test results is
  if(is.null(test_type)){
    if("Dg" %in% colnames(knn_res$statistics)){
      test_type = "GeoSP"
    }else{
      if("Dgeo" %in% colnames(knn_res$statistics)){
        test_type = "GenSP"
      }else{
        stop("cannot recognize the data format of `knn_res`. is it the output of `detect_outlier_in_GeoSpace` or `detect_outlier_in_GeneticSpace`?\n")
      }
    }
  }

  if(test_type == "GeoSP"){
    if(is.null(gen_coord)){stop("`get_knn_pvalue` requires `gen_coord` since you provide the output from `detect_outlier_in_GeoSpace`")}
    if(nrow(gen_coord) != nrow(knn_res$knn_index)){stop("the sample sizes of `knn_res` and `gen_coord` are unequal. please check if you provide correct inputs.")}
  }
  if(test_type == "GenSP"){
    if(is.null(geo_coord)){stop("`get_knn_pvalue` requires `geo_coord` since you provide the output from `detect_outlier_in_GeneticSpace`")}
    if(nrow(geo_coord) != nrow(knn_res$knn_index)){stop("the sample sizes of `knn_res` and `geo_coord` are unequal. please check if you provide correct inputs.")}
    if(ncol(geo_coord) != 2){stop("please ensure you provide correct format of geographical coordinates: a dataframe with longitude in the 1st column and latitude in the 2nd column.")}
  }

  # duplicate p values to make a matrix
  k <- ncol(knn_res$knn_index)
  knn.p <- matrix(rep(knn_res$statistics$p.value, times = k), nrow = nrow(knn_res$statistics), ncol = k)

  return(knn.p)
} # get_knn_pvalue_v2 end
