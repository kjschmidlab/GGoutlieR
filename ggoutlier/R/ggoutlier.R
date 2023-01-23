#' Identify samples geographically remote from K genetically nearest neighbors
#' @param geo_coord a two column matrix or data.frame. the first column is longitude and the second one is latitude.
#' @param gen_coord a matrix of "coordinates in a genetic space". Users can provide ancestry coefficients or eigenvectors for calculation. If, for example, ancestry coefficients are given, each column corresponds to an ancestral population. Samples are ordered in rows as in `geo_coord`. Users have to provide `pgdM` if `gen_coord` is not given.
#' @param pgdM a pairwise genetic distance matrix. Users can provide a customized genetic distance matrix with this argument. Samples are ordered in rows and columns as in the rows of `geo_coord`. The default of `pgdM` is `NULL`. If `pgdM` is not provided, a genetic distance matrix will be calculated from `gen_coord`.
#' @param k_geneticKNN number of the nearest neighbor in a genetic space. the default is `NULL`.
#' @param k_geoKNN number of the nearest neighbor in a geographical space. the default is `NULL`.
#' @param klim if `k = NULL`, an optimal k will be searched between the first and second value of `klim`
#' @param s a scalar of geographical distance. The default `s=100` scales the distance to a unit of 1 kilometer.
#' @param plot_dir the path to save plots
#' @param w_power a value controlling the power of distance weight in KNN prediction. For example, if `w_power=2`, the weight of KNN is 1/d^2/sum(1/d^2).
#' @param p_thres a significance level
#' @param n number of samples to draw from the null distribution (to obtain the range of x axis to make a curve plot of null distribution)
#' @param multi_stages logic. a multi-stage test will be performed if is `TRUE` (the default is `TRUE`).
#' @param maxIter maximal iteration number of multi-stage KNN test.
#' @param keep_all_stg_res logic. results from all iterations of the multi-stage test will be retained if it is`TRUE`. (the default is `FALSE`)
#' @param warning_minR2 the prediction accuracy of KNN is evaluated as R^2 to assess the violation of isolation-by-distance expectation. if any R^2 is larger than `warning_minR2`, a warning message will be reported at the end of your analysis.
#' @return a list including five items. `statistics` is a `data.frame` consisting of the D_geo values, p values and a column of logic values showing if a sample is an outlier or not. `threshold` is a `data.frame` recording the significance threshold. `gamma_parameter` is a vector recording the parameter of the heuristic Gamma distribution. `knn_index` and `knn_name` are a `data.frame` recording the K nearest neighbors of each sample.
#'
#'
#'
#' @examples
#' @export

ggoutlier <- function(geo_coord,
                      gen_coord,
                      pgdM = NULL,
                      k_geneticKNN = NULL,
                      k_geoKNN = NULL,
                      klim = c(3,50),
                      plot_dir = ".",
                      w_power = 2,
                      p_thres = 0.05,
                      n = 10^6,
                      s = 100,
                      min_nn_dist = NULL,
                      multi_stages = T,
                      maxIter=NULL,
                      keep_all_stg_res = F,
                      warning_minR2 = 0.9,
                      cpu = 1,
                      geneticKNN_output = NULL,
                      geoKNN_output = NULL,
                      method = c("geneticKNN", "geoKNN", "composite"),
                      verbose = TRUE

    ){
  method <- match.arg(method)
  if(method == "geneticKNN"){
    out <-
      ggoutlier_geneticKNN(geo_coord = geo_coord,
                           gen_coord = gen_coord,
                           pgdM = pgdM,
                           k = k_geneticKNN,
                           klim = klim,
                           plot_dir = plot_dir,
                           w_power = w_power,
                           p_thres = p_thres,
                           n = n,
                           s = s,
                           multi_stages = multi_stages,
                           maxIter = maxIter,
                           keep_all_stg_res = F,
                           warning_minR2 = warning_minR2,
                           cpu = cpu,
                           verbose = verbose
      )
  }
  if(method == "geoKNN"){
    out <-
      ggoutlier_geoKNN(geo_coord = geo_coord,
                       gen_coord = gen_coord,
                       min_nn_dist = min_nn_dist,
                       k = k_geoKNN,
                       klim = klim,
                       s = s,
                       plot_dir = plot_dir,
                       w_power = w_power,
                       p_thres = p_thres,
                       n = n,
                       multi_stages = multi_stages,
                       maxIter=maxIter,
                       keep_all_stg_res = keep_all_stg_res,
                       warning_minR2 = warning_minR2,
                       cpu = cpu,
                       verbose = verbose
      )
  }
  if(method == "composite"){
    if(!multi_stages){warning("`composite` method is recommended to be done with the multi-stage test (`multi_stages = TRUE`)")}
    out <-
      ggoutlier_compositeKNN(geo_coord = geo_coord,
                             gen_coord = gen_coord,
                             pgdM = pgdM,
                             k_geneticKNN = k_geneticKNN,
                             k_geoKNN = k_geoKNN,
                             klim = klim,
                             plot_dir = plot_dir,
                             w_power = w_power,
                             p_thres = p_thres,
                             n = n,
                             s = s,
                             min_nn_dist = min_nn_dist,
                             multi_stages = multi_stages,
                             maxIter=maxIter,
                             keep_all_stg_res = keep_all_stg_res,
                             warning_minR2 = warning_minR2,
                             cpu = cpu,
                             geneticKNN_output = geneticKNN_output,
                             geoKNN_output = geoKNN_output,
                             verbose = verbose
      )
  }
  return(out)
} # ggoutlier end



