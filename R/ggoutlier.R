#' Identify outliers with unusual geo-genetic patterns
#' @description This function is used to identify outliers with unusual geo-genetic patterns using the KNN approach.
#' @param geo_coord matrix or data.frame with two columns. The first column is longitude and the second one is latitude.
#' @param gen_coord matrix. A matrix of "coordinates in a genetic space". Users can provide ancestry coefficients or eigenvectors for calculation. If, for example, ancestry coefficients are given, each column corresponds to an ancestral population. Samples are ordered in rows as in `geo_coord`.
#' @param pgdM matrix. A pairwise genetic distance matrix. Users can provide a customized genetic distance matrix with this argument. Samples are ordered in rows and columns as in the rows of `geo_coord`. The default of `pgdM` is `NULL`. If `pgdM` is not provided, a genetic distance matrix will be calculated from `gen_coord`.
#' @param method string. The method to run `GGoutlieR`. It can be "composite", "geneticKNN", or "geoKNN".
#' @param K integer. Number of the nearest neighbors. If `K` is not `NULL`, the value will pass to `k_geneticKNN` and `k_geoKNN`.
#' @param k_geneticKNN integer. Number of the nearest neighbors in a genetic space. The default is `NULL`. The `ggoutlier` will search the optimal K if `k_geneticKNN = NULL`.
#' @param k_geoKNN integer. Number of the nearest neighbors in a geographical space. the default is `NULL`. The `ggoutlier` will search the optimal K if `k_geoKNN = NULL`.
#' @param klim vector. A range of K to search for the optimal number of nearest neighbors. The default is `klim = c(3, 50)`
#' @param make_fig logic. If `make_fig = TRUE`, plots for diagnosing GGoutlieR analysis will be generated and saved to `plot_dir`. The default is `FALSE`
#' @param plot_dir string. The path to save plots
#' @param w_geo numeric. A value controlling the power of distance weight in geographical KNN prediction.
#' @param w_genetic numeric. A value controlling the power of distance weight in genetic KNN prediction.
#' @param p_thres numeric. A significance level
#' @param s integer. A scalar of geographical distance. The default `s=100` scales the distance to a unit of 0.1 kilometer.
#' @param min_nn_dist numeric. A minimal geographical distance for searching KNNs. Neighbors of a focal sample within this distance will be excluded from the KNN searching procedure.
#' @param cpu integer. Number of CPUs to use for searching the optimal K.
#' @param geneticKNN_output output of `ggoutlier_geneticKNN`. Users can use this argument if running `ggoutlier_geneticKNN` in advance.
#' @param geoKNN_output output of `ggoutlier_geoKNN`. Users can use this argument if running `ggoutlier_geoKNN` in advance.
#' @param verbose logic. If `verbose = FALSE`, `ggoutlier` will suppress printout messages.
#' @param multi_stages logic. A multi-stage test will be performed if is `TRUE` (the default is `TRUE`).
#' @param maxIter numeric. Maximal iteration number of multi-stage KNN test.
#' @param keep_all_stg_res logic. Results from all iterations of the multi-stage test will be retained if it is`TRUE`. (the default is `FALSE`)
#' @param warning_minR2 numeric. The prediction accuracy of KNN is evaluated as R^2 to assess the violation of isolation-by-distance expectation. If any R^2 is larger than `warning_minR2`, a warning message will be reported at the end of your analysis.
#' @return a list including five items. `statistics` is a `data.frame` consisting of the D_geo values, p values and a column of logic values showing if a sample is an outlier or not. `threshold` is a `data.frame` recording the significance threshold. `gamma_parameter` is a vector recording the parameter of the heuristic Gamma distribution. `knn_index` and `knn_name` are a `data.frame` recording the K nearest neighbors of each sample.
#' @examples
#' library(GGoutlieR)
#' data("ipk_anc_coef") # get ancestry coefficients
#' data("ipk_geo_coord") # get geographical coordinates
#'
#' ## To reduce computational time, a random subset
#' ## of the IPK barley landrace collection is used
#' ## in this minimal example.
#' ## The argument setting `multi_stages = TRUE`
#' ## and a larger range of `klim` are recommended
#' ## in your analysis.
#'
#' indx <- sample(1:nrow(ipk_geo_coord), size = 100)
#' mini_example <- ggoutlier(geo_coord = ipk_geo_coord[indx,],
#'                              gen_coord = ipk_anc_coef[indx,],
#'                              klim = c(3,6),
#'                              p_thres = 0.01,
#'                              cpu = 2,
#'                              method = "composite",
#'                              verbose = FALSE,
#'                              min_nn_dist = 1000,
#'                              multi_stages = FALSE)
#' summary_ggoutlier(mini_example)
#'
#' # DON'T RUN: this analysis will take a few minutes
#' # full_example <- ggoutlier(geo_coord = ipk_geo_coord,
#' #                           gen_coord = ipk_anc_coef,
#' #                           klim = c(3,6),
#' #                           p_thres = 0.01,
#' #                           cpu = 2,
#' #                           method = "composite",
#' #                           verbose = FALSE,
#' #                           min_nn_dist = 1000,
#' #                           multi_stages = FALSE)
#'
#' @export

ggoutlier <- function(geo_coord,
                      gen_coord,
                      pgdM = NULL,
                      method = c("geneticKNN", "geoKNN", "composite"),
                      K = NULL,
                      k_geneticKNN = NULL,
                      k_geoKNN = NULL,
                      klim = c(3,50),
                      make_fig = FALSE,
                      plot_dir = ".",
                      w_geo = 1,
                      w_genetic = 2,
                      p_thres = 0.05,
                      s = 100,
                      min_nn_dist = 100,
                      cpu = 1,
                      geneticKNN_output = NULL,
                      geoKNN_output = NULL,
                      verbose = TRUE,
                      multi_stages = TRUE,
                      maxIter=NULL,
                      keep_all_stg_res = FALSE,
                      warning_minR2 = 0.9


    ){
  method <- match.arg(method)

  if(!is.null(K)){
      k_geneticKNN <- K
      k_geoKNN <- K
    }
  if(method == "geneticKNN"){
    out <-
      ggoutlier_geneticKNN(geo_coord = geo_coord,
                           gen_coord = gen_coord,
                           pgdM = pgdM,
                           k = k_geneticKNN,
                           klim = klim,
                           plot_dir = plot_dir,
                           w_power = w_genetic,
                           p_thres = p_thres,
                           n = 10^6,
                           s = s,
                           multi_stages = multi_stages,
                           maxIter = maxIter,
                           keep_all_stg_res = F,
                           warning_minR2 = warning_minR2,
                           cpu = cpu,
                           verbose = verbose,
                           make_fig = make_fig
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
                       w_power = w_geo,
                       p_thres = p_thres,
                       n = 10^6,
                       multi_stages = multi_stages,
                       maxIter=maxIter,
                       keep_all_stg_res = keep_all_stg_res,
                       warning_minR2 = warning_minR2,
                       cpu = cpu,
                       verbose = verbose,
                       make_fig = make_fig
      )
  }
  if(method == "composite"){
    if(!multi_stages){message("`composite` method is recommended to be done with the multi-stage test (`multi_stages = TRUE`)")}
    out <-
      ggoutlier_compositeKNN(geo_coord = geo_coord,
                             gen_coord = gen_coord,
                             pgdM = pgdM,
                             k_geneticKNN = k_geneticKNN,
                             k_geoKNN = k_geoKNN,
                             klim = klim,
                             plot_dir = plot_dir,
                             w_geo = w_geo,
                             w_genetic = w_genetic,
                             p_thres = p_thres,
                             n = 10^6,
                             s = s,
                             min_nn_dist = min_nn_dist,
                             multi_stages = multi_stages,
                             maxIter=maxIter,
                             keep_all_stg_res = keep_all_stg_res,
                             warning_minR2 = warning_minR2,
                             cpu = cpu,
                             geneticKNN_output = geneticKNN_output,
                             geoKNN_output = geoKNN_output,
                             verbose = verbose,
                             make_fig = make_fig
      )
  }
  return(out)
} # ggoutlier end



