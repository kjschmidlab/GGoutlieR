#' Identify samples geographically remote from K genetically nearest neighbors
#' @param geo_coord a two column matrix or data.frame. the first column is longitude and the second one is latitude.
#' @param gen_coord a matrix of "coordinates in a genetic space". Users can provide ancestry coefficients or eigenvectors for calculation. If, for example, ancestry coefficients are given, each column corresponds to an ancestral population. Samples are ordered in rows as in `geo_coord`. Users have to provide `pgdM` if `gen_coord` is not given.
#' @param pgdM a pairwise genetic distance matrix. Users can provide a customized genetic distance matrix with this argument. Samples are ordered in rows and columns as in the rows of `geo_coord`. The default of `pgdM` is `NULL`. If `pgdM` is not provided, a genetic distance matrix will be calculated from `gen_coord`.
#' @param k_geneticKNN number of the nearest neighbor in a genetic space. the default is `NULL`.
#' @param k_geo number of the nearest neighbor in a geographical space. the default is `NULL`.
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
# STATUS: NOT FINISHED
ggoutlier_compositeKNN <- function(geo_coord,
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
                                   geoKNN_output = NULL
                                    ){
  require(geosphere) # for calculating geographical distances
  require(stats4) # package to perform maximum likelihood estimation
  require(FastKNN) # search KNN with an arbitrary distance matrix
  if(cpu > 1){
    require(foreach)
    require(doParallel)
    max_cores=detectCores()
    if(cpu >= max_cores){
        warning(paste0("\n The maximum number of CPUs is ", max_cores, ". Set `cpu` to ", max_cores-1," \n"))
        cpu <- max_cores -1 # reserve max-1 cpus if users request too many cpus
    }
    if(cpu > 1){
      cat(paste0("\n Parallelize computation using ", cpu, " cores \n"))
      cl <- makeCluster(cpu)
      do_par <- TRUE
    }
  } else {
    if(cpu < 1){stop("`cpu` has to be at least 1")}
    do_par <- FALSE
    cat("\n Computation using 1 core. You can parallelize computation by setting a higher value to `cpu` argument \n")
  }

  # check inputs
  if(all(is.null(gen_coord))){
    stop("`gen_coord` has to be provided!")
  }
  if(all(!is.null(gen_coord) & !is.null(pgdM))){
    warning("Both `gen_coord` and `pgdM` are provided. `pgdM` will be used instead of calculating genetic distances from `gen_coord`.\n")
  }
  if(ncol(geo_coord) != 2){stop("Please ensure the `geo_coord` having two columns!")}
  if(!is.null(gen_coord)){
    if(nrow(geo_coord) != nrow(gen_coord)){stop("`geo_coord` and `gen_coord` has different sample size!")}
  }else{
    if(!is.null(pgdM)){
      if(nrow(geo_coord) != nrow(pgdM) | nrow(geo_coord) != ncol(pgdM)){stop("`geo_coord` and `pgdM` has different sample size!")}
    }
  }

  if(is.null(rownames(geo_coord))){
    rownames(geo_coord) <- paste("sample", 1:nrow(geo_coord), sep = "")
  }else{
    if(any(is.na(rownames(geo_coord)))){
      cat("Some samples in `geo_coord` have IDs but some do not. Arbitratry IDs are assigned to all samples...")
      rownames(geo_coord) <- paste("sample", 1:nrow(geo_coord), sep = "")
    }
  }
  # check missing values
  if(any(apply(geo_coord, 1, function(x){any(is.na(x))}))){
    n_geomis <- sum(apply(geo_coord, 1, function(x){any(is.na(x))}))
    warning(paste0(n_geomis," samples are missing in your `geo_coord`, they are removed from the analysis"))
  }
  if(any(apply(gen_coord, 1, function(x){any(is.na(x))}))){
    n_genmis <- sum(apply(geo_coord, 1, function(x){any(is.na(x))}))
    warning(paste0(n_genmis," samples are missing in your `gen_coord`, they are removed from the analysis"))
  }
  to.rm <- apply(geo_coord, 1, function(x){any(is.na(x))}) | apply(gen_coord, 1, function(x){any(is.na(x))})
  if(any(to.rm)){
    geo_coord <- geo_coord[!to.rm,]
    gen_coord <- gen_coord[!to.rm,]
  }

  #--------------------------------
  #--------------------------------
  # Data pre-treatment
  #--------------------------------
  #--------------------------------
  ## calculate geographical distance
  geo.dM <- distm(x = geo_coord)/s
  ## handle samples with identical geographical coordinates
  if(any(geo.dM[lower.tri(geo.dM)] == 0)){
    cat("Find samples with identical geographical coordinates.\n")
    if(any(min_nn_dist == 0 , is.null(min_nn_dist))){
      cat("Add one unit of distance to individual pairs\n")
      geo.dM <- geo.dM + 1
      diag(geo.dM) <- 0
    }else{
      cat(paste0("\n\nIgnore neighbors within ",min_nn_dist, " unit(s) of distance (the default unit of distance is km)\n"))
    }
  }

  # calculate genetic distance using ancestry coefficients
  if(is.null(pgdM)){
    pgdM <- as.matrix(dist(gen_coord, method = "euclidean"))
  }
  # handle samples with identical genotypes
  if(any(pgdM[lower.tri(pgdM)] == 0)){
    cat("Find samples with zero genetic distances.\n")
    cat("Add 10^-6 unit of distance to individual pairs\n")
    pgdM <- pgdM + 10^-6
    diag(pgdM) <- 0
  }

  ## collect data if `ggoutlier_geneticKNN` and `ggoutlier_geoKNN` are run beforehand

  ## parse `ggoutlier_geneticKNN` output
  if(!is.null(geneticKNN_output)){
    if(attributes(geneticKNN_output, "model") != "ggoutlier_geneticKNN"){
      stop("The input for `geneticKNN_output` is incorrect. Please check if you give an output of `geneticKNN_output` or not.")
    }
    if(nrow(geneticKNN_output$knn_index) != nrow(geo_coord)){
      stop(paste0("The sample size of `geneticKNN_output` is ", nrow(geneticKNN_output$knn_index), " It does not match the sample size of the other inputs"))
    }
    if(is.null(k_geneticKNN)){
      k_geneticKNN <- ncol(geneticKNN_output$knn_index) # get optimal K
    }
    gamma_par_geneticKNN <- geneticKNN_output$gamma_parameter # get parameters of gamma distribution
  }else{gamma_par_geneticKNN <- NULL}

  ## parse `ggoutlier_geoKNN` output
  if(!is.null(geoKNN_output)){
    if(attributes(geoKNN_output, "model") != "ggoutlier_geoKNN"){
      stop("The input for `geoKNN_output` is incorrect. Please check if you give an output of `ggoutlier_geoKNN` or not.")
    }
    if(nrow(geneticKNN_output$knn_index) != nrow(geo_coord)){
      stop(paste0("The sample size of `geoKNN_output` is ", nrow(geoKNN_output$knn_index), " It does not match the sample size of the other inputs"))
    }
    if(is.null(k_geoKNN)){
      k_geoKNN <- ncol(geoKNN_output$knn_index) # get optimal K
    }
    gamma_par_geoKNN <- geoKNN_output$gamma_parameter # get parameters of gamma distribution
  }else{gamma_par_geoKNN <- NULL}

  ##---------------------------------------------------
  ##---------------------------------------------------
  ## KNN in genetic space: find the optimal K
  ##---------------------------------------------------
  ##---------------------------------------------------
  ## NOTE: this step would be skipped if `geneticKNN_output` is available
  if(is.null(k_geneticKNN)){
    cat(paste("\n `k_geneticKNN` is NULL; searching for optimal k between", klim[1], "and", klim[2],"\nthis process can take time...\n"))
    all.D = find_optimalK_geneticKNN(geo_coord = geo_coord,
                                     pgdM = pgdM,
                                     w_power = w_power,
                                     klim = klim,
                                     do_par = do_par,
                                     s = s,
                                     cl = cl)
    opt.k = c(klim[1]:klim[2])[which.min(all.D)]
    k.sel.plot <- paste(plot_dir, "/KNN_Dgeo_optimal_k_selection.pdf", sep = "")
    pdf(k.sel.plot, width = 5, height = 4)
    par(mar=c(4,6,1,1))
    plot(x = klim[1]:klim[2], y = all.D, xlab="K", ylab=expression(sum(D["geo,i"], i==1, n)))
    abline(v = opt.k)
    legend("top",legend = paste("optimal k =", opt.k), pch="", bty = "n",cex = 1.2)
    invisible(dev.off())

    k_geneticKNN = opt.k
    cat(paste("\n The optimal k_geneticKNN is ",opt.k,". Its figure is saved at ", k.sel.plot," \n", sep = ""))

  }
  if(is.null(k_geneticKNN)){stop("k_geneticKNN is NULL!")}

  ##---------------------------------------------------
  ##---------------------------------------------------
  ## KNN in genetic space: prediction with the optimal K (or K given by users)
  ##---------------------------------------------------
  ##---------------------------------------------------
  cat("\nSearching K nearest neighbors (in a genetic space)...\n")
  knn.indx_geneticKNN <- find_gen_knn(pgdM, k=k_geneticKNN)
  pred.geo_coord <- pred_geo_coord_knn(geo_coord = geo_coord,
                                       pgdM = pgdM,
                                       knn.indx = knn.indx_geneticKNN,
                                       w_power = w_power)
  # calculate Dgeo statistic
  cat(paste("\n\n D geo is scaled to a unit of",s,"meters \n"))
  Dgeo <- cal_Dgeo(pred.geo_coord = pred.geo_coord, geo_coord = geo_coord, scalar = s)

  # rescale Dgeo
  # -> maximum likelihood estimation would run into an error if the values of Dgeo is too large
  # -> rescaling does not influence our outlier test
  maxD <- max(Dgeo)
  if(maxD > 20){
    tmps <- 1
    while (maxD > 20) {
      tmps <- tmps * 10
      maxD <- maxD / 10
    }
    orig_s <- s
    s <- s * tmps # new scalar
    cat(paste0("\n\n GGoutlieR adjusts the given scalar `s` value from `s=", orig_s, "` to `s=", s, "` to prevent an error in the maximum likelihood estimation process\n"))
    cat(paste0("\n\n D geo is re-scaled to a unit of ",s," meters \n\n"))
  }
  Dgeo <- cal_Dgeo(pred.geo_coord = pred.geo_coord, geo_coord = geo_coord, scalar = s)
  if(any(Dgeo == 0)){
    tmp.indx <- Dgeo == 0
    cat(paste("\n\n",sum(tmp.indx),"samples have D geo = 0. Zeros are replaced with 10^-8 to prevent the error in maximum likelihood estimation."))
    Dgeo[tmp.indx] <- 10^-8
  }

  predR2 <- diag(cor(pred.geo_coord, geo_coord))^2
  if(min(predR2) >= warning_minR2){
    warning(paste0("\n\n\n The lowest prediction accuracy (R^2) of KNN among two geographical dimensions (according to the given `geo_coord`) is ", round(min(predR2), digits = 4),
                   ". Maybe only few extreme outliers in your samples. \nYou could manually check the results with `plot_GGoutlieR` and adjust its `p_thres` argument to see which threshold is more appropriate. \n\n\n\n"))
  }


  ##---------------------------------------------------
  ##---------------------------------------------------
  ## KNN in geographical space: find the optimal K
  ##---------------------------------------------------
  ##---------------------------------------------------
  ## NOTE: this step would be skipped if `geoKNN_output` is available

  if(is.null(k_geoKNN)){
    # automatically select k if k=NULL
    cat(paste("\n\n `k_geoKNN` is NULL; searching for optimal k between", klim[1], "and", klim[2],"\nthis process can take a lot of time...\n"))
    all.D <- find_optimalK_geoKNN(geo_coord = geo_coord,
                                  gen_coord = gen_coord,
                                  geo.dM = geo.dM,
                                  w_power = w_power,
                                  klim = klim,
                                  do_par = do_par,
                                  min_nn_dist = min_nn_dist,
                                  cl = cl)

    # make a figure for K searching procedure
    opt.k = c(klim[1]:klim[2])[which.min(all.D)]
    k_geoKNN = opt.k # replace k with the optimal k
    k.sel.plot <- paste(plot_dir, "/KNN_Dg_optimal_k_selection.pdf", sep = "")
    pdf(k.sel.plot, width = 5, height = 4)
    par(mar=c(4,6,1,1))
    plot(x = klim[1]:klim[2], y = all.D, xlab="K", ylab=expression(sum(D[i], i==1, n)))
    abline(v = opt.k)
    legend("top",legend = paste("optimal k =", opt.k), pch="", bty = "n",cex = 1.2)
    invisible(dev.off())

    cat(paste("\n The optimal k_geoKNN is ",opt.k,". Its figure is saved at ", k.sel.plot," \n", sep = ""))
  }
  if(is.null(k_geoKNN)){stop("k_geoKNN is NULL!")}


  ##---------------------------------------------------
  ##---------------------------------------------------
  ## KNN in geographical space: prediction with the optimal K (or K given by users)
  ##---------------------------------------------------
  ##---------------------------------------------------
  cat("\nSearching K nearest neighbors...\n")
  knn.indx_geoKNN <- find_geo_knn(geo.dM = geo.dM, k=k_geoKNN, min_nn_dist=min_nn_dist)
  pred.q <- pred_q_knn(geo_coord = geo_coord, gen_coord = gen_coord, geo.dM =  geo.dM, knn.indx_geoKNN, w_power = w_power)
  # calculate Dg statistic
  Dg <- cal_Dg(pred.q, gen_coord)
  if(any(Dg == 0)){
    tmp.indx <- Dg == 0
    cat(paste("\n\n", sum(tmp.indx),"samples have Dg = 0. Zeros are replaced with 10^-8 to avoid the error in maximum likelihood estimation."))
    Dg[tmp.indx] <- 10^-8
  }

  # evaluate prediction accuracy as R^2
  predR2 <- diag(cor(pred.q, gen_coord))^2
  if(min(predR2) >= warning_minR2){
    warning(paste0("\n\n\nThe lowest prediction accuracy (R^2) of KNN among ", ncol(gen_coord),
                   " genetic dimensions (according to the given `gen_coord`) is ", round(min(predR2), digits = 4),
                   ". Maybe only few extreme outliers in your samples. \nYou could manually check the results with `plot_GGNet_map` and adjust its `p_thres` argument to see which threshold is more appropriate.\n\n\n\n"))
  }


  ##-------------------------------------------------------------------
  ##-------------------------------------------------------------------
  ## Get null distribution
  ##-------------------------------------------------------------------
  ##-------------------------------------------------------------------
  ### Define negative log likelihood function
  negLL_geneticKNN <- function(a_geneticKNN, b_geneticKNN){
    -sum(dgamma(Dgeo, shape = a_geneticKNN,
                      rate = b_geneticKNN, log = T))
  }
  negLL_geoKNN <- function(a_geoKNN, b_geoKNN){
    -sum(dgamma(Dg, shape = a_geoKNN,
                  rate = b_geoKNN, log = T))
  }


  ### setup initial values
  if(!is.null(gamma_par_geneticKNN)){
    initial.a_geneticKNN <- gamma_par_geneticKNN[1]
    initial.b_geneticKNN <- gamma_par_geneticKNN[2]
  }else{
    initial.a_geneticKNN <- (mean(Dgeo)^2)/var(Dgeo)
    initial.b_geneticKNN <- mean(Dgeo)/var(Dgeo)
  }
  if(!is.null(gamma_par_geoKNN)){
    initial.a_geoKNN <- gamma_par_geoKNN[1]
    initial.b_geoKNN <- gamma_par_geoKNN[2]
  }else{
    initial.a_geoKNN <- (mean(Dgeo)^2)/var(Dgeo)
    initial.b_geoKNN <- mean(Dgeo)/var(Dgeo)
  }

  ## Maximum likelihood estimation
  mle.res_geneticKNN <- mle(minuslogl = negLL_geneticKNN,
                            start = list(a_geneticKNN = initial.a_geneticKNN,
                                         b_geneticKNN = initial.b_geneticKNN),
                            lower = 10^-8, upper = 10^8, method = "L-BFGS-B")
  mle.res_geoKNN <- mle(minuslogl = negLL_geoKNN,
                        start = list(a_geoKNN = initial.a_geoKNN,
                                     b_geoKNN = initial.b_geoKNN),
                        lower = 10^-8, upper = 10^8, method = "L-BFGS-B")
  current.a_geneticKNN <- unname(mle.res_geneticKNN@coef["a_geneticKNN"])
  current.b_geneticKNN <- unname(mle.res_geneticKNN@coef["b_geneticKNN"])
  current.a_geoKNN <- unname(mle.res_geoKNN@coef["a_geoKNN"])
  current.b_geoKNN <- unname(mle.res_geoKNN@coef["b_geoKNN"])


  #--------------------------------------------------
  # make a figure for the null Gamma distribution
  plot_nullDistr <- function(current.a,
                             current.b,
                             n,
                             p_thres,
                             plot_dir,
                             prefix, # prefix of the output plot PDF file
                             plot_main # plot title
                             ){
    null.fun <- function(x){dgamma(x , shape = current.a, rate = current.b)}
    gamma.thres <- qgamma(1-p_thres, shape = current.a, rate = current.b)
    null.distr <- rgamma(n , shape = current.a, rate = current.b)
    null.plot <- paste0(plot_dir, "/",prefix,"_null_distribution.pdf")
    cat(paste("\nThe plot of null distribution is saved at ", null.plot," \n", sep = ""))
    pdf(null.plot, width = 5, height = 4)
    curve(null.fun, from = 0, to = max(null.distr), add = F, col = "blue",
          ylab = "Density",
          xlab = bquote(Gamma ~ '(' ~ alpha ~'='~.(round(current.a, digits = 3))~','~beta~'='~.(round(current.b, digits = 3))~')' )
          ,bty = "n", main = plot_main)
    abline(v = gamma.thres, col = "red", lty = 2)
    text(x = gamma.thres, y = par("usr")[4], labels = paste("p =", p_thres), xpd=NA)
    par(mgp = c(3,0,0))
    axis(side = 1 ,at = round(gamma.thres, digits = 3), line = 0.3, tck = 0.02, font = 2)
    invisible(dev.off())
  }
  plot_nullDistr(current.a = current.a_geneticKNN,
                 current.b = current.b_geneticKNN,
                 n = n,
                 p_thres = p_thres,
                 plot_dir = plot_dir,
                 prefix = "KNN_Dgeo",
                 plot_main = expression("Null distribution of D"[geo]))

  plot_nullDistr(current.a = current.a_geoKNN,
                 current.b = current.b_geoKNN,
                 n = n,
                 p_thres = p_thres,
                 plot_dir = plot_dir,
                 prefix = "KNN_Dg",
                 plot_main = expression("Null distribution of D"[g]))




  #----------------------------
  # make output (ggoutlier_geneticKNN)
  gamma.thres_geneticKNN <- qgamma(1-p_thres, shape = current.a_geneticKNN, rate = current.b_geneticKNN)
  sig.indx_geneticKNN <- Dgeo > gamma.thres_geneticKNN
  p.value_geneticKNN <- 1 - pgamma(Dgeo, shape = current.a_geneticKNN, rate = current.b_geneticKNN)
  out_geneticKNN <- data.frame("Dgeo" = Dgeo,
                               "p.value" = p.value_geneticKNN,
                               "significant" = sig.indx_geneticKNN)
  rownames(out_geneticKNN) <- rownames(geo_coord)
  thres_geneticKNN <- data.frame(pvalue = c(0.05,0.01,0.005,0.001),statistic = qgamma(1 - c(0.05,0.01,0.005,0.001), shape = current.a_geneticKNN, rate = current.b_geneticKNN))
  gamma.par_geneticKNN <- c(current.a_geneticKNN, current.b_geneticKNN)
  names(gamma.par_geneticKNN) <- c("alpha", "beta")
  rownames(knn.indx_geneticKNN) <- rownames(geo_coord)
  knn.name_geneticKNN <- apply(knn.indx_geneticKNN,2, function(x){rownames(geo_coord)[x]})
  res.out_geneticKNN <- list(out_geneticKNN,
                            thres_geneticKNN,
                            gamma.par_geneticKNN,
                            knn.indx_geneticKNN,
                            knn.name_geneticKNN, s) # s is the scalar of geographical distance
  names(res.out_geneticKNN) <- c("statistics","threshold","gamma_parameter", "knn_index", "knn_name", "scalar")

  #--------------------------------------------
  # make output (ggoutlier_geoKNN)
  gamma.thres_geoKNN <- qgamma(1-p_thres, shape = current.a_geoKNN, rate = current.b_geoKNN)
  sig.indx_geoKNN <- Dg > gamma.thres_geoKNN
  p.value_geoKNN <- 1 - pgamma(Dg, shape = current.a_geoKNN, rate = current.b_geoKNN)
  out_geoKNN <- data.frame("Dg" = Dg,
                           "p.value" = p.value_geoKNN,
                           "significant" = sig.indx_geoKNN)
  rownames(out_geoKNN) <- rownames(geo_coord)
  thres_geoKNN <- data.frame(pvalue = c(0.05,0.01,0.005,0.001),statistic = qgamma(1 - c(0.05,0.01,0.005,0.001), shape = current.a_geoKNN, rate = current.b_geoKNN))
  gamma.par_geoKNN <- c(current.a_geoKNN, current.b_geoKNN)
  names(gamma.par_geoKNN) <- c("alpha", "beta")
  rownames(knn.indx_geoKNN) <- rownames(geo_coord)
  knn.name_geoKNN <- apply(knn.indx_geoKNN,2, function(x){rownames(geo_coord)[x]})
  res.out_geoKNN <- list(out_geoKNN,
                         thres_geoKNN,
                         gamma.par_geoKNN,
                         knn.indx_geoKNN,
                         knn.name_geoKNN)
  names(res.out_geoKNN) <- c("statistics","threshold","gamma_parameter", "knn_index", "knn_name")

  if(!multi_stages){
    attr(res.out_geneticKNN, "model") <- "ggoutlier_geneticKNN"
    attr(res.out_geoKNN, "model") <- "ggoutlier_geoKNN"
    return(list("geneticKNN_result" = res.out_geneticKNN,
                "geoKNN_result" = res.out_geoKNN))
  }else{
    #-------------------------------------------------------------
    # multi-stage test
    cat(paste0("\n\nStart the multi-stage `geneticKNN` test process with k=",k_geneticKNN,
                   " using a null Gamma distribution with shape=", round(current.a_geneticKNN, digits = 3),
                   " and rate=",round(current.b_geneticKNN, digits = 3),
                   " (the parameters of Gamma distribution were determined by MLE)\n"))
    cat(paste0("Start the multi-stage `geoKNN` test process with k=",k_geoKNN,
               " using a null Gamma distribution with shape=", round(current.a_geoKNN, digits = 3),
               " and rate=",round(current.b_geoKNN, digits = 3),
               " (the parameters of Gamma distribution were determined by MLE)\n\n"))

    i=1
    run_geneticKNN <- min(res.out_geneticKNN$statistics$p.value) <= min(res.out_geoKNN$statistics$p.value)

    if(run_geneticKNN){
      to_keep <- res.out_geneticKNN$statistics$p.value > min(res.out_geneticKNN$statistics$p.value)
    }else{
      to_keep <- res.out_geoKNN$statistics$p.value > min(res.out_geoKNN$statistics$p.value)
    }

    tmp.geo.dM <- geo.dM[to_keep, to_keep]
    tmp.pgdM <- pgdM[to_keep, to_keep]
    tmp.geo_coord <- geo_coord[to_keep,]
    tmp.gen_coord <- gen_coord[to_keep,]

    res.Iters_geneticKNN <- list(res.out_geneticKNN)
    res.Iters_geoKNN <- list(res.out_geoKNN)
    # if `maxIter` is NULL -> let it equal to 50% of sample size
    if(is.null(maxIter)){maxIter <- round(nrow(gen_coord) * 0.5)}
    while (i <= maxIter) {
      if(i > 1){
        tmp.pgdM <- tmp.pgdM[to_keep, to_keep]
        tmp.geo.dM <- tmp.geo.dM[to_keep, to_keep]
        tmp.geo_coord <- tmp.geo_coord[to_keep,]
        tmp.gen_coord <- tmp.gen_coord[to_keep,]
      }

      cat(paste0("Iteration ", i,"\r"), append = F)

      # find KNN
      tmp.knn.indx_geneticKNN <- find_gen_knn(tmp.pgdM, k=k_geneticKNN)
      tmp.knn.indx_geoKNN <- find_geo_knn(geo.dM = tmp.geo.dM, k=k_geoKNN, min_nn_dist=min_nn_dist)

      # KNN prediction
      tmp.pred.geo_coord <- pred_geo_coord_knn(geo_coord = tmp.geo_coord,
                            pgdM = tmp.pgdM,
                            knn.indx = tmp.knn.indx_geneticKNN,
                            w_power = w_power)
      tmp.pred.q <- pred_q_knn(geo_coord = tmp.geo_coord,
                               gen_coord = tmp.gen_coord,
                               geo.dM = tmp.geo.dM,
                               knn.indx = tmp.knn.indx_geoKNN,
                               w_power = w_power)
      # calculate Dgeo statistic and p values
      tmp.Dgeo <- cal_Dgeo(pred.geo_coord = tmp.pred.geo_coord,
                           geo_coord = tmp.geo_coord,
                           scalar = s)
      tmp.p.value_geneticKNN <- 1 - pgamma(tmp.Dgeo, shape = current.a_geneticKNN, rate = current.b_geneticKNN)


      # calculate Dg statistic and p values
      tmp.Dg <- cal_Dg(tmp.pred.q, tmp.gen_coord)
      tmp.p.value_geoKNN <- 1 - pgamma(tmp.Dg, shape = current.a_geoKNN, rate = current.b_geoKNN)

      # get `to_keep` flags ->
      #      if the most significant p value from `geneticKNN`,
      #      create `to_keep` flags according to `geneticKNN` output.
      #      if not, create `to_keep` flags according to `geoKNN` output.
      run_geneticKNN <- min(tmp.p.value_geneticKNN) <= min(tmp.p.value_geoKNN)
      if(run_geneticKNN){
        to_keep <- tmp.p.value_geneticKNN > min(tmp.p.value_geneticKNN)
      }else{
        to_keep <- tmp.p.value_geoKNN > min(tmp.p.value_geoKNN)
      }

      if(all(tmp.p.value_geneticKNN > p_thres) & all(tmp.p.value_geoKNN > p_thres)){
        cat("\nNo new significant sample is identified. Multi-stage testing process ends...\n")
        break
      }else{
        # save the results from the current iteration
        ## make output of geneticKNN
        tmp.out_geneticKNN <- data.frame(tmp.Dgeo, tmp.p.value_geneticKNN, significant = tmp.p.value_geneticKNN < p_thres)
        colnames(tmp.out_geneticKNN) <- c("Dgeo", "p.value", "significant")
        rownames(tmp.out_geneticKNN) <- rownames(tmp.geo_coord)
        rownames(tmp.knn.indx_geneticKNN) <- rownames(tmp.geo_coord)
        tmp.knn.name_geneticKNN <- apply(tmp.knn.indx_geneticKNN,2,
                              function(x){rownames(tmp.geo_coord)[x]})
        tmp.res_geneticKNN <- list(tmp.out_geneticKNN,
                                   thres_geneticKNN,
                                   gamma.par_geneticKNN,
                                   tmp.knn.indx_geneticKNN,
                                   tmp.knn.name_geneticKNN, s)
        names(tmp.res_geneticKNN) <- c("statistics","threshold","gamma_parameter", "knn_index", "knn_name", "scalar")
        res.Iters_geneticKNN <- c(res.Iters_geneticKNN, list(tmp.res_geneticKNN))

        ## make output of geoKNN
        tmp.out_geoKNN <- data.frame(tmp.Dg, tmp.p.value_geoKNN, significant = tmp.p.value_geoKNN < p_thres)
        colnames(tmp.out_geoKNN) <- c("Dg", "p.value", "significant")
        rownames(tmp.out_geoKNN) <- rownames(tmp.geo_coord)
        rownames(tmp.knn.indx_geoKNN) <- rownames(tmp.geo_coord)
        tmp.knn.name_geoKNN <- apply(tmp.knn.indx_geoKNN,2, function(x){rownames(tmp.geo_coord)[x]})
        tmp.res_geoKNN <- list(tmp.out_geoKNN, thres_geoKNN, gamma.par_geoKNN, tmp.knn.indx_geoKNN, tmp.knn.name_geoKNN)
        names(tmp.res_geoKNN) <- c("statistics","threshold","gamma_parameter", "knn_index", "knn_name")
        res.Iters_geoKNN <- c(res.Iters_geoKNN, list(tmp.res_geoKNN))
      }
      i = i+1

    } # while loop end

    ## collect results of all iterations
    ## NOTE: In each iteration of 'multi-stage test', the sample with the most significant p values will be exclude from the KNN searching procedure in the next iteration
    ##       The loop here sequentially collect the outputs from each iteration and update the data.frame `collapse_res`

    ### processing geneticKNN output
    collapse_res_geneticKNN <- res.Iters_geneticKNN[[1]]
    if(length(res.Iters_geneticKNN) > 1){
      for(i in 2:length(res.Iters_geneticKNN)){
        tmp <- res.Iters_geneticKNN[[i]]
        tmp$statistics <- tmp$statistics[match(rownames(collapse_res_geneticKNN$statistics),
                                               rownames(tmp$statistics)),]
        tmp$knn_name <- tmp$knn_name[match(rownames(collapse_res_geneticKNN$statistics),
                                           rownames(tmp$knn_index)),]

        tmp.indx <- which(!is.na(tmp$statistics$Dgeo))
        if(length(tmp.indx)>0){
          collapse_res_geneticKNN$statistics[tmp.indx,] <- tmp$statistics[tmp.indx,]
          collapse_res_geneticKNN$knn_name[tmp.indx,] <- tmp$knn_name[tmp.indx,]
          for(j in tmp.indx){
            collapse_res_geneticKNN$knn_index[j,] <- match(collapse_res_geneticKNN$knn_name[j,], rownames(collapse_res_geneticKNN$statistics))
          }
        }
      }
    }

    ### processing geoKNN output
    collapse_res_geoKNN <- res.Iters_geoKNN[[1]]
    if(length(res.Iters_geoKNN) > 1){
      for(i in 2:length(res.Iters_geoKNN)){
        tmp <- res.Iters_geoKNN[[i]]
        tmp$statistics <- tmp$statistics[match(rownames(collapse_res_geoKNN$statistics),
                                               rownames(tmp$statistics)),]
        tmp$knn_name <- tmp$knn_name[match(rownames(collapse_res_geoKNN$statistics),
                                           rownames(tmp$knn_index)),]
        tmp.indx <- which(!is.na(tmp$statistics$Dg))
        if(length(tmp.indx)>0){
          collapse_res_geoKNN$statistics[tmp.indx,] <- tmp$statistics[tmp.indx,]
          collapse_res_geoKNN$knn_name[tmp.indx,] <- tmp$knn_name[tmp.indx,]
          for(j in tmp.indx){
            collapse_res_geoKNN$knn_index[j,] <- match(collapse_res_geoKNN$knn_name[j,], rownames(collapse_res_geoKNN$statistics))
          }
        }
      }
    }

    # make a figure comparing the results of single stage and multi-stage tests
    ## plot of geneticKNN
    logp.plot <- paste0(plot_dir, "/geneticKNN_test_multi_stage_Log10P_comparison.pdf")
    cat(paste("\n\n\nThe plot for comparing -logP between single-stage and multi-stage KNN tests is saved at ", logp.plot," \n", sep = ""))
    pdf(logp.plot, width = 4, height = 4.2)
    plot(-log10(res.out_geneticKNN$statistics$p.value),
         -log10(collapse_res_geneticKNN$statistics$p.value),
         xlab = expression("-log"[10]~"(p) of single-stage KNN test"),
         ylab = expression("-log"[10]~"(p) of multi-stage KNN test"),
         main = "KNN in Genetic space")
    invisible(dev.off())

    ## plot of geoKNN
    logp.plot <- paste0(plot_dir, "/geoKNN_test_multi_stage_Log10P_comparison.pdf")
    cat(paste("\n\n\nThe plot of comparing -logP between single-stage and multi-stage KNN tests is saved at ", logp.plot," \n", sep = ""))
    pdf(logp.plot, width = 4, height = 4.2)
    plot(-log10(res.out_geoKNN$statistics$p.value),
         -log10(collapse_res_geoKNN$statistics$p.value),
         xlab = expression("-log"[10]~"(p) of single-stage KNN test"),
         ylab = expression("-log"[10]~"(p) of multi-stage KNN test"),
         main = "KNN in Geographical space")
    invisible(dev.off())

    # reture multi-stage test outputs
    if(keep_all_stg_res){
      ## geneticKNN output
      names(res.Iters_geneticKNN) <- paste0("Iter_", 1:length(res.Iters_geneticKNN))
      res.out_geneticKNN <- c(res.Iters_geneticKNN, collapse_res_geneticKNN)
      out_geneticKNN <- c(res.Iters_geneticKNN, combined_result = list(collapse_res_geneticKNN))
      attr(out_geneticKNN, "model") <- "ggoutlier_geneticKNN"

      ## geoKNN output
      names(res.Iters_geoKNN) <- paste0("Iter_", 1:length(res.Iters_geoKNN))
      res.out_geoKNN <- c(res.Iters_geoKNN, collapse_res_geoKNN)
      out_geoKNN <- c(res.Iters_geoKNN, combined_result = list(collapse_res_geoKNN))
      attr(out_geoKNN, "model") <- "ggoutlier_geoKNN"

      return(list("geneticKNN_result" = out_geneticKNN,
                  "geoKNN_result" = out_geoKNN))
    }else{
      attr(collapse_res_geneticKNN, "model") <- "ggoutlier_geneticKNN"
      attr(collapse_res_geoKNN, "model") <- "ggoutlier_geoKNN"
      return(list("geneticKNN_result" = collapse_res_geneticKNN,
                  "geoKNN_result" = collapse_res_geoKNN))
    }
  } # multi-stage test end
} # ggoutlier_compositeKNN end
