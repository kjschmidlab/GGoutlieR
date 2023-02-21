#' GGoutlieR with the geographical KNN approach
#' @description identify samples genetically different from K nearest geographical neighbors (geographical KNN)
#' @param geo_coord matrix or data.frame with two columns. The first column is longitude and the second one is latitude.
#' @param gen_coord matrix. A matrix of "coordinates in a genetic space". Users can provide ancestry coefficients or eigenvectors for calculation. If, for example, ancestry coefficients are given, each column corresponds to an ancestral population. Samples are ordered in rows as in `geo_coord`.
#' @param k integer. Number of the nearest neighbors.
#' @param k_geneticKNN integer. Number of the nearest neighbors in a genetic space. The default is `NULL`. The `ggoutlier` will search the optimal K if `k_geneticKNN = NULL`.
#' @param k_geoKNN integer. Number of the nearest neighbors in a geographical space. the default is `NULL`. The `ggoutlier` will search the optimal K if `k_geoKNN = NULL`.
#' @param klim vector. A range of K to search for the optimal number of nearest neighbors. The default is `klim = c(3, 50)`
#' @param plot_dir string. The path to save plots
#' @param w_power numeric. A value controlling the power of distance weight in geographical KNN prediction.
#' @param p_thres numeric. A significance level
#' @param s integer. A scalar of geographical distance. The default `s=100` scales the distance to a unit of 0.1 kilometer.
#' @param cpu integer. Number of CPUs to use for searching the optimal K.
#' @param verbose logic. If `verbose = FALSE`, `ggoutlier` will suppress printout messages.
#' @param multi_stages logic. A multi-stage test will be performed if is `TRUE` (the default is `TRUE`).
#' @param maxIter numeric. Maximal iteration number of multi-stage KNN test.
#' @param keep_all_stg_res logic. Results from all iterations of the multi-stage test will be retained if it is`TRUE`. (the default is `FALSE`)
#' @param warning_minR2 numeric. The prediction accuracy of KNN is evaluated as R^2 to assess the violation of isolation-by-distance expectation. If any R^2 is larger than `warning_minR2`, a warning message will be reported at the end of your analysis.
#' @return a list including five items. `statistics` is a `data.frame` consisting of the D_g values, p values and a column of logic values showing if a sample is an outlier or not. `threshold` is a `data.frame` recording the significance threshold. `gamma_parameter` is a vector recording the parameter of the heuristic Gamma distribution. `knn_index` and `knn_name` are a `data.frame` recording the K nearest neighbors of each sample.
#' @export
ggoutlier_geoKNN <- function(geo_coord,
                             gen_coord,
                             min_nn_dist = 100,
                             k = NULL,
                             klim = c(3,50),
                             s = 100,
                             plot_dir = ".",
                             w_power = 1,
                             p_thres = 0.05,
                             n = 10^6,
                             multi_stages = TRUE,
                             maxIter=NULL,
                             keep_all_stg_res = FALSE,
                             warning_minR2 = 0.9,
                             cpu = 1,
                             verbose = TRUE
){
  required_pkgs <- c("geosphere", # for calculating geographical distances
                     "stats4", # package to perform maximum likelihood estimation
                     "FastKNN", # KNN algorithm using a given distance matrix (other packages do not take arbitrary distance matrices)
                     "foreach", "doParallel",
                     "iterators","parallel")
  invisible(lapply(required_pkgs, FUN=function(x){suppressPackageStartupMessages(library(x, verbose = FALSE, character.only = TRUE))}))

  if(cpu > 1){
    max_cores=detectCores()
    if(cpu >= max_cores){
      warning(paste0("\n The maximum number of CPUs is ", max_cores, ". Set `cpu` to ", max_cores-1," \n"))
      cpu <- max_cores -1 # reserve max-1 cpus if users request too many cpus
    }
    if(cpu > 1){
      if(verbose) cat(paste0("\n Parallelize computation using ", cpu, " cores \n"))
      cl <- makeCluster(cpu)
      do_par <- TRUE
    }
  } else {
    if(cpu < 1){stop("`cpu` has to be at least 1")}
    do_par <- FALSE
    cl <- NULL
    if(verbose) cat("\n Computation using 1 core. you can parallelize computation by setting a higher value to `cpu` argument \n")
  }
  # check inputs
  if(ncol(geo_coord) != 2){stop("Please ensure the `geo_coord` having two columns!")}
  if(nrow(geo_coord) != nrow(gen_coord)){stop("`geo_coord` and `gen_coord` has different sample size!")}

  if(is.null(rownames(geo_coord))){
    rownames(geo_coord) <- paste("sample", 1:nrow(geo_coord), sep = "")
  }else{
    if(any(is.na(rownames(geo_coord)))){
      if(verbose) cat("Some samples in `geo_coord` have IDs but some do not. Arbitratry IDs are assigned to all samples...")
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
  # Data pre-treatment
  ## calculate geographical distance
  geo.dM <- distm(x = geo_coord)/s
  ## handle samples with identical geographical coordinates
  if(any(geo.dM[lower.tri(geo.dM)] == 0)){
    if(verbose) cat("Find samples with identical geographical coordinates.\n")
    if(any(min_nn_dist == 0 , is.null(min_nn_dist))){
      if(verbose) cat("Add one unit of distance to individual pairs\n")
      geo.dM <- geo.dM + 1
      diag(geo.dM) <- 0
    }else{
      if(verbose) cat("Ignore neighbors whose pairwise distance is less than `min_nn_dist` when searching KNNs.\n")
      orig_min_nn_dist <- min_nn_dist
      min_nn_dist <- orig_min_nn_dist/s
      if(verbose) cat(paste0("\n\nIgnore neighbors within ",min_nn_dist, " unit(s) of distance (the default unit of distance is m; `min_nn_dist` is adjusted automatically according to `s`)\n"))
    }
  }

  #---------------------------------------
  # Search the optimal K for KNN
  if(is.null(k)){
    # automatically select k if k=NULL
    if(verbose) cat(paste("\n\n `k` is NULL; searching for optimal k between", klim[1], "and", klim[2],"\nthis process can take a lot of time...\n"))
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
    k = opt.k # replace k with the optimal k
    k.sel.plot <- paste(plot_dir, "/KNN_Dg_optimal_k_selection.pdf", sep = "")
    pdf(k.sel.plot, width = 5, height = 4)
    par(mar=c(4,6,1,1))
    plot(x = klim[1]:klim[2], y = all.D, xlab="K", ylab=expression(sum(D["genetic,i"], i==1, n)))
    abline(v = opt.k)
    legend("top",legend = paste("optimal k =", opt.k), pch="", bty = "n",cex = 1.2)
    dev.off()


    if(verbose) cat(paste("\n The optimal k is ",opt.k,". Its figure is saved at ", k.sel.plot," \n", sep = ""))

  }
  if(is.null(k)){stop("k is NULL!")}
  if(verbose) cat("\nSearching K nearest neighbors...\n")

  #-----------------------------------
  # KNN prediction with the optimal K (or K given by users)
  knn.indx <- find_geo_knn(geo.dM = geo.dM, k=k, min_nn_dist=min_nn_dist)
    pred.q <- pred_q_knn(geo_coord = geo_coord, gen_coord = gen_coord, geo.dM =  geo.dM, knn.indx, w_power = w_power)
  # calculate Dg statistic
  Dg <- cal_Dg(pred.q, gen_coord)
  if(any(Dg == 0)){
    tmp.indx <- Dg == 0
    if(verbose) cat(paste("\n\n", sum(tmp.indx),"samples have Dg = 0. Zeros are replaced with 10^-8 to avoid the error in maximum likelihood estimation."))
    Dg[tmp.indx] <- 10^-8
  }

  # evaluate prediction accuracy as R^2
  predR2 <- diag(cor(pred.q, gen_coord))^2
  if(min(predR2) >= warning_minR2){
    warning(paste0("\n\n\nThe lowest prediction accuracy (R^2) of KNN among ", ncol(gen_coord),
                   " genetic dimensions (according to the given `gen_coord`) is ", round(min(predR2), digits = 4),
                   ". Maybe only few extreme outliers in your samples. \nYou could manually check the results with `plot_GGNet_map` and adjust its `p_thres` argument to see which threshold is more appropriate.\n\n\n\n"))
  }
  #----------------------------------------------------------
  # Get null distribution
  ## Maximum likelihood assuming Gamma distribution
  ### Define negative log likelihood function
  negLL <- function(a, b){
    -sum(dgamma(Dg, shape = a, rate = b, log = T))
  }
  ### Loop until parameters converging
  initial.a <- (mean(Dg)^2)/var(Dg)
  initial.b <- mean(Dg)/var(Dg)
  #if(verbose) cat("\n\n Using maximum likelihood estimation to infer the null Gamma distribution...\n")
  mle.res <- mle(minuslogl = negLL, start = list(a = initial.a, b = initial.b),
                 lower = 10^-8, upper = 10^8, method = "L-BFGS-B")
  current.a <- unname(mle.res@coef["a"])
  current.b <- unname(mle.res@coef["b"])

  # get significance threshold
  null.fun <- function(x){dgamma(x , shape = current.a, rate = current.b)}
  gamma.thres <- qgamma(1-p_thres, shape = current.a, rate = current.b)
  null.distr <- rgamma(n , shape = current.a, rate = current.b)
  null.plot <- paste(plot_dir, "/KNN_Dg_null_distribution.pdf", sep = "")
  if(verbose) cat(paste("\nThe plot of null distribution is saved at ", null.plot," \n", sep = ""))

  ## make a plot for the null distribution
  pdf(null.plot, width = 5, height = 4)
  curve(null.fun, from = 0, to = max(null.distr), add = F, col = "blue",
        ylab = "Density",
        xlab = bquote(Gamma ~ '(' ~ alpha ~'='~.(round(current.a, digits = 3))~','~beta~'='~.(round(current.b, digits = 3))~')' )
        ,bty = "n", main = expression("Null distribution of D"[g]))
  abline(v = gamma.thres, col = "red", lty = 2)
  text(x = gamma.thres, y = par("usr")[4], labels = paste("p =", p_thres), xpd=NA)
  par(mgp = c(3,0,0))
  axis(side = 1 ,at = round(gamma.thres, digits = 3), line = 0.3, tck = 0.02, font = 2)
  dev.off()

  #--------------------------------------------
  # Return results
  sig.indx <- Dg > gamma.thres
  p.value <- 1 - pgamma(Dg, shape = current.a, rate = current.b)
  out <- data.frame(Dg, p.value, significant = sig.indx)
  rownames(out) <- rownames(geo_coord)
  thres <- data.frame(pvalue = c(0.05,0.01,0.005,0.001),statistic = qgamma(1 - c(0.05,0.01,0.005,0.001), shape = current.a, rate = current.b))
  gamma.par <- c(current.a, current.b)
  names(gamma.par) <- c("alpha", "beta")
  rownames(knn.indx) <- rownames(geo_coord)
  knn.name <- apply(knn.indx,2, function(x){rownames(geo_coord)[x]})
  res.out <- list(out, thres, gamma.par, knn.indx, knn.name)
  names(res.out) <- c("statistics","threshold","gamma_parameter", "knn_index", "knn_name")

  ## arguments used in GGoutlieR
  arguments <- c(s, w_power, k, min_nn_dist, multi_stages)
  names(arguments) <- c("scalar", "geo_weight_power", "K", "min_neighbor_dist", "multi_stage_test")

  if(!multi_stages){
    attr(res.out, "model") <- "ggoutlier_geoKNN"
    ## save arguments used in GGoutlieR
    attr(res.out, "arguments") <- arguments
    return(res.out)
  }else{
    #---------------------
    # multi-stage test
    if(verbose) cat(paste0("\n\nStart multi-stage KNN test process with k=",k,
                   " and using the null Gamma distribution with shape=", round(current.a, digits = 3),
                   " and rate=",round(current.b, digits = 3),
                   " (the parameters of Gamma distribution were determined by MLE)\n\n"))

    i = 1
    to_keep <- res.out$statistics$p.value > min(res.out$statistics$p.value)
    tmp.gen_coord <- gen_coord[to_keep,]
    tmp.geo.dM <- geo.dM[to_keep, to_keep]
    tmp.geo_coord <- geo_coord[to_keep,]

    res.Iters <- list(res.out)
    # if `maxIter` is NULL -> let it equal to 50% of sample size
    if(is.null(maxIter)){maxIter <- round(nrow(gen_coord) * 0.5)}
    while (i <= maxIter) {
      if(i > 1){
        tmp.gen_coord <- tmp.gen_coord[to_keep,]
        tmp.geo.dM <- tmp.geo.dM[to_keep, to_keep]
        tmp.geo_coord <- tmp.geo_coord[to_keep,]
      }

      if(verbose) cat(paste0("Iteration ", i,"\r"), append = F)

      # find KNN
      tmp.knn.indx <- find_geo_knn(geo.dM = tmp.geo.dM, k=k, min_nn_dist=min_nn_dist)
      # KNN prediction
      tmp.pred.q <- pred_q_knn(geo_coord = tmp.geo_coord,
                              gen_coord = tmp.gen_coord,
                              geo.dM = tmp.geo.dM,
                              knn.indx = tmp.knn.indx,
                              w_power = w_power)
      # calculate Dg statistic
      tmp.Dg <- cal_Dg(tmp.pred.q, tmp.gen_coord)
      tmp.p.value <- 1 - pgamma(tmp.Dg, shape = current.a, rate = current.b)
      to_keep <- tmp.p.value > min(tmp.p.value)


      if(all(tmp.p.value > p_thres)){
        if(verbose) cat("\nNo new significant sample is identified. Multi-stage testing process ends...\n")
        break
      }else{
        out <- data.frame(tmp.Dg, tmp.p.value, significant = tmp.p.value < p_thres)
        colnames(out) <- c("Dg", "p.value", "significant")
        rownames(out) <- rownames(tmp.geo_coord)
        thres <- data.frame(pvalue = c(0.05,0.01,0.005,0.001),statistic = qgamma(1 - c(0.05,0.01,0.005,0.001), shape = current.a, rate = current.b))
        gamma.par <- c(current.a, current.b)
        names(gamma.par) <- c("alpha", "beta")
        rownames(tmp.knn.indx) <- rownames(tmp.geo_coord)
        tmp.knn.name <- apply(tmp.knn.indx,2, function(x){rownames(tmp.geo_coord)[x]})
        tmp.res <- list(out, thres, gamma.par, tmp.knn.indx, tmp.knn.name)
        names(tmp.res) <- c("statistics","threshold","gamma_parameter", "knn_index", "knn_name")
        res.Iters <- c(res.Iters, list(tmp.res))
      }
      i = i+1

    } # while loop end

    ## collect results of all iterations
    collapse_res <- res.Iters[[1]]
    if(length(res.Iters) > 1){
      for(i in 2:length(res.Iters)){
        tmp <- res.Iters[[i]]
        tmp$statistics <- tmp$statistics[match(rownames(collapse_res$statistics),
                                               rownames(tmp$statistics)),]
        tmp$knn_name <- tmp$knn_name[match(rownames(collapse_res$statistics),
                                           rownames(tmp$knn_index)),]
        tmp.indx <- which(!is.na(tmp$statistics$Dg))
        if(length(tmp.indx)>0){
          collapse_res$statistics[tmp.indx,] <- tmp$statistics[tmp.indx,]
          collapse_res$knn_name[tmp.indx,] <- tmp$knn_name[tmp.indx,]
          for(j in tmp.indx){
            collapse_res$knn_index[j,] <- match(collapse_res$knn_name[j,], rownames(collapse_res$statistics))
          }
        }
      }
    }

    # make a figure comparing the results of single stage and multi-stage tests
    logp.plot <- paste0(plot_dir, "/geoKNN_test_multi_stage_Log10P_comparison.pdf")
    if(verbose) cat(paste("\n\n\nThe plot of comparing -logP between single-stage and multi-stage KNN tests is saved at ", logp.plot," \n", sep = ""))
    pdf(logp.plot, width = 4, height = 4.2)
    plot(-log10(res.out$statistics$p.value),
         -log10(collapse_res$statistics$p.value),
         xlab = expression("-log"[10]~"(p) of single-stage KNN test"),
         ylab = expression("-log"[10]~"(p) of multi-stage KNN test"),
         main = "KNN in Geographical space")
    dev.off()

    if(keep_all_stg_res){
      names(res.Iters) <- paste0("Iter_", 1:length(res.Iters))
      out <- c(res.Iters, combined_result = list(collapse_res))
      attr(out, "model") <- "ggoutlier_geoKNN"
      ## save arguments used in GGoutlieR
      attr(out, "arguments") <- arguments
      return(out)
    }else{
      attr(collapse_res, "model") <- "ggoutlier_geoKNN"
      ## save arguments used in GGoutlieR
      attr(collapse_res, "arguments") <- arguments
      return(collapse_res)
    }
  } # multi-stage test end
} # ggoutlier_geoKNN end

