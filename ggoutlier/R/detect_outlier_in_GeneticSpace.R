#' Identify samples geographically remote from K genetically nearest neighbors
#' @param geo_coord a two column matrix or data.frame. the first column is longitude and the second one is latitude.
#' @param gen_coord a matrix of "coordinates in a genetic space". Users can provide ancestry coefficients or eigenvectors for calculation. If, for example, ancestry coefficients are given, each column corresponds to an ancestral population. Samples are ordered in rows as in `geo_coord`. Users have to provide `pgdM` if `gen_coord` is not given.
#' @param pgdM a pairwise genetic distance matrix. Users can provide a customized genetic distance matrix with this argument. Samples are ordered in rows and columns as in the rows of `geo_coord`. The default of `pgdM` is `NULL`. If `pgdM` is not provided, a genetic distance matrix will be calculated from `gen_coord`.
#' @param k number of the nearest neighbor. the default is `NULL`.
#' @param klim if `k = NULL`, an optimal k will be searched between the first and second value of `klim`
#' @param s a scalar of geographical distance. The default `s=1000` scales the distance to a unit of 1 kilometer.
#' @param plot_dir the path to save plots
#' @param w_power a value controlling the power of distance weight in KNN prediction. For example, if `w_power=2`, the weight of KNN is 1/d^2/sum(1/d^2).
#' @param p_thres a significance level
#' @param n number of samples to draw from the null distribution (to obtain the range of x axis to make a curve plot of null distribution)
#' @param multi_stages logic. a multi-stage test will be performed if is `TRUE` (the default is `TRUE`).
#' @param maxIter maximal iteration number of multi-stage KNN test.
#' @param keep_all_stg_res logic. results from all iterations of the multi-stage test will be retained if it is`TRUE`. (the default is `FALSE`)
#' @param warning_minR2 the prediction accuracy of KNN is evaluated as R^2 to assess the violation of isolation-by-distance expectation. if any R^2 is less than `warning_minR2`, a warning message will be reported at the end of your analysis.
#' @return a list including five items. `statistics` is a `data.frame` consisting of the D_g values, p values and logic values showing if a sample is an outlier or not. `threshold` is a `data.frame` recording the significance threshold. `gamma_parameter` is a vector recording the parameter of the heuristic Gamma distribution. `knn_index` and `knn_name` are a `data.frame` recording the K nearest neighbors of each sample.
#' @export
#'
#'
#' @examples
#' # example 1 -> using ancestry coefficients:
#'  geo_coord <- read.table("./data/georef_1661ind_geo_coord_for_locator.txt", header = T, stringsAsFactors = F)
#'  rownames(geo_coord) <- geo_coord[,1]
#'  geo_coord <- geo_coord[,-1]
#'  anc.coef <- t(as.matrix(read.csv("./data/alstructure_Q_hat_1661inds.csv", header = F, stringsAsFactors = F)))
#'  sus2 <- ggoutlier_geneticspace(geo_coord = geo_coord, gen_coord = anc.coef,
#'                                     plot_dir = "./fig", p_thres = 0.05, w_power = 2,
#'                                     k = NULL, klim = c(3,100), keep_all_stg_res = F)
#'  example 2 -> using eigenvectors:
#'  pc <- read.table("data/gbs_georef_1661inds_PCA.eigenvec", stringsAsFactors = F)
#'  rownames(pc) <- gsub(pc[,2], pattern = "^0_", replacement = "")
#'  pc <- apply(pc[,-c(1:2)], 2, function(x){scale(x)}) # removing FID and IID and normalizing data
#'  sus2 <- ggoutlier_geneticspace(geo_coord = geo_coord, gen_coord = pc[,1:5],
#'                                        plot_dir = "./fig", p_thres = 0.05, w_power = 2,
#'                                        k = NULL, klim = c(3,100), keep_all_stg_res = F)

ggoutlier_geneticspace <- function(geo_coord, gen_coord = NULL, pgdM = NULL,
                                           k = NULL,
                                           klim = c(3,100),
                                           plot_dir = ".",
                                           w_power = 2,
                                           p_thres = 0.05,
                                           n = 10^6,
                                           s = 1000,
                                           multi_stages = T,
                                           maxIter=100,
                                           keep_all_stg_res = F,
                                           warning_minR2 = 0.9,
                                           cpu = 1
                                    ){
  require(geosphere) # for calculating geographical distances
  require(stats4) # package to perform maximum likelihood estimation
  require(FastKNN) # search KNN with an arbitrary distance matrix
  if(cpu > 1){
    require(foreach)
    require(doParallel)
    max_cores=detectCores()
    if(cpu >= max_cores){}
    message()
  } else {
    if(cpu < 1){stop("`cpu` has to be at least 1")}
    message()
  }

  # check inputs
  if(all(is.null(gen_coord) & is.null(pgdM))){
    stop("Either `gen_coord` or `pgdM` has to be provided!")
  }
  if(all(!is.null(gen_coord) & !is.null(pgdM))){
    warning("Both `gen_coord` and `pgdM` are provided. `pgdM` will be used instead of calculating genetic distances from `gen_coord`.\n")
  }
  if(ncol(geo_coord) != 2){stop("Please ensure the `geo_coord` having two columns!")}
  if(!is.null(gen_coord)){
    if(nrow(geo_coord) != nrow(gen_coord)){stop("`geo_coord` and `gen_coord` has different sample size!")}
  }else{
    if(nrow(geo_coord) != nrow(pgdM) | nrow(geo_coord) != ncol(pgdM)){stop("`geo_coord` and `pgdM` has different sample size!")}
  }

  if(is.null(rownames(geo_coord))){
    rownames(geo_coord) <- paste("sample", 1:nrow(geo_coord), sep = "")
  }else{
    if(any(is.na(rownames(geo_coord)))){
      message("Some samples in `geo_coord` have IDs but some do not. Arbitratry IDs are assigned to all samples...")
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

  # calculate genetic distance using ancestry coefficients
  if(is.null(pgdM)){
    pgdM <- as.matrix(dist(gen_coord, method = "euclidean"))
  }
  # handle samples with identical genotypes
  if(any(pgdM[lower.tri(pgdM)] == 0)){
    message("Find samples with zero genetic distances.\n")
    message("Add 10^-6 unit of distance to individual pairs\n")
    pgdM <- pgdM + 10^-6
    diag(pgdM) <- 0
  }

  # search KNN
  ## a function to find KNN
  find_gen_knn <- function(pgdM, k){
    indx <- 1:ncol(pgdM)
    res <- matrix(NA, nrow = nrow(pgdM), ncol = k)
    for(i in indx){
      #tmp.d <- cbind(indx, pgdM[,i])
      #tmp.d <- tmp.d[order(tmp.d[,2], decreasing = F),]
      # get the k nearest neighbors but with more than 0 unit dist
      #res[i,] <- tmp.d[tmp.d[,2] > 0,][1:k,1]
      #res[i,] <- tmp.d[tmp.d[,1] != i,][1:k,1]
      res[i,] <- FastKNN::k.nearest.neighbors(i = i, distance_matrix = pgdM, k = k)
    }
    return(res)
  } # find_gen_knn end

  ## a function to predict with KNN
  pred_geo_coord_knn <- function(geo_coord, pgdM, knn.indx){
    res <- matrix(NA, nrow = nrow(geo_coord), ncol = ncol(geo_coord))
    for(j in 1:nrow(pgdM)){
      tmp.indx <- knn.indx[j,]
      tmp.geo_coord <- geo_coord[tmp.indx,]
      tmp.d <- (pgdM[tmp.indx,j]) ^ w_power # new approach
      w <- (1/tmp.d)/(sum(1/tmp.d))
      res[j,] <- apply(tmp.geo_coord, 2, function(x){weighted.mean(x, w)})
    }
    return(res)
  } # pred_geo_coord_knn end

  ## a function to calculate Dg
  cal_Dgeo <- function(pred.geo_coord, geo_coord, scalar){
    geodist <- (sapply(1:nrow(geo_coord),function(a){distm(x = geo_coord[a,], y = pred.geo_coord[a,])/scalar}))
    return(geodist)
  } # cal_Dgeo end


  if(is.null(k)){
    # automatically select k if k=NULL
    message(paste("searching for optimal k between", klim[1], "and", klim[2],"\nthis process can take a lot of time...\n"))
    k = min(klim)
    all.D <- c()

    while(k <= max(klim)){
      message(paste0("calculating Dgeo with k = ",k,"\r"), appendLF = F)
      # find KNN
      knn.indx <- find_gen_knn(pgdM, k=k)
      # KNN predicttion
      pred.geo_coord <- pred_geo_coord_knn(geo_coord, pgdM, knn.indx)

      # calculate Dg statistic
      Dgeo <- cal_Dgeo(pred.geo_coord = pred.geo_coord, geo_coord = geo_coord, scalar = s)
      all.D <- c(all.D, sum(Dgeo))
      k=k+1
    }
    opt.k = c(klim[1]:klim[2])[which.min(all.D)]
    k.sel.plot <- paste(plot_dir, "/KNN_Dgeo_optimal_k_selection.pdf", sep = "")
    pdf(k.sel.plot, width = 5, height = 4)
    par(mar=c(4,6,1,1))
    plot(x = klim[1]:klim[2], y = all.D, xlab="K", ylab=expression(sum(D["geo,i"], i==1, n)))
    abline(v = opt.k)
    legend("top",legend = paste("optimal k =", opt.k), pch="", bty = "n",cex = 1.2)
    dev.off()

    k = opt.k
    message(paste("the optimal k is ",opt.k,". the figure is saved at ", k.sel.plot," \n", sep = ""))

  }
  if(is.null(k)){stop("k is NULL!")}
  message("searching K nearest neighbors...\n")
  knn.indx <- find_gen_knn(pgdM, k=k)
  # KNN predicttion
  pred.geo_coord <- pred_geo_coord_knn(geo_coord, pgdM, knn.indx)
  # calculate Dgeo statistic
  message(paste("D geo is scaled to a unit of",s,"meters \n"))
  Dgeo <- cal_Dgeo(pred.geo_coord = pred.geo_coord, geo_coord = geo_coord, scalar = s)
  if(any(Dgeo == 0)){
    tmp.indx <- Dgeo == 0
    message(paste(sum(tmp.indx),"samples have D geo = 0. Zeros are replaced with 10^-8 to avoid the error in maximum likelihood estimation."))
    Dgeo[tmp.indx] <- 10^-8
  }

  predR2 <- diag(cor(pred.geo_coord, geo_coord))^2
  if(min(predR2) >= warning_minR2){
    warning(paste0("\n\n\nThe lowest prediction accuracy (R^2) of KNN among two geographical dimensions (according to the given `geo_coord`) is ", round(min(predR2), digits = 4),
                   ". Maybe only few extreme outliers in your samples. \nYou could manually check the results with `plot_GGNet_map` and adjust its `p_thres` argument to see which threshold is more appropriate.\n\n\n\n"))
  }
  # get a null distribution
  ## Maximum likelihood assuming Gamma distribution
  ### Define negative log likelihood function
  negLL <- function(a, b){
    -sum(dgamma(Dgeo, shape = a, rate = b, log = T))
  }

  ### Loop until paramters converging
  initial.a <- (mean(Dgeo)^2)/var(Dgeo)
  initial.b <- mean(Dgeo)/var(Dgeo)
  message("Using maximum likelihood estimation to infer the null Gamma distribution...\n")
  mle.res <- mle(minuslogl = negLL, start = list(a = initial.a, b = initial.b),
                 lower = 10^-8, upper = 10^8, method = "L-BFGS-B")
  current.a <- unname(mle.res@coef["a"])
  current.b <- unname(mle.res@coef["b"])


  # get significance threshold
  null.fun <- function(x){dgamma(x , shape = current.a, rate = current.b)}
  gamma.thres <- qgamma(1-p_thres, shape = current.a, rate = current.b)
  null.distr <- rgamma(n , shape = current.a, rate = current.b)
  null.plot <- paste(plot_dir, "/KNN_Dgeo_null_distribution.pdf", sep = "")
  message(paste("the plot of null distribution is saved at ", null.plot," \n", sep = ""))
  pdf(null.plot, width = 5, height = 4)
  curve(null.fun, from = 0, to = max(null.distr), add = F, col = "blue",
        ylab = "Density",
        xlab = bquote(Gamma ~ '(' ~ alpha ~'='~.(round(current.a, digits = 3))~','~beta~'='~.(round(current.b, digits = 3))~')' )
        ,bty = "n", main = expression("Null distribution of D"[geo]))
  abline(v = gamma.thres, col = "red", lty = 2)
  text(x = gamma.thres, y = par("usr")[4], labels = paste("p =", p_thres), xpd=NA)
  par(mgp = c(3,0,0))
  axis(side = 1 ,at = round(gamma.thres, digits = 3), line = 0.3, tck = 0.02, font = 2)
  dev.off()

  sig.indx <- Dgeo > gamma.thres
  message("calculating p values...\n")
  p.value <- 1 - pgamma(Dgeo, shape = current.a, rate = current.b)


  # return results
  out <- data.frame(Dgeo, p.value, significant = sig.indx)
  rownames(out) <- rownames(geo_coord)
  thres <- data.frame(pvalue = c(0.05,0.01,0.005,0.001),statistic = qgamma(1 - c(0.05,0.01,0.005,0.001), shape = current.a, rate = current.b))
  gamma.par <- c(current.a, current.b)
  names(gamma.par) <- c("alpha", "beta")
  rownames(knn.indx) <- rownames(geo_coord)
  knn.name <- apply(knn.indx,2, function(x){rownames(geo_coord)[x]})
  res.out <- list(out, thres, gamma.par, knn.indx, knn.name, s)
  names(res.out) <- c("statistics","threshold","gamma_parameter", "knn_index", "knn_name", "scalar")


  if(!multi_stages){
    return(res.out)
  }else{
    # multi-stage testing
    message(paste0("Start multi-stage KNN test process with k=",k,
                   " and using the null Gamma distribution with shape=", round(current.a, digits = 3),
                   " and rate=",round(current.b, digits = 3),
                   " (the parameters of Gamma distribution were determined by MLE)\n"))

    i=1
    declining_rate = 0.8
    to_keep <- res.out$statistics$p.value >= p_thres*(1-declining_rate^i)
    tmp.pgdM <- pgdM[to_keep, to_keep]
    tmp.geo_coord <- geo_coord[to_keep,]

    res.Iters <- list(res.out)
    while (i <= maxIter) {
      if(i > 1){
        tmp.pgdM <- tmp.pgdM[to_keep, to_keep]
        tmp.geo_coord <- tmp.geo_coord[to_keep,]
      }

      message(paste0("Iteration ", i,"\r"), appendLF = F)


      # find KNN
      tmp.knn.indx <- find_gen_knn(tmp.pgdM, k=k)
      # KNN predicttion
      tmp.pred.geo_coord <- pred_geo_coord_knn(tmp.geo_coord, tmp.pgdM, tmp.knn.indx)
      # calculate Dg statistic
      tmp.Dgeo <- cal_Dgeo(pred.geo_coord = tmp.pred.geo_coord, geo_coord = tmp.geo_coord, scalar = s)
      tmp.p.value <- 1 - pgamma(tmp.Dgeo, shape = current.a, rate = current.b)
      to_keep <- tmp.p.value >= p_thres*(1-declining_rate^i)


      if(all(tmp.p.value > p_thres)){
        message("\nNo new significant sample is identified. Multi-stage testing process ends...\n")
        break
      }else{
        out <- data.frame(tmp.Dgeo, tmp.p.value, significant = tmp.p.value < p_thres)
        colnames(out) <- c("Dgeo", "p.value", "significant")
        rownames(out) <- rownames(tmp.geo_coord)
        thres <- data.frame(pvalue = c(0.05,0.01,0.005,0.001),statistic = qgamma(1 - c(0.05,0.01,0.005,0.001), shape = current.a, rate = current.b))
        gamma.par <- c(current.a, current.b)
        names(gamma.par) <- c("alpha", "beta")
        rownames(tmp.knn.indx) <- rownames(tmp.geo_coord)
        tmp.knn.name <- apply(tmp.knn.indx,2, function(x){rownames(tmp.geo_coord)[x]})
        tmp.res <- list(out, thres, gamma.par,  tmp.knn.indx, tmp.knn.name, s)
        names(tmp.res) <- c("statistics","threshold","gamma_parameter", "knn_index", "knn_name", "scalar")
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
        #tmp.indx <- which(tmp$statistics$Dgeo > collapse_res$statistics$Dgeo)
        tmp.indx <- which(!is.na(tmp$statistics$Dgeo))
        if(length(tmp.indx)>0){
          collapse_res$statistics[tmp.indx,] <- tmp$statistics[tmp.indx,]
          collapse_res$knn_name[tmp.indx,] <- tmp$knn_name[tmp.indx,]
          for(j in tmp.indx){
            collapse_res$knn_index[j,] <- match(collapse_res$knn_name[j,], rownames(collapse_res$statistics))
          }
        }
      }
    }

    logp.plot <- paste0(plot_dir, "/KNN_GenSP_test_multi_stage_Log10P_comparison.pdf")
    message(paste("the plot of comparing -logP between single-stage and multi-stage KNN tests is saved at ", logp.plot," \n", sep = ""))
    pdf(logp.plot, width = 4, height = 4.2)
    plot(-log10(res.out$statistics$p.value),
         -log10(collapse_res$statistics$p.value),
         xlab = expression("-log"[10]~"(p) of single-stage KNN test"),
         ylab = expression("-log"[10]~"(p) of multi-stage KNN test"),
         main = "KNN in Genetic space")
    dev.off()

    if(keep_all_stg_res){
      names(res.Iters) <- paste0("Iter_", 1:length(res.Iters))
      res.out <- c(res.Iters, collapse_res)
      return(c(res.Iters, combined_result = list(collapse_res)))
    }else{
      return(collapse_res)
    }
  } # multi-stage testing end

} # ggoutlier_geneticspace end
