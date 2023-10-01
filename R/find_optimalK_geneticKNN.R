find_optimalK_geneticKNN <- function(geo_coord,
                                     pgdM,
                                     w_power,
                                     klim,
                                     do_par,
                                     s,
                                     cl
                                     ){
  # automatically select k if k=NULL

  k = min(klim)
  all.D <- c()

  ## parallel computation (if `cpu` > 1)
  if(do_par){
    doParallel::registerDoParallel(cl)
    kindx <- seq(from = min(klim), to = max(klim), by = 1)
    parallel::clusterExport(cl = cl,
                            unclass(lsf.str(envir = asNamespace("GGoutlieR"),
                                            all = T)),
                            envir = as.environment(asNamespace("GGoutlieR"))
    )

    # NOTE: abolish progress bar to meet the requirement of CRAN submission
    ## setup a progress bar for foreach
    #pb <- txtProgressBar(max = max(klim), min = min(klim), style = 3)
    #progress <- function(n){setTxtProgressBar(pb, n)}
    #opts <- list(progress = progress)
    #clusterExport(cl, "opts", envir = environment())
    #all.D <- foreach(k = kindx, .packages=c('geosphere','FastKNN'), .combine="c", .options.snow = opts) %dopar% {

    all.D <- foreach(k = kindx, .packages=c('sf','FastKNN'), .combine="c") %dopar% {
      knn.indx <- find_gen_knn(pgdM, k=k)
      # KNN prediction
      pred.geo_coord <- pred_geo_coord_knn(geo_coord = geo_coord,
                                           pgdM = pgdM,
                                           knn.indx = knn.indx,
                                           w_power = w_power)
      # calculate Dgeo statistic
      return(sum(cal_Dgeo(pred.geo_coord = pred.geo_coord, geo_coord = geo_coord, scalar = s)))
    }
    stopCluster(cl)
  } else {
    ## computation with a single cpu

    ## setup a progress bar
    pb <- txtProgressBar(max = max(klim), min = min(klim), style = 3)

    while(k <= max(klim)){
      # find KNN
      knn.indx <- find_gen_knn(pgdM, k=k)
      # KNN prediction
      pred.geo_coord <- pred_geo_coord_knn(geo_coord = geo_coord,
                                           pgdM = pgdM,
                                           knn.indx = knn.indx,
                                           w_power = w_power)
      # calculate Dg statistic
      Dgeo <- cal_Dgeo(pred.geo_coord = pred.geo_coord, geo_coord = geo_coord, scalar = s)
      all.D <- c(all.D, sum(Dgeo))

      setTxtProgressBar(pb, k)
      k=k+1
    }
  }

  return(all.D)
} # find_optimalK_geneitcKNN end
