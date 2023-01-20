find_optimalK_geoKNN <- function(geo_coord,
                                 gen_coord,
                                     geo.dM,
                                     w_power,
                                     klim,
                                     do_par,
                                     min_nn_dist,
                                     cl
){
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
    all.D <- foreach(k = kindx, .packages='FNN', .combine="c") %dopar% {
      # find KNN
      knn.indx <- find_geo_knn(geo.dM = geo.dM,
                               k=k,
                               min_nn_dist=min_nn_dist)
      # KNN prediction
      pred.q <- pred_q_knn(geo_coord = geo_coord,
                           gen_coord = gen_coord,
                           geo.dM =  geo.dM,
                           knn.indx = knn.indx,
                           w_power = w_power)
      # calculate Dg statistic
      return(sum(cal_Dg(pred.q = pred.q, gen_coord = gen_coord)))
    }
    stopCluster(cl)
  } else {
    ## computation with a single cpu
    while(k <= max(klim)){
      cat(paste0("Calculating Dg with k = ",k,"\r"), append = F)
      # find KNN
      knn.indx <- find_geo_knn(geo.dM = geo.dM, k=k, min_nn_dist=min_nn_dist)
      # KNN prediction
      pred.q <- pred_q_knn(geo_coord = geo_coord, gen_coord = gen_coord, geo.dM =  geo.dM, knn.indx, w_power = w_power)

      # calculate Dg statistic
      Dg <- cal_Dg(pred.q, gen_coord)
      all.D <- c(all.D, sum(Dg))
      k=k+1
    }
  }
  return(all.D)
} # find_optimalK_geoKNN
