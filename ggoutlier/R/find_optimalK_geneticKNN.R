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
    all.D <- foreach(k = kindx, .packages=c('geosphere','FNN'), .combine="c") %dopar% {
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
    while(k <= max(klim)){
      cat(paste0("Calculating Dgeo with k = ",k,"\r"), append = F)
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
      k=k+1
    }
  }

  return(all.D)
} # find_optimalK_geneitcKNN end
