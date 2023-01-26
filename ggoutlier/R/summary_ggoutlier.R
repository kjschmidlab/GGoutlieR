#'
#' @details
#' Get a summary table from the `ggoutlier` output
#' @param  ggoutlier_res output from the function `ggoutlier`

#' @export

summary_ggoutlier <- function(ggoutlier_res){
  method <- attributes(ggoutlier_res)$model
  if(method == "composite"){
    ## get significant outliers from genetic space KNN result
    sig_index <- ggoutlier_res$geneticKNN_result$statistics$significant
    sig_geneticKNN <- data.frame(
      "ID" = rownames(ggoutlier_res$geneticKNN_result$statistics)[sig_index],
      "method" = "geneticKNN",
      "p.value" = ggoutlier_res$geneticKNN_result$statistics[sig_index,"p.value"])
    sig_geneticKNN <- sig_geneticKNN[order(sig_geneticKNN$p.value, decreasing = F),]

    ## get significant outliers from geographical space KNN result
    sig_index <- ggoutlier_res$geoKNN_result$statistics$significant
    sig_geoKNN <- data.frame(
      "ID" = rownames(ggoutlier_res$geoKNN_result$statistics)[sig_index],
      "method" = "geoKNN",
      "p.value" = ggoutlier_res$geoKNN_result$statistics[sig_index,"p.value"])
    sig_geoKNN <- sig_geoKNN[order(sig_geoKNN$p.value, decreasing = F),]

    prop_sig <- c(mean(ggoutlier_res$geneticKNN_result$statistics$significant),
                  mean(ggoutlier_res$geoKNN_result$statistics$significant))
    names(prop_sig) <- c("geneticKNN", "geoKNN")
    summary_out <- rbind.data.frame(sig_geneticKNN, sig_geoKNN)
    summary_out <- summary_out[order(summary_out$p.value, decreasing = F),]
  }
  if(method == "ggoutlier_geoKNN"){
    ## get significant outliers from genetic space KNN result
    sig_index <- ggoutlier_res$statistics$significant
    sig_geneticKNN <- data.frame(
      "ID" = rownames(ggoutlier_res$statistics)[sig_index],
      "p.value" = ggoutlier_res$statistics[sig_index,"p.value"])
    sig_geneticKNN <- sig_geneticKNN[order(sig_geneticKNN$p.value, decreasing = F),]
    prop_sig <- c(mean(ggoutlier_res$geneticKNN_result$statistics$significant),
                  mean(ggoutlier_res$geoKNN_result$statistics$significant))
    names(prop_sig) <- c("geneticKNN")
    summary_out <- sig_geneticKNN
  }
  if(method == "ggoutlier_geoKNN"){
    ## get significant outliers from geographical space KNN result
    sig_index <- ggoutlier_res$geoKNN_result$statistics$significant
    sig_geoKNN <- data.frame(
      "ID" = rownames(ggoutlier_res$geoKNN_result$statistics)[sig_index],
      "p.value" = ggoutlier_res$geoKNN_result$statistics[sig_index,"p.value"])
    sig_geoKNN <- sig_geoKNN[order(sig_geoKNN$p.value, decreasing = F),]

    prop_sig <- c(mean(ggoutlier_res$geoKNN_result$statistics$significant))
    names(prop_sig) <- c("geoKNN")
    summary_out <- sig_geoKNN
  }
  rownames(summary_out) <- NULL
  sig_ID_list <- unique(summary_out$ID)
  attr(summary_out, "proportion_of_outliers") = prop_sig
  attr(summary_out, "significant_samples" = sig_ID_list)
  class(summary_out) <- "ggoutlier_summary"
  return(summary_out)
} # summary_ggoutlier end
