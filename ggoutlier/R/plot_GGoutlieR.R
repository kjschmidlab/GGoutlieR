#' Project GGNet on a geographical map
#' @details
#' Red links on the map denote individual pairs that are genetically similar but geographically remote.
#' The color depth and thickness of red links are proportional to -log10(p) based on the empirical Gamma distribution obtained from `detect_outlier_in_GeneticSpace`.
#' Blue links on the map denote individual pairs that are genetically different but geographically close.
#' The color depth and thickness of blue links are proportional to -log10(p) based on the empirical Gamma distribution obtained from `detect_outlier_in_GeoSpace`
#' @param  geo_coord a two column matrix or data.frame. the first column is longitude and the second one is latitude.
#' @param  GGNet_adjm the output of `get_GGNet_adjacency_matrix`
#' @param  anc_coef a matrix of ancestry coefficients with samples ordered by rows. Ancestry coefficients are used to make pie charts on a geographical map. This argument is optional.
#' @param  pie_color colors of pie charts. colors are automatically assigned if `pie_color = NULL` (which is the default). This argument is optional.
#' @param  p_thres a value of significant level. only samples with p values less than `p_thres` are mapped on a geographical map if `p_thres` is provided (the default is `NULL`). This argument is optional.
#' @param  GeoSP_knn_res an output from `detect_outlier_in_GeoSpace`
#' @param  GenSP_knn_res an output from `detect_outlier_in_GeneticSpace`
#' @param  color_res the resolution of color scale
#' @param  dot_cex the size of dots denoting the positions of samples on a geographical map.
#' @param  map_type the type of plot to draw. it can be `"geographic_knn"`, `"genetic_knn"` and `"both"`.
#' @param  geo_xlim values controlling longitude boundaries to select outliers to present on a geographical map. the default is `geo_xlim = c(-180,180)`.
#' @param  geo_ylim values controlling latitude boundaries to select outliers to present on a geographical map. the default is `geo_ylim = c(-90,90)`.
#' @param  plot_xlim values controlling longitude boundaries of a map.
#' @param  plot_ylim values controlling latitude boundaries of a map.
#' @param  only_edges_in_xylim logic. only the edges with starting points within the given `geo_xlim` and `geo_ylim` will display on a geographical map. the default is `TRUE`.
#' @param  pie_r_scale a scale controlling the radius of pie charts
#' @param  red_alpha a value controlling the transparency of red lines. the default is 0.8
#' @param  map_resolution the resolution of the geographical map. See details in the manual of `rworldmap::getMap()`
#' @param  show_knn_pie logic. If `TRUE`, the ancestry coefficients of K nearest neighbors of significant samples will display on the map. The default is `FALSE`.
#' @param  which_sample a string vector. If users want to only
#' @export

#------------------------------------------------------------------
# status: finished (2022-05-03)
#------------------------------------------------------------------

plot_GGNet_map <- function(geo_coord,
                           GGNet_adjm,
                           anc_coef = NULL,
                           pie_color = NULL,
                           p_thres = NULL,
                           GeoSP_knn_res = NULL,
                           GenSP_knn_res = NULL,
                           color_res = 10,
                           dot_cex = 0.4,
                           map_type = c("geographic_knn", "genetic_knn", "both"),
                           geo_xlim = c(-180,180),
                           geo_ylim = c(-90,90),
                           plot_xlim = NULL,
                           plot_ylim = NULL,
                           only_edges_in_xylim = TRUE,
                           pie_r_scale = 1,
                           red_alpha = 0.8,
                           map_resolution = "high",
                           show_knn_pie = FALSE,
                           which_sample = NULL

){
  require(rworldmap)
  require(scales)
  require(plotrix)
  require(mapplots)
  require(RColorBrewer)
  require(rworldxtra)
  # check input
  map_type = match.arg(map_type)
  if(!is.null(anc_coef)|!is.null(p_thres)){
    if(map_type == "both" | map_type == "geographic_knn"){
      if(is.null(GeoSP_knn_res)){
        stop("As `plot_GGNet_map` is designed to map only outliers on a geographical map when `p_thres`` is given, please provide `GeoSP_knn_res` (an output of `detect_outlier_in_GeoSpace`)")
      }else{
        if(nrow(GeoSP_knn_res$statistics) != nrow(geo_coord)){
          stop("the sample size in `GeoSP_knn_res` and `geo_coord` does not match!")
        }
      }
    }
    if(map_type == "both" | map_type == "genetic_knn"){
      if(is.null(GenSP_knn_res)){
        stop("As `plot_GGNet_map` is designed to map only outliers on a geographical map when `p_thres`` is given, please provide `GenSP_knn_res` (an output of `detect_outlier_in_GeneticSpace`) if you set `map_type = 'both'` or `map_type = 'GenSP'`")
      }
      if(nrow(GenSP_knn_res$statistics) != nrow(geo_coord)){
        stop("the sample size in `GenSP_knn_res` and `geo_coord` does not match!")
      }
    }

    # assign colors if pie_color is NULL
    if(!is.null(anc_coef) & is.null(pie_color)){
      # check if sample sizes match
      if(nrow(anc_coef) != nrow(geo_coord)){stop("the sample size in `anc_coef` and `geo_coord` does not match!")}
      if(ncol(anc_coef) <= 12){
        pie_color = brewer.pal(ncol(anc_coef), "Set3")
      }else{
        # color generation: https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        # shuffle colors
        pie_color = c(matrix(c(matrix(col_vector[1:72], ncol = 8, byrow = T))[c(matrix(c(1:36, rev(37:72)), nrow = 3))], ncol = 4,byrow = T))
        pie_color = pie_color[1:ncol(anc_coef)]
      }
    }
    if(!is.null(anc_coef)){message("ancestry coefficients (`anc_coef`) are mapped as pie charts")}
  }
  if(!is.null(anc_coef) & is.null(p_thres)){
    warning("`p_thres` is `NULL`. ancestry coefficients will NOT be projected to your geographical map. the pie charts of outliers could be added to the map if you give `p_thres`.")
  }


  geomap <- getMap(resolution = map_resolution)

  # get edge colors
  edge.col <- get_GGNet_edge_col(GGNet_adjm, color_res = color_res)

  #******************************
  # edges of geographical space (GeoSP) KNN outlier likelihood
  if(map_type %in% c("both", "geographic_knn")){
    geosp.indx <- which(GGNet_adjm$GeoSP_pvalue > 0, arr.ind = T)
    geosp.pvalue <- GGNet_adjm$GeoSP_pvalue[GGNet_adjm$GeoSP_pvalue > 0]
    geosp.col <- edge.col$GeoSP_col[GGNet_adjm$GeoSP_pvalue > 0]
    if(!is.null(p_thres)){
      sig.tag <-geosp.indx[,1] %in% which(GeoSP_knn_res$statistics$p.value < p_thres)
      geosp.df <- data.frame(geosp.indx, geosp.pvalue, geosp.col, sig.tag, stringsAsFactors=FALSE)
    }else{
      geosp.df <- data.frame(geosp.indx, geosp.pvalue, geosp.col, stringsAsFactors=FALSE)
    }
    ## sort geosp.df -> make deeper colors on the top of the figure
    geosp.df <- geosp.df[order(geosp.df$geosp.pvalue, decreasing = T),]

    ## keep specific samples if which_sample is not NULL
    if(!is.null(which_sample) & all(nrow(geosp.df) > 0)){
      tmp <- apply(sapply(which_sample, function(x){
        grepl(pattern = x, rownames(geosp.df))
      }), 1, any)
      geosp.df <- geosp.df[tmp,]
    }
    ## remove edges out of the given margins
    if(only_edges_in_xylim){
      geo_xlim_max <- max(geo_xlim)
      geo_xlim_min <- min(geo_xlim)
      geo_ylim_max <- max(geo_ylim)
      geo_ylim_min <- min(geo_ylim)
      tmp1 <- apply(geo_coord[geosp.df$row,],1, function(x){x[1] < geo_xlim_max & x[1] > geo_xlim_min & x[2] < geo_ylim_max & x[2] > geo_ylim_min})
      #tmp2 <- apply(geo_coord[geosp.df$col,],1, function(x){x[1] < geo_xlim_max & x[1] > geo_xlim_min & x[2] < geo_ylim_max & x[2] > geo_ylim_min})

      if(!is.null(p_thres)){
        tmp3 <- geosp.df$row %in% which(GeoSP_knn_res$statistics$p.value < p_thres)
        to_keep <- mapply(a = tmp1, b = tmp3, function(a,b,c){return(a&b)}) # if starting points are out of the given range and also not outliers -> remove
      }else{
        #to_keep <- apply(cbind(tmp1, tmp2), 1, any)
        to_keep <- tmp1
      }
      geosp.df <- geosp.df[which(to_keep),]
    }
  }

  ## set plot_xlim and plot_ylim if NULL
  if(is.null(plot_xlim)){plot_xlim <- geo_xlim}
  if(is.null(plot_ylim)){plot_ylim <- geo_ylim}
  par(mar = c(1,1,1,1))
  #******************************
  # edges of genetic space (GenSP) KNN outlier likelihood
  if(map_type %in% c("both", "genetic_knn")){
    gensp.indx <- which(GGNet_adjm$GenSP_pvalue > 0, arr.ind = T)
    gensp.pvalue <- GGNet_adjm$GenSP_pvalue[GGNet_adjm$GenSP_pvalue > 0]
    gensp.col <- edge.col$GenSP_col[GGNet_adjm$GenSP_pvalue > 0]

    if(!is.null(p_thres)){
      # sig.tag shows the edges attach to the outliers
      sig.tag <-gensp.indx[,1] %in% which(GenSP_knn_res$statistics$p.value < p_thres)
      gensp.df <- data.frame(gensp.indx, gensp.pvalue, gensp.col, sig.tag, stringsAsFactors=FALSE)
    }else{
      gensp.df <- data.frame(gensp.indx, gensp.pvalue, gensp.col, stringsAsFactors=FALSE)
    }
    ## sort gensp.df -> make deeper colors on the top of the figure
    gensp.df <- gensp.df[order(gensp.df$gensp.pvalue, decreasing = T),]

    ## keep specific samples if which_sample is not NULL
    if(!is.null(which_sample) & all(nrow(gensp.df) > 0)){
      tmp <- apply(sapply(which_sample, function(x){
        grepl(pattern = x, rownames(gensp.df))
      }), 1, any)
      gensp.df <- gensp.df[tmp,]
    }

    ## remove edges out of the given margins
    if(only_edges_in_xylim){
      geo_xlim_max <- max(geo_xlim)
      geo_xlim_min <- min(geo_xlim)
      geo_ylim_max <- max(geo_ylim)
      geo_ylim_min <- min(geo_ylim)
      tmp1 <- apply(geo_coord[gensp.df$row,],1, function(x){x[1] < geo_xlim_max & x[1] > geo_xlim_min & x[2] < geo_ylim_max & x[2] > geo_ylim_min})
      #tmp2 <- apply(geo_coord[gensp.df$col,],1, function(x){x[1] < geo_xlim_max & x[1] > geo_xlim_min & x[2] < geo_ylim_max & x[2] > geo_ylim_min})
      if(!is.null(p_thres)){
        tmp3 <- gensp.df$row %in% which(GenSP_knn_res$statistics$p.value < p_thres)
        to_keep <- mapply(a = tmp1, b = tmp3, function(a,b){return(a&b)}) # if starting points are out of the given range and also not outliers -> remove
      }else{
        #to_keep <- apply(cbind(tmp1, tmp2), 1, any)
        to_keep <- tmp1
      }
      gensp.df <- gensp.df[which(to_keep),]
    }
  }


  plot(geomap, xlim = plot_xlim, ylim = plot_ylim, col = NA, border = "white")
  if(map_type == "both" | map_type == "geographic_knn"){
    if(is.null(p_thres)){
      segments(x0 = geo_coord[geosp.df$row,1],
               y0 = geo_coord[geosp.df$row,2],
               x1 = geo_coord[geosp.df$col,1],
               y1 = geo_coord[geosp.df$col,2],
               col = alpha(geosp.df$geosp.col, alpha = 0.9),
               lwd = -log10(geosp.df$geosp.pvalue)
      )
    }else{
      # sig.row corresponds to the rows (the samples under testing)
      # sig.col corresponds to the columns (the K nearest neighbors of the samples under testing)
      sig.row <- geosp.df$row[geosp.df$sig.tag]
      sig.col <- geosp.df$col[geosp.df$sig.tag]
      sig.color <- geosp.df$geosp.col[geosp.df$sig.tag]
      sig.likelihood <- geosp.df$geosp.pvalue[geosp.df$sig.tag]
      segments(x0 = geo_coord[sig.row,1],
               y0 = geo_coord[sig.row,2],
               x1 = geo_coord[sig.col,1],
               y1 = geo_coord[sig.col,2],
               col = alpha(sig.color, alpha = 0.9),
               lwd = -log10(sig.likelihood)
      )
    }
  }
  if(map_type == "both" | map_type == "genetic_knn"){
    # adjust transparency according to the type of figures
    if(map_type == "both"){
      tmp.alpha = red_alpha * 0.75
    }else{
      tmp.alpha = red_alpha
    }
    if(is.null(p_thres)){
      segments(x0 = geo_coord[gensp.df$row,1],
               y0 = geo_coord[gensp.df$row,2],
               x1 = geo_coord[gensp.df$col,1],
               y1 = geo_coord[gensp.df$col,2],
               col = alpha(gensp.df$gensp.col, alpha = tmp.alpha),
               lwd = -log10(gensp.df$gensp.pvalue)
      )
    }else{
      # sig.row corresponds to the rows (the samples under testing)
      # sig.col corresponds to the columns (the K nearest neighbors of the samples under testing)
      sig.row <- gensp.df$row[gensp.df$sig.tag]
      sig.col <- gensp.df$col[gensp.df$sig.tag]
      sig.color <- gensp.df$gensp.col[gensp.df$sig.tag]
      sig.likelihood <- gensp.df$gensp.pvalue[gensp.df$sig.tag]
      segments(x0 = geo_coord[sig.row,1],
               y0 = geo_coord[sig.row,2],
               x1 = geo_coord[sig.col,1],
               y1 = geo_coord[sig.col,2],
               col = alpha(sig.color, alpha = tmp.alpha),
               lwd = -log10(sig.likelihood)
      )
    }
  }

  plot(geomap, xlim = plot_xlim, ylim = plot_ylim, add=T)
  if(is.null(anc_coef)){points(x = geo_coord[,1], y = geo_coord[,2], pch = 16, cex = dot_cex)}
  if(!is.null(anc_coef)){
    points(x = geo_coord[,1], y = geo_coord[,2], pch = 16, cex = dot_cex)
    #}else{
    if(!is.null(p_thres) & map_type %in% c("both","geographic_knn")){
      sig.row <- unique(geosp.df$row[geosp.df$sig.tag])
      if(length(sig.row) == 0){
          warning("there is no outlier in the region selected according to `geo_xlim` and `geo_ylim`")
      }else{
        pie.r = abs(diff(par("usr")[1:2]) )*0.005
        for(i in 1:length(sig.row)){
          add.pie(z = round(anc_coef[sig.row[i],]*10^5),
                  x = geo_coord[sig.row[i],1],
                  y = geo_coord[sig.row[i],2],
                  col = pie_color,
                  labels = NA,
                  radius = pie.r*pie_r_scale)
        }
        # add pie charts of KNNs if show_knn_pie is TRUE
        if(show_knn_pie){
          for(i in 1:length(sig.row)){
            for(j in GeoSP_knn_res$knn_index[sig.row[i],]){
              add.pie(z = round(anc_coef[j,]*10^5),
                      x = geo_coord[j,1],
                      y = geo_coord[j,2],
                      col = pie_color,
                      labels = NA,
                      radius = pie.r*pie_r_scale)
            }
          }
        }
      }
    }
    if(!is.null(p_thres)  & map_type %in% c("both","genetic_knn")){
      sig.row <- unique(gensp.df$row[gensp.df$sig.tag])
      if(length(sig.row) == 0){
        warning("there is no outlier in the region selected according to `geo_xlim` and `geo_ylim`")
      }else{
        pie.r = abs(diff(par("usr")[1:2]) )*0.005
        for(i in 1:length(sig.row)){
          add.pie(z = round(anc_coef[sig.row[i],]*10^5),
                  x = geo_coord[sig.row[i],1],
                  y = geo_coord[sig.row[i],2],
                  col = pie_color,
                  labels = NA,
                  radius = pie.r*pie_r_scale)
        }
        # add pie charts of KNNs if show_knn_pie is TRUE
        if(show_knn_pie){
          for(i in 1:length(sig.row)){
            for(j in GenSP_knn_res$knn_index[sig.row[i],]){
              add.pie(z = round(anc_coef[j,]*10^5),
                      x = geo_coord[j,1],
                      y = geo_coord[j,2],
                      col = pie_color,
                      labels = NA,
                      radius = pie.r*pie_r_scale)
            }
          }
        }
      }
    }
  }
  if(map_type == "genetic_knn"){
    color.legend( xl =plot_xlim[2] - abs(plot_xlim[2] - plot_xlim[1])*0.15,
                  xr = plot_xlim[2],
                  yb = plot_ylim[2] - abs(plot_ylim[2] - plot_ylim[1])*0.07,
                  yt = plot_ylim[2] - abs(plot_ylim[2] - plot_ylim[1])*0.05 , # the coordinates
                  legend = c(0, round(max(-log10(gensp.pvalue)))) ,
                  gradient="x",
                  rect.col=edge.col$GenSP_colkey, align="rb")
    text(x = plot_xlim[2] - abs(plot_xlim[2] - plot_xlim[1])*0.075,
         y = plot_ylim[2],
         labels = expression(-log[10](p)), font = 2, cex = 1.1
    )
  }

  if(map_type == "geographic_knn"){
    color.legend( xl =plot_xlim[2] - abs(plot_xlim[2] - plot_xlim[1])*0.15,
                  xr = plot_xlim[2],
                  yb = plot_ylim[2] - abs(plot_ylim[2] - plot_ylim[1])*0.07,
                  yt = plot_ylim[2] - abs(plot_ylim[2] - plot_ylim[1])*0.05 , # the coordinates
                  legend = c(0, round(max(-log10(geosp.pvalue)))) ,
                  gradient="x",
                  rect.col=edge.col$GeoSP_colkey, align="rb")
    text(x = plot_xlim[2] - abs(plot_xlim[2] - plot_xlim[1])*0.075,
         y = plot_ylim[2],
         labels = expression(-log[10](p)), font = 2, cex = 1.1
    )
  }
  if(map_type == "both"){
    par(xpd = NA)
    color.legend( xl =plot_xlim[2] - abs(plot_xlim[2] - plot_xlim[1])*0.15,
                  xr = plot_xlim[2],
                  yb = plot_ylim[2] - abs(plot_ylim[2] - plot_ylim[1])*0.01,
                  yt = plot_ylim[2] - abs(plot_ylim[2] - plot_ylim[1])*(-0.02) , # the coordinates
                  legend = c(0, round(max(-log10(geosp.pvalue)))) ,
                  gradient="x",
                  rect.col=edge.col$GeoSP_colkey, align="rb")
    color.legend( xl =plot_xlim[2] - abs(plot_xlim[2] - plot_xlim[1])*0.15,
                  xr = plot_xlim[2],
                  yb = plot_ylim[2] - abs(plot_ylim[2] - plot_ylim[1])*0.10,
                  yt = plot_ylim[2] - abs(plot_ylim[2] - plot_ylim[1])*0.07, # the coordinates
                  legend = c(0, round(max(-log10(gensp.pvalue)))) ,
                  gradient="x",
                  rect.col=edge.col$GenSP_colkey, align="rb")

    text(x = plot_xlim[2] - abs(plot_xlim[2] - plot_xlim[1])*0.075,
         y = plot_ylim[2] + abs(plot_ylim[2] - plot_ylim[1])*0.1,
         labels = expression(-log[10](p)), font = 2, cex = 1.1
    )
    par(xpd = FALSE)
  }
} # plot_GGNet_map end


#-------------------------------------------------------------
# assign colors to edges
# `get_GGNet_edge_col` is a function used in the `plot_GGNet_map`
#-------------------------------------------------------------
# Argument:
# GGNet_adjm: an output of `get_GGNet_adjacency_matrix
# color_res: the resolution of color scale
get_GGNet_edge_col <- function(GGNet_adjm, color_res = 10){
  require(dichromat)
  GeoSP.colfun <- colorRampPalette(c("white", "dodgerblue"))
  GeoSP.col <- GeoSP.colfun(color_res)
  GeoSP_col <- rep(NA, length(GGNet_adjm$GeoSP_pvalue))
  geosp.indx <- c(GGNet_adjm$GeoSP_pvalue) > 0
  GeoSP_col[geosp.indx] <- as.character(GeoSP.col[as.numeric(cut(-log10(GGNet_adjm$GeoSP_pvalue[geosp.indx]), breaks = color_res))])

  GenSP.colfun <- colorRampPalette(c("white", "orangered"))
  GenSP.col <- GenSP.colfun(color_res)
  GenSP_col <- rep(NA, length(GGNet_adjm$GenSP_pvalue))
  gensp.indx <- c(GGNet_adjm$GenSP_pvalue) > 0
  GenSP_col[gensp.indx] <- as.character(GenSP.col[as.numeric(cut(-log10(GGNet_adjm$GenSP_pvalue[gensp.indx]), breaks = color_res))])

  return(list(GeoSP_col = c(GeoSP_col),
              GenSP_col = c(GenSP_col),
              GeoSP_colkey = c(GeoSP.col),
              GenSP_colkey = c(GenSP.col)
  ))
} # get_GGNet_edge_col end
