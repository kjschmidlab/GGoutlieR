#' Visualize GGoutlieR results on a geographical map
#' @description Visualize geo-genetic patterns of outliers with their K nearest neighbors
#' @details
#' Red links on the map denote individual pairs that are genetically similar but geographically remote.
#' The color depth and thickness of red links are proportional to -log10(p) based on the empirical Gamma distribution obtained from `detect_outlier_in_GeneticSpace`.
#' Blue links on the map denote individual pairs that are genetically different but geographically close.
#' The color depth and thickness of blue links are proportional to -log10(p) based on the empirical Gamma distribution obtained from `detect_outlier_in_GeoSpace`
#' @param  ggoutlier_res output of `ggoutlier`
#' @param  geo_coord matrix or data.frame with two columns. The first column is longitude and the second one is latitude.
#' @param  anc_coef matrix. A matrix of ancestry coefficients with samples ordered by rows. Ancestry coefficients are used to make pie charts on a geographical map. This argument is optional.
#' @param  gen_coord matrix. A matrix of "coordinates in a genetic space". It should be identical to the `gen_coord` used for running `ggoutlier`
#' @param  pie_color string. Colors of pie charts. colors are automatically assigned if `pie_color = NULL` (which is the default). This argument is optional.
#' @param  map_color string. Colors of map contours. The default is `map_color = "black"`
#' @param  p_thres numeric. A value of significant level. Only outliers (p values less than `p_thres`) are mapped on a geographical map if `p_thres` is provided (the default is `NULL`). This argument is optional.
#' @param  color_res integer. The resolution of color scale.
#' @param  dot_cex numeric. The size of dots denoting the positions of samples on a geographical map.
#' @param  map_type string. The type of plot to draw. It can be `"geographic_knn"`, `"genetic_knn"` and `"both"`.
#' @param  select_xlim vector. Values controlling longitude boundaries of a window to select outliers to present on a geographical map. The default is `select_xlim = c(-180,180)`.
#' @param  select_ylim vector. Values controlling latitude boundaries of a window to select outliers to present on a geographical map. The default is `select_ylim = c(-90,90)`.
#' @param  plot_xlim vector. Values controlling longitude boundaries of a map.
#' @param  plot_ylim vector. Values controlling latitude boundaries of a map.
#' @param  only_edges_in_xylim logic. only the edges with starting points within the given `select_xlim` and `select_ylim` will display on a geographical map. If `FALSE`, the edges out of the given boundaries will be removed from your plot. The default is `TRUE`.
#' @param  pie_r_scale numeric. A scale controlling the radius of pie charts
#' @param  map_resolution a character string. The resolution of the geographical map. See details of the `scale` argument in the manual of `rnaturalearth::ne_countries()`. The default is `map_resolution = "medium"`
#' @param  show_knn_pie logic. If `TRUE`, the ancestry coefficients of K nearest neighbors of significant samples will display on the map. The default is `FALSE`.
#' @param  which_sample a string vector of sample ID(s). If users want to only show specific sample(s)
#' @param  add_benchmark_graph logic. If `TRUE`, a benchmark graph with only pie charts of ancestry coefficients for comparison with the outlier graph.
#' @param  adjust_p_value_projection logic. If `TRUE`, the function will perform KNN prediction by forcing K=1 and compute new p-values for visualization.
#' @param  linewidth_range numeric. A vector of two values. It is used to control the minimal and maximal width of KNN network on the geographical map.
#' @param  plot_labels character. A string of labels for plots. The default is `plot_labels = "auto"`. This parameter is for `cowplot::plot_grid`.
#'
#' @returns ggplot object. The plot is geographical map(s) with colored lines showing sample pairs with unusual geo-genetic associations.
#' @examples
#' library(GGoutlieR)
#' data("ipk_anc_coef") # get ancestry coefficients
#' data("ipk_geo_coord") # get geographical coordinates
#' data(ggoutlier_example) # get an example output of ggoutlier
#' \dontrun{
#' plot_ggoutlier(ggoutlier_res = ggoutlier_example,
#'                gen_coord = ipk_anc_coef,
#'                geo_coord = ipk_geo_coord,
#'                p_thres = 0.025,
#'                map_type = "both",
#'                select_xlim = c(-20,140),
#'                select_ylim = c(10,62),
#'                plot_xlim = c(-20,140),
#'                plot_ylim = c(10,62),
#'                pie_r_scale = 1.8,
#'                map_resolution = "medium")
#'                }
#'
#' @export

#------------------------------------------------------------------
# status: finished (need to update user manual)
#------------------------------------------------------------------

plot_ggoutlier <- function(ggoutlier_res,
                           geo_coord,
                           anc_coef = NULL,
                           gen_coord = NULL,
                           pie_color = NULL,
                           map_color = "black",
                           p_thres = NULL,
                           color_res = 10,
                           dot_cex = NULL,
                           map_type = c("geographic_knn", "genetic_knn", "both"),
                           select_xlim = c(-180,180),
                           select_ylim = c(-90,90),
                           plot_xlim = NULL,
                           plot_ylim = NULL,
                           only_edges_in_xylim = TRUE,
                           pie_r_scale = 1,
                           map_resolution = "medium",
                           show_knn_pie = FALSE,
                           which_sample = NULL,
                           add_benchmark_graph = TRUE,
                           adjust_p_value_projection = FALSE,
                           linewidth_range = c(0.5,3),
                           plot_labels = "auto"


){
  required_pkgs <- c(#"rworldmap",
                     #"mapplots",
                     #"rworldxtra","dichromat", "sp",
                     "RColorBrewer","scales",
                     "rnaturalearth",
                     "rnaturalearthdata", "sf",
                     "ggplot2", "cowplot",
                     "ggforce", "rlang", "ggfun",
                     "stats", "tidyr", "dplyr", "utils")
  invisible(lapply(required_pkgs, FUN=function(x){suppressMessages(library(x, verbose = FALSE, character.only = TRUE))}))

  # use on.exit to prevent changes in users' pars
  oldpar <- par(no.readonly = TRUE) # save original par of users
  on.exit(oldpar)

  # extract data
  if(attributes(ggoutlier_res)$model == "composite"){
    GeoSP_knn_res <- ggoutlier_res$geoKNN_result
    GenSP_knn_res <- ggoutlier_res$geneticKNN_result
  }
  if(attributes(ggoutlier_res)$model == "ggoutlier_geoKNN"){
    GeoSP_knn_res <- ggoutlier_res
    GenSP_knn_res <- NULL
  }
  if(attributes(ggoutlier_res)$model == "ggoutlier_geneticKNN"){
    GeoSP_knn_res <- NULL
    GenSP_knn_res <- ggoutlier_res
  }

  # check input
  if(is.null(anc_coef) & is.null(gen_coord)){
    stop("please provide at least either `anc_coef` or `gen_coord`")
  }else{
    if(is.null(gen_coord)){gen_coord <- anc_coef}
    if(is.null(anc_coef) & !is.null(gen_coord)){
      if(all(abs(apply(gen_coord,1 , sum) - 1) < 10^-5)){ # double check if gen_coord is as expectation of ancestry coefficients
        message("assign `gen_coord` as ancestry coefficients (`anc_coef`)")
        anc_coef <- gen_coord
      }else{
        message("\n ancestry coefficients are missing.\n")
      }
    } # if is.null(anc_coef) end
  }

  map_type = match.arg(map_type)
  if(!is.null(anc_coef) & !is.null(p_thres)){
    if(map_type == "both" | map_type == "geographic_knn"){
      if(is.null(GeoSP_knn_res)){
        stop("`plot_ggoutlier` is designed to map outliers on a geographical map, please provide `GeoSP_knn_res` (an output of `detect_outlier_in_GeoSpace`) if you set `map_type = both` or `map_type = geographic_knn`")
      }else{
        if(nrow(GeoSP_knn_res$statistics) != nrow(geo_coord)){
          stop("the sample size in `GeoSP_knn_res` and `geo_coord` does not match!")
        }
      }
    }
    if(map_type == "both" | map_type == "genetic_knn"){
      if(is.null(GenSP_knn_res)){
        stop("`plot_ggoutlier` is designed to map outliers on a geographical map, please provide `GenSP_knn_res` (an output of `detect_outlier_in_GeneticSpace`) if you set `map_type = 'both'` or `map_type = 'genetic_knn'`")
      }
      if(nrow(GenSP_knn_res$statistics) != nrow(geo_coord)){
        stop("the sample size in `GenSP_knn_res` and `geo_coord` does not match!")
      }
    }
    if(is.null(dot_cex)){
      dot_cex <- min(linewidth_range)
    }

    # assign colors if pie_color is NULL
    if(!is.null(anc_coef) & is.null(pie_color)){
      # check if sample sizes match
      if(nrow(anc_coef) != nrow(geo_coord)){stop("the sample size in `anc_coef` and `geo_coord` does not match!")}
      if(ncol(anc_coef) <= 12){
        pie_color = RColorBrewer::brewer.pal(ncol(anc_coef), "Set3")
      }else{
        # color generation: https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
        qual_col_pals = RColorBrewer::brewer.pal.info[brewer.pal.info$category == 'qual',]
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



  ## geographical data: data.frame to sf format
  geo_data_df <- data.frame(ID = rownames(geo_coord),
                            x = geo_coord$x,
                            y = geo_coord$y)

  geo_data_sf <- sf::st_as_sf(geo_data_df,
                              coords = c("x", "y"), crs = 4326)
  # get world geographical map
  geomap <- rnaturalearth::ne_countries(scale = map_resolution, returnclass = "sf")
  #geomap <- getMap(resolution = map_resolution)

  # get an adjancency matrix 'GGNet_adjm'
  if(all(dim(gen_coord) == dim(anc_coef))){
    if(all(round(gen_coord, digits = 5) == round(anc_coef, digits = 5))){
      gen_coord_eq_AncCoeff = TRUE
    }else{gen_coord_eq_AncCoeff = FALSE}
  }

  GGNet_adjm <- get_GGNet_adjacency_matrix(ggoutlier_res = ggoutlier_res,
                                           gen_coord = gen_coord,
                                           geo_coord = geo_coord,
                                           mutual = FALSE,
                                           adjust_p_value = adjust_p_value_projection
                                           )

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
      select_xlim_max <- max(select_xlim)
      select_xlim_min <- min(select_xlim)
      select_ylim_max <- max(select_ylim)
      select_ylim_min <- min(select_ylim)
      tmp1 <- apply(geo_coord[geosp.df$row,],1, function(x){x[1] < select_xlim_max & x[1] > select_xlim_min & x[2] < select_ylim_max & x[2] > select_ylim_min})
      #tmp2 <- apply(geo_coord[geosp.df$col,],1, function(x){x[1] < select_xlim_max & x[1] > select_xlim_min & x[2] < select_ylim_max & x[2] > select_ylim_min})

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
  if(is.null(plot_xlim)){plot_xlim <- select_xlim}
  if(is.null(plot_ylim)){plot_ylim <- select_ylim}
  #par(mar = c(1,1,1,1))
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
      select_xlim_max <- max(select_xlim)
      select_xlim_min <- min(select_xlim)
      select_ylim_max <- max(select_ylim)
      select_ylim_min <- min(select_ylim)
      tmp1 <- apply(geo_coord[gensp.df$row,],1, function(x){x[1] < select_xlim_max & x[1] > select_xlim_min & x[2] < select_ylim_max & x[2] > select_ylim_min})
      #tmp2 <- apply(geo_coord[gensp.df$col,],1, function(x){x[1] < select_xlim_max & x[1] > select_xlim_min & x[2] < select_ylim_max & x[2] > select_ylim_min})
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
  if(!is.null(anc_coef)){
    ## prepare a data frame for making pie charts on ggplot
    anc_coef_for_pie <- data.frame(geo_coord,
                                   pie_r_scale = pie_r_scale,
                                   anc_coef
    )
    cols <- colnames(anc_coef_for_pie)[4:ncol(anc_coef_for_pie)]
    cols2 <- enquo(cols)
    anc_ggplot_df <- gather(anc_coef_for_pie, "type", "value", !!cols2)
    anc_ggplot_df$type <- factor(anc_ggplot_df$type, levels = cols) # set legend order based on order of "cols"
  }

  if(add_benchmark_graph & !is.null(anc_coef)){
    #if(vertical_plots){
      #par(mfrow = c(2,1), mar = c(1,1,1,1))
    #}else{
      #par(mfrow = c(1,2), mar = c(1,1,1,1))
    #}


    ## make benchmark plot
    geomap_plot_benchmark <-
      ggplot2::ggplot(data = geomap) +
      ggplot2::geom_sf(color = map_color, show.legend = NA) +
      coord_sf(xlim = plot_xlim, ylim = plot_ylim, expand = TRUE) +

      scale_fill_manual(values=pie_color) +
      geom_arc_bar(aes(x0 = x, y0 = y,r0=0,
                       fill = type, r=pie_r_scale,
                       amount=value),
                   data=anc_ggplot_df, stat='pie',
                   col = "black",
                   show.legend = FALSE) +
      theme(axis.title.x = element_text(color = "white"),
            axis.title.y = element_text(color = "white"),
            # hide legend but keep blank space
            legend.text = element_text(color = "white"),
            legend.title = element_text(color = "white"),
            legend.key = element_rect(fill = NA)
            ) +
      guides(color = guide_legend(override.aes =
                                    list(fill = NA,
                                         color = c("white"),
                                         alpha = 0)))
      #print(geomap_plot_benchmark)


    #save.image("tmp_test_sf.Rdata")


    #sp::plot(geomap, xlim = plot_xlim, ylim = plot_ylim, border = map_color)
    #pie.r = abs(diff(par("usr")[1:2]) )*0.005
    #for(i in 1:nrow(anc_coef)){
    #  add.pie(z = round(anc_coef[i,]*10^5),
    #          x = geo_coord[i,1],
    #          y = geo_coord[i,2],
    #          col = pie_color,
    #          labels = NA,
    #          radius = pie.r*pie_r_scale)
    #}
  } # make benchmark plot end

  ##---------------------------------------------------------
  ## prepare ggplot input for the GGoutlieR network plot(s)
  ##---------------------------------------------------------

  #sp::plot(geomap, xlim = plot_xlim, ylim = plot_ylim, col = NA, border = "white")
  if(map_type == "both" | map_type == "geographic_knn"){
    if(is.null(p_thres)){

      ## create a data.frame for drawing segments on ggplot
      ## x0 and y0 correspond to the rows (the samples under testing)
      ## x1 and y1 correspond to the columns (the K nearest neighbors of the samples under testing)
      link_geosp_df <-
        data.frame(
          x0 = geo_coord[geosp.df$row,1],
          y0 = geo_coord[geosp.df$row,2],
          x1 = geo_coord[geosp.df$col,1],
          y1 = geo_coord[geosp.df$col,2],
          logp = round(-log10(geosp.df$geosp.pvalue), digits = 2)
        )

    }else{

      ## create a data.frame for drawing segments on ggplot
      ## x0 and y0 correspond to the rows (the samples under testing)
      ## x1 and y1 correspond to the columns (the K nearest neighbors of the samples under testing)
      link_geosp_df <-
        data.frame(
          x0 = geo_coord[geosp.df$row[geosp.df$sig.tag],1],
          y0 = geo_coord[geosp.df$row[geosp.df$sig.tag],2],
          x1 = geo_coord[geosp.df$col[geosp.df$sig.tag],1],
          y1 = geo_coord[geosp.df$col[geosp.df$sig.tag],2],
          logp = round(-log10(geosp.df$geosp.pvalue[geosp.df$sig.tag]), digits = 2)
        )
    }

  } # if map_type == "both" or "geographic_knn" end
  if(map_type == "both" | map_type == "genetic_knn"){
    # adjust transparency according to the type of figures
    if(is.null(p_thres)){
      ## create a data.frame for drawing segments on ggplot
      ## x0 and y0 correspond to the rows (the samples under testing)
      ## x1 and y1 correspond to the columns (the K nearest neighbors of the samples under testing)
      link_gensp_df <-
        data.frame(
          x0 = geo_coord[gensp.df$row,1],
          y0 = geo_coord[gensp.df$row,2],
          x1 = geo_coord[gensp.df$col,1],
          y1 = geo_coord[gensp.df$col,2],
          #col = alpha(gensp.df$gensp.col, alpha = tmp.alpha),
          logp = round(-log10(gensp.df$gensp.pvalue), digits = 2)
        )

    }else{
      ## create a data.frame for drawing segments on ggplot
      ## x0 and y0 correspond to the rows (the samples under testing)
      ## x1 and y1 correspond to the columns (the K nearest neighbors of the samples under testing)
      link_gensp_df <-
        data.frame(
          x0 = geo_coord[gensp.df$row[gensp.df$sig.tag],1],
          y0 = geo_coord[gensp.df$row[gensp.df$sig.tag],2],
          x1 = geo_coord[gensp.df$col[gensp.df$sig.tag],1],
          y1 = geo_coord[gensp.df$col[gensp.df$sig.tag],2],
          #col = alpha(gensp.df$gensp.col[gensp.df$sig.tag], alpha = tmp.alpha),
          logp = round(-log10(gensp.df$gensp.pvalue[gensp.df$sig.tag]), digits = 2)
        )

      # sig.row corresponds to the rows (the samples under testing)
      # sig.col corresponds to the columns (the K nearest neighbors of the samples under testing)
    }
  }

  ##---------------------------------------------
  ## make ggplot
  ## make a blank map
  geomap_plot <-
    ggplot2::ggplot(data = geomap) +
    ggplot2::geom_sf(color = map_color) +
    coord_sf(xlim = plot_xlim, ylim = plot_ylim, expand = TRUE) +
    theme(axis.title.x = element_text(color = "white"),
          axis.title.y = element_text(color = "white")
    )
  ## add links of geographic KNN to geographical map
  if(map_type == "both" | map_type == "geographic_knn"){
    geomap_plot_GeoSpKNN <-
      geomap_plot +
      geom_segment(aes(x=x0,
                       y=y0,
                       xend=x1,
                       yend=y1,
                       alpha = logp,
                       linewidth = logp
      ),
      data=link_geosp_df,
      lineend = "round",
      show.legend = NA,
      inherit.aes=FALSE, color = "blue") +
      labs(linewidth =bquote("-"~log[10]~"(p)"),
           alpha = bquote("-"~log[10]~"(p)"))+
      scale_linewidth(range = linewidth_range) +
      geom_point(mapping = aes(x = x, y = y), data = geo_coord, size = dot_cex)

    ## add pies
    if(!is.null(p_thres)){
      sig.row <- unique(geosp.df$row[geosp.df$sig.tag])
      if(length(sig.row) == 0){
        warning("there is no outlier in the region selected according to `select_xlim` and `select_ylim`")
      }else{

        if(show_knn_pie){# add pie charts of KNNs if show_knn_pie is TRUE
          sig_knn_index <- unique(unlist(lapply(sig.row, function(i){GeoSP_knn_res$knn_index[i,]})))
          sig.row <- c(sig.row, sig_knn_index)
        } # show_knn_pie end

        sig_anc_ggplot_df <- gather(anc_coef_for_pie[sig.row,], "type", "value", !!cols2)
        sig_anc_ggplot_df$type <- factor(sig_anc_ggplot_df$type, levels = cols) # set legend order based on order of "cols"
        geomap_plot_GeoSpKNN <-
          geomap_plot_GeoSpKNN +
          scale_fill_manual(values=pie_color) +
          geom_arc_bar(aes(x0 = x, y0 = y,r0=0,
                           fill = type, r=pie_r_scale,
                           amount=value),
                       data=sig_anc_ggplot_df, stat='pie',
                       col = "black", show.legend = FALSE)
      } # if no significant end
    } # add pie Geographic KNN end
  } # make geomap_plot_GeoSpKNN end
  if(map_type == "both" | map_type == "genetic_knn"){
    geomap_plot_GeneticKNN <-
      geomap_plot +
      geom_segment(aes(x=x0,
                       y=y0,
                       xend=x1,
                       yend=y1,
                       alpha = logp,
                       linewidth = logp
      ),
      data=link_gensp_df,
      lineend = "round",
      show.legend = NA,
      inherit.aes=FALSE, color = "red") +
      labs(linewidth =bquote("-"~log[10]~"(p)"),
           alpha = bquote("-"~log[10]~"(p)"))+
      scale_linewidth(range = linewidth_range) +
      geom_point(mapping = aes(x = x, y = y), data = geo_coord, size = dot_cex)

    ## add pies
    if(!is.null(p_thres)){
      sig.row <- unique(gensp.df$row[gensp.df$sig.tag])
      if(length(sig.row) == 0){
        warning("there is no outlier in the region selected according to `select_xlim` and `select_ylim`")
      }else{
        if(show_knn_pie){# add pie charts of KNNs if show_knn_pie is TRUE
          sig_knn_index <- unique(unlist(lapply(sig.row, function(i){GenSP_knn_res$knn_index[i,]})))
          sig.row <- c(sig.row, sig_knn_index)
        } # show_knn_pie end

        sig_anc_ggplot_df <- gather(anc_coef_for_pie[sig.row,], "type", "value", !!cols2)
        sig_anc_ggplot_df$type <- factor(sig_anc_ggplot_df$type, levels = cols) # set legend order based on order of "cols"
        geomap_plot_GeneticKNN <-
          geomap_plot_GeneticKNN +
          scale_fill_manual(values=pie_color) +
          geom_arc_bar(aes(x0 = x, y0 = y,r0=0,
                           fill = type, r=pie_r_scale,
                           amount=value),
                       data=sig_anc_ggplot_df, stat='pie',
                       col = "black", show.legend = FALSE)
      } # if no significant end
    } # add pie Genetic KNN end
  } # make geomap_plot_GeneticKNN end

  ##-----------------------------------
  ## output plot
  ## -> combine plots using `cowplot`
  if(map_type == "both"){
    if(add_benchmark_graph){
      plot_out <-
        cowplot::plot_grid(
          geomap_plot_benchmark,
          geomap_plot_GeneticKNN,
          geomap_plot_GeoSpKNN,
          labels = plot_labels, ncol = 1
       )
    }else{
      plot_out <-
        cowplot::plot_grid(
          geomap_plot_GeneticKNN,
          geomap_plot_GeoSpKNN,
          labels = plot_labels, ncol = 1
        )
    }
    return(plot_out)
  } # if both end
  if(map_type == "genetic_knn"){
    if(add_benchmark_graph){
      plot_out <-
        cowplot::plot_grid(
          geomap_plot_benchmark,
          geomap_plot_GeneticKNN,
          labels = plot_labels, ncol = 1
        )
      return(plot_out)
    }else{
      return(geomap_plot_GeneticKNN)
    }
  } # if genetic_knn end
  if(map_type == "geographic_knn"){
    if(add_benchmark_graph){
      plot_out <-
        cowplot::plot_grid(
          geomap_plot_benchmark,
          geomap_plot_GeoSpKNN,
          labels = plot_labels, ncol = 1
        )
      return(plot_out)
    }else{
      return(geomap_plot_GeoSpKNN)
    }
  }
} # plot_ggoutlier end


#-------------------------------------------------------------
# assign colors to edges
# `get_GGNet_edge_col` is a function used in the `plot_ggoutlier`
#-------------------------------------------------------------
# Argument:
# GGNet_adjm: an output of `get_GGNet_adjacency_matrix
# color_res: the resolution of color scale
get_GGNet_edge_col <- function(GGNet_adjm, color_res = 10){

  GeoSP.colfun <- colorRampPalette(c("white", "dodgerblue"))
  GeoSP.col <- GeoSP.colfun(color_res)
  if(!is.null(GGNet_adjm$GeoSP_pvalue)){
    GeoSP_col <- rep(NA, length(GGNet_adjm$GeoSP_pvalue))
    geosp.indx <- c(GGNet_adjm$GeoSP_pvalue) > 0
    GeoSP_col[geosp.indx] <- as.character(GeoSP.col[as.numeric(cut(-log10(GGNet_adjm$GeoSP_pvalue[geosp.indx]), breaks = color_res))])
  }

  GenSP.colfun <- colorRampPalette(c("white", "orangered"))
  GenSP.col <- GenSP.colfun(color_res)
  if(!is.null(GGNet_adjm$GenSP_pvalue)){
    GenSP_col <- rep(NA, length(GGNet_adjm$GenSP_pvalue))
    gensp.indx <- c(GGNet_adjm$GenSP_pvalue) > 0
    GenSP_col[gensp.indx] <- as.character(GenSP.col[as.numeric(cut(-log10(GGNet_adjm$GenSP_pvalue[gensp.indx]), breaks = color_res))])
  }

  # if only geographical KNN result available
  if(!is.null(GGNet_adjm$GeoSP_pvalue) & is.null(GGNet_adjm$GenSP_pvalue)){
    return(list(GeoSP_col = c(GeoSP_col),
                GeoSP_colkey = c(GeoSP.col)
    ))
  }
  # if only genetic KNN result available
  if(is.null(GGNet_adjm$GeoSP_pvalue) & !is.null(GGNet_adjm$GenSP_pvalue)){
    return(list(GenSP_col = c(GenSP_col),
                GenSP_colkey = c(GenSP.col)
    ))
  }

  # if both geographical KNN result and genetic KNN result
  if(!is.null(GGNet_adjm$GeoSP_pvalue) & !is.null(GGNet_adjm$GenSP_pvalue)){
    return(list(GeoSP_col = c(GeoSP_col),
                GenSP_col = c(GenSP_col),
                GeoSP_colkey = c(GeoSP.col),
                GenSP_colkey = c(GenSP.col)
    ))
  }

} # get_GGNet_edge_col end
