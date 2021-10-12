#' Beta value error plot with CpG coordinates
#' 
#' @description Plots methylation beta values in the form of error bars. For
#' every group of samples, the mean methylation value is plotted together
#' with the confidence interval. Depending on the function parameters, CpG sites 
#' represented in the plot can be plotted equidistantly or their position on
#' the x-axis can represent their true genomic coordinates. In the latter case,
#' additional information can be included in the plot such as position of genes
#' and exons.
#' 
#' @param beta_values Data frame with methylation beta values. CpG sites should
#' be listed in the rows with the row names being the Illumina CpG identifiers.
#' Samples should be listed in the columns.
#' @param sample_groups Data frame which defines the sample groups. It has to 
#' have two columns, the first column with sample IDs and the second column 
#' with the corresponding group.
#' @param cg_list List of CpG sites whose methylation values should be presented
#' in the plot. If provided, CpGs will be plotted equidistantly on the x-axis in
#' the given order.
#' @param chr Chromosome on which the plotted CpG sites are found. In the format
#' `1`:`22` or `"X"` or `"Y"`.
#' @param first_coord The first coordinate to be represented in the plot.
#' @param last_coord The last coordinate to be represented in the plot.
#' @param add_lines Whether to connect the error bars. If `TRUE`, mean value
#' points will be connected with the lines. Default to `TRUE`.
#' @param plot_cg Whether to plot the position of all CpG sites on x-axis.
#' Defaults to `FALSE`.
#' @param plot_cgi Whether to plot the position of CpG island on x-axis.
#' Defaults to `FALSE`.
#' @param plot_gene Whether to plot the position of genes/transcripts on x-axis.
#' Defaults to `FALSE`.
#' @param plot_exon Whether to plot the position of exons on x-axis. Defaults to
#' `FALSE`.
#' @param transcript_types Genes/transcripts to be plotted. `"NM"` stands for
#' protein-coding and `"NR"` stands for non-protein coding transcripts.
#' By default, both types of transcripts are included.
#' @param col_groups Color used for the error bars. It can be given as a list
#' of colors of the same length as the number of groups. Alternatively, it
#' can be given as a data frame with sample groups in the first column and
#' the corresponding color in the second column.
#' @param col_cg Color of the other CpG sites. Used only if `plot_cg` is `TRUE`.
#' Default color is `grey77`.
#' @param col_gene Color of the genes/exons. Used only if `plot_genes` or
#' `plot_exons` is `TRUE`. Default color is `"cyan4"`.
#' @param col_cgi Color of the CpG islands. Used only if `plot_cgi` is `TRUE`.
#' Default color is `"grey77"`.
#' @param title Plot title. By default, omitted.
#' @param x_label X-axis label. By default, omitted.
#' @param y_label Y-axis label. By default, omitted.
#' @param text_cex Numeric character expansion factor, used for all displayed 
#' text. Defaults to `1`.
#' @param beta_min Minimum beta value to be displayed. Has to be in `[0,1>` range.
#' Defaults to `0`.
#' @param beta_max Maximum beta value to be displayed. Has to be in `<0,1]` range.
#' Defaults to `1`.
#' @param error_size Number of standard deviations which are to be plotted
#' below and above the mean. Defaults to `1`.
#' @param cg_names Whether to add names of the plotted CpG sites below the
#' x-axis. Defaults to `FALSE`.
#' @param cg_names_expand Expand factor for the `cg_names`. If the labels are 
#' cluttered, the expand factor could be set to separate the labels. Expand 
#' factor defines the minimum distance between two consecutive labels and is
#' expressed as a percentage of the distance between the first and the last
#' label. Optimal separation is usually achieved for values close to `1`.
#' For the default `NA` value, labels are not expanded.
#' @param plot_legend Whether to plot the legend. Default to `FALSE`.
#' @param legend_cex Numeric expansion factor, used for all elements of the
#' legend. Defaults to `1`.
#' @param plot_outside Whether error plot should be plotted outside of the
#' plotting area. Plotting area is determined with `beta_min` and `beta_max`. 
#' Defaults to `TRUE`.
#' 
#' @return Vector with names of plotted CpG sites, in the same order as in the
#' plot.
#'
#' @details Plots error bars and confidence intervals for beta values of given 
#' CpGs. If a list of CpGs is given, they are plotted equidistantly. 
#' If a range of coordinates and a chromosome are given, then all CpGs inside 
#' the range are plotted and positioned on x-axis based on their chromosomal 
#' coordinate (from 5' end to 3' end). 
#' If the latter case, additional information can be displayed on
#' the x-axis such as positions of unmeasured CpGs, positions of CpG islands
#' and/or positions of genes and exons. If we are specifically interested in
#' protein-coding or non-protein-coding genes, we can accordingly set the values
#' of parameter `transcript_types`. The color for each genomic element can be 
#' changed with `col_cg`, `col_cgi` and `col_gene` parameters.
#' 
#' For error bars, the default plotted confidence interval represents distance of
#' one standard deviation from the mean. This can be changed with `error_size` 
#' parameter. If we are interested only in certain beta values, we can restrict the
#' y-axis with `beta_min` and `beta_max` parameters.
#' Error bars are separately plotted for each group of samples and the color of each
#' group is set with `col_groups` parameter. Groups of samples are defined in
#' `sample_groups`. To better observe the methylation patterns, lines which 
#' connect mean methylation values can be added with `add_lines`.
#' 
#' Plotting area can be customized with `title`, `x_label` and `y_label`. 
#' Below the x-axis, names of CpG sites can be displayed setting `cg_names` to
#' `TRUE`. In cg names are cluttered, we can expand them setting `cg_names_expand`
#' factor to an appropriate numerical value. Usually this value is a number close
#' to `1`. An increase of the expand factor will result in more separate labels.
#' For the optimal separation, the best strategy is to try several different
#' values until the right value is found. The size of all displayed text is 
#' regulated simultaneously with `text_cex` parameter.
#' 
#' The legend can be added to the plot with `plot_legend` parameter. There are 
#' two parts of the legend. Color legend for sample groups and the legend which
#' explains the meaning of plotted genomic elements. The size of the whole legend
#' is regulated with `legend_cex` parameter.
#' 
#' In some cases, error bars have values outside the `beta_min`-`beta_max` range.
#' If `plot_outside` is `FALSE`, these values will be trimmed. Otherwise, they will
#' be plotted. If `beta_min` is `0` and `beta_max` is `1`, there still may be values
#' outside of the range. This happens for methylation values which are very close 
#' to `0` or `1` since the confidence interval is wider than the real methylation
#' values. In this case, it is recommended for `plot_outside` to be set to `TRUE`.
#' 
#' 
#' @import graphics
#' 
#' @export
plotCI_meCpG <- function(beta_values, sample_groups,
                         cg_list=NULL, chr=NULL, first_coord=NULL, last_coord=NULL,
                         add_lines=TRUE, plot_cg=FALSE, plot_cgi=FALSE,
                         plot_gene=FALSE, plot_exon=FALSE, transcript_types=c("NM", "NR"),
                         col_groups=NULL, col_cg="grey77", col_gene="cyan4", col_cgi="grey77",
                         title="", x_label="", y_label="", text_cex=1,
                         beta_min=0, beta_max=1, error_size=1, 
                         cg_names=FALSE, cg_names_expand=NA,
                         plot_legend=FALSE, legend_cex=1, plot_outside=TRUE){
  
  data("CpG_annoHM450K", package="meNet")
  data("UCSC_genes", package="meNet")
  data("UCSC_CGI", package="meNet")
  
  cg_in_plot <- c() # just to hold a space
  
  # check validity of given parameters
  ####################################
  # check if CpGs are provided
  if(is.null(cg_list) & (is.null(chr)|is.null(first_coord)|is.null(last_coord))){
    stop('CpGs not specified.')
  }
  
  # check logical values
  if(!is.logical(add_lines)){
    stop('"add_lines" has to have a logical value.')
  }
  if(!is.logical(plot_cg)){
    stop('"plot_cg" has to have a logical value.')
  }
  if(!is.logical(plot_gene)){
    stop('"plot_gene" has to have a logical value.')
  }
  if(!is.logical(plot_exon)){
    stop('"plot_exon" has to have a logical value.')
  }
  if(!is.logical(plot_cgi)){
    stop('"plot_cgi" has to have a logical value.')
  }
  if(!is.logical(cg_names)){
    stop('"cg_names" has to have a logical value.')
  }
  if(!is.logical(plot_legend)){
    stop('"plot_legend" has to have a logical value.')
  }
  if(!is.logical(plot_outside)){
    stop('"plot_outside" has to have a logical value.')
  }
  
  # check if cex values are numeric
  if(!is.numeric(text_cex)){
    stop('"text_cex" has to have a numeric value.')
  }
  if(!is.numeric(legend_cex)){
    stop('"legend_cex" has to have a numeric value.')
  }
  
  # check for cg_names_expand
  if(!is.na(cg_names_expand)&!is.numeric(cg_names_expand)){
    stop('"cg_names_expand" has to either have a numeric value or be set to NA.')
  }else if(is.numeric(cg_names_expand)&cg_names_expand<0){
    stop('"cg_names_expand" has to have a non-negative value.')
  }
  
  # transcript types
  transcript_types <- intersect(transcript_types, c("NM", "NR"))
  if(length(transcript_types)==0 & (plot_exon|plot_gene)){
    stop('"transcript_types" incorrectly specified.')
  }
  
  # check titles/labels
  if(!is.character(title)){
    stop('"title" has to be of character type.')
  }
  if(!is.character(x_label)){
    stop('"x_label" has to be of character type.')
  }
  if(!is.character(y_label)){
    stop('"y_label" has to be of character type.')
  }
  
  # check min/max beta values
  if(!is.numeric(beta_min)){
    warning('"beta_min" has to have a numeric value. Given value is discarded.')
    y_min <- 0
  }else if(beta_min<0|beta_min>=1){
    warning('"beta_min" out of bounds. Given value is discarded.')
    y_min <- 0
  }else{
    y_min <- beta_min
  }
  #
  if(!is.numeric(beta_max)){
    warning('"beta_max" has to have a numeric value. Given value is discarded.')
    y_max <- 1
  }else if(beta_max<=0|beta_max>1){
    warning('"beta_max" out of bounds. Given value is discarded.')
    y_max <- 1
  }else{
    y_max <- beta_max
  }
  
  # check error_size
  if(!is.numeric(error_size)){
    warning('"error_size" has to have a numeric value. Given value is discarded.')
    error_size <- 1
  }else if(error_size<0){
    warning('"error_size" has to be a non-negative number. Given value is discarded.')
    error_size <- 1
  }
  
  # beta values
  if(!is.data.frame(beta_values)){
    stop('"beta_values" has to be a data frame.')
  }
  if(nrow(beta_values)==0|ncol(beta_values)==0){
    stop('"beta_values" missing data.')
  }
  if(!all(rownames(beta_values)%in%CpG_annoHM450K[,"IlmnID"])){
    stop('Row names of "beta_values" are not CpGs measured by Illumina HM450K microarray.')
  }
  
  # sample groups
  if(!is.data.frame(sample_groups)){
    stop('"sample_groups" has to be a data frame.')
  }else if(ncol(sample_groups)<2 ){
    stop('"sample_groups" has to have two columns, the first with sample IDs and the second with the corresponding group.')
  }else if(!all(colnames(beta_values)%in%sample_groups[,1])){
    stop('First column of "sample_groups" has to contain all column names of "beta_values"')
  }
  
  # remove samples not present in "beta_values"
  sample_groups <- sample_groups[sample_groups[,1]%in%colnames(beta_values),]
  # remove possible factors
  sample_groups <- data.frame(lapply(sample_groups, as.character), stringsAsFactors=FALSE)
  
  # groups of samples
  all_groups <- unique(sample_groups[,2])
  
  # colors
  # cg, gene and cgi colors
  if(!.is_color(col_cg)|length(col_cg)>1){
    warning('"col_cg" has an incorrect value. Value set to default.')
    col_cgi <- "grey77"
  }
  if(!.is_color(col_gene)|length(col_gene)>1){
    warning('"col_gene" has an incorrect value. Value set to default.')
    col_gene <- "cyan4"
  }
  if(!.is_color(col_cgi)|length(col_cgi)>1){
    warning('"col_cgi" has an incorrect value. Value set to default.')
    col_cgi <- "grey77"
  }
  
  # colors for groups
  col_df <- data.frame(group=all_groups)
  # can be data frame or vector; if wrongly specified the same as NULL
  if(!is.null(col_groups)){
    if(is.data.frame(col_groups)){
      if(all_groups%in%col_groups[,1] & ncol(col_groups)>1 & all(.is_color(col_groups[,2]))){
        colnames(col_groups)[1:2] <- c("group", "color")
        col_df <- plyr::join(col_df, col_groups, by="group", type="left", match="first")
      }else{
        warning('"col_groups" has to have two columns with the first column containing group names and the second corresponding colors.')
        col_groups <- NULL
      }
    }else if(length(col_groups)!=length(all_groups)){
      warning('"col_groups" has to be of the same length as the number of groups.')
      col_groups <- NULL
    }else if(!all(.is_color(col_groups))){
      warning('"col_groups" contains undefined colors.')
      col_groups <- NULL
    }else{
      col_df$color <- col_groups
    }
  }
  if(is.null(col_groups)){
    # we specify random color
    col_df$color <- viridisLite::viridis(nrow(col_df), alpha = 0.8, begin = 0.1, end = 0.9, direction = -1, option = "C")
  }
  
  
  # prepare for plotting
  ######################
  if(!is.null(cg_list)){
    
    if(!is.null(chr)|!is.null(first_coord)|!is.null(last_coord)){
      warning('Only "cg_list" was used to specify the list of CpGs.')
    }
    #cg_in_plot <- intersect(cg_list, rownames(beta_values))
    cg_in_plot <- cg_list[cg_list%in%rownames(beta_values)]
    if(length(cg_in_plot)==0){
      stop("Given CpGs don't have measured beta values.")
    }
    # define plotting variables
    plot_cg <- FALSE
    plot_cgi <- FALSE
    plot_gene <- FALSE
    plot_exon <- FALSE
    x_min <- 0.5
    x_coords <- 1:length(cg_in_plot)
    x_max <- length(cg_in_plot)+0.5
    
  }else{
    
    # find CpGs
    CpG_annoHM450K <- CpG_annoHM450K[(CpG_annoHM450K$CHR==chr)&(CpG_annoHM450K$MAPINFO>=first_coord)&(CpG_annoHM450K$MAPINFO<=last_coord),]
    if(nrow(CpG_annoHM450K)==0){
      stop("No CpG found in the specified interval. Check chromosome name, the first and the last coordinate.")
    }
    # check CpGs have measured beta values
    CpG_annoHM450K <- CpG_annoHM450K[CpG_annoHM450K$IlmnID%in%rownames(beta_values),]
    if(nrow(CpG_annoHM450K)==0){
      stop("No CpG with measured beta values found in the specified interval.")
    }
    # sort CpGs
    CpG_annoHM450K <- CpG_annoHM450K[with(CpG_annoHM450K, order(CHR, MAPINFO)),]
    cg_in_plot <- CpG_annoHM450K$IlmnID
    
    x_coords <- CpG_annoHM450K$MAPINFO
    x_min <- first_coord
    x_max <- last_coord
    
    if(plot_cg){
      cg_all_coords <- find_CpG(chromosome=paste0("chr",chr), min_coord=first_coord, max_coord=last_coord)
      if(length(cg_all_coords)==0){
        plot_cg <- FALSE
      }
    }
    
    if(plot_cgi){
      # UCSC_CGI had chromosome format: 1:22, X, Y
      UCSC_CGI <- UCSC_CGI[UCSC_CGI$chr==chr,]
      # match coordinates
      UCSC_CGI <- UCSC_CGI[UCSC_CGI$start<last_coord & UCSC_CGI$end>first_coord,]
      
      if(nrow(UCSC_CGI)>0){
        cgi_starts <- sapply(UCSC_CGI$start, function(x) max(x,first_coord))
        cgi_ends <- sapply(UCSC_CGI$end, function(x) min(x,last_coord))
      }else{
        plot_cgi <- FALSE
      }
    }
    
    if(plot_gene | plot_exon){
      # keep only relevant genes
      # chr: different notation so we need to add "chr" before the number
      UCSC_genes <- UCSC_genes[UCSC_genes$chrom==paste0("chr",chr),]
      # coord: txStart 
      UCSC_genes <- UCSC_genes[UCSC_genes$txStart<last_coord & UCSC_genes$txEnd>first_coord,]
      # only correct transcript types
      UCSC_genes <- UCSC_genes[grepl(paste0("^", transcript_types, collapse="|"), UCSC_genes[,1]),]
      #
      
      if(nrow(UCSC_genes)>0){
        if(plot_gene){
          # where to plot genes:
          gene_starts <- sapply(UCSC_genes$txStart, function(x) max(x,first_coord))
          gene_ends <- sapply(UCSC_genes$txEnd, function(x) min(x,last_coord))
        }
        if(plot_exon){
          exon_starts <- unlist(sapply(UCSC_genes$exonStarts, function(x) as.integer(strsplit(x, ',')[[1]])),use.name=FALSE)
          exon_ends <- unlist(sapply(UCSC_genes$exonEnds, function(x) as.integer(strsplit(x, ',')[[1]])),use.name=FALSE)
          exons_to_keep <- exon_starts<last_coord & exon_ends>first_coord
          exon_starts <- exon_starts[exons_to_keep]
          exon_ends <- exon_ends[exons_to_keep]
          # where to plot them:
          exon_starts <- sapply(exon_starts, function(x) max(x,first_coord))
          exon_ends <- sapply(exon_ends, function(x) min(x,last_coord))
        }
      }else{ #no genes or exons found
        plot_gene <- FALSE
        plot_exon <- FALSE
      }
    }
  }
  
  # cg_in_plot is already in the right order
  beta_values <- beta_values[cg_in_plot,]

  # mean and sd: list of vectors
  cg_mean_perGroup <- lapply(all_groups, function(x) apply(beta_values[,sample_groups[sample_groups[,2]==x,1]], 1, mean))
  names(cg_mean_perGroup) <- all_groups
  cg_sd_perGroup <- lapply(all_groups, function(x) apply(beta_values[,sample_groups[sample_groups[,2]==x,1]], 1, sd))
  names(cg_sd_perGroup) <- all_groups
  
  # y-coordinates of the error bars
  cg_errorCoord_lower <-  lapply(all_groups, function(x) cg_mean_perGroup[[x]]-error_size*cg_sd_perGroup[[x]])
  names(cg_errorCoord_lower) <- all_groups
  cg_errorCoord_upper <-  lapply(all_groups, function(x) cg_mean_perGroup[[x]]+error_size*cg_sd_perGroup[[x]])
  names(cg_errorCoord_upper) <- all_groups
  
  # adjust values based on y_min and y_max
  min_plotted_beta <- min(sapply(all_groups, function(x) min(cg_errorCoord_lower[[x]])))
  max_plotted_beta <- max(sapply(all_groups, function(x) max(cg_errorCoord_upper[[x]])))
  if(min_plotted_beta<y_min|max_plotted_beta>y_max){
    if(plot_outside){
      warning('Plotting outside of bounds.') #sometimes makes sense (for beta values close to 0/1)
    }else{
      cg_errorCoord_lower <-  lapply(all_groups, function(x) cg_errorCoord_lower[[x]]*(cg_errorCoord_lower[[x]]>=y_min)+y_min*(cg_errorCoord_lower[[x]]<y_min))
      cg_errorCoord_upper <-  lapply(all_groups, function(x) cg_errorCoord_upper[[x]]*(cg_errorCoord_upper[[x]]<=y_max)+y_max*(cg_errorCoord_upper[[x]]>y_max))
    }
  }
  
  
  # PLOTTING
  ##############################
  
  # device size
  device_width <- grDevices::dev.size(units="px")[1]
  device_height <- grDevices::dev.size(units="px")[2]
  # universal unit is one thousandth of the smaller dimension
  universal_unit <- min(device_width, device_height)/1000 # one thousandth of screen in pixels
  pixel_uu <-1/universal_unit # one pixel in universal units
  
  # plot dimensions
  #################
  plottingArea_width <- x_max - x_min
  plottingArea_height <- y_max - y_min
 
  plot_x_min <- x_min
  plot_x_max <- x_max
  plot_y_min <- y_min
  plot_y_max <- y_max
  
  # plotting parameters: expressed in universal units
  plotting_parameters <- list(x_axis_distance=30, x_tick_length=15,
                              x_label_space=70*text_cex,
                              y_axis_distance=20, y_tick_length=15,
                              y_number_space=30*text_cex, # numbers are adjusted to the left border
                              y_label_space=70*text_cex,
                              title_space=100*text_cex, legend_space=170*legend_cex,
                              cg_names_space=100*text_cex,
                              exon_width=7, cgi_width=5,
                              error_cap_length=8, 
                              margin=30)
  # plotting parameters for the legend
  if(plot_legend){
    legend_parameters <- list(left_space=10, between_line_height=20,
                              line_width=55, line_text_space=12,
                              polygon_height=6, cg_height=8,
                              group_line_lwd=1.5,
                              point_size_multiplier=1,
                              text_size=0.5)
    legend_parameters <- sapply(legend_parameters, function(x) x*legend_cex, simplify=FALSE, USE.NAMES=TRUE)
  }
  
  # if legend is plotted, check if there are some long names
  if(plot_legend){
    groupName_maxWidth <- max(sapply(all_groups, function(x) strwidth(x, units="inches")))
    legendString_optimalWidth <- strwidth("Gene w/ exons", units="inches")
    # find maximal string width within the legend
    legendString_maxWidth <- 0
    if(plot_gene&plot_exon){
      legendString_maxWidth <- strwidth("Gene w/ exons", units="inches")
    }else if(plot_cgi){
      legendString_maxWidth <- strwidth("CpG island", units="inches")
    }else if(plot_cg){
      legendString_maxWidth <- strwidth("CpG sites", units="inches")
    }else if(plot_gene|plot_exon){
      legendString_maxWidth <- strwidth("Gene", units="inches")
    }
    legendString_maxWidth <- max(legendString_maxWidth, groupName_maxWidth)
    
    # space for text within the legend
    legend_text_percentage <- (plotting_parameters[["legend_space"]]-legend_parameters[["left_space"]]-legend_parameters[["line_width"]]-legend_parameters[["line_text_space"]])/plotting_parameters[["legend_space"]]
    # additional name size of maxSize should doule the text length: legend_text_percentage=maxSize/?
    plotting_parameters[["legend_space"]] <- (1+legend_text_percentage*(legendString_maxWidth/legendString_optimalWidth-1))*plotting_parameters[["legend_space"]]
  }
  
  # text size: these numbers are multiplied by text_cex as parameters of graphics::text()
  text_size_parameters <- list(title=1.3, label=0.9, number=0.75, cg_names=0.6)
  text_size_parameters <- sapply(text_size_parameters, function(x) x*text_cex, simplify=FALSE, USE.NAMES=TRUE)
  # size of point that depicts the mean within error bars
  point_size <- 0.8
  
  # the remaining space is reserved for the plotting area
  plottingArea_width_uu <- device_width/universal_unit - 
                           (plotting_parameters[["y_axis_distance"]] +
                            plotting_parameters[["y_tick_length"]] +
                            plotting_parameters[["y_number_space"]] +
                            plotting_parameters[["y_label_space"]]*(y_label!="") +
                            plotting_parameters[["legend_space"]]*plot_legend +
                            plotting_parameters[["margin"]]*2)
  
  plottingArea_height_uu <- device_height/universal_unit - 
                            (plotting_parameters[["x_axis_distance"]] +
                             plotting_parameters[["x_tick_length"]] +
                             plotting_parameters[["cg_names_space"]]*(cg_names) +
                             plotting_parameters[["x_label_space"]]*(x_label!="") +
                             plotting_parameters[["title_space"]]*(title!="") +
                             plotting_parameters[["margin"]]*2)
  
  # conversion from universal units to window coordinates
  dx <- plottingArea_width/plottingArea_width_uu
  dy <- plottingArea_height/plottingArea_height_uu
  
  # space for axes and ticks
  plot_x_min <- plot_x_min - (plotting_parameters[["y_axis_distance"]] +
                              plotting_parameters[["y_tick_length"]] +
                              plotting_parameters[["y_number_space"]] +
                              plotting_parameters[["y_label_space"]]*(y_label!="") +
                              plotting_parameters[["margin"]])*dx
  plot_x_max <- plot_x_max + (plotting_parameters[["legend_space"]]*plot_legend +
                              plotting_parameters[["margin"]])*dx
  
  plot_y_min <- plot_y_min - (plotting_parameters[["x_axis_distance"]] +
                              plotting_parameters[["x_tick_length"]] +
                              plotting_parameters[["cg_names_space"]]*(cg_names) +
                              plotting_parameters[["x_label_space"]]*(x_label!="") +
                              plotting_parameters[["margin"]])*dy
  plot_y_max <- plot_y_max + (plotting_parameters[["title_space"]]*(title!="") +
                              plotting_parameters[["margin"]])*dy
  
  # create plot
  #############
  par(mar=c(0,0,0,0))
  par(lend=2) # to make lines have squared ends
  plot.new()
  plot.window(xlim=c(plot_x_min, plot_x_max), ylim=c(plot_y_min, plot_y_max), xaxs="i", yaxs="i") # xaxs='i' and yaxs="i inhibit expansion
  
  # axes
  ######
  # cg_names for x-axis; axis and ticks should be added later to be plotted above cgis/genes
  if(cg_names){
    if(is.numeric(cg_names_expand)){
      cg_names_coords <- .new_position_for_label(x_coords, cg_names_expand)
    }else{
      cg_names_coords <- x_coords
    }
    text(x=cg_names_coords, y=y_min-(plotting_parameters[["x_axis_distance"]]+plotting_parameters[["x_tick_length"]]+5*pixel_uu)*dy,
         labels=cg_in_plot, adj=c(1,0.5), srt=90, cex=text_size_parameters[["cg_names"]])
  }
  
  # y-axis with ticks and numbers
  segments(x0=x_min-plotting_parameters[["y_axis_distance"]]*dx, y0=y_min, y1=y_max)
  y_axis_numbers <- c(0.0,0.2,0.4,0.6,0.8,1.0)
  y_axis_numbers <- y_axis_numbers[y_axis_numbers>=y_min&y_axis_numbers<=y_max]
  if(!(y_min %in% y_axis_numbers)){
    y_axis_numbers <- c(y_min, y_axis_numbers)
  }
  if(!(y_max %in% y_axis_numbers)){
    y_axis_numbers <- c(y_axis_numbers, y_max)
  }
  segments(x0=x_min-plotting_parameters[["y_axis_distance"]]*dx, y0=y_axis_numbers, 
           x1=x_min-(plotting_parameters[["y_axis_distance"]]+plotting_parameters[["y_tick_length"]])*dx)
  text(x=x_min-(plotting_parameters[["y_axis_distance"]]+plotting_parameters[["y_tick_length"]]+plotting_parameters[["y_number_space"]])*dx, 
       y=y_axis_numbers, labels=y_axis_numbers, adj=c(0.5,1), srt=90, cex=text_size_parameters[["number"]])
  
  # title and axis labels
  #######################
  if(y_label!=""){ #just right to the left margin
    text(x=plot_x_min+plotting_parameters[["margin"]]*dx,
         y=mean(c(y_min,y_max)), labels=y_label, adj=c(0.5,1), srt=90, cex=text_size_parameters[["label"]])
  }
  if(x_label!=""){ #just above the margin
    text(x=mean(c(x_min,x_max)), 
         y=plot_y_min+plotting_parameters[["margin"]]*dy,
         labels=x_label, adj=c(0.5,0), cex=text_size_parameters[["label"]])
  }
  if(title!=""){ #just below the margin
    text(x=mean(c(x_min,x_max)), 
         y=plot_y_max-plotting_parameters[["margin"]]*dy,
         labels=title, adj=c(0.5,1), cex=text_size_parameters[["title"]], font=2)
  }
  
  
  # all cpgs/cgis/genes/exons
  ######################
  # CpG islands will be plotted above the x-axis
  if(plot_cgi){
    for(i in 1:length(cgi_starts)){
      polygon(x=c(cgi_starts[i], cgi_ends[i], cgi_ends[i], cgi_starts[i]),
              y=y_min-plotting_parameters[["x_axis_distance"]]*dy+c(0,0,1,1)*plotting_parameters[["cgi_width"]]*dy+pixel_uu*dy,
              border=col_cgi, col=col_cgi)
    }
  }
  # genes will be plotted below the x-axis from the first to the last transcription site
  if(plot_gene){
    segments(x0=gene_starts, y0=y_min-plotting_parameters[["x_axis_distance"]]*dy-pixel_uu*dy,
             x1=gene_ends, col=col_gene)
  }
  if(plot_exon){
    for(i in 1:length(exon_starts)){
      polygon(x=c(exon_starts[i], exon_ends[i], exon_ends[i], exon_starts[i]),
              y=y_min-plotting_parameters[["x_axis_distance"]]*dy-pixel_uu*dy+
                c(-1,-1,0,0)*plotting_parameters[["exon_width"]]*dy,
              border=col_gene, col=col_gene)
    }
  }
  
  # other CpGs are plotted above the x-axis
  if(plot_cg){
    segments(x0=cg_all_coords, y0=y_min-plotting_parameters[["x_axis_distance"]]*dy+pixel_uu*dy, 
             y1=y_min-(plotting_parameters[["x_axis_distance"]]-plotting_parameters[["x_tick_length"]])*dy,
             col=col_cg)
  }
  
  # x-axis with ticks
  # CpGs: ticks go both below and above the x-axis
  segments(x0=x_min, y0=y_min-plotting_parameters[["x_axis_distance"]]*dy, x1=x_max)
  segments(x0=x_coords, y0=y_min-(plotting_parameters[["x_axis_distance"]]-plotting_parameters[["x_tick_length"]])*dy, 
           y1=y_min-(plotting_parameters[["x_axis_distance"]]+plotting_parameters[["x_tick_length"]])*dy)
  
  
  # error bars
  ############
  # lines
  if(add_lines){
    n_points <- length(x_coords)
    for(i in 1:length(all_groups)){
      segments(x0=x_coords[1:(n_points-1)], x1=x_coords[2:n_points],
               y0=cg_mean_perGroup[[i]][1:(n_points-1)], y1=cg_mean_perGroup[[i]][2:n_points],
               col=col_df[i,"color"])
    }
  }
  # mean points and error bars
  for(i in 1:length(all_groups)){
    segments(x0=x_coords, y0=cg_errorCoord_lower[[i]], y1=cg_errorCoord_upper[[i]], 
             col=col_df[i,"color"])
    points(x=x_coords, y=cg_mean_perGroup[[i]], pch=16, col=col_df[i,"color"], cex=point_size)
    #caps
    segments(x0=x_coords-plotting_parameters[["error_cap_length"]]*dx/2,
             x1=x_coords+plotting_parameters[["error_cap_length"]]*dx/2,
             y0=cg_errorCoord_lower[[i]], col=col_df[i,"color"])
    segments(x0=x_coords-plotting_parameters[["error_cap_length"]]*dx/2,
             x1=x_coords+plotting_parameters[["error_cap_length"]]*dx/2,
             y0=cg_errorCoord_upper[[i]], col=col_df[i,"color"])
  }
  
  
  # legend
  ########
  if(plot_legend){
    # plotting parameters
    legend_x0 <- x_max+legend_parameters[["left_space"]]*dx
    legend_x1 <- legend_x0 + legend_parameters[["line_width"]]*dx
    
    # groups
    segments(y0=y_max-0:(length(all_groups)-1)*legend_parameters[["between_line_height"]]*dy,
             x0=legend_x0, x1=legend_x1, col=scales::alpha(col_df$color,1), 
             lwd=legend_parameters[["group_line_lwd"]])
    points(x=rep((legend_x0+legend_x1)/2, length(all_groups)), 
           y=y_max-0:(length(all_groups)-1)*legend_parameters[["between_line_height"]]*dy, 
           pch=16, col=scales::alpha(col_df$color,1), 
           cex=point_size*legend_parameters[["point_size_multiplier"]])
    # text
    text(x=legend_x1+legend_parameters[["line_text_space"]]*dx,
         y=y_max-0:(length(all_groups)-1)*legend_parameters[["between_line_height"]]*dy,
         labels=all_groups, adj=c(0,0.5), cex=legend_parameters[["text_size"]])
    
    # cg/cgi/gene/exon
    current_position_y <- y_max-length(all_groups)*legend_parameters[["between_line_height"]]*dy
    if(plot_cg){
      # we take actual positions of all CpGs of the same length as line width in legend
      cg_to_use <- (cg_all_coords>first_coord)&(cg_all_coords<=first_coord+legend_parameters[["line_width"]]*dx)
      cg_legend_positions <- cg_all_coords[cg_to_use]-first_coord
      # we want the first interval which has more than the average number of CpGs (this interval must exist)
      number_of_intervals <- ceiling((x_max-x_min)/(legend_parameters[["line_width"]]*dx))
      average_number_of_cg <- length(cg_all_coords)/number_of_intervals
      segment_index <- 1
      while(length(cg_legend_positions)<average_number_of_cg){
        cg_to_use <- (cg_all_coords>first_coord+segment_index*legend_parameters[["line_width"]]*dx)&(cg_all_coords<=first_coord+(segment_index+1)*legend_parameters[["line_width"]]*dx)
        cg_legend_positions <- cg_all_coords[cg_to_use]-(first_coord+segment_index*legend_parameters[["line_width"]]*dx)
        segment_index <- segment_index + 1
        # at some point an interval has to have at least mean number of CpGs within it
      }
      segments(x0=legend_x0+cg_legend_positions,
               y0=current_position_y+0.5*legend_parameters[["cg_height"]]*dy,
               y1=current_position_y-0.5*legend_parameters[["cg_height"]]*dy, 
               col=scales::alpha(col_cg, 1))
      # plot black line which represents x-axis
      segments(x0=legend_x0, x1=legend_x1, 
               y0=current_position_y-0.5*legend_parameters[["cg_height"]]*dy)
      #text
      text(x=legend_x1+legend_parameters[["line_text_space"]]*dx,
           y=current_position_y,
           labels="CpG sites", adj=c(0,0.5), cex=legend_parameters[["text_size"]])
      #
      current_position_y <- current_position_y-legend_parameters[["between_line_height"]]*dy
    }
    
    if(plot_cgi){
      polygon(x=legend_x0+c(0.05,0.95,0.95,0.05)*legend_parameters[["line_width"]]*dx,
              y=current_position_y+0.5*c(1,1,-1,-1)*legend_parameters[["polygon_height"]]*dy,
              border=scales::alpha(col_cgi,1), col=scales::alpha(col_cgi,1))
      # plot black line which represents x-axis
      segments(x0=legend_x0, x1=legend_x1, 
               y0=current_position_y-0.5*legend_parameters[["polygon_height"]]*dy)
      #text
      text(x=legend_x1+legend_parameters[["line_text_space"]]*dx,
           y=current_position_y,
           labels="CpG island", adj=c(0,0.5), cex=legend_parameters[["text_size"]])
      #
      current_position_y <- current_position_y-legend_parameters[["between_line_height"]]*dy
    }
    
    if(plot_gene|plot_exon){
      if(!plot_exon){
        segments(x0=legend_x0, x1=legend_x1, y0=current_position_y, col=scales::alpha(col_gene,1))
        #text
        text(x=legend_x1+legend_parameters[["line_text_space"]]*dx,
             y=current_position_y,
             labels="Gene", adj=c(0,0.5), cex=legend_parameters[["text_size"]])
      }else if(!plot_gene){
        polygon(x=c(legend_x0, legend_x1, legend_x1, legend_x0),
                y=current_position_y+0.5*c(1,1,-1,-1)*legend_parameters[["polygon_height"]]*dy,
                border=scales::alpha(col_cgi,1), col=scales::alpha(col_gene,1))
        #text
        text(x=legend_x1+legend_parameters[["line_text_space"]]*dx,
             y=current_position_y,
             labels="Exon", adj=c(0,0.5), cex=legend_parameters[["text_size"]])
        #
      }else{
        segments(x0=legend_x0, x1=legend_x1, 
                 y0=current_position_y+0.5*legend_parameters[["polygon_height"]]*dy, 
                 col=scales::alpha(col_gene,1))
        polygon(x=legend_x0+c(0.05,0.37,0.37,0.05)*legend_parameters[["line_width"]]*dx,
                y=current_position_y+0.5*c(1,1,-1,-1)*legend_parameters[["polygon_height"]]*dy,
                border=scales::alpha(col_gene,1), col=scales::alpha(col_gene,1))
        polygon(x=legend_x0+c(0.63,0.95,0.95,0.63)*legend_parameters[["line_width"]]*dx,
                y=current_position_y+0.5*c(1,1,-1,-1)*legend_parameters[["polygon_height"]]*dy,
                border=scales::alpha(col_gene,1), col=scales::alpha(col_gene,1))
        #text
        text(x=legend_x1+legend_parameters[["line_text_space"]]*dx,
             y=current_position_y,
             labels="Gene w/ exons", adj=c(0,0.5), cex=legend_parameters[["text_size"]])
      }
    }
  }
  
  # return names of plotted CpGs in the order of plotting
  return(cg_in_plot)
}
