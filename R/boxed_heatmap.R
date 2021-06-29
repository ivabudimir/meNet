#' Heat map without the use of `layout` function
#'
#' @description Function which draws heat map with box around it and without any legend.
#' The function calls `image` and because of its simplicity it doesn't use `layout` so it can be used
#' inside the `layout` function together with other plots.
#'
#' @param cor_matrix Correlation matrix.
#' @param clrs_list A list of colors passed to `graphics::colorRampPalette`. 
#' If three colors are given, first corresponds to negative correlation, second 
#' to zero correlation and third to positive correlation.
#' Defaults to `c("blue", "white", "red")`.
#' @param title Title for the heat map. Passed to the `graphics::title` function.
#' @param color_title Color of the title. Passed to the `graphics::title` function.
#' @param cex_title Size of the title. Passed to the `graphics::title` function.
#' @param line_title Line at which title is added to the plot. Passed to the `graphics::title`
#' function.
#' @param margin_add Additional size for the margins. `c(0.5,0.5,0.5,0.5)+margin_add`
#' is passed as the `mar` parameter to `graphics::par`. Defaults to 1.5.
#' @param box_width Width of the box around the heat map. 
#' Passed to the `graphics::box` function.
#'
#' @import graphics
#'
#' @export
boxed_heatmap <- function(cor_matrix, clrs_list=c("blue", "white", "red"), title="", color_title="black", cex_title=1.3, line_title=0.5, margin_add=1.5, box_width=1){
  my_palette <- grDevices::colorRampPalette(c("blue", "white", "red"))(n = 200) #to have exact white
  cor_matrix <- t(apply(cor_matrix, 2, rev)) #rotate for 90 degrees in clockwise direction because image plots starting form bottom left corner
  par(mar=c(0.5,0.5,0.5,0.5)+margin_add)
  image(cor_matrix, axes=FALSE, col=my_palette, breaks=seq(-1,1,length=201))
  box(lwd=box_width)
  title(title, col.main=color_title, cex.main=cex_title, line=line_title)
}
