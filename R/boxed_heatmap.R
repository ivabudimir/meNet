# it doesn't call "layout" function so it can be used in the layout
#'@export
boxed_heatmap <- function(cor_matrix, title="", color_title="black", cex_title=1.3, margin_add=1.5, box_width=1, line_title=0.5){
  my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 200) #to have exact white
  cor_matrix <- t(apply(cor_matrix, 2, rev)) #rotate for 90 degrees in clockwise direction because image plots starting form bottom left corner
  par(mar=c(0.5,0.5,0.5,0.5)+margin_add)
  image(cor_matrix, axes=FALSE, col=my_palette, breaks=seq(-1,1,length=201))
  box(lwd=box_width)
  title(title, col.main=color_title, cex.main=cex_title, line=line_title)
}