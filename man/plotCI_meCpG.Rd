% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotCI_meCpG.R
\name{plotCI_meCpG}
\alias{plotCI_meCpG}
\title{Beta value error plot with CpG coordinates}
\usage{
plotCI_meCpG(
  beta_values,
  sample_groups,
  cg_list = NULL,
  chr = NULL,
  first_coord = NULL,
  last_coord = NULL,
  add_lines = TRUE,
  plot_cg = FALSE,
  plot_cgi = FALSE,
  plot_gene = FALSE,
  plot_exon = FALSE,
  transcript_types = c("NM", "NR"),
  col_groups = NULL,
  col_cg = "grey77",
  col_gene = "cyan4",
  col_cgi = "grey77",
  title = "",
  x_label = "",
  y_label = "",
  text_cex = 1,
  beta_min = 0,
  beta_max = 1,
  error_size = 1,
  cg_names = FALSE,
  cg_names_expand = NA,
  plot_legend = FALSE,
  legend_cex = 1,
  plot_outside = TRUE
)
}
\arguments{
\item{beta_values}{Data frame with methylation beta values. CpG sites should
be listed in the rows with the row names being the Illumina CpG identifiers.
Samples should be listed in the columns.}

\item{sample_groups}{Data frame which defines the sample groups. It has to 
have two columns, the first column with sample IDs and the second column 
with the corresponding group.}

\item{cg_list}{List of CpG sites whose methylation values should be presented
in the plot. If provided, CpGs will be plotted equidistantly on the x-axis in
the given order.}

\item{chr}{Chromosome on which the plotted CpG sites are found. In the format
`1`:`22` or `"X"` or `"Y"`.}

\item{first_coord}{The first coordinate to be represented in the plot.}

\item{last_coord}{The last coordinate to be represented in the plot.}

\item{add_lines}{Whether to connect the error bars. If `TRUE`, mean value
points will be connected with the lines. Default to `TRUE`.}

\item{plot_cg}{Whether to plot the position of all CpG sites on x-axis.
Defaults to `FALSE`.}

\item{plot_cgi}{Whether to plot the position of CpG island on x-axis.
Defaults to `FALSE`.}

\item{plot_gene}{Whether to plot the position of genes/transcripts on x-axis.
Defaults to `FALSE`.}

\item{plot_exon}{Whether to plot the position of exons on x-axis. Defaults to
`FALSE`.}

\item{transcript_types}{Genes/transcripts to be plotted. `"NM"` stands for
protein-coding and `"NR"` stands for non-protein coding transcripts.
By default, both types of transcripts are included.}

\item{col_groups}{Color used for the error bars. It can be given as a list
of colors of the same length as the number of groups. Alternatively, it
can be given as a data frame with sample groups in the first column and
the corresponding color in the second column.}

\item{col_cg}{Color of the other CpG sites. Used only if `plot_cg` is `TRUE`.
Default color is `grey77`.}

\item{col_gene}{Color of the genes/exons. Used only if `plot_genes` or
`plot_exons` is `TRUE`. Default color is `"cyan4"`.}

\item{col_cgi}{Color of the CpG islands. Used only if `plot_cgi` is `TRUE`.
Default color is `"grey77"`.}

\item{title}{Plot title. By default, omitted.}

\item{x_label}{X-axis label. By default, omitted.}

\item{y_label}{Y-axis label. By default, omitted.}

\item{text_cex}{Numeric character expansion factor, used for all displayed 
text. Defaults to `1`.}

\item{beta_min}{Minimum beta value to be displayed. Has to be in `[0,1>` range.
Defaults to `0`.}

\item{beta_max}{Maximum beta value to be displayed. Has to be in `<0,1]` range.
Defaults to `1`.}

\item{error_size}{Number of standard deviations which are to be plotted
below and above the mean. Defaults to `1`.}

\item{cg_names}{Whether to add names of the plotted CpG sites below the
x-axis. Defaults to `FALSE`.}

\item{cg_names_expand}{Expand factor for the `cg_names`. If the labels are 
cluttered, the expand factor could be set to separate the labels. Expand 
factor defines the minimum distance between two consecutive labels and is
expressed as a percentage of the distance between the first and the last
label. Optimal separation is usually achieved for values close to `1`.
For the default `NA` value, labels are not expanded.}

\item{plot_legend}{Whether to plot the legend. Default to `FALSE`.}

\item{legend_cex}{Numeric expansion factor, used for all elements of the
legend. Defaults to `1`.}

\item{plot_outside}{Whether error plot should be plotted outside of the
plotting area. Plotting area is determined with `beta_min` and `beta_max`. 
Defaults to `TRUE`.}
}
\value{
Vector with names of plotted CpG sites, in the same order as in the
plot.
}
\description{
Plots methylation beta values in the form of error bars. For
every group of samples, the mean methylation value is plotted together
with the confidence interval. Depending on the function parameters, CpG sites 
represented in the plot can be plotted equidistantly or their position on
the x-axis can represent their true genomic coordinates. In the latter case,
additional information can be included in the plot such as position of genes
and exons.
}
\details{
Plots error bars and confidence intervals for beta values of given 
CpGs. If a list of CpGs is given, they are plotted equidistantly. 
If a range of coordinates and a chromosome are given, then all CpGs inside 
the range are plotted and positioned on x-axis based on their chromosomal 
coordinate (from 5' end to 3' end). 
If the latter case, additional information can be displayed on
the x-axis such as positions of unmeasured CpGs, positions of CpG islands
and/or positions of genes and exons. If we are specifically interested in
protein-coding or non-protein-coding genes, we can accordingly set the values
of parameter `transcript_types`. The color for each genomic element can be 
changed with `col_cg`, `col_cgi` and `col_gene` parameters.

For error bars, the default plotted confidence interval represents distance of
one standard deviation from the mean. This can be changed with `error_size` 
parameter. If we are interested only in certain beta values, we can restrict the
y-axis with `beta_min` and `beta_max` parameters.
Error bars are separately plotted for each group of samples and the color of each
group is set with `col_groups` parameter. Groups of samples are defined in
`sample_groups`. To better observe the methylation patterns, lines which 
connect mean methylation values can be added with `add_lines`.

Plotting area can be customized with `title`, `x_label` and `y_label`. 
Below the x-axis, names of CpG sites can be displayed setting `cg_names` to
`TRUE`. In cg names are cluttered, we can expand them setting `cg_names_expand`
factor to an appropriate numerical value. Usually this value is a number close
to `1`. An increase of the expand factor will result in more separate labels.
For the optimal separation, the best strategy is to try several different
values until the right value is found. The size of all displayed text is 
regulated simultaneously with `text_cex` parameter.

The legend can be added to the plot with `plot_legend` parameter. There are 
two parts of the legend. Color legend for sample groups and the legend which
explains the meaning of plotted genomic elements. The size of the whole legend
is regulated with `legend_cex` parameter.

In some cases, error bars have values outside the `beta_min`-`beta_max` range.
If `plot_outside` is `FALSE`, these values will be trimmed. Otherwise, they will
be plotted. If `beta_min` is `0` and `beta_max` is `1`, there still may be values
outside of the range. This happens for methylation values which are very close 
to `0` or `1` since the confidence interval is wider than the real methylation
values. In this case, it is recommended for `plot_outside` to be set to `TRUE`.
}
