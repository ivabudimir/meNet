% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_multiplex.R
\name{plot_multiplex}
\alias{plot_multiplex}
\title{Multiplex plot}
\usage{
plot_multiplex(
  first_layer,
  second_layer,
  layout_layer = "first",
  module_df = NULL,
  first_weighted = TRUE,
  second_weighted = TRUE,
  first_normalization_fun = .max_normalization,
  second_normalization_fun = .neg_max_normalization,
  include_negative = TRUE,
  node_sort_fun = .coord_sort,
  layout_fun = igraph::layout_with_fr,
  layout_weighted = FALSE,
  layout_weight_transformation_fun = .max_normalization,
  main_title = NULL,
  first_title = "Correlation",
  second_title = "CpG island",
  node_col = "black",
  node_size = 2,
  pos_edge_col = "grey70",
  neg_edge_col = "blue",
  plot_module = FALSE,
  module_border_col = "black",
  module_border_expand = 10,
  plot_layout_graph = FALSE
)
}
\arguments{
\item{first_layer}{The first layer of the multiplex as `igraph` object.}

\item{second_layer}{The second layer of the multiplex as `igraph` object.}

\item{layout_layer}{Defines the plotting method used to determines the layout
of nodes. Can take values `"first"`, `"second"`, `"module"`, `"aggregated"` or 
`"coord"`. See 'Details'.}

\item{module_df}{Data frame with community structure. The first column 
contains names of nodes and the second column contains corresponding
community. It is mandatory parameter if `layout_layer="module"` or
`plot_module=TRUE`.}

\item{first_weighted}{Whether the first layer is weighted. If `FALSE`,
all weights are set to `1`. Default to `TRUE`.}

\item{second_weighted}{Whether the second layer is weighted. If `FALSE`,
all weights are set to `1`. Default to `TRUE`.}

\item{first_normalization_fun}{Normalization function applied on the weights
of the first layer, if the layer is weighted. If `NULL`, no normalization
is applied. Defaults to `max_normalization`.}

\item{second_normalization_fun}{Normalization function applied on the weights
of the second layer, if the layer is weighted. If `NULL`, no normalization
is applied. Defaults to `neg_max_normalization`.}

\item{node_sort_fun}{Sort function used on the nodes of the network.
Default to `.coord_sort`, which sorts CpGs based on their genomic location.}

\item{layout_fun}{Algorithm for the layout of the nodes used on 
the layout graph.
In choosing of the function, attention should be paid whether weights of 
the layout graph represent "strength" or "cost".
Default to `igraph::layout_with_fr`.}

\item{layout_weighted}{Whether the edge weights of the layout graph are
used by `layout_fun` to determine the layout of nodes.
Default to `FALSE`.}

\item{layout_weight_transformation_fun}{Normalization function applied on the 
weights of the layout graph, if the graph is weighted. If `NULL`, no normalization
is applied. Defaults to `max_normalization`.}

\item{main_title}{Main title of the plot. By default, there is no main title.}

\item{first_title}{Title of the first subplot with the first layer of 
multiplex. Default value is `"Correlation"`.}

\item{second_title}{Title of the second subplot with the second layer of 
multiplex. Default value is `"CpG island"`.}

\item{node_col}{Color of the nodes. Color name if all nodes have the same
color or data frame with node names in the first column and corresponding 
color in the second column. If some nodes are missing from the given data
frame, their color is set to `"black"`.
Defaults value of `node_col` is `"black"`.}

\item{node_size}{Size of the nodes. Defaults to `2`.}

\item{pos_edge_col}{Color of positive edges. Default value is `"grey70"`.}

\item{neg_edge_col}{Color of negative edges. Default value is `"blue"`.}

\item{plot_module}{Whether to plot the border around the communities defined
in `module_df`.
Defaults to `FALSE`.}

\item{module_border_col}{Color of the border around the communities. Passed
to `plot.igraph` as `mark.border` parameter.
Defaults to `"black"`.}

\item{module_border_expand}{Size of the border around the communities. Passed
to `plot.igraph` as `mark.expand` parameter.
Defaults to `10`.}

\item{plot_layout_graph}{Whether to add a plot of layout graph. If `TRUE`, third
plot is added to the right of the multiplex with the same layout of nodes and
with edges from the layout graph. Default to `FALSE`.}
}
\value{
A vector of removed nodes.
}
\description{
Plots a multiplex network from two `igraph` layers. It is a 
corresponding plotting function for `meMultiplex` function
which plots two layers next to each other preserving the location of 
shared nodes.
As in `meMultiplex`, the first layer can contain negative edge weights while 
the second layer must have only non-negative weights.
}
\details{
Plots two-layer multiplex from two `igraph` layers. Layers are plotted
next to each other with the same layout of nodes. Plotted edges represent
edges found in a certain layer of multiplex. Parameters `main_title`,
`first_title`, `second_title`, `node_col`, `node_size`, `pos_edge_col` and
`neg_edge_col` help in customization of the plot.

Nodes which are isolated in both layers are removed from the multiplex. Nodes
which are found only in one layer are added to the other layer as isolated 
nodes. This equalizes the set of nodes between two layers.

If `module_df` is given, borders can be plotted around the communities of
nodes setting `plot_module` to `TRUE`.
If the same node belong to multiple communities, only the first 
membership will be used in plotting. For customization, parameters
`module_border_col` and `module_border_expand` can be used.

To determine the layout of nodes, new layout graph is constructed and used
in the `layout_fun`. The nodes of the layout graph are the same as in
the first and the second layer. 
If `layout_weighted`, the weight of edges are used by
`layout_fun`. There are several option for the construction of the layout
graph, depending on the value of `layout_layer`. If `layout_layer` has value
\itemize{
   \item{`"first"`}{layout graph is the first layer with removed negative
   edges;}
   \item{`"second"`}{Layout graph is the second layer.}
   \item{`"module"`}{Layout graph is constructed from `module_df`. In layout
   graph edges exist between nodes which are found in the same community.
   Weight of an edge is the number of communities in which a selected pair
   of nodes is found. Weight can be larger than `1` only if the same node
   can be a member of multiple communities in `module_df`.}
   \item{`"aggregated"`}{Layout graph is obtained aggregating the first and
   the second layer. Weight of an edge in layout graph is the sum of weights
   of the same edge in the first and the second layer.}
   \item{`"coord"`}{Layout graph is disconnected graph with no edges. This
   only makes sense if `layout_fun` uses circular layout since then
   nodes will be ordered by `node_sort_fun`.}
}
Structure of the layout graph can be plotted as third layer setting
`plot_layout_graph` to `TRUE`. It if not advisable to use this in the final
plot, but this feature can be useful to find the appropriate `layout_fun`.

If layers are weighted, function searches for the edge attributes with
names`"Cor"` and `"Dist"`. If the attribute doesn't exist, it takes the first
edge attribute instead. The function handles weight of 
the first layer as correlations and weights of the second layer as distances. 
Function checks if all edge weights of `second_layer` are positive.
}
