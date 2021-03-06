% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meMultiplex_communities.R
\name{meMultiplex_communities}
\alias{meMultiplex_communities}
\title{Multiplex community structure}
\usage{
meMultiplex_communities(
  cor_layer,
  supplementary_layer,
  physical_nodes = FALSE,
  layer = 1,
  folder = "./meNet/",
  file_basename = "meNet_infomap",
  infomap_call = "infomap",
  cor_weighted = TRUE,
  supp_weighted = TRUE,
  cor_normalization_fun = max_normalization,
  supp_normalization_fun = neg_max_normalization,
  inter_cor_supp = NULL,
  inter_supp_cor = inter_cor_supp,
  infomap_seed = NULL,
  relaxation_rate = 0.15,
  delete_files = FALSE
)
}
\arguments{
\item{cor_layer}{Correlation layer of the multiplex as igraph object.}

\item{supplementary_layer}{Supplementary layer of the multiples as igraph
object.}

\item{physical_nodes}{Should physical or state nodes be used for community
structure nodes. Default to `FALSE`, i.e. to state nodes.}

\item{layer}{The layer from which the community structure is obtained.
Only used for state nodes.}

\item{folder}{Folder in which the files will be saved. It is automatically
created if not already present. Default value is `"./meNet/"`.}

\item{file_basename}{Base name of the creates files. For different files,
different suffix is added to the base name. 
Default value is `"meNet_infomap"`.}

\item{infomap_call}{Path to the `Infomap` on user's system.
Defaults to `"infomap"`.}

\item{cor_weighted}{Whether the correlation layer is weighted. Passed to
`meMultiplex`. Defaults to `TRUE`.}

\item{supp_weighted}{Whether the supplementary layer is weighted. Passed to
`meMultiplex`. Defaults to `TRUE`.}

\item{cor_normalization_fun}{Normalization function for the correlation layer.
Passed to `meMultiplex`. Defaults to `max_normalization`.}

\item{supp_normalization_fun}{Normalization function for the supplementary
layer. Passed to `meMultiplex`. Defaults to `neg_max_normalization`.}

\item{inter_cor_supp}{Weight of the inter-layer edges. Passed to
`meMultiplex`. Defaults to `NULL`.}

\item{inter_supp_cor}{Weight of the inter-layer edges. Passed to
`meMultiplex`. By default, the values is equal to `inter_cor_supp`.}

\item{infomap_seed}{`seed` parameter in `Infomap` call.
By default, no seed is set.}

\item{relaxation_rate}{`multilayer-relax-rate` parameter in `Infomap` call.
Default to `0.15`.}

\item{delete_files}{Should created files be automatically deleted
from the user's system.
Default to `FALSE`. Changing the parameter to `TRUE` should be done with 
caution since it will allow the function to delete files from user's system.}
}
\value{
Data frame with two columns. The first column contains names of nodes
and the second column contains corresponding `Infomap` community.
}
\description{
For two `igraph` layers with shared nodes, the function finds `Infomap`
communities of the constructed multiplex. Function interprets the first layer 
as correlation layer. The second layer, if weighted, should have only 
positive weights and can be interpreted as distance layer.

`Infomap` has to be locally installed for this function to work. 
Help with installation and alternative ways to run the `Infomap` algorithm 
are explained in 'Details'.
}
\details{
Function finds communities of nodes (CpGs) performing a sequence of three  
action. Firstly, `meMultiplex` function is called and the structure of
multiplex if written to a file. The file is used in a call of to
`Infomap` which performs the clustering of nodes. The clusters are saved to
a file. In the end, function `read_infomap_communities` reads the community
structure of the multiplex from the two files and returns a tidy data frame
with nodes and their corresponding communities.

For multiplex reconstruction, parameters `cor_layer`, `supplementary_layer`,
`cor_weighted`, `supp_weighted`, `cor_normalization_fun`, `supp_normalization_fun`,
`inter_cor_supp` and `inter_supp_cor` are passed to the `meMultiplex` with
`output_type` set to `"infomap"`.

`Infomap` is called from the command line using the path given in
`infomap_call` with flags `-i multilayer -o clu`.
If specified, `infomap_seed` is added as `--seed` flag and `relaxation-rate`
is added as `--multilayer-relax-rate` flag.
For help with the installation of the `Infomap` algorithm, visit
`Infomap` website: https://www.mapequation.org/infomap/.
If the installation fails, the output of the `meMultiplex` can be manually
uploaded online on the `Infomap` website and the resulting '.clu` files can 
be used as input to `read_infomap_community` function. For futher details,
refer to `Infomap` paper \insertCite{infomap}{meNet}.

To read multiplex communities, function `read_infomap_communities` is called
with the location of two files from the previous steps. Parameters 
`physical_nodes` and `layer` are passed to the function.

In the first two steps of the function algorithm, files are created and
saved in the `folder` using `file_basename`. These files are by-product of
the function and can be deleted if the user doesn't want to examine their
details. If `delete_files` is `TRUE`, these files will be automatically 
deleted from the user's system. If the folder contains only these files,
it will also be deleted.
}
\references{
\insertAllCited{}
}
