% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_infomap_communities.R
\name{read_infomap_communities}
\alias{read_infomap_communities}
\title{Multiplex community structure from file}
\usage{
read_infomap_communities(
  multiplex_file,
  clu_file,
  physical_nodes = FALSE,
  layer = NULL
)
}
\arguments{
\item{multiplex_file}{File with the `Infomap` structure of the multiplex.}

\item{clu_file}{File with the communities obtained as output of the 
`Infomap` algorithm.}

\item{physical_nodes}{Whether the `clu_file` contains physical of state
nodes. Default to `FALSE`.}

\item{layer}{The layer from which the community structure is obtained.
Only used if `clu_file` contains state nodes.}
}
\value{
Data frame with two columns. The first column contains names of nodes
and the second column contains corresponding `Infomap` community.
}
\description{
Reads `Infomap` community structure for multiplex from two files.
The first file contains `Infomap`-style multiplex structure and the second
file is an output from the `Infomap` algorithm containing the communities.
Function his returns data frame with nodes and their corresponding community.
}
\details{
The first file can be obtained with `meMultiplex` function with the
output type set to `"Infomap"`. The second file is a result of `Infomap` 
algorithm for community detection.
`Infomap` algorithm finds community structure of the multiplex using
`multiplex_file` and outputs it in `clu_file`. 
This function unites the two files in one tidy data frame with multiplex 
nodes and their corresponding community.

`Infomap` returns two different types of communities, communities based on 
physical nodes (normal nodes) and communities based on state nodes. State 
nodes are abstract nodes created by `Infomap` to label physical nodes which 
exist in multiple layers. For example, for a physical node which is shared 
between two layers, two state nodes are created. Thus community structure of 
the state nodes contains two community structures, one existing in each layer 
of multiplex and both od them influenced by the whole multiplex structure.
If `physical_nodes` is `FALSE`, `clu_file` must contain state node output.
In this case `layer` must be specified and the community structure is taken
from that layer.
If `physcal_nodes` is `TRUE`, `clu_file` must contain physical node output.
This output overlaps two state node community structures existing in different
layers and thus the same physical node can exist in multiple communities.
For more information, visit `Infomap` website \insertCite{infomap}{meNet}.
}
\references{
\insertAllCited{}
}
