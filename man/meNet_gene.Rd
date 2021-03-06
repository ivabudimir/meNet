% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meNet_gene.R
\name{meNet_gene}
\alias{meNet_gene}
\title{CpG network: edges are based on gene affiliation}
\usage{
meNet_gene(
  cg_list = NULL,
  gene_list = NULL,
  weighted = TRUE,
  link_method = "twoLyr_clust",
  cor_matrix = NULL,
  data = NULL,
  cor_normalization_fun = max_normalization,
  dist_normalization_fun = neg_max_normalization,
  cor_threshold = 0.2,
  neg_cor_threshold = NULL,
  cor_stDev = NULL,
  cor_alpha = NULL,
  n_repetitions = 1000,
  alternative = "two_sided",
  infomap_call = "infomap",
  folder = "./meNet/",
  file_basename = "meNet_CGI_infomap",
  relaxation_rate = 0.15,
  cgGene_meta = data("CpG_genes", package = "meNet"),
  cgGene_meta_cols = list(cg_id = "IlmnID", cg_coord = "MAPINFO", gene_id =
    "UCSC_RefGene_Name", gene_region = "UCSC_RefGene_Group"),
  gene_regions = c("Promoter", "Body", "3'UTR"),
  expand_cg_list = FALSE,
  normalization_fun = NULL,
  save_all_files = FALSE,
  delete_files = FALSE
)
}
\arguments{
\item{cg_list}{}
}
\value{

}
\description{
For every present gene, function "meNet_singleGene" is called 
and the resulting network is a union of the single-gene networks. 
(Optional) weights are chromosomal distances.
}
\details{
if both cg_list and gene_list are given, we merge both lists;
    in the human genome" by J. Sandoval et al.
}
