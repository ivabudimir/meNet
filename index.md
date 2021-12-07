
<!-- README.md is generated from README.Rmd -->

## meNet

Cytosine methylation in the humane genome is an important and
well-studied epigenetic mark that has the potential to regulate gene
expression. Since the process of methylation predominantly occurs at CG
dinucleotide sequences, the target of studies are these so-called CpG
sites. One of the available tools, Infinium 450k assay measures the
level of methylation of \~450,000 CpG sites spread across the genome.
Studies suggest that CpG sites are not independent, but rather that the
methylation profile is guided by the complicated network structure of
CpG sites. **meNet** package utilizes prior biological knowledge and
correlation analysis for reconstruction of the methylation network.

## Installation

``` r
devtools::install_github("ivabudimir/meNet")
```

For installation of vignettes together with the package, add
“build\_vignettes=TRUE”.

## License

meNet is open source software, licensed under
[GPL-3](https://github.com/ivabudimir/visProteomics/blob/master/LICENSE).
