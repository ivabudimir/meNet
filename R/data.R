#' Unsequences genome
#'
#' Data frame with unsequenced genomic regions for human genome assembly hg19.
#' Unsequences nucleotides are expressed as intervals since the unsequenced
#' regions span over larger portions of the genome.
#' 
#' @format A data frame with 362 rows and 4 variables:
#' \describe{
#'    \item{`chr`}{Chromosome containing unsequenced interval.}
#'    \item{`start`}{Chromosomal coordinate of the first nuclotide within the interval.}
#'    \item{`end`}{Chromosomal coordinate of the last nuclotide within the interval.}
#'    \item{`length`}{Length of the interval, expressed in base pairs.}
#' }.
#'
"N_ranges" 



#' Description of HM450K CpG sites
#'
#' Description file for CpGs  measured by Illumina 450K microarray.
#' It is Illumina `“HumanMethylation450 15017482 v1-2.csv”` manifest file with
#' reduced number of columns.
#'
#' @format A data frame with 485512 rows and variables:
#' \describe{
#'    \item{IlmnID}{Unique CpG identifier from the Illumina CG database.}
#'    \item{CHR}{Chromosome containing the CpG.}
#'    \item{MAPINFO}{Chromosomal coordinates of the CpG.}
#'    \item{Strand}{The Forward (F) or Reverse (R) designation of the Design 
#'    Strand.}
#'    \item{UCSC_RefGene_Name}{Target gene name(s), from the UCSC database. 
#'    *Note: multiple listings of the same gene name indicate splice variants.}
#'    \item{UCSC_RefGene_Accession}{The UCSC accession number(s) of the target 
#'    transcript(s). Accession numbers are given in the same order as the target
#'    gene transcripts.}
#'    \item{UCSC_RefGene_Group}{Gene region feature category describing the CpG 
#'    position, from UCSC. Features listed in the same order as the target gene 
#'    transcripts.}
#'    \item{UCSC_CpG_Islands_Name}{Chromosomal coordinates of the CpG island 
#'    from UCSC}
#'    \item{Relation_to_UCSC_CpG_Island}{The location of the CpG relative to the
#'    CpG island.}
#' }
#'
#'@source \url{https://support.illumina.com/downloads/humanmethylation450_15017482_v1-2_product_files.html}
#'
"CpG_annoHM450K"



#' Genes associated with HM450K CpG sites
#' 
#' Data frame based on the `CpG_anno450K`, Illumina Infinium HumanMethylation450
#' manifest file, describing the relationship between CpG sites and
#' genes/transcripts. If CpG is associated with multiple genes/transcripts,
#' it is present in multiple rows.
#' 
#' @format A data frame 687137 rows and 6 variables:
#' \describe{
#'    \item{IlmnID}{Unique CpG identifier from Illumina CG database.}
#'    \item{CHR}{Chromosome containing the CpG.}
#'    \item{MAPINFO}{Chromosomal coordinate of the CpG.}
#'    \item{UCSC_RefGene_Name}{Target gene name, from the UCSC database.}
#'    \item{UCSC_RefGene_Accession}{The UCSC accession number of the target transcript.}
#'    \item{UCSC_RefGene_Accession}{Gene region describing the CpG position.}
#' }
#'  
"CpG_genes"



#' UCSC CpG islands
#' 
#' Data frame with all UCSC CpG islands (CGI).
#' 
#' @format A data frame with 27718 rows and 6 variables:
#' \describe{
#'    \item{CpG_island_name}{Name of the CGI.}
#'    \item{chr}{Chromosome containing the GCI.}
#'    \item{start}{Chromosomal coordinate of the start of the CGI.}
#'    \item{end}{Chromosomal coordinate of the end of the CGI.}
#'    \item{length}{Length of the CGI, in base pairs.}
#'    \item{UCSC_CpG_Islands_Name}{Name of the CGI which matches 
#'    `UCSC_CpG_Islands_Name` used in Illumina 450K manifest file,
#'    or `meNet::CpG_anno450K` dataset.}
#' }
#' 
"UCSC_CGI"



#' UCSC transcripts
#' 
#' Data frame with description of both protein-coding and non-protein-coding 
#' transcripts including positions of exons. Transcript name starts with
#' `"NM"` for protein-coding and with `"NR"` for non-protein coding transcripts.
#' It is downloaded as `“UCSC RefSeq (refGene)”` table from the UCSC Table
#' Browser for assembly `“Feb. 2009 (GRCh37/hg19)”` and group 
#' `“Genes and Gene Predictions”`. 
#' 
#' @format A data frame with 78288 rows and variables:
#' \describe{
#'    \item{name}{Name of gene/transcript.}
#'    \item{chrom}{Chromosome containing the gene.}
#'    \item{strand}{Genome strand on which the gene is located.}
#'    \item{txStart}{Transcription start position (or end position for minus strand).}
#'    \item{txEnd}{Transcription end position (or start position for minus strand).}
#'    \item{cdsStart}{Coding region start position (or end position for minus strand).}
#'    \item{cdsEnd}{Coding region end position (or start position for minus strand).}
#'    \item{exonCount}{Number of exons.}
#'    \item{exonStarts}{Exon start positions (or end positions for minus strand).}
#'    \item{exonEnds}{Exon end position (or start positions for minus strand).}
#' }
#'
#' @source \url{http://genome.ucsc.edu/cgi-bin/hgTables}
"UCSC_genes"



#' Bacalini HM450K methylation data
#'
#' Beta values measured by Illumina HM450K microarray for 87 samples obtained by
#' Bacalini et al. \insertCite{bacalini2015identification}{meNet}. Data set has
#' family-based structure with the aim to study the effect of Down syndrome on 
#' the methylation signature. 29 Down syndrome patients were included in the study
#' together with a mother and one unaffected sibling. `meNet::Bacalini_betaValues`
#' data frame is a subset of the original data set with only 1000 CpG sites
#' from chromosome 21.
#' 
#' @format A data frame with 1000 rows and 87 columns. Row names are unique CpG
#' identifiers from the Illumina CG database. Column names are unique sample
#' identifiers further described in `meNet::Bacalini_sampleSheet`.
#'
#' @source Original data set was downloaded from the Gene Expression Omnibus (GEO)
#' data repository as a data set with accession number GSE52588
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52588}.
#' 
#' @references
#'       \insertAllCited{}
#'       
"Bacalini_betaValues"



#' Bacalini sample sheet
#' 
#' Description file for samples in Bacalini data set \insertCite{bacalini2015identification}{meNet}.
#' 87 samples have a family-based structure with a Down syndrome patient, a mother
#' and one unaffected sibling coming from each of 29 included families. Thus, each
#' sample belongs to one of three groups: Down syndrome patients (DSP), Down
#' syndrome mothers (DSM) or Down syndrome siblings (DSS). `meNet::Bacalini_sampleSheet`
#' is a subset of the original sample sheet which includes more details for every
#' samples.
#' 
#' @format A data frame with 87 rows and columns:
#' \describe{
#'    \item{gsm}{Unique identifier for every sample.}
#'    \item{Group}{A group to which the sample belongs. One of DSP, DSM or DSS.}
#' }
#'
#' @source Original sample sheet was downloaded from the Gene Expression Omnibus (GEO)
#' data repository as part of the data set with accession number GSE52588
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52588}.
#' 
#' @references
#'       \insertAllCited{}
#'       
"Bacalini_sampleSheet"



#' Johansson HM450K methylation data
#'
#' Methylation values based on the data set obtained by Johansson et al.
#' \insertCite{johansson2013continuous}{meNet}. Methylation beta values for 728
#' healthy samples were obtained with Illumina HM450K microarray. In the step of
#' data preparation and preprocessing, problematic CpGs were removed and the
#' original beta values were adjusted for age, sex and cell counts.
#' `meNet::Johansson_residuals` data frame contains the residuals of the 
#' preprocessing for 1000 CpG sites from chromosome 21.
#'
#' @format A data frame with 1000 rows and 728 columns. Row names are unique CpG
#' identifiers from the Illumina CG database.
#' 
#' @source Original data set was downloaded from the Gene Expression Omnibus (GEO)
#' data repository as a data set with accession number GSE87571
#' \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE87571}.
#' 
#' @references
#'       \insertAllCited{}
"Johansson_residuals"