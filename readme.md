# GeneSetMatch

## About
GeneSetMatch is a novel methodology and associated software tool which aims to streamline a variety of statistical analyses for large -omic data sets. GeneSetMatch allows the user to access a menu of appropriate differential expression analyses, dimensionality reductions, transformations and a novel informative visualization of the output in a single, user-friendly software tool.

## Repository layout

* Data/ holds several bulk and single0-cell RNA-Seq example datasets. These datasets are used for case-study analyses in the associated publication (link and citation forthcoming with publication).
* R/ holds all wrapper functions used in the application
* inst/analysis holds a variety of markdowns and scripts outlining relevant preliminary and in-depth example analytical protocols that GeneSetMatch can be utilized for. Sub-folders include high resolution .pdf files of the corresponding heatmaps generated by respective analyses.
* inst/shiny holds the code for execution of the app (in development)
* man/, as in all R packages, houses .Rd files and will be auto-generated after running devtools::document()

## Installing the package

To install the package, simply run
```
devtools::install_github("rmoffitt/GeneSetMatch")
```
## Workflow

GeneSetMatch offers a variety of analytical pipelines. The individual workflows will differ slightly based on the input data (bulk vs single-cell), desired statistical analyses (DE vs NMF), preparation of data for pathway enrichemnt step (trasformations & filtering) and subsequent visualization (appropriate heatmap selection). 

<img src="images/GSM Workflow new.jpg" width="100%"/>



### Data Structure

To allow smooth implementation, please parse your input data into the following format: a list containing three elements;
* a) a numeric matrix of raw expression data where rows are indexed by gene IDs and columns are indexed by experimental samples
* b) a data frame with information corresponding to all samples, biological replicates, or any other relevant data to the downstream analysis
* c) a data frame ideally containing ENSEMBL, SYMBOL or ENTREZ IDs which correspond to the indexed gene data in the expression matrix. If some gene IDs types are not available, the GeneSetMatch package offers a variety of conversion functions that aid in interspecies gene ID conversions, as well as animal-model-to-human homologues.

<img src="images/Data Structure.png" width="100%"/>



### Differential Expression Analysis 

GeneSetMatch offers an extensive menu of DE methods, all of which can be performed without additional manipulation assuming the input data frames is structured appropriately. Available methods exist in the “All_DE_tools_for_concordance_analysis.R” script and will be automatically loaded in with the build of the GeneSetMatch package. Each DE tool expects the same input to perform its analysis;
* - count data (i.e ‘full_data$ex’)
* - biological condition (i.e ‘full_data$sampInfo$condition’).

Each method will return a data frame with gene IDs as row labels and a variety of statistics as column labels. The user has the option to choose whichever statistic is of most interest, perhaps run several DE analyses in parallel and combine the results into a comprehensive matrix. This step allows for significant flexibility, and as long as the resulting data frame is shaped in a ‘gene vs statistic’ matrix, all of the downstream tools will accept the input.  

#### Example code for two DeSEQ2 computations on a bulk RNA-Seq dataset compared by KO status
<img src="images/Example DE.jpg" width="100%"/>

#### Resulting dataframe of the example code
<img src="images/DE dataframe.jpg" width="100%"/>



### Dimentional Reduction Analysis

Because bulk RNA-Seq may contain a combination of signals, deconvolution of the data can aid in subtype discovery. Since gene expression data is innately non-negative, non-negative matrix factorization is a compatible choice. Adapted from Brunet et al, we implemented NMF - an unsupervised clustering method which generates a predetermined number of metagenes and metasamples ready for deeper investigation. For purposes of the downstream GSM visualization, the k parameter must be set to 3.

#### Example code for rank = 3 NMF deconvolutionan of a bulk RNA-Seq rat expression matrix
<img src="images/Example NMF.jpg" width="100%"/>

#### Resulting unmixed matrix
<img src="images/NMF LLB.jpg" width="100%"/>



### Pathway Enrichment Analysis

For pathway enrichment step, we designed a modified fgsea function "GSMGSEA.R" which accepts the outputs of preceding analyses. For this step to run smoothly, the function expects
* theGoFile - a .gmt file containing a reference geneset or database. Unless declared otherwise, GeneSetMatch will generate .gmt files for the entirety of MSigDB via the “convert_msigdb_obj_to_gmt.R” function. Please make sure the .gmt file is made for the correct model species. 
* gene_list - a named vector of genes rank ordered by their corresponding numeric values. Gene names must be ENTREZ IDs.
Of important note - appropriate transformations need to be applied to NMF vectors ("transformedPM.R") and single-cell RNA-Seq lists ("FindAllMarkers.R") to ensure compatibility with GSMGSEA. Please refer to ins/analysis/LLB or inst/analysis/TAB MUR respectively for example Markdowns outlining the necessary steps.

The GSMGSEA function automatically structures the results into a layered list, designed to be fed directly into the heatmap functions without any further manipulation. 



### Visualizations

Deeply informative visualizations of GSEA outputs remain challenging. To approach this problem, we propose two novel heatmaps, one for results of a DE analysis and another for visualization of dimensionality reduction analysis. The results is a clustered heatmap that highlights individual genes in the context of their relevant gene sets. Importantly, our methods convey both the magnitude and direction of differential changes and display a high-level view of the degree of redundancy or heterogeneity of the GSEA result. 
Both GeneSetMatch heatmap functions accept the direct output of the previous pathway enrichment analysis step; 
* GSM.Heatmap(gsmgsea_result)
* GSM.Multilevel.Heatmap(gsmgsea_result)


#### DE Heatmap

In the simpler case of differential expression, where genes ‘go’ either up or down, the leading-edge genes are mapped back to their respective pathways in a tabular matrix format. For each .gmt file, the top 30 enriched pathways and their top 30 respective leading-edge genes are constructed into a matrix X, populated by the original DE statistic. This matrix is pseudo-colored in a red-to-blue fashion indicating the direction of differential expression.

##### Example of a DE Heatmap for KO vs WT analysis on bulk RNA-Seq PDAC Mouse data
<img src="images/Ex DE Heatmap.jpg" width="100%"/>

#### Multimodal Heatmap

In the more complex case of a deconvolution analysis, we generate our multi-modal heatmaps by consolidating three sets of GSMGSEA results as generated from each NMF factor vector. The result is a multi-colored heatmap in which the y-axis is populated by enriched gene sets, and the x-axis is populated by leading edge genes as identified by GSMGSEA. The leading edge genes from each sample vector are clustered together by pathway membership (ES), color coded by sample (R/G/B respectively for each NMF vector) and levels of normalized expression (an a-priori normalized statistic generated by the NMF deconvolution step), and clustered by our distance metric. Names of leading-edge genes are translated from ENTREZIDs to SYMBOL for easier interpretation.

##### Example of a multimodal heatmap for an NMF deconvolved mouse liver/lung/brain single-cell RNA-Seq dataset
<img src="images/Ex NMF heatmap.jpg" width="100%"/>


## License

All source code in this repository, with extensions `.R` and `.Rmd`, is available under an MIT license. All are welcome to use and repurpose this code.

**The MIT License (MIT)**

Copyright (c) 2016 Richard A. Moffitt

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
