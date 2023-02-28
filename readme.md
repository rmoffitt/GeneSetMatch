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

To allow smooth implementation, please parse your input data into the following format: a list containing three elements;   *a) a numeric matrix of raw expression data where rows are indexed by gene IDs and columns are indexed by experimental samples
  *b) a data frame with information corresponding to all samples, biological replicates, or any other relevant data to the downstream analysis
  *c) a data frame ideally containing ENSEMBL, SYMBOL or ENTREZ IDs which correspond to the indexed gene data in the expression matrix. If some gene IDs types are not available, the GeneSetMatch package offers a variety of conversion functions that aid in interspecies gene ID conversions, as well as animal-model-to-human homologues.

<img src="images/Data Structure.png" width="100%"/>



### Differential Expression Analysis 

tools for DE (image)
example script (image)



### Dimentional Reduction Analysis
NMF
example script (image)



### Tranfromation and filtering


### Pathway Enrichment Analysis
make sure it's the correct species when maing gmt files
example script (image)
talk about structure (image)



### Visualizations



#### DE Heatmap

#### Multimodal Heatmap




## License

All source code in this repository, with extensions `.R` and `.Rmd`, is available under an MIT license. All are welcome to use and repurpose this code.

**The MIT License (MIT)**

Copyright (c) 2016 Richard A. Moffitt

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
