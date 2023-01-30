# GeneSetMatch

## General Layout 
* **data/** houses only .rda or .RData file types. These files are pre-loaded on build for easy access (think package::example)
* **R/** houses only .R files, be it wrapper functions or initial parsing scripts
* **inst/analysis/** houses .rmds for analysis
* **inst/extdata** houses raw data files (.csv,.txt) that are parsed via scripts in .R to generate the .RData objects
* **man/** houses .Rd files and will be auto-generated after running devtools::document(), this is where your "man"uals go for your locally written scripts (see included examples)
