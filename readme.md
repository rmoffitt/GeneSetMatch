# Making a new R project in the Moffitt Lab? Use this repo as a template to ensure you have the right architecture!
### *Note: You should never need to pull this repository, leave it on github*

1. From the template_repo page, click "Use this template"
    * Alternatively, from the new repository page, click "From a template and select 'lthealy/template_repo'"
2. Name your repository
3. **SET IT TO PRIVATE**
3. Clone it to your Rserver, run devtools::document() to generate a namespace and .Rd files for prewritten functions
4. Open the DESCRIPTION file and change your package name to match your repository name
5. Build and Install and get going!

## General Layout 
* **data/** houses only .rda or .RData file types. These files are pre-loaded on build for easy access (think package::example)
* **R/** houses only .R files, be it wrapper functions or initial parsing scripts
* **inst/analysis/** houses .rmds for analysis
* **inst/extdata** houses raw data files (.csv,.txt) that are parsed via scripts in .R to generate the .RData objects
* **man/** houses .Rd files and will be auto-generated after running devtools::document(), this is where your "man"uals go for your locally written scripts (see included examples)
