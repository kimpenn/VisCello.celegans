VisCello for C.elegans Embryogenesis
================
Qin Zhu, Kim Lab, University of Pennsylvania


About VisCello.celegans
------------------------

This is a tool for interactive visualization of C. elegans embryogenesis single cell data published with Packer, J.S., Zhu, Q., et al., 2019.

It has entire C. elegans embryogensis data **built in**, with features specifically designed for showing analysis results in the C. elegans paper, such as co-visualization of umap, lineage annotation and lineage tree.

For using VisCello for other single cell data visualization, please check https://github.com/qinzhu/VisCello.

Screenshot:

[![Alt text](inst/app/www/screenshot.png?raw=true "VisCello screenshot")](https://cello.shinyapps.io/celegans/)

Online version
------------------------

Link: https://cello.shinyapps.io/celegans/

Bugs you found with the online tool please post to this github repo.

See also:
--------------------------------------

* C. elegans L2 data: https://github.com/qinzhu/Celegans.L2.Cello

* C. elegans Tintori et al. data (up to 16 cell stage): https://github.com/qinzhu/Celegans.Tintori.Cello

Installation
--------------------------------------

Please follow protocol listed below to install:

1. Get the package by downloading from Zenodo: https://zenodo.org/record/8397398 (NOTE we no longer provide full package via github). Unzip it (keep the folder structure so that the unzipped folder name is VisCello.celegans, and inside the folder you have folders including inst/, man/, R/).


2. Now launch R:

    ```
    R   #Launch R, Windows user open R/Rstudio, setwd() to parent folder of VisCello. 
    ```
    Run code below inside R:
    ``` r
    install.packages("devtools") 
    devtools::install_github("cysouw/qlcMatrix") # Due to recent removal of this package from CRAN
    devtools::install_local("VisCello.celegans", force=T)
    # Temporarily needed, install a version of "tidytree" package to avoid a bug in newer version
    packageurl <- "https://cran.r-project.org/src/contrib/Archive/tidytree/tidytree_0.2.4.tar.gz"
    install.packages(packageurl, repos=NULL, type="source")
    ```

    Now VisCello is ready to go! To launch Viscello, in R:

    ``` r
    library(VisCello.celegans)
    cello()
    ```
  

FAQ
-------------------------

* Q: Got installation error when running devtools::install_local: `package 'AnnotationDbi' was built under R version xxx, lazy loading failed for package xx`.

  A: Run following code in R:
  
    ```
    BiocManager::install(c("celegans.db", "GO.db", "DO.db"))
    install.packages("VisCello.celegans", repos = NULL, type = "source")
    ```
  

* Q: Can I use VisCello to explore my own data?
    
  A: Yes. Please use the general version: https://github.com/qinzhu/VisCello and kindly cite VisCello if you use it for publication.
  

Cite VisCello
-------------------------

Packer, J. S., Q. Zhu, C. Huynh, P. Sivaramakrishnan, E. Preston, H. Dueck, D. Stefanik, K. Tan, C. Trapnell, J. Kim, R. H. Waterston and J. I. Murray (2019). A lineage-resolved molecular atlas of C. elegans embryogenesis at single-cell resolution. Science: eaax1971.

Q. Zhu, J. I. Murray, K. Tan, J. Kim, qinzhu/VisCello: VisCello v1.0.0 (2019; https://zenodo.org/record/3262313).


