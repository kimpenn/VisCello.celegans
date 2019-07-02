VisCello for C.elegans Embryogenesis
================
Qin Zhu, Kim Lab, University of Pennsylvania


About VisCello.celegans
------------------------

This is a tool for interactive visualization of C. elegans embryogenesis single cell data published with Packer, J.S., Zhu, Q., et al., 2019.

It has entire C. elegans embryogensis data **built in**, with features specifically designed for showing analysis results in the C. elegans paper, such as co-visualization of umap, lineage annotation and lineage tree.

For using VisCello for other single cell data visualization, please check https://github.com/qinzhu/VisCello.

Screenshot:

![](inst/app/www/screenshot.png?raw=true "VisCello screenshot")

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

Due to large data file currently being hosted on Git LFS, you cannot use devtools::install_github to install this package. 
Please follow protocol listed below to install:

0. Windows and Linux: install git from https://git-scm.com/book/en/v2/Getting-Started-Installing-Git. 

    MACOS and Linux: Install git-lfs from https://git-lfs.github.com/ (See FAQ for example installation using brew or conda).
    
    Any system: Install R (>3.4 required) and latest bioconductor (code below, copy and paste inside R):
    
    ```
    if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
    BiocManager::install()
    # For first time installer, if prompted "Update all/some/none? [a/s/n]", type a and press return.
    ```

1. In terminal copy paste followng line by line with return:

    ```
    git lfs clone https://github.com/qinzhu/VisCello.celegans.git
    R   #Launch R, Windows user open R/Rstudio, setwd() to parent folder of VisCello. 
    ```

2. Now inside R and do the following (line by line):

    ``` r
    install.packages("devtools") 
    devtools::install_local("VisCello.celegans")
    ```

    Now VisCello is ready to go! To launch Viscello, in R:

    ``` r
    library(VisCello.celegans)
    cello()
    ```
  
  If you have also installed VisCello, please make sure you do not load VisCello in the same session, as some functions are re-used and may lead to conflicts.

FAQ
-------------------------

* Q: I see error in R complaining: "Error in readRDS("data/eset.rds"): unknown input format", what should I do?
    
    A: Install git-lfs from Install git-lfs from https://git-lfs.github.com/ and then go to step 1. Example code using brew or conda:
    
    ```
    brew install git-lfs
    git lfs install
    git-lfs clone https://github.com/qinzhu/VisCello.celegans.git
    ```
    
    ```
    conda install -c conda-forge git-lfs
    git lfs install
    git-lfs clone https://github.com/qinzhu/VisCello.celegans.git
    ```

* Q: Can I use VisCello to explore my own data?
    
  A: Yes. Please use the general version: https://github.com/qinzhu/VisCello and kindly cite VisCello if you use it for publication.
  

Cite VisCello
-------------------------

Q. Zhu, J. I. Murray, K. Tan, J. Kim, qinzhu/VisCello: VisCello v1.0.0 (2019; https://zenodo.org/record/3262313).

Packer, J.S., Zhu, Q., Huynh, C., Sivaramakrishnan, P., Preston, E., Dueck, H., Stefanik, D., Tan, K., Trapnell, C., Kim, J. and Waterston, R.H., 2019. A lineage-resolved molecular atlas of C. elegans embryogenesis at single cell resolution. BioRxiv, p.565549.
