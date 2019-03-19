VisCello for C.elegans Embryogenesis
================
Qin Zhu, Kim Lab, University of Pennsylvania

Online version
------------------------

Link: https://cello.shinyapps.io/celegans/

Bugs you found with the online tool please post to this github repo.

Installation
--------------------------------------

Due to large data file currently being hosted on Git LFS, you cannot use devtools::install_github to install this package. 
Please follow protocol listed below to install:

0. Windows and Linux: install git from https://git-scm.com/book/en/v2/Getting-Started-Installing-Git. 

    MACOS and Linux: Install git-lfs from https://git-lfs.github.com/ (See FAQ for example installation using brew).
    
    Any system: Install R (>3.4 required) and latest bioconductor (code below, copy and paste inside R):
    
    ```
    if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
    BiocManager::install()
    ```

1. In terminal copy paste followng line by line with return:

    ```
    git lfs clone https://github.com/qinzhu/VisCello.git
    R   #Launch R, Windows user open R/Rstudio, setwd() to parent folder of VisCello. 
    ```

2. Now inside R and do the following (line by line):

    ``` r
    install.packages("devtools") 
    devtools::install_local("VisCello")
    ```

    Now VisCello is ready to go! To launch Viscello, in R:

    ``` r
    library(VisCello)
    cello()
    ```

FAQ
-------------------------

* Q: I see error in R complaining: "Error in readRDS("data/eset.rds"): unknown input format", what should I do?
    
    A: Install git-lfs from Install git-lfs from https://git-lfs.github.com/ and then go to step 1. Example code using brew:
    
    ```
    brew install git-lfs
    git lfs install
    git-lfs clone https://github.com/qinzhu/VisCello.git
    ```

Cite VisCello
-------------------------

Packer, J.S., Zhu, Q., Huynh, C., Sivaramakrishnan, P., Preston, E., Dueck, H., Stefanik, D., Tan, K., Trapnell, C., Kim, J. and Waterston, R.H., 2019. A lineage-resolved molecular atlas of C. elegans embryogenesis at single cell resolution. BioRxiv, p.565549.
