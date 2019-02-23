VisCello for C.elegans Embryogenesis
================
Qin Zhu, Kim Lab, University of Pennsylvania

Installation
------------------------

Compile (install) the VisCello package
--------------------------------------

Open R and do the following:

``` r
setwd("/Users/yourname/Download/") # Replace the path with path to PARENT folder of VisCello
install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("monocle", update = F, ask=F)
BiocManager::install("clusterProfiler", update = F, ask=F)
BiocManager::install("org.Mm.eg.db", update = F, ask=F)
BiocManager::install("org.Hs.eg.db", update = F, ask=F)
devtools::install_local("VisCello.base")
```

Now VisCello is ready to go! To launch Viscello, in R:

``` r
library(VisCello.base)
viscello()
```

Host VisCello on a server
-------------------------

To host VisCello on a server, you need either ShinyServer (<https://www.rstudio.com/products/shiny/shiny-server/>) or use the shinyapps.io service (<https://www.shinyapps.io/>). Start Viscello on server by calling above code.

Reference
---------

Paul, Franziska, Ya’ara Arkin, Amir Giladi, Diego Adhemar Jaitin, Ephraim Kenigsberg, Hadas Keren-Shaul, Deborah Winter, et al. 2015. “Transcriptional Heterogeneity and Lineage Commitment in Myeloid Progenitors.” *Cell* 163 (7). Elsevier: 1663–77.
