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

In terminal copy paste followng line by line with return:

```
git clone https://github.com/qinzhu/VisCello.git
R
```

Now inside R and do the following (line by line):

``` r
install.packages("devtools") 
devtools::install_local("VisCello")
```

Now VisCello is ready to go! To launch Viscello, in R:

``` r
library(VisCello)
cello()
```

Cite VisCello
-------------------------

Paper in submission at this moment.


