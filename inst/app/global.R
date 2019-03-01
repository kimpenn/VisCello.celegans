

if(Sys.getenv('SHINY_PORT') == "") options(shiny.maxRequestSize=10000*1024^2)

# For DEPLOY only

# deps<-c(
# "shinyBS",
# "shinycssloaders",
# "shinyWidgets",
# "htmlwidgets",
# "shinyBS",
# "R.utils",
# "purrr",
# "ggplot2",
# "RColorBrewer",
# "reshape2",
# "gridExtra",
# "shinythemes",
# "plotly",
# "heatmaply",
# "pheatmap",
# "DT",
# "dplyr",
# "openxlsx",
# "knitr",
# "monocle",
# "clusterProfiler",
# "celegans.db",
# "ggraph",
# "tidygraph",
# "data.table")
# 
# library(shinyBS);library(shinycssloaders);library(shinyWidgets);library(htmlwidgets);library(shinyBS);library(R.utils);library(purrr);library(ggplot2);library(RColorBrewer);library(reshape2);library(gridExtra);library(shinythemes);library(plotly);library(heatmaply);library(pheatmap);library(DT);library(dplyr);library(openxlsx);library(knitr);library(monocle);library(clusterProfiler);library(celegans.db);library(ggraph);library(tidygraph);library(data.table)

# for(d in deps) {
#     library(d,character.only = TRUE)
# }



Cello <- setClass("Cello",
                  slots = c(
                      name = "character", # The name of cvis
                      idx = "numeric", # The index of global cds object
                      proj = "list", # The projections as a list of data frame
                      pmeta = "data.frame", # The local meta data
                      notes = "character" # Other information to display to the user
                  )
)

# Load data
data_names=list.files('data/', pattern="*.rda", full.names = TRUE)
lapply(data_names,base::load, env = .GlobalEnv)

eset <- readRDS("data/eset.rds")
ct_tbl <-  readRDS("data/s6_tbl.rds")
lin_tbl <-  readRDS("data/s7_tbl.rds")
#lapply(list.files("src/", pattern = "\\.(r|R)$", recursive = F, full.names = TRUE), function(x){source(file = x)})


### Which meta data to show, and in what order ###
ctype_cols_advanced <- pmeta_attr$meta_id
names(ctype_cols_advanced) <- pmeta_attr$meta_name
ctype_cols_basic <- ctype_cols_advanced[c("Cell type (broad)", "Cell type + cell subtype", "Lineage", "Gene expression", "Embryo time")]
elin_cols_basic <- ctype_cols_advanced[c("Lineage", "150min early lineage","250min early lineage","Gene expression", "Embryo time")]
elin_cols_advanced <- c(elin_cols_basic, ctype_cols_advanced[!ctype_cols_advanced %in% elin_cols_basic])

bp_colorBy_choices <- ctype_cols_advanced[c("Cell type (broad)", "Cell subtype","Lineage", "Embryo time bin")]

de_meta_options <- ctype_cols_advanced[c("Cell type (broad)", "Cell subtype", "Lineage", "Cluster")]

numeric_palettes <- numeric_color_opt()
names(numeric_palettes) <- numeric_palettes

image_palettes <- numeric_palettes[numeric_palettes %in% c("RdYlBu", "RdBu", "viridis", "magma", "plasma", "inferno")]
heatmap_palettes <- image_palettes

image_colorBy_choices <- graph_genes

max_pc_show <- 10


