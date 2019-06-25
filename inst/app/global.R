

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
clist <- readRDS("data/clist.rds")
elist <- readRDS("data/elist.rds")
ct_tbl <-  readRDS("data/s6_tbl.rds")
lin_tbl <-  readRDS("data/s7_tbl.rds")
tree_tbl <- as_tibble(readRDS("data/lineage_tree_tbl.rds"))
lin_sc_expr <- readRDS("data/lin_sc_expr_190602.rds")
lin.expanded.list <- readRDS("data/lin_expanded_list_0602.rds")
avail_nodes <- readRDS("data/avail_nodes.rds")
cell_type_markers <- read.xlsx("data/Supplementary_Tables_190611.xlsx",sheet=1, startRow=4)
lineage_markers <-  read.xlsx("data/Supplementary_Tables_190611.xlsx",sheet=4, startRow=7)
lineage_markers$IsTerminal <- NULL; lineage_markers$Redundant<-NULL
#lapply(list.files("src/", pattern = "\\.(r|R)$", recursive = F, full.names = TRUE), function(x){source(file = x)})
gene_tbl <- data.table("Gene names:" = rownames(eset))


# "Meta data" of cell meta data, decide which meta to show, which color to apply to different meta

colorby_order <- c(
    "Cell type/subtype" = "plot.cell.type",
    "Cell type (broad)" = "cell.type",
    "Cell subtype" = "cell.subtype",
    "Lineage" = "lineage",
    
    "Gene expression" = "gene.expr", # This will not be in any meta column
    "UMI count" = "n.umi",
    "Number of expressed Genes" = "num.genes.expressed",
    "Size factor" = "Size_Factor",
    "Likely doublets/debris" = "to.filter",
    
    "Batch" = "batch",
    "Time point" = "time.point",
    "Embryo time" = "raw.embryo.time",
    "Embryo time bin" = "raw.embryo.time.bin",
    
    "Embryo time (smoothed)" = "embryo.time",
    "Embryo time bin (smoothed)" = "embryo.time.bin",
    "Cluster" = "Cluster" # Note only this meta data is taken from local cvis object, 'Cluster' is now allowed in global meta colnames.
)

pmeta_attr <- data.frame(meta_id = colorby_order, meta_name = names(colorby_order), stringsAsFactors=FALSE)

pmeta_attr$is_numeric <- sapply(as.character(pmeta_attr$meta_id), function(x) {
    if(x %in% colnames(pData(eset))) {
        is.numeric(pData(eset)[[x]])
    } else if(x %in% colnames(clist[[1]]@pmeta)) {
        is.numeric(clist[[1]]@pmeta[[x]])
    } else if(x == "gene.expr") {
        T
    } else {
        NA
    }
})


# optional, default_pal of each meta
pmeta_attr$dpal <- ifelse(pmeta_attr$is_numeric, "viridis", "Set1")
pmeta_attr$dpal[grepl("time", pmeta_attr$meta_id)] <- "BlueGreenRed"
pmeta_attr$dscale <- ifelse(pmeta_attr$is_numeric, "log10", NA)
pmeta_attr$dscale[which(pmeta_attr$meta_id %in% c("Size_Factor", "raw.embryo.time", "embryo.time"))] <- "identity"

### Which meta data to show, and in what order ###
ctype_cols_advanced <- pmeta_attr$meta_id
names(ctype_cols_advanced) <- pmeta_attr$meta_name
ctype_cols_basic <- ctype_cols_advanced[c("Cell type/subtype", "Cell type (broad)", "Lineage", "Gene expression", "Embryo time")]
elin_cols_basic <- ctype_cols_advanced[c("Lineage","Cell type/subtype", "Cell type (broad)", "Gene expression", "Embryo time")]
elin_cols_advanced <- c(elin_cols_basic, ctype_cols_advanced[!ctype_cols_advanced %in% elin_cols_basic])
bp_colorBy_choices <- ctype_cols_advanced[c("Cell type/subtype", "Cell type (broad)", "Cell subtype","Lineage", "Embryo time bin")]
de_meta_options <- ctype_cols_advanced[c("Cell type/subtype", "Cell type (broad)", "Cell subtype", "Lineage", "Cluster")]

# Lineage subset default show
elin_sets_basic <- c(names(elist)[c(1:11)])

numeric_palettes <- numeric_color_opt()
names(numeric_palettes) <- numeric_palettes

image_palettes <- numeric_palettes[numeric_palettes %in% c("RdYlBu", "RdBu", "viridis", "magma", "plasma", "inferno")]
heatmap_palettes <- image_palettes




image_colorBy_choices <- graph_genes

max_pc_show <- 10


