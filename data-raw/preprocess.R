




clist <- readRDS("data-raw/clist.rds")
elist <- readRDS("data-raw/elist.rds")
all_cds <- readRDS("data-raw/all_cds.rds")
g_all <- readRDS("data-raw/g_all.rds")
time_umap_list <- readRDS("data-raw/time_umap_list.rds")
time_deg_list <- readRDS("data-raw/time_deg_list.rds")
tf_tbl <- read.csv("data-raw/TFcisBP_dataCSV839.csv")
cell_type_markers <- read.csv("data-raw/Celltypemarkers.csv")
graph_genes <- readRDS("data-raw/graph_genes.rds")

names(clist[[1]]@proj)[2] <- "UMAP-2D [Paper]"

max_pc <- 10
clist_cvis <- lapply(1:length(clist), function(i) {
    x <- clist[[i]]
    cvis <- new("cvis")
    cvis@name <- names(clist)[i]
    cvis@idx <- x@idx
    cvis@proj <- x@proj
    # Also reduce PCA projection dimension to reduce size
    if(!is.null(x@proj[["PCA"]]))    cvis@proj[["PCA"]] <- x@proj[["PCA"]][,1:max_pc]
    local_pmeta <- data.frame(Cluster = x@cluster)
    rownames(local_pmeta) <- colnames(all_cds)[cvis@idx]
    #print(identical(rownames(local_pmeta), rownames(x@proj[[1]])))
    cvis@pmeta <- local_pmeta
    #cvis@notes <- "some notes"
    return(cvis)
})
names(clist_cvis) <- names(clist)
clist <- clist_cvis


names(clist)[which(names(clist) == "Early Embryo Germline Rectum")] <- "Early Embryo, Germline and Rectum"
names(clist)[which(names(clist) == "Non Ciliated Neurons")] <- "Non-Ciliated Neurons"
names(clist) <- tools::toTitleCase(names(clist))

max_pc <- 10
keep_elist <- c(
    "Time 150min Lineage"="time150_disp1_bc_regress", 
    "Time 200min Lineage"="time200_disp1_bc_regress", 
    "Time 250min Lineage"="time250_disp1_bc_regress", 
    "Time 300min Lineage"="time300_disp1_bc_regress",
    "Time 250min AB Neuron + Glia"="ABneuronglia[time250]", 
    "Time 250min AB Pharynx"="ABpharynx[time250]", 
    "Time 250min AB Hypodermis"="AB_HYP[time250]", 
    "Time 250min MSxa"="MSXa[time250]", 
    "Time 250min MSxp"="MSXp[time250]", 
    "Time 250min MSxa + MSxp"="MSXa+MSXp[time250]",
    "Time 250min CD BWM"="CD_BWM[time250]")
elist <- elist[keep_elist]
elist_cvis <- lapply(1:length(elist), function(i) {
    x <- elist[[i]]
    cvis <- new("cvis")
    cvis@name <- names(elist)[i]
    cvis@idx <- x@idx
    cvis@proj <- x@proj
    # Also reduce PCA projection dimension to reduce size
    if(!is.null(x@proj[["PCA"]]))    cvis@proj[["PCA"]] <- x@proj[["PCA"]][,1:max_pc]
    local_pmeta <- data.frame(Cluster = x@cluster)
    rownames(local_pmeta) <- colnames(all_cds)[cvis@idx]
    #print(identical(rownames(local_pmeta), rownames(x@proj[[1]])))
    cvis@pmeta <- local_pmeta
    #cvis@notes <- "some notes"
    return(cvis)
})
names(elist_cvis) <- names(keep_elist)
elist <- elist_cvis




sexpr <- exprs(all_cds)
any(duplicated(fData(all_cds)$gene_short_name))
# If duplicate exists, make.names
rownames(sexpr) <- fData(all_cds)$gene_short_name
all_cds@assayData$exprs <- sexpr
rownames(fData(all_cds)) <- fData(all_cds)$gene_short_name
sexpr_nmlog <- all_cds@auxOrderingData$normalize_expr_data
rownames(sexpr_nmlog) <- fData(all_cds)$gene_short_name
all_cds@assayData$normalize_expr_data <- sexpr_nmlog
all_cds@auxOrderingData$normalize_expr_data <- NULL


# convert all non-numeric meta into factor, and replace _ with space to make names look nicer
for (x in as.character(pmeta_attr$meta_id)) {
    if(grepl("time.bin",x)) next
    if(x %in% colnames(pData(all_cds))) {
        if(!is.numeric(pData(all_cds)[[x]])) {
            print(x)
            unique_levels<- as.character(unique(pData(all_cds)[[x]]))
            if("unannotated" %in% unique_levels) unique_levels <- c(unique_levels[unique_levels!="unannotated"], "unannotated")
            replace_levels <- gsub("_", " ", unique_levels, fixed=TRUE)
            names(replace_levels) <- unique_levels
            pData(all_cds)[[x]] <- replace_levels[as.character(pData(all_cds)[[x]])]
            pData(all_cds)[[x]] <- factor(pData(all_cds)[[x]], levels = replace_levels)
        }
    } 
}

pData(all_cds)$embryo.time.bin[which(pData(all_cds)$embryo.time.bin == "unannotated")] <- NA
pData(all_cds)$raw.embryo.time.bin[which(pData(all_cds)$embryo.time.bin == "unannotated")] <- NA
first_num <- as.numeric(stringi::stri_extract_first_regex(unique(pData(all_cds)$embryo.time.bin), "[0-9]+"))
bin_level <- unique(pData(all_cds)$embryo.time.bin)[order(first_num)]
pData(all_cds)$embryo.time.bin <- factor(pData(all_cds)$embryo.time.bin, levels = bin_level[!is.na(bin_level)])
pData(all_cds)$raw.embryo.time.bin <- factor(pData(all_cds)$raw.embryo.time.bin, levels = bin_level[!is.na(bin_level)])
# Now construct a "meta attribute dataframe" to store the info of what each column is about, and which ones you want to show
# This will be used in the 'Color By' option




# Load new t250 lineage annotaion
meta_new <- read.csv("data-raw/t250_metadata_revised.csv")
meta_old <- pData(all_cds)
meta_old$t250.lineages <- as.character(meta_old$t250.lineages)
meta_old$t250.lineages[match(meta_new$X, rownames(meta_old))] <- as.character(meta_new$t250.lineages.new)

unique_levels<- unique(as.character(meta_old$t250.lineages))
if("unannotated" %in% unique_levels) unique_levels <- c(unique_levels[unique_levels!="unannotated"], "unannotated")
meta_old$t250.lineages <- factor(meta_old$t250.lineages, levels = unique_levels)
pData(all_cds) <- meta_old


# Make a new cds to clean up stuff?

pd <- new("AnnotatedDataFrame", data = pData(all_cds))
fd <- new("AnnotatedDataFrame", data = fData(all_cds))
cds <- newCellDataSet(exprs(all_cds), phenoData = pd, featureData = fd)
cds@assayData$normalize_expr_data <- all_cds@assayData$normalize_expr_data
all_cds <- cds


cello <- new("ExpressionSet", assayData = assayDataNew("environment", exprs=exprs(all_cds)),
             phenoData =  pd,
             featureData = fd)
cello@assayData$normalize_expr_data <- all_cds@assayData$normalize_expr_data



# Graph plot

graph_genes <- graph_genes[order(graph_genes)]
g_meta_list<- readRDS("data-raw/image_gene_meta_list.rds")
g_meta_list <- lapply(g_meta_list, function(x) {
    colnames(x) <- c("Gene", "Series", "Allele", "Strain", "Reporter type", "Data source")
    x$Gene <- NULL
    x[["Data source"]][which(x[["Data source"]] == "EPIC")] <- "EPiC"
    x <- x[,c("Series", "Reporter type", "Strain", "Allele", "Data source")]
    return(x)
})




# Load Sup Table S1 and modify
ct_marker_tbl <- read.csv("data-raw/Table S1_ marker genes for terminal cell types - Sheet1.csv")
ct_marker_tbl$X <- NULL
cell_type_markers <- ct_marker_tbl
cur_umap_names<-as.character(unique(cell_type_markers$UMAP))
names(cur_umap_names) <- cur_umap_names
cur_umap_names<-tools::toTitleCase(cur_umap_names)
cur_umap_names[which(!cur_umap_names %in% names(clist))]
cur_umap_names["Hypodermis and seam cells"] <- "Hypodermis and Seam"
cur_umap_names["Early embryo, germline, and rectum"] <- "Early Embryo, Germline and Rectum"
cur_umap_names[which(!cur_umap_names %in% names(clist))]

cell_type_markers$UMAP <- cur_umap_names[as.character(cell_type_markers$UMAP)]


# Load Sup Table S5 and modify
ct_marker_s5 <- read.csv("data-raw/Table S5_ marker genes for pre-terminal lineages - t250.csv")
lineage_markers <- ct_marker_s5
replace_names <- names(keep_elist)
keep_elist2 <- keep_elist
keep_elist2[which(keep_elist2 == "MSXp[time250]")] <- "MSxp[time250]"
names(replace_names) <- keep_elist2
lineage_markers$UMAP <- as.character(lineage_markers$UMAP)
lineage_markers$UMAP <- replace_names[as.character(lineage_markers$UMAP)]
lineage_markers<-lineage_markers[lineage_markers$Lineage.Name != "Unknown MSxa descendants", ]
write.csv(lineage_markers, "data-raw/lineage_markers_cleaned.csv")

x <- as.character(lineage_markers$Markers)
genes<-trimws(unlist(strsplit(x, ",")), which = "both")
genes[which(!genes %in% gene_tbl$`Gene names:`)]


# update lin250

# Rename ABalaapppp to ABalapxpap (adding these cells to the existing ABalapxpap cells)
# Rename ABalapaapp to ABalaapppp/ABalapaapp
# Rename ABalaapppa/ABalapaapa to ABalapxpaa
# Rename ABalapxpaa to ABalapxppp
# Rename ABarpapp/ABplaaap to ABarpappp/ABplaaapp
# Rename ABarpapaa/ABplaaaaa to ABarpappa/ABplaaapa

name_change <- c(
    "ABalapxpap" = "ABalapxpap",
    "ABalapaapp" = "ABalaapppp/ABalapaapp",
    "ABalaapppa/ABalapaapa" = "ABalapxpaa",
    "ABalapxpaa" = "ABalapxppp",
    "ABarpapp/ABplaaap" = "ABarpappp/ABplaaapp",
    "ABarpapaa/ABplaaaaa" = "ABarpappa/ABplaaapa"
)

pData(all_cds)$t250.lineages <- as.character(pData(all_cds)$t250.lineages)
for(i in 1:length(name_change)) {
    old_name <- names(name_change)[i]
    new_name <- name_change[i]
    pData(all_cds)$t250.lineages[which(pData(all_cds)$t250.lineages == old_name)] <- new_name
}

# Rename the cells in this attached metadata file to ABalaapppa/ABalapaapa 
Actual.RIA.NB.ABalaapppa.and.ABalapaapa <- read.csv("~/Documents/LYNCH/Celegans/VisCello/data-raw/Actual RIA NB ABalaapppa and ABalapaapa.txt", row.names=1, stringsAsFactors=FALSE)
pData(all_cds)$t250.lineages[which(rownames(pData(all_cds)) %in% rownames(Actual.RIA.NB.ABalaapppa.and.ABalapaapa))] <- "ABalaapppa/ABalapaapa"

# Rename the cells in this attached metadata file to ABarpapaa/ABplaaaaa
Actual.ABplaaaaa.lineage <- read.csv("~/Documents/LYNCH/Celegans/VisCello/data-raw/Actual ABplaaaaa lineage.txt", row.names=1, stringsAsFactors=FALSE)
pData(all_cds)$t250.lineages[which(rownames(pData(all_cds)) %in% rownames(Actual.ABplaaaaa.lineage))] <- "ABarpapaa/ABplaaaaa"


# Convert all_cds monocle to ExpressionSet to reduce memory use
fmeta <- fData(all_cds)
fmeta <- fmeta[, c(1,2)]
colnames(fmeta) <- c("id", "symbol")
eset <- new("ExpressionSet",
            assayData = assayDataNew( "environment", exprs=exprs(all_cds), norm_exprs = all_cds@assayData$normalize_expr_data),
            phenoData =  new("AnnotatedDataFrame", data = pData(all_cds)),
            featureData = new("AnnotatedDataFrame", data = fmeta))
eset$Size_Factor <- all_cds$Size_Factor[match(colnames(eset), colnames(all_cds))]
saveRDS(eset, paste0("inst/app/data/eset.rds"))
saveRDS(eset, paste0("data-raw/eset.rds"))


# Update annotation to latest
eset <- readRDS(paste0("data-raw/eset.rds"))
cds_190301 <- readRDS("data-raw/for-qin.cds.all.bg.corrected.unfiltered.rds")
identical(colnames(eset), colnames(cds_190301)) # TRUE - good to go

# Remove annotations not to be shown
remove_ids <- c("mm.lineage", "temp.ABala.250")
pData(eset)<-pData(eset)[, -which(colnames(pData(eset)) %in% remove_ids)]

pData(eset)$combine_lineage <- pData(cds_190301)$lineage
pData(eset)$combine_lineage[which(is.na(pData(eset)$combine_lineage))] <- "unannotated"
saveRDS(eset, paste0("inst/app/data/eset.rds"))

pmeta_attr <- rbind(pmeta_attr, c(meta_id = "combine_lineage", meta_name = "Lineage", is_numeric=F, dpal = "Set1", dscale = NA))
save(pmeta_attr, file = paste0("inst/app/data/pmeta_attr.rda"))

# Additional change from master eset
Actual.ABalapxpaa <- read.csv("~/Documents/LYNCH/Celegans/VisCello/data-raw/Actual ABalapxpaa.txt")
pData(eset)$t250.lineages[which(rownames(pData(eset)) %in% Actual.ABalapxpaa$X)] <- "ABalapxpaa"

pData(eset)$t250.lineages[which(pData(eset)$t250.lineages == "ABalaapppp")] <- "IL1/IL2/OLQ/OLL neuroblasts"
pData(eset)$t250.lineages[which(pData(eset)$t250.lineages == "Unknown MSxa descendants")] <- "unannotated"
saveRDS(eset, paste0("inst/app/data/eset.rds"))


# saveRDS(eset, paste0("inst/app/data/eset1.rds"))
# save(eset, file=paste0("inst/app/data/eset1.RData"))

# system.time(eset <- readRDS(paste0("inst/app/data/eset1.rds")))
# system.time(load(paste0("inst/app/data/eset1.RData")))









Cello <- setClass("Cello",
                  slots = c(
                      name = "character", # The name of cvis
                      idx = "numeric", # The index of global cds object
                      proj = "list", # The projections as a list of data frame
                      pmeta = "data.frame", # The local meta data
                      notes = "character" # Other information to display to the user
                  )
)

clist_cello <- lapply(1:length(clist), function(i) {
    x <- clist[[i]]
    cvis <- new("Cello")
    cvis@name <- names(clist)[i]
    cvis@idx <- x@idx
    cvis@proj <- x@proj
    # Also reduce PCA projection dimension to reduce size
    if(!is.null(x@proj[["PCA"]]))    cvis@proj[["PCA"]] <- x@proj[["PCA"]]
    #print(identical(rownames(local_pmeta), rownames(x@proj[[1]])))
    cvis@pmeta <- x@pmeta
    #cvis@notes <- "some notes"
    return(cvis)
})
names(clist_cello) <- names(clist)
clist <- clist_cello
base::save(clist, file="data/clist.rda")

elist_cello <- lapply(1:length(elist), function(i) {
    x <- elist[[i]]
    cvis <- new("Cello")
    cvis@name <- names(elist)[i]
    cvis@idx <- x@idx
    cvis@proj <- x@proj
    # Also reduce PCA projection dimension to reduce size
    if(!is.null(x@proj[["PCA"]]))    cvis@proj[["PCA"]] <- x@proj[["PCA"]]
    #print(identical(rownames(local_pmeta), rownames(x@proj[[1]])))
    cvis@pmeta <- x@pmeta
    #cvis@notes <- "some notes"
    return(cvis)
})
names(elist_cello) <- names(elist)
elist <- elist_cello
base::save(elist, file="data/elist.rda")


### Final test and load data ### 
# usethis::use_data(clist, overwrite = T)
# usethis::use_data(elist, overwrite = T)
# usethis::use_data(all_cds, overwrite = T)
# usethis::use_data(pmeta_attr, overwrite = T)
# 
# usethis::use_data(g_all, overwrite = T)
# usethis::use_data(g_meta_list, overwrite = T)
# 
# usethis::use_data(tf_tbl, overwrite = T)
# usethis::use_data(cell_type_markers, overwrite = T)
# usethis::use_data(lineage_markers, overwrite = T)
# usethis::use_data(graph_genes, overwrite = T)
# 
# tools::add_datalist("../VisCello")

devtools::load_all()




# S3 and S6

s6_tbl<-read.table(file = 'data-raw/Table_S6.tsv', sep = '\t', header = TRUE)
saveRDS(s6_tbl, "data-raw/s6_tbl.rds")
saveRDS(s6_tbl, "inst/app/data/s6_tbl.rds")


s7_tbl<-read.table(file = 'data-raw/Table_S7.tsv', sep = '\t', header = TRUE)
saveRDS(s7_tbl, "data-raw/s7_tbl.rds")
saveRDS(s7_tbl, "inst/app/data/s7_tbl.rds")



# S6 map to graph


image_lin <- as.data.frame(g_all)$name
orig_lin <- as.character(levels(s6_tbl$lineage))

# Match name and modify if necessary
sexpr_lin <- strsplit(orig_lin, "/")
lin_split <- lapply(sexpr_lin, function(x) {strsplit(x, "x")})
sexpr_lin[which(!sexpr_lin %in% image_lin)]

# Change x to proper 2 lin

first_pass<-sapply(lin_split, function(x){
    unlist(lapply(x, function(l) {
        if(length(l) ==2) {
            if(l[1] == "ABp") {
                return(c(paste(l, collapse = "l"), paste(l, collapse = "r")))
            } else if(l[1] == "MS") {
                return(c(paste(l, collapse = "a"), paste(l, collapse = "p")))
            } else if(l[1] == "C") {
                return(c(paste(l, collapse = "a"), paste(l, collapse = "p")))
            } else if(l[1] == "D") {
                return(c(paste(l, collapse = "a"), paste(l, collapse = "p")))
            } else if(l[1] %in% c("ABalap", "ABarpp","ABaraap")){
                return(c(paste(l, collapse = "a"), paste(l, collapse = "p")))
            } 
        } else if(length(l) == 1) {
            return(l)
        }
        return(NA)
    }))
})

names(first_pass) <- orig_lin
which(sapply(first_pass, function(x) {any(is.na(x))}))


# Manual fix
first_pass["ABpraaa_and_ABalppp"] <- list(c("ABpraaa", "ABalppp"))
#first_pass["MSxppaa_or_MSxppaax_before_branching"] <- list(c("MSappaa", "MSpppaa", "MSappaaa", "MSappaap", "MSpppaaa", "MSpppaap"))
first_pass["MSxppaa_or_MSxppaax_before_branching"] <- NA # Seems some are duplicated
first_pass["MSxpppxa/MSxpappx"] <- list(c("MSapppaa", "MSappppa", "MSppppaa","MSpppppa", "MSapappa","MSapappp","MSppappa", "MSppappp"))
first_pass["MSxpppxp"] <- list(c("MSapppap","MSappppp","MSppppap","MSpppppp"))
first_pass["MSxxx"] <- list(c("MSaaa","MSaap","MSapa","MSapp","MSpaa","MSpap","MSppa","MSppp"))

unlist(first_pass)[which(duplicated(unlist(first_pass)))]
# MAP expression to graph, for duplicated cell, use the one with shorter name (more specific in annotation)
graph_tbl <- as.data.frame(g_all)
elin_idx<-sapply(graph_tbl$name, function(x) {
    idx<-which(sapply(first_pass, function(y){
        x %in% y
    }))
    if(length(idx) > 1) {
        idx <- idx[which.min(sapply(first_pass[idx], length))]
    }
    return(idx)
})

elin_match<-sapply(elin_idx, function(i)orig_lin[i])

usethis::use_data(elin_match, overwrite = T)

e_tpm <- sapply(e_idx, function(i){
    if(length(i)){
        s6_tbl$bootstrap.median.tpm[i]
    } else {
        NA
    }
})

g_all <- g_all %>% activate("nodes") %>% 
    mutate(scExprTPM = e_tpm)

plot_col <- "scExprTPM"
t_cut <- 80
g<-g_all %>% activate("nodes") %>% 
    mutate(text.size = ifelse(time > t_cut, 0, 1/log(time))) %>%
    mutate(name = ifelse(time > t_cut, "", name)) %>%
    filter(!(time > 200 &  is.na(!!as.name(plot_col))))
range(as.data.frame(g)$text.size)
g1<-plotGraph(g, color.by=plot_col, pal="gg_color_hue", label="name", type = "numeric",border.size=.3, background="white") + 
    theme(text=element_text(family = "Helvetica", size=5),
          axis.ticks.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(-10,-10,-10,-10),
          plot.margin = unit(c(.3,.5,.3,.3), "cm"))





## Post submission update: lncRNA quantification
eset <- readRDS("inst/app/data/eset.rds")
source("data-raw/scripts/Load_10x_data.R")

sample_dirs <- c(
    "Waterston 300 minutes"= "/kimdata/zhuqin/celegans/lncQuant/Waterston_300_min_ncRNA/",
    "Waterston 400 minutes"="/kimdata/zhuqin/celegans/lncQuant/Waterston_400_min_ncRNA/",
    "Waterston 500 minutes batch 1" = "/kimdata/zhuqin/celegans/lncQuant/Waterston_500_min_batch_1_ncRNA/",
    "Waterston 500 minutes batch 2" = "/kimdata/zhuqin/celegans/lncQuant/Waterston_500_min_batch_2_ncRNA/",
    "Murray r17"="/kimdata/zhuqin/celegans/lncQuant/Murray_r17_ncRNA/",
    "Murray b01" = "/kimdata/zhuqin/celegans/lncQuant/Murray_b01_ncRNA/",
    "Murray b02"="/kimdata/zhuqin/celegans/lncQuant/Murray_b02_ncRNA/"
)

pmeta <- pData(eset)
pmeta$barcode <- paste0(sapply(strsplit(rownames(pData(eset)), "-"), function(x)x[[1]]), "-1")


Counter=0
for(i in 1:length(sample_dirs)) {
    dirName <- sample_dirs[i]
    batch <- names(sample_dirs)[i]
    message(dirName)
    raw_matrix <- load_cellranger_matrix(dirName, barcode_filtered = F)
    barcode_in <- pmeta$barcode[which(pmeta$batch == batch)]
    names(barcode_in) <- rownames(pmeta)[which(pmeta$batch == batch)]
    ThisMatrix = exprs(raw_matrix)[, match(barcode_in, colnames(raw_matrix))]
    colnames(ThisMatrix)<-names(barcode_in)
    print(dim(ThisMatrix))
    if(Counter==0){
        Counter=Counter+1
        CombinedMatrix = ThisMatrix
    }else{
        Counter=Counter+1
        CombinedMatrix = cbind(CombinedMatrix,ThisMatrix)
    }        
}

CombinedMatrix <- CombinedMatrix[, match(colnames(eset), colnames(CombinedMatrix))]
saveRDS(CombinedMatrix, "data-raw/rds/lncRNA_combined_mtx.rds")


# Also load annotation to add gene name

library(rtracklayer)
lncs <- import("../ce.WS260.ncRNA/genes/genes.gtf")

tid <- sapply(rownames(CombinedMatrix), function(x) {
    print(x)
    cur_lnc <- lncs[which(lncs$gene_id == x)]
    cur_tid <- cur_lnc$transcript_id[which(!is.na(cur_lnc$transcript_id))]
    return(cur_tid[1])
})
fmeta <- data.frame(id = rownames(CombinedMatrix), symbol = tid)
sum(duplicated(fmeta$symbol))
rownames(fmeta) <- fmeta$symbol
rownames(CombinedMatrix) <- fmeta$symbol

# Normalization together with mRNA

all_cds <- readRDS(paste0("/kimdata/zhuqin/celegans/data_analysis/CAtlas/data/all_cds.rds"))
identical(colnames(all_cds), colnames(CombinedMatrix))

lnc_cds<- newCellDataSet(CombinedMatrix,
            phenoData = new("AnnotatedDataFrame", data = pData(eset)),
            featureData = new("AnnotatedDataFrame", data = fmeta))
identical(rownames(pData(lnc_cds)), rownames(pData(eset)))
pData(lnc_cds)$Size_Factor <- pData(eset)$Size_Factor
source("data-raw/scripts/compute_dimR.R")
FM <- normalize_expr_data2(lnc_cds, "log", 1, use_order_gene = F)
saveRDS(FM, paste0("data-raw/rds/lnc_normalized_FM.rds"))

# Archive eset before change
# saveRDS(eset, paste0("data-raw/rds/eset_archive0403.rds"))


combine_expr <- Matrix(rbind(exprs(eset), exprs(lnc_cds)), sparse = T)
combind_norm <- Matrix(rbind(eset@assayData$norm_exprs, FM), sparse = T)
combined_fmeta <- rbind(fData(eset), fData(lnc_cds))

eset <- new("ExpressionSet",
            assayData = assayDataNew( "environment", exprs=combine_expr, norm_exprs = combind_norm),
            phenoData =  new("AnnotatedDataFrame", data = pData(eset)),
            featureData = new("AnnotatedDataFrame", data = combined_fmeta))

saveRDS(eset, paste0("inst/app/data/eset.rds"))

expr_cnt <- rowSums(exprs(eset))





# Post submission update: t150 lineages
eset <- readRDS("inst/app/data/eset.rds")
load("data-raw/postsub_update/State-2019-03-02_with_t150_updates.rda")
new_t150_df <- r_data$cmeta$df
nrow(new_t150_df)

pData(eset) <- readRDS("data-raw/pData_archive/pData-190317.rds")
identical(rownames(new_t150_df), rownames(pData(eset)))

#saveRDS(pData(eset), "data-raw/pData_archive/pData-190317.rds")
pData(eset)$t150.lineages <- new_t150_df$t150.lineages

# Wipe annotation before t250
early_cidx <- which(pData(eset)$raw.embryo.time <= 250)
pData(eset)$combine_lineage[early_cidx] <- "unannotated"

# Write over combined lineages
pData(eset)$combine_lineage[which(pData(eset)$t150.lineages!="unannotated")] <- pData(eset)$t150.lineages[which(pData(eset)$t150.lineages!="unannotated")]

#saveRDS(eset, "inst/app/data/eset.rds")

# Post submission update: NewLineageAnnots_post_t150_updates4.rda
load("data-raw/postsub_update/NewLineageAnnots_post_t150_updates4.rda")
new_lin_df <- r_data$cmeta$df
new_lin_df$NewLineage[which(new_lin_df$NewLineage == "ALM_PLM cluster")] <- "unannotated"
nrow(new_lin_df)

identical(rownames(new_lin_df), colnames(eset))
update_idx <- which(new_lin_df$NewLineage != "unannotated")
# Don't update t150 annotation
ignore_idx <- which(pData(eset)$t150.lineages != "unannotated")
length(ignore_idx)
update_idx <- setdiff(update_idx, ignore_idx)
length(update_idx)

pData(eset)$t150.lineages <- NULL
pData(eset)$t250.lineages <- NULL
pData(eset)$combine_lineage[update_idx] <- new_lin_df$NewLineage[update_idx]

# Check if any of the old lineages that should be removed are still present
missing_lin_markers <- read.csv("data-raw/postsub_update/missing_lin_markers.csv", row.names = 1)
pData(eset)$combine_lineage[which(pData(eset)$combine_lineage %in% missing_lin_markers$Lineage.Name[missing_lin_markers$Operation == "RM"])] <- "unannotated"

# Reverse order of r and l
lin_split <- strsplit(pData(eset)$combine_lineage, "/")
lin_split_len<-sapply(lin_split, length)
#lin_split[which(lin_split_len > 2)]
lin_sorted <- lapply(lin_split, sort)
lin_sorted <- sapply(lin_sorted, function(x) {
    paste0(x, collapse = "/")
})

correct_idx <- which(lin_sorted != pData(eset)$combine_lineage)
lin_correct <- data.frame(Orig_name = pData(eset)$combine_lineage[correct_idx], New_name = lin_sorted[correct_idx])
lin_correct <- lin_correct[!duplicated(lin_correct), ]
lin_correct <- lin_correct[order(lin_correct$New_name),]
write.csv(lin_correct, paste0("data-raw/postsub_update/lineage_name_resorted.csv"))

pData(eset)$combine_lineage <- lin_sorted
sum(pData(eset)$combine_lineage != "unannotated")
saveRDS(pData(eset), "data-raw/pData_archive/pData-190417_updated.rds")
saveRDS(eset, "inst/app/data/eset.rds")




# Update lineage marker sheet 

load("data-raw/postsub_update/NewLineageAnnots_post_t150_updates4.rda")
# Previous lineage
lineage_markers <- read.xlsx("data-raw/postsub_update/LineageMarkers_addPrevious.xlsx")
lineage_markers$Alternate.name <- NULL
lineage_markers$UMAP.group <- NULL
# Count cells for the annoatated lineage
lineage_markers$Lineage.Name <- sapply(lapply(strsplit(lineage_markers$Lineage.Name, "/"),sort), function(x)paste0(x, collapse = "/"))

lineage_markers$count <- table(pData(eset)$combine_lineage)[as.character(lineage_markers$Lineage.Name)]
lineage_markers[which(is.na(lineage_markers$count) & !lineage_markers$Cells.produced %in% c("x", "death")),]

unique(lineage_markers$UMAP)[!unique(lineage_markers$UMAP) %in% c(names(r_data$usr$clist), names(r_data$usr$elist))]

# saveRDS(elist, "data-raw/pData_archive/elist_190418_archive.rds")
# saveRDS(clist, "data-raw/pData_archive/clist_190418_archive.rds")
cur_elist <- readRDS("data-raw/pData_archive/elist_190418_archive.rds")
cur_clist <- readRDS("data-raw/pData_archive/clist_190418_archive.rds")
new_elist_name <- unique(lineage_markers$UMAP)[unique(lineage_markers$UMAP) %in% c(names(r_data$usr$clist), names(r_data$usr$elist))]
new_elist_name <- setdiff(new_elist_name, c(names(cur_elist), names(cur_clist)))
elist <- c(cur_elist, c(r_data$usr$clist, r_data$usr$elist)[new_elist_name])
save(elist, file="inst/app/data/elist.rda")


# Update cells produced
library(stringi)
library(tibble)
cell_list <- read.csv("data-raw/celegans_cell_list.csv")
cell_list <- cell_list %>% mutate(lineage_name =stri_replace_all_fixed(cell_list$Lineage.Name, " ", ""))

Lineage.Name.expanded <- strsplit(lineage_markers$Lineage.Name, "/")

lin_expand <-sapply(Lineage.Name.expanded, function(x){
    paste0(sort(unlist(lapply(x, function(l) {
        dts<- grep(paste0("^", stri_replace_all_fixed(l, "x", ".")), cell_list$lineage_name, value = T)
        unique(strtrim(dts, nchar(l)))
    }))), collapse = "/")
})
lineage_markers<- add_column(lineage_markers,  Lineage.Name.expanded= lin_expand, .after = 1)
lineage_markers$Lineage.Name.expanded[lineage_markers$Lineage.Name.expanded == ""] <- lineage_markers$Lineage.Name[lineage_markers$Lineage.Name.expanded == ""]
cell_list$Cell <- as.character(cell_list$Cell)
dt_cells <- sapply(Lineage.Name.expanded, function(x){
    paste0(unlist(lapply(x, function(l) {
        dts<- sort(grep(paste0("^", stri_replace_all_fixed(l, "x", ".")), cell_list$lineage_name, value = T))
        dts <- unique(cell_list$Cell[match(dts, cell_list$lineage_name)])
    })), collapse = "/")
})

lineage_markers$Cells.produced <- dt_cells

# Also combine columns note and trajectory leads to

lineage_markers$Notes <- ifelse(!is.na(lineage_markers[["Trajectory.leads.to:"]]), 
                                    ifelse(is.na(lineage_markers$Notes), 
                                           paste0("trajectory leads to: ", lineage_markers[["Trajectory.leads.to:"]]), 
                                           paste0(lineage_markers$Notes, "; trajectory leads to: ", lineage_markers[["Trajectory.leads.to:"]])),
                                    lineage_markers$Notes)

lineage_markers[["Trajectory.leads.to:"]] <- NULL
lineage_markers$scRNA.cell.count <- lineage_markers$count
lineage_markers$count <- NULL
# Reorder columns
col_order <- c("Lineage.Name", "Lineage.Name.expanded", "Depth", "Cells.produced", "UMAP","Markers", "New.markers", "scRNA.cell.count","Notes") 
lineage_markers <- lineage_markers[,col_order]
# Reorder by lineage
lineage_markers <- lineage_markers[order(lineage_markers$Lineage.Name), ]

lineage_markers$UMAP[which(!lineage_markers$UMAP %in% c(names(clist), names(elist)))] <- ""
write.xlsx(lineage_markers, "data-raw/postsub_update/LineageMarkers_JM_QZcleaned.xlsx")
save(lineage_markers, file = "inst/app/data/lineage_markers.rda")


load("data-raw/postsub_update/NewLineageAnnots_post_t150_updates4.rda")

replace_names <- names(keep_elist)
keep_elist2 <- keep_elist
keep_elist2[which(keep_elist2 == "MSXp[time250]")] <- "MSxp[time250]"
names(replace_names) <- keep_elist2
lineage_markers$UMAP <- as.character(lineage_markers$UMAP)
lineage_markers$UMAP <- replace_names[as.character(lineage_markers$UMAP)]
lineage_markers<-lineage_markers[lineage_markers$Lineage.Name != "Unknown MSxa descendants", ]
write.csv(lineage_markers, "data-raw/lineage_markers_cleaned.csv")

x <- as.character(lineage_markers$Markers)
genes<-trimws(unlist(strsplit(x, ",")), which = "both")
genes[which(!genes %in% gene_tbl$`Gene names:`)]




