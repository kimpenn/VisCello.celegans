


#' @export
substrRight <- function(x, n){
    substr(x, nchar(x)-n+1, nchar(x))
}

#' @export
# Computationally check if daugters of all parents have the right p-a order
lineage_tree_flip <- function(treeplot, silence = T) {
    p1 <- treeplot
    for(p in unique(p1$data$parent)) {
        cur_d <- which(p1$data$parent == p)
        cur_d_name <- p1$data$label[cur_d]
        cur_d_n <- p1$data$node[cur_d]
        if(length(cur_d) == 2) {
            cur_d_end <- substrRight(cur_d_name,1)
            cur_d_y <- p1$data$y[cur_d]
            names(cur_d_y) <- cur_d_end
            end_paste <- paste0(cur_d_end, collapse = ",")
            if(end_paste%in% c("a,p", "p,a", "d,v", "v,d", "l,r", "r,l")) {
                if(end_paste%in% c("a,p", "p,a")){
                    if(cur_d_y["p"] < cur_d_y["a"]) {
                        if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
                        p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
                    }
                }
                if(end_paste%in% c("d,v", "v,d")){
                    if(cur_d_y["v"] < cur_d_y["d"]) {
                        if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
                        p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
                    }
                }
                if(end_paste%in% c("l,r", "r,l")){
                    if(cur_d_y["r"] < cur_d_y["l"]) {
                        if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
                        p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
                    }
                }
            } else {
                cur_d_y <- p1$data$y[cur_d]
                names(cur_d_y) <- cur_d_name
                if(all(c("Z2", "Z3") %in% cur_d_name)) {
                    if(cur_d_y["Z2"] < cur_d_y["Z3"]) {
                        if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
                        p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
                    }
                }
                if(all(c("P2", "EMS") %in% cur_d_name)) {
                    if(cur_d_y["P2"] < cur_d_y["EMS"]) {
                        if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
                        p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
                    }
                }
                if(all(c("E", "MS") %in% cur_d_name)) {
                    if(cur_d_y["E"] < cur_d_y["MS"]) {
                        if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
                        p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
                    }
                }
                if(all(c("C","P3") %in% cur_d_name)) {
                    if(cur_d_y["P3"] < cur_d_y["C"]) {
                        if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
                        p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
                    }
                }
                if(all(c("D", "P4") %in% cur_d_name)) {
                    if(cur_d_y["P4"] < cur_d_y["D"]) {
                        if(!silence) message(paste0("Flipped: ", paste0(cur_d_name, collapse = ",")))
                        p1 <- p1 %>% flip(cur_d_n[1], cur_d_n[2])
                    }
                }
            }
        }
    }
    return(p1)
}




#' @export
find_child_root <- function(root, tree_tbl) {
    find_child <- function(root, tree_tbl) {
        dts <- tree_tbl$to[which(tree_tbl$from == root)]
        return(c(dts, lapply(dts, find_child, tree_tbl)))
    }
    
    res<-find_child(root, tree_tbl)
    res <- unlist(res, recursive = T)
    return(res)
}




#' @export
# Function to get lowest common ancester
get_ancestor <- function(df, node){
    anc <- get_parent(df, node)
    i <- 1
    while (i <= length(anc)) {
        anc <- c(anc, get_parent(df, anc[i]))
        i <- i + 1
    }
    return(anc)
}


#' @export
get_parent <- function (df, node) 
{
    parent_id <- df$from[df$to == node]
    parent_id[parent_id != node]
}

#' @export
get_level <- function (df, node, root = "P0") 
{
    lev = 0
    cur_n <- node
    parent  <- node
    while(parent != root) {
        parent <- get_parent(df, cur_n)
        cur_n <- parent
        lev <- lev+1
    }
    return(lev)
}

