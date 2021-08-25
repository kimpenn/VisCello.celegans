


#' @export
get_tree_color <- function(ctree, pal) {
    tcolor <- colorRampPalette(pal)(length(as_tibble(ctree)$label)) ## (n)
    names(tcolor) <- sort(as_tibble(ctree)$label)
    return(tcolor)
}

#' @export
make_lineage_ggtree <- function(in_tree = NULL, root = "P0", time.cut = 300, 
                                color.annot = "label", branch.length='lifetime', 
                                tip.lab = F, tip.lab.size = 2, tip.lab.align = F,tip.lab.angle = 0,
                                node.lab = F, node.lab.size = 4, node.lab.timecut = 80,node.lab.angle = 90,
                                xmax = 450, xlim = NULL,
                                color.pal = NULL, edge.size = 1) {
    in_tree<-in_tree %>% filter(br_time <=time.cut)
    plot_tree <- in_tree %>% filter(to %in% find_child_root(root, in_tree))
    correct_idx <- which(plot_tree$d_time > xmax)
    plot_tree$d_time[correct_idx] = xmax
    plot_tree$lifetime[correct_idx] <- (plot_tree$d_time - plot_tree$br_time)[correct_idx]
    plot_tree <- as.treedata(plot_tree)
    out_tree <- ggtree(plot_tree, branch.length=branch.length, aes_string(color = color.annot), size = edge.size, ladderize = F)
    out_tree <- lineage_tree_flip(out_tree, silence = T)
    
    if(!is.null(color.pal)) {
        out_tree<-out_tree + scale_color_manual(values = color.pal, na.value= "grey")
    }
    if(tip.lab) {
        out_tree<-out_tree + 
            geom_tiplab(size = tip.lab.size,color = "black", align = tip.lab.align, angle = tip.lab.angle)
    }
    if(node.lab) {
        out_tree<-out_tree + geom_text(aes(x=branch, label=label), data = out_tree$data[out_tree$data$br_time < node.lab.timecut,], size = node.lab.size, angle=node.lab.angle, color = "black")
    }
    if(!is.null(xlim)) {
        out_tree<-out_tree + 
            scale_x_continuous(lim = xlim)
    }
    return(out_tree)
}


#' @export
make_tree_dimr <- function(proj = NULL, left_tree = NULL, right_tree = NULL, top_tree = NULL, 
                        colorBy = "combine_lineage", shared.col = "combine_lineage",
                        label.time.cut = 80, label.type = "label",
                           pt.size = 1, edge.size = 1, label.size = 3,
                           tip.size = 0,
                           tree.color = NULL, tree.h.scale = 1/4, shift.x.scale = 1/7, shift.y.scale = 1/10,
                           plot.link = NULL, link.size = .5, link.alpha = .5, link.x.shift = 30, link.x.shift.absolute = F,
                           return_coords = F) {
    if(is.null(proj) ) {
        stop("projection must be specified")
    }
    
    colnames(proj)[c(1,2)] <- c("x","y")
    proj_x_range <- range(proj$x)
    proj_y_range <- range(proj$y)
    label_function <- if(label.type == "label") geom_label else if(label.type == "text") geom_text else NULL
    
    lr_x_length <- (proj_x_range[2] - proj_x_range[1]) * tree.h.scale
    tp_y_length <- (proj_y_range[2] - proj_y_range[1]) * tree.h.scale
    lr_x_adjust <- (proj_x_range[2] - proj_x_range[1]) * shift.x.scale
    lr_y_adjust <- (proj_y_range[2] - proj_y_range[1]) * shift.y.scale
    left_x_range <- c(proj_x_range[1] - lr_x_length -  lr_x_adjust, proj_x_range[1] -  lr_x_adjust)
    right_x_range <-  c(proj_x_range[2] + lr_x_adjust, proj_x_range[2] +  lr_x_adjust  +  lr_x_length)
    top_y_range <- c(proj_y_range[2] + lr_y_adjust, proj_y_range[2] + lr_y_adjust+ tp_y_length)
    
    pp <- left_tree 
    pp <- plotProj(proj, group.by = colorBy, pal = tree.color, legend = F) + theme_void()
    
    if(!is.null(left_tree)) {
        d1 <- left_tree$data
        d1$branch <-  scales::rescale(d1$branch, to = left_x_range)
        d1$x <- scales::rescale(d1$x, to = left_x_range)
        d1$y <- scales::rescale(d1$y, to = proj_y_range)
        pp <- pp+
            geom_tree(data=d1, size = edge.size)
        if(!is.null(label_function)) {
            pp <- pp+label_function(aes(x=branch, y = y, label=label), data =d1[d1$br_time < label.time.cut | d1$branch.length == 0,], size = label.size, angle = 90,color = "black")}
        if(tip.size) {
            pp <- pp + geom_tiplab(data=d1, size = tip.size, color = "black", align = F, hjust = -.1) 
        }
    }
    if(!is.null(right_tree)) {
        d2 <- right_tree$data
        d2$x <- scales::rescale(-d2$x, to = right_x_range)
        d2$branch <-  scales::rescale(-d2$branch, to = right_x_range)
        d2$y <- scales::rescale(d2$y, to = rev(proj_y_range))
        pp <- pp+
            geom_tree(data=d2, size = edge.size)
        if(!is.null(label_function)) {
            pp <- pp+label_function(aes(x=branch, y = y, label=label), data =d2[d2$br_time < label.time.cut | d2$branch.length == 0,], size = label.size, angle = -90,color = "black")}
        if(tip.size) {
            pp <- pp + geom_tiplab(data=d2, size = tip.size, color = "black", align = F, hjust = 1.1) 
        }
    }
    if(!is.null(top_tree)) { 
        d3 <- top_tree$data
        d3$y <-  -top_tree$data$x
        d3$x <- top_tree$data$y
        d3$y <- scales::rescale(d3$y, to = top_y_range)
        d3$x <- scales::rescale(d3$x, to = proj_x_range)
        pp <- pp + geom_segment(data=d3, aes(x    = x,
                                             xend = x,
                                             y    = y[parent],
                                             yend = y), size = edge.size )+
            geom_segment(data=d3, aes(x    = x[parent],
                                      xend = x,
                                      y    = y[parent],
                                      yend = y[parent]), size = edge.size )
        if(!is.null(label_function)) {
            pp <- pp+label_function(aes(x=x, y=y+tp_y_length/20, label=label), data =d3[d3$br_time < label.time.cut | d3$branch.length == 0,], size = label.size, angle = 0,color = "black")}
        if(tip.size) {
            pp <- pp + geom_tiplab(data=d3, size = tip.size, color = "black", align = F, angle = -90, hjust = -.1) 
        }
    }
    pp <- pp
    if(!is.null(tree.color)) {
        pp <- pp+ 
            scale_color_manual(values = tree.color, na.value = "grey")
    }
    #assign("pp", pp, env=.GlobalEnv)
    if(!is.null(plot.link)){
        proj <- proj[, c("x","y", shared.col)]
        proj_center <- proj %>% group_by_at(shared.col) %>% summarize_at(c("x", "y"), median)
        use_col <- c("x","y",shared.col)
        df_bind <- list(NULL,NULL,NULL)
        
        if(!is.null(left_tree)) {
            df_bind[[1]]<-d1[,use_col]
            df_bind[[1]]$source <- "left"
            if(link.x.shift.absolute) {
                df_bind[[1]]$x <- link.x.shift
            } else {
                df_bind[[1]]$x <- df_bind[[1]]$x + link.x.shift
            }
        }
        
        if(!is.null(right_tree)) {
            df_bind[[2]]<-d2[,use_col]
            df_bind[[2]]$source <- "right"
            df_bind[[2]]$x <- df_bind[[2]]$x - link.x.shift
        }
        if(!is.null(top_tree)) {
            df_bind[[3]]<-d3[,use_col]
            df_bind[[3]]$source <- "top"
        }
        dd_tree <- bind_rows(df_bind) 
        dd_tree$x2 <- proj_center$x[match(dd_tree[[shared.col]], proj_center[[shared.col]])]
        dd_tree$y2 <- proj_center$y[match(dd_tree[[shared.col]], proj_center[[shared.col]])]
        if(all(plot.link %in% c("left", "right", "top"))) {
            sub_tree<- dd_tree %>% filter(source %in% plot.link)
        } else {
            sub_tree <- dd_tree
        }
        pp<-pp + 
            geom_segment(aes(x = x, y=y, xend = x2, yend =y2), data=sub_tree, color='grey', alpha = link.alpha, size = link.size)
    }
    if(return_coords) {
        return(list(
            plot = pp,
            coords = list(
                proj = proj,
                d1 = if(!is.null(left_tree)) d1 else NULL,
                d2 = if(!is.null(right_tree)) d2 else NULL,
                d3 = if(!is.null(top_tree)) d3 else NULL
            )
        ))
    } else {
        return(pp)
    }
}















