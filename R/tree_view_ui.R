
#' @export
tree_ui <- function(id) {
    ns <- NS(id)
    gtbl <- gene_tbl
    colnames(gtbl) <- "sc-expression:"
    tagList(
        wellPanel(
            fluidRow(
                column(3, selectizeInput(ns("lin_group"), "Color by:", choices = c("Birth time"="br_time", gtbl), selected = "br_time")),
                column(3, selectizeInput(ns("tree_root"), "Tree root", choices = data.table(Node = avail_nodes[1:101]))),
                column(3, selectInput(ns("tree_style"), "Tree style:", choices = c("Rectangular" = "rect", "Slanted" = "slanted", "Arc bar" = "arc-bar"))),
                column(3, selectInput(ns("lin_pal"), "Palette", choices=image_palettes[3:6]))
            ),
            fluidRow(
                column(3, colourpicker::colourInput(ns("na_color"), "Umapped node color", value = "#d3d3d3")),
                column(3, numericInput(ns("max_time"), "Birth time cut:", value = 500, step = 50, min = 10, max = 1000)),
                column(3, numericInput(ns("plot_height"), "Plot height:", value = 5000, step = 100, min = 100, max = 6000)),
                column(3, selectInput(ns("tree_label_col"), "Tree Label", choices = c("Lineage+cell" = "lin_cell", "Lineage+cell+description" = "lin_cell_des", "Lineage" = "label", "Cell name" = "Cell", "Cell description" = "Description")))
            ),
            fluidRow(
                column(6),
                column(3, checkboxInput(ns("filter_na"), tags$b("Filter unmapped leaves"), value = F)),
                column(3, downloadButton(ns("download_tree_plot"), "Download tree", class = "btn-primary btn_rightAlign"))
            )
        ),
        tags$b("Lineage tree can be colored by expression of gene X as determined from single cell data (choose gene from 'Color by' drop down). "),
        fluidRow(
            column(12, uiOutput(ns("tree_view_ui")))
        )
    )
}


#' @export
tree_server <- function(input, output, session){
    output$tree_view_ui <- renderUI({
        ns <- session$ns
        plotOutput(ns("tree_view"), width = "100%", height = paste0(input$plot_height, "px"))
    })
    
    lin_tree_plot <- reactive({
        req(input$lin_group %in% c("br_time", gene_tbl[[1]]))
        if(input$lin_group == "br_time") {
            cby <- "br_time"
            cby_name <- "Birth time"
        } else {
            curg <- input$lin_group
            tree_tbl$sc_expr <- lin_sc_expr[curg,][match(tree_tbl$to, colnames(lin_sc_expr))]
            cby <- "sc_expr"
            cby_name <- curg
        }
        
        if(input$tree_root %in% avail_nodes) {
            cur_tree <- tree_tbl %>% filter(to %in% find_child_root(input$tree_root, tree_tbl))
            cur_root <- input$tree_root
        } else {
            cur_tree <- tree_tbl
            cur_root <- "P0"
        }
        
        assign("cur_tree", cur_tree, env =.GlobalEnv)
        
        cur_tree <- cur_tree %>% filter(br_time <= input$max_time)
        tBy <- input$tree_label_col
        if(input$tree_label_col == "lin_cell") {
            cur_tree$lin_cell <- paste0(cur_tree$to, ifelse(is.na(cur_tree$Cell), "", paste0("; ", cur_tree$Cell)))
        } else if(input$tree_label_col == "lin_cell_des") {
            cur_tree$lin_cell_des <- paste0(cur_tree$to, ifelse(is.na(cur_tree$Cell), "", paste0("; ", cur_tree$Cell)), ifelse(is.na(cur_tree$Description), "", paste0("; ", cur_tree$Description)))
        } 
        
        tree_td <- as.treedata(cur_tree)
        if(input$filter_na && cby != "br_time") {
            tip_to_drop <- fortify(tree_td)%>% filter(isTip & is.na(sc_expr))
            tree_td<-as.treedata(cur_tree %>% filter(!to %in% tip_to_drop$label))
            tip_to_drop2 <- fortify(tree_td)%>% filter(isTip & is.na(sc_expr)) # this is unfortunately the only working approach, as treeio has bugs in drop.tip.treedata
            tree_td<-as.treedata(cur_tree %>% filter(!to %in% c(tip_to_drop$label, tip_to_drop2$label)))
        }
        
        if(input$tree_style != "arc-bar") {
            p1<-ggtree(tree_td, branch.length='lifetime', aes_string(color = cby), size = 2, ladderize = F, layout = input$tree_style)
            p1 <- lineage_tree_flip(p1, silence = T)
            assign("p1", p1, env = .GlobalEnv)
            p_final <- p1 + 
                geom_tiplab(aes_string(label = tBy), size = 3,color = "black", align = F) + 
                geom_text(aes(label=label), data = p1$data[p1$data$br_time < 80,], size = 5, angle = 90,color = "black") + 
                scale_color_gradientn(colours = get_numeric_color(input$lin_pal), na.value=input$na_color) + 
                labs(color = cby_name)+
                scale_y_continuous(expand = c(0, 3))+
                theme(axis.text=element_text(size=12))  +
                theme(legend.position="top")
            if(cur_root == "P0") {
                p_final <- p_final + theme_tree2(legend.position = "top") + 
                    scale_x_continuous(sec.axis = dup_axis())
                    #scale_x_continuous(sec.axis = dup_axis(), lim = c(0,input$max_time + 100))
            } else { # time does not match axis so for now drop them
                p_final <- p_final + xlab(NULL) + ylab(NULL)  #+ 
                    #scale_x_continuous(lim = c(0,input$max_time + 100))
            }
            return(p_final)
        } else {
            # ggraph
            tbl_graph <- as_tbl_graph(tree_td@phylo) %>% activate(nodes) %>% left_join(cur_tree, by = c("name" = "to")) %>%
                mutate(text.size = 1) %>%
                mutate(label = ifelse(br_time > 80, "", name))
            plotGraph(tbl_graph, color.by=cby, pal="viridis", label="label", type = "numeric",border.size=.3, background="white", na_color = input$na_color, legend.title = cby_name) + 
                geom_node_text(aes(x = x*1.035, y=y*1.035, filter = leaf, label = name,angle = -((-node_angle(x, y)+90)%%180)+90), hjust='outward', size = 3, show.legend = FALSE) + 
                theme(text=element_text(family = "Helvetica", size=10),
                      axis.ticks.x=element_blank(),
                      axis.text.x=element_blank(),
                      axis.ticks.y=element_blank(),
                      axis.text.y=element_blank(),
                      legend.text=element_text(color="black"),
                      legend.title = element_text(colour="black"),
                      legend.margin=margin(0,0,0,0),
                      legend.box.margin=margin(-10,-10,-10,-10),
                      plot.margin = unit(c(.3,.5,.3,.3), "cm"))                                                              
        }
    })
    
    output$tree_view <- renderPlot({
        req(lin_tree_plot())
        lin_tree_plot()
    })
    
    output$download_tree_plot <- downloadHandler(
        filename = function(format = "pdf") {
            paste('Lineage_tree-', input$lin_group, "_", Sys.Date(), ".pdf", sep='')
        },
        content = function(con, format = input$plotf) {
            fn_dev<-"pdf"
            req(lin_tree_plot())
            ggsave(con, plot = lin_tree_plot(), device = fn_dev, width = 8, height = input$plot_height/100, limitsize = FALSE)
            shut_device <- dev.list()[which(names(dev.list()) != "quartz_off_screen")]
            if(length(shut_device)) dev.off(which = shut_device) # Make sure ggsave does not change graphic device
        }
    )
}








































