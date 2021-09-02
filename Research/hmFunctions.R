#' getHeatAnno returns a heatmap annotation object where you input the list of variables you would like
#'
#' @param colorValues the variable name paired with color for each sample by variable
#' @param anno_df the dataframe that ou are annotating, these are the filtered columns of the varColumns
getHeatAnno <- function(colorValues,anno_df){
    ha <- ComplexHeatmap::HeatmapAnnotation(
        df = anno_df,
        col = colorValues,
        na_col = "white",
        annotation_legend_param = list(
            title_gp = gpar(
                fontsize = 14,
                lineheight = 2,
                fontface = "bold"
            ),
            labels_gp = gpar(fontsize = 14, lineheight = 2),
            grid_height = unit(0.8, "cm"),
            grid_width = unit(0.8, "cm"),
            gt_render = list(padding = unit(c(2, 2, 2, 2), "cm"))
        ),
        annotation_name_gp = gpar(
            fontsize = 14,
            fontface = "bold",
            lineheight = 2,
            gt_render = list(padding = unit(c(2, 2, 2, 2), "cm"))
        ),
        annotation_name_side = "left"
    )
    par(c(1.5, 1.5, 1.5, 2.5) + 0.5)
    return(ha)
}

#' assignColors returns a ComplexHeatmap annotation object with a list of variable columns you have
#'
#' @param targets your dataframe of sample IDs, and variables i.e. Age, Sex, Tumor, etc.
#' @param varColumns string name of the columns you want to create colors for in the complexheatmap within targets
#' @param col_vect NULL is default, these are all the colors you want to use in you heatmap annotations
assignColors <- function(targets, varColumns = c("Type", "Grade"), col_vect = NULL) {
    if (is.null(col_vect)) {col_vect <- glasbey()}
    dimnames(targets)[[2]]
    dat <- targets[,varColumns]
    anno_df <- data.frame(dat)
    vars2Color <- as.list(lapply(dat, unique))
    colorValues <-lapply(vars2Color, function(x) {x = (col_vect)[1:(length(x))]})
    for (x in 1:length(vars2Color)) {
        for (varNum in 1:length(vars2Color[x])) {
            names(colorValues[x][[1]]) = c(vars2Color[x][[1]])
        }
    }
    if(all(names(colorValues)==names(anno_df))){
        return(getHeatAnno(colorValues,anno_df))
    }
}

## Define the Heatmap ----------------------
col_fun <- circlize::colorRamp2(c(0, 0.25, 0.5, 0.75, 1), c("darkblue","deepskyblue", "white", "tomato","red"))

drawHeatMap <- function(yourHeatMap) {
    return(
        draw(
            yourHeatMap,
            merge_legend = T,
            ht_gap = unit(3, "cm"),
            heatmap_legend_side = "left",
            annotation_legend_side = "left",
            padding = unit(c(20, 20, 20, 20), "mm")
        )
    )
}

getHeatMap <- function(betaRanges, titleValue, ha){
    titleOfPlot <- paste("Heatmap of",titleValue,sep = " ")
    hmTopNumbers <- Heatmap(
        betaRanges,
        col = col_fun,  ## Define the color scale
        cluster_columns = T,  ## Cluster the columns
        #cluster_rows = rowcluster,
        show_column_names = T,  ## Show the Column Names (which is sample #)
        column_names_gp = gpar(fontsize = 16),  ## Column Name Size
        show_row_names = F,  ## Show Row names (which is probes)
        row_names_side = "left",
        row_title_side = "left",
        row_names_gp = gpar(fontsize = 10),
        row_title_gp = gpar(fontsize = 16, fontface = "bold"),
        show_row_dend = F,
        show_column_dend = T,
        use_raster=T,
        show_heatmap_legend = T,
        top_annotation = ha,
        column_title = titleOfPlot,
        column_title_gp = gpar(fontsize = 24,fontface = "bold"),
        raster_device = "png",
        raster_quality = 3,
        heatmap_legend_param = list(
            title = "Beta Value",
            #legend_height = unit(5, "in"),
            labels_gp = gpar( fontsize = 14),
            title_gp = gpar(fontsize = 18, fontface = "bold"),
            legend_direction = "vertical",
            heatmap_legend_side = "right", annotation_legend_side = "right",
            legend_height =  unit(2.5, "in")
        )
    )
    return(drawHeatMap(hmTopNumbers))
}
