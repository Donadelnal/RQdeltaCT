

#' @title read_Ct_wide
#'
#' @description
#' This function imports Ct data in a wide-format table with sample names given in columns.
#'
#' @details
#' This function needs two files to import: a wide-format table with Ct values and an additional file with group names
#' (see parameters path.Ct.file and path.design.file for further details on tables structure).
#' Both files are merged to return a long-format table ready for analysis.
#' All parameters must be specified; there are no default values.
#'
#' @param path.Ct.file Path to wide-format table in .txt or csv. format with Ct values.
#' This table must contain gene names in the first column, and sample names in the first row (genes by rows and samples by columns).
#' @param path.design.file Path to table in .txt  or csv. file that contains two columns: column named "Sample" with names of samples
#' and column named "Group" with names of groups assigned to samples. The names of samples in this file must
#' correspond to the names of columns in the file with Ct values.
#' @param sep Character of a field separator in both imported files.
#' @param dec Character used for decimal points in Ct values.
#'
#' @return Data frame in long format ready to analysis.
#' @export
#'
#' @examples
#' path.Ct.file <- system.file("extdata",
#'                             "data_Ct_wide.txt",
#'                             package = "RQdeltaCT")
#' path.design.file <- system.file("extdata",
#'                                 "data_design.txt",
#'                                 package = "RQdeltaCT")
#'
#' library(tidyverse)
#' data.Ct <- read_Ct_wide(path.Ct.file = path.Ct.file,
#'                         path.design.file = path.design.file,
#'                         sep ="\t",
#'                         dec = ".")
#' str(data.Ct)
#'
#' @importFrom utils read.csv
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @import tidyverse
#'
read_Ct_wide <- function(path.Ct.file,
                         path.design.file,
                         sep,
                         dec) {
  data_wide <- read.csv(path.Ct.file,
                        header = TRUE,
                        sep = sep,
                        dec = dec)

  data_wide_design <- read.csv(path.design.file,
                               header = TRUE,
                               sep = sep)

  colnames(data_wide)[1] <- "Gene"
  data_wide <- mutate(data_wide,
                      across(everything(),
                             as.character))
  data_slim <- pivot_longer(data_wide, -Gene,
                            names_to = "Sample",
                            values_to = "Ct")
  data_slim[, "Group"] <- NA

  for (x in 1:nrow(data_wide_design)) {
    index <- which(data_slim$Sample == data_wide_design$Sample[x])
    data_slim$Group[index] <- data_wide_design$Group[x]
  }
  return(data_slim)
}






#' @title read_Ct_long
#'
#' @description
#' Imports a long-format table with Ct values.
#'
#' @param path Path to a .txt or csv. file with long-type table with Ct values. This table must contain at least 4 columns, separately for
#' sample names, gene names, Ct values and group names (these columns will be imported by this function).
#' Imported table could also contain a column with a flag information,
#' which could be optionally imported (see add.col.Flag and col.Flag parameters).
#' @param sep Character of a field separator in imported file.
#' @param dec Character used for decimal points in Ct values.
#' @param skip Integer: number of lines of the data file to skip before beginning to read data. Default to 0.
#' @param column.Sample Integer: number of column with sample names.
#' @param column.Gene Integer: number of column with gene names.
#' @param column.Ct Integer: number of column with Ct values.
#' @param column.Group Integer: number of column with group names.
#' @param add.column.Flag Logical: if data contains a column with flag information which should also be imported,
#' this parameter should be set to TRUE. Default to FALSE.
#' @param column.Flag Integer: number of column with flag information. Should be specified if add.col.Flag = TRUE.
#' This column should contain a character-type values (e.g. "Undetermined" and "OK"), however,
#' other types of values are allowed (e.g. numeric), but must be converted to character or factor
#' after importing data (see examples).
#'
#' @return Data.frame in long format ready for analysis.
#' @export
#'
#' @examples
#' path <- system.file("extdata",
#'                     "data_Ct_long.txt",
#'                     package = "RQdeltaCT")
#'
#' library(tidyverse)
#' data.Ct <- read_Ct_long(path = path,
#'                         sep = "\t",
#'                         dec = ".",
#'                         skip = 0,
#'                         add.column.Flag = TRUE,
#'                         column.Sample = 1,
#'                         column.Gene = 2,
#'                         column.Ct = 5,
#'                         column.Group = 9,
#'                         column.Flag = 4)
#' str(data.Ct)
#'
#' data.Ct <- mutate(data.Ct,
#'                   Flag = ifelse(Flag < 1, "Undetermined", "OK"))
#' str(data.Ct)
#'
#' @importFrom utils read.csv
#'
read_Ct_long <- function(path,
                         sep,
                         dec,
                         skip = 0,
                         column.Sample,
                         column.Gene,
                         column.Ct,
                         column.Group,
                         add.column.Flag = FALSE,
                         column.Flag) {
  data <- read.csv(
    path,
    header = TRUE,
    sep = sep,
    dec = dec,
    skip = skip
  )

  if (add.column.Flag == FALSE) {
    data <- data[, c(column.Sample,
                     column.Gene,
                     column.Ct,
                     column.Group)]
    colnames(data) <- c("Sample", "Gene", "Ct", "Group")

  } else {
    data <- data[, c(column.Sample,
                     column.Gene,
                     column.Ct,
                     column.Group,
                     column.Flag)]
    colnames(data) <- c("Sample", "Gene", "Ct", "Group", "Flag")
  }

  return(data)
}







#' @title control_Ct_barplot_sample
#'
#' @description
#' Sample-wide quality control of raw Ct values by illustrating the numbers of Ct values labelled as reliable or not by using reliability criteria (see function parameters).
#'
#' @details
#' This function labels Ct values as reliable or not using given reliability criteria, counts them, and presents them graphically.
#' Results are useful to identify samples with low numbers of reliable Ct values. This function does not perform data filtering.
#'
#' @param data Object returned from read_Ct_long() or read_Ct_wide() function,
#' or data frame containing column named "Sample" with sample names, column named "Gene" with gene names,
#' column named "Ct" with raw Ct values, and column named "Group" with group names.
#' Optionally, data frame could contain column named "Flag" with flag information
#' (e.g. "Undetermined" and "OK"), which will be used for reliability assessment.
#' @param flag.Ct Character of a flag used for undetermined Ct values. Default to "Undetermined".
#' @param maxCt Numeric, a maximum of Ct value allowed. Default to 35.
#' @param flag Character of a flag used in the Flag column for values which are unreliable. Default to "Undetermined".
#' @param colors Character vector length of two, containing colors for Ct values that were labelled as reliable (first element of vector) or not (second element of vector).
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param x.axis.title Character: title of x axis. Default to "".
#' @param y.axis.title Character: title of y axis. Default to "Number".
#' @param legend.title Character: title of legend. Default to "Reliable Ct value?".
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param legend.title.size Integer: font size of legend title.  Default to 12.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be set to "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "Ct_control_barplot_for_samples".
#'
#' @return List containing plot and table with counts of reliable and no reliable Ct values in samples.
#' Additional information about returned table is also printed to help the user to properly interpret returned table.
#' Plot is also displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' sample.Ct.control <- control_Ct_barplot_sample(data.Ct)
#' sample.Ct.control[[2]]
#'
#' @importFrom dplyr mutate arrange filter rename desc
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr %>%
#' @importFrom stats reorder
#' @importFrom ggplot2 ggplot geom_bar coord_flip scale_fill_manual xlab ylab labs theme_classic theme element_text scale_x_discrete ggsave
#' @import ggplot2
#' @import tidyverse
#'
control_Ct_barplot_sample <- function(data,
                                      flag.Ct = "Undetermined",
                                      maxCt = 35,
                                      flag = "Undetermined",
                                      colors = c("#66c2a5", "#fc8d62"),
                                      x.axis.title = "",
                                      y.axis.title = "Number",
                                      axis.title.size = 11,
                                      axis.text.size = 10,
                                      plot.title = "",
                                      plot.title.size = 14,
                                      legend.title = "Reliable Ct value?",
                                      legend.title.size = 11,
                                      legend.text.size = 11,
                                      legend.position = "top",
                                      save.to.tiff = FALSE,
                                      dpi = 600,
                                      width = 15,
                                      height = 15,
                                      name.tiff = "Ct_control_barplot_for_samples") {
  data$Ct[data$Ct == flag.Ct] <- 41
  data$Ct <- as.numeric(data$Ct)

  if (sum(colnames(data) %in% "Flag") > 0) {
    data <-
      mutate(data, Reliable = ifelse(Ct > maxCt |
                                       Flag == flag, yes = "No",  no = "Yes"))

  } else {
    data <-
      mutate(data, Reliable = ifelse(Ct > maxCt , yes = "No",  no = "Yes"))
  }

  bar <- as.data.frame(table(data$Reliable, data$Sample))
  order <- arrange(filter(bar, Var1 == "Yes"), desc(Freq))$Var2

  barplot.samples <-
    ggplot(bar, aes(
      x = reorder(Var2, desc(Freq)),
      y = Freq,
      fill = Var1
    )) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(breaks = c("Yes", "No"),
                      values = c("Yes" = colors[1], "No" = colors[2])) +
    xlab(x.axis.title) +
    ylab(y.axis.title) +
    labs(fill = legend.title, title = plot.title) +
    theme_classic() +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, color = 'black')) +
    theme(axis.title = element_text(size = axis.title.size, color = 'black')) +
    theme(legend.title = element_text(size = legend.title.size, colour =
                                        "black")) +
    theme(legend.text = element_text(size = legend.text.size, colour = "black")) +
    theme(plot.title = element_text(size = plot.title.size)) +
    scale_x_discrete(limits = order)

  print(barplot.samples)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      barplot.samples,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  } else {

  }

  tab <- table(data$Reliable, data$Sample)
  tab <- tab %>%
    as.data.frame() %>%
    pivot_wider(names_from = Var1, values_from = Freq)

  if (sum(colnames(tab) %in% "No") > 0) {

    tab <- tab %>%
      arrange(desc(No)) %>%
      mutate(Not.reliable.fraction = No / (No + Yes)) %>%
      rename(Sample = Var2,
             Not.reliable = No,
             Reliable = Yes)
  } else {

    tab <- tab %>%
      arrange(Yes) %>%
      rename(Sample = Var2,
             Reliable = Yes)
    message("Note: All Ct values were labeled as 'Reliable'.")
  }

  return(list(barplot.samples, tab))
}









#' @title control_Ct_barplot_gene
#'
#' @description
#' Gene-wide quality control of raw Ct values across groups by illustrating numbers of Ct values labelled as reliable or not by using reliability criteria (see function parameters).
#'
#' @details
#' This function does not perform data filtering, but only counts Ct values labelled as reliable or not and presents them graphically.
#' Could be useful to identify genes with low number of reliable Ct values.
#'
#' @param data Object returned from read_Ct_long() or read_Ct_wide() function,
#' or data frame containing column named "Sample" with sample names, column named "Gene" with gene names,
#' column named "Ct" with raw Ct values, column named "Group" with group names.
#' Optionally, data frame can contain column named "Flag" with flag information (e.g. "Undetermined" and "OK"),
#' which will be used for reliability assessment.
#' @param flag.Ct Character of a flag used for undetermined Ct values. Default to "Undetermined".
#' @param maxCt Numeric, a maximum of Ct value allowed. Default to 35.
#' @param flag Character of a flag used in the Flag column for values which are unreliable. Default to "Undetermined".
#' @param colors Character vector length of two, containing colors for Ct values which are labelled as reliable (first element of vector) or not (second element of vector).
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param x.axis.title Character: title of x axis. Default to "".
#' @param y.axis.title Character: title of y axis. Default to "Number".
#' @param legend.title Character: title of legend. Default to "Reliable Ct value?".
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param legend.title.size Integer: font size of legend title.  Default to 12.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file.  Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "Ct_control_barplot_for_genes".
#'
#' @return List containing plot and table with counts of reliable and non reliable Ct values in genes.
#' Additional information about returned table is also printed to help the user to properly interpret returned table.
#' Plot is also displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' gene.Ct.control <- control_Ct_barplot_gene(data.Ct)
#' gene.Ct.control[[2]]
#'
#' @importFrom dplyr mutate arrange filter rename select desc
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom magrittr %>%
#' @importFrom tidyselect any_of
#' @importFrom stats reorder
#' @importFrom ggplot2 ggplot geom_bar coord_flip scale_fill_manual xlab ylab labs theme_classic theme element_text scale_x_discrete facet_wrap ggsave
#' @import ggplot2
#' @import tidyverse
#'
control_Ct_barplot_gene <- function(data,
                                    flag.Ct = "Undetermined",
                                    maxCt = 35,
                                    flag = "Undetermined",
                                    colors = c("#66c2a5", "#fc8d62"),
                                    x.axis.title = "",
                                    y.axis.title = "Number",
                                    axis.title.size = 11,
                                    axis.text.size = 10,
                                    legend.title = "Reliable Ct value?",
                                    legend.title.size = 11,
                                    legend.text.size = 11,
                                    legend.position = "top",
                                    plot.title = "",
                                    plot.title.size = 14,
                                    save.to.tiff = FALSE,
                                    dpi = 600,
                                    width = 15,
                                    height = 15,
                                    name.tiff = "Ct_control_barplot_for_genes") {
  data$Ct[data$Ct == flag.Ct] <- 41
  data$Ct <- as.numeric(data$Ct)

  if (sum(colnames(data) %in% "Flag") > 0) {
    data <-
      mutate(data, Reliable = ifelse(Ct > maxCt |
                                       Flag == flag, yes = "No",  no = "Yes"))
  } else {
    data <-
      mutate(data, Reliable = ifelse(Ct > maxCt , yes = "No",  no = "Yes"))
  }

  bar <-
    as.data.frame(table(data$Reliable, data$Gene, data$Group))
  order <- arrange(filter(bar, Var1 == "Yes"), Freq)$Var2

  barplot.genes <-
    ggplot(bar, aes(
      x = reorder(Var2, desc(Freq)),
      y = Freq,
      fill = Var1
    )) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(breaks = c("Yes", "No"),
                      values = c("Yes" = colors[1], "No" = colors[2])) +
    xlab(x.axis.title) + ylab(y.axis.title) +
    labs(fill = legend.title, title = plot.title) +
    theme_classic() + theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, color = 'black')) +
    theme(axis.title = element_text(size = axis.title.size, color = 'black')) +
    theme(legend.title = element_text(size = legend.title.size, colour =
                                        "black")) +
    theme(legend.text = element_text(size = legend.text.size, colour =
                                       "black")) +
    theme(plot.title = element_text(size = plot.title.size)) +
    scale_x_discrete(limits = rev(unique(bar$Var2))) +
    facet_wrap(vars(Var3))

  print(barplot.genes)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      barplot.genes,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }

  tab <- table(data$Reliable, data$Gene, data$Group)
  tab <- tab %>%
    as.data.frame() %>%
    pivot_wider(names_from = Var1, values_from = Freq)

  if (sum(colnames(tab) %in% "No") > 0) {

    tab <- tab %>%
      arrange(desc(No)) %>%
      mutate(Not.reliable.fraction = No / (No + Yes)) %>%
      rename(Gene = Var2,
             Group = Var3,
             Not.reliable = No,
             Reliable = Yes)
  } else {

    tab <- tab %>%
      arrange(Yes) %>%
      rename(Gene = Var2,
             Group = Var3,
             Reliable = Yes)
    message("Note: All Ct values were labeled as 'Reliable'.")
  }

  return(list(barplot.genes, tab))
}







#' @title control_heatmap
#'
#' @description
#' This function creates heatmap with number of replicates for each gene across samples.
#'
#' @param data Object returned from read_Ct_long() or read_Ct_wide() function,
#' or data frame containing at least column named "Sample" with sample names, column named "Gene" with gene names.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all genes.
#' @param display.numbers Logical: if TRUE (default), numeric values will be printed to the cells.
#' @param show.legend Logical: if TRUE, legend will be shown. Default to FALSE.
#' @param dist.row,dist.col Character: name of method used for calculation of distances between rows or columns, derived from stats::dist() function,
#' must be one of "euclidean" (default) , "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param clust.method Character: name of used method for agglomeration, derived from stats::hclust() function,
#' must be one of "ward.D", "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".
#' @param colors Vector with colors used to fill created heatmap.
#' @param show.colnames,show.rownames Logical: of TRUE, names of columns (sample names) and rows (gene names) will be shown. Both default to TRUE.
#' @param border.color Character: color of cell borders on heatmap. If set to NA (default) no border will be drawn.
#' @param fontsize Numeric: global fontsize of heatmap. Default to 10.
#' @param fontsize.col,fontsize.row Numeric: fontsize of colnames and rownames. Default to 10.
#' @param angle.col Integer: angle of the column labels, one of the 0, 45, 90, 270, and 315.
#' @param cellwidth,cellheight Numeric: width and height of individual cell. Both default to NA.
#' These parameters are useful in situations where margins are too small and the plot is cropped (column names and annotation legend are sometimes partially hidden).
#' Specification of this parameter allows to adjust size of the plot and solve this problem.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "heatmap_results".
#'
#' @return Heatmap with hierarchical clustering, displayed on the graphic device (if save.to.tiff = FALSE)
#' or saved to .tiff file (if save.to.tiff = TRUE).
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(pheatmap)
#' data(data.Ct)
#' # Vector of colors to fill the heatmap can be specified to fit the user needings:
#' colors <- c("#4575B4","#FFFFBF","#C32B23")
#' control_heatmap(data.Ct,
#'                 sel.Gene = "all",
#'                 colors = colors,
#'                 show.colnames = TRUE,
#'                 show.rownames = TRUE,
#'                 fontsize = 11,
#'                 fontsize.row = 11)
#'
#' @importFrom stats hclust dist
#' @importFrom grDevices tiff
#' @importFrom dplyr select ungroup
#' @importFrom tidyselect any_of
#' @importFrom pheatmap pheatmap
#' @import tidyverse
#'
control_heatmap <- function(data,
                            sel.Gene = "all",
                            display.numbers = TRUE,
                            dist.row = "euclidean",
                            dist.col = "euclidean",
                            clust.method = "average",
                            colors = c(
                              "#4575B4",
                              "#FFFFBF",
                              "#C32B23"),
                            show.colnames = TRUE,
                            show.rownames = TRUE,
                            border.color = NA,
                            fontsize = 10,
                            fontsize.col = 10,
                            fontsize.row = 10,
                            angle.col = 0,
                            cellwidth = NA,
                            cellheight = NA,
                            show.legend = FALSE,
                            save.to.tiff = FALSE,
                            dpi = 600,
                            width = 15,
                            height = 15,
                            name.tiff = "heatmap_results") {
  if (sel.Gene[1] == "all") {
    data <- data

  } else {
    data <- select(data, Group, Sample, any_of(sel.Gene))
  }

  data <- as.data.frame(data)
  data<- table(data$Sample, data$Gene)

  if (save.to.tiff == TRUE) {

    device <- .Device

    if (device == "RStudioGD") {
      dev.off()
    }

    tiff(
      paste(name.tiff, ".tiff", sep = ""),
      res = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )

    pheatmap(
      data,
      display_numbers = display.numbers,
      number_format = "%.0f",
      clustering_method = clust.method,
      clustering_distance_cols = dist.col,
      clustering_distance_rows = dist.row,
      show_colnames = show.colnames,
      show_rownames = show.rownames,
      border_color = border.color,
      fontsize = fontsize,
      fontsize_col = fontsize.col,
      fontsize_row = fontsize.row,
      color = colors,
      angle_col = angle.col,
      cellwidth = cellwidth,
      cellheight = cellheight,
      legend = show.legend
    )

    dev.off()

  } else {

    pheatmap(
      data,
      display_numbers = display.numbers,
      number_format = "%.0f",
      clustering_method = clust.method,
      clustering_distance_cols = dist.col,
      clustering_distance_rows = dist.row,
      show_colnames = show.colnames,
      show_rownames = show.rownames,
      border_color = border.color,
      fontsize = fontsize,
      fontsize_col = fontsize.col,
      fontsize_row = fontsize.row,
      color = colors,
      angle_col = angle.col,
      cellwidth = cellwidth,
      cellheight = cellheight,
      legend = show.legend
    )
  }
}






#' @title filter_Ct
#'
#' @description
#' This function filters Ct data according to the used filtering criteria (see parameters).
#'
#' @param data Object returned from read_Ct_long() or read_Ct_wide() function,
#' or data frame containing column named "Sample" with sample names, column named "Gene" with gene names,
#' column named "Ct" with raw Ct values, column named "Group" with group names.
#' Optionally, data frame can contain column named "Flag" with flag information (e.g. "Undetermined" and "OK"),
#' which will be used for filtering.
#' @param flag.Ct Character of a flag used for undetermined Ct values, default to "Undetermined".
#' @param maxCt Numeric, a maximum of Ct value allowed.
#' @param flag Character: flag used in Flag column for values which should be filtered out, default to "Undetermined".
#' @param remove.Gene Character: vector with names of genes to remove from data
#' @param remove.Sample Character: vector with names of samples to remove from data
#' @param remove.Group Character: vector with names of groups to remove from data
#'
#' @return Data frame with filtered data.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#'
#' dim(data.Ct)
#' dim(data.CtF)
#'
#' @importFrom dplyr filter
#' @import tidyverse
#'
filter_Ct <- function(data,
                      flag.Ct = "Undetermined",
                      maxCt = 35,
                      flag = c("Undetermined"),
                      remove.Gene = c(""),
                      remove.Sample = c(""),
                      remove.Group = c("")) {
  data <- filter(data, Ct != flag.Ct)
  data$Ct <- as.numeric(data$Ct)
  data <- filter(
    data,
    Ct <= maxCt,
    !Gene %in% remove.Gene,
    !Sample %in% remove.Sample,
    !Group %in% remove.Group
  )

  if (sum(colnames(data) %in% "Flag") > 0) {
    data <- filter(data, !Flag %in% flag)
  }

  return(data)
}





#' @title make_Ct_ready
#'
#' @description
#' This function collapses technical replicates (if present in data) by means and (optionally) imputes missing data by means calculated separately for each group.
#' This function also prepares Ct data for further steps of analysis.
#'
#' @param data Data object returned from read_Ct_long(), read_Ct_wide() or filter_Ct() function,
#' or data frame containing column named "Sample" with sample names, column named "Gene" with gene names,
#' column named "Ct" with raw Ct values (must be numeric), column named "Group" with group names. Presence of any other columns is allowed,
#' but they will not be used by this function.
#' @param imput.by.mean.within.groups Logical: if TRUE, missing values will be imputed by means calculated separately for each group.
#' This parameter can influence results, thus to draw more user attention on this parameter, no default value was set.
#' @param save.to.txt Logical: if TRUE, returned data will be saved to .txt file. Default to FALSE.
#' @param name.txt Character: name of saved .txt file, without ".txt" name of extension. Default to "Ct_ready".
#'
#' @return Data frame with prepared data. Information about number and percentage of missing values is also printed.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#'data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#'head(data.CtF.ready)
#'
#' @importFrom utils write.table
#' @importFrom dplyr select summarise group_by mutate across
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect where
#' @import tidyverse
#'
make_Ct_ready <- function(data,
                          imput.by.mean.within.groups,
                          save.to.txt = FALSE,
                          name.txt = "Ct_ready") {
  data <- data %>%
    group_by(Group, Gene, Sample) %>%
    summarise(mean = mean(Ct, na.rm = TRUE), .groups = "keep") %>%
    as.data.frame()

  data_wide <- data %>%
    select(Group, Sample, Gene, mean) %>%
    pivot_wider(names_from = Gene, values_from = mean)

  nas <- sum(is.na(data_wide))
  percentage <-
    sum(is.na(data_wide)) / ((ncol(data_wide) - 2) * nrow(data_wide))

  if (imput.by.mean.within.groups == TRUE) {
    data_wide_imp <- data_wide %>%
      group_by(Group) %>%
      mutate(across(where(is.numeric), ~ replace(., is.na(.), mean(., na.rm = TRUE))))

    message(
      "The data contain ",
      nas,
      " missing values that constitute ",
      round(percentage * 100, 5),
      " percent of the total data.\nMissing values were imputed using means within compared groups.\n"
    )

    if (save.to.txt == TRUE) {
      write.table(as.data.frame(data_wide_imp),
                  paste(name.txt, ".txt", sep = ""))
    }
    return(data_wide_imp)

  } else {
    message(
      "The data contain ",
      nas,
      " missing values that constitute ",
      round(percentage * 100, 5),
      " percent of the total data."
    )

    if (save.to.txt == TRUE) {
      write.table(as.data.frame(data_wide),
                  paste(name.txt, ".txt", sep = ""))
    }
    return(data_wide)
  }
}













#' @title norm_finder
#'
#' @description
#' This function calculates stability scores using NormFinder algorithm (https://www.moma.dk/software/normfinder)
#' published in https://aacrjournals.org/cancerres/article/64/15/5245/511517/Normalization-of-Real-Time-Quantitative-Reverse.
#' This function is internally used by other RQdeltaCT package function, find_ref_gene(); therefore norm_finder() function does not need to be used separately.
#'
#' @param data Object returned from make_Ct_ready() functions.
#' @param candidates Character: vector of names of genes - candidates for reference gene.
#' @param save.to.txt Logical: if TRUE, returned table with results will be saved to .txt file. Default to FALSE.
#' @param name.txt Character: name of saved .txt file, without ".txt" name of extension. Default to "norm_finder_results".
#'
#' @return Table with calculated stability score; the lowest value the best candidate for reference gene.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' reference.stability.nF <- norm_finder(data.CtF.ready,
#'                                       candidates = c("Gene4",
#'                                                      "Gene8",
#'                                                      "Gene10",
#'                                                      "Gene16",
#'                                                      "Gene17",
#'                                                      "Gene18"))
#' @importFrom dplyr filter select
#' @importFrom utils write.table
#' @importFrom tidyr pivot_longer pivot_wider
#' @import tidyverse
#'
norm_finder <- function(data,
                        candidates,
                        save.to.txt = FALSE,
                        name.txt = "NormFinder_results") {
  if (sum(duplicated(data$Sample)) != 0) {
    data$Sample <- paste(data$Sample, data$Group, sep = "_")
  }
  data.t <- data %>%
    pivot_longer(!c(Sample, Group),
                 names_to = "Gene",
                 values_to = "Ct") %>%
    filter(Gene %in% candidates) %>%
    pivot_wider(
      id_cols = !c(Group),
      names_from = "Sample",
      values_from = "Ct"
    )

  last.row <- c("Group", data$Group)
  dat0 <- rbind(as.data.frame(data.t), last.row)
  rownames(dat0) <- dat0[, "Gene"]
  dat0 <- select(dat0, -Gene)

  ntotal = dim(dat0)[2]
  k0 = dim(dat0)[1]
  ngenes = k0 - 1
  genenames = rownames(dat0)[-k0]
  grId = dat0[k0,]
  dat0 = dat0[-k0,]

  dat = matrix(as.numeric(unlist(dat0)), ngenes, ntotal)

  samplenames = colnames(dat0)
  grId = factor(unlist(grId))
  groupnames = levels(grId)
  ngr = length(levels(grId))
  nsamples = rep(0, ngr)

  for (group in 1:ngr) {
    nsamples[group] = sum(grId == groupnames[group])
  }

  MakeStab = function(da) {
    ngenes = dim(da)[1]
    sampleavg = apply(da, 2, mean)
    genegroupavg = matrix(0, ngenes, ngr)

    for (group in 1:ngr) {
      genegroupavg[, group] = apply(da[, grId == groupnames[group]], 1, mean)
    }

    groupavg = rep(0, ngr)

    for (group in 1:ngr) {
      groupavg[group] = mean(da[, grId == groupnames[group]])
    }

    GGvar = matrix(0, ngenes, ngr)

    for (group in 1:ngr) {
      grset = (grId == groupnames[group])
      a = rep(0, ngenes)

      for (gene in 1:ngenes) {
        a[gene] = sum((da[gene, grset] - genegroupavg[gene, group] -
                         sampleavg[grset] + groupavg[group]) ^ 2) / (nsamples[group] -
                                                                       1)
      }
      GGvar[, group] = (a - sum(a) / (ngenes * ngenes - ngenes)) / (1 -
                                                                      2 / ngenes)
    }

    genegroupMinvar = matrix(0, ngenes, ngr)

    for (group in 1:ngr) {
      grset = (grId == groupnames[group])
      z = da[, grset]

      for (gene in 1:ngenes) {
        varpair = rep(0, ngenes)

        for (gene1 in 1:ngenes) {
          varpair[gene1] = var(z[gene,] - z[gene1,])
        }
        genegroupMinvar[gene, group] = min(varpair[-gene]) / 4
      }
    }

    GGvar = ifelse(GGvar < 0, genegroupMinvar, GGvar)

    dif = genegroupavg
    difgeneavg = apply(dif, 1, mean)
    difgroupavg = apply(dif, 2, mean)
    difavg = mean(dif)

    for (gene in 1:ngenes) {
      for (group in 1:ngr) {
        dif[gene, group] = dif[gene, group] - difgeneavg[gene] - difgroupavg[group] +
          difavg
      }
    }

    nsampMatrix = matrix(rep(nsamples, ngenes), ngenes, ngr, byrow = T)
    vardif = GGvar / nsampMatrix
    gamma = sum(dif * dif) / ((ngr - 1) * (ngenes - 1)) - sum(vardif) /
      (ngenes * ngr)
    gamma = ifelse(gamma < 0, 0, gamma)

    difnew = dif * gamma / (gamma + vardif)
    varnew = vardif + gamma * vardif / (gamma + vardif)
    Ostab0 = abs(difnew) + sqrt(varnew)
    Ostab = apply(Ostab0, 1, mean)

    mud = rep(0, ngenes)
    for (gene in 1:ngenes) {
      mud[gene] = 2 * max(abs(dif[gene,]))
    }

    genevar = rep(0, ngenes)
    for (gene in 1:ngenes) {
      genevar[gene] = sum((nsamples - 1) * GGvar[gene,]) / (sum(nsamples) - ngr)
    }
    Gsd = sqrt(genevar)

    return(cbind(mud, Gsd, Ostab, rep(gamma, ngenes), GGvar, dif))
  }

  res = MakeStab(dat)
  ord = order(res[, 3])
  FinalRes = data.frame("Stability" = round(res[ord, 3], 2),
                        row.names = genenames[ord])

  if (save.to.txt == TRUE) {
    write.table(as.data.frame(FinalRes), paste(name.txt, ".txt", sep = ""))
  }

  return(FinalRes)
}






#' @title find_ref_gene
#'
#' @description
#' This function assess gene expression stability by calculation the following parameters: minimum, maximum, standard deviation, variance,
#' and stability measures from NormFinder and geNorm algorithms. It also presents C values graphically on a line plot. This function is helpful to select the best
#' reference gene for normalization of Ct values.
#'
#' @param data Object returned from make_Ct_ready() functions.
#' @param groups Character vector with names of groups used for analysis. If all groups should be included to the analysis, groups parameter should be set to "all".
#' @param candidates Character: vector of names of genes - candidates for gene reference.
#' @param colors Character: vector of colors for genes, the number of colors should be equal to the number of candidate genes.
#' @param norm.finder.score Logical: if TRUE, NormFinder stability score will be calculated. Default to TRUE.
#' @param genorm.score Logical: if TRUE, geNorm stability score will be calculated. Default to TRUE.
#' @param line.width Numeric: width of lines drawn in the plot. Default to 1.
#' @param angle Integer: value of angle in which names of genes are displayed. Default to 0.
#' @param x.axis.title Character: title of x axis. Default to "".
#' @param y.axis.title Character: title of y axis.  Default to "Ct".
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param legend.title Character: title of legend. Default to "".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "Ct_reference_gene_selection".
#' @param save.to.txt Logical: if TRUE, returned data will be saved to .txt file. Default to FALSE.
#' @param name.txt Character: name of saved .txt file, without ".txt" name of extension. Default to "ref_gene_results".
#'
#' @return List containing an object with plot and a table with calculated parameters. Created plot is also displayed on the graphic device.
#'
#' @export
#'
#' @examples
#'library(ctrlGene)
#'library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                        remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                        remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' ref <- find_ref_gene(data.CtF.ready,
#'                      groups = c("Disease","Control"),
#'                      candidates = c("Gene4", "Gene8","Gene10","Gene16","Gene17", "Gene18"),
#'                      col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "#1F77B4", "black"),
#'                      norm.finder.score = TRUE,
#'                      genorm.score = TRUE)
#' ref[[2]]
#'
#' @importFrom dplyr mutate arrange filter rename select group_by summarise ungroup full_join join_by
#' @importFrom tidyr pivot_wider pivot_longer
#' @importFrom magrittr %>%
#' @importFrom tidyselect any_of
#' @importFrom stats reorder sd var
#' @importFrom ctrlGene geNorm
#' @importFrom ggplot2 ggplot geom_line guides scale_color_manual xlab ylab labs theme_classic theme element_text scale_x_discrete facet_wrap ggsave
#' @import ggplot2
#' @import tidyverse
#' @importFrom utils write.table
#'
find_ref_gene <- function(data,
                          groups,
                          candidates,
                          colors,
                          norm.finder.score = TRUE,
                          genorm.score = TRUE,
                          line.width = 1,
                          angle = 0,
                          x.axis.title = "",
                          y.axis.title = "Ct",
                          axis.title.size = 11,
                          axis.text.size = 10,
                          legend.title = "",
                          legend.title.size = 11,
                          legend.text.size = 11,
                          legend.position = "top",
                          plot.title = "",
                          plot.title.size = 14,
                          save.to.tiff = FALSE,
                          dpi = 600,
                          width = 15,
                          height = 15,
                          name.tiff = "Ct_reference_gene_selection",
                          save.to.txt = FALSE,
                          name.txt = "reference_gene_selection") {
  if (groups[1] == "all") {
    data <- data

  } else {
    data <- filter(data, Group %in% groups)
  }

  if (sum(duplicated(data$Sample)) != 0) {
    data$Sample <- paste(data$Sample, data$Group, sep = "_")
  }

  ref <- data %>%
    pivot_longer(
      cols = -c(Group, Sample),
      names_to = "Gene",
      values_to = "Ct"
    ) %>%
    filter(Gene %in% candidates)

  ref_plot <-
    ggplot(ref, aes(
      x = Sample,
      y = Ct,
      color = Gene,
      group = Gene
    )) +
    geom_line(linewidth = line.width) +
    scale_color_manual(values = c(colors)) +
    guides(x =  guide_axis(angle = angle)) +
    theme_bw() +
    xlab(x.axis.title) +
    ylab(y.axis.title) +
    labs(color = legend.title, title = plot.title) +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, color = 'black')) +
    theme(axis.title = element_text(size = axis.title.size, color = 'black')) +
    theme(legend.title = element_text(size = legend.title.size, colour =
                                        "black")) +
    theme(legend.text = element_text(size = legend.text.size, colour = "black")) +
    theme(plot.title = element_text(size = plot.title.size))

  print(ref_plot)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      ref_plot,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }

  ref_var <- ref %>%
    group_by(Gene) %>%
    summarise(
      min = min(Ct),
      max = max(Ct),
      sd = sd(Ct, na.rm = TRUE),
      var = var(Ct, na.rm = TRUE),
      .groups = "keep"
    ) %>%
    as.data.frame()

  if (norm.finder.score == TRUE) {
    reference.stability.nF <- norm_finder(data, candidates = candidates)
    colnames(reference.stability.nF) <- "NormFinder_score"
    reference.stability.nF$Gene <- rownames(reference.stability.nF)
    ref_var <- ref_var %>%
      full_join(reference.stability.nF, by = join_by(Gene))
  }
  if (genorm.score == TRUE) {
    data <- data %>%
      ungroup() %>%
      select(-Group) %>%
      select(any_of(c("Sample", candidates))) %>%
      as.data.frame()

    rownames(data) <- data[, "Sample"]
    data <- select(data,-Sample)
    data <- as.matrix(data)

    reference.stability.gF <- geNorm(
      data,
      genes = data.frame(Gene = character(0), geNorm_score = numeric(0)),
      ctVal = TRUE
    )
    colnames(reference.stability.gF) <- c("Gene", "geNorm_score")
    ref_var <- ref_var %>%
      full_join(reference.stability.gF, by = join_by(Gene))

  }

if (save.to.txt == TRUE) {
    write.table(as.data.frame(ref_var), paste(name.txt, ".txt", sep = ""))
  }

  return(list(ref_plot, as.data.frame(ref_var)))
}







#' @title delta_Ct
#'
#' @description
#' This function calculates delta Ct (dCt) values by subtracting Ct values of reference
#' gene or genes from Ct values of remaining genes. Obtained dCt values can be further
#' transformed by using 2^-dCt formula (if transform == TRUE).
#'
#' @param data Data object returned from make_Ct_ready function.
#' @param ref Character vector with name of one or more reference genes.
#' @param normalise Logical: if TRUE, data normalisation will be done using reference gene/genes
#' (provided in ref parameter).
#' @param transform Logical: if TRUE, calculated dCt values will be transformed using 2^-dCt formula. Default to FALSE.
#' @param save.to.txt Logical: if TRUE, returned data will be saved to .txt file. Default to FALSE.
#' @param name.txt Character: name of saved .txt file, without ".txt" name of extension. Default to "data_dCt".
#'
#' @return Data frame with dCt values.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' head(data.dCt)
#' data.dCt.exp <- delta_Ct(data.CtF.ready, ref = "Gene8", transform = TRUE)
#' head(data.dCt.exp)
#'
#' @importFrom utils write.table
#' @importFrom tidyselect all_of
#' @importFrom dplyr mutate_at select vars
#' @import tidyverse
#'
delta_Ct <- function(data,
                     ref,
                     normalise = TRUE,
                     transform = FALSE,
                     save.to.txt = FALSE,
                     name.txt = "data_dCt") {
  if (normalise == TRUE) {

  if (length(ref) > 1) {
    data$calibrator <- rowMeans(data[, colnames(data) %in% ref])
    dCt <-  mutate_at(data,
                      vars(-c(
                        "Group", "Sample", all_of(ref), calibrator
                      )),
                      list(dCt = ~ . - calibrator))
  } else {
    dCt <- mutate_at(data,
                     vars(-c("Group", "Sample", all_of(ref))),
                     list(dCt = ~ . - .data[[ref]]))
  }

  dCt <- select(dCt, Group, Sample, ends_with("dCt"))
  colnames(dCt) <- sub("_dCt*", "", colnames(dCt))
  }

  if (transform == TRUE) {
  exp <- function(x) {
    x <- 2 ^ -x
  }

  dCt <- dCt %>%
    mutate_at(vars(-c("Group", "Sample")), exp)
  }

  if (save.to.txt == TRUE) {
    write.table(as.data.frame(dCt), paste(name.txt, ".txt", sep = ""))
  }

  return(dCt)
}









#' @title control_boxplot_sample
#'
#' @description
#' Boxplot that illustrate distribution of data in each sample. This function is helpful to identify outlier samples.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions. Also, a data frame with
#' fold change values returned from RQ_dCt() and RQ_ddCt() functions run in a pairwise mode can be used.
#' @param sel.Sample Character vector with names of samples to include, or "all" (default) to use all samples.
#' @param pairwise.FCh Logical: If fold change values returned from RQ_dCt() and RQ_ddCt() functions
#' in a pairwise approach are used as data, this parameter should be set to TRUE, otherwise to FALSE (default).
#' @param coef Numeric: how many times of interquartile range should be used to determine range point for whiskers. Default to 1.5.
#' @param colors Character vector containing colors for compared groups. Numbers of colors must be equal to number of groups. Default to c("#66c2a5", "#fc8d62").
#' If pairwise.FCh parameter is set to TRUE, one color is required (data contain no groups).
#' @param x.axis.title Character: title of x axis. Default to "Sample".
#' @param y.axis.title Character: title of y axis. Default to "value".
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 12.
#' @param legend.title Character: title of legend. Default to "Group".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_boxplot_samples".
#'
#' @return Object with a boxplot illustrating distribution of data in each sample. Created plot is also displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' control.boxplot.sample <- control_boxplot_sample(data.dCt)
#'
#' @importFrom dplyr filter select
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_boxplot coord_flip scale_color_manual xlab ylab labs theme_classic theme element_text ggsave
#' @import ggplot2
#' @import tidyverse
#'
control_boxplot_sample <- function(data,
                                   sel.Sample = "all",
                                   pairwise.FCh = FALSE,
                                   coef = 1.5,
                                   colors = c("#66c2a5", "#fc8d62"),
                                   x.axis.title = "Sample",
                                   y.axis.title = "value",
                                   axis.title.size = 11,
                                   axis.text.size = 12,
                                   legend.title = "Group",
                                   legend.title.size = 11,
                                   legend.text.size = 11,
                                   legend.position = "right",
                                   plot.title = "",
                                   plot.title.size = 14,
                                   save.to.tiff = FALSE,
                                   dpi = 600,
                                   width = 15,
                                   height = 15,
                                   name.tiff = "control_boxplot_samples") {
  if (sel.Sample[1] == "all") {
    data <- data

  } else {
    data <- filter(data, Sample %in% sel.Sample)
  }

  if (pairwise.FCh == TRUE) {
    data <- data %>%
      select(Sample, Gene, FCh)

    box_control_sample <- ggplot(data, aes(x = Sample, y = FCh)) +
      geom_boxplot(coef = 1.5, fill = colors[1]) +
      scale_x_discrete(limits = rev(unique(data$Sample))) +
      labs(title = plot.title)

  } else {
    data <-
      pivot_longer(data,
                   !c(Sample, Group),
                   names_to = "Gene" ,
                   values_to = "value")

    box_control_sample <-
      ggplot(data, aes(x = Sample, y = value, color = Group)) +
      geom_boxplot(coef = coef) +
      scale_x_discrete(limits = rev(unique(data$Sample))) +
      scale_color_manual(values = c(colors)) +
      labs(color = legend.title, title = plot.title)
  }

  box_control_sample <- box_control_sample +
    coord_flip() +
    xlab(x.axis.title) + ylab(y.axis.title) +
    theme_classic() +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
    theme(legend.title = element_text(size = legend.title.size, colour =
                                        "black")) +
    theme(legend.text = element_text(size = legend.text.size, colour =
                                       "black")) +
    theme(plot.title = element_text(size = plot.title.size))

  print(box_control_sample)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      box_control_sample,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }
  return(box_control_sample)
}






#' @title control_boxplot_gene
#'
#' @description
#' This function creates boxplot that illustrate distribution of data in each gene.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions. Also, a data frame with
#' fold change values returned from RQ_dCt() and RQ_ddCt() functions run in a pairwise mode can be used.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all genes.
#' @param pairwise.FCh Logical: If fold change values returned from RQ_dCt() and RQ_ddCt() functions
#' in a pairwise approach are used as data, this parameter should be set to TRUE, otherwise to FALSE (default).
#' @param by.group Logical: if TRUE, distributions will be drawn by compared groups of samples.
#' @param coef Numeric: how many times of interquartile range should be used to determine range point for whiskers. Default to 1.5.
#' @param colors Character vector containing colors for groups, length of one (when by.group = FALSE) or equal to the number of groups (when by.group = TRUE).
#' @param x.axis.title Character: title of x axis. Default to "Gene".
#' @param y.axis.title Character: title of y axis. Default to "value".
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 12.
#' @param legend.title Character: title of legend. Default to "Group".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_boxplot_genes".
#'
#' @return Object with boxplot illustrating distribution of data for each gene. Created plot is also displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' control.boxplot.gene <- control_boxplot_gene(data.dCt)
#'
#' @importFrom dplyr filter select
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_boxplot coord_flip scale_color_manual xlab ylab labs theme_classic theme element_text ggsave scale_x_discrete
#' @import ggplot2
#' @import tidyverse
#'
control_boxplot_gene <- function(data,
                                 sel.Gene = "all",
                                 pairwise.FCh = FALSE,
                                 coef = 1.5,
                                 by.group = TRUE,
                                 colors = c("#66c2a5", "#fc8d62"),
                                 axis.title.size = 11,
                                 axis.text.size = 12,
                                 x.axis.title = "Gene",
                                 y.axis.title = "value",
                                 legend.title = "Group",
                                 legend.title.size = 11,
                                 legend.text.size = 11,
                                 legend.position = "right",
                                 plot.title = "",
                                 plot.title.size = 14,
                                 save.to.tiff = FALSE,
                                 dpi = 600,
                                 width = 15,
                                 height = 15,
                                 name.tiff = "control_boxplot_genes") {
  if (sel.Gene[1] == "all") {
    data <- data

  } else {
    if (pairwise.FCh == TRUE) {
      data <- filter(data, Gene %in% sel.Gene)

    } else {
      data <- select(data, Group, Sample, any_of(sel.Gene))
    }
  }



  if (pairwise.FCh == TRUE) {
    data <- data %>%
      select(Sample, Gene, FCh)

  } else {
    data <-
      pivot_longer(data,
                   !c(Sample, Group),
                   names_to = "Gene" ,
                   values_to = "value")
  }

  if (pairwise.FCh == TRUE) {
    box_control_genes <- ggplot(data, aes(x = Gene, y = FCh)) +
      geom_boxplot(coef = coef, fill = colors[1]) +
      labs(title = plot.title)

  } else {
    if (by.group == TRUE) {
      box_control_genes <-
        ggplot(data, aes(x = Gene, y = value, color = Group)) +
        geom_boxplot(coef = coef) +
        scale_color_manual(values = c(colors)) +
        labs(color = legend.title, title = plot.title)

    } else {
      box_control_genes <- ggplot(data, aes(x = Gene, y = value)) +
        geom_boxplot(coef = coef, fill = colors[1]) +
        labs(title = plot.title)
    }
  }

  box_control_genes <- box_control_genes +
    scale_x_discrete(limits = rev(unique(data$Gene))) +
    coord_flip() +
    xlab(x.axis.title) +
    ylab(y.axis.title) +
    theme_classic() +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
    theme(legend.title = element_text(size = legend.title.size, colour =
                                        "black")) +
    theme(legend.text = element_text(size = legend.text.size, colour = "black")) +
    theme(plot.title = element_text(size = plot.title.size))


  print(box_control_genes)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      box_control_genes,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }
  return(box_control_genes)
}






#' @title control_cluster_sample
#'
#' @description
#' This function performs hierarchical clustering of samples based on the data, useful to identify outlier samples.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' Also, table with fold change values returned from RQ_dCt() and RQ_ddCt() functions can be used,
#' but in such situation a pairwise.FCh parameter must be set to TRUE.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all genes.
#' @param method.dist Character: name of method used for calculation of distances, derived from stats::dist() function,
#' must be one of "euclidean" (default) , "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param method.clust Character: name of used method for agglomeration, derived from stats::hclust() function,
#' must be one of "ward.D", "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".
#' @param pairwise.FCh Logical: If fold change values returned from RQ_dCt() and RQ_ddCt() functions
#' in a pairwise approach are used as data, this parameter should be set to TRUE, otherwise to FALSE (default).
#' @param x.axis.title Character: title of x axis. Default to "Samples".
#' @param y.axis.title Character: title of y axis. Default to "Height".
#' @param label.size Numeric: size of text labels. Default to 0.8.
#' @param plot.title Character: title of plot. Default to "".
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_clust_samples".
#'
#' @return Plot with hierarchical clustering of samples, displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                        remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' control_cluster_sample(data.dCt)
#'
#' @importFrom stats hclust dist
#' @importFrom grDevices tiff
#' @importFrom dplyr select ungroup
#' @importFrom tidyr pivot_wider
#' @import tidyverse
#'
control_cluster_sample <- function(data,
                                   sel.Gene = "all",
                                   method.dist = "euclidean",
                                   method.clust = "average",
                                   pairwise.FCh = FALSE,
                                   x.axis.title = "Samples",
                                   y.axis.title = "Height",
                                   label.size = 0.8,
                                   plot.title = "",
                                   save.to.tiff = FALSE,
                                   dpi = 600,
                                   width = 15,
                                   height = 15,
                                   name.tiff = "control_clust_samples") {
  if (sel.Gene[1] == "all") {
    data <- data

  } else {
    if (pairwise.FCh == TRUE) {
      data <- filter(data, Gene %in% sel.Gene)

    } else {
      data <- select(data, Group, Sample, any_of(sel.Gene))
    }
  }

  if (pairwise.FCh == TRUE) {
    data <- data %>%
      select(Sample, Gene, FCh) %>%
      pivot_wider(names_from = Gene, values_from = FCh) %>%
      ungroup()

    cluster <- hclust(dist(select(data,-Sample),
                           method = method.dist),
                      method = method.clust)
  } else{
    data <- ungroup(data)
    cluster <- hclust(dist(select(data,-Group,-Sample),
                           method = method.dist),
                      method = method.clust)
  }

  cluster$labels <- data$Sample
  plot(
    cluster,
    xlab = x.axis.title,
    ylab = y.axis.title,
    main = plot.title,
    cex = label.size
  )

  if (save.to.tiff == TRUE) {
    tiff(
      paste(name.tiff, ".tiff", sep = ""),
      res = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )

    plot(
      cluster,
      xlab  = x.axis.title,
      ylab = y.axis.title,
      main = plot.title,
      cex = label.size
    )
    dev.off()
  }
}






#' @title control_cluster_gene
#'
#' @description
#' This function performs hierarchical clustering of genes based on the data, useful to gain insight into similarity in expression of analyzed genes.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' Also, table with fold change values returned from RQ_dCt() and RQ_ddCt() functions can be used,
#' but in such situation a pairwise.FCh parameter must be set to TRUE.
#' @param sel.Sample Character vector with names of samples to include, or "all" (default) to use all samples.
#' @param method.dist Character: name of method used for calculation of distances, derived from stats::dist() function,
#' must be one of "euclidean" (default) , "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param method.clust Character: name of used method for agglomeration, derived from stats::hclust() function,
#' must be one of "ward.D", "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".
#' @param pairwise.FCh Logical: If fold change values returned from RQ_dCt() and RQ_ddCt() functions
#' in a pairwise approach are used as data, this parameter should be set to TRUE, otherwise to FALSE (default).
#' @param x.axis.title Character: title of x axis. Default to "Genes".
#' @param y.axis.title Character: title of y axis. Default to "Height".
#' @param plot.title Character: title of plot. Default to "".
#' @param label.size Numeric: size of text labels. Default to 0.8.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_clust_genes".
#'
#' @return Plot with hierarchical clustering of genes, displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' control_cluster_gene(data.dCt)
#'
#' @importFrom stats hclust dist
#' @importFrom grDevices tiff
#' @importFrom dplyr select ungroup filter
#' @importFrom tidyr pivot_wider
#' @import tidyverse
#'
control_cluster_gene <- function (data,
                                  method.dist = "euclidean",
                                  sel.Sample = "all",
                                  method.clust = "average",
                                  pairwise.FCh = FALSE,
                                  x.axis.title = "Genes",
                                  y.axis.title = "Height",
                                  label.size = 0.8,
                                  plot.title = "",
                                  save.to.tiff = FALSE,
                                  dpi = 600,
                                  width = 15,
                                  height = 15,
                                  name.tiff = "control_clust_genes")
{
  if (sel.Sample[1] == "all") {
    data <- data

  } else {
    data <- filter(data, Sample %in% sel.Sample)
  }

  if (pairwise.FCh == TRUE) {
    data_t <- data %>%
      select(Sample, Gene, FCh) %>%
      pivot_wider(names_from = Sample, values_from = FCh) %>%
      ungroup()

    cluster <- hclust(dist(select(data_t, -Gene),
                           method = method.dist),
                      method = method.clust)
    cluster$labels <- data_t$Gene

  } else{
    data_t <- ungroup(data)
    data_t <- t(select(data_t, -Group, -Sample))
    colnames(data_t) <- data$Sample

    cluster <- hclust(dist(as.data.frame(data_t),
                           method = method.dist),
                      method = method.clust)
  }

  plot(
    cluster,
    xlab = x.axis.title,
    ylab = y.axis.title,
    main = plot.title,
    cex = label.size
  )

  if (save.to.tiff == TRUE) {
    tiff(
      paste(name.tiff, ".tiff", sep = ""),
      res = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
    plot(
      cluster,
      xlab = x.axis.title,
      ylab = y.axis.title,
      main = plot.title,
      cex = label.size
    )
    dev.off()
  }
}



#' @title control_pca_sample
#'
#' @description
#' This function performs principal component analysis (PCA) for samples and generate plot that illustrate spatial arrangement
#' of samples based on the two first components. This plot is useful to identify outlier samples.
#' PCA analysis can not deal with missing values, thus all samples with at least one missing value are removed from data before analysis.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' Also, table with fold change values returned from RQ_dCt() and RQ_ddCt() functions can be used,
#' but in such situation a pairwise.FCh parameter must be set to TRUE.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all genes.
#' @param pairwise.FCh Logical: If fold change values returned from RQ_dCt() and RQ_ddCt() functions
#' in a pairwise approach are used as data, this parameter should be set to TRUE, otherwise to FALSE (default).
#' @param point.size Numeric: size of points. Default to 4.
#' @param point.shape Integer: shape of points. Default to 19.
#' @param alpha Numeric: transparency of points, a value between 0 and 1. Default to 0.7.
#' @param label.size Numeric: size of points labels (names of samples). Default to 3.
#' @param hjust Numeric: horizontal position of points labels. Default to 0.
#' @param vjust Numeric: vertical position of points labels.  Default to -1.
#' @param colors Character vector containing colors for compared groups. The number of colors must be equal to the number of groups.
#' Default to c("#66c2a5", "#fc8d62").
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param legend.title Character: title of legend. Default to "Group".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_pca_samples".
#'
#' @return Object with plot. The plot is also displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' control_pca_sample(data.dCt)
#'
#' @importFrom stats na.omit prcomp
#' @importFrom dplyr select filter
#' @importFrom tidyr pivot_wider
#' @importFrom ggplot2 ggplot geom_point geom_text scale_color_manual xlab ylab labs theme_classic theme element_text ggsave
#' @import ggplot2
#' @import tidyverse
#'
control_pca_sample <- function(data,
                               sel.Gene = "all",
                               pairwise.FCh = FALSE,
                               point.size = 4,
                               point.shape = 19,
                               alpha = 0.7,
                               colors = c("#66c2a5", "#fc8d62"),
                               label.size = 3,
                               hjust = 0.5,
                               vjust = -1,
                               axis.title.size = 11,
                               axis.text.size = 10,
                               legend.text.size = 11,
                               legend.title = "Group",
                               legend.title.size = 11,
                               legend.position = "right",
                               plot.title = "",
                               plot.title.size = 14,
                               save.to.tiff = FALSE,
                               dpi = 600,
                               width = 15,
                               height = 15,
                               name.tiff = "control_pca_samples") {
  if (sel.Gene[1] == "all") {
    data <- data

  } else {
    if (pairwise.FCh == TRUE) {
      data <- filter(data, Gene %in% sel.Gene)

    } else {
      data <- select(data, Group, Sample, any_of(sel.Gene))
    }
  }


  if (pairwise.FCh == TRUE) {
    data <- data %>%
      select(Sample, Gene, FCh) %>%
      pivot_wider(names_from = Gene, values_from = FCh) %>%
      as.data.frame()
    rownames(data) <- data$Sample
    data <- na.omit(data)
    pca <- prcomp(select(data,-Sample), scale = TRUE)
    var_pca1 <- summary(pca)$importance[2, ][1]
    var_pca2 <- summary(pca)$importance[2, ][2]
    pca_comp <- as.data.frame(pca$x)
    pca_comp$Sample <- data$Sample

    control_pca <-
      ggplot(pca_comp, aes(x = PC1, y = PC2, label = Sample)) +
      geom_point(
        size = point.size,
        shape = point.shape,
        alpha = alpha,
        col = colors[1]
      ) +
      labs(title = plot.title) +
      theme_bw() +
      labs(
        x = paste("PC1: ", round(var_pca1 * 100, 2), "% variance explained", sep = ""),
        y = paste("PC2: ", round(var_pca2 * 100, 2), "% variance explained", sep = "")
      ) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
      theme(plot.title = element_text(size = plot.title.size)) +
      geom_text(
        aes(label = Sample),
        hjust = hjust,
        vjust = vjust,
        size = label.size
      )

  } else {
    data <- as.data.frame(data)
    data <- unite(data, Rownames, Sample, Group, remove = FALSE)
    rownames(data) <- data$Rownames
    data <- na.omit(data)
    pca <-
      prcomp(select(data,-Rownames,-Sample,-Group), scale = TRUE)
    var_pca1 <- summary(pca)$importance[2, ][1]
    var_pca2 <- summary(pca)$importance[2, ][2]
    pca_comp <- as.data.frame(pca$x)
    pca_comp$Sample <- data$Sample
    pca_comp$Group <- data$Group

    control_pca <-
      ggplot(pca_comp, aes(
        x = PC1,
        y = PC2,
        label = Sample,
        color = Group
      )) +
      geom_point(size = point.size,
                 shape = point.shape,
                 alpha = alpha) +
      scale_color_manual(values = c(colors)) +
      labs(colour = legend.title, title = plot.title) +
      theme_bw() +
      labs(
        x = paste("PC1: ", round(var_pca1 * 100, 2), "% variance explained", sep = ""),
        y = paste("PC2: ", round(var_pca2 * 100, 2), "% variance explained", sep = "")
      ) +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
      theme(legend.text = element_text(size = legend.text.size, colour = "black")) +
      theme(legend.title = element_text(size = legend.title.size, colour =
                                          "black")) +
      theme(plot.title = element_text(size = plot.title.size)) +
      geom_text(
        aes(label = Sample),
        hjust = hjust,
        vjust = vjust,
        size = label.size
      )
  }

  print(control_pca)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      control_pca,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }
  return(control_pca)
}








#' @title control_pca_gene
#'
#' @description
#' This function performs principal component analysis (PCA) for genes and generates plot illustrating spatial arrangement of genes using two first PCA components.
#' This plot allows to gain insight into similarity in expression of analyzed genes. PCA analysis can not deal with missing values,
#' thus all genes with at least one missing value are removed from data before analysis.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' Also, table with fold change values returned from RQ_dCt() and RQ_ddCt() functions can be used,
#' but in such situation a pairwise.FCh parameter must be set to TRUE.
#' @param sel.Sample Character vector with names of samples to include, or "all" (default) to use all samples.
#' @param pairwise.FCh Logical: If fold change values returned from RQ_dCt() and RQ_ddCt() functions
#' in a pairwise approach are used as data, this parameter should be set to TRUE, otherwise to FALSE (default).
#' @param point.size Numeric: size of points. Default to 4.
#' @param point.shape Integer: shape of points. Default to 19.
#' @param alpha Numeric: transparency of points, a value between 0 and 1. Default to 0.7.
#' @param label.size Numeric: size of points labels (names of samples). Default to 3.
#' @param hjust Numeric: horizontal position of points labels. Default to 0.
#' @param vjust Numeric: vertical position of points labels. Default to -1.
#' @param color Character: color used for points.
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param legend.title Character: title of legend. Default to "Group".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_pca_genes".
#'
#' @return Object with plot. The plot is also displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' control_pca_gene(data.dCt)
#'
#' @importFrom stats na.omit prcomp
#' @importFrom dplyr select filter
#' @importFrom ggplot2 ggplot geom_point geom_text xlab ylab labs theme_classic theme element_text ggsave
#' @import ggplot2
#' @import tidyverse
#'
control_pca_gene <- function(data,
                             sel.Sample = "all",
                             pairwise.FCh = FALSE,
                             point.size = 4,
                             point.shape = 19,
                             alpha = 0.7,
                             label.size = 3,
                             hjust = 0.5,
                             vjust = -1,
                             color = "black",
                             axis.title.size = 11,
                             axis.text.size = 10,
                             legend.text.size = 11,
                             legend.title = "Group",
                             legend.title.size = 11,
                             legend.position = "right",
                             plot.title = "",
                             plot.title.size = 14,
                             save.to.tiff = FALSE,
                             dpi = 600,
                             width = 15,
                             height = 15,
                             name.tiff = "control_pca_genes") {
  if (sel.Sample[1] == "all") {
    data <- data

  } else {
    data <- filter(data, Sample %in% sel.Sample)
  }

  if (pairwise.FCh == TRUE) {
    data <- data %>%
      select(Sample, Gene, FCh) %>%
      pivot_wider(names_from = Sample, values_from = FCh) %>%
      as.data.frame()
    rownames(data) <- data$Gene
    data <- na.omit(data)
    pca <- prcomp(select(data,-Gene), scale = TRUE)
    var_pca1 <- summary(pca)$importance[2, ][1]
    var_pca2 <- summary(pca)$importance[2, ][2]
    pca_comp <- as.data.frame(pca$x)
    pca_comp$Gene <- data$Gene

    control_pca <-
      ggplot(pca_comp, aes(x = PC1, y = PC2, label = Gene)) +
      geom_point(
        size = point.size,
        shape = point.shape,
        alpha = alpha,
        col = color
      ) +
      labs(title = plot.title) +
      theme_bw() +
      labs(
        x = paste("PC1: ", round(var_pca1 * 100, 2), "% variance explained", sep = ""),
        y = paste("PC2: ", round(var_pca2 * 100, 2), "% variance explained", sep = "")
      ) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
      theme(plot.title = element_text(size = plot.title.size)) +
      geom_text(
        aes(label = Gene),
        hjust = hjust,
        vjust = vjust,
        size = label.size
      )

  } else {
    data <- ungroup(data)
    data <- as.data.frame(data)
    data_t <- t(select(data,-Group,-Sample))
    colnames(data_t) <- data$Sample
    data_t <- na.omit(data_t)
    pca <- prcomp(data_t, scale = TRUE)
    var_pca1 <- summary(pca)$importance[2, ][1]
    var_pca2 <- summary(pca)$importance[2, ][2]
    pca_comp <- as.data.frame(pca$x)
    pca_comp$Gene <- rownames(pca_comp)

    control_pca <-
      ggplot(pca_comp, aes(x = PC1, y = PC2, label = Gene)) +
      geom_point(
        size = point.size,
        shape = point.shape,
        alpha = alpha,
        col = color
      ) +
      labs(title = plot.title) +
      theme_bw() +
      labs(
        x = paste("PC1: ", round(var_pca1 * 100, 2), "% variance explained", sep = ""),
        y = paste("PC2: ", round(var_pca2 * 100, 2), "% variance explained", sep = "")
      ) +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
      theme(plot.title = element_text(size = plot.title.size)) +
      geom_text(
        aes(label = Gene),
        hjust = hjust,
        vjust = vjust,
        size = label.size
      )
  }

  print(control_pca)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      control_pca,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }
  return(control_pca)
}





#' @title corr_gene
#'
#' @description
#' This function performs correlation analysis of genes based on the data, useful to gain insight into relationships between analyzed genes.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all genes.
#' @param type Character: type of displayed matrix, must be one of the 'full' (full matrix), 'upper' (upper triangular, default) or 'lower' (lower triangular).
#' @param add.coef If correlation coefficients should be add to the plot, specify color of coefficients (default to "black").
#' If NULL, correlation coefficients will not be printed.
#' @param method Character: type of correlations to compute, can be "pearson" (default) for Pearson's correlation coefficients
#' or "spearman" for Spearman's rank correlation coefficients.
#' @param order Character: method used for ordering the correlation matrix, inherited from corrplot::corrplot() function.
#' Must be one of the "original" (original order), "AOE" (angular order of the eigenvectors), "FPC" (first principal component order),
#' "hclust" (hierarchical clustering order, default), or "alphabet" (alphabetical order).
#' @param hclust.method Character: name of method used for hclust agglomeration, must be one of "ward", ward.D",
#' "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".
#' @param size Numeric: size of variable names and numbers in legend. Default to 0.6.
#' @param coef.size Numeric: size of correlation coefficients. Default to 0.6.
#' @param p.adjust.method Character: p value correction method for multiple testing, one of the "holm", "hochberg", "hommel",
#' "bonferroni", "BH" (default), "BY","fdr", or "none". See documentation for stats::p.adjust() function for details.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "corr_genes".
#' @param save.to.txt Logical: if TRUE, correlation results (sorted by absolute values of correlation coefficients in descending order)
#' will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension.. Default to "corr_genes".
#'
#' @return Plot illustrating correlation matrix (displayed on the graphic device) and
#' data frame with computed correlation coefficients and p values.
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(Hmisc)
#' library(corrplot)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' corr.genes <- corr_gene(data.dCt)
#' head(corr.genes)
#'
#' @importFrom stats p.adjust
#' @importFrom utils write.table
#' @importFrom dplyr select arrange
#' @importFrom Hmisc rcorr
#' @importFrom corrplot corrplot
#' @importFrom grDevices tiff
#' @import corrplot
#' @import tidyverse
#'
corr_gene <- function(data,
                      sel.Gene = "all",
                      type = "upper",
                      method = "pearson",
                      add.coef = "black",
                      order = "hclust",
                      hclust.method = "average",
                      size = 0.6,
                      coef.size = 0.6,
                      p.adjust.method = "BH",
                      save.to.tiff = FALSE,
                      dpi = 600,
                      width = 15,
                      height = 15,
                      name.tiff = "corr_genes",
                      save.to.txt = FALSE,
                      name.txt = "corr_genes") {
  if (sel.Gene[1] == "all") {
    data <- data

  } else {
    data <- select(data, Group, Sample, any_of(sel.Gene))
  }

  data <- as.data.frame(data)
  data <- select(data,-Group,-Sample)
  res_cor <- rcorr(as.matrix(data), type = method)

  if (order == "hclust") {
    corrplot(
      res_cor$r,
      type = type,
      addCoef.col = add.coef,
      tl.cex = size,
      cl.cex = size,
      tl.col = "black",
      number.cex = coef.size,
      order = order,
      hclust.method = hclust.method
    )

    if (save.to.tiff == TRUE) {
      tiff(
        paste(name.tiff, ".tiff", sep = ""),
        res = dpi,
        width = width,
        height = height,
        units = "cm",
        compression = "lzw"
      )
      corrplot(
        res_cor$r,
        type = type,
        tl.cex = size,
        tl.col = "black",
        cl.cex = size,
        order = order,
        hclust.method = hclust.method,
        addCoef.col = add.coef,
        number.cex = coef.size
      )
      dev.off()
    }
  } else {
    corrplot(
      res_cor$r,
      type = type,
      tl.cex = size,
      tl.col = "black",
      cl.cex = size,
      order = order,
      addCoef.col = add.coef,
      number.cex = coef.size
    )

    if (save.to.tiff == TRUE) {
      tiff(
        paste(name.tiff, ".tiff", sep = ""),
        res = dpi,
        width = width,
        height = height,
        units = "cm",
        compression = "lzw"
      )
      corrplot(
        res_cor$r,
        type = type,
        tl.cex = size,
        tl.col = "black",
        cl.cex = size,
        order = order,
        addCoef.col = "black",
        number.cex = coef.size
      )
      dev.off()
    }
  }
  corr <- upper.tri(res_cor$r)
  corr.data <- data.frame(
    row = rownames(res_cor$r)[row(res_cor$r)[corr]],
    column = rownames(res_cor$r)[col(res_cor$r)[corr]],
    cor  = (res_cor$r)[corr],
    p = res_cor$P[corr]
  )
  corr.data$p.adj <-
    p.adjust(corr.data$p, method = p.adjust.method)
  corr.data.sort <- arrange(corr.data,-abs(cor))

  if (save.to.txt == TRUE) {
    write.table(as.data.frame(corr.data.sort), paste(name.txt, ".txt", sep = ""))
  }
  return(as.data.frame(corr.data.sort))
}





#' @title corr_sample
#'
#' @description
#' This function performs correlation analysis of samples based on the data. Results are useful to gain insight into relationships between analyzed samples.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' @param sel.Sample Character vector with names of samples to include, or "all" (default) to use all samples.
#' @param type Character: type of displayed matrix, must be one of the 'full' (full matrix), 'upper' (upper triangular, default) or 'lower' (lower triangular).
#' @param add.coef If correlation coefficients should be add to the plot, specify color of coefficients (default to "black").
#' If NULL, correlation coefficients will not be printed.
#' @param method Character: type of correlations to compute, can be "pearson" (default) for Pearson's correlation coefficients
#' or "spearman" for Spearman's rank correlation coefficients.
#' @param order Character: method used for ordering the correlation matrix, inherited from corrplot::corrplot() function.
#' Must be one of the "original" (original order), "AOE" (angular order of the eigenvectors), "FPC" (first principal component order),
#' "hclust" (hierarchical clustering order, default), or "alphabet" (alphabetical order).
#' @param hclust.method Character: name of method used for hclust agglomeration, must be one of "ward", ward.D",
#' "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".
#' @param size Numeric: size of variable names and numbers in legend. Default to 0.6.
#' @param coef.size Numeric: size of correlation coefficients. Default to 0.6.
#' @param p.adjust.method Character: p value correction method for multiple testing, one of the "holm", "hochberg", "hommel",
#' "bonferroni", "BH" (default), "BY","fdr", or "none". See documentation for stats::p.adjust() function for details.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "corr_samples".
#' @param save.to.txt Logical: if TRUE, correlation results (sorted by absolute values of correlation coefficients in descending order)
#' will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension.. Default to "corr_samples".
#'
#' @return Plot illustrating correlation matrix (displayed on the graphic device) and data.frame with computed correlation coefficients and p values.
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(Hmisc)
#' library(corrplot)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' corr.samples <- corr_sample(data.CtF.ready)
#' head(corr.samples)
#'
#' @importFrom stats p.adjust
#' @importFrom utils write.table
#' @importFrom dplyr select arrange filter
#' @importFrom Hmisc rcorr
#' @importFrom corrplot corrplot
#' @importFrom grDevices tiff
#' @import corrplot
#' @import tidyverse
#'
corr_sample <- function(data,
                        sel.Sample = "all",
                        type = "upper",
                        method = "pearson",
                        add.coef = "black",
                        order = "hclust",
                        hclust.method = "average",
                        size = 0.6,
                        coef.size = 0.6,
                        p.adjust.method = "BH",
                        save.to.tiff = FALSE,
                        dpi = 600,
                        width = 15,
                        height = 15,
                        name.tiff = "corr_samples",
                        save.to.txt = FALSE,
                        name.txt = "corr_samples") {
  if (sel.Sample[1] == "all") {
    data <- data

  } else {
    data <- filter(data, Sample %in% sel.Sample)
  }

  data <- as.data.frame(data)
  data_t <- t(select(data,-Group,-Sample))
  colnames(data_t) <- data$Sample
  res_cor <- rcorr(as.matrix(data_t), type = method)

  if (order == "hclust") {
    corrplot(
      res_cor$r,
      type = type,
      addCoef.col = add.coef,
      tl.cex = size,
      cl.cex = size,
      tl.col = "black",
      number.cex = coef.size,
      order = order,
      hclust.method = hclust.method
    )

    if (save.to.tiff == TRUE) {
      tiff(
        paste(name.tiff, ".tiff", sep = ""),
        res = dpi,
        width = width,
        height = height,
        units = "cm",
        compression = "lzw"
      )
      corrplot(
        res_cor$r,
        type = type,
        tl.cex = size,
        tl.col = "black",
        cl.cex = size,
        order = order,
        hclust.method = hclust.method,
        addCoef.col = add.coef,
        number.cex = coef.size
      )
      dev.off()
    }
  } else {
    corrplot(
      res_cor$r,
      type = type,
      tl.cex = size,
      tl.col = "black",
      cl.cex = size,
      order = order,
      addCoef.col = add.coef,
      number.cex = coef.size
    )

    if (save.to.tiff == TRUE) {
      tiff(
        paste(name.tiff, ".tiff", sep = ""),
        res = dpi,
        width = width,
        height = height,
        units = "cm",
        compression = "lzw"
      )
      corrplot(
        res_cor$r,
        type = type,
        tl.cex = size,
        tl.col = "black",
        cl.cex = size,
        order = order,
        addCoef.col = "black",
        number.cex = coef.size
      )
      dev.off()
    }
  }

  corr <- upper.tri(res_cor$r)
  corr.data <- data.frame(
    row = rownames(res_cor$r)[row(res_cor$r)[corr]],
    column = rownames(res_cor$r)[col(res_cor$r)[corr]],
    cor  = (res_cor$r)[corr],
    p = res_cor$P[corr]
  )
  corr.data$p.adj <-
    p.adjust(corr.data$p, method = p.adjust.method)
  corr.data.sort <- arrange(corr.data,-abs(cor))

  if (save.to.txt == TRUE) {
    write.table(as.data.frame(corr.data.sort), paste(name.txt, ".txt", sep = ""))
  }

  return(as.data.frame(corr.data.sort))
}







#' @title single_pair_gene
#'
#' @description
#' This function generates scatter plot with linear regression line for two specified genes.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' @param x,y Characters: names of genes to use.
#' @param by.group Logical: if TRUE (default), relationships will be shown separately for compared groups.
#' @param point.size Numeric: size of points. Default to 4.
#' @param point.shape Integer: shape of points. Default to 19.
#' @param point.alpha Numeric: transparency of points, a value between 0 and 1. Default to 0.7.
#' @param colors Character vector containing colors for compared groups. The number of colors must be equal to the number of groups.
#' Default to c("#66c2a5", "#fc8d62").
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param legend.title Character: title of legend. Default to "Group".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param labels Logical: if TRUE (default), a regression statistics will be added to the plot using ggpimsc::stat_poly_eq() function.
#' @param label Character: specifies regression statistics to add, names specified by ggpimsc::stat_poly_eq() function must be used.
#' Default to c("eq", "R2", "p") for regression equation, coefficient of determination and p value, respectively.
#' @param label.position.x,label.position.y  Numeric: coordinates for position of regression statistics. If by.group = TRUE,
#' a vector length of groups number can be provided to avoid overlapping. See description of label.x and label.y parameters from ggpimsc::stat_poly_eq() function.
#' @param small.p,small.r Logical, if TRUE, p character in p value label and r character in coefficient of determination label will be lowercase. Default to FALSE.
#' @param rr.digits,p.digits Integer: number of digits after the decimal point in coefficient of determination and p value in labels.
#' Default to 2 for rr.digits and 3 for p.digits.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "Gene1_Gene2_single_pair_plot".
#'
#' @return Object with plot. Created plot is also displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(ggpmisc)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' single_pair_gene(data.dCt, "Gene16", "Gene17")
#'
#' @importFrom ggplot2 ggplot geom_point geom_smooth scale_color_manual xlab ylab labs theme_classic theme element_text ggsave
#' @importFrom ggpmisc stat_poly_eq use_label
#' @import ggplot2
#' @import ggpmisc
#'
single_pair_gene <- function(data,
                             x,
                             y,
                             by.group = TRUE,
                             point.size = 4,
                             point.shape = 19,
                             point.alpha = 0.7,
                             colors = c("#66c2a5", "#fc8d62"),
                             axis.title.size = 11,
                             axis.text.size = 10,
                             legend.title = "Group",
                             legend.title.size = 11,
                             legend.text.size = 11,
                             legend.position = "right",
                             plot.title = "",
                             plot.title.size = 14,
                             labels = TRUE,
                             label = c("eq", "R2", "p"),
                             label.position.x = c(1, 1),
                             label.position.y = c(1, 0.95),
                             small.p = FALSE,
                             small.r = FALSE,
                             p.digits = 3,
                             rr.digits = 2,
                             save.to.tiff = FALSE,
                             dpi = 600,
                             width = 15,
                             height = 15,
                             name.tiff = "_single_pair_plot") {
  if (by.group == TRUE) {
    single_pair <-
      ggplot(data, aes(x = .data[[x]], y = .data[[y]], color = Group)) +
      geom_point(size = point.size,
                 shape = point.shape,
                 alpha = point.alpha) +
      geom_smooth(method = 'lm', se = FALSE) +
      scale_color_manual(values = c(colors)) +
      xlab(x) +
      ylab(y) +
      labs(color = legend.title, title = plot.title) +
      theme_classic() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
      theme(legend.text = element_text(size = legend.text.size, colour =
                                         "black")) +
      theme(legend.title = element_text(size = legend.title.size, colour =
                                          "black")) +
      theme(plot.title = element_text(size = plot.title.size))

  } else {
    single_pair <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
      geom_point(
        size = point.size,
        shape = point.shape,
        alpha = point.alpha,
        col = colors[1]
      ) +
      geom_smooth(method = 'lm', se = FALSE) +
      xlab(x) + ylab(y) +
      labs(title = plot.title) +
      theme_classic() +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
      theme(plot.title = element_text(size = plot.title.size))
  }
  if (labels == TRUE) {
    single_pair <- single_pair +
      stat_poly_eq(
        use_label(label),
        label.y = c(label.position.y),
        label.x = c(label.position.x),
        small.p = small.p,
        small.r = small.r,
        p.digits = p.digits,
        rr.digits = rr.digits
      )
  }

  print(single_pair)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(x, "_", y, name.tiff, ".tiff", sep = ""),
      single_pair,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }
  return(single_pair)
}







#' @title single_pair_sample
#'
#' @description
#' This function generates scatter plot with linear regression line for two specified samples.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' @param pairwise.data Logical: set to TRUE if a pairwise data is used, default to FALSE.
#' @param by.group Logical: If TRUE (only for pairwise analysis), relationships will be shown separately for groups.
#' Default to FALSE.
#' @param x,y Characters: names of samples to use.
#' @param point.size Numeric: size of points. Default to 4.
#' @param point.shape Integer: shape of points. Default to 19.
#' @param point.alpha Numeric: transparency of points, a value between 0 and 1. Default to 0.7.
#' @param colors Character vector containing colors for compared groups. Use only if pairwise.data = TRUE.
#' The number of colors must be equal to the number of groups. Default to c("#66c2a5", "#fc8d62").
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param legend.title Character: title of legend. Default to "Group".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param labels Logical: if TRUE (default), a regression statistics will be added to the plot using ggpimsc::stat_poly_eq() function.
#' @param label Character: specifies regression statistics to add, names specified by ggpimsc::stat_poly_eq() function must be used.
#' Default to c("eq", "R2", "p") for regression equation, coefficient of determination and p value, respectively.
#' @param label.position.x,label.position.y  Numeric: coordinates for position of regression statistics.
#' See description of label.x and label.y parameters from ggpimsc::stat_poly_eq() function.
#' @param small.p,small.r Logical, if TRUE, p character in p value label and r character in coefficient of determination label will be lowercase. Default to FALSE.
#' @param rr.digits,p.digits Integer: number of digits after the decimal point in coefficient of determination and p value in labels.
#' Default to 2 for rr.digits and 3 for p.digits.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "Sample1_Sample2_single_pair_plot".
#'
#' @return Object with plot. Created plot is also displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(ggpmisc)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' single_pair_sample(data.dCt, "Disease6", "Control17")
#'
#' @importFrom dplyr select
#' @importFrom ggplot2 ggplot geom_point geom_smooth xlab ylab labs theme_classic theme element_text ggsave
#' @importFrom ggpmisc stat_poly_eq use_label
#' @import ggplot2
#' @import ggpmisc
#'
single_pair_sample <- function(data,
                               x,
                               y,
                               pairwise.data = FALSE,
                               by.group = FALSE,
                               point.size = 4,
                               point.shape = 19,
                               point.alpha = 0.7,
                               colors = c("#66c2a5", "#fc8d62"),
                               axis.title.size = 11,
                               axis.text.size = 10,
                               legend.title = "Group",
                               legend.title.size = 11,
                               legend.text.size = 11,
                               legend.position = "right",
                               plot.title = "",
                               plot.title.size = 14,
                               labels = TRUE,
                               label = c("eq", "R2", "p"),
                               label.position.x = 1,
                               label.position.y = 1,
                               small.p = FALSE,
                               small.r = FALSE,
                               p.digits = 3,
                               rr.digits = 2,
                               save.to.tiff = FALSE,
                               dpi = 600,
                               width = 15,
                               height = 15,
                               name.tiff = "samples_single_pair_plot") {
  if (pairwise.data == TRUE) {
    data <- as.data.frame(data)
    data <-
      pivot_longer(data,
                   !c(Sample, Group),
                   names_to = "Gene" ,
                   values_to = "value")
    data_t <-
      pivot_wider(data, names_from = Sample, values_from = value)

  } else {
    data <- as.data.frame(data)
    data_t <- t(select(data,-Group,-Sample))
    colnames(data_t) <- data$Sample
  }

  if (by.group == FALSE) {
    single_pair_t <-
      ggplot(as.data.frame(data_t), aes(x = .data[[x]], y = .data[[y]])) +
      geom_point(
        size = point.size,
        shape = point.shape,
        alpha = point.alpha,
        color = colors[1]
      ) +
      geom_smooth(method = 'lm', se = FALSE) +
      xlab(x) +
      ylab(y) +
      labs(title = plot.title) +
      theme_classic() +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
      theme(plot.title = element_text(size = plot.title.size))

  } else {
    single_pair_t <-
      ggplot(data_t, aes(x = .data[[x]], y = .data[[y]], color = Group)) +
      geom_point(size = point.size,
                 shape = point.shape,
                 alpha = point.alpha) +
      geom_smooth(method = 'lm', se = FALSE) +
      scale_color_manual(values = c(colors)) +
      xlab(x) +
      ylab(y) +
      labs(color = legend.title, title = plot.title) +
      theme_classic() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
      theme(legend.text = element_text(size = legend.text.size, colour =
                                         "black")) +
      theme(legend.title = element_text(size = legend.title.size, colour =
                                          "black")) +
      theme(plot.title = element_text(size = plot.title.size))
  }

  if (labels == TRUE) {
    single_pair_t <- single_pair_t +
      stat_poly_eq(
        use_label(label),
        label.y = c(label.position.y),
        label.x = c(label.position.x),
        small.p = small.p,
        small.r = small.r,
        p.digits = p.digits,
        rr.digits = rr.digits
      )
  }

  print(single_pair_t)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(x, "_", y, name.tiff, ".tiff", sep = ""),
      single_pair_t,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }
  return(single_pair_t)
}





#' @title filter_transformed_data
#'
#' @description
#' This function filters transformed Ct data (filtered Ct, delta Ct, and 2^-dCt data)
#' according to the used filtering criteria (see parameters).
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' @param remove.Gene Character: vector with names of genes to remove from data.
#' @param remove.Sample Character: vector with names of samples to remove from data.
#' @param remove.Group Character: vector with names of groups to remove from data.
#'
#' @return Data frame with filtered data.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' data.dCtF <- filter_transformed_data(data.dCt, remove.Sample = c("Control11"))
#'
#' dim(data.dCt)
#' dim(data.dCtF)
#'
#'
#' @importFrom dplyr filter select
#' @importFrom tidyselect any_of
#' @importFrom magrittr %>%
#' @import tidyverse
#'
filter_transformed_data <- function(data,
                                    remove.Gene = c(""),
                                    remove.Sample = c(""),
                                    remove.Group = c("")) {
  data <- filter(data,!Sample %in% remove.Sample,!Group %in% remove.Group)

  data <- select(data,-any_of(remove.Gene))

  return(data)
}





#' @title RQ_dCt
#'
#' @description
#' This function performs relative quantification of gene expression using 2^-dCt method.
#'
#' @details
#' This function calculates:
#' 1. Means (returned in columns with the "_mean" pattern) and standard deviations (returned in columns with the "_sd" pattern)
#' of 2^-dCt transformed dCt values for each analyzed gene across compared groups.
#' 2. P values of normality test (Shapiro_Wilk test) performed on 2^-dCt values across compared groups (returned in
#' columns with the "_norm_p" pattern).
#' 3. Fold change values (returned in "FCh" column) calculated for each gene by dividing  mean of 2^-dCt values in study group
#' by mean of 2^-dCt values in reference group.
#' 4. Statistics (returned in column with the "_test_stat" pattern) and p values (returned in column with "_test_p" pattern) of
#' differences in 2^-dCt values between study group and reference group using Student's t test and Mann-Whitney U test.
#' 5. P values adjusted for multiple testing using a selected method.
#' @param data Data object returned from delta_Ct() or filter_transformed_data() function.
#' @param group.study Character: name of study group (group of interest).
#' @param group.ref Character: name of reference group.
#' @param do.tests Logical: if TRUE, statistical significance of differences between compared groups will be calculated using Student's t test and Mann-Whitney U test. Default to TRUE.
#' @param pairwise Logical: if TRUE, a pairwise analysis will be performed (see details). Default to FALSE.
#' @param alternative Character: alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param p.adjust.method Character: p value correction method for multiple testing, one of the "holm", "hochberg", "hommel",
#' "bonferroni", "BH" (default), "BY","fdr", or "none". See documentation for stats::p.adjust() function for details.
#' @param save.to.txt Logical: if TRUE, returned table with results will be saved to .txt file. Default to FALSE.
#' @param name.txt Character: name of saved .txt file, without ".txt" name of extension. Default to "results_dCt".
#
#' @return Data frame with results (if pairwise = FALSE) or, if pairwise = TRUE, the list object
#' with two elements: a table with the results and the second table with fold change values calculated individually for each sample.
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(coin)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                      remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                      remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8", transform = TRUE)
#' results.dCt <- RQ_dCt(data.dCt,
#'                            group.study = "Disease",
#'                            group.ref = "Control")
#' head(results.dCt)
#'
#' @importFrom stats sd shapiro.test t.test p.adjust
#' @importFrom coin wilcox_test wilcoxsign_test pvalue statistic
#' @importFrom utils write.table
#' @importFrom dplyr filter select rename_with full_join group_by summarise mutate
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tidyselect all_of ends_with everything
#' @importFrom magrittr %>%
#' @import tidyverse
#'
RQ_dCt <- function(data,
                          group.study,
                          group.ref,
                          do.tests = TRUE,
                          pairwise = FALSE,
                          alternative = "two.sided",
                          p.adjust.method = "BH",
                          save.to.txt = FALSE,
                          name.txt = "results_dCt") {
   data_slim <- data %>%
    filter(Group == group.study | Group == group.ref) %>%
    pivot_longer(
      cols = -c(Group, Sample),
      names_to = "Gene",
      values_to = "value"
    )

  data_mean <- data_slim %>%
    group_by(Group, Gene) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = value) %>%
    rename_with(~ paste0(.x, "_mean", recycle0 = TRUE), all_of(c(group.study, group.ref)))


  data_sd <- data_slim %>%
    group_by(Group, Gene) %>%
    summarise(value_sd = sd(value, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = value_sd) %>%
    rename_with(~ paste0(.x, "_sd", recycle0 = TRUE), all_of(c(group.study, group.ref)))

  data_mean_sd <- full_join(data_mean, data_sd, by = c("Gene"))


  if (pairwise == FALSE) {
    data_FCh <- data_slim %>%
      group_by(Group, Gene) %>%
      summarise(value = mean(value, na.rm = TRUE), .groups = "keep") %>%
      pivot_wider(names_from = Group, values_from = value) %>%
      mutate(FCh = .data[[group.study]] / .data[[group.ref]]) %>%
      mutate(log10FCh = log10(FCh)) %>%
      select(Gene, FCh, log10FCh)

  } else {
    data_FCh <- data %>%
      filter(Group == group.study | Group == group.ref) %>%
      pivot_longer(
        cols = -c(Group, Sample),
        names_to = "Gene",
        values_to = "value"
      ) %>%
      pivot_wider(names_from = Group, values_from = value) %>%
      mutate(FCh = .data[[group.study]] / .data[[group.ref]])

    data_FCh_mean <- data_FCh %>%
      group_by(Gene) %>%
      summarise(FCh = mean(FCh, na.rm = TRUE), .groups = "keep") %>%
      mutate(log10FCh = log10(FCh))

    data_FCh_sd <- data_FCh %>%
      group_by(Gene) %>%
      summarise(FCh_sd = sd(FCh, na.rm = TRUE), .groups = "keep")
  }

  if (do.tests == TRUE) {
    data_norm <- data_slim %>%
      group_by(Group, Gene) %>%
      summarise(shap_wilka_p = shapiro.test(value)$p.value,
                .groups = "keep") %>%
      pivot_wider(names_from = Group, values_from = shap_wilka_p) %>%
      rename_with(~ paste0(.x, "_norm_p", recycle0 = TRUE), all_of(c(group.study, group.ref)))

    data_mean_sd_norm <-
      full_join(data_mean_sd, data_norm, by = c("Gene"))

    data_slim$Group <- as.factor(data_slim$Group)

    if (pairwise == FALSE) {
      data_tests <- data_slim %>%
        group_by(Gene) %>%
        summarise(
          t_test_p = t.test(value ~ Group, alternative = alternative)$p.value,
          t_test_stat = t.test(value ~ Group, alternative = alternative)$statistic,
          MW_test_p = coin::pvalue(wilcox_test(value ~ Group, alternative = alternative)),
          MW_test_stat = coin::statistic(wilcox_test(value ~ Group, alternative = alternative)),
          .groups = "keep"
        )

    } else {
      data_tests <- data_slim %>%
        pivot_wider(names_from = "Group", values_from = "value") %>%
        group_by(Gene) %>%
        summarise(
          t_test_p = t.test(
            .data[[group.study]],
            .data[[group.ref]],
            alternative = alternative,
            paired = TRUE
          )$p.value,
          t_test_stat = t.test(
            .data[[group.study]],
            .data[[group.ref]],
            alternative = alternative,
            paired = TRUE
          )$statistic,
          MW_test_p = coin::pvalue(
            wilcoxsign_test(.data[[group.study]] ~ .data[[group.ref]], alternative = alternative)
          ),
          MW_test_stat = coin::statistic(
            wilcoxsign_test(.data[[group.study]] ~ .data[[group.ref]], alternative = alternative)
          ),
          .groups = "keep"
        )
    }
    data_tests$t_test_p_adj <-
      p.adjust(data_tests$t_test_p, method = p.adjust.method)
    data_tests$MW_test_p_adj <-
      p.adjust(data_tests$MW_test_p, method = p.adjust.method)

    if (pairwise == TRUE) {
      data_mean_sd_norm_FChmean <-
        full_join(data_mean_sd_norm, data_FCh_mean, by = c("Gene"))
      data_mean_sd_norm_FChmean_FChsd <-
        full_join(data_mean_sd_norm_FChmean, data_FCh_sd, by = c("Gene"))
      data_results <-
        full_join(data_mean_sd_norm_FChmean_FChsd,
                  data_tests,
                  by = c("Gene"))

      if (save.to.txt == TRUE) {
        write.table(as.data.frame(data_results), paste(name.txt, ".txt", sep = ""))
      } else {}

      return(list(data_results, data_FCh))

    } else {
      data_mean_sd_FCh <-
        full_join(data_mean_sd_norm, data_FCh, by = c("Gene"))
      data_results <-
        full_join(data_mean_sd_FCh, data_tests, by = c("Gene"))

      if (save.to.txt == TRUE) {
        write.table(as.data.frame(data_results), paste(name.txt, ".txt", sep = ""))
      } else {}

      return(data_results)
    }

  } else {
    if (pairwise == TRUE) {
      data_mean_sd_FChmean <-
        full_join(data_mean_sd, data_FCh_mean, by = c("Gene"))
      data_results <-
        full_join(data_mean_sd_FChmean, data_FCh_sd, by = c("Gene"))

      if (save.to.txt == TRUE) {
        write.table(as.data.frame(data_results), paste(name.txt, ".txt", sep = ""))
      } else {}

      return(list(data_results, data_FCh))

    } else {
      data_results <- full_join(data_mean_sd, data_FCh, by = c("Gene"))

      if (save.to.txt == TRUE) {
        write.table(as.data.frame(data_results), paste(name.txt, ".txt", sep = ""))
      } else {}

      return(data_results)
    }
  }
}





#' @title results_boxplot
#'
#' @description
#' This function creates boxplot that illustrate distribution of data for selected genes.
#' It is similar to control_boxplot_gene() function; however, some new options are added,
#' including gene selection, faceting, adding mean labels to boxes, and adding statistical significance labels.
#' This function can be used to present results for finally selected genes.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' @param coef Numeric: how many times of interquartile range should be used to determine range point for whiskers. Default to 1.5.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all genes.
#' @param by.group Logical: if TRUE (default), distributions will be drawn by compared groups of samples.
#' @param signif.show Logical: if TRUE, labels for statistical significance will be added to the plot. Default to FALSE.
#' @param signif.labels Character vector with statistical significance labels (e.g. "ns.","***", etc.). The number
#' of elements must be equal to the number of genes used for plotting. All elements in the vector must be different; therefore,
#' symmetrically white spaces to repeated labels must be add to the same labels, e.g. "ns.", " ns. ", "  ns.  ".
#' @param signif.length Numeric: length of horizontal bars under statistical significance labels, values from 0 to 1.
#' @param signif.dist Numeric: distance between the highest value and statistical significance labels.
#' Can be in y axis units (if faceting = FALSE) or fraction of y axis value reached by the highest value (if faceting = TRUE).
#' @param faceting Logical: if TRUE (default), plot will be drawn with facets with free scales.
#' @param facet.row,facet.col Integer: number of rows and columns to arrange facets.
#' @param angle Integer: value of angle in which names of genes are displayed. Default to 0.
#' @param y.exp.low,y.exp.up Numeric: space between data on the plot and lower and upper axis. Useful to add extra space for statistical significance labels when faceting = TRUE.
#' @param rotate Logical: if TRUE, boxplots will be arranged horizontally. Default to FALSE.
#' @param add.mean Logical: if TRUE, mean points will be added to boxes as squares. Default to TRUE.
#' @param add.mean.size Numeric: size of squares indicating means. Default to 2.
#' @param add.mean.color Character: color of squares indicating means. Default to "black".
#' @param colors Character vector length of one (when by.group = FALSE) or more (when by.group = TRUE), containing colors for groups.
#' The number of colors must be equal to the number of groups. Default to c("#66c2a5", "#fc8d62").
#' @param x.axis.title Character: title of x axis. Default to "Gene".
#' @param y.axis.title Character: title of y axis. Default to "value".
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param legend.title Character: title of legend. Default to "Group".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "results_boxplot".
#'
#' @return Object with plot. Created plot will be also displayed on the graphic device.
#' @export
#'
#' @examples
#' library(ggsignif)
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' data.dCtF <- filter_transformed_data(data.dCt, remove.Sample = c("Control11"))
#' results_boxplot(data.dCtF,
#'                 sel.Gene = c("Gene1","Gene16","Gene19","Gene20"),
#'                 signif.labels = c("****","*","***"," * "),
#'                 angle = 30,
#'                 signif.dist = 1.05,
#'                 facet.row = 1,
#'                 facet.col = 4,
#'                 y.exp.up = 0.1,
#'                 y.axis.title = bquote(~2^-dCt))
#'
#' @importFrom dplyr select filter group_by summarise
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot geom_boxplot scale_fill_manual coord_flip guides xlab ylab labs theme_classic theme element_text ggsave scale_y_continuous expansion facet_wrap element_blank stat_summary
#' @importFrom ggsignif geom_signif
#' @import ggplot2
#' @import tidyverse
#'
results_boxplot <- function(data,
                            coef = 1.5,
                            sel.Gene = "all",
                            by.group = TRUE,
                            signif.show = FALSE,
                            signif.labels,
                            signif.length = 0.2,
                            signif.dist = 0.2,
                            faceting = TRUE,
                            facet.row,
                            facet.col,
                            y.exp.low = 0.1,
                            y.exp.up = 0.2,
                            angle = 0,
                            rotate = FALSE,
                            add.mean = TRUE,
                            add.mean.size = 2,
                            add.mean.color = "black",
                            colors = c("#66c2a5", "#fc8d62"),
                            x.axis.title = "",
                            y.axis.title = "value",
                            axis.title.size = 11,
                            axis.text.size = 10,
                            legend.text.size = 11,
                            legend.title = "Group",
                            legend.title.size = 11,
                            legend.position = "top",
                            plot.title = "",
                            plot.title.size = 14,
                            save.to.tiff = FALSE,
                            dpi = 600,
                            width = 15,
                            height = 15,
                            name.tiff = "results_boxplot") {
  data <-
    pivot_longer(data,
                 !c(Sample, Group),
                 names_to = "Gene" ,
                 values_to = "value")

  if (sel.Gene[1] == "all") {
    data <- data

  } else {
    data <- filter(data, Gene %in% sel.Gene)
  }

  if (by.group == TRUE) {
    box_results <- ggplot(data, aes(x = Gene, y = value)) +
      geom_boxplot(aes(fill = Group), coef = coef) +
      scale_fill_manual(values = c(colors)) +
      labs(fill = legend.title, title = plot.title)

    if (signif.show == TRUE) {
      label.height <- data %>%
        group_by(Gene) %>%
        summarise(height = max(value), .groups = "keep")

      data.label <-
        data.frame(matrix(nrow = length(unique(
          label.height$Gene
        )), ncol = 4))
      rownames(data.label) <- label.height$Gene
      colnames(data.label) <- c("x", "xend", "y", "annotation")


      if (faceting == TRUE) {
        data.label$x <- rep(1 - signif.length, nrow(data.label))
        data.label$xend <- rep(1 + signif.length, nrow(data.label))
        data.label$y <- label.height$height * signif.dist
      } else {
        data.label$x <- (1:nrow(data.label)) - signif.length
        data.label$xend <- (1:nrow(data.label)) + signif.length
        data.label$y <- label.height$height + signif.dist

      }
      data.label$annotation <- signif.labels
      data.label$Gene <- rownames(data.label)

      box_results <- box_results +
        geom_signif(
          stat = "identity",
          data = data.label,
          aes(
            x = x,
            xend = xend,
            y = y,
            yend = y,
            annotation = annotation
          ),
          color = "black",
          manual = TRUE
        ) +
        scale_y_continuous(expand = expansion(mult = c(y.exp.low, y.exp.up)))

    }

    if (faceting == TRUE) {
      box_results <- box_results +
        facet_wrap(vars(Gene),
                   scales = "free",
                   nrow = facet.row,
                   ncol = facet.col)
    }
  } else {
    box_results <- ggplot(data, aes(x = Gene, y = value)) +
      geom_boxplot(coef = coef, fill = colors[1]) +
      theme(axis.text.x = element_text(size = axis.text.size, colour = "black")) +
      labs(title = plot.title)

    if (faceting == TRUE) {
      box_results <- box_results +
        facet_wrap(vars(Gene),
                   scales = "free",
                   nrow = facet.row,
                   ncol = facet.col)
    }
  }

  box_results <- box_results +
    guides(x =  guide_axis(angle = angle)) +
    xlab(x.axis.title) +
    ylab(y.axis.title) +
    theme_bw() +
    theme(legend.position = legend.position) +
    theme(axis.text.y = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
    theme(legend.text = element_text(size = legend.text.size, colour = "black")) +
    theme(legend.title = element_text(size = legend.title.size, colour = "black")) +
    theme(plot.title = element_text(size = plot.title.size)) +
    theme(panel.grid.major.x = element_blank())


  if (rotate == TRUE) {
    box_results <- box_results +
      coord_flip()
  }

  if (add.mean == TRUE & by.group == TRUE) {
    box_results <- box_results +
      stat_summary(
        aes(group = Group),
        fun = mean,
        position = position_dodge(width = .75),
        geom = "point",
        shape = 15,
        size = add.mean.size,
        color = add.mean.color
      )
  }
  if (add.mean == TRUE & by.group == FALSE) {
    box_results <- box_results +
      stat_summary(
        fun = mean,
        position = position_dodge(width = .75),
        geom = "point",
        shape = 15,
        size = add.mean.size,
        color = add.mean.color
      )
  }


  if (faceting == TRUE) {
    box_results <- box_results +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }


  print(box_results)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      box_results,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }
  return(box_results)
}







#' @title results_barplot
#'
#' @description
#' This function creates a barplot that illustrate mean and sd values of genes.
#' Faceting and adding custom labels of statistical significance are available.
#' This function is useful to present results for finally selected genes.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all genes.
#' @param bar.width Numeric: width of bars.
#' @param signif.show Logical: if TRUE, labels for statistical significance will be added to the plot. Default to FALSE.
#' @param signif.labels Character vector with statistical significance labels (e.g. "ns.","***", etc.). The number
#' of elements must be equal to the number of genes used for plotting. All elements in the vector must be different; therefore,
#' symmetrically white spaces to repeated labels must be add to the same labels, e.g. "ns.", " ns. ", "  ns.  ".
#' @param signif.length Numeric: length of horizontal bars under statistical significance labels, values from 0 to 1.
#' @param signif.dist Numeric: distance between errorbar and statistical significance labels.
#' Can be in y axis units (if faceting = FALSE) or fraction of y axis value reached by errorbar (mean + sd value) (if faceting = TRUE).
#' @param faceting Logical: if TRUE (default), plot will be drawn with facets with free scales.
#' @param facet.row,facet.col Integer: number of rows and columns to arrange facets.
#' @param y.exp.low,y.exp.up Numeric: space between data on the plot and lower or upper axis. Useful to add extra space for statistical significance labels when faceting = TRUE.
#' @param angle Integer: value of angle in which names of genes are displayed. Default to 0.
#' @param rotate Logical: if TRUE, boxplots will be arranged horizontally. Default to FALSE.
#' @param colors Character vector length of one (when by.group = FALSE) or two (when by.group = TRUE), containing colors for groups.
#' The numbers of colors must be equal to the number of groups. Default to c("#66c2a5", "#fc8d62").
#' @param x.axis.title Character: title of x axis. Default to "Gene".
#' @param y.axis.title character: title of y axis. Default to "value".
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param legend.title Character: title of legend. Default to "Group".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "results_barplot".
#'
#' @return Object with plot. Created plot will be also displayed on graphic device.
#' @export
#'
#' @examples
#' library(ggsignif)
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' data.dCtF <- filter_transformed_data(data.dCt, remove.Sample = c("Control11"))
#' results_barplot(data.dCtF,
#'                 sel.Gene = c("Gene1","Gene16","Gene19","Gene20"),
#'                 signif.labels = c("****","*","***"," * "),
#'                 angle = 30,
#'                 signif.dist = 1.05,
#'                 facet.row = 1,
#'                 facet.col = 4,
#'                 y.exp.up = 0.1,
#'                 y.axis.title = bquote(~2^-dCt))
#'
#' @importFrom dplyr select filter group_by summarise mutate
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot geom_col geom_errorbar scale_fill_manual coord_flip guides xlab ylab labs theme_classic theme element_text ggsave scale_y_continuous expansion facet_wrap element_blank stat_summary
#' @importFrom ggsignif geom_signif
#' @import ggplot2
#' @import tidyverse
#'
results_barplot <- function(data,
                            sel.Gene = "all",
                            bar.width = 0.8,
                            signif.show = FALSE,
                            signif.labels,
                            signif.length = 0.2,
                            signif.dist = 0.2,
                            faceting = FALSE,
                            facet.row,
                            facet.col,
                            y.exp.low = 0.1,
                            y.exp.up = 0.2,
                            angle = 0,
                            rotate = FALSE,
                            colors = c("#66c2a5", "#fc8d62"),
                            x.axis.title = "",
                            y.axis.title = "value",
                            axis.title.size = 11,
                            axis.text.size = 10,
                            legend.text.size = 11,
                            legend.title = "Group",
                            legend.title.size = 11,
                            legend.position = "top",
                            plot.title = "",
                            plot.title.size = 14,
                            save.to.tiff = FALSE,
                            dpi = 600,
                            width = 15,
                            height = 15,
                            name.tiff = "results_barplot") {
  data <-
    pivot_longer(data,
                 !c(Sample, Group),
                 names_to = "Gene" ,
                 values_to = "value")

  if (sel.Gene[1] == "all") {
    data <- data

  } else {
    data <- filter(data, Gene %in% sel.Gene)
  }

  data.mean <- data %>%
    group_by(Group, Gene) %>%
    summarise(mean = mean(value, na.rm = TRUE), .groups = "keep")

  data.sd <- data %>%
    group_by(Group, Gene) %>%
    summarise(sd = sd(value, na.rm = TRUE), .groups = "keep")

  data.mean$sd <- data.sd$sd

  bar_results <- ggplot(data.mean, aes(x = Gene, y = mean)) +
    geom_errorbar(
      aes(
        group = Group,
        y = mean,
        ymin = ifelse(mean < 0, mean - abs(sd), mean),
        ymax = ifelse(mean > 0, mean + abs(sd), mean)
      ),
      width = .2,
      position = position_dodge(.9)
    ) +
    geom_col(
      aes(fill = Group, group = Group),
      position = position_dodge(.9),
      width = bar.width,
      color = "black"
    ) +
    scale_fill_manual(values = c(colors)) +
    guides(x =  guide_axis(angle = angle)) +
    xlab(x.axis.title) +
    ylab(y.axis.title) +
    labs(fill = legend.title, title = plot.title) +
    theme_bw() +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
    theme(legend.text = element_text(size = legend.text.size, colour =
                                       "black")) +
    theme(legend.title = element_text(size = legend.title.size, colour =
                                        "black")) +
    theme(plot.title = element_text(size = plot.title.size)) +
    theme(panel.grid.major.x = element_blank())

  if (faceting == TRUE) {
    bar_results <- bar_results +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
      facet_wrap(vars(Gene),
                 scales = "free",
                 nrow = facet.row,
                 ncol = facet.col)
  }

  if (signif.show == TRUE) {
    label.height <- data.mean %>%
      mutate(max = mean + sd) %>%
      group_by(Gene) %>%
      summarise(height = max(max, na.rm = TRUE), .groups = "keep")

    data.label <-
      data.frame(matrix(nrow = length(unique(data.mean$Gene)), ncol = 4))
    rownames(data.label) <- unique(data.mean$Gene)
    colnames(data.label) <- c("x", "xend", "y", "annotation")

    if (faceting == TRUE) {
      data.label$x <-
        rep(1 - signif.length, length(unique(data.mean$Gene)))
      data.label$xend <-
        rep(1 + signif.length, length(unique(data.mean$Gene)))
      data.label$y <- label.height$height * signif.dist

    } else {
      data.label$x <- (1:length(unique(data.mean$Gene))) - signif.length
      data.label$xend <-
        (1:length(unique(data.mean$Gene))) + signif.length
      data.label$y <- label.height$height + signif.dist
    }

    data.label$annotation <- signif.labels
    data.label$Gene <- rownames(data.label)

    bar_results <- bar_results +
      geom_signif(
        stat = "identity",
        data = data.label,
        aes(
          x = x,
          xend = xend,
          y = y,
          yend = y,
          annotation = annotation
        ),
        color = "black",
        manual = TRUE
      ) +
      scale_y_continuous(expand = expansion(mult = c(y.exp.low, y.exp.up)))
  }

  if (rotate == TRUE) {
    bar_results <- bar_results +
      coord_flip()
  }

  print(bar_results)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      bar_results,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }
  return(bar_results)
}




#' @title RQ_ddCt
#'
#' @description
#' This function performs relative quantification of gene expression using 2^-ddCt method.
#'
#' @details
#' This function performs:
#' 1. Calculation of means (returned in columns with the "_mean" pattern) and standard deviations (returned
#' in columns with the "_sd" pattern) of delta Ct values for analyzed genes across compared groups.
#' 2. Normality tests (Shapiro_Wilk test) of delta Ct values of analyzed genes across compared groups and returned p values are
#' stored in columns with the "_norm_p" pattern.
#' 3. Calculation of differences in mean delta Ct values of genes between compared groups, obtaining delta delta Ct values (returned  in"ddCt" column).
#' 4. Calculation of fold change values (returned in "FCh" column) for each gene by exponentiation of ddCt values using 2^-ddCt formula.
#' 5. Statistical testing of differences between study group and reference group using Student's t test and Mann-Whitney U test.
#' Resulted statistics (in column with "_test_stat" pattern) and p values (in column with "_test_p" pattern) are returned.
#' 6. P values adjusted for multiple testing using a selected method.
#' @param data Data object returned from delta_Ct() function.
#' @param group.study Character: name of study group (group of interest).
#' @param group.ref Character: name of reference group.
#' @param do.tests Logical: if TRUE, statistical significance of delta delta Ct values between compared groups
#' is calculated using Student's t test and Mann-Whitney U test. Default to TRUE.
#' @param pairwise Logical: set to TRUE if a pairwise analysis should be done.
#' @param alternative Character: alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param p.adjust.method Character: p value correction method for multiple testing, one of the "holm", "hochberg", "hommel",
#' "bonferroni", "BH" (default), "BY","fdr", or "none". See documentation for stats::p.adjust() function for details.
#' @param save.to.txt Logical: if TRUE, returned table with results is saved to .txt file. Default to FALSE.
#' @param name.txt Character: name of saved .txt file, without ".txt" name of extension. Default to "results_ddCt".
#
#' @return Data frame with relative quantification results.
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(coin)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' results.ddCt <- RQ_ddCt(data.dCt, "Disease", "Control")
#' head(results.ddCt)
#'
#' @importFrom stats sd shapiro.test t.test p.adjust Pair
#' @importFrom coin wilcox_test wilcoxsign_test pvalue statistic
#' @importFrom utils write.table
#' @importFrom dplyr filter select rename_with full_join group_by summarise mutate
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tidyselect all_of ends_with everything
#' @importFrom magrittr %>%
#' @import tidyverse
#'
RQ_ddCt <- function(data,
                    group.study,
                    group.ref,
                    do.tests = TRUE,
                    pairwise = FALSE,
                    alternative = "two.sided",
                    p.adjust.method = "BH",
                    save.to.txt = FALSE,
                    name.txt = "ddCt_RQ_results") {
  data_slim <- data %>%
    filter(Group == group.study | Group == group.ref) %>%
    pivot_longer(
      cols = -c(Group, Sample),
      names_to = "Gene",
      values_to = "dCt"
    )

  data_slim$Group <- as.factor(data_slim$Group)

  data_mean <- data_slim %>%
    group_by(Group, Gene) %>%
    summarise(value = mean(dCt, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = value) %>%
    rename_with(~ paste0(.x, "_mean", recycle0 = TRUE), all_of(c(group.study, group.ref)))

  data_sd <- data_slim %>%
    group_by(Group, Gene) %>%
    summarise(value_sd = sd(dCt, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = value_sd) %>%
    rename_with(~ paste0(.x, "_sd", recycle0 = TRUE), all_of(c(group.study, group.ref)))

  data_mean_sd <- full_join(data_mean, data_sd, by = c("Gene"))

  data_slim$Group <- as.factor(data_slim$Group)

  if (pairwise == FALSE) {
    data_ddCt <- data_slim %>%
      group_by(Group, Gene) %>%
      summarise(ddCt = mean(dCt, na.rm = TRUE), .groups = "keep") %>%
      pivot_wider(names_from = Group, values_from = ddCt) %>%
      mutate(ddCt = .data[[group.study]] - .data[[group.ref]]) %>%
      mutate(FCh = 2 ^ -ddCt) %>%
      mutate(log10FCh = log10(FCh)) %>%
      select(Gene, ddCt, FCh, log10FCh)

    data_ddCt_sd <- data_slim %>%
      group_by(Group, Gene) %>%
      summarise(dCt_sd = sd(dCt, na.rm = TRUE), .groups = "keep") %>%
      pivot_wider(names_from = Group, values_from = dCt_sd) %>%
      rename_with( ~ paste0(.x, "_sd", recycle0 = TRUE), all_of(c(group.study, group.ref)))

  } else {
    data_FCh <- data %>%
      filter(Group == group.study | Group == group.ref) %>%
      pivot_longer(
        cols = -c(Group, Sample),
        names_to = "Gene",
        values_to = "value"
      ) %>%
      pivot_wider(names_from = Group, values_from = value) %>%
      mutate(ddCt = .data[[group.study]] - .data[[group.ref]]) %>%
      mutate(FCh = 2 ^ -ddCt) %>%
      select(-ddCt)

    data_FCh_mean <- data_FCh %>%
      group_by(Gene) %>%
      summarise(FCh = mean(FCh, na.rm = TRUE), .groups = "keep") %>%
      mutate(log10FCh = log10(FCh))

    data_FCh_sd <- data_FCh %>%
      group_by(Gene) %>%
      summarise(FCh_sd = sd(FCh, na.rm = TRUE), .groups = "keep")

  }

  if (do.tests == TRUE) {
    data_norm <- data_slim %>%
      group_by(Group, Gene) %>%
      summarise(shap_wilka_p = shapiro.test(dCt)$p.value,
                .groups = "keep") %>%
      pivot_wider(names_from = Group, values_from = shap_wilka_p) %>%
      rename_with(~ paste0(.x, "_norm_p", recycle0 = TRUE), all_of(c(group.study, group.ref)))

    data_mean_sd_norm <-
      full_join(data_mean_sd, data_norm, by = c("Gene"))


    if (pairwise == FALSE) {
      data_tests <- data_slim %>%
        group_by(Gene) %>%
        summarise(
          t_test_p = t.test(dCt ~ Group, alternative = alternative)$p.value,
          t_test_stat = t.test(dCt ~ Group, alternative = alternative)$statistic,
          MW_test_p = coin::pvalue(wilcox_test(dCt ~ Group, alternative = alternative)),
          MW_test_stat = coin::statistic(wilcox_test(dCt ~ Group, alternative = alternative)),
          .groups = "keep"
        )
    } else {
      data_tests <- data_slim %>%
        pivot_wider(names_from = "Group", values_from = "dCt") %>%
        group_by(Gene) %>%
        summarise(
          t_test_p = t.test(
            .data[[group.study]],
            .data[[group.ref]],
            alternative = alternative,
            paired = TRUE
          )$p.value,
          t_test_stat = t.test(
            .data[[group.study]],
            .data[[group.ref]],
            alternative = alternative,
            paired = TRUE
          )$statistic,
          MW_test_p = coin::pvalue(
            wilcoxsign_test(.data[[group.study]] ~ .data[[group.ref]], alternative = alternative)
          ),
          MW_test_stat = coin::statistic(
            wilcoxsign_test(.data[[group.study]] ~ .data[[group.ref]], alternative = alternative)
          ),
          .groups = "keep"
        )
    }

    data_tests$t_test_p_adj <- p.adjust(data_tests$t_test_p,
                                        method = p.adjust.method)
    data_tests$MW_test_p_adj <- p.adjust(data_tests$MW_test_p,
                                         method = p.adjust.method)

    if (pairwise == TRUE) {
      data_mean_sd_norm_FChmean <-
        full_join(data_mean_sd_norm, data_FCh_mean, by = c("Gene"))
      data_mean_sd_norm_FChmean_FChsd <-
        full_join(data_mean_sd_norm_FChmean, data_FCh_sd, by = c("Gene"))
      data_results <-
        full_join(data_mean_sd_norm_FChmean_FChsd,
                  data_tests,
                  by = c("Gene"))

if (save.to.txt == TRUE) {
    write.table(as.data.frame(data_results), paste(name.txt, ".txt", sep = ""))
  } else {}


      return(list(data_results, data_FCh))

    } else {
      data_mean_sd_FCh <-
        full_join(data_mean_sd_norm, data_ddCt, by = c("Gene"))
      data_results <-
        full_join(data_mean_sd_FCh, data_tests, by = c("Gene"))

if (save.to.txt == TRUE) {
    write.table(as.data.frame(data_results), paste(name.txt, ".txt", sep = ""))
  } else {}


      return(data_results)
    }

  } else {
    if (pairwise == TRUE) {
      data_mean_sd_FChmean <-
        full_join(data_mean_sd, data_FCh_mean, by = c("Gene"))
      data_results <-
        full_join(data_mean_sd_FChmean, data_FCh_sd, by = c("Gene"))

if (save.to.txt == TRUE) {
    write.table(as.data.frame(data_results), paste(name.txt, ".txt", sep = ""))
  } else {}


      return(list(data_results, data_FCh))

    } else {
      data_results <- full_join(data_mean_sd, data_ddCt, by = c("Gene"))

if (save.to.txt == TRUE) {
    write.table(as.data.frame(data_results), paste(name.txt, ".txt", sep = ""))
  } else {}

      return(data_results)
    }
  }
}





#' @title FCh_plot
#'
#' @description
#' This function creates barplot that illustrate fold change values with indicating of significance by different colors of bars.
#'
#' @param data Object returned from RQ_dCt() or RQ_ddCt() functions.
#' @param use.p Logical: if TRUE, p value threshold will be used to label gene as significant.
#' @param mode Character: which p value should be used? One of the "t" (p values from Student's t test), "t.adj" (adjusted p values from Student's t test),
#' "mw" (p values from Mann-Whitney U test),"mw.adj" (adjusted p values from Mann-Whitney U test),
#' "depends" (if data in both compared groups were considered as derived from normal distribution
#' (p value from Shapiro_Wilk test > 0.05) - p values from Student's t test will be used for significance assignment,
#' otherwise p values from Mann-Whitney U test will be used), "depends.adj"
#' (if data in both compared groups were considered as derived from normal distribution
#' (p value from Shapiro_Wilk test > 0.05) - adjusted p values from Student's t test will be used for significance assignment,
#' otherwise adjusted p values from Mann-Whitney U test will be used), and "user" that can be used
#' the user intend to use another p values, e.g. obtained from other statistical test. In such situation, before run FCh_plot function, the user should prepare
#' data frame object named "user" that contains two columns, the first of them with Gene names and the second with p values.
#' The order of columns must be kept as described.
#' @param p.threshold Numeric: threshold of p values for statistical significance. Default to 0.05.
#' @param use.FCh Logical: if TRUE, the criterion of fold change will be also used for significance assignment of genes.
#' @param FCh.threshold Numeric: threshold of fold change values used for significance assignment of genes.
#' If is set to 2 (default), genes with 2-fold changed expression (increased or decreased) between groups will be assigned as significant.
#' @param use.sd Logical: if TRUE, errorbars with standard deviations will be added to the plot.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all names of genes.
#' @param bar.width numeric: width of bars.
#' @param signif.show Logical: if TRUE, labels for statistical significance will be added to the plot. Default to FALSE.
#' @param signif.labels Character vector with statistical significance labels (e.g. "ns.","***", etc.). The number
#' of elements must be equal to the number of genes used for plotting. All elements in the vector must be different; therefore,
#' symmetrically white spaces to repeated labels must be add to the same labels, e.g. "ns.", " ns. ", "  ns.  ".
#' @param signif.length Numeric: length of horizontal bars under statistical significance labels, values from 0 to 1.
#' @param signif.dist Numeric: distance between errorbar and statistical significance labels.
#' @param y.exp.low,y.exp.up Numeric: space between data on the plot and lower or upper axis. Useful to add extra space for statistical significance labels when faceting = TRUE.
#' @param angle Integer: value of angle in which names of genes are displayed. Default to 0.
#' @param rotate Logical: if TRUE, bars will be arranged horizontally. Deafault to FALSE.
#' @param colors Character vector length of one (when use.p = FALSE) or two (when use.p = TRUE), containing colors for significant and no significant genes.
#' @param border.color Character: color for borders of bars.
#' @param x.axis.title Character: title of x axis. Default to "Gene".
#' @param y.axis.title Character: title of y axis. Default to "value".
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param legend.title Character: title of legend. Default to "Group".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be one of "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "FCh_plot".
#'
#' @return List containing object with barplot and data frame with results. Created plot is also displayed on graphic device.
#' @export
#'
#' @examples
#' library(ggsignif)
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' data.dCtF <- filter_transformed_data(data.dCt, remove.Sample = c("Control11"))
#' results.ddCt <- RQ_ddCt(data.dCtF, "Disease", "Control")
#'
#' signif.labels <- c("****",
#'                    "**",
#'                    "ns.",
#'                    " ns. ",
#'                    "  ns.  ",
#'                    "   ns.   ",
#'                    "    ns.    ",
#'                    "     ns.     ",
#'                    "      ns.      ",
#'                    "       ns.       ",
#'                    "        ns.        ",
#'                    "         ns.         ",
#'                    "          ns.          ",
#'                    "***")
#' FCh.plot <- FCh_plot(results.ddCt,
#'                    mode = "depends",
#'                    use.FCh = TRUE,
#'                    FCh.threshold = 2.5,
#'                    signif.labels = signif.labels,
#'                    angle = 30)
#' head(FCh.plot[[2]])
#'
#' # with user p values - calculated using stats::wilcox.test() function:
#' user <- data.dCt %>%
#' pivot_longer(cols = -c(Group, Sample), names_to = "Gene", values_to = "dCt") %>%
#'   group_by(Gene) %>%
#'   summarise(MW_test_p = wilcox.test(dCt ~ Group)$p.value, .groups = "keep")
#'
#' FCh.plot <- FCh_plot(results.ddCt,
#'                    mode = "user",
#'                    use.FCh = TRUE,
#'                    FCh.threshold = 2,
#'                    signif.labels = signif.labels,
#'                    angle = 30)
#' head(FCh.plot[[2]])
#
#' @importFrom dplyr select filter group_by summarise mutate ungroup desc
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot geom_col geom_errorbar scale_fill_manual coord_flip guides xlab ylab labs theme_classic theme element_text ggsave scale_y_continuous expansion facet_wrap element_blank stat_summary
#' @importFrom ggsignif geom_signif
#' @import ggplot2
#' @import tidyverse
#'
FCh_plot <- function(data,
                    use.p = TRUE,
                    mode,
                    p.threshold = 0.05,
                    use.FCh = FALSE,
                    FCh.threshold = 2,
                    use.sd = FALSE,
                    sel.Gene = "all",
                    bar.width = 0.8,
                    signif.show = FALSE,
                    signif.labels,
                    signif.length = 0.2,
                    signif.dist = 0.5,
                    y.exp.low = 0.1,
                    y.exp.up = 0.1,
                    angle = 0,
                    rotate = FALSE,
                    colors = c("#66c2a5", "#fc8d62"),
                    border.color = "black",
                    x.axis.title = "",
                    y.axis.title = "log10(Fold change)",
                    axis.title.size = 11,
                    axis.text.size = 10,
                    legend.text.size = 11,
                    legend.title = "Selected as significant?",
                    legend.title.size = 11,
                    legend.position = "top",
                    plot.title = "",
                    plot.title.size = 14,
                    dpi = 600,
                    width = 15,
                    height = 15,
                    save.to.tiff = FALSE,
                    name.tiff = "FCh_plot") {
  if (sel.Gene[1] != "all") {
    data <- filter(data, Gene %in% sel.Gene)
  }

  if (use.sd == TRUE){
  data <- mutate(data, factor = FCh_sd / FCh)
}

  if (use.p == TRUE) {
    if (mode == "t") {
      data$p.used <- data$t_test_p
    }
    if (mode == "t.adj") {
      data$p.used <- data$t_test_p_adj
    }
    if (mode == "mw") {
      data$p.used <- data$MW_test_p
    }
    if (mode == "mw.adj") {
      data$p.used <- data$MW_test_p_adj
    }
    if (mode == "depends") {
      data <- ungroup(data)
      vars <- colnames(select(data, ends_with("norm_p")))
      data <-
        mutate(
          data,
          test.for.comparison = ifelse(
            .data[[vars[[1]]]] >= 0.05 &
              .data[[vars[[2]]]] >= 0.05,
            "t.student's.test",
            "Mann-Whitney.test"
          )
        )
      data <-
        mutate(data,
               p.used = ifelse(
                 test.for.comparison == "t.student's.test",
                 t_test_p,
                 MW_test_p
               ))
    }
    if (mode == "depends.adj") {
      data <- ungroup(data)
      vars <- colnames(select(data, ends_with("norm_p")))
      data <-
        mutate(
          data,
          test.for.comparison = ifelse(
            .data[[vars[[1]]]] >= 0.05 &
              .data[[vars[[2]]]] >= 0.05,
            "t.student's.test",
            "Mann-Whitney.test"
          )
        )
      data <-
        mutate(
          data,
          p.used = ifelse(
            test.for.comparison == "t.student's.test",
            t_test_p_adj,
            MW_test_p_adj
          )
        )
    }
    if (mode == "user") {
      colnames(user) <- c("Gene", "p.used")
      data <- full_join(data, user, by = c("Gene"))
    }
    if (use.FCh == TRUE) {
      data <-
        mutate(data,
               `Selected as significant?` = ifelse(
                 p.used > p.threshold,
                 yes = "No",
                 no = ifelse(abs(log10(FCh)) <  log10(FCh.threshold), "No", "Yes")
               ))
      data$`Selected as significant?` <-
        factor(data$`Selected as significant?`, levels = c("Yes", "No"))

    } else {
      data <-
        mutate(data,
               `Selected as significant?` = ifelse(p.used > p.threshold, yes = "No",  no = "Yes"))
      data$`Selected as significant?` <-
        factor(data$`Selected as significant?`, levels = c("Yes", "No"))
    }

    RQ <-
      ggplot(data, aes(x = reorder(Gene,-FCh), y = log10FCh)) +
      geom_col(
        aes(fill = `Selected as significant?`, group = `Selected as significant?`),
        width = bar.width,
        color = border.color
      ) +
      scale_fill_manual(values = c(colors)) +
      labs(fill = legend.title, title = plot.title)

  } else {
    RQ <- ggplot(data, aes(x = reorder(Gene,-FCh), y = log10FCh)) +
      geom_col(width = bar.width,
               fill = colors[1],
               color = border.color) +
      labs(title = plot.title)
  }

  if (use.sd == TRUE) {
    RQ <- RQ +
      geom_errorbar(aes(
        y = log10FCh,
        ymin = ifelse(log10FCh < 0, log10FCh + (log10FCh * factor), log10FCh),
        ymax = ifelse(log10FCh > 0, log10FCh + (log10FCh * factor), log10FCh),
        width = .2
      ))
  }

  RQ <- RQ +
    guides(x =  guide_axis(angle = angle)) +
    xlab(x.axis.title) +
    ylab(y.axis.title) +
    theme_bw() +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
    theme(legend.text = element_text(size = legend.text.size, colour = "black")) +
    theme(legend.title = element_text(size = legend.title.size, colour = "black")) +
    theme(plot.title = element_text(size = plot.title.size)) +
    theme(panel.grid.major.x = element_blank()) +
    geom_hline(yintercept = 0, linewidth = 0.4)


  if (signif.show == TRUE) {
    data.label <- data.frame(matrix(nrow = nrow(data), ncol = 4))
    colnames(data.label) <- c("x", "xend", "y", "annotation")
    data <- arrange(data, desc(FCh))
    data.label$x <- (1:nrow(data)) - signif.length
    data.label$xend <- (1:nrow(data)) + signif.length

    if (use.sd == TRUE) {
      data.label <- mutate(
        data.label,
        y = ifelse(
          data$log10FCh > 0,
          data$log10FCh + (data$log10FCh * data$factor) + signif.dist,
          data$log10FCh + (data$log10FCh * data$factor) - signif.dist
        )
      )
    } else {
      data.label <- mutate(data.label,
                           y = ifelse(
                             data$log10FCh > 0,
                             data$log10FCh + signif.dist,
                             data$log10FCh - signif.dist
                           ))
    }
    data.label$annotation <- signif.labels

    RQ <- RQ +
      geom_signif(
        stat = "identity",
        data = data.label,
        aes(
          x = x,
          xend = xend,
          y = y,
          yend = y,
          annotation = annotation
        ),
        color = "black",
        manual = TRUE
      ) +
      scale_y_continuous(expand = expansion(mult = c(y.exp.low, y.exp.up)))

  }

  if (rotate == TRUE) {
    RQ <- RQ +
      coord_flip()
  }

  print(RQ)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      RQ,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }

  if (use.sd == TRUE){
    data <- select(data, -factor)
  }

  return(list(RQ, data))
}





#' @title ROCh
#'
#' @description
#' This function is designed to perform Receiver Operating Characteristic (ROC) analysis based on the gene expression data.
#' This kind of analysis is useful to further examine performance of samples classification into two groups.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all genes.
#' @param groups Character vector length of two with names of two compared groups.
#' @param panels.row,panels.col Integer: number of rows and columns to arrange panels with plots.
#' @param text.size Numeric: size of text on the plot. Default to 1.1.
#' @param print.auc Logical: if TRUE, AUC values with confidence interval will be added to the plot. Default to TRUE.
#' @param print.auc.size Numeric: size of AUC text on the plot. Default to 0.8.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "ROC_plot".
#' @param save.to.txt Logical: if TRUE, returned table with results will be saved to .txt file. Default to FALSE.
#' @param name.txt Character: name of saved .txt file, without ".txt" name of extension. Default to "ROC_results".
#'
#' @return Data frame with ROC parameters including AUC, threshold, specificity, sensitivity, accuracy,
#' positive predictive value, negative predictive value, and Youden's J statistic.
#' Plot with ROC curves can be saved to .tiff file and opened from the working directory (will not be displayed on the graphic device).
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(pROC)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' data.dCtF <- filter_transformed_data(data.dCt, remove.Sample = c("Control11"))
#' roc_parameters <- ROCh(data.dCtF, sel.Gene = c("Gene1","Gene16","Gene19","Gene20"),
#'                         groups = c("Disease","Control"),
#'                         panels.row = 2,
#'                         panels.col = 2)
#'
#' @importFrom dplyr select filter
#' @importFrom utils write.table
#' @importFrom pROC roc coords plot.roc
#' @importFrom grDevices tiff dev.off
#' @importFrom graphics par
#' @import tidyverse
#'
ROCh <- function(data,
                 sel.Gene = "all",
                 groups,
                 panels.row,
                 panels.col,
                 text.size = 1.1,
                 print.auc = TRUE,
                 print.auc.size = 0.8,
                 save.to.tiff = FALSE,
                 dpi = 600,
                 width = 15,
                 height = 15,
                 name.tiff = "ROC_plot",
                 save.to.txt = FALSE,
                 name.txt = "ROC_results") {

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  data <- filter(data, Group %in% groups)

  if (sel.Gene[1] != "all") {
    data <- data[, colnames(data) %in% c("Group", "Sample", sel.Gene)]
  }

  roc_param <- as.data.frame(matrix(nrow = ncol(data) - 2, ncol = 9))
  colnames(roc_param) <-
    c(
      "Gene",
      "Threshold",
      "Specificity",
      "Sensitivity",
      "Accuracy",
      "ppv",
      "npv",
      "youden",
      "AUC"
    )
  roc_param$Gene <- colnames(data)[-c(1:2)]

  for (x in 1:nrow(roc_param)) {
    myproc <- roc(
      response = data$Group,
      predictor = as.data.frame(data)[, x + 2],
      levels = c(groups),
      smooth = FALSE,
      auc = TRUE,
      plot = FALSE,
      ci = TRUE,
      of = "auc",
      quiet = TRUE
    )

    parameters <- coords(
      myproc,
      "best",
      ret = c(
        "threshold",
        "specificity",
        "sensitivity",
        "accuracy",
        "ppv",
        "npv",
        "youden"
      )
    )
    roc_param[x, 2:8] <- parameters
    roc_param[x, 9] <- myproc$auc
    roc_param[x, 1] <- colnames(data)[x + 2]

    if (nrow(parameters) > 1) {
      message(
        '\nImportant: ',
        colnames(data)[x + 2],
        ' has more than 1 threshold value for calculated Youden J statistic.\n'
      )
    } else {
    }
  }
  if (save.to.tiff == TRUE) {
    tiff(
      paste(name.tiff, ".tiff", sep = ""),
      res = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )

    par(mfrow = c(panels.row, panels.col))

    for (x in 1:nrow(roc_param)) {
      myproc <- roc(
        response = data$Group,
        predictor = as.data.frame(data)[, x + 2],
        levels = c(groups),
        smooth = FALSE,
        auc = TRUE,
        plot = FALSE,
        ci = TRUE,
        of = "auc",
        quiet = TRUE
      )

      plot.roc(
        myproc,
        main = roc_param$Gene[x],
        smooth = FALSE,
        cex.axis = text.size,
        cex.lab = text.size,
        identity.lwd = 2,
        plot = TRUE,
        percent = TRUE,
        print.auc = print.auc,
        print.auc.x = 0.85,
        print.auc.y = 0.1,
        print.auc.cex = print.auc.size
      )
    }
    dev.off()
  }

  if (save.to.txt == TRUE) {
    write.table(as.data.frame(roc_param), paste(name.txt, ".txt", sep = ""))
  }

  return(roc_param)
}






#' @title log_reg
#'
#' @description
#' This function performs logistic regression analysis, computes odd ratio values, and presents them graphically.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all genes.
#' @param group.study Character: name of study group (group of interest).
#' @param group.ref Character: name of reference group.
#' @param increment Numeric or "mean": the change of expression for which odds ratio values are calculated.
#' If increment = 1, odds ratio values regards for a one-unit increase in gene expression
#' (more suitable where data are not transformed using 2^-value formula). If increment = "mean", odds ratio values are
#' calculated for the situation where gene expression increases by mean of gene expression levels in all samples
#' (more suitable where data were transformed using 2^-value formula).
#' @param centerline Numeric: position of vertical centerline on the plot. Default to 1.
#' @param ci Numeric: confidence level used for computation of confidence interval. Default to 0.95.
#' @param log.axis Logical: if TRUE, axis with odds ratio values will be in log10 scale. Default to FALSE.
#' @param x.axis.title Character: title of x axis. Default to "Gene".
#' @param y.axis.title Character: title of y axis. Default to "value".
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param legend.title Character: title of legend. Default to "Group".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be one of "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "OR_plot".
#' @param save.to.txt Logical: if TRUE, returned table with results will be saved to .txt file. Default to FALSE.
#' @param name.txt Character: name of saved .txt file, without ".txt" name of extension. Default to "OR_results".
#'
#' @return A list that contains an object with plot and data frame with calculated results. Created plot is also displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(oddsratio)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                        remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready,
#'                      ref = "Gene8")
#' log.reg.results <- log_reg(data.dCt,
#'                            sel.Gene = c("Gene1","Gene16","Gene19","Gene20"),
#'                            group.study = "Disease",
#'                            group.ref = "Control",
#'                            increment = 1)
#'
#' @importFrom dplyr filter mutate ungroup
#' @importFrom stats coef glm setNames
#' @importFrom utils write.table
#' @importFrom ggplot2 ggplot geom_vline geom_point geom_errorbarh scale_color_continuous geom_text coord_flip scale_x_log10 guides xlab ylab labs theme_classic theme element_text ggsave scale_y_continuous expansion facet_wrap element_blank stat_summary
#' @import ggplot2
#' @import tidyverse
#' @import oddsratio
#'
log_reg <- function(data,
                    sel.Gene = "all",
                    group.study,
                    group.ref,
                    increment,
                    centerline = 1,
                    ci = 0.95,
                    log.axis = FALSE,
                    x.axis.title = "Odds ratio",
                    y.axis.title = "",
                    axis.title.size = 11,
                    axis.text.size = 10,
                    legend.title = "p value",
                    legend.text.size = 11,
                    legend.title.size = 11,
                    legend.position = "right",
                    plot.title = "",
                    plot.title.size = 14,
                    save.to.tiff = FALSE,
                    dpi = 600,
                    width = 15,
                    height = 15,
                    name.tiff = "OR_plot",
                    save.to.txt = FALSE,
                    name.txt = "OR_results") {
  data <- filter(data, Group %in% c(group.study, group.ref))

  if (sel.Gene[1] != "all") {
    data <- data[, colnames(data) %in% c("Group", "Sample", sel.Gene)]
    genes <-
      colnames(data)[!colnames(data) %in% c("Group", "Sample")]

  } else {
    genes <- colnames(data)[!colnames(data) %in% c("Group", "Sample")]
  }

  data <- mutate(data, Group_num = ifelse(Group == group.ref, 0, 1))
  n.genes <- ncol(data) - 3

  data.CI <- as.data.frame(matrix(ncol = 9, nrow = n.genes))
  colnames(data.CI) <-
    c(
      "Gene",
      "oddsratio",
      "CI_low",
      "CI_high",
      "Increment",
      "Intercept",
      "coeficient",
      "p_intercept",
      "p_coef"
    )


  if (increment == "mean") {
    for (x in 1:n.genes) {
      data.m <- data %>%
        ungroup() %>%
        select(all_of(c(genes[x], "Group_num")))
      m <- glm(data.m$Group_num ~ ., data = data.m, family = binomial)
      l <-
        setNames(as.list(mean(as.data.frame(data.m[, genes[x]])[, 1])), genes[x])
      or <- or_glm(
        data = data.m,
        model = m,
        incr = l,
        ci = 0.95
      )

      data.CI[x, 1:5] <- or
      data.CI[x, 6:7] <- m$coefficients
      data.CI[x, 8:9] <- coef(summary(m))[, 4]
    }
  } else{
    for (x in 1:n.genes) {
      data.m <- data %>%
        ungroup() %>%
        select(all_of(c(genes[x], "Group_num")))
      m <-
        glm(data.m$Group_num ~ ., data = data.m, family = binomial)
      l <- setNames(as.list(increment), genes[x])
      or <- or_glm(
        data = data.m,
        model = m,
        incr = l,
        ci = 0.95
      )

      data.CI[x, 1:5] <- or
      data.CI[x, 6:7] <- m$coefficients
      data.CI[x, 8:9] <- coef(summary(m))[, 4]
    }

  }

  od_df <- data.frame(
    yAxis = 1:nrow(data.CI),
    boxOdds = data.CI$oddsratio,
    boxCILow = data.CI$CI_low,
    boxCIHigh = data.CI$CI_high,
    boxLabels = data.CI$Gene,
    p = data.CI$p_coef
  )

  odd.ratio <-
    ggplot(od_df, aes(x = boxOdds, y = boxLabels, label = boxOdds)) +
    geom_vline(aes(xintercept = centerline),
               linewidth = .25,
               linetype = "dashed") +
    geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow),
                   linewidth = .5,
                   height = .2) +
    geom_point(aes(color = p), size = 3.5) +
    scale_color_continuous(type = "viridis") +
    geom_text(
      aes(label = boxOdds),
      hjust = 0.5,
      vjust = -1,
      size = 3
    ) +
    xlab(x.axis.title) + ylab(y.axis.title) +
    labs(color = legend.title, title = plot.title) +
    theme_bw() +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
    theme(legend.text = element_text(size = legend.text.size, colour = "black")) +
    theme(legend.title = element_text(size = legend.title.size, colour =
                                        "black")) +
    theme(plot.title = element_text(size = plot.title.size))

  if (log.axis == TRUE) {
    odd.ratio <- odd.ratio +
      scale_x_log10()
  }

  print(odd.ratio)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      odd.ratio,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }

  if (save.to.txt == TRUE) {
    write.table(as.data.frame(data.CI), paste(name.txt, ".txt", sep = ""))
  }

  return(list(odd.ratio, data.CI))
}









#' @title results_heatmap
#'
#' @description
#' This function creatse heatmap with hierarchical clustering.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all genes.
#' @param dist.row,dist.col Character: name of method used for calculation of distances between rows or columns, derived from stats::dist() function,
#' must be one of "euclidean" (default) , "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param clust.method Character: name of used method for agglomeration, derived from stats::hclust() function,
#' must be one of "ward.D", "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".
#' @param colors Vector with colors used to fill created heatmap.
#' @param col.groups A named list with colors for groups annotation (see example).
#' @param show.colnames,show.rownames Logical: of TRUE, names of columns (sample names) and rows (gene names) will be shown. Both default to TRUE.
#' @param border.color Character: color of cell borders on heatmap. If set to NA (default) no border will be drawn.
#' @param fontsize Numeric: global fontsize of heatmap. Default to 10.
#' @param fontsize.col,fontsize.row Numeric: fontsize of colnames and rownames. Default to 10.
#' @param angle.col Integer: angle of the column labels, one of the 0, 45, 90, 270, and 315.
#' @param cellwidth,cellheight Numeric: width and height of individual cell. Both default to NA.
#' These parameters are useful in situations where margins are too small and the plot is cropped (column names and annotation legend are sometimes partially hidden).
#' Specification of this parameter allows to adjust size of the plot and solve this problem.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "heatmap_results".
#'
#' @return Heatmap with hierarchical clustering, displayed on the graphic device (if save.to.tiff = FALSE)
#' or saved to .tiff file (if save.to.tiff = TRUE).
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(pheatmap)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                        remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' # Remember to firstly create named list with colors for groups annotation:
#' colors.for.groups = list("Group" = c("Disease"="firebrick1","Control"="green3"))
#' # Vector of colors to fill the heatmap can be also specified to fit the user needings:
#' colors <- c("navy","navy","#313695","#4575B4","#74ADD1","#ABD9E9",
#'             "#E0F3F8","#FFFFBF","#FEE090","#FDAE61","#F46D43",
#'             "#D73027","#C32B23","#A50026","#8B0000",
#'             "#7E0202","#000000")
#' results_heatmap(data.dCt,
#'                 sel.Gene = "all",
#'                 col.groups = colors.for.groups,
#'                 colors = colors,
#'                 show.colnames = FALSE,
#'                 show.rownames = TRUE,
#'                 fontsize = 11,
#'                 fontsize.row = 11)
#'
#' @importFrom stats hclust dist
#' @importFrom grDevices tiff colorRampPalette
#' @importFrom dplyr select ungroup
#' @importFrom tidyselect any_of
#' @importFrom pheatmap pheatmap
#' @import tidyverse
#'
results_heatmap <- function(data,
                            sel.Gene = "all",
                            dist.row = "euclidean",
                            dist.col = "euclidean",
                            clust.method = "average",
                            col.groups,
                            colors = c(
                              "navy",
                              "#313695",
                              "#4575B4",
                              "#74ADD1",
                              "#ABD9E9",
                              "#E0F3F8",
                              "#FFFFBF",
                              "#FEE090",
                              "#FDAE61",
                              "#F46D43",
                              "#D73027",
                              "#C32B23",
                              "#A50026",
                              "#8B0000",
                              "#7E0202",
                              "#000000"
                            ),
                            show.colnames = TRUE,
                            show.rownames = TRUE,
                            border.color = NA,
                            fontsize = 10,
                            fontsize.col = 10,
                            fontsize.row = 10,
                            angle.col = 0,
                            cellwidth = NA,
                            cellheight = NA,
                            save.to.tiff = FALSE,
                            dpi = 600,
                            width = 15,
                            height = 15,
                            name.tiff = "heatmap_results") {
  if (sel.Gene[1] == "all") {
    data <- data

  } else {
    data <- select(data, Group, Sample, any_of(sel.Gene))
  }

  data <- as.data.frame(data)
  labels_col <- data$Sample
  data <- unite(data, Rownames, Sample, Group, remove = FALSE)
  rownames(data) <- data$Rownames
  mat <-
    data[,!colnames(data) %in% c("Group", "Sample", "Rownames")]
  mat <- t(mat)
  mat <- mat - rowMeans(mat)
  df <- as.data.frame(data$Group)
  rownames(df) <- colnames(mat)
  colnames(df) <- c("Group")
  colors.to.fill <- colorRampPalette(colors)(255)

  if (save.to.tiff == TRUE) {
    device <- .Device

    if (device == "RStudioGD") {
      dev.off()
    }

    tiff(
      paste(name.tiff, ".tiff", sep = ""),
      res = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )

    pheatmap(
      mat,
      annotation_col = df,
      annotation_colors = col.groups,
      clustering_method = clust.method,
      clustering_distance_cols = dist.col,
      clustering_distance_rows = dist.row,
      show_colnames = show.colnames,
      show_rownames = show.rownames,
      border_color = border.color,
      fontsize = fontsize,
      fontsize_col = fontsize.col,
      fontsize_row = fontsize.row,
      color = colors.to.fill,
      angle_col = angle.col,
      cellwidth = cellwidth,
      cellheight = cellheight
    )
    dev.off()

  } else {
    pheatmap(
      mat,
      annotation_col = df,
      annotation_colors = col.groups,
      clustering_method = clust.method,
      clustering_distance_cols = dist.col,
      clustering_distance_rows = dist.row,
      show_colnames = show.colnames,
      show_rownames = show.rownames,
      border_color = border.color,
      fontsize = fontsize,
      fontsize_col = fontsize.col,
      fontsize_row = fontsize.row,
      color = colors.to.fill,
      angle_col = angle.col,
      cellwidth = cellwidth,
      cellheight = cellheight,
      labels_col = labels_col
    )
  }
}






#' @title results_volcano
#'
#' @description
#' This function creates volcano plot that illustrate the arrangement of genes based on
#' fold change values and p values. Significant genes can be pointed out using
#' specified p value and fold change thresholds, and highlighted on the plot by color and (optionally) isolated by thresholds lines.
#'
#' @param data Object returned from RQ_dCt() or RQ_ddCt() functions.
#' @param mode Character: which p value should be used? One of the "t" (p values from Student's t test), "t.adj" (adjusted p values from Student's t test),
#' "mw" (p values from Mann-Whitney U test),"mw.adj" (adjusted p values from Mann-Whitney U test),
#' "depends" (if data in both compared groups were considered as derived from normal distribution
#' (p value from Shapiro_Wilk test > 0.05) - p values from Student's t test will be used for significance assignment,
#' otherwise p values from Mann-Whitney U test will be used), "depends.adj"
#' (if data in both compared groups were considered as derived from normal distribution
#' (p value from Shapiro_Wilk test > 0.05) - adjusted p values from Student's t test will be used for significance assignment,
#' otherwise adjusted p values from Mann-Whitney U test will be used), and "user" that can be used
#' the user intend to use another p values, e.g. obtained from other statistical test. In such situation, before run FCh_plot function, the user should prepare
#' data frame object named "user" that contains two columns, the first of them with Gene names and the second with p values.
#' The order of columns must be kept as described.
#' @param p.threshold Numeric: threshold of p values for statistical significance.
#' @param FCh.threshold Numeric: threshold of fold change values used for significance assignment of genes.
#' If is set to 2 (default), genes with 2-fold changed expression (increased or decreased) between groups will be assigned as significant.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all names of genes.
#' @param point.size Numeric: size of points. Default to 4.
#' @param point.shape Integer: shape of points. Default to 19..
#' @param alpha Numeric: transparency of points, a value between 0 and 1. Default to 0.7.
#' @param colors Character vector length of two, containing colors for significant and no significant genes.
#' @param add.thr.lines Logical, if TRUE (default), threshold lines will be added to the plot.
#' @param linewidth Numeric: width of the added threshold lines. Default to 0.25.
#' @param linetype Character: type of the added threshold lines. One of the "solid", "dashed" (default), "dotted",
#' "dotdash", "longdash", "twodash", and "blank".
#' @param x.axis.title Character: title of x axis. Default to "Gene".
#' @param y.axis.title Character: title of y axis. Default to "value".
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param legend.title Character: title of legend. Default to "Group".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be one of "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff Character: name of saved .tiff file, without ".tiff" name of extension. Default to "FCh_plot".
#'
#' @return List containing object with barplot and data frame with results. Created plot is also displayed on graphic device.
#' @export
#'
#' @examples
#' library(ggsignif)
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                        remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                        remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' data.dCtF <- filter_transformed_data(data.dCt, remove.Sample = c("Control11"))
#' results.dCt <- RQ_dCt(data.dCtF, "Disease", "Control")
#'
#' RQ.volcano <- results_volcano(data = results.dCt,
#'                          mode = "depends",
#'                          p.threshold = 0.05,
#'                          FCh.threshold = 2)
#' head(RQ.volcano[[2]])
#
#' @importFrom dplyr select filter group_by summarise mutate ungroup desc
#' @importFrom tidyr pivot_longer
#' @importFrom magrittr %>%
#' @importFrom ggplot2 ggplot geom_point scale_color_manual xlab ylab labs theme_classic theme element_text ggsave element_blank geom_vline geom_hline
#' @import ggplot2
#' @import tidyverse
#'
results_volcano <- function(data,
                       mode,
                       p.threshold = 0.05,
                       FCh.threshold = 2,
                       sel.Gene = "all",
                       point.size = 4,
                       point.shape = 19,
                       alpha = 0.7,
                       colors = c("#66c2a5", "#fc8d62"),
                       add.thr.lines = TRUE,
                       linewidth = 0.25,
                       linetype = "dashed",
                       x.axis.title = "log2(fold change)",
                       y.axis.title = "-log10(p value)",
                       axis.title.size = 11,
                       axis.text.size = 10,
                       legend.text.size = 11,
                       legend.title = "Selected as significant?",
                       legend.title.size = 11,
                       legend.position = "top",
                       plot.title = "",
                       plot.title.size = 14,
                       dpi = 600,
                       width = 15,
                       height = 15,
                       save.to.tiff = FALSE,
                       name.tiff = "results_volcano") {
  if (sel.Gene[1] != "all") {
    data <- filter(data, Gene %in% sel.Gene)
  }

  if (mode == "t") {
    data$p.used <- data$t_test_p
  }
  if (mode == "t.adj") {
    data$p.used <- data$t_test_p_adj
  }
  if (mode == "mw") {
    data$p.used <- data$MW_test_p
  }
  if (mode == "mw.adj") {
    data$p.used <- data$MW_test_p_adj
  }
  if (mode == "depends") {
    data <- ungroup(data)
    vars <- colnames(select(data, ends_with("norm_p")))
    data <-
      mutate(
        data,
        test.for.comparison = ifelse(
          .data[[vars[[1]]]] >= 0.05 &
            .data[[vars[[2]]]] >= 0.05,
          "t.student's.test",
          "Mann-Whitney.test"
        )
      )
    data <-
      mutate(data,
             p.used = ifelse(
               test.for.comparison == "t.student's.test",
               t_test_p,
               MW_test_p
             ))
  }
  if (mode == "depends.adj") {
    data <- ungroup(data)
    vars <- colnames(select(data, ends_with("norm_p")))
    data <-
      mutate(
        data,
        test.for.comparison = ifelse(
          .data[[vars[[1]]]] >= 0.05 &
            .data[[vars[[2]]]] >= 0.05,
          "t.student's.test",
          "Mann-Whitney.test"
        )
      )
    data <-
      mutate(
        data,
        p.used = ifelse(
          test.for.comparison == "t.student's.test",
          t_test_p_adj,
          MW_test_p_adj
        )
      )
  }
  if (mode == "user") {
    colnames(user) <- c("Gene", "p.used")
    data <- full_join(data, user, by = c("Gene"))
  }

  data <-
    mutate(data,
           `Selected as significant?` = ifelse(
             p.used > p.threshold,
             yes = "No",
             no = ifelse(abs(log2(FCh)) <  log2(FCh.threshold), "No", "Yes")
           ))
  data$`Selected as significant?` <-
    factor(data$`Selected as significant?`, levels = c("Yes", "No"))

  Vol <-
    ggplot(data,
           aes(
             x = log2(FCh),
             y = -log10(p.used),
             color = `Selected as significant?`
           )) +
    geom_point(size = point.size,
               shape = point.shape,
               alpha = alpha) +
    scale_color_manual(values = c(colors)) +
    labs(colour = legend.title, title = plot.title) +
    xlab(x.axis.title) +
    ylab(y.axis.title) +
    theme_bw() +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
    theme(legend.text = element_text(size = legend.text.size, colour = "black")) +
    theme(legend.title = element_text(size = legend.title.size, colour =
                                        "black")) +
    theme(plot.title = element_text(size = plot.title.size)) +
    theme(panel.grid.major.x = element_blank())

  if (add.thr.lines == TRUE) {
    Vol <- Vol +
      geom_vline(aes(xintercept = log2(FCh.threshold)),
                 linewidth = linewidth,
                 linetype = linetype) +
      geom_vline(aes(xintercept = -log2(FCh.threshold)),
                 linewidth = linewidth,
                 linetype = linetype) +
      geom_hline(aes(yintercept = -log10(p.threshold)),
                 linewidth = linewidth,
                 linetype = linetype)
  }

  print(Vol)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      Vol,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }
  return(list(Vol, data))
}







#' @title pca_kmeans
#'
#' @description
#' This function performs principal component analysis (PCA) together with k means analysis for samples,
#' and generate plot that illustrate spatial arrangement of samples based on the two first components and with assignment to k means clusters.
#' PCA analysis can not deal with missing values, thus all samples with at least one missing value are removed from data before analysis.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' @param do.k.means Logical: if TRUE (default), k means analysis will be performed.
#' @param k.clust Integer: number of clusters for k means analysis. Default to 2.
#' @param clust.names Character vector with names of clusters, must be equal to the number of clusters specified in the k.clust parameter.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all genes.
#' @param point.size Numeric: size of points. Default to 4.
#' @param point.shape Integer: shape of points. If do.k.means = TRUE, the number of provided values must be equal to the
#' number of cluster (k.clust). Default to c(19, 17).
#' @param alpha Numeric: transparency of points, a value between 0 and 1. Default to 0.7.
#' @param add.sample.labels Logical: if TRUE, points labels (names of samples) will be added. Default to FALSE.
#' @param label.size Numeric: size of points labels (names of samples). Default to 3.
#' @param hjust Numeric: horizontal position of points labels. Default to 0.
#' @param vjust Numeric: vertical position of points labels.  Default to -1.
#' @param point.color Character vector containing colors for compared groups. The number of colors must be equal to the number of groups. Default to c("#66c2a5", "#fc8d62").
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param legend.title.group Character: title of legend for groups. Default to "Group".
#' @param legend.title.cluster Character: title of legend for k means clusters. Default to "Clusters".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "pca_and_kmeans".
#'
#' @return A list containing object with plot and, if do.k.means = TRUE, a confusion matrix that show classification performance of k means method.
#' Created plot is also displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' pca_kmeans(data.dCt, sel.Gene = c("Gene1","Gene16","Gene19","Gene20"))
#'
#' @importFrom stats na.omit prcomp kmeans
#' @importFrom dplyr select
#' @importFrom tidyr unite
#' @importFrom ggplot2 ggplot geom_point geom_text scale_color_manual xlab ylab labs theme_classic theme element_text ggsave
#' @import ggplot2
#' @import tidyverse
#'
pca_kmeans <- function(data,
                       do.k.means = TRUE,
                       k.clust = 2,
                       clust.names = c("Cluster1", "Cluster2"),
                       sel.Gene = "all",
                       point.size = 4,
                       point.shape = c(19, 17),
                       alpha = 0.7,
                       point.color = c("#66c2a5", "#fc8d62"),
                       add.sample.labels = FALSE,
                       label.size = 3,
                       hjust = 0,
                       vjust = -1,
                       axis.title.size = 11,
                       axis.text.size = 10,
                       legend.text.size = 11,
                       legend.title.group = "Group",
                       legend.title.cluster = "Cluster",
                       legend.title.size = 11,
                       legend.position = "right",
                       plot.title = "",
                       plot.title.size = 14,
                       save.to.tiff = FALSE,
                       dpi = 600,
                       width = 15,
                       height = 15,
                       name.tiff = "pca_and_kmeans") {
  if (sel.Gene[1] == "all") {
    data <- data

  } else {
    data <- select(data, Group, Sample, any_of(sel.Gene))
  }

  data <- as.data.frame(data)
  data <- unite(data, Rownames, Sample, Group, remove = FALSE)
  rownames(data) <- data$Rownames
  data <- na.omit(data)
  pca <-
    prcomp(select(data,-Rownames,-Sample,-Group), scale = TRUE)
  var_pca1 <- summary(pca)$importance[2, ][1]
  var_pca2 <- summary(pca)$importance[2, ][2]
  pca_comp <- as.data.frame(pca$x)
  pca_comp$Sample <- data$Sample
  pca_comp$Group <- data$Group

  if (do.k.means == TRUE) {
    km <- kmeans(pca_comp[, 1:2], k.clust)
    pca_comp$k.means.cluster <- km$cluster
    pca_comp$k.means.cluster <- as.factor(pca_comp$k.means.cluster)
    levels(pca_comp$k.means.cluster) <- clust.names
    conf.matrix <- table(pca_comp$Group, pca_comp$k.means.cluster)
    pca.plot <-
      ggplot(pca_comp,
             aes(
               x = PC1,
               y = PC2,
               label = Sample,
               color = Group,
               shape = k.means.cluster
             )) +
      geom_point(size = point.size, alpha = alpha) +
      scale_shape_manual(values = c(point.shape)) +
      labs(shape = legend.title.cluster)

  } else {
    pca.plot <-
      ggplot(pca_comp, aes(
        x = PC1,
        y = PC2,
        label = Sample,
        color = Group
      )) +
      geom_point(size = point.size,
                 shape = point.shape,
                 alpha = alpha)
  }


  pca.plot <- pca.plot +
    scale_color_manual(values = c(point.color)) +
    labs(colour = legend.title.group, title = plot.title) +
    theme_bw() +
    labs(
      x = paste("PC1: ", round(var_pca1 * 100, 2), "% variance explained", sep = ""),
      y = paste("PC2: ", round(var_pca2 * 100, 2), "% variance explained", sep = "")
    ) +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
    theme(legend.text = element_text(size = legend.text.size, colour = "black")) +
    theme(legend.title = element_text(size = legend.title.size, colour =
                                        "black")) +
    theme(plot.title = element_text(size = plot.title.size))



  if (add.sample.labels == TRUE) {
    pca.plot <- pca.plot +
      geom_text(
        aes(label = Sample),
        hjust = hjust,
        vjust = vjust,
        size = label.size
      )
  }

  print(pca.plot)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      pca.plot,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }

  if (do.k.means == TRUE) {
    return(list(pca.plot, conf.matrix))

  } else{
    return(pca.plot)
  }
}







#' @title parallel_plot
#'
#' @description
#' This function illustrates expression values in a pairwise samples as series of lines connected across each axis.
#' This function can be used only for a pairwise data.
#'
#' @param data Object returned from make_Ct_ready() or delta_Ct() functions.
#' @param sel.Gene Character vector with names of genes to include, or "all" (default) to use all genes.
#' @param scale Character: a scale used for data presentation, one of the passed to ggparcoord() function.
#' Generally, scaling is not required since variables are the same units. Default to "globalminmax (no scaling).
#' @param order Character: method for groups ordering, one of the used in the ggparcoord() function. Default to 'anyClass'.
#' Must either be a vector of column indices (obligatory if only one gene is plotted), starting from 3 (e.g., for two groups it can be c(3,4) or c(4,3)),
#' or one of 'skewness', 'allClass', 'anyClass' (default), as well as scagnostic measures available in the 'scagnostics'
#' package (must be loaded): 'Outlying', `Skewed', 'Clumpy', 'Sparse', 'Striated', 'Convex', 'Skinny', 'Stringy', 'Monotonic'.
#' @param alpha Numeric: transparency of lines, a value between 0 and 1. Default to 0.7.
#' @param custom.colors Logical: if custom vector colors for genes is provided (and passed to the colors parameter),
#' it should be set to TRUE. For default colors, use custom.colors = FALSE (default).
#' @param colors Character vector containing custom colors for genes.
#' The number of colors must be equal to the number of presented genes. Must be provided if custom.colors = TRUE.
#' @param linewidth Numeric: width of lines. Default to 1.
#' @param show.points Logical: if TRUE (default), points will be also shown.
#' @param x.axis.title Character: title of x axis. Default to "".
#' @param y.axis.title Character: title of y axis. Default to "value".
#' @param axis.title.size Integer: font size of axis titles. Default to 11.
#' @param axis.text.size Integer: font size of axis text. Default to 10.
#' @param legend.title Character: title of legend for groups. Default to "Gene".
#' @param legend.title.size Integer: font size of legend title. Default to 11.
#' @param legend.text.size Integer: font size of legend text. Default to 11.
#' @param legend.position Position of the legend, can be "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param plot.title Character: title of plot. Default to "".
#' @param plot.title.size Integer: font size of plot title. Default to 14.
#' @param save.to.tiff Logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi Integer: resolution of saved .tiff file. Default to 600.
#' @param width Numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height Numeric: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "pca_and_kmeans".
#'
#' @return An object with plot. Created plot is also displayed on the graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(GGally)
#' data(data.Ct.pairwise)
#' data.CtF.pairwise <- filter_Ct(data = data.Ct.pairwise,
#'                                flag.Ct = "Undetermined",
#'                                maxCt = 35,
#'                                flag = c("Undetermined"),
#'                                remove.Gene = c("Gene9", "Gene2","Gene5", "Gene11","Gene1"))
#' data.CtF.ready.pairwise <- make_Ct_ready(data = data.CtF.pairwise,
#'                                          imput.by.mean.within.groups = TRUE)
#' data.dCt.pairwise <- delta_Ct(data = data.CtF.ready.pairwise,
#'                               ref = "Gene4")
#' parallel.plot <- parallel_plot(data = data.dCt.pairwise,
#'                                sel.Gene = c("Gene8","Gene19"))
#'
#' @importFrom stats na.omit prcomp kmeans
#' @importFrom dplyr select
#' @importFrom tidyr unite
#' @importFrom ggplot2 ggplot geom_point geom_text scale_color_manual xlab ylab labs theme_classic theme element_text ggsave
#' @import ggplot2
#' @import tidyverse
#' @import GGally
#'
parallel_plot <- function(data,
                          sel.Gene = "all",
                          scale = "globalminmax",
                          alpha = 0.7,
                          custom.colors = FALSE,
                          order = "anyClass",
                          colors,
                          linewidth = 1,
                          show.points = TRUE,
                          x.axis.title = "",
                          y.axis.title = "value",
                          axis.title.size = 11,
                          axis.text.size = 10,
                          legend.text.size = 11,
                          legend.title = "Gene",
                          legend.title.size = 11,
                          legend.position = "top",
                          plot.title = "",
                          plot.title.size = 14,
                          save.to.tiff = FALSE,
                          dpi = 600,
                          width = 15,
                          height = 15,
                          name.tiff = "parallel_plot") {
  if (sel.Gene[1] == "all") {
    data <- data

  } else {
    data <- select(data, Group, Sample, any_of(sel.Gene))
  }

  data <- as.data.frame(data)
  data <-
    pivot_longer(data,
                 !c(Sample, Group),
                 names_to = "Gene" ,
                 values_to = "value")
  data <- pivot_wider(data, names_from = Group, values_from = value)
  data$Gene <- as.factor(data$Gene)

  lines <- ggparcoord(
    data,
    columns = 3:length(colnames(data)),
    groupColumn = 2,
    order = order,
    scale = scale,
    showPoints = show.points,
    title = plot.title,
    alphaLines = alpha,
    mapping = ggplot2::aes(linewidth = linewidth)
  ) +
    ggplot2::scale_linewidth_identity()

  if (custom.colors == TRUE) {
    lines <- lines +
      scale_color_manual(values = colors)
  }

  lines <- lines +
    labs(colour = legend.title) +
    theme_bw() +
    xlab(x.axis.title) +
    ylab(y.axis.title) +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour = "black")) +
    theme(legend.text = element_text(size = legend.text.size, colour = "black")) +
    theme(legend.title = element_text(size = legend.title.size, colour = "black")) +
    theme(plot.title = element_text(size = plot.title.size))

  print(lines)

  if (save.to.tiff == TRUE) {
    ggsave(
      paste(name.tiff, ".tiff", sep = ""),
      lines,
      dpi = dpi,
      width = width,
      height = height,
      units = "cm",
      compression = "lzw"
    )
  }
  return(lines)
}






if (getRversion() >= "2.15.1")
  utils::globalVariables(
    c(
      "Group",
      "Gene",
      "Sample",
      "Ct",
      "dCt",
      "ddCt",
      "dCt_sd",
      "shap_wilka_p",
      "value",
      "value_sd",
      "test.for.comparison",
      "t_test_p",
      "MW_test_p",
      "p.used",
      "FCh",
      "Selected as significant?",
      "x",
      "y",
      "xend",
      "annotation",
      "Flag",
      "Var1",
      "Freq",
      "Var2",
      "Var3",
      "No",
      "Yes",
      "PC1",
      "PC2",
      "cor",
      "binomial",
      "boxOdds",
      "boxLabels",
      "boxCIHigh",
      "boxCILow",
      "p",
      "k.means.cluster",
      "t_test_p_adj",
      "MW_test_p_adj",
      "log10FCh",
      "FCh_sd",
      "Rownames",
      "calibrator"
    )
  )
