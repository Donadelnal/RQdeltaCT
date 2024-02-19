
#' @title read_Ct_wide
#'
#' @description
#' Enables to import Ct dataset in a wide-format table with sample names given in columns.
#'
#' @details
#' This function needs two files to import: a wide-format table with Ct values and design file (see parameters path.Ct.file and path.design.file for further details regarding tables structure).
#' Subsequently merges both files to return long-format table ready for analysis.
#' All parameters must be specified, there is no default values.
#'
#' @param path.Ct.file path to wide-format table in .txt format containing Ct values, gene names in the first column, and
#' sample names in the first row. In other words, this table should contain genes by rows and samples by columns.
#' @param path.design.file path to .txt file with two columns: column named "Sample" with names of samples
#' and column named "Group" with names of groups assigned to samples. Names of samples in this file
#' should correspond to the names of columns in file with Ct values.
#' @param sep character of a field separator in both imported files.
#' @param dec character used for decimal points in Ct values.
#'
#' @return Data frame in long format ready to analysis.
#' @export
#'
#' @examples
#' path.Ct.file <- system.file("extdata", "data_Ct_wide.txt", package = "RQdeltaCT")
#' path.design.file <- system.file("extdata", "data_design.txt", package = "RQdeltaCT")
#'
#' library(tidyverse)
#' data.Ct <- read_Ct_wide(path.Ct.file = path.Ct.file,
#'                    path.design.file = path.design.file,
#'                    sep ="\t",
#'                    dec = ".")
#' str(data.Ct)
#'
#' @importFrom utils read.csv
#' @importFrom tidyr pivot_longer
#' @importFrom base which
#' @import tidyverse
#'
read_Ct_wide <- function(path.Ct.file,
                         path.design.file,
                         sep,
                         dec){

  data_wide <- read.csv(path.Ct.file,
                        header = TRUE,
                        sep = sep,
                        dec = dec)

  data_wide_design <- read.csv(path.design.file,
                               header = TRUE,
                               sep = sep)

  colnames(data_wide)[1] <- "Gene"
  data_wide <- mutate(data_wide, across(everything(), as.character))
  data_slim <- pivot_longer(data_wide, -Gene, names_to = "Sample", values_to = "Ct")
  data_slim[ ,"Group"] <- NA

  for (x in 1:nrow(data_wide_design)) {
    index <- which(data_slim$Sample == data_wide_design$Sample[x])
    data_slim$Group[index] <- data_wide_design$Group[x]
  }
  return(data_slim)
}





#' @title read_Ct_long
#'
#' @description
#' Imports long-format table with Ct values.
#'
#' @param path path to a .txt file with long-type table of Ct values. This table should contain at least  4 columns, with
#' sample names, gene names, Ct values and group names (those columns will be imported by this function).
#' Imported table could also contain a column with flag information, which could be optionally imported (see add.col.Flag and col.Flag parameters).
#'
#' @param sep character of a field separator in imported file.
#' @param dec character used for decimal points in Ct values.
#' @param skip integer: number of lines of the data file to skip before beginning to read data. Default to 0.
#' @param col.Sample integer: number of column with sample names.
#' @param col.Gene integer: number of column with gene names.
#' @param col.Ct integer: number of column with Ct values.
#' @param col.Group integer: number of column with group names.
#' @param add.col.Flag logical: if data contains a column with flag information which should be also imported, this parameters should be set to TRUE. Default to FALSE.
#' @param col.Flag integer: number of column with flag information. Should be specified if add.col.Flag = TRUE.
#' This column should contain a character-type values (ex. "Undetermined" and "OK"), however,
#' other types of values are allowed (ex. numeric), but must be converted to character or factor after importing data (see examples).
#'
#' @return Data.frame in long format ready to analysis.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "data_Ct_long.txt", package = "RQdeltaCT")
#'
#' library(tidyverse)
#' data.Ct <- read_Ct_long(path = path, sep = "\t",dec = ".",skip = 0,
#'                         add.column.Flag = TRUE, column.Sample = 1, column.Gene = 2,
#'                         column.Ct = 5, column.Group = 9, column.Flag = 4)
#' str(data.Ct)
#'
#' data.Ct <- mutate(data.Ct, Flag = ifelse(Flag < 1, "Undetermined", "OK"))
#' str(data.Ct)
#'
#' @importFrom utils read.csv
#' @importFrom base colnames
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
                         column.Flag){

  data <- read.csv(path,
                   header = TRUE,
                   sep = sep,
                   dec = dec,
                   skip = skip)

  if (add.column.Flag == FALSE){
  data <- data[ ,c(column.Sample, column.Gene, column.Ct, column.Group)]
  colnames(data) <- c("Sample", "Gene", "Ct", "Group")
  }

  if (add.column.Flag == TRUE){
    data <- data[ ,c(column.Sample, column.Gene, column.Ct, column.Group, column.Flag)]
    colnames(data) <- c("Sample", "Gene", "Ct", "Group", "Flag")
  }

  return(data)
}







#' @title control_Ct_barplot_sample
#'
#' @description
#' Sample-wide control of raw Ct values by illustrating numbers of Ct values labeled as reliable or not by using reliability criteria (see function parameters).
#'
#' @details
#' This function does not perform data filtering, but only numbers Ct values labeled as reliable or not and presents them graphically.
#' Results could be useful to identify samples with low number of reliable Ct values.
#'
#' @param data object returned from read_Ct_long() or read_Ct_wide() function,
#' or data frame containing column named "Sample" with sample names, column named "Gene" with gene names,
#' column named "Ct" with raw Ct values, and column named "Group" with group names.
#' Optionally, data frame could contain column named "Flag" with flag information (ex. "Undetermined" and "OK"), which will be used for reliability assessment.
#' @param flag.Ct character of a flag used for undetermined Ct values. Default to "Undetermined".
#' @param maxCt numeric, a maximum of Ct value allowed. Default to 35.
#' @param flag character of a flag used in Flag column for values which are unreliable. Default to "Undetermined".
#' @param colors character vector length of two, containing colors for Ct values which were labeled as reliable (first element of vector) or not (second element of vector).
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 10.
#' @param x.axis.title character: title of x axis. Default to "".
#' @param y.axis.title character: title of y axis. Default to "Number".
#' @param legend.title character: title of legend. Default to "Reliable Ct value?".
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param legend.title.size integer: font size of legend title.  Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file.  Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "Ct_control_barplot_for_samples".
#'
#' @return List containing plot and table with numbers of reliable and not reliable Ct values in samples.
#' Additional information about returned table is also printed, it could help user to properly interpret returned table.
#' Plot will be displayed on graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' sample.Ct.control <- control_Ct_barplot_sample(data.Ct)
#' sample.Ct.control[[2]]
#'
#' @importFrom base as.numeric as.data.frame cat print table paste sum colnames list
#' @importFrom dplyr mutate arrange filter
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
                                      save.to.tiff = FALSE, dpi = 600, width = 15, height = 15,
                                      name.tiff = "Ct_control_barplot_for_samples"){

  data$Ct[data$Ct == flag.Ct] <- 100
  data$Ct <- as.numeric(data$Ct)

  if(sum(colnames(data) %in% "Flag") > 0){
    data <- mutate(data, Reliable = ifelse(Ct > maxCt | Flag == flag, yes = "No",  no = "Yes"))

  } else {
  data <- mutate(data, Reliable = ifelse(Ct > maxCt , yes = "No",  no = "Yes"))
  }

  bar <- as.data.frame(table(data$Reliable, data$Sample))
  order <- arrange(filter(bar, Var1 == "Yes"), Freq)$Var2
  cat("Returned table contains numbers of Ct values labeled as reliable or not in each sample, as well as fraction of unreliable Ct values in each sample.\n")

  barplot.samples <- ggplot(bar, aes(x = reorder(Var2, desc(Freq)), y = Freq, fill = Var1)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(breaks = c("Yes", "No"), values = c("Yes" = colors[1],"No" = colors[2])) +
    xlab(x.axis.title) +
    ylab(y.axis.title) +
    labs(fill = legend.title, title = plot.title) +
    theme_classic() +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, color = 'black')) +
    theme(axis.title = element_text(size = axis.title.size, color = 'black')) +
    theme(legend.title = element_text(size = legend.title.size, colour="black")) +
    theme(legend.text = element_text(size = legend.text.size, colour="black")) +
  	theme(plot.title = element_text(size = plot.title.size)) +
    scale_x_discrete(limits = order)

  print(barplot.samples)

    if (save.to.tiff == TRUE){
      ggsave(paste(name.tiff, ".tiff", sep = ""), barplot.samples, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
    } else {}

  tab <- table(data$Reliable, data$Sample)
  tab <- tab %>%
    as.data.frame() %>%
    pivot_wider(names_from = Var1, values_from = Freq) %>%
    arrange(desc(No)) %>%
    mutate(Not.reliale.fraction = No/(No+Yes)) %>%
    rename(Sample = Var2, Not.reliable = No, Reliable = Yes)

  return(list(barplot.samples, tab))
}







#' @title control_Ct_barplot_gene
#'
#' @description
#' Gene-wide control of raw Ct values across groups by illustrating numbers of Ct values labeled as reliable or not by using reliability criteria (see function parameters).
#'
#' @details
#' This function does not perform data filtering, but only numbers Ct values labeled as reliable or not and presents them graphically.
#' Could be useful to identify genes with low number of reliable Ct values.
#'
#' @param data object returned from read_Ct_long() or read_Ct_wide() function,
#' or data frame containing column named "Sample" with sample names, column named "Gene" with gene names,
#' column named "Ct" with raw Ct values, column named "Group" with group names.
#' Optionally, data frame could contain column named "Flag" with flag information (ex. "Undetermined" and "OK"),
#' which will be used for reliability assessment.
#' @param flag.Ct character of a flag used for undetermined Ct values. Default to "Undetermined".
#' @param maxCt numeric, a maximum of Ct value allowed. Default to 35.
#' @param flag character of a flag used in Flag column for values which are unreliable. Default to "Undetermined".
#' @param colors character vector length of two, containing colors for Ct values which were labeled as reliable (first element of vector) or not (second element of vector).
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 10.
#' @param x.axis.title character: title of x axis. Default to "".
#' @param y.axis.title character: title of y axis. Default to "Number".
#' @param legend.title character: title of legend. Default to "Reliable Ct value?".
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param legend.title.size integer: font size of legend title.  Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file.  Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension.  Default to "Ct_control_barplot_for_genes".
#'
#' @return List containing plot and table with numbers of reliable and not reliable Ct values in genes.
#' Additional information about returned table is also printed, it could help user to properly interpret returned table.
#' Plot will be displayed on graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' gene.Ct.control <- control_Ct_barplot_gene(data.Ct)
#' gene.Ct.control[[2]]
#'
#' @importFrom base as.numeric as.data.frame cat print table paste sum colnames
#' @importFrom dplyr mutate arrange filter
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
                                      dpi = 600, width = 15, height = 15,
                                      name.tiff = "Ct_control_barplot_for_genes"){

  data$Ct[data$Ct == flag.Ct] <- 100
  data$Ct <- as.numeric(data$Ct)

  if(sum(colnames(data) %in% "Flag") > 0){
    data <- mutate(data, Reliable = ifelse(Ct > maxCt | Flag == flag, yes = "No",  no = "Yes"))
  } else {
    data <- mutate(data, Reliable = ifelse(Ct > maxCt , yes = "No",  no = "Yes"))
  }

    bar <- as.data.frame(table(data$Reliable, data$Gene, data$Group))
    order <- arrange(filter(bar, Var1 == "Yes"), Freq)$Var2
    cat("Returned table contains numbers of Ct values labeled as reliable or not in each gene, as well as fraction of unreliable Ct values in each gene.\n")

    barplot.genes <- ggplot(bar, aes(x = reorder(Var2, desc(Freq)), y = Freq, fill = Var1)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(breaks = c("Yes", "No"), values = c("Yes" = colors[1],"No" = colors[2])) +
      xlab(x.axis.title) + ylab(y.axis.title) +
      labs(fill = legend.title, title = plot.title) +
      theme_classic() + theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, color = 'black')) +
      theme(axis.title = element_text(size = axis.title.size, color = 'black')) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
	    theme(plot.title = element_text(size = plot.title.size)) +
      scale_x_discrete(limits = rev(unique(bar$Var2))) +
      facet_wrap(vars(Var3))

    print(barplot.genes)

    if (save.to.tiff == TRUE){
      ggsave(paste(name.tiff, ".tiff", sep = ""), barplot.genes, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
    }

    tab <- table(data$Reliable, data$Gene, data$Group)
    tab <- tab %>%
      as.data.frame() %>%
      pivot_wider(names_from = Var1, values_from = Freq) %>%
      arrange(desc(No)) %>%
      mutate(Not.reliable.fraction = No/(No+Yes)) %>%
      rename(Gene = Var2, Group = Var3, Not.reliable = No, Reliable = Yes)

    return(list(barplot.genes, tab))
}






#' @title filter_Ct
#'
#' @description
#' Filters Ct data according to the used filtering criteria (see parameters).
#'
#' @param data object returned from read_Ct_long() or read_Ct_wide() function,
#' or data frame containing column named "Sample" with sample names, column named "Gene" with gene names,
#' column named "Ct" with raw Ct values, column named "Group" with group names.
#' Optionally, data frame could contain column named "Flag" with flag information (ex. "Undetermined" and "OK"),
#' which will be used for filtering.
#'
#' @param flag.Ct character of a flag used for undetermined Ct values, default to "Undetermined".
#' @param maxCt numeric, a maximum of Ct value allowed.
#' @param flag character: flag used in Flag column for values which should be filtered out, default to "Undetermined".
#' @param remove.Gene character: vector with names of genes which should be removed from data
#' @param remove.Sample character: vector with names of samples which should be removed from data
#' @param remove.Group character: vector with names of groups which should be removed from data
#'
#' @return Data.frame with filtered data.
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
#' @importFrom base as.numeric sum colnames
#' @importFrom dplyr filter
#' @import tidyverse
#'
filter_Ct <- function(data,
                      flag.Ct = "Undetermined",
                      maxCt = 35,
                      flag = c("Undetermined"),
                      remove.Gene = c(""),
                      remove.Sample = c(""),
                      remove.Group = c("")){

  data <- filter(data, Ct != flag.Ct)
  data$Ct <- as.numeric(data$Ct)
  data <- filter(data, Ct <= maxCt,
                        !Gene %in% remove.Gene,
                        !Sample %in% remove.Sample,
                        !Group %in% remove.Group)

  if(sum(colnames(data) %in% "Flag") > 0){
    data <- filter(data, !Flag %in% flag)
  }

  return(data)
}




#' @title make_Ct_ready
#'
#' @description
#' This function collapses technical replicates (if present in data) by means counts and imputes missing data by means within groups (if so indicated).
#' These actions also prepare Ct data for control functions.
#'
#' @param data data object returned from read_Ct_long(), read_Ct_wide() or filter_Ct() function,
#' or data frame containing column named "Sample" with sample names, column named "Gene" with gene names,
#' column named "Ct" with raw Ct values (must be numeric), column named "Group" with group names. Any other columns could exist, but will not be used by this function.
#' @param imput.by.mean.within.groups logical: if TRUE, missing values will be imputed by means within groups. This parameter could influence results, thus to draw more user attention on this parameter, no default value was set.
#' @param save.to.txt logical: if TRUE, returned data will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension. Default to "Ct_ready".
#'
#' @return Data.frame with prepared data and printed information about number and percentage of missing values.
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
#' @importFrom base mean
#' @importFrom dplyr select summarise
#' @import tidyverse
#'
make_Ct_ready <- function(data,
                   imput.by.mean.within.groups,
                   save.to.txt = FALSE,
                   name.txt = "Ct_ready"){
  data <- data %>%
    group_by(Group, Gene, Sample) %>%
    summarise(mean = mean(Ct, na.rm = TRUE), .groups = "keep") %>%
    as.data.frame()

  data_wide <- data %>%
    select(Group, Sample, Gene, mean) %>%
    pivot_wider(names_from = Gene, values_from = mean)

  nas <- sum(is.na(data_wide))
  percentage <- sum(is.na(data_wide))/((ncol(data_wide)-2)*nrow(data_wide))

  if (imput.by.mean.within.groups == TRUE){
    data_wide_imp <- data_wide %>%
      group_by(Group) %>%
      mutate(across(where(is.numeric), ~ replace(., is.na(.), mean(., na.rm = TRUE))))

    cat("Data contained", nas, "missing values that constitute", round(percentage*100, 5), "percent of the total data.\n Missing values were imputed using means within compared groups.\n")

    if (save.to.txt == TRUE){
      write.table(as.data.frame(data_wide_imp), paste(name.txt,".txt", sep = ""))
    }
    return(data_wide_imp)

  } else {

    cat("Data contains", nas, "missing values that constitute", round(percentage*100, 5), "percent of the total data.")

    if (save.to.txt == TRUE){
      write.table(as.data.frame(data_wide), paste(name.txt,".txt", sep = ""))
    }
    return(data_wide)
  }
}





#' @title exp_Ct_dCt
#'
#' @description
#' This function exponentiates Ct or delta Ct (dCt) values by using formula 2^(-Ct) or 2^(-dCt), respectively.
#'
#' @param data data object returned from make_Ct_ready() or delta_Ct() functions.
#' @param save.to.txt logical: if TRUE, returned data will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension. Default to "data_exp_Ct_dCt".
#'
#' @return Data frame with exponentiated Ct or dCt values.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.Ct.exp <- exp_Ct_dCt(data.CtF.ready)
#' head(data.Ct.exp)
#'
#' @importFrom base as.data.frame paste
#' @importFrom utils write.table
#' @importFrom dplyr mutate_at
#' @import tidyverse
#'
exp_Ct_dCt <- function(data,
                    save.to.txt = FALSE,
                    name.txt = "data_exp_Ct_dCt"){

  exp <- function(x){
    x <- 2^-x
  }

  data_exp <- data %>%
    mutate_at(vars(-c("Group","Sample")), exp)

  if (save.to.txt == TRUE){
    write.table(as.data.frame(data_exp), paste(name.txt,".txt", sep = ""))
  }
  return(data_exp)
}




#' @title RQ_exp_Ct_dCt
#'
#' @description
#' Performs relative quantification of gene expression using 2^(-Ct) and 2^(-dCt) methods.
#'
#' @details
#' This function calculates:
#' * Means (returned in columns with "_mean" pattern) and standard deviations (returned in columns with "_sd" pattern) of exponentiated Ct or dCt values of analyzed genes across compared groups.
#' * Normality tests (Shapiro_Wilk test) of exponentiated Ct or dCt values of analyzed genes across compared groups and returned p values in columns with "_norm_p" pattern.
#' * Fold Change values (return in "FCh" column) together with log10 Fold change values (return in "log10FCh" column).
#'   Fold change values were calculated for each gene by dividing  mean of exponentiated Ct od dCt values in study group by mean of exponentiated Ct or dCt values in reference group.
#' * Statistical testing of differences in exponentiated Ct or dCt values between study group and reference group.
#'   Student's t test and Mann-Whitney U test are implemented and resulted statistics (in column with "_test_stat" pattern) and p values (in column with "_test_p" pattern) are returned.
#'
#' @param data data object returned from exp_Ct_dCt() function.
#' @param group.study character: name of study group (group of interest).
#' @param group.ref character: name of reference group.
#' @param do.tests logical: if TRUE, statistical significance of differences between compared groups will be calculated using Student's t test and Mann-Whitney U test. Default to TRUE.
#' @param save.to.txt logical: if TRUE, returned table with results will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension. Default to "RQ_exp_results".
#
#' @return Data frame with transformed Ct values and printed information about number and percentage of missing values.
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(coin)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF)
#' data.Ct.exp <- exp_Ct_dCt(data.CtF.ready)
#' RQ.Ct.exp <- RQ_exp_Ct_dCt(data.Ct.exp, group.study = "Disease", group.ref = "Control")
#' head(RQ.Ct.exp)
#'
#' @importFrom base as.data.frame as.factor mean
#' @importFrom stats sd shapiro.test t.test
#' @importFrom coin wilcox_test pvalue statistic
#' @importFrom utils write.table
#' @importFrom dplyr filter select rename_with full_join
#' @import tidyverse
#'
RQ_exp_Ct_dCt <- function(data,
                     group.study,
                     group.ref,
                     do.tests = TRUE,
                     save.to.txt = FALSE,
                     name.txt = "RQ_exp_results"){

  data_slim <- data %>%
    filter(Group == group.study | Group == group.ref) %>%
    pivot_longer(cols = -c(Group, Sample), names_to = "Gene", values_to = "value")

  data_slim$Group <- as.factor(data_slim$Group)

  data_FCh <- data_slim %>%
    group_by(Group, Gene) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = value) %>%
    mutate(FCh = .data[[group.study]]/.data[[group.ref]]) %>%
    as.data.frame()

  data_FCh_sd <- data_slim %>%
    group_by(Group, Gene) %>%
    summarise(value_sd = sd(value, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = value_sd) %>%
    rename_with(~paste0(.x, "_sd", recycle0 = TRUE), all_of(c(group.study, group.ref)))

  if (do.tests == TRUE){

    data_FCh_norm <- data_slim %>%
      group_by(Group, Gene) %>%
      summarise(shap_wilka_p = shapiro.test(value)$p.value, .groups = "keep") %>%
      pivot_wider(names_from = Group, values_from = shap_wilka_p) %>%
      rename_with(~paste0(.x, "_norm_p", recycle0 = TRUE), all_of(c(group.study, group.ref))) %>%
      full_join(data_FCh, by = c("Gene")) %>%
      rename_with(~paste0(.x, "_mean", recycle0 = TRUE), all_of(c(group.study, group.ref)))

    data_FCh_tests <- data_slim %>%
      group_by(Gene) %>%
      summarise(t_test_p = t.test(value ~ Group, alternative = "two.sided")$p.value,
                t_test_stat = t.test(value ~ Group, alternative = "two.sided")$statistic,
                MW_test_p = pvalue(wilcox_test(value ~ Group)),
                MW_test_stat = statistic(wilcox_test(value ~ Group)), .groups = "keep")
    data_FCh_norm_tests <- full_join(data_FCh_norm, data_FCh_tests, by = c("Gene"))
    data_FCh_results <- full_join(data_FCh_norm_tests, data_FCh_sd, by = c("Gene"))
    data_FCh_results <- select(data_FCh_results, Gene, ends_with("_mean"), ends_with("_sd"), everything())

    return(data_FCh_results)

  } else {

    data_FCh_results <- full_join(data_FCh, data_FCh_sd, by = c("Gene"))
    data_FCh_results <- rename_with(data_FCh_results, ~paste0(.x, "_mean", recycle0 = TRUE), all_of(c(group.study, group.ref)))
    data_FCh_results <- select(data_FCh_results, Gene, ends_with("_mean"), ends_with("_sd"), everything())

    return(data_FCh_results)
  }
  if (save.to.txt == TRUE){
    write.table(data_FCh_results, paste(name.txt,".txt", sep = ""))
  }
}






#' @title norm_finder
#'
#' @description
#' This function calculate stability scores using a code adapted from the [NormFinder algorithm](https://www.moma.dk/software/normfinder) published in [this article](https://aacrjournals.org/cancerres/article/64/15/5245/511517/Normalization-of-Real-Time-Quantitative-Reverse).
#' This function is internally used by `RQdeltaCT::find_ref_gene()` function and do not need to be use separately.
#'
#' @param data object returned from make_Ct_ready() functions.
#' @param candidates character: vector of names of genes - candidates for reference gene.
#' @param save.to.txt logical: if TRUE, returned table with results will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension. Default to "RQ_exp_results".
#'
#' @return Table with calculated stability scores, the lowest value the best candidate for reference gene.
#'
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
#'                                    candidates = c("Gene4", "Gene8","Gene10","Gene16","Gene17", "Gene18"))
#'
#' @importFrom dplyr filter select
#' @import tidyverse
#'
norm_finder <- function(data,
                        candidates,
                        save.to.txt = FALSE,
                        name.txt = "NormFinder_results"){

  data.t <- data %>%
    pivot_longer(!c(Sample, Group), names_to = "Gene", values_to = "Ct") %>%
    filter(Gene %in% candidates) %>%
    pivot_wider(id_cols = !c(Group), names_from = "Sample", values_from = "Ct")

  last.row <- c("Group", data$Group)
  dat0 <- rbind(as.data.frame(data.t), last.row)
  rownames(dat0) <- dat0[,"Gene"]
  dat0 <- select(dat0, -Gene)

  ntotal = dim(dat0)[2]
  k0 = dim(dat0)[1]
  ngenes=k0-1
  genenames=rownames(dat0)[-k0]
  grId=dat0[k0,]
  dat0=dat0[-k0,]

  dat=matrix(as.numeric(unlist(dat0)),ngenes,ntotal)

  samplenames=colnames(dat0)
  grId=factor(unlist(grId))
  groupnames=levels(grId)
  ngr=length(levels(grId))
  nsamples=rep(0,ngr)

  for (group in 1:ngr){nsamples[group]=sum(grId==groupnames[group])}

  MakeStab=function(da){
    ngenes=dim(da)[1]
    sampleavg=apply(da,2,mean)
    genegroupavg=matrix(0,ngenes,ngr)

    for (group in 1:ngr){
      genegroupavg[,group]=apply(da[,grId==groupnames[group]],1,mean)}

    groupavg=rep(0,ngr)

    for (group in 1:ngr){groupavg[group]=mean(da[,grId==groupnames[group]])}

    GGvar=matrix(0,ngenes,ngr)

    for (group in 1:ngr){
      grset=(grId==groupnames[group])
      a=rep(0,ngenes)

      for (gene in 1:ngenes){
        a[gene]=sum((da[gene,grset]-genegroupavg[gene,group]-
                       sampleavg[grset]+groupavg[group])^2)/(nsamples[group]-1)
      }
      GGvar[,group]=(a-sum(a)/(ngenes*ngenes-ngenes))/(1-2/ngenes)
    }

    genegroupMinvar=matrix(0,ngenes,ngr)

    for (group in 1:ngr){
      grset=(grId==groupnames[group])
      z=da[,grset]

      for (gene in 1:ngenes){
        varpair=rep(0,ngenes)

        for (gene1 in 1:ngenes){varpair[gene1]=var(z[gene,]-z[gene1,])}
        genegroupMinvar[gene,group]=min(varpair[-gene])/4
      }
    }

    GGvar=ifelse(GGvar<0,genegroupMinvar,GGvar)

    dif=genegroupavg
    difgeneavg=apply(dif,1,mean)
    difgroupavg=apply(dif,2,mean)
    difavg=mean(dif)

    for (gene in 1:ngenes){
      for (group in 1:ngr){
        dif[gene,group]=dif[gene,group]-difgeneavg[gene]-difgroupavg[group]+difavg
      }
    }

    nsampMatrix=matrix(rep(nsamples,ngenes),ngenes,ngr,byrow=T)
    vardif=GGvar/nsampMatrix
    gamma=sum(dif*dif)/((ngr-1)*(ngenes-1))-sum(vardif)/(ngenes*ngr)
    gamma=ifelse(gamma<0,0,gamma)

    difnew=dif*gamma/(gamma+vardif)
    varnew=vardif+gamma*vardif/(gamma+vardif)
    Ostab0=abs(difnew)+sqrt(varnew)
    Ostab=apply(Ostab0,1,mean)

    mud=rep(0,ngenes)
    for (gene in 1:ngenes){
      mud[gene]=2*max(abs(dif[gene,]))
    }

    genevar=rep(0,ngenes)
    for (gene in 1:ngenes){
      genevar[gene]=sum((nsamples-1)*GGvar[gene,])/(sum(nsamples)-ngr)
    }
    Gsd=sqrt(genevar)

    return(cbind(mud,Gsd,Ostab,rep(gamma,ngenes),GGvar,dif))
  }

  res = MakeStab(dat)
  ord = order(res[,3])
  FinalRes = data.frame("Stability" = round(res[ord,3],2),
                        row.names = genenames[ord])

  if (save.to.txt == TRUE){
    write.table(FinalRes, paste(name.txt,".txt", sep = ""))
  }

  return(FinalRes)
}





#' @title find_ref_gene
#'
#' @description
#' This function draw a line plot and calculate parameters useful to assess gene expression stability:
#'  minimum, maximum, standard deviation, variance, colinearity coefficient (VIF),
#'  and stability measures from NormFinder and geNorm algorithms. This function could be helpful to select the best
#' reference gene for normalization of Ct values.
#'
#' @param data object returned from make_Ct_ready() functions,
#' @param candidates character: vector of names of genes - candidates for gene reference.
#' @param groups character vector length of two with names of compared groups
#' @param colors character: vector of colors for genes, number of elements should be equal to number of candidate genes (elements in `candidates` vector).
#' @param vif.score logical: if TRUE, VIF colinearity coefficient will be calculated.
#' @param norm.finder.score logical: if TRUE, NormFinder stability score will be calculated.
#' @param genorm.score logical: if TRUE, geNorm stability score will be calculated.
#' @param line.width numeric: width of lines drawn in the plot. Default to 1.
#' @param angle integer: value of angle in which names of genes should be displayed. Default to 0.
#' @param x.axis.title character: title of x axis. Default to "".
#' @param y.axis.title character: title of y axis.  Default to "Ct".
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 10.
#' @param legend.title character: title of legend. Default to "".
#' @param legend.title.size integer: font size of legend title. Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "Ct_reference_gene_selection".
#'
#' @return List containing plot object and table with calculated parameters. Created plot is displayed on graphic device.
#'
#' @export
#'
#' @examples
#'library(car)
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
#'                      vif.score = TRUE,
#'                      norm.finder.score = TRUE,
#'                      genorm.score = TRUE)
#' ref[[2]]
#'
#' @importFrom base mean print min max as.data.frame
#' @importFrom stats sd var
#' @importFrom dplyr filter
#' @importFrom car vif
#' @import ctrlGene
#' @import ggplot2
#' @import tidyverse
#'
find_ref_gene <- function(data,
                          groups,
                          candidates,
                          colors,
                          vif.score = FALSE,
                          norm.finder.score = FALSE,
                          genorm.score = FALSE,
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
                          dpi = 600, width = 15, height = 15,
                          name.tiff = "Ct_reference_gene_selection"){

  ref <- data %>%
      #filter(Group == groups[1] | Group == groups[2]) %>%
      pivot_longer(cols = -c(Group, Sample), names_to = "Gene", values_to = "Ct") %>%
      filter(Gene %in% candidates)

  ref_plot <- ggplot(ref, aes(x = Sample, y = Ct, color = Gene, group = Gene)) +
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
    theme(legend.title = element_text(size = legend.title.size, colour="black")) +
    theme(legend.text = element_text(size = legend.text.size, colour="black")) +
	  theme(plot.title = element_text(size = plot.title.size))

  print(ref_plot)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff, ".tiff", sep = ""), ref_plot, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }

  ref_var <- ref %>%
    group_by(Gene) %>%
    summarise(min = min(Ct),
              max = max(Ct),
              sd = sd(Ct, na.rm = TRUE),
              var = var(Ct, na.rm = TRUE), .groups = "keep") %>%
    as.data.frame()

  if (vif.score == TRUE){
  ref_lm <- data %>%
    filter(Group %in% groups[1]) %>%
    ungroup()

  ref_lm$dum <- sample(1:nrow(ref_lm), nrow(ref_lm), replace = FALSE)
  model <- lm(dum ~ ., data = select(ref_lm, -Sample, -Group))
  vif <- vif(model)
  vif_sel <- vif[names(vif) %in% candidates]
  ref_var$VIF <- vif_sel
  colnames(ref_var)[colnames(ref_var) == "VIF"] = paste(groups[1], "_VIF", sep = "")

  ref_lm2 <- data %>%
    filter(Group %in% groups[2]) %>%
    ungroup()

  ref_lm2$dum <- c(1:nrow(ref_lm2))
  model2 <- lm(dum ~ ., data = select(ref_lm2, -Sample, -Group))
  vif2 <- vif(model2)
  vif_sel2 <- vif2[names(vif2) %in% candidates]
  ref_var$VIF2 <- vif_sel2
  colnames(ref_var)[colnames(ref_var) == "VIF2"] = paste(groups[2], "_VIF", sep = "")
  }

if (norm.finder.score == TRUE){
  reference.stability.nF <- norm_finder(data, candidates = candidates)
  colnames(reference.stability.nF) <- "NormFinder_score"
  reference.stability.nF$Gene <- rownames(reference.stability.nF)
  ref_var <- ref_var %>%
    full_join(reference.stability.nF, by = join_by(Gene))
}
  if (genorm.score == TRUE){
  data <- data %>%
    ungroup() %>%
    select(-Group) %>%
    select(any_of(c("Sample", candidates))) %>%
    as.data.frame()

  rownames(data) <- data[,"Sample"]
  data <- select(data, -Sample)
  data <- as.matrix(data)

  reference.stability.gF <- geNorm(data,
         genes = data.frame(Gene = character(0), geNorm_score = numeric(0)),
         ctVal = TRUE)
  colnames(reference.stability.gF) <- c("Gene", "geNorm_score")
  ref_var <- ref_var %>%
    full_join(reference.stability.gF, by = join_by(Gene))

  }

  return(list(ref_plot, as.data.frame(ref_var)))
}







#' @title delta_Ct
#'
#' @description
#' This function calculates delta Ct (dCt) values by subtracting Ct values of reference gene from Ct values of other genes.
#'
#' @param data data object returned from make_Ct_ready function,
#' @param ref character: name of reference gene.
#' @param save.to.txt logical: if TRUE, returned data will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension. Default to "data_dCt".
#'
#' @return Data.frame with dCt values and printed information about number and percentage of missing values.
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
#'
#' @importFrom base as.data.frame colnames paste
#' @importFrom utils write.table
#' @importFrom dplyr mutate_at select
#' @import tidyverse
#'
delta_Ct <- function(data,
                     ref,
					           save.to.txt = FALSE,
                     name.txt = "data_dCt"){

  dCt <- mutate_at(data,
                   vars(-c("Group", "Sample", all_of(ref))),
                   list(dCt = ~ . - .data[[ref]]))
  dCt <- select(dCt, Group, Sample, ends_with("dCt"))
  colnames(dCt) <- sub("_dCt*", "", colnames(dCt))

  if (save.to.txt == TRUE){
    write.table(as.data.frame(dCt), paste(name.txt,".txt", sep = ""))
  }

  return(dCt)
}







#' @title control_boxplot_sample
#'
#' @description
#' Boxplot illustrating distribution of data in each sample. Could be useful to identify outlier samples.
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param sel.Sample character vector with names of samples to include, or "all" (default) to use all samples.
#' @param coef numeric: how many times of interquartile range should be used to indicate the most extend data point for whiskers. Default to 1.5.
#' @param colors character vector containing colors for compared groups. Numbers of colors must be equal to number of groups. Default to c("#66c2a5", "#fc8d62").
#' @param x.axis.title character: title of x axis. Default to "Sample".
#' @param y.axis.title character: title of y axis. Default to "value".
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 12.
#' @param legend.title character: title of legend. Default to "Group".
#' @param legend.title.size integer: font size of legend title. Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_boxplot_samples".
#'
#' @return Object with boxplot illustrating distribution of data in each sample. Created plot is also displayed on graphic device.
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
#' @importFrom base print
#' @import ggplot2
#' @import tidyverse
#'
control_boxplot_sample <- function(data,
                                   sel.Sample = "all",
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
                                   dpi = 600, width = 15, height = 15,
                                   name.tiff = "control_boxplot_samples"){
  if (sel.Sample[1] == "all"){
    data <- data

  } else {

    data <- filter(data, Sample %in% sel.Sample)
  }

    data <- pivot_longer(data, !c(Sample, Group), names_to = "Gene" , values_to = "value")

    box_control_sample <- ggplot(data, aes(x = Sample, y = value, color = Group)) +
      geom_boxplot(coef = coef) +
      scale_x_discrete(limits = rev(unique(data$Sample))) +
      scale_color_manual(values = c(colors)) +
      coord_flip() +
      xlab(x.axis.title) + ylab(y.axis.title) +
      labs(color = legend.title, title = plot.title) +
      theme_classic() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
	    theme(plot.title = element_text(size = plot.title.size))

    print(box_control_sample)

	if (save.to.tiff == TRUE){
      ggsave(paste(name.tiff,".tiff", sep = ""), box_control_sample, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
	}
    return(box_control_sample)
}






#' @title control_boxplot_gene
#'
#' @description
#' This function creates boxplot illustrating distribution of data in each gene. Could be useful to compare expression of analyzed genes.
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param sel.Gene character vector with names of genes to include, or "all" (default) to use all names of genes.
#' @param by.group logical: if TRUE, distributions will be drawn by compared groups of samples.
#' @param coef numeric: how many times of interquartile range should be used to indicate the most extend data point for whiskers.
#' @param colors character vector containing colors for groups, length of one (when by.group = FALSE) or equal to number of groups (when by.group = TRUE).
#' @param coef numeric: how many times of interquartile range should be used to indicate the most extend data point for whiskers. Default to 1.5.
#' @param x.axis.title character: title of x axis. Default to "Gene".
#' @param y.axis.title character: title of y axis. Default to "value".
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 12.
#' @param legend.title character: title of legend. Default to "Group".
#' @param legend.title.size integer: font size of legend title. Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_boxplot_genes".
#'
#' @return Object with boxplot illustrating distribution of data for each gene. Created plot will be also displayed on graphic device.
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
#' @importFrom base print
#' @import ggplot2
#' @import tidyverse
#'
control_boxplot_gene <- function(data,
                                   sel.Gene = "all",
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
                                   dpi = 600, width = 15, height = 15,
                                   name.tiff = "control_boxplot_genes"){
  if (sel.Gene[1] == "all"){
    data <- data

  } else {

    data <- select(data, Group, Sample, any_of(sel.Gene))
  }

  data <- pivot_longer(data, !c(Sample, Group), names_to = "Gene" , values_to = "value")

  if (by.group == TRUE){
    box_control_genes <- ggplot(data, aes(x = Gene, y = value, color = Group)) +
      geom_boxplot(coef = coef) +
      scale_color_manual(values = c(colors)) +
      labs(color = legend.title, title = plot.title)

  } else {

    box_control_genes <- ggplot(data, aes(x = Gene, y = value)) +
      geom_boxplot(coef = coef, fill = colors[1]) +
      labs(title = plot.title)
  }

  box_control_genes <- box_control_genes +
    scale_x_discrete(limits = rev(unique(data$Gene))) +
    coord_flip() +
    xlab(x.axis.title) +
    ylab(y.axis.title) +
    theme_classic() +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour="black")) +
    theme(legend.title = element_text(size = legend.title.size, colour="black")) +
    theme(legend.text = element_text(size = legend.text.size, colour="black")) +
    theme(plot.title = element_text(size = plot.title.size))


  print(box_control_genes)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), box_control_genes, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
  return(box_control_genes)
}





#' @title control_cluster_sample
#'
#' @description
#' Performs hierarchical clustering of samples based on the data. Could be useful to identify outlier samples.
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param sel.Sample character vector with names of samples to include, or "all" (default) to use all samples.
#' @param method.dist character: name of method used for calculation of distances, derived from stats::dist() function, should be one of "euclidean" (default) , "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param method.clust character: name of used method for agglomeration, derived from stats::hclust() function, should be one of "ward.D", "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".
#' @param x.axis.title character: title of x axis. Default to "Samples".
#' @param y.axis.title character: title of y axis. Default to "Height".
#' @param label.size numeric: size of text labels. Default to 1.
#' @param plot.title character: title of plot. Default to "".
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_clust_samples".
#'
#' @return Plot with hierarchical clustering of samples, displayed on graphic device.
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
#' control_cluster_sample(data.dCt)
#'
#' @importFrom stats hclust dist
#' @importFrom base plot
#' @importFrom dplyr select
#' @import tidyverse
#'
control_cluster_sample <- function(data,
                                   sel.Sample = "all",
                                   method.dist = "euclidean",
                                   method.clust = "average",
                                   x.axis.title = "Samples",
                                   y.axis.title = "Height",
                                   label.size = 1,
                                   plot.title = "",
                                   save.to.tiff = FALSE,
                                   dpi = 600, width = 15, height = 15,
                                   name.tiff = "control_clust_samples"){
  if (sel.Sample[1] == "all"){
    data <- data

  } else {

    data <- filter(data, Sample %in% sel.Sample)
  }

  data <- ungroup(data)
  cluster <- hclust(dist(select(data, -Group, -Sample), method = method.dist), method = method.clust)
  cluster$labels <- data$Sample
  plot(cluster, xlab = x.axis.title, ylab = y.axis.title, main = plot.title, cex = label.size)

  if (save.to.tiff == TRUE){
    tiff(paste(name.tiff, ".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
    plot(cluster, xlab  = x.axis.title, ylab = y.axis.title, main = plot.title, cex = label.size)
    dev.off()
  }
}






#' @title control_cluster_gene
#'
#' @description
#' Performs hierarchical clustering of genes based on the data. Could be useful to gain insight into similarity in expression of analyzed genes.
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param sel.Gene character vector with names of genes to include, or "all" (default) to use all names of genes.
#' @param method.dist character: name of method used for calculation of distances, derived from stats::dist() function, should be one of "euclidean" (default) , "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param method.clust character: name of used method for agglomeration, derived from stats::hclust() function, should be one of "ward.D", "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".
#' @param x.axis.title character: title of x axis. Default to "Genes".
#' @param y.axis.title character: title of y axis. Default to "Height".
#' @param plot.title character: title of plot. Default to "".
#' @param label.size numeric: size of text labels. Defaault to 1.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_clust_genes".
#'
#' @return Plot with hierarchical clustering of genes, displayed on graphic device.
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
#' @importFrom base plot t colnames
#' @importFrom dplyr select
#' @import tidyverse
#'
control_cluster_gene <- function (data,
                                    method.dist = "euclidean",
                                    sel.Gene = "all",
                                    method.clust = "average",
                                    x.axis.title = "Genes",
                                    y.axis.title = "Height",
                                    label.size = 1,
                                    plot.title = "",
                                    save.to.tiff = FALSE,
                                    dpi = 600, width = 15, height = 15,
                                    name.tiff = "control_clust_genes")
{
  if (sel.Gene[1] == "all"){
    data <- data

  } else {

    data <- select(data, Group, Sample, any_of(sel.Gene))
  }

  data_t <- ungroup(data)
  data_t <- t(select(data_t, -Group, -Sample))
  colnames(data_t) <- data$Sample
  cluster <- hclust(dist(as.data.frame(data_t), method = method.dist), method = method.clust)
  plot(cluster, xlab = x.axis.title, ylab = y.axis.title, main = plot.title, cex = label.size)

  if (save.to.tiff == TRUE) {
    tiff(paste(name.tiff, ".tiff", sep = ""), res = dpi,
         width = width, height = height, units = "cm", compression = "lzw")
    plot(cluster, xlab = x.axis.title, ylab = y.axis.title,
         main = plot.title, cex = label.size)
    dev.off()
  }
}





#' @title control_pca_sample
#'
#' @description
#' Performs principal component analysis (PCA) for samples and generate plot illustrating spatial arrangement of samples using two first components. Could be useful to identify outlier samples.
#'     IMPORTANT: PCA analysis can not deal with missing values, thus all samples with at least one missing value are removed from data before analysis. It is recommended to run this function on data after imputation of missing values.
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param sel.Sample character vector with names of samples to include, or "all" (default) to use all samples.
#' @param point.size numeric: size of points. Default to 4.
#' @param point.shape integer: shape of points. Default to 19.
#' @param alpha numeric: transparency of points, a value between 0 and 1. Default to 0.7.
#' @param label.size numeric: size of points labels (names of samples). Default to 3.
#' @param hjust numeric: horizontal position of points labels. Default to 0.
#' @param vjust numeric: vertical position of points labels.  Default to -1.
#' @param colors character vector containing colors for compared groups. Numbers of colors must be equal to number of groups. Default to c("#66c2a5", "#fc8d62").
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 10.
#' @param legend.title character: title of legend. Default to "Group".
#' @param legend.title.size integer: font size of legend title. Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_pca_samples".
#'
#' @return Object with plot illustrating spatial arrangement of samples according to coordinates of two first components obtained from principal component analysis (PCA). The plot will be also displayed in graphic device.
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
#' @importFrom base print as.data.frame rownames summary paste round
#' @importFrom stats na.omit prcomp
#' @import ggplot2
#' @import tidyverse
#'
control_pca_sample <- function(data,
                               sel.Sample = "all",
                               point.size = 4,
                               point.shape = 19,
                               alpha = 0.7,
                               colors = c("#66c2a5", "#fc8d62"),
                               label.size = 3,
                               hjust = 0,
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
                               dpi = 600, width = 15, height = 15,
                               name.tiff = "control_pca_samples"){
  if (sel.Sample[1] == "all"){
    data <- data

  } else {

    data <- filter(data, Sample %in% sel.Sample)
  }

  data <- as.data.frame(data)
  rownames(data) <- data$Sample
  data <- na.omit(data)
  pca <- prcomp(select(data, -Sample, -Group), scale = TRUE)
  var_pca1 <- summary(pca)$importance[2,][1]
  var_pca2 <- summary(pca)$importance[2,][2]
  pca_comp <- as.data.frame(pca$x)
  pca_comp$Sample <- data$Sample
  pca_comp$Group <- data$Group

  control_pca <- ggplot(pca_comp, aes(x = PC1, y = PC2, label = Sample, color = Group)) +
    geom_point(size = point.size, shape = point.shape, alpha = alpha) +
    scale_color_manual(values = c(colors)) +
    labs(colour = legend.title, title = plot.title) +
    theme_bw() +
    labs(x = paste("PC1: ", round(var_pca1*100,2), "% variance explained", sep = ""),
         y = paste("PC2: ", round(var_pca2*100,2), "% variance explained", sep = "")) +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour="black")) +
    theme(legend.text = element_text(size = legend.text.size, colour="black")) +
    theme(legend.title = element_text(size = legend.title.size, colour="black")) +
    theme(plot.title = element_text(size = plot.title.size)) +
    geom_text(aes(label = Sample), hjust = hjust, vjust = vjust, size = label.size)

  print(control_pca)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), control_pca, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
  return(control_pca)
}









#' @title control_pca_gene
#'
#' @description
#' Performs principal component analysis (PCA) for genes and generate plot illustrating spatial arrangement of genes using 2 components. Could be useful to gain insight into similarity in expression of analyzed genes.
#' IMPORTANT: PCA analysis can not deal with missing values, thus all genes with at least one missing value are removed from data before analysis. It is recommended to run this function on data after imputation of missing values.
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param sel.Gene character vector with names of genes to include, or "all" (default) to use all names of genes.
#' @param point.size numeric: size of points. Default to 4.
#' @param point.shape integer: shape of points. Default to 19.
#' @param alpha numeric: transparency of points, a value between 0 and 1. Default to 0.7.
#' @param label.size numeric: size of points labels (names of samples). Default to 3.
#' @param hjust numeric: horizontal position of points labels. Default to 0.
#' @param vjust numeric: vertical position of points labels.  Default to -1.
#' @param color character: color used for points.
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 10.
#' @param legend.title character: title of legend. Default to "Group".
#' @param legend.title.size integer: font size of legend title. Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_pca_genes".
#'
#' @return Object with plot illustrating spatial arrangement of genes according to coordinates of 2 components obtained from principal component analysis (PCA). The plot will be also displayed in graphic device.
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
#' @importFrom base print as.data.frame rownames summary paste round t colnames
#' @importFrom stats na.omit prcomp
#' @import ggplot2
#' @import tidyverse
#'
control_pca_gene <- function(data,
                               sel.Gene = "all",
                               point.size = 4,
                               point.shape = 19,
                               alpha = 0.7,
                               label.size = 3, hjust = 0, vjust = -1,
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
                               dpi = 600, width = 15, height = 15,
                              name.tiff = "control_pca_genes"){
  if (sel.Gene[1] == "all"){
    data <- data

  } else {

    data <- select(data, Group, Sample, any_of(sel.Gene))
  }

  data <- ungroup(data)
  data <- as.data.frame(data)
  data_t <- t(select(data, -Group, -Sample))
  colnames(data_t) <- data$Sample
  data_t <- na.omit(data_t)
  pca <- prcomp(data_t, scale = TRUE)
  var_pca1 <- summary(pca)$importance[2,][1]
  var_pca2 <- summary(pca)$importance[2,][2]
  pca_comp <- as.data.frame(pca$x)
  pca_comp$Gene <- rownames(pca_comp)

  control_pca <- ggplot(pca_comp, aes(x = PC1, y = PC2, label = Gene)) +
    geom_point(size = point.size, shape = point.shape, alpha = alpha, col = color) +
    labs(title = plot.title) +
    theme_bw() +
    labs(x = paste("PC1: ", round(var_pca1*100,2), "% variance explained", sep = ""),
         y = paste("PC2: ", round(var_pca2*100,2), "% variance explained", sep = "")) +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour="black")) +
    theme(plot.title = element_text(size = plot.title.size)) +
    geom_text(aes(label = Gene), hjust = hjust, vjust = vjust, size = label.size)

  print(control_pca)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), control_pca, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
  return(control_pca)
}





#' @title corr_gene
#'
#' @description
#' Performs correlation analysis of genes based on the data. Could be useful to gain insight into relationships between analyzed genes.
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param sel.Gene character vector with names of genes to include, or "all" (default) to use all names of genes.
#' @param add.coef if coefficients should be add to the plot, specify color of coefficients (default to "black"), otherwise set to NULL.
#' @param method character: type of correlations to compute, specify "pearson" (default) for Pearson's correlation coefficients
#' or "spearman" for Spearman's rank correlation coefficients.
#' @param order character: method used for ordering the correlation matrix (see documentation for corrplot::corrplot() function),
#' one of the "original" (original order), "AOE" (angular order of the eigenvectors),
#' "FPC" (first principal component order), "hclust" (hierarchical clustering order, default), or "alphabet" (alphabetical order).
#' @param hclust.method character: name of used method for hclust agglomeration, should be one of "ward", ward.D", "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".
#' @param size numeric: size of variable names and numbers in legend. Default to 0.6.
#' @param coef.size numeric: size of correlation coefficients. Default to 0.6.
#' @param p.adjust.method character: p value correction method for multiple testing, one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", or "none". See documentation for stats::p.adjust() function for details.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "corr_genes".
#' @param save.to.txt logical: if TRUE, correlation results (sorted by descending absolute values of correlation coefficients) will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension.. Default to "corr_genes".
#'
#' @return Plot illustrating correlation matrix (displayed in graphic device) and data.frame with computed correlation coefficients and p values.
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
#' @importFrom base as.data.frame paste upper.tri rownames
#' @importFrom stats p.adjust
#' @importFrom utils write.table
#' @importFrom Hmisc rcorr
#' @import corrplot
#' @import tidyverse
#'
corr_gene <- function(data,sel.Gene = "all",
                                method = "pearson",
                                add.coef = "black",
                                order = "hclust",
                                hclust.method = "average",
                                size = 0.6,
                                coef.size = 0.6,
                                p.adjust.method = "BH",
                                save.to.tiff = FALSE,
                                dpi = 600, width = 15, height = 15,
                                name.tiff = "corr_genes",
                                save.to.txt = FALSE,
                                name.txt = "corr_genes"){
  if (sel.Gene[1] == "all"){
    data <- data

  } else {

    data <- select(data, Group, Sample, any_of(sel.Gene))
  }

  data <- as.data.frame(data)
  data <- select(data, -Group, -Sample)
  res_cor <- rcorr(as.matrix(data), type = method)

  if (order == "hclust"){
    corrplot(res_cor$r,
             type = "upper",
             addCoef.col = add.coef,
             tl.cex = size,
             cl.cex = size,
             tl.col = "black",
             number.cex = coef.size,
             order = order,
             hclust.method = hclust.method)

    if (save.to.tiff == TRUE){
      tiff(paste(name.tiff,".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
      corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
               order = order,
               hclust.method = hclust.method,
               addCoef.col = add.coef,
               number.cex = coef.size)
      dev.off()
    }
  } else {

    corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
             order = order,
             addCoef.col = add.coef,
             number.cex = coef.size)

    if (save.to.tiff == TRUE){
      tiff(paste(name.tiff,".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
      corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
               order = order,
               addCoef.col = "black",
               number.cex = coef.size)
      dev.off()
    }
  }
    corr <- upper.tri(res_cor$r)
    corr.data <- data.frame(
      row = rownames(res_cor$r)[row(res_cor$r)[corr]],
      column = rownames(res_cor$r)[col(res_cor$r)[corr]],
      cor  =(res_cor$r)[corr],
      p = res_cor$P[corr]
    )
    corr.data$p.adj <- p.adjust(corr.data$p, method = p.adjust.method)
    corr.data.sort <- arrange(corr.data, -abs(cor))

     if (save.to.txt == TRUE){

    write.table(corr.data.sort, paste(name.txt, ".txt", sep = ""))
  }
  return(as.data.frame(corr.data.sort))
}





#' @title corr_sample
#'
#' @description
#' Performs correlation analysis of samples based on the data. Could be useful to gain insight into relationships between analyzed samples.
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param sel.Sample character vector with names of samples to include, or "all" (default) to use all samples.
#' @param add.coef if coefficients should be add to the plot, specify color of coefficients (default to "black"), otherwise set to NULL.
#' @param method character: type of correlations to compute, specify "pearson" (default) for Pearson's correlation coefficients
#' or "spearman" for Spearman's rank correlation coefficients.
#' @param order character: method used for ordering the correlation matrix (see documentation for corrplot::corrplot() function),
#' one of the "original" (original order), "AOE" (angular order of the eigenvectors),
#' "FPC" (first principal component order), "hclust" (hierarchical clustering order, default), or "alphabet" (alphabetical order).
#' @param hclust.method character: name of used method for hclust agglomeration, should be one of "ward", ward.D", "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".
#' @param size numeric: size of variable names and numbers in legend. Default to 0.6.
#' @param coef.size numeric: size of correlation coefficients. Default to 0.6.
#' @param p.adjust.method character: p value correction method for multiple testing, one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", or "none". See documentation for stats::p.adjust() function for details.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "corr_samples".
#' @param save.to.txt logical: if TRUE, correlation results (sorted by descending absolute values of correlation coefficients) will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension.. Default to "corr_samples".
#'
#' @return Plot illustrating correlation matrix (displayed in graphic device) and data.frame with computed correlation coefficients and p values.
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
#' @importFrom base as.data.frame paste upper.tri rownames
#' @importFrom stats p.adjust
#' @importFrom utils write.table
#' @importFrom Hmisc rcorr
#' @import corrplot
#' @import tidyverse
#'
corr_sample <- function(data, sel.Sample = "all",
                                method = "pearson",
                                add.coef = "black",
                                order = "hclust",
                                hclust.method = "average",
                                size = 0.6,
                                coef.size = 0.6,
                                p.adjust.method = "BH",
                                save.to.tiff = FALSE,
                                dpi = 600, width = 15, height = 15,
                                name.tiff = "corr_samples",
                                save.to.txt = FALSE,
                                name.txt = "corr_samples"){
  if (sel.Sample[1] == "all"){
    data <- data

  } else {
    data <- filter(data, Sample %in% sel.Sample)
  }

  data <- as.data.frame(data)
  data_t <- t(select(data, -Group, -Sample))
  colnames(data_t) <- data$Sample
  res_cor <- rcorr(as.matrix(data_t), type = method)

  if (order == "hclust"){
    corrplot(res_cor$r,
             type = "upper",
             addCoef.col = add.coef,
             tl.cex = size,
             cl.cex = size,
             tl.col = "black",
             number.cex = coef.size,
             order = order,
             hclust.method = hclust.method)

    if (save.to.tiff == TRUE){
      tiff(paste(name.tiff,".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
      corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
               order = order,
               hclust.method = hclust.method,
               addCoef.col = add.coef,
               number.cex = coef.size)
      dev.off()
    }
  } else {

    corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
             order = order,
             addCoef.col = add.coef,
             number.cex = coef.size)

    if (save.to.tiff == TRUE){
      tiff(paste(name.tiff,".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
      corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
               order = order,
               addCoef.col = "black",
               number.cex = coef.size)
      dev.off()
    }
  }

    corr <- upper.tri(res_cor$r)
    corr.data <- data.frame(
      row = rownames(res_cor$r)[row(res_cor$r)[corr]],
      column = rownames(res_cor$r)[col(res_cor$r)[corr]],
      cor  =(res_cor$r)[corr],
      p = res_cor$P[corr]
    )
    corr.data$p.adj <- p.adjust(corr.data$p, method = p.adjust.method)
    corr.data.sort <- arrange(corr.data, -abs(cor))

    if (save.to.txt == TRUE){
    write.table(corr.data.sort, paste(name.txt, ".txt", sep = ""))
    }

  return(as.data.frame(corr.data.sort))
}







#' @title single_pair_gene
#'
#' @description
#' Generate scatter plot with linear regression line for two specified genes. Could be useful to assess linear relationship between these genes.
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param x,y characters: names of genes to use.
#' @param by.group logical: if TRUE (default), relationships will be shown separately for compared groups.
#' @param point.size numeric: size of points. Default to 4.
#' @param point.shape integer: shape of points. Default to 19.
#' @param point.alpha numeric: transparency of points, a value between 0 and 1. Default to 0.7.
#' @param colors character vector containing colors for compared groups. Numbers of colors must be equal to number of groups. Default to c("#66c2a5", "#fc8d62").
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 10.
#' @param legend.title character: title of legend. Default to "Group".
#' @param legend.title.size integer: font size of legend title. Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for legend.position parameter in ggplot2::theme() function.
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param labels logical: if TRUE (default), a regression statistics will be added to the plot using ggpimsc::stat_poly_eq() function.
#' @param label character: which regression statistics should be drawn, use names specified by ggpimsc::stat_poly_eq() function. Default to c("eq", "R2", "p") for regression equation, coefficient of determination and p value, respectively.
#' @param label.position.x,label.position.y  numeric: coordinates for position of regression statistics. If by.group = TRUE, two values could be provided for each of these parameters (for different regression lines). See description of label.x and label.y parameters from ggpimsc::stat_poly_eq() function.
#' @param small.p,small.r logical, if TRUE, p in p value label and r in coefficient of determination label, will be lowercase. Default to FALSE.
#' @param rr.digits,p.digits integer: number of digits after the decimal point in coefficient of determination and p value in labels. Default to 2 for rr.digits and 3 for p.digits.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "samples_single_plot".
#'
#' @return Object with plot illustrating spatial arrangement of samples according to coordinates of 2 components obtained from principal component analysis (PCA). Created plot will be also displayed on graphic device.
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
#' @importFrom base print paste
#' @import ggplot2
#' @import ggpmisc
#'
single_pair_gene <- function(data, x, y, by.group = TRUE,
                            point.size = 4, point.shape = 19, point.alpha = 0.7,
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
                            label.position.x = c(1,1),
                            label.position.y = c(1,0.95),
                            small.p = FALSE,
                            small.r = FALSE,
                            p.digits = 3,
                            rr.digits = 2,
                            save.to.tiff = FALSE,
                            dpi = 600, width = 15, height = 15,
                            name.tiff = "genes_single_pair_plot"){

  if (by.group == TRUE){
    single_pair <- ggplot(data, aes(x = .data[[x]], y = .data[[y]], color = Group)) +
      geom_point(size = point.size, shape = point.shape, alpha = point.alpha) +
      geom_smooth(method='lm', se = FALSE) +
      scale_color_manual(values = c(colors)) +
      xlab(x) +
      ylab(y) +
      labs(color = legend.title, title = plot.title) +
      theme_classic() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(plot.title = element_text(size = plot.title.size))

    if (labels == TRUE){
      single_pair <- single_pair +
        stat_poly_eq(use_label(label),
                     label.y = c(label.position.y),
                     label.x = c(label.position.x),
                     small.p = small.p,
                     small.r = small.r,
                     p.digits = p.digits,
                     rr.digits = rr.digits)
    }
  } else {

    single_pair <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
      geom_point(size = point.size, shape = point.shape, alpha = point.alpha, col = colors[1]) +
      geom_smooth(method='lm', se = FALSE) +
      stat_poly_eq(use_label(label),
                   label.y = c(label.position.y),
                   label.x = c(label.position.x)) +
      xlab(x) + ylab(y) +
      labs(title = plot.title) +
      theme_classic() +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(plot.title = element_text(size = plot.title.size))

    if (labels == TRUE){
      single_pair <- single_pair +
        stat_poly_eq(use_label(label),
                     label.y = c(label.position.y),
                     label.x = c(label.position.x),
                     small.p = small.p,
                     small.r = small.r,
                     p.digits = p.digits,
                     rr.digits = rr.digits)
      }
    }

  print(single_pair)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), single_pair, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
  return(single_pair)
}






#' @title single_pair_sample
#'
#' @description
#' Generate scatter plot with linear regression line for two specified samples Could be useful to assess linear relationship between these samples.
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param x,y characters: names of genes to use.
#' @param point.size numeric: size of points.
#' @param point.shape integer: shape of points.
#' @param point.alpha numeric: transparency of points, a value between 0 and 1.
#' @param color character: color used for points.
#' @param axis.title.size integer: font size of axis titles.
#' @param axis.text.size integer: font size of axis text.
#' @param plot.title character: title of plot.
#' @param plot.title.size integer: font size of plot title.
#' @param labels logical: if TRUE, a regression statistics will be added to the plot using ggpimsc::stat_poly_eq() function.
#' @param label character: which regression statistics should be drawn, use names specified by ggpimsc::stat_poly_eq() function. Default to c("eq", "R2", "p") for regression equation, coefficient of determination and p value, respectively.
#' @param label.position.x,label.position.y  numeric: coordinates for position of regression statistics. If by.group = TRUE, two values could be provided for each of these parameters (for different regression lines). See description of label.x and label.y parameters from ggpimsc::stat_poly_eq() function.
#' @param small.p,small.r logical, if TRUE, p in p value label and r in coefficient of determination label, will be lowercase.
#' @param rr.digits,p.digits integer: number of digits after the decimal point in coefficient of determination and p value in labels.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file.
#' @param dpi integer: resolution of saved .tiff file.
#' @param width numeric: width (in cm) of saved .tiff file.
#' @param height integer: height (in cm) of saved .tiff file.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension.
#'
#' @return Plot illustrating spatial arrangement of samples according to coordinates of 2 components obtained from principal component analysis (PCA).
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
#' @importFrom base print paste t colnames
#' @importFrom dplyr select
#' @import ggplot2
#' @import ggpmisc
#'
single_pair_sample <- function(data, x, y,
                                   point.size = 4, point.shape = 19, point.alpha = 0.7,
                                   color = "black",
                                   axis.title.size = 11,
                                   axis.text.size = 10,
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
                                   dpi = 600, width = 15, height = 15,
                                   name.tiff = "samples_single_pair_plot"){

  data <- as.data.frame(data)
  data_t <- t(select(data, -Group, -Sample))
  colnames(data_t) <- data$Sample

  single_pair_t <- ggplot(as.data.frame(data_t), aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(size = point.size, shape = point.shape, alpha = point.alpha, color = color) +
    geom_smooth(method='lm', se = FALSE) +
    xlab(x) +
    ylab(y) +
    labs(title = plot.title) +
    theme_classic() +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour="black")) +
    theme(plot.title = element_text(size = plot.title.size))

  if (labels == TRUE){
    single_pair_t <- single_pair_t +
      stat_poly_eq(use_label(label),
                   label.y = c(label.position.y),
                   label.x = c(label.position.x),
                   small.p = small.p,
                   small.r = small.r,
                   p.digits = p.digits,
                   rr.digits = rr.digits)
  }

  print(single_pair_t)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), single_pair_t, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
  return(single_pair_t)
}







#' @title filter_transformed_data
#'
#' @description
#' Filters transformed Ct data (2^(-Ct), delta Ct, and 2^(-dCt) data) according to the used filtering criteria (see parameters).
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param remove.Gene character: vector with names of genes which should be removed from data.
#' @param remove.Sample character: vector with names of samples which should be removed from data.
#' @param remove.Group character: vector with names of groups which should be removed from data.
#'
#' @return Data.frame with filtered data.
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
#' data.dCt.exp <- exp_Ct_dCt(data.dCt)
#' data.dCt.expF <- filter_transformed_data(data.dCt.exp, remove.Sample = c("Control11"))
#'
#' dim(data.dCt.exp)
#' dim(data.dCt.expF)
#'
#'
#' @importFrom dplyr filter select
#' @import tidyverse
#'
filter_transformed_data <- function(data,
                                  remove.Gene = c(""),
                                  remove.Sample = c(""),
                                  remove.Group = c("")){

  data <- filter(data,
                 !Sample %in% remove.Sample,
                 !Group %in% remove.Group)

  data <- select(data, -any_of(remove.Gene))

  return(data)
}












#' @title results_boxplot
#'
#' @description
#' This function creates boxplot illustrating distribution of data in selected genes.
#'     It is similar to control_boxplot_gene() function; however, some new options are added,
#'     including gene selection, faceting, and adding mean labels to boxes.
#'     This, this function could be useful to present results for finally selected genes.
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param coef numeric: how many times of interquartile range should be used to indicate the most extend data point for whiskers. Default to 1.5.
#' @param sel.Gene character vector with names of genes to include, or "all" (default) to use all names of genes.
#' @param by.group logical: if TRUE (default), distributions will be drawn by compared groups of samples.
#' @param signif.show logical: if TRUE, labels for statistical significance will be added to the plot.
#' @param signif.labels character vector with statistical significance labels (ex. "ns.","***", etc.) with number
#' of elements equal to nimber of genes used for plotting.
#' The ggsignif package was used for convenient adding labels, but there is one tricky point:
#' the same elements of labels can not be handled by used package and must be different.
#' It could be achieved by adding symmetrically white spaces to repeated labels, ex. "ns.", " ns. ", "  ns.  ".
#' @param signif.length numeric: length of horizontal bars, values from 0 to 1.
#' @param signif.dist numeric: distance between errorbar and significance label.
#' Could be in y axis units (if `faceting` = TRUE) or fraction of y axis value reached by errorbar (mean + sd value) (if `faceting` = TRUE).
#' @param faceting logical: if TRUE (default), plot will be drawn with facets with free scales using ggplot2::facet_wrap() function (see its documentation for more details).
#' @param facet.row,facet.col integer: number of rows and columns to arrange facets.
#' @param angle integer: value of angle in which names of genes should be displayed. Default to 0.
#' @param y.exp.low,y.exp.up numeric: space between data on the plot and lower or upper axis. Useful to add extra space for statistical significance labels when `faceting` = TRUE.
#' @param rotate logical: if TRUE, boxplots will be arranged horizontally. Deafault to FALSE.
#' @param add.mean logical: if TRUE, means will be added to boxes as squares. Default to TRUE.
#' @param add.mean.size numeric: size of squares indicating means. Default to 2.
#' @param add.mean.color character: color of squares indicating means. Default to "black".
#' @param colors character vector length of one (when by.group = FALSE) or more (when by.group = TRUE), containing colors for groups. Numbers of colors must be equal to number of groups. Default to c("#66c2a5", "#fc8d62").
#' @param x.axis.title character: title of x axis. Default to "Gene".
#' @param y.axis.title character: title of y axis. Default to "value".
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 10.
#' @param legend.title character: title of legend. Default to "Group".
#' @param legend.title.size integer: font size of legend title. Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "results_boxplot".
#'
#' @return Object with boxplot illustrating distribution of data for selected genes. Created plot will be also displayed on graphic device.
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
#' data.dCt.exp <- exp_delta_Ct(data.dCt)
#' data.dCt.expF <- filter_transformed_data(data.dCt.exp, remove.Sample = c("Control11"))
#' results_boxplot(data.dCt.exp,
#'                 sel.Gene = c("Gene1","Gene16","Gene19","Gene20"),
#'                 signif.labels = c("****","*","***"," * "),
#'                 angle = 30,
#'                 signif.dist = 1.05,
#'                 facet.row = 1,
#'                 facet.col = 4,
#'                 y.exp.up = 0.1,
#'                 y.axis.title = bquote(~2^-dCt))
#'
#' @importFrom base print paste
#' @importFrom dplyr filter
#' @import ggsignif
#' @import ggplot2
#' @import tidyverse
#'
results_boxplot <- function(data,
                            coef = 1.5,
                            sel.Gene = "all",
                            by.group = TRUE,
                            signif.show = TRUE,
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
                            dpi = 600, width = 15, height = 15,
                            name.tiff = "results_boxplot"){

  data <- pivot_longer(data, !c(Sample, Group), names_to = "Gene" , values_to = "value")

  if (sel.Gene[1] == "all"){
    data <- data

  } else {

    data <- filter(data, Gene %in% sel.Gene)
  }

  if (by.group == TRUE){

    box_results <- ggplot(data, aes(x = Gene, y = value)) +
      geom_boxplot(aes(fill = Group), coef = 1.5) +
      scale_fill_manual(values = c(colors)) +
      labs(fill = legend.title, title = plot.title)

    if (signif.show == TRUE){
      label.height <- data %>%
        group_by(Gene) %>%
        summarise(height = max(value), .groups = "keep")

      data.label <- data.frame(matrix(nrow = length(unique(label.height$Gene)), ncol = 4))
      rownames(data.label) <- label.height$Gene
      colnames(data.label) <- c("x", "xend", "y", "annotation")


      if (faceting == TRUE){

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
        geom_signif(stat = "identity",
                    data = data.label,
                    aes(x = x,
                        xend = xend,
                        y = y,
                        yend = y,
                        annotation = annotation),
                    color = "black",
                    manual = TRUE) +
        scale_y_continuous(expand = expansion(mult = c(y.exp.low, y.exp.up)))

      }

    if (faceting == TRUE) {
      box_results <- box_results +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(vars(Gene), scales = "free", nrow = facet.row, ncol = facet.col)
    }
      } else {

    box_results <- ggplot(data, aes(x = Gene, y = value)) +
      geom_boxplot(coef = coef, fill = colors[1]) +
      labs(title = plot.title)

    if (faceting == TRUE) {
      box_results <- box_results +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(vars(Gene), scales = "free", nrow = facet.row, ncol = facet.col)
    }
  }

  box_results <- box_results +
    guides(x =  guide_axis(angle = angle)) +
    xlab(x.axis.title) +
    ylab(y.axis.title) +
    theme_bw() +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour="black")) +
    theme(legend.text = element_text(size = legend.text.size, colour="black")) +
    theme(legend.title = element_text(size = legend.title.size, colour="black")) +
    theme(plot.title = element_text(size = plot.title.size)) +
    theme(panel.grid.major.x = element_blank())


  if (rotate == TRUE){
    box_results <- box_results +
      coord_flip()
  }

  if (add.mean == TRUE & by.group == TRUE){
    box_results <- box_results +
      stat_summary(aes(group = Group),
                   fun = mean,
                   position = position_dodge(width = .75),
                   geom = "point",
                   shape = 15,
                   size = add.mean.size,
                   color = add.mean.color)
  }
  if (add.mean == TRUE & by.group == FALSE){
    box_results <- box_results +
      stat_summary(fun = mean,
                   position = position_dodge(width = .75),
                   geom = "point",
                   shape = 15,
                   size = add.mean.size,
                   color = add.mean.color)
  }

  print(box_results)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), box_results, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
  return(box_results)
}







#' @title results_barplot
#'
#' @description
#' This function creates a barplot illustrating mean and sd values of each gene.
#' Faceting and adding custom labels of statistical significance is possible.
#' This function could be useful to present results for finally selected genes.
#'
#' @param data object returned from make_Ct_ready(), exp_Ct_dCt() or delta_Ct() functions.
#' @param sel.Gene character vector with names of genes to include, or "all" (default) to use all names of genes.
#' @param bar.width numeric: width of bars.
#' @param signif.show logical: if TRUE, labels for statistical significance will be added to the plot.
#' @param signif.labels character vector with statistical significance labels (ex. "ns.","***", etc.) with number
#' of elements equal to nimber of genes used for plotting.
#' The ggsignif package was used for convenient adding labels, but there is one tricky point:
#' the same elements of labels can not be handled by used package and must be different.
#' It could be achieved by adding symmetrically white spaces to repeated labels, ex. "ns.", " ns. ", "  ns.  ".
#' @param signif.length numeric: length of horizontal bars, values from 0 to 1.
#' @param signif.dist numeric: distance between errorbar and significance label.
#' Could be in y axis units (if `faceting` = TRUE) or fraction of y axis value reached by errorbar (mean + sd value) (if `faceting` = TRUE).
#' @param faceting logical: if TRUE (default), plot will be drawn with facets with free scales using ggplot2::facet_wrap() function (see its documentation for more details).
#' @param facet.row,facet.col integer: number of rows and columns to arrange facets.
#' @param y.exp.low,y.exp.up numeric: space between data on the plot and lower or upper axis. Useful to add extra space for statistical significance labels when `faceting` = TRUE.
#' @param angle integer: value of angle in which names of genes should be displayed. Default to 0.
#' @param rotate logical: if TRUE, boxplots will be arranged horizontally. Deafault to FALSE.
#' @param colors character vector length of one (when by.group = FALSE) or two (when by.group = TRUE), containing colors for groups. Numbers of colors must be equal to number of groups. Default to c("#66c2a5", "#fc8d62").
#' @param x.axis.title character: title of x axis. Default to "Gene".
#' @param y.axis.title character: title of y axis. Default to "value".
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 10.
#' @param legend.title character: title of legend. Default to "Group".
#' @param legend.title.size integer: font size of legend title. Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "results_boxplot".
#'
#' @return Object with boxplot illustrating distribution of data for selected genes. Created plot will be also displayed on graphic device.
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
#' data.dCt.exp <- exp_delta_Ct(data.dCt)
#' data.dCt.expF <- filter_transformed_data(data.dCt.exp, remove.Sample = c("Control11"))
#' results_barplot(data.dCt.exp,
#'                 sel.Gene = c("Gene1","Gene16","Gene19","Gene20"),
#'                 signif.labels = c("****","*","***"," * "),
#'                 angle = 30,
#'                 signif.dist = 1.05,
#'                 facet.row = 1,
#'                 facet.col = 4,
#'                 y.exp.up = 0.1,
#'                 y.axis.title = bquote(~2^-dCt))
#'
#' @importFrom base print paste
#' @importFrom dplyr filter
#' @import ggsignif
#' @import ggplot2
#' @import tidyverse
#'
results_barplot <- function(data,
                            sel.Gene = "all",
                            bar.width = 0.8,
                            signif.show = TRUE,
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
                            dpi = 600, width = 15, height = 15,
                            name.tiff = "results_barplot"){

  data <- pivot_longer(data, !c(Sample, Group), names_to = "Gene" , values_to = "value")

  if (sel.Gene[1] == "all"){
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
      geom_errorbar(aes(group = Group, ymin=mean, ymax=mean+sd), width=.2,
                    position=position_dodge(.9)) +
      geom_col(aes(fill = Group, group = Group), position=position_dodge(.9), width = bar.width, color = "black") +
      scale_fill_manual(values = c(colors)) +
      guides(x =  guide_axis(angle = angle)) +
      xlab(x.axis.title) +
      ylab(y.axis.title) +
      labs(fill = legend.title, title = plot.title) +
      theme_bw() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(plot.title = element_text(size = plot.title.size)) +
      theme(panel.grid.major.x = element_blank())

    if (faceting == TRUE) {
      bar_results <- bar_results +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(vars(Gene), scales = "free", nrow = facet.row, ncol = facet.col)
    }

    if (signif.show == TRUE) {
      label.height <- data.mean %>%
        mutate(max = mean + sd) %>%
        group_by(Gene) %>%
        summarise(height = max(max, na.rm = TRUE), .groups = "keep")

      data.label <- data.frame(matrix(nrow = length(unique(data.mean$Gene)), ncol = 4))
      rownames(data.label) <- unique(data.mean$Gene)
      colnames(data.label) <- c("x", "xend", "y", "annotation")

      if (faceting == TRUE) {
        data.label$x <- rep(1 - signif.length, length(unique(data.mean$Gene)))
        data.label$xend <- rep(1 + signif.length, length(unique(data.mean$Gene)))
        data.label$y <- label.height$height * signif.dist

        } else {

      data.label$x <- (1:length(unique(data.mean$Gene))) - signif.length
      data.label$xend <- (1:length(unique(data.mean$Gene))) + signif.length
      data.label$y <- label.height$height + signif.dist
    }

      data.label$annotation <- signif.labels
      data.label$Gene <- rownames(data.label)

      bar_results <- bar_results +
        geom_signif(stat = "identity",
                    data = data.label,
                    aes(x = x,
                      xend = xend,
                      y = y,
                      yend = y,
                      annotation = annotation),
                   color = "black",
                   manual = TRUE) +
        scale_y_continuous(expand = expansion(mult = c(y.exp.low, y.exp.up)))
    }

  if (rotate == TRUE){
    bar_results <- bar_results +
      coord_flip()
  }

  print(bar_results)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), bar_results, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
  return(bar_results)
}






#' @title RQ_ddCt
#'
#' @description
#' Performs relative quantification of gene expression using 2^(-ddCt) method.
#'
#' @details
#' This function calculates:
#' * Means (return in columns with "_mean" pattern) and standard deviations (return in columns with "_sd" pattern) of delta Ct values of analyzed genes across compared groups.
#' * Normality tests (Shapiro_Wilk test) of delta Ct values of analyzed genes across compared groups and returned p values in columns with "_norm_p" pattern.
#' * Differences in mean delta Ct values of genes between compared groups, obtaining delta delta Ct values (in returned "ddCt" column).
#' * Fold change values (return in "FCh" column) together with log10 Fold change values (return in "log10FCh" column).
#'   Fold change values are calculated for each gene by exponentiating ddCt values using 2^-ddCt formula.
#' * Statistical testing of delta delta Ct values (differences between study group and reference group).
#'   Student's t test and Mann-Whitney U test are implemented and resulted statistics (in column with "_test_stat" pattern) and p values (in column with "_test_p" pattern) are returned.
#'
#' @param data data object returned from delta_Ct() function.
#' @param group.study character: name of study group (group of interest).
#' @param group.ref character: name of reference group.
#' @param do.tests logical: if TRUE, statistical significance of delta delta Ct values between compared groups will be calculated using Student's t test and Mann-Whitney U test. Default to TRUE.
#' @param save.to.txt logical: if TRUE, returned table with results will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension. Default to "RQ_expCt_results".
#
#' @return Data.frame with relative quantification results.
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
#' RQ.ddCt <- RQ_ddCt(data.dCt, "Disease", "Control")
#' head(RQ.ddCt)
#'
#' @importFrom base as.data.frame as.factor mean
#' @importFrom stats sd shapiro.test t.test
#' @importFrom coin wilcox_test pvalue statistic
#' @importFrom utils write.table
#' @importFrom dplyr filter select rename_with full_join
#' @import tidyverse
#'
RQ_ddCt <- function(data,
                    group.study,
                    group.ref,
                    do.tests = TRUE,
                    save.to.txt = FALSE,
                    name.txt = "ddCt_RQ_results"){

  data_slim <- data %>%
    filter(Group == group.study | Group == group.ref) %>%
    pivot_longer(cols = -c(Group, Sample), names_to = "Gene", values_to = "dCt")

  data_slim$Group <- as.factor(data_slim$Group)

  data_ddCt <- data_slim %>%
    group_by(Group, Gene) %>%
    summarise(ddCt = mean(dCt, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = ddCt) %>%
    mutate(ddCt = .data[[group.study]] - .data[[group.ref]]) %>%
    mutate(FCh = 2^-ddCt) %>%
    rename_with(~paste0(.x, "_mean", recycle0 = TRUE), all_of(c(group.study, group.ref)))

  data_ddCt_sd <- data_slim %>%
    group_by(Group, Gene) %>%
    summarise(dCt_sd = sd(dCt, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = dCt_sd) %>%
    rename_with(~paste0(.x, "_sd", recycle0 = TRUE), all_of(c(group.study, group.ref)))

  if (do.tests == TRUE){

    data_ddCt_norm <- data_slim %>%
      group_by(Group, Gene) %>%
      summarise(shap_wilka_p = shapiro.test(dCt)$p.value, .groups = "keep") %>%
      pivot_wider(names_from = Group, values_from = shap_wilka_p) %>%
      full_join(data_ddCt, by = c("Gene")) %>%
      rename_with(~paste0(.x, "_norm_p", recycle0 = TRUE), all_of(c(group.study, group.ref)))

    data_ddCt_tests <- data_slim %>%
      group_by(Gene) %>%
      summarise(t_test_p = t.test(dCt ~ Group, alternative = "two.sided")$p.value,
                t_test_stat = t.test(dCt ~ Group, alternative = "two.sided")$statistic,
                MW_test_p = pvalue(wilcox_test(dCt ~ Group)),
                MW_test_stat = statistic(wilcox_test(dCt ~ Group)), .groups = "keep")

    data_ddCt_norm_tests <- full_join(data_ddCt_norm, data_ddCt_tests, by = c("Gene"))
    data_ddCt_results <- full_join(data_ddCt_norm_tests, data_ddCt_sd, by = c("Gene"))
    data_ddCt_results <- select(data_ddCt_results, Gene, ends_with("_mean"), ends_with("_sd"), everything())

    return(data_ddCt_results)

    } else{

      data_ddCt_results <- full_join(data_ddCt, data_ddCt_sd, by = c("Gene"))
      data_ddCt_results <- rename_with(data_ddCt_results, ~paste0(.x, "_mean", recycle0 = TRUE), all_of(c(group.study, group.ref)))
      data_ddCt_results <- select(data_ddCt_results, Gene, ends_with("_mean"), ends_with("_sd"), everything())

    return(data_ddCt_results)
  }
  if (save.to.txt == TRUE){
    write.table(data_ddCt_norm_tests, paste(name.txt,".txt", sep = ""))
  }
}




#' @title RQ_plot
#'
#' @description
#' This function creates barplot illustrating fold change (when 2^-Ct or 2^-dCt methods are used) or RQ (when 2^-ddCt method is used) values.
#' with indicating of statistical significance.
#'
#' @param data object returned from RQ_exp_Ct_dCt() or RQ_ddCt() functions.
#' @param use.p logical: if TRUE, bars of statistically significant genes will be distinguished by colors.
#' @param mode character: which p value should be used? One of the "t" (p values from Student's t test),
#' "mw" (p values from Mann-Whitney U test), "depends" (if data in both compared groups were considered as derived from normal distribution (p value from Shapiro_Wilk test > 0.05) - p
#' values from Student's t test will be used for significance assessment, otherwise p values from Mann-Whitney U test will be used for significance assessment).
#' There is one more option, if user intend to use another p values, ex. obtained from other statistical test,
#' a mode parameter could be set to "user". In this situation, before run RQ_plot function, user should to prepare
#' data.frame object names "user"with two columns, one named "Gene" with Gene names and second with p values. The order of columns must be kept as described.
#' @param p.threshold numeric: threshold of p values for statistical significance.
#' @param use.log10FCh logical: if TRUE, the criterion of fold change will be also used for significance assessment of genes.
#' @param log10FCh.threshold numeric: threshold of log10 fold change values used for significance assessment of genes.
#' @param sel.Gene character vector with names of genes to include, or "all" (default) to use all names of genes.
#' @param bar.width numeric: width of bars.
#' @param signif.show logical: if TRUE, labels for statistical significance will be added to the plot.
#' @param signif.labels character vector with statistical significance labels (ex. "ns.","***", etc.) with number
#' of elements equal to nimber of genes used for plotting.
#' The ggsignif package was used for convenient adding labels, but there is one tricky point:
#' the same elements of labels can not be handled by used package and must be different.
#' It could be achieved by adding symmetrically white spaces to repeated labels, ex. "ns.", " ns. ", "  ns.  ".
#' @param signif.length numeric: length of horizontal bars, values from 0 to 1.
#' @param signif.dist numeric: distance between errorbar and significance label.
#' Could be in y axis units (if `faceting` = TRUE) or fraction of y axis value reached by errorbar (mean + sd value) (if `faceting` = TRUE).
#' @param y.exp.low,y.exp.up numeric: space between data on the plot and lower or upper axis. Useful to add extra space for statistical significance labels when `faceting` = TRUE.
#' @param angle integer: value of angle in which names of genes should be displayed. Default to 0.
#' @param rotate logical: if TRUE, bars will be arranged horizontally. Deafault to FALSE.
#' @param colors character vector length of one (when use.p = FALSE) or two (when use.p = TRUE), containing colors for significant and no significant genes.
#' @param x.axis.title character: title of x axis. Default to "Gene".
#' @param y.axis.title character: title of y axis. Default to "value".
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 10.
#' @param legend.title character: title of legend. Default to "Group".
#' @param legend.title.size integer: font size of legend title. Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "results_boxplot".
#'
#' @return List containing object with barplot and dataframe with results. Created plot will be also displayed on graphic device.
#' @export
#'
#' @examples
#' library(signif)
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.ready <- make_Ct_ready(data.CtF, imput.by.mean.within.groups = TRUE)
#' data.dCt <- delta_Ct(data.CtF.ready, ref = "Gene8")
#' data.dCt.exp <- exp_delta_Ct(data.dCt)
#' data.dCt.expF <- filter_transformed_data(data.dCt.exp, remove.Sample = c("Control11"))
#' RQ.ddCt <- RQ_ddCt(data.dCt.expF, "Disease", "Control")
#'
#' signif.labels <- c("****","**","ns."," ns. ","  ns.  ","   ns.   ","    ns.    ","     ns.     ","      ns.      ","       ns.       ",        "ns.        ","         ns.         ","          ns.          ","***")
#' RQ.plot <- RQ_plot(RQ.ddCt,
#'                    mode = "depends",
#'                    use.log10FCh = TRUE,
#'                    log10FCh.threshold = 0.30103,
#'                    signif.labels = signif.labels,
#'                    angle = 30)
#' head(RQ.plot[[2]])
#'
#' # with user p values - in this example used p values are calculated using stats::wilcox.test() function:
#' user <- data.dCt %>%
#' pivot_longer(cols = -c(Group, Sample), names_to = "Gene", values_to = "dCt") %>%
#'   group_by(Gene) %>%
#'   summarise(MW_test_p = wilcox.test(dCt ~ Group)$p.value, .groups = "keep")
#'
#' RQ.plot <- RQ_plot(RQ.ddCt,
#'                    mode = "user",
#'                    use.log10FCh = TRUE,
#'                    log10FCh.threshold = 0.30103,
#'                    signif.labels = signif.labels,
#'                    angle = 30)
#' head(RQ.plot[[2]])
#
#' @importFrom base print paste colnames factor
#' @importFrom dplyr filter
#' @importFrom stats reorder
#' @import ggsignif
#' @import ggplot2
#' @import tidyverse
#'
RQ_plot <- function(data,
                    use.p = TRUE,
                    mode,
                    p.threshold = 0.05,
                    use.log10FCh = FALSE,
                    log10FCh.threshold = 0,
                    sel.Gene = "all",
                    bar.width = 0.8,
                    signif.show = TRUE,
                    signif.labels,
                    signif.length = 0.2,
                    signif.dist = 0.1,
                    y.exp.low = 0.1,
                    y.exp.up = 0.1,
                    angle = 0,
                    rotate = FALSE,
                    colors = c("#66c2a5", "#fc8d62"),
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
                    dpi = 600, width = 15, height = 15,
                    save.to.tiff = FALSE,
                    name.tiff = "RQ_plot"){

  if (sel.Gene[1] != "all"){
    data <- filter(data, Gene %in% sel.Gene)
  }

  if (use.p == TRUE){
    if (mode == "t"){
      data$p.used <- data$t_test_p
    }
    if (mode == "mw"){
      data$p.used <- data$MW_test_p
    }
    if (mode == "depends"){
      data <- ungroup(data)
      vars <- colnames(select(data, ends_with("norm_p")))
      data <- mutate(data, test.for.comparison = ifelse(.data[[vars[[1]]]] >= 0.05 & .data[[vars[[2]]]] >= 0.05, "t.student's.test", "Mann-Whitney.test"))
      data <- mutate(data, p.used = ifelse(test.for.comparison == "t.student's.test", t_test_p, MW_test_p))
    }
    if (mode == "user"){
      colnames(user) <- c("Gene","p.used")
      data <- full_join(data, user, by = c("Gene"))
    }
    if (use.log10FCh == TRUE){
      data <- mutate(data, `Selected as significant?` = ifelse(p.used > p.threshold, yes = "No",
                                                                 no = ifelse(abs(log10(FCh)) <  log10FCh.threshold, "No", "Yes")))
      data$`Selected as significant?` <- factor(data$`Selected as significant?`, levels = c("Yes", "No"))

          } else {

    data <- mutate(data, `Selected as significant?` = ifelse(p.used > p.threshold, yes = "No (p > 0.05)",  no = "Yes (p <= 0.05)"))
    data$`Selected as significant?` <- factor(data$`Selected as significant?`, levels = c("Yes (p <= 0.05)", "No (p > 0.05)"))
          }

    RQ <- ggplot(data, aes(x = reorder(Gene, -FCh), y = log10(FCh))) +
      geom_col(aes(fill = `Selected as significant?`, group = `Selected as significant?`), width = bar.width) +
      scale_fill_manual(values = c(colors)) +
      labs(fill = legend.title, title = plot.title)

    } else {

      RQ <- ggplot(data, aes(x = reorder(Gene, -FCh), y = log10(FCh))) +
        geom_col(width = bar.width, fill = colors[1]) +
        labs(title = plot.title)
  }


  RQ <- RQ +
    guides(x =  guide_axis(angle = angle)) +
    xlab(x.axis.title) +
    ylab(y.axis.title) +
    theme_bw() +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour="black")) +
    theme(legend.text = element_text(size = legend.text.size, colour="black")) +
    theme(legend.title = element_text(size = legend.title.size, colour="black")) +
    theme(plot.title = element_text(size = plot.title.size)) +
    theme(panel.grid.major.x = element_blank()) +
    geom_hline(yintercept = 0, linewidth = 0.4)


  if (signif.show == TRUE) {

    data.label <- data.frame(matrix(nrow = nrow(data), ncol = 4))
    colnames(data.label) <- c("x", "xend", "y", "annotation")
    data <- arrange(data, desc(FCh))
    data.label$x <- (1:nrow(data)) - signif.length
    data.label$xend <- (1:nrow(data)) + signif.length

    data.label <- mutate(data.label, y = ifelse(log10(data$FCh) > 0,
                                                log10(data$FCh) + signif.dist,
                                                log10(data$FCh) - signif.dist))

    data.label$annotation <- signif.labels

    RQ <- RQ +
      geom_signif(stat = "identity",
                  data = data.label,
                  aes(x = x,
                      xend = xend,
                      y = y,
                      yend = y,
                      annotation = annotation),
                  color = "black",
                  manual = TRUE) +
      scale_y_continuous(expand = expansion(mult = c(y.exp.low, y.exp.up)))

  }

  if (rotate == TRUE){
    RQ <- RQ +
      coord_flip()
  }

  print(RQ)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), RQ, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }

  return(list(RQ, data))
}







#' @title ROCh
#'
#' @description
#' This function performs Receiver Operating Characteristic (ROC) analysis of samples of genes based on the expression data.
#' This analysis could be useful to further examine performance of samples classification into groups using gene expression data.
#'
#' @param data object returned from exp_Ct_dCt() or delta_Ct() functions.
#' @param sel.Gene character vector with names of genes to include, or "all" (default) to use all names of genes.
#' @param groups character vector length of two with names of compared groups
#' @param panels.row,panels.col integer: number of rows and columns to arrange panels with plots.
#' @param text.size numeric: size of text on the plot. Default to 1.1.
#' @param print.auc logical: if TRUE, AUC values will be added to the plot. Default to TRUE.
#' @param print.auc.size numeric: size of AUC text on the plot. Default to 0.8.
#' @param plot.title character: title of the plot.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "ROC_plot".
#' @param save.to.txt logical: if TRUE, returned table with results will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension. Default to "ROC_results".
#'
#' @return Data frame with ROC parameters including AUC, threshold, specificity, sensitivity, accuracy,
#' positive predictive value, negative predictive value, and Youden's J statistic.
#' Plot with ROC curves could be saved to .tiff file (will not be displayed on the graphic device).
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
#' data.dCt.expF <- filter_transformed_data(data.dCt.exp, remove.Sample = c("Control11"))
#' roc_parameters <- ROCh(data.dCt, sel.Gene = c("Gene1","Gene16","Gene19","Gene20"),
#'                        groups = c("Disease","Control"),
#'                        panels.row = 2,
#'                        panels.col = 2)
#'
#' @importFrom base plot colnames as.data.frame ncol matrix
#' @importFrom dplyr select filter
#' @importFrom utils write.table
#' @import tidyverse
#' @import pROC
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
                dpi = 600, width = 15, height = 15,
                name.tiff = "ROC_plot",
                save.to.txt = FALSE,
                name.txt = "ROC_results"){

  data <- filter(data, Group %in% groups)

  if (sel.Gene[1] != "all"){
    data <- data[, colnames(data) %in% c("Group", "Sample", sel.Gene)]
  }

  roc_param <- as.data.frame(matrix(nrow = ncol(data)-2, ncol = 9))
  colnames(roc_param) <- c("Gene","Threshold", "Specificity", "Sensitivity", "Accuracy", "ppv", "npv", "youden", "AUC")
  roc_param$Gene <- colnames(data)[-c(1:2)]

  for (x in 1:nrow(roc_param)){
    myproc <- roc(response = data$Group, predictor = as.data.frame(data)[ ,x+2], levels = c(groups),
                  smooth = FALSE, auc = TRUE, plot=FALSE, ci=TRUE, of = "auc", quiet = TRUE)
    parameters <- coords(myproc, "best", ret = c("threshold", "specificity", "sensitivity","accuracy", "ppv", "npv", "youden"))
    roc_param[x,2:8] <- parameters
    roc_param[x,9] <- myproc$auc
    roc_param[x,1] <- colnames(data)[x+2]

    if (nrow(parameters) > 1){
      cat('Warning: ',colnames(data)[x+2],'has more than 1 threshold value for calculated Youden’s J statistic.\n')
    } else {}
  }

  if (save.to.tiff == TRUE){
    tiff(paste(name.tiff,".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
    par(mfrow = c(panels.row, panels.col))

    for (x in 1:nrow(roc_param)){
      myproc <- roc(response = data$Group, predictor = as.data.frame(data)[ ,x+2], levels = c(groups),
                    smooth = FALSE, auc = TRUE, plot=FALSE, ci=TRUE, of = "auc", quiet = TRUE)
      plot.roc(myproc, main = roc_param$Gene[x],
               smooth = FALSE, cex.axis = text.size, cex.lab = text.size, identity.lwd = 2,
               plot = TRUE, percent = TRUE, print.auc = print.auc, print.auc.x = 0.85, print.auc.y = 0.1, print.auc.cex = print.auc.size)
    }
    dev.off()
  }

  if (save.to.txt == TRUE){
    write.table(roc_param, paste(name.txt,".txt", sep = ""))
  }

  return(roc_param)
}







#' @title log_reg
#'
#' @description
#' This function performs logistic regression analysis, computes odd ratio values and presents them graphically.
#'
#' @param data object returned from exp_Ct_dCt() or delta_Ct() functions.
#' @param sel.Gene character vector with names of genes to include, or "all" (default) to use all names of genes.
#' @param group.study character: name of study group (group of interest).
#' @param group.ref character: name of reference group.
#' @param centerline numeric: position of vertical centerline on the plot, if logaxis = TRUE, centerline should be set to 0, otherwise to 1 (default).
#' @param ci numeric: confidence level used for computation of confidence interval. Default to 0.95.
#' @param log.axis logical: if TRUE, axis with odds ratio values will be in log10 scale. If axis is in logarithmic scale, centerline parameter should be set to 0. Default ot FALSE.
#' @param x.axis.title character: title of x axis. Default to "Gene".
#' @param y.axis.title character: title of y axis. Default to "value".
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 10.
#' @param legend.title character: title of legend. Default to "Group".
#' @param legend.title.size integer: font size of legend title. Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top" (default), "right", "bottom", "left", or "none" (no legend).
#' See description for legend.position in ggplot2::theme() function.
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "OR_plot".
#' @param save.to.txt logical: if TRUE, returned table with results will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension. Default to "OR_results".
#'
#' @return List containing object with barplot and dataframe with results. Created plot will be also displayed on graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(oddsratio)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Gene = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' log.reg.results <- log_reg(data.dCt,
#'                             sel.Gene = c("Gene1","Gene16","Gene19","Gene20"),
#'                             group.study = "Disease",
#'                             group.ref = "Control")
#
#' @importFrom base print paste colnames ncol lapply as.data.frame as.vrctor names summary data.frame
#' @importFrom dplyr filter
#' @importFrom stats coef
#' @importFrom utils write.table
#' @import ggplot2
#' @import tidyverse
#' @import oddsratio
#'
log_reg <- function(data,
                    sel.Gene = "all",
                    group.study,
                    group.ref,
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
                    dpi = 600, width = 15, height = 15,
                    name.tiff = "OR_plot",
                    save.to.txt = FALSE,
                    name.txt = "OR_results"){

  data <- filter(data, Group %in% c(group.study, group.ref))

  if (sel.Gene[1] != "all"){
    data <- data[, colnames(data) %in% c("Group", "Sample", sel.Gene)]
  }

  data <- mutate(data, Group_num = ifelse(Group == group.study, 0, 1))
  n.genes <- ncol(data)-3
  list.models <- lapply(data[3:(n.genes+2)], function(x) glm(data$Group_num ~ x, data = data, family = binomial))
  list.CI <- lapply(names(list.models)[1:n.genes], function(x) or_glm(data = data,
                                                           model = list.models[[x]],
                                                           incr = list(x = 1),
                                                           ci = ci))
  data.CI <- as.data.frame(matrix(ncol = 8, nrow = n.genes))
  colnames(data.CI) <- c("Gene", "oddsratio", "CI_low", "CI_high", "Intercept", "coeficient","p_intercept","p_coef")

  for (x in 1:n.genes){
    data.CI$Gene <- names(list.models)
    data.CI[x,2:4] <- as.vector(list.CI)[[x]][2:4]
    data.CI[x,5:6] <- list.models[[x]]$coefficients
    data.CI[x,7:8] <- coef(summary(list.models[[x]]))[,4]
  }
  od_df <- data.frame(yAxis = 1:nrow(data.CI),
                        boxOdds = data.CI$oddsratio,
                        boxCILow = data.CI$CI_low,
                        boxCIHigh = data.CI$CI_high,
                        boxLabels = data.CI$Gene,
                        p = data.CI$p_coef)

  odd.ratio <- ggplot(od_df, aes(x = boxOdds, y = boxLabels, label = boxOdds)) +
    geom_vline(aes(xintercept = centerline), linewidth = .25, linetype = "dashed") +
    geom_errorbarh(aes(xmax = boxCIHigh, xmin = boxCILow), linewidth = .5, height = .2) +
    geom_point(aes(color = p), size = 3.5) +
    scale_color_continuous(type = "viridis") +
    geom_text(aes(label = boxOdds), hjust=0.5, vjust = -1, size = 3) +
    xlab(x.axis.title) + ylab(y.axis.title) +
    labs(color = legend.title, title = plot.title) +
    theme_bw() +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour="black")) +
    theme(legend.text = element_text(size = legend.text.size, colour="black")) +
    theme(legend.title = element_text(size = legend.title.size, colour="black")) +
    theme(plot.title = element_text(size = plot.title.size))

  if (log.axis == TRUE){
    odd.ratio <- odd.ratio +
      scale_x_log10()
  }

  print(odd.ratio)

   if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), odd.ratio, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }

  if (save.to.txt == TRUE){
    write.table(data.CI, paste(name.txt,".txt", sep = ""))
  }

  return(list(odd.ratio, data.CI))
}
