
#' @title read_Ct_wide
#'
#' @description
#' Enables to import Ct dataset in a wide-format table with sample names given in columns.
#'
#' @details
#' This function needs two tables to import: a wide-format table with Ct values and design file (see parameters `path.Ct.file` and `path.design.file` for further details regarding tables structure).
#' Subsequently merges both files to return long-format table ready for analysis.
#' All parameters must be specified, there is no default values.
#'
#' @param path.Ct.file path to wide-format table in .txt format containing Ct values, target names in the first column, and
#' sample names in the first row. In other words, this table should contain targets by rows and samples by columns.
#' @param path.design.file path to .txt file with two columns: column named "Sample" with names of samples
#' and column named "Group" with names of groups assigned to samples. Names of samples in this file
#' should correspond to the names of columns in file with Ct values.
#' @param sep character of a field separator in both imported files.
#' @param dec character used for decimal points in Ct values.
#'
#' @return Data.frame in long format ready to analysis.
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

  colnames(data_wide)[1] <- "Target"
  data_wide <- mutate(data_wide, across(everything(), as.character))
  data_slim <- pivot_longer(data_wide, -Target, names_to = "Sample", values_to = "Ct")
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
#' sample names, target names, Ct values and group names (those columns will be imported by this function).
#' Imported table could also contain a column with flag information, which could be optionally imported (see add.col.Flag and col.Flag parameters).
#'
#' @param sep character of a field separator in imported file.
#' @param dec character used for decimal points in Ct values.
#' @param skip integer: number of lines of the data file to skip before beginning to read data. Default to 0.
#' @param col.Sample integer: number of column with sample names.
#' @param col.Target integer: number of column with target names.
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
#'                         add.column.Flag = TRUE, column.Sample = 1, column.Target = 2,
#'                         column.Ct = 5, column.Group = 9, column.Flag = 4)
#' str(data.Ct)
#'
#' data.Ct <- mutate(Ct, Flag = ifelse(Flag < 1, "Undetermined", "OK"))
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
                         column.Target,
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
  data <- data[ ,c(column.Sample, column.Target, column.Ct, column.Group)]
  colnames(data) <- c("Sample", "Target", "Ct", "Group")
  }

  if (add.column.Flag == TRUE){
    data <- data[ ,c(column.Sample, column.Target, column.Ct, column.Group, column.Flag)]
    colnames(data) <- c("Sample", "Target", "Ct", "Group", "Flag")
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
#' or data frame containing column named "Sample" with sample names, column named "Target" with target names,
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
                                      axis.title.size = 12,
                                      axis.text.size = 10,
                                      plot.title = "",
						                          plot.title.size = 14,
                                      legend.title = "Reliable Ct value?",
                                      legend.title.size = 12,
                                      legend.text.size = 12,
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
  cat("Returned table contains number of targets retained (Yes) or not retained (No) in samples after reliability assessment.")

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

  return(list(barplot.samples, tab))
}




#' @title control_Ct_barplot_target
#'
#' @description
#' Target-wide control of raw Ct values across groups by illustrating numbers of Ct values labeled as reliable or not by using reliability criteria (see function parameters).
#'
#' @details
#' This function does not perform data filtering, but only numbers Ct values labeled as reliable or not and presents them graphically.
#' Could be useful to identify targets with low number of reliable Ct values.
#'
#' @param data object returned from read_Ct_long() or read_Ct_wide() function,
#' or data frame containing column named "Sample" with sample names, column named "Target" with target names,
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
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension.  Default to "Ct_control_barplot_for_targets".
#'
#' @return List containing plot and table with numbers of reliable and not reliable Ct values in targets.
#' Additional information about returned table is also printed, it could help user to properly interpret returned table.
#' Plot will be displayed on graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' target.Ct.control <- control_Ct_barplot_target(data.Ct)
#' target.Ct.control[[2]]
#'
#' @importFrom base as.numeric as.data.frame cat print table paste sum colnames
#' @importFrom dplyr mutate arrange filter
#' @import ggplot2
#' @import tidyverse
#'
control_Ct_barplot_target <- function(data,
                                      flag.Ct = "Undetermined",
                                      maxCt = 35,
                                      flag = "Undetermined",
                                      colors = c("#66c2a5", "#fc8d62"),
                                      x.axis.title = "",
                                      y.axis.title = "Number",
                                      axis.title.size = 12,
                                      axis.text.size = 10,
                                      legend.title = "Reliable Ct value?",
                                      legend.title.size = 12,
                                      legend.text.size = 12,
                                      legend.position = "top",
                                      plot.title = "",
						              	          plot.title.size = 14,
                                      save.to.tiff = FALSE,
                                      dpi = 600, width = 15, height = 15,
                                      name.tiff = "Ct_control_barplot_for_targets"){

  data$Ct[data$Ct == flag.Ct] <- 100
  data$Ct <- as.numeric(data$Ct)

  if(sum(colnames(data) %in% "Flag") > 0){
    data <- mutate(data, Reliable = ifelse(Ct > maxCt | Flag == flag, yes = "No",  no = "Yes"))
  } else {
    data <- mutate(data, Reliable = ifelse(Ct > maxCt , yes = "No",  no = "Yes"))
  }

    bar <- as.data.frame(table(data$Reliable, data$Target, data$Group))
    order <- arrange(filter(bar, Var1 == "Yes"), Freq)$Var2
    cat("Returned table contains number of samples retained (Yes) or not retained (No) for each target after filtering.")

    barplot.targets <- ggplot(bar, aes(x = reorder(Var2, desc(Freq)), y = Freq, fill = Var1)) +
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

    print(barplot.targets)

    if (save.to.tiff == TRUE){
      ggsave(paste(name.tiff, ".tiff", sep = ""), barplot.targets, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
    }

    tab <- table(data$Reliable, data$Target, data$Group)

    return(list(barplot.targets, tab))

}




#' @title samples_to_remove
#'
#' @description
#' Indicates samples with the amount of unreliable Ct values higher than specified fraction.
#' Could be useful for identification of samples with specified fraction of unreliable Ct values.
#' Samples with high amount of unreliable data could be considered to remove from the dataset.
#'
#' @param data object returned from control_Ct_barplot_sample() function.
#' @param fraction numeric: a threshold fraction (ranged from 0 to 1),
#' samples with the higher fraction of Ct values labeled as unreliable will be returned. Default to 0.5.
#' @return Vector with names of samples.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' sample.Ct.control <- control_Ct_barplot_sample(data.Ct)
#' samples.to.remove <- samples_to_remove(sample.Ct.control[[2]], 0.4)
#' samples.to.remove
#'
#' @importFrom base as.data.frame as.vector
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr mutate filter
#' @import tidyverse
#'
samples_to_remove <- function(data, fraction = 0.5){

  data1 <- data %>%
    as.data.frame() %>%
    pivot_wider(names_from = Var1, values_from = Freq) %>%
    mutate(rem = ifelse(No >= ((No + Yes) * fraction), "remove", "leave"))

  to.remove <- filter(data1, rem == "remove")$Var2

  return(as.vector(to.remove))
}



#' @title targets_to_remove
#'
#' @description
#' Indicates targets with the amount of unreliable Ct values higher than specified fraction in at least one of the compared groups.
#' Could be useful for identification of targets with specified fraction of unreliable Ct values.
#' Targets with high amount of unreliable data could be considered to remove from the data.
#'
#' @param data object returned from control_Ct_barplot_target() function.
#' @param fraction numeric: a threshold fraction (ranged from 0 to 1), targets with the higher fraction of Ct values labeled as unreliable in at least one of the compared groups will be returned. Default to 0.5.
#' @param groups character vector length of two, containing names of compared groups.
#'
#' @return Vector with names of targets.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' target.Ct.control <- control_Ct_barplot_target(data.Ct)
#' targets.to.remove <- targets_to_remove(target.Ct.control[[2]], fraction = 0.5, groups = c("Disease", "Control"))
#' targets.to.remove
#'
#' @importFrom base as.data.frame as.vector unique
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr mutate filter
#' @import tidyverse
#'
targets_to_remove <- function(data, fraction = 0.5, groups){

  data1 <- data %>%
    as.data.frame() %>%
    filter(Var3 == groups[1]) %>%
    pivot_wider(names_from = Var1, values_from = Freq) %>%
    mutate(rem = ifelse(No >= ((No + Yes) * fraction), "remove", "leave"))

  to.remove1 <- filter(data1, rem == "remove")$Var2

  data2 <- data %>%
    as.data.frame() %>%
    filter(Var3 == groups[2]) %>%
    pivot_wider(names_from = Var1, values_from = Freq) %>%
    mutate(rem = ifelse(No >= ((No + Yes) * fraction), "remove", "leave"))

  to.remove2 <- filter(data2, rem == "remove")$Var2

  return(unique(c(as.vector(to.remove1), as.vector(to.remove2))))
}




#' @title filter_Ct
#'
#' @description
#' Filters Ct data according to the used filtering criteria (see parameters).
#'
#' @param data object returned from read_Ct_long() or read_Ct_wide() function,
#' or data frame containing column named "Sample" with sample names, column named "Target" with target names,
#' column named "Ct" with raw Ct values, column named "Group" with group names.
#' Optionally, data frame could contain column named "Flag" with flag information (ex. "Undetermined" and "OK"),
#' which will be used for filtering.
#'
#' @param flag.Ct character of a flag used for undetermined Ct values, default to "Undetermined".
#' @param maxCt numeric, a maximum of Ct value allowed.
#' @param flag character: flag used in Flag column for values which should be filtered out, default to "Undetermined".
#' @param remove.Target character: vector with names of targets which should be removed from data
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
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
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
                      remove.Target = c(""),
                      remove.Sample = c(""),
                      remove.Group = c("")){

  data <- filter(data, Ct != flag.Ct)
  data$Ct <- as.numeric(data$Ct)
  data <- filter(data, Ct <= maxCt,
                        !Target %in% remove.Target,
                        !Sample %in% remove.Sample,
                        !Group %in% remove.Group)

  if(sum(colnames(data) %in% "Flag") > 0){
    data <- filter(data, !Flag %in% flag)
  }

  return(data)
}






#' @title Ct_for_control
#'
#' @description
#' This function collapses technical replicates (if present in data) by means counts and imputes missing data by means within groups (if so indicated).
#' These actions also prepare Ct data for control functions.
#'
#' @param data data object returned from read_Ct_long(), read_Ct_wide() or filter_Ct() function,
#' or data frame containing column named "Sample" with sample names, column named "Target" with target names,
#' column named "Ct" with raw Ct values (must be numeric), column named "Group" with group names. Any other columns could exist, but will not be used by this function.
#' @param imput.by.mean.within.groups logical: if TRUE, missing values will be imputed by means within groups. Default to TRUE.
#' @param save.to.txt logical: if TRUE, returned data will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension. Default to "data.exp.Ct".
#'
#' @return Data.frame with prepared data and printed information about number and percentage of missing values.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#'data.CtF.prepared <- Ct_for_control(data.CtF)
#'head(data.CtF.prepared)
#'
#' @importFrom base mean
#' @importFrom dplyr select summarise
#' @import tidyverse
#'
Ct_for_control <- function(data,
                   imput.by.mean.within.groups = TRUE,
                   save.to.txt = FALSE,
                   name.txt = "data.exp.Ct"){
  data <- data %>%
    group_by(Group, Target, Sample) %>%
    summarise(mean = mean(Ct, na.rm = TRUE), .groups = "keep") %>%
    as.data.frame()

  data_wide <- data %>%
    select(Group, Sample, Target, mean) %>%
    pivot_wider(names_from = Target, values_from = mean)

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







#' @title exp_Ct
#'
#' @description
#' This function collapses technical replicates (if present in data) by means, counts and imputes missing data by means within groups (if so indicated),
#' and finally exponentiates Ct values by using formula 2^(-Ct).
#'
#' @param data data object returned from read_Ct_long(), read_Ct_wide() or filter_Ct() function,
#' or data frame containing column named "Sample" with sample names, column named "Target" with target names,
#' column named "Ct" with raw Ct values (must be numeric), column named "Group" with group names. Any other columns could exist, but will not be used by this function.
#' @param imput.by.mean.within.groups logical: if TRUE, missing values will be imputed by means within groups. Default to TRUE.
#' @param save.to.txt logical: if TRUE, returned data will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension. Default to "data.exp.Ct".
#'
#' @return Data.frame with exponentiated Ct values and printed information about number and percentage of missing values.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.Ct.exp <- exp_Ct(data.CtF, imput.by.mean.within.groups = TRUE)
#' head(data.Ct.exp)
#'
#' @importFrom base as.data.frame as.vector sum is.na ncol nrow paste cat
#' @importFrom utils write.table
#' @importFrom dplyr mutate filter select
#' @import tidyverse
#'
exp_Ct <- function(data,
                  imput.by.mean.within.groups = TRUE,
                  save.to.txt = FALSE,
                  name.txt = "data.exp.Ct"){
  data <- data %>%
    group_by(Group, Target, Sample) %>%
    summarise(mean = base::mean(Ct, na.rm = TRUE), .groups = "keep") %>%
    as.data.frame()

  calcRQ <- function(x){
    x <- 2^-x
  }

  data_wide <- data %>%
    select(Group, Sample, Target, mean) %>%
    pivot_wider(names_from = Target, values_from = mean)

  nas <- sum(is.na(data_wide))
  percentage <- sum(is.na(data_wide))/((ncol(data_wide)-2)*nrow(data_wide))

  if (imput.by.mean.within.groups == TRUE){
    data_wide_imp_exp <- data_wide %>%
      group_by(Group) %>%
      mutate(across(where(is.numeric), ~ replace(., is.na(.), mean(., na.rm = TRUE)))) %>%
      mutate_at(vars(-c("Group","Sample")), calcRQ)

    cat("Data contained", nas, "missing values that constitute", round(percentage*100, 5), "percent of the total data.\n Missing values were imputed using means within compared groups.\n")

    if (save.to.txt == TRUE){
      write.table(as.data.frame(data_wide_imp), paste(name.txt,".txt", sep = ""))
    }
    return(data_wide_imp_exp)

  } else {

    data_wide_exp <- data_wide %>%
      mutate_at(vars(-c("Group","Sample")), calcRQ)

    cat("Data contains", nas, "missing values that constitute", round(percentage*100, 5), "percent of the total data.")

    if (save.to.txt == TRUE){
      write.table(as.data.frame(data_wide_exp), paste(name.txt,".txt", sep = ""))
    }
    return(data_wide_exp)
  }
}





#' @title RQ_exp_Ct
#'
#' @description
#' Performs relative quantification of gene expression using 2^(-Ct) and 2^(-dCt) methods.
#'
#' @details
#' This function calculates:
#' * Means (returned in columns with "_mean" pattern) and standard deviations (returned in columns with "_sd" pattern) of exponentiated Ct or dCt values of analyzed targets across compared groups.
#' * Normality tests (Shapiro_Wilk test) of exponentiated Ct or dCt values of analyzed targets across compared groups and returned p values in columns with "_norm_p" pattern.
#' * Fold Change values (return in "FCh" column) together with log10 Fold change values (return in "log10FCh" column).
#'   Fold change values were calculated for each target by dividing  mean of exponentiated Ct od dCt values in study group by mean of exponentiated Ct or dCt values in reference group.
#' * Statistical testing of differences in exponentiated Ct or dCt values between study group and reference group.
#'   Student's t test and Mann-Whitney U test are implemented and resulted statistics (in column with "_test_stat" pattern) and p values (in column with "_test_p" pattern) are returned.
#'
#' @param data data object returned from exp_Ct() or exp_delta_Ct() function.
#' @param group.study character: name of study group (group of interest).
#' @param group.ref character: name of reference group.
#' @param do.tests logical: if TRUE, statistical significance of differences in exponentiated Ct values between compared groups will be calculated using Student's t test and Mann-Whitney U test. Default to TRUE.
#' @param save.to.txt logical: if TRUE, returned table with results will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension. Default to "RQ_expCt_results".
#
#' @return Data.frame with transformed Ct values and printed information about number and percentage of missing values.
#' @export
#'
#' @examples
#' library(tidyverse)
#' library(coin)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.Ct.exp <- exp_Ct(data.CtF, imput.by.mean.within.groups = TRUE)
#' RQ.Ct.exp <- RQ_exp_Ct(data.Ct.exp, group.study = "Disease", group.ref = "Control", do.tests = TRUE)
#' head(RQ.Ct.exp)
#'
#' @importFrom base as.data.frame as.factor mean
#' @importFrom stats sd shapiro.test t.test
#' @importFrom coin wilcox_test pvalue statistic
#' @importFrom utils write.table
#' @importFrom dplyr filter select rename_with full_join
#' @import tidyverse
#'
RQ_exp_Ct <- function(data,
                     group.study,
                     group.ref,
                     do.tests = TRUE,
                     save.to.txt = FALSE,
                     name.txt = "RQ_expCt_results"){

  data_slim <- data %>%
    filter(Group == group.study | Group == group.ref) %>%
    pivot_longer(cols = -c(Group, Sample), names_to = "Target", values_to = "value")

  data_slim$Group <- as.factor(data_slim$Group)

  data_FCh <- data_slim %>%
    group_by(Group, Target) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = expCt) %>%
    mutate(FCh = .data[[group.study]]/.data[[group.ref]]) %>%
    mutate(log10FCh = log10(FCh)) %>%
    as.data.frame()

  data_FCh_sd <- data_slim %>%
    group_by(Group, Target) %>%
    summarise(value_sd = sd(value, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = value_sd) %>%
    rename_with(~paste0(.x, "_sd", recycle0 = TRUE), all_of(c(group.study, group.ref)))

  if (do.tests == TRUE){

    data_FCh_norm <- data_slim %>%
      group_by(Group, Target) %>%
      summarise(shap_wilka_p = shapiro.test(value)$p.value, .groups = "keep") %>%
      pivot_wider(names_from = Group, values_from = shap_wilka_p) %>%
      rename_with(~paste0(.x, "_norm_p", recycle0 = TRUE), all_of(c(group.study, group.ref))) %>%
      full_join(data_FCh, by = c("Target")) %>%
      rename_with(~paste0(.x, "_mean", recycle0 = TRUE), all_of(c(group.study, group.ref)))

    data_FCh_tests <- data_slim %>%
      group_by(Target) %>%
      summarise(t_test_p = t.test(value ~ Group, alternative = "two.sided")$p.value,
                t_test_stat = t.test(value ~ Group, alternative = "two.sided")$statistic,
                MW_test_p = pvalue(wilcox_test(value ~ Group)),
                MW_test_stat = statistic(wilcox_test(value ~ Group)), .groups = "keep")
    data_FCh_norm_tests <- full_join(data_FCh_norm, data_FCh_tests, by = c("Target"))
    data_FCh_results <- full_join(data_FCh_norm_tests, data_FCh_sd, by = c("Target"))
    data_FCh_results <- select(data_FCh_results, Target, ends_with("_mean"), ends_with("_sd"), everything())

    return(data_FCh_results)

  } else {

    data_FCh_results <- full_join(data_FCh, data_FCh_sd, by = c("Target"))
    data_FCh_results <- rename_with(data_FCh_results, ~paste0(.x, "_mean", recycle0 = TRUE), all_of(c(group.study, group.ref)))
    data_FCh_results <- select(data_FCh_results, Target, ends_with("_mean"), ends_with("_sd"), everything())

    return(data_FCh_results)
  }
  if (save.to.txt == TRUE){
    write.table(data_FCh_results, paste(name.txt,".txt", sep = ""))
  }
}




#' @title select_ref_gene
#'
#' @description
#' This function draw line plot and calculate statistics (minimum, maximum, standard deviation,
#' variance and colinearity coefficient VIF) that could be helpful to select the best
#' reference gene for normalization of Ct values.
#'
#' @param data object returned from read_Ct_long(), read_Ct_wide() or filter_Ct() functions,
#' or data frame containing column named "Sample" with sample names, column named "Target" with target names,
#' column named "Ct" with raw Ct values, column named "Group" with group names.
#'     Any other columns could exist, but will not be used by this function.
#'
#' @param candidates character: vector of names of targets - candidates for gene reference.
#' @param colors character: vector of colors for targets, number of elements should be equal to number of candidate genes (elements in `candidates` vector).
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
#' @return List containing plot object and table with calculated statistics.
#'     Additional information about returned table is also printed, it could help user to properly interpret returned table.
#'     Created plot is displayed on graphic device.
#' @export
#'
#' @examples
#' library(car)
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' ref <-  select_ref_gene(data.CtF,
#'                         candidates = c("Gene4", "Gene8","Gene10","Gene16","Gene17", "Gene18"),
#'                         col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "#1F77B4", "black"))
#' ref[[2]]
#'
#' @importFrom base mean print min max as.data.frame
#' @importFrom stats sd var
#' @importFrom dplyr filter
#' @importFrom car vif
#' @import ggplot2
#' @import tidyverse
#'
select_ref_gene <- function(data,
                            candidates,
                            colors,
                            line.width = 1,
                            angle = 0,
                            x.axis.title = "",
                            y.axis.title = "Ct",
                            axis.title.size = 12,
                            axis.text.size = 10,
                            legend.title = "",
                            legend.title.size = 12,
                            legend.text.size = 12,
                            legend.position = "top",
                            plot.title = "",
				          	        plot.title.size = 14,
				          	        save.to.tiff = FALSE,
                            dpi = 600, width = 15, height = 15,
                            name.tiff = "Ct_reference_gene_selection"){

  ref <- data %>%
    group_by(Group, Target, Sample) %>%
    summarise(mean = mean(Ct, na.rm = TRUE),  .groups = "keep") %>%
    filter(Target %in% candidates)

  ref_plot <- ggplot(ref, aes(x = Sample, y = mean, color = Target, group = Target)) +
    geom_line(linewidth = line.width) +
    scale_color_manual(values = c(colors)) +
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

  if (angle != 0){
    ref_plot <- ref_plot +
      guides(x =  guide_axis(angle = angle))
  }

  print(ref_plot)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff, ".tiff", sep = ""), ref_plot, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }

  ref_var <- ref %>%
    group_by(Target) %>%
    summarise(min = min(mean),
              max = max(mean),
              sd = sd(mean, na.rm = TRUE),
              var = var(mean, na.rm = TRUE), .groups = "keep") %>%
    as.data.frame()

  ref_vif <- data %>%
    group_by(Group, Target, Sample) %>%
    summarise(mean = mean(Ct, na.rm = TRUE), .groups = "keep") %>%
    ungroup()

  ref_lm <- pivot_wider(ref_vif, names_from = Target, values_from = mean)

  ref_lm$dum <- c(1:nrow(ref_lm))
  model <- lm(dum ~ ., data = select(ref_lm, -Sample, -Group))
  vif <- vif(model)
  vif_sel <- vif[names(vif) %in% candidates]
  ref_var$VIF <- vif_sel
  ref_var$VIF.Target <- names(vif_sel)

  return(list(ref_plot, ref_var))
}





#' @title delta_Ct
#'
#' @description
#' This function collapses technical replicates (if present in data) by means,
#' counts and imputes missing data by means within groups (if so indicated),
#' and finally calculates delta Ct (dCt) values by subtracting Ct values of reference gene from Ct values of other genes.
#'
#' @param data data object returned from read_Ct_long(), read_Ct_wide() or filter_Ct() function,
#' or data frame containing column named "Sample" with sample names, column named "Target" with target names,
#' column named "Ct" with raw Ct values, column named "Group" with group names.
#'     Any other columns could exist, but will not be used by this function.
#' @param imput.by.mean.within.groups logical: if TRUE, missing values will be imputed by means within groups. Default to TRUE.
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
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' head(data.dCt)
#'
#' @importFrom base as.data.frame as.vector sum is.na ncol nrow paste cat
#' @importFrom utils write.table
#' @importFrom dplyr mutate filter select
#' @import tidyverse
#'
delta_Ct <- function(data,
                    imput.by.mean.within.groups = TRUE,
                    ref,
					          save.to.txt = FALSE,
                    name.txt = "data_dCt"){
  data <- data %>%
    group_by(Group, Target, Sample) %>%
    summarise(mean = base::mean(Ct, na.rm = TRUE)) %>%
    as.data.frame()
  data_wide <- data %>%
    select(Group, Sample, Target, mean) %>%
    pivot_wider(names_from = Target, values_from = mean)

  nas <- sum(is.na(data_wide))
  percentage <- sum(is.na(data_wide))/((ncol(data_wide)-2)*nrow(data_wide))

  if (imput.by.mean.within.groups == TRUE){

  data_wide_imp <- data_wide %>%
      group_by(Group) %>%
    mutate(across(where(is.numeric), ~ replace(., is.na(.), mean(., na.rm = TRUE))))

  dCt <- mutate_at(data_wide_imp,
                   vars(-c("Group", "Sample", ref)),
                   list(dCt = ~ . - .data[[ref]]))
  dCt <- select(dCt, Group, Sample, ends_with("dCt"))
  colnames(dCt) <- sub("_dCt*", "", colnames(dCt))
  cat("Data contained", nas, "missing values that constitute", round(percentage*100, 5), "percent of the total data.\n Missing values were imputed using means within compared groups.")
  if (save.to.txt == TRUE){
    write.table(as.data.frame(dCt), paste(name.txt,".txt", sep = ""))
  }
  return(dCt)

  } else {

    dCt <- mutate_at(data_wide,
                     vars(-c("Group", "Sample", ref)),
                     list(dCt = ~ . - .data[[ref]]))
    dCt <- select(dCt, Group, Sample, ends_with("dCt"))
    colnames(dCt) <- sub("_dCt*", "", colnames(dCt))
    cat("Data contains", nas, "missing values that constitute", round(percentage*100, 5), "percent of the total data.")
    if (save.to.txt == TRUE){
      write.table(as.data.frame(dCt), paste(name.txt,".txt", sep = ""))
    }
    return(dCt)
  }
}




#' @title control_boxplot_sample
#'
#' @description
#' Boxplot illustrating distribution of data in each sample. Could be useful to identify outlier samples.
#'
#' @param data object returned from delta_Ct(), exp_Ct() or exp_delta_Ct() function.
#' @param coef numeric: how many times of interquartile range should be used to indicate the most extend data point for whiskers. Default to 1.5.
#' @param colors character vector length of two, containing colors for compared groups. Default to c("#66c2a5", "#fc8d62").
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
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' control.boxplot.sample <- control_boxplot_sample(data.dCt)
#'
#' @importFrom base print
#' @import ggplot2
#' @import tidyverse
#'
control_boxplot_sample <- function(data, coef = 1.5,
                                   colors = c("#66c2a5", "#fc8d62"),
							                     x.axis.title = "Sample",
							                     y.axis.title = "value",
                                   axis.title.size = 12,
                                   axis.text.size = 12,
                                   legend.title = "Group",
                                   legend.title.size = 12,
                                   legend.text.size = 12,
                                   legend.position = "right",
                                   plot.title = "",
				                  			   plot.title.size = 14,
							                     save.to.tiff = FALSE,
                                   dpi = 600, width = 15, height = 15,
                                   name.tiff = "control_boxplot_samples"){

    data <- pivot_longer(data, !c(Sample, Group), names_to = "Target" , values_to = "value")

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




#' @title control_boxplot_target
#'
#' @description
#' This function creates boxplot illustrating distribution of data in each target. Could be useful to compare expression of analyzed targets.
#'
#' @param data object returned from delta_Ct(), exp_Ct() or exp_delta_Ct() function.
#' @param by.group logical: if TRUE, distributions will be drawn by compared groups of samples.
#' @param coef numeric: how many times of interquartile range should be used to indicate the most extend data point for whiskers.
#' @param colors character vector length of one (when by.group = FALSE) or two (when by.group = TRUE), containing colors for groups.
#' @param coef numeric: how many times of interquartile range should be used to indicate the most extend data point for whiskers. Default to 1.5.
#' @param x.axis.title character: title of x axis. Default to "Target".
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
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_boxplot_targets".
#'
#' @return Object with boxplot illustrating distribution of data for each target. Created plot will be also displayed on graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' control.boxplot.target <- control_boxplot_target(data.dCt)
#'
#' @importFrom base print
#' @import ggplot2
#' @import tidyverse
#'
control_boxplot_target <- function(data, coef = 1.5,
                                   by.group = TRUE,
                                   colors = c("#66c2a5", "#fc8d62"),
                                   axis.title.size = 12,
                                   axis.text.size = 12,
                                   x.axis.title = "Target",
                                   y.axis.title = "value",
                                   legend.title = "Group",
                                   legend.title.size = 12,
                                   legend.text.size = 12,
                                   legend.position = "right",
                                   plot.title = "",
                                   plot.title.size = 14,
                                   save.to.tiff = FALSE,
                                   dpi = 600, width = 15, height = 15,
                                   name.tiff = "control_boxplot_targets"){

  data <- pivot_longer(data, !c(Sample, Group), names_to = "Target" , values_to = "value")

  if (by.group == TRUE){
    box_control_targets <- ggplot(data, aes(x = Target, y = value, color = Group)) +
      geom_boxplot(coef = coef) +
      scale_x_discrete(limits = rev(unique(data$Target))) +
      scale_color_manual(values = c(colors)) +
      coord_flip() +
      xlab(x.axis.title) +
      ylab(y.axis.title) +
      labs(color = legend.title, title = plot.title) +
      theme_classic() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
      theme(plot.title = element_text(size = plot.title.size))

  } else {

    box_control_targets <- ggplot(data, aes(x = Target, y = value)) +
      geom_boxplot(coef = coef, fill = colors[1]) +
      scale_x_discrete(limits = rev(unique(data$Target))) +
      coord_flip() +
      xlab(x.axis.title) +
      ylab(y.axis.title) +
      labs(title = plot.title) +
      theme_classic() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
      theme(plot.title = element_text(size = plot.title.size))
  }

  print(box_control_targets)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), box_control_targets, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
  return(box_control_targets)
}





#' @title control_cluster_sample
#'
#' @description
#' Performs hierarchical clustering of samples based on the data. Could be useful to identify outlier samples.
#'
#' @param data object returned from delta_Ct(), exp_Ct() or exp_delta_Ct() function.
#' @param method.dist character: name of method used for calculation of distances, derived from stats::dist() function, should be one of "euclidean" (default) , "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param method.clust character: name of used method for agglomeration, derived from stats::hclust() function, should be one of "ward.D", "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".
#' @param x.axis.title character: title of x axis. Default to "Samples".
#' @param y.axis.title character: title of y axis. Default to "Height".
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
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' control_cluster_sample(data.dCt)
#'
#' @importFrom stats hclust dist
#' @importFrom base plot
#' @importFrom dplyr select
#' @import tidyverse
#'
control_cluster_sample <- function(data,
                                   method.dist = "euclidean",
                                   method.clust = "average",
                                   x.axis.title = "Samples",
                                   y.axis.title = "Height",
                                   plot.title = "",
                                   save.to.tiff = FALSE,
                                   dpi = 600, width = 15, height = 15,
                                   name.tiff = "control_clust_samples"){
  data <- ungroup(data)
  cluster <- hclust(dist(select(data, -Group, -Sample), method = method.dist), method = method.clust)
  cluster$labels <- data$Sample
  plot(cluster, xlab = x.axis.title, ylab = y.axis.title, main = plot.title)

  if (save.to.tiff == TRUE){
    tiff(paste(name.tiff, ".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
    plot(cluster, xlab  = x.axis.title, ylab = y.axis.title, main = plot.title)
    dev.off()
  }
}




#' @title control_cluster_target
#'
#' @description
#' Performs hierarchical clustering of targets based on the data. Could be useful to gain insight into similarity in expression of analyzed targets.
#'
#' @param data object returned from delta_Ct(), exp_Ct() or exp_delta_Ct() function.
#' @param method.dist character: name of method used for calculation of distances, derived from stats::dist() function, should be one of "euclidean" (default) , "maximum", "manhattan", "canberra", "binary" or "minkowski".
#' @param method.clust character: name of used method for agglomeration, derived from stats::hclust() function, should be one of "ward.D", "ward.D2", "single", "complete", "average" (default), "mcquitty", "median" or "centroid".
#' @param x.axis.title character: title of x axis. Default to "Targets".
#' @param y.axis.title character: title of y axis. Default to "Height".
#' @param plot.title character: title of plot. Default to "".
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_clust_targets".
#'
#' @return Plot with hierarchical clustering of targets, displayed on graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' control_cluster_target(data.dCt)
#'
#' @importFrom stats hclust dist
#' @importFrom base plot t colnames
#' @importFrom dplyr select
#' @import tidyverse
#'
control_cluster_target <- function (data,
                                    method.dist = "euclidean",
                                    method.clust = "average",
                                    x.axis.title = "Targets",
                                    y.axis.title = "Height",
                                    plot.title = "",
                                    save.to.tiff = FALSE,
                                    dpi = 600, width = 15, height = 15,
                                    name.tiff = "control_clust_targets")
{
  data_t <- ungroup(data)
  data_t <- t(select(data_t, -Group, -Sample))
  colnames(data_t) <- data$Sample
  cluster <- hclust(dist(as.data.frame(data_t), method = method.dist), method = method.clust)
  plot(cluster, xlab = x.axis.title, ylab = y.axis.title, main = plot.title)

  if (save.to.tiff == TRUE) {
    tiff(paste(name.tiff, ".tiff", sep = ""), res = dpi,
         width = width, height = height, units = "cm", compression = "lzw")
    plot(cluster, xlab = x.axis.title, ylab = y.axis.title,
         main = plot.title)
    dev.off()
  }
}




#' @title control_pca_sample
#'
#' @description
#' Performs principal component analysis (PCA) for samples and generate plot illustrating spatial arrangement of samples using two first components. Could be useful to identify outlier samples.
#'     IMPORTANT: PCA analysis can not deal with missing values, thus all samples with at least one missing value are removed from data before analysis. It is recommended to run this function on data after imputation of missing values.
#'
#' @param data object returned from delta_Ct(), exp_Ct() or exp_delta_Ct() function.
#' @param point.size numeric: size of points. Default to 4.
#' @param point.shape integer: shape of points. Default to 19.
#' @param alpha numeric: transparency of points, a value between 0 and 1. Default to 0.7.
#' @param label.size numeric: size of points labels (names of samples). Default to 3.
#' @param hjust numeric: horizontal position of points labels. Default to 0.
#' @param vjust numeric: vertical position of points labels.  Default to -1.
#' @param colors character vector length of two, containing colors for compared groups.
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
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.prepared <- Ct_for_control(data.CtF)
#' control_pca_sample(data.CtF.prepared)
#'
#' @importFrom base print as.data.frame rownames summary paste round
#' @importFrom stats na.omit prcomp
#' @import ggplot2
#' @import tidyverse
#'
control_pca_sample <- function(data,
                               point.size = 4,
                               point.shape = 19,
                               alpha = 0.7,
                               colors = c("#66c2a5", "#fc8d62"),
                               label.size = 3,
                               hjust = 0,
                               vjust = -1,
                               axis.title.size = 12,
                               axis.text.size = 10,
                               legend.text.size = 12,
                               legend.title = "Group",
                               legend.title.size = 12,
                               legend.position = "right",
                               plot.title = "",
                               plot.title.size = 14,
                               save.to.tiff = FALSE,
                               dpi = 600, width = 15, height = 15,
                               name.tiff = "control_pca_samples"){

  data <- as.data.frame(data)
  rownames(data) <- data$Sample
  data <- na.omit(data)
  pca <- prcomp(select(data, -Sample, -Group))
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



#' @title control_pca_target
#'
#' @description
#' Performs principal component analysis (PCA) for targets and generate plot illustrating spatial arrangement of targets using 2 components. Could be useful to gain insight into similarity in expression of analyzed targets.
#' IMPORTANT: PCA analysis can not deal with missing values, thus all targets with at least one missing value are removed from data before analysis. It is recommended to run this function on data after imputation of missing values.
#'
#' @param data object returned from `delta_Ct()`, `exp_Ct()` or `exp_delta_Ct()` function.
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
#' See description for `legend.position` parameter in [ggplot2::theme()] function.
#' @param plot.title character: title of plot. Default to "".
#' @param plot.title.size integer: font size of plot title. Default to 14.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file. Default to FALSE.
#' @param dpi integer: resolution of saved .tiff file. Default to 600.
#' @param width numeric: width (in cm) of saved .tiff file. Default to 15.
#' @param height integer: height (in cm) of saved .tiff file. Default to 15.
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_pca_targets".
#'
#' @return Object with plot illustrating spatial arrangement of targets according to coordinates of 2 components obtained from principal component analysis (PCA). The plot will be also displayed in graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.CtF.prepared <- Ct_for_control(data.CtF)
#' control_pca_target(data.CtF.prepared)
#'
#' @importFrom base print as.data.frame rownames summary paste round t colnames
#' @importFrom stats na.omit prcomp
#' @import ggplot2
#' @import tidyverse
#'
control_pca_target <- function(data, point.size = 4, point.shape = 19, alpha = 0.7,
                               label.size = 3, hjust = 0, vjust = -1,
                               color = "black",
                               axis.title.size = 12,
                               axis.text.size = 10,
                               legend.text.size = 12,
                               legend.title = "Group",
                               legend.title.size = 12,
                               legend.position = "right",
                               plot.title = "",
                               plot.title.size = 14,
                               save.to.tiff = FALSE,
                               dpi = 600, width = 15, height = 15,
                              name.tiff = "control_pca_targets"){

  data <- ungroup(data)
  data <- as.data.frame(data)
  data_t <- t(select(data, -Group, -Sample))
  colnames(data_t) <- data$Sample
  data_t <- na.omit(data_t)
  pca <- prcomp(data_t)
  var_pca1 <- summary(pca)$importance[2,][1]
  var_pca2 <- summary(pca)$importance[2,][2]
  pca_comp <- as.data.frame(pca$x)
  pca_comp$Target <- rownames(pca_comp)

  control_pca <- ggplot(pca_comp, aes(x = PC1, y = PC2, label = Target)) +
    geom_point(size = point.size, shape = point.shape, alpha = alpha, col = color) +
    labs(title = plot.title) +
    theme_bw() +
    labs(x = paste("PC1: ", round(var_pca1*100,2), "% variance explained", sep = ""),
         y = paste("PC2: ", round(var_pca2*100,2), "% variance explained", sep = "")) +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour="black")) +
    theme(plot.title = element_text(size = plot.title.size)) +
    geom_text(aes(label = Target), hjust = hjust, vjust = vjust, size = label.size)

  print(control_pca)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), control_pca, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
  return(control_pca)
}





#' @title control_corr_target
#'
#' @description
#' Performs correlation analysis of targets based on the data. Could be useful to gain insight into relationships between analyzed targets.
#'
#' @param data object returned from delta_Ct(), exp_Ct() or exp_delta_Ct() function.
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
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_corr_targets".
#' @param save.to.txt logical: if TRUE, correlation results (sorted by descending absolute values of correlation coefficients) will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension.. Default to "control_corr_targets".
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
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' corr.targets <- control_corr_target(data.dCt)
#' head(corr.targets)
#'
#' @importFrom base as.data.frame paste upper.tri rownames
#' @importFrom stats p.adjust
#' @importFrom utils write.table
#' @importFrom Hmisc rcorr
#' @import corrplot
#' @import tidyverse
#'
control_corr_target <- function(data,
                                method = "pearson",
                                add.coef = "black",
                                order = "hclust",
                                hclust.method = "average",
                                size = 0.6,
                                coef.size = 0.6,
                                p.adjust.method = "BH",
                                save.to.tiff = FALSE,
                                dpi = 600, width = 15, height = 15,
                                name.tiff = "control_corr_targets",
                                save.to.txt = FALSE,
                                name.txt = "control_corr_targets"){

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
  return(corr.data.sort)
}





#' @title control_corr_sample
#'
#' @description
#' Performs correlation analysis of samples based on the data. Could be useful to gain insight into relationships between analyzed samples.
#'
#' @param data object returned from delta_Ct(), exp_Ct() or exp_delta_Ct() function.
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
#' @param name.tiff character: name of saved .tiff file, without ".tiff" name of extension. Default to "control_corr_samples".
#' @param save.to.txt logical: if TRUE, correlation results (sorted by descending absolute values of correlation coefficients) will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension.. Default to "control_corr_samples".
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
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' corr.samples <- control_corr_sample(data.dCt)
#' head(corr.samples)
#'
#' @importFrom base as.data.frame paste upper.tri rownames
#' @importFrom stats p.adjust
#' @importFrom utils write.table
#' @importFrom Hmisc rcorr
#' @import corrplot
#' @import tidyverse
#'
control_corr_sample <- function(data,
                                method = "pearson",
                                add.coef = "black",
                                order = "hclust",
                                hclust.method = "average",
                                size = 0.6,
                                coef.size = 0.6,
                                p.adjust.method = "BH",
                                save.to.tiff = FALSE,
                                dpi = 600, width = 15, height = 15,
                                name.tiff = "control_corr_samples",
                                save.to.txt = FALSE,
                                name.txt = "control_corr_samples"){

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

  return(corr.data.sort)
}







#' @title single_pair_target
#'
#' @description
#' Generate scatter plot with linear regression line for two specified targets. Could be useful to assess linear relationship between these targets.
#'
#' @param data object returned from delta_Ct(), exp_Ct() or exp_delta_Ct() function.
#' @param x,y characters: names of targets to use.
#' @param by.group logical: if TRUE (default), relationships will be shown separately for compared groups.
#' @param point.size numeric: size of points. Default to 4.
#' @param point.shape integer: shape of points. Default to 19.
#' @param point.alpha numeric: transparency of points, a value between 0 and 1. Default to 0.7.
#' @param colors character: color for points (if by.group = FALSE) or vector of two colors for points, used to distinguish groups (if by.group = TRUE).
#' @param line.alpha numeric: transparency of regression line, a value between 0 and 1. Default to 0.7.
#' @param axis.title.size integer: font size of axis titles. Default to 12.
#' @param axis.text.size integer: font size of axis text. Default to 10.
#' @param legend.title character: title of legend. Default to "Group".
#' @param legend.title.size integer: font size of legend title. Default to 12.
#' @param legend.text.size integer: font size of legend text. Default to 12.
#' @param legend.position position of the legend, one of "top", "right" (default), "bottom", "left", or "none" (no legend).
#' See description for `legend.position` parameter in ggplot2::theme() function.
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
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' single_pair_target(data.dCt, "Gene16, "Gene17)
#'
#' @importFrom base print paste
#' @import ggplot2
#' @import ggpmisc
#'
single_pair_target <- function(data, x, y, by.group = TRUE,
                            point.size = 4, point.shape = 19, point.alpha = 0.7,
                            colors = c("#66c2a5", "#fc8d62"),
                            line.alpha = 0.7,
                            axis.title.size = 12,
                            axis.text.size = 10,
                            legend.title = "Group",
                            legend.title.size = 12,
                            legend.text.size = 12,
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
                            name.tiff = "targetss_single_plot"){

  if (by.group == TRUE){
    single_pair <- ggplot(data, aes(x = .data[[x]], y = .data[[y]], color = Group)) +
      geom_point(size = point.size, shape = point.shape, alpha = point.alpha) +
      geom_smooth(method='lm', se = FALSE, alpha = line.alpha) +
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
      single_pair <- corr_control +
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
      geom_smooth(method='lm', se = FALSE, alpha = line.alpha) +
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
#' @param data object returned from delta_Ct(), exp_Ct() or exp_delta_Ct() function.
#' @param x,y characters: names of targets to use.
#' @param point.size numeric: size of points.
#' @param point.shape integer: shape of points.
#' @param point.alpha numeric: transparency of points, a value between 0 and 1.
#' @param color character: color used for points.
#' @param line.alpha numeric: transparency of regression line, a value between 0 and 1.
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
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
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
                                   line.alpha = 0.7,
                                   axis.title.size = 12,
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
                                   name.tiff = "samples_single_plot"){

  data <- as.data.frame(data)
  data_t <- t(select(data, -Group, -Sample))
  colnames(data_t) <- data$Sample

  single_pair_t <- ggplot(as.data.frame(data_t), aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(size = point.size, shape = point.shape, alpha = point.alpha, color = color) +
    geom_smooth(method='lm', se = FALSE, alpha = alpha) +
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





#' @title filter_transformed_Ct
#'
#' @description
#' Filters transformed Ct data (2^(-Ct), delta Ct, and 2^(-dCt) data) according to the used filtering criteria (see parameters).
#'
#' @param data object returned from delta_Ct(), exp_Ct() or exp_delta_Ct() functions.
#' @param remove.Target character: vector with names of targets which should be removed from data.
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
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' data.dCt.exp <- exp_delta_Ct(data.dCt)
#' data.dCt.expF <- filter_transformed_Ct(data.dCt.exp, remove.Sample = c("Control11"))0
#'
#' dim(data.dCt.exp)
#' dim(data.dCt.expF)
#'
#'
#' @importFrom dplyr filter select
#' @import tidyverse
#'
filter_transformed_Ct <- function(data,
                                  remove.Target = c(""),
                                  remove.Sample = c(""),
                                  remove.Group = c("")){

  data <- filter(data,
                 !Sample %in% remove.Sample,
                 !Group %in% remove.Group)

  data <- select(data, -any_of(remove.Target))

  return(data)
}





#' @title exp_delta_Ct
#'
#' @description
#' This function exponentiates delta Ct (dCt) values by using formula 2^(-dCt).
#'
#' @param data data object returned from delta_Ct() function.
#' @param save.to.txt logical: if TRUE, returned data will be saved to .txt file. Default to FALSE.
#' @param name.txt character: name of saved .txt file, without ".txt" name of extension. Default to "data.exp.dCt".
#'
#' @return Data.frame with exponentiated dCt values.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' data.dCt.exp <- exp_delta_Ct(data.dCt)
#' head(data.dCt.exp)
#'
#' @importFrom base as.data.frame paste
#' @importFrom utils write.table
#' @importFrom dplyr mutate_at
#' @import tidyverse
#'
exp_dCt <- function(data,
                    save.to.txt = FALSE,
                    name.txt = "data.exp.dCt"){

  calcRQ <- function(x){
    x <- 2^-x
  }

  data_exp <- data %>%
    mutate_at(vars(-c("Group","Sample")), calcRQ)

  if (save.to.txt == TRUE){
    write.table(as.data.frame(data_exp), paste(name.txt,".txt", sep = ""))
  }
  return(data_exp)
}




#' @title results_boxplot
#'
#' @description
#' This function creates boxplot illustrating distribution of data in selected targets.
#'     It is similar to control_boxplot_target() function; however, some new options are added,
#'     including target selection, faceting, and adding mean labels to boxes.
#'     This, this function could be useful to present results for finally selected targets.
#'
#' @param data object returned from delta_Ct(), exp_Ct() or exp_delta_Ct() function.
#' @param coef numeric: how many times of interquartile range should be used to indicate the most extend data point for whiskers. Default to 1.5.
#' @param sel.Target character vector with names of targets to include, or "all" (default) to use all names of targets.
#' @param by.group logical: if TRUE (default), distributions will be drawn by compared groups of samples.
#' @param faceting logical: if TRUE (default), plot will be drawn with facets with free scales using ggplot2::facet_wrap() function (see its documentation for more details).
#' @param facet.row,facet.col integer: number of rows and column to arrange facets.
#' @param angle integer: value of angle in which names of genes should be displayed. Default to 0.
#' @param rotate logical: if TRUE, boxplots will be arranged horizontally. Deafault to FALSE.
#' @param add.mean logical: if TRUE, means will be added to boxes as squares. Default to TRUE.
#' @param add.mean.size numeric: size of squares indicating means. Default to 2.
#' @param add.mean.color character: color of squares indicating means. Default to "black".
#' @param colors character vector length of one (when by.group = FALSE) or two (when by.group = TRUE), containing colors for groups.
#' @param x.axis.title character: title of x axis. Default to "Target".
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
#' @return Object with boxplot illustrating distribution of data for selected targets. Created plot will be also displayed on graphic device.
#' @export
#'
#' @examples
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' data.dCt.exp <- exp_delta_Ct(data.dCt)
#' data.dCt.expF <- filter_transformed_Ct(data.dCt.exp, remove.Sample = c("Control11"))
#' results_boxplott(data.dCt,
#'                  sel.Target = c("Gene1","Gene12","Gene16","Gene19","Gene20")
#'                  facet.row = 3,
#'                  facet.col = 5)
#'
#' @importFrom base print paste
#' @importFrom dplyr filter
#' @import ggplot2
#' @import tidyverse
#'
results_boxplot <- function(data,
                            coef = 1.5,
                            sel.Target = "all",
                            by.group = TRUE,
                            faceting = TRUE,
                            facet.row,
                            facet.col,
                            angle = 0,
                            rotate = FALSE,
                            add.mean = TRUE,
                            add.mean.size = 2,
                            add.mean.color = "black",
                            colors = c("#66c2a5", "#fc8d62"),
                            x.axis.title = "Sample",
                            y.axis.title = "value",
                            axis.title.size = 12,
                            axis.text.size = 10,
                            legend.text.size = 12,
                            legend.title = "Group",
                            legend.title.size = 12,
                            legend.position = "top",
                            plot.title = "",
                            plot.title.size = 14,
                            save.to.tiff = FALSE,
                            dpi = 600, width = 15, height = 15,
                            name.tiff = "results_boxplot"){

  data <- pivot_longer(data, !c(Sample, Group), names_to = "Target" , values_to = "value")

  if (sel.Target[1] == "all"){
    data <- data

  } else {

    data <- filter(data, Target %in% sel.Target)
  }

  if (by.group == TRUE){

    box_results <- ggplot(data, aes(x = Target, y = value, fill = Group)) +
      geom_boxplot(coef = coef) +
      scale_fill_manual(values = c(colors)) +
      xlab(x.axis.title) +
      ylab(y.axis.title) +
      labs(fill = legend.title, title = plot.title) +
      theme_bw() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(panel.grid.major.x = element_blank())

    if (faceting == TRUE) {
      box_results <- box_results +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(vars(Target), scales = "free", nrow = facet.row, ncol = facet.col)
    }

  } else {

    box_results <- ggplot(data, aes(x = Target, y = value)) +
      geom_boxplot(coef = coef, fill = colors[1]) +
      xlab(x.axis.title) + ylab(y.axis.title) +
      labs(fill = legend.title, title = plot.title) +
      theme_bw() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(panel.grid.major.x = element_blank())

    if (faceting == TRUE) {
      box_results <- box_results +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
        facet_wrap(vars(Target), scales = "free", nrow = facet.row, ncol = facet.col)
    }

  }

  if (angle != 0){
    box_results <- box_results +
      guides(x =  guide_axis(angle = angle))
  }

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





#' @title RQ_exp_ddCt
#'
#' @description
#' Performs relative quantification of gene expression using 2^(-ddCt) method.
#'
#' @details
#' This function calculates:
#' * Means (return in columns with "_mean" pattern) and standard deviations (return in columns with "_sd" pattern) of delta Ct values of analyzed targets across compared groups.
#' * Normality tests (Shapiro_Wilk test) of delta Ct values of analyzed targets across compared groups and returned p values in columns with "_norm_p" pattern.
#' * Differences in mean delta Ct values of targets between compared groups, obtaining delta delta Ct values (in returned "ddCt" column).
#' * RQ values (return in "RQ" column) together with log10 RQ values (return in "log10RQ" column).
#'   RQ values are calculated for each target by exponentiating ddCt values using 2^-ddCt formula.
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
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' RQ.ddCt.exp <- RQ_exp_ddCt(data, "Disease", "Control")
#' head(RQ.ddCt.exp)
#'
#' @importFrom base as.data.frame as.factor mean
#' @importFrom stats sd shapiro.test t.test
#' @importFrom coin wilcox_test pvalue statistic
#' @importFrom utils write.table
#' @importFrom dplyr filter select rename_with full_join
#' @import tidyverse
#'
RQ_exp_ddCt <- function(data,
                    group.study,
                    group.ref,
                    do.tests = TRUE,
                    save.to.txt = FALSE,
                    name.txt = "ddCt_RQ_results"){

  data_slim <- data %>%
    filter(Group == group.study | Group == group.ref) %>%
    pivot_longer(cols = -c(Group, Sample), names_to = "Target", values_to = "dCt")

  data_slim$Group <- as.factor(data_slim$Group)

  data_ddCt <- data_slim %>%
    group_by(Group, Target) %>%
    summarise(ddCt = mean(dCt, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = ddCt) %>%
    mutate(ddCt = .data[[group.study]] - .data[[group.ref]]) %>%
    mutate(RQ = 2^-ddCt) %>%
    rename_with(~paste0(.x, "_mean", recycle0 = TRUE), all_of(c(group.study, group.ref)))

  if (do.tests == TRUE){

    data_ddCt_norm <- data_slim %>%
      group_by(Group, Target) %>%
      summarise(shap_wilka_p = shapiro.test(dCt)$p.value, .groups = "keep") %>%
      pivot_wider(names_from = Group, values_from = shap_wilka_p) %>%
      full_join(data_ddCt, by = c("Target")) %>%
      rename_with(~paste0(.x, "_norm_p", recycle0 = TRUE), all_of(c(group.study, group.ref)))

    data_ddCt_tests <- data_slim %>%
      group_by(Target) %>%
      summarise(t.test.p = t.test(dCt ~ Group, alternative = "two.sided")$p.value,
                t.test.stat = t.test(dCt ~ Group, alternative = "two.sided")$statistic,
                MW.test.p = pvalue(wilcox_test(dCt ~ Group)),
                MW.test.stat = statistic(wilcox_test(dCt ~ Group)), .groups = "keep")

    data_ddCt_norm_tests <- full_join(data_ddCt_norm, data_ddCt_tests, by = c("Target"))

    return(data_ddCt_norm_tests)

    } else{

    return(data_ddCt)
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
#' @param data object returned from RQ_exp_ddCt() or RQ_exp_Ct() functions.
#' @param use.p logical: if TRUE, bars of statistically significant genes will be distinguished by colors.
#' @param mode character: which p value should be used? One of the "t" (p values from Student's t test),
#' "mw" (p values from Mann-Whitney U test), "depends" (if data in both compared groups were considered as derived from normal distribution (p value from Shapiro_Wilk test > 0.05) - p
#' values from Student's t test will be used for significance assessment, otherwise p values from Mann-Whitney U test will be used for significance assessment).
#' There is one more option, if user intend to use another p values, ex. obtained from other statistical test,
#' a mode parameter could be set to "user". In this situation, before run RQ_plot function, user should to prepare
#' data.frame object names "user"with two columns, one named "Target" with Target names and second with p values. The order of columns must be kept as described.
#' @param p.threshold numeric: threshold of p values for statistical significance.
#' @param use.log10FCh logical: if TRUE, the criterion of fold change will be also used for significance assessment of genes.
#' @param log10FCh.threshold numeric: threshold of log10 fold change values used for significance assessment of genes.
#' @param sel.Target character vector with names of targets to include, or "all" (default) to use all names of targets.
#' @param bar.width numeric: width of bars..
#' @param angle integer: value of angle in which names of genes should be displayed. Default to 0.
#' @param rotate logical: if TRUE, bars will be arranged horizontally. Deafault to FALSE.
#' @param colors character vector length of one (when use.p = FALSE) or two (when use.p = TRUE), containing colors for significant and no significant genes.
#' @param x.axis.title character: title of x axis. Default to "Target".
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
#' library(tidyverse)
#' data(data.Ct)
#' data.CtF <- filter_Ct(data.Ct,
#'                       remove.Target = c("Gene2","Gene5","Gene6","Gene9","Gene11"),
#'                       remove.Sample = c("Control08","Control16","Control22"))
#' data.dCt <- delta_Ct(data.CtF,
#'                      imput.by.mean.within.groups = TRUE,
#'                      ref = "Gene8")
#' data.dCt.exp <- exp_delta_Ct(data.dCt)
#' data.dCt.expF <- filter_transformed_Ct(data.dCt.exp, remove.Sample = c("Control11"))
#' RQ.ddCt.exp <- RQ_exp_ddCt(data, "Disease", "Control")
#' RQ.plot <- RQ_plot(RQ.ddCt.exp, mode = "depends", Target.sel = c("Gene1","Gene12","Gene16","Gene19","Gene20"))
#' head(RQ.plot[[2]])
#'
#' @importFrom base print paste colnames factor
#' @importFrom dplyr filter
#' @importFrom stats reorder
#' @import ggplot2
#' @import tidyverse
#'
RQ_plot <- function(data,
                    use.p = TRUE,
                    mode,
                    p.threshold = 0.05,
                    use.log10FCh = FALSE,
                    log10FCh.threshold = 0,
                    Target.sel = "all",
                    bar.width = 0.8,
                    angle = 0,
                    rotate = FALSE,
                    colors = c("#66c2a5", "#fc8d62"),
                    x.axis.title = "Target",
                    y.axis.title = "log10(Fold change)",
                    axis.title.size = 12,
                    axis.text.size = 10,
                    legend.text.size = 12,
                    legend.title = "Selected as significant?",
                    legend.title.size = 12,
                    legend.position = "top",
                    plot.title = "",
                    dpi = 600, width = 15, height = 15,
                    save.to.tiff = FALSE,
                    name.tiff = "RQ_plot"){

  if (Target.sel[1] != "all"){
    data <- filter(data, Target %in% Target.sel)
  }

  if(sum(colnames(data) %in% "RQ") > 0){
    data <- rename(data, FCh = RQ)
  }

  if (use.p == TRUE){
    if (mode == "t"){
      data$p.used <- data$t.test.p
    }
    if (mode == "mw"){
      data$p.used <- data$MW.test.p
    }
    if (mode == "depends"){
      vars <- colnames(select(data, ends_with("norm_p")))
      data <- mutate(data, test.for.comparison = ifelse(.data[[vars[[1]]]] >= 0.05 & .data[[vars[[2]]]] >= 0.05, "t.student's.test", "Mann-Whitney.test"))
      data <- mutate(data, p.used = ifelse(test.for.comparison == "t.student's.test", t.test.p, MW.test.p))
    }
    if (mode == "user"){
      colnames(user) <- c("Target","p.used")
      data <- full_join(data, user, by = c("Target"))
    }
    if (use.log10FCh == TRUE){
      data <- mutate(data, `Selected as significant?` = ifelse(p.used > p.threshold, yes = "No",
                                                                 no = ifelse(abs(log10(FCh)) <  log10FCh.threshold, "No", "Yes")))
      data$`Selected as significant?` <- factor(data$`Selected as significant?`, levels = c("Yes", "No"))

          } else {

    data <- mutate(data, `Selected as significant?` = ifelse(p.used > p.threshold, yes = "No (p > 0.05)",  no = "Yes (p <= 0.05)"))
    data$`Selected as significant?` <- factor(data$`Selected as significant?`, levels = c("Yes (p <= 0.05)", "No (p > 0.05)"))
          }

    RQ <- ggplot(data, aes(x = reorder(Target, -FCh), y = log10(FCh), fill = `Selected as significant?`)) +
      geom_col(width = bar.width) +
      scale_fill_manual(values = c(colors)) +
      xlab(x.axis.title) +
      ylab(y.axis.title) +
      labs(fill = legend.title, title = plot.title) +
      theme_bw() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(panel.grid.major.x = element_blank()) +
      geom_hline(yintercept = 0, linewidth = 0.4)

    } else {

      RQ <- ggplot(data, aes(x = reorder(Target, -FCh), y = log10(FCh))) +
      geom_col(width = bar.width, fill = colors[1]) +
      xlab(x.axis.title) + ylab(y.axis.title) +
      labs(title = plot.title) +
      theme_bw() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(panel.grid.major.x = element_blank())
  }

  if (angle != 0){
    RQ <- RQ +
      guides(x =  guide_axis(angle = angle))
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













ROC <- function(data,
                Target.sel = "all",
                Groups,
                save.to.tiff = FALSE, plot.title = "",
                dpi = 600, width = 15, height = 15,
                panels.row,
                panels.col,
                name.tiff = "dCt_ROC_plot",
                text.size = 1.1,
                print.auc = TRUE,
                print.auc.cex = 0.8,
                save.to.txt = FALSE,
                name.txt = "dCt_ROC_results"){
  data <- filter(data, Group %in% Groups)
  if (Target.sel[1] != "all"){
    data <- data[, colnames(data) %in% c("Group", "Sample", Target.sel)]
  }
  roc_param <- as.data.frame(matrix(nrow = ncol(data)-2, ncol = 10))
  colnames(roc_param) <- c("Target","Threshold", "Specificity", "Sensitivity", "Accuracy", "ppv", "npv", "youden", "AUC", "Target.control")
  roc_param$Target <- colnames(data)[-c(1:2)]
  for (x in 1:nrow(roc_param)){
    myproc <- roc(response = data$Group, predictor = as.data.frame(data)[ ,x+2], levels = c(Groups),
                  smooth = FALSE, auc = TRUE, plot=FALSE, ci=TRUE, of = "auc", quiet = TRUE)
    parameters <- coords(myproc, "best", ret = c("threshold", "specificity", "sensitivity","accuracy", "ppv", "npv", "youden"))
    roc_param[x,2:8] <- parameters
    roc_param[x,9] <- myproc$auc
    roc_param[x,10] <- colnames(data)[x+2]
    if (nrow(parameters) > 1){
      cat('Warning: ',colnames(data)[x+2],'has more than 1 threshold value for calculated Youden’s J statistic.\n')
    } else {}
  }
  if (save.to.tiff == TRUE){
    tiff(paste(name.tiff,".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
    par(mfrow = c(panels.row, panels.col))
    for (x in 1:nrow(roc_param)){
      myproc <- roc(response = data$Group, predictor = as.data.frame(data)[ ,x+2], levels = c(Groups),
                    smooth = FALSE, auc = TRUE, plot=FALSE, ci=TRUE, of = "auc", quiet = TRUE)
      plot.roc(myproc, main = roc_param$Target[x],
               smooth = FALSE, cex.axis = text.size, cex.lab = text.size, identity.lwd = 2,
               plot = TRUE, percent = TRUE, print.auc = print.auc, print.auc.x = 0.85, print.auc.y = 0.1, print.auc.cex = print.auc.cex)
    }
    dev.off()
  }
  if (save.to.txt == TRUE){
    write.table(roc_param, paste(name.txt,".txt", sep = ""))
  }
  return(roc_param)
}





log_reg <- function(data, remove.Inf.NA = FALSE,
                    Target.sel = "all",
                    group.study,
                    group.ref,
                    ci = 0.95,
                    axis.title.size = 12,
                    axis.text.size = 10,
                    legend.text.size = 12,
                    x.axis.title = "Odds ratio",
                    y.axis.title = "",
                    legend.title = "p value",
                    legend.title.size = 12,
                    legend.position = "right",
                    plot.title = "",
                    dpi = 600, width = 15, height = 15,
                    save.to.tiff = FALSE,
                    name.tiff = "dCt_OR_plot",
                    save.to.txt = FALSE,
                    name.txt = "dCt_OR_results",
                    log.axis = FALSE){
  data <- filter(data, Group %in% c(group.study, group.ref))
  if (Target.sel[1] != "all"){
    data <- data[, colnames(data) %in% c("Group", "Sample", Target.sel)]
  }
  data <- mutate(data, Group_num = ifelse(Group == group.study, 0, 1))
  n.targets <- ncol(data)-3
  list.models <- lapply(data[3:(n.targets+2)], function(x) glm(data$Group_num ~ x, data = data, family = binomial))
  list.CI <- lapply(names(list.models)[1:n.targets], function(x) or_glm(data = data,
                                                           model = list.models[[x]],
                                                           incr = list(x = 1),
                                                           ci = ci))
  data.CI <- as.data.frame(matrix(ncol = 8, nrow = n.targets))
  colnames(data.CI) <- c("Target", "oddsratio", "CI_low", "CI_high", "Intercept", "coeficient","p_intercept","p_coef")
  for (x in 1:n.targets){
    data.CI$Target <- names(list.models)
    data.CI[x,2:4] <- as.vector(list.CI)[[x]][2:4]
    data.CI[x,5:6] <- list.models[[x]]$coefficients
    data.CI[x,7:8] <- coef(summary(list.models[[x]]))[,4]
  }
  od_df <- data.frame(yAxis = 1:nrow(data.CI),
                        boxOdds = data.CI$oddsratio,
                        boxCILow = data.CI$CI_low,
                        boxCIHigh = data.CI$CI_high,
                        boxLabels = data.CI$Target,
                        p = data.CI$p_coef)

  if (remove.Inf.NA == TRUE){
    od_df[sapply(od_df, is.infinite)] <- NA
    od_df <- na.omit(od_df)
  }

  odd.ratio <- ggplot(od_df, aes(x = boxOdds, y = boxLabels, label = boxOdds)) +
    geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
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
    theme(legend.title = element_text(size = legend.title.size, colour="black"))
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
  return(data.CI)
}

