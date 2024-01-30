
#' @title read_Ct_wide
#'
#' @description
#' Imports two tables: first is a wide-type table with Ct values and design file. Merges both files to return long table ready for analysis.
#' All parameters of this function have not default values and schould be specified by user.
#'
#'
#' @param path.Ct.file path to wide-format table in .txt format containing Ct values, target names in the first column, and
#' sample names in the first row. In other words, this table should contain targets by rows and samples by columns.
#' @param path.design.file path to .txt file with two columns: column named "Sample" with names of samples
#' and column named "Group" with names of groups assigned to samples. Names of samples in this file
#' should correspond to the names of columns in file with Ct values (Ct.file).
#' @param sep character of a field separator in both imported files.
#' @param dec character used for decimal points in Ct values.
#'
#' @return data.frame in long format ready to analysis.
#' @export
#'
#' @examples
#' path.Ct.file <- system.file("extdata", "data_Ct_wide.txt", package = "RQdeltaCT")
#' path.design.file <- system.file("extdata", "data_design.txt", package = "RQdeltaCT")
#' data.Ct <- read_Ct_wide(path.Ct.file = path.Ct.file,
#'                    path.design.file = path.design.file,
#'                    sep ="\t",
#'                    dec = ".")
#'
#' @importFrom utils read.csv
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom base which
#' @import tidyverse
#'
read_Ct_wide <- function(path.Ct.file, path.design.file,  sep, dec){
  data_wide <- read.csv(Ct.file,
                        header = TRUE,
                        sep = sep,
                        dec = dec)

  data_wide_design <- read.csv(design.file,
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
#' Imports long-type table with Ct values.
#'
#' @param path path to a .txt file with long-type table of Ct values. This table should contain at least  4 columns, with
#' sample names, target names, Ct values and group names (those columns will be imported by this function).
#' Imported table could also contain a column with flag information, which could be optionally imported (see add.col.Flag and col.Flag parameters).
#'
#' @param sep character of a field separator in both imported files.
#' @param dec character used for decimal points in Ct values.
#' @param skip integer: number of lines of the data file to skip before beginning to read data.
#' @param col.Sample integer: number of column with sample names.
#' @param col.Target integer: number of column with target names.
#' @param col.Ct integer: number of column with Ct values.
#' @param col.Group integer: number of column with group names.
#' @param add.col.Flag logical: if dataset contains a column with flag informatio which should be also imported, this parameters should be set to TRUE.
#' @param col.Flag integer: number of column with flag information.
#' This column should contain a character-type values (ex. "Undetermined" and "OK"), however, other types of values (ex. numeric), could be converted after importing dataset (see examples).
#'
#' @return data.frame in long format ready to analysis.
#' @export
#'
#' @examples
#' path <- system.file("extdata", "data_Ct_long.txt", package = "RQdeltaCT")
#' data.Ct <- read_Ct_long(path = path, sep = "\t",dec = ".",skip = 0,
#'                    add.col.Flag = TRUE, col.Sample = 1, col.Target = 2,
#'                    col.Ct = 5, col.Group = 9, col.Flag = 4)
#'
#' @importFrom utils read.csv
#'
read_Ct_long <- function(path, sep, dec, skip = 0, column.Sample, column.Target, column.Ct, column.Group, add.column.Flag = FALSE, column.Flag){
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
#' This function does not perform data filtering, but only numbers Ct values labeled as reliable or not and presents them graphically.
#' This function Could be useful to identify samples with low number of reliable Ct values.
#'
#' @param data object returned from read_Ct_long or read_Ct_wide function, or data frame containing column named "Sample" with sample names, column named "Target" with target names,
#' column named "Ct" with raw Ct values, column named "Group" with group names. Optionally, data frame could contain column named "Flag" with flag information (ex. "Undetermined" and "OK"),
#' which will be used for relibility assessment.
#' This parameter has no default value (must be provided).
#' @param flag.Ct character of a flag used for undetermined Ct values, default to "Undetermined".
#' @param maxCt numeric, a maximum of Ct value allowed.
#' @param flag character of a flag used in Flag column for values which are unreliable, default to "Undetermined".
#' @param colors character vector length of two, containing colors for Ct values which were labeled as reliable (first element of vector) or not (second element of vector).
#' @param axis.title.size integer: font size of axis titles.
#' @param axis.text.size integer: font size of axis text.
#' @param x.axis.title character: title of x axis.
#' @param y.axis.title character: title of y axis.
#' @param legend.title character: title of legend.
#' @param plot.title character: title of plot.
#' @param legend.title.size integer: font size of legend title.
#' @param legend.text.size integer: font size of legend text.
#' @param legend.position position of the legend, one of "top", "right", "bottom", "left".
#' @param save.to.tiff logical: if plot should be saved as .tiff file, set to TRUE.
#' @param dpi integer: resolution of saved .tiff file.
#' @param width numeric: width (in cm) of saved .tiff file.
#' @param height integer: height (in cm) of saved .tiff file.
#' @param name.tiff character: name of saved .tiff file.
#'
#' @return data.frame in long format ready to analysis.
#' @export
#'
#' @examples
#' sample.Ct.control <- control_Ct_barplot_sample(data.Ct)
#'
#' @importFrom base as.numeric as.data.frame cat print table paste sum colnames
#' @importFrom dplyr mutate arrange filter
#' @import ggplot2
#' @import tidyverse
#'
control_Ct_barplot_sample <- function(data, flag.Ct = "Undetermined", maxCt = 35, flag = "Undetermined",
                           colors = c("#66c2a5", "#fc8d62"),
                           x.axis.title = "",
                           y.axis.title = "Number",
                           axis.title.size = 12,
                           axis.text.size = 10,
                           plot.title = "",
                           legend.title = "Reliable?",
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
  cat("Returned object contains number of targets retained (Yes) or not retained (No) in samples after reliability assessment.")
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
    scale_x_discrete(limits = order)

  print(barplot.samples)

    if (save.to.tiff == TRUE){
      ggsave(paste(name.tiff, ".tiff", sep = ""), barplot.samples, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
    } else {}

    return(table(data$Reliable, data$Sample))
}




#' @title control_Ct_barplot_target
#'
#' @description
#' Target-wide control of raw Ct values across groups by illustrating numbers of Ct values labeled as reliable or not by using reliability criteria (see function parameters).
#' This function does not perform data filtering, but only numbers Ct values labeled as reliable or not and presents them graphically.
#' This function Could be useful to identify targets with low number of reliable Ct values.
#'
#' @param data object returned from read_Ct_long or read_Ct_wide function, or data frame containing column named "Sample" with sample names, column named "Target" with target names,
#' column named "Ct" with raw Ct values, column named "Group" with group names. Optionally, data frame could contain column named "Flag" with flag information (ex. "Undetermined" and "OK"),
#' which will be used for reliability assessment.
#' This parameter has no default value (must be provided).
#' @param flag.Ct character of a flag used for undetermined Ct values, default to "Undetermined".
#' @param maxCt numeric, a maximum of Ct value allowed.
#' @param flag character of a flag used in Flag column for values which are unreliable, default to "Undetermined".
#' @param colors character vector length of two, containing colors for Ct values which were labeled as reliable (first element of vector) or not (second element of vector).
#' @param axis.title.size integer: font size of axis titles.
#' @param axis.text.size integer: font size of axis text.
#' @param x.axis.title character: title of x axis.
#' @param y.axis.title character: title of y axis.
#' @param legend.title character: title of legend.
#' @param plot.title character: title of plot.
#' @param legend.title.size integer: font size of legend title.
#' @param legend.text.size integer: font size of legend text.
#' @param legend.position position of the legend, one of "top", "right", "bottom", "left".
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file.
#' @param dpi integer: resolution of saved .tiff file.
#' @param width numeric: width (in cm) of saved .tiff file.
#' @param height integer: height (in cm) of saved .tiff file.
#' @param name.tiff character: name of saved .tiff file.
#'
#' @return data.frame in long format ready to analysis.
#' @export
#'
#' @examples
#' target.Ct.control <- control_Ct_barplot_target(data.Ct)
#'
#' @importFrom base as.numeric as.data.frame cat print table paste sum colnames
#' @importFrom dplyr mutate arrange filter
#' @import ggplot2
#' @import tidyverse
#'
control_Ct_barplot_target <- function(data, flag.Ct = "Undetermined", maxCt = 35, flag = "Undetermined",
                               colors = c("#66c2a5", "#fc8d62"),
                               x.axis.title = "",
                               y.axis.title = "Number",
                               axis.title.size = 12,
                               axis.text.size = 10,
                               legend.title = "Reliable?",
                               legend.title.size = 12,
                               legend.text.size = 12,
                               legend.position = "top",
                               plot.title = "",
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
    cat("Returned object contains number of samples retained (Yes) or not retained (No) for each target after filtering.")
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
      scale_x_discrete(limits = rev(unique(bar$Var2))) +
      facet_wrap(vars(Var3))

    print(barplot.targets)

    if (save.to.tiff == TRUE){
      ggsave(paste(name.tiff, ".tiff", sep = ""), barplot.targets, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
    }
    return(table(data$Reliable, data$Target, data$Group))

}




#' @title samples_to_remove
#'
#' @description
#' Indicates samples with the amount of unreliable Ct values higher than specified fraction.
#' Could be useful for identification of samples with specified fraction of unreliable Ct values.
#' Samples with high amount of unreliable data could be considered to remove from the dataset.
#'
#' @param data object returned from control_Ct_barplot_sample function.
#' @param fraction numeric from 0 to 1: a threshold fraction, samples with the higher fraction of Ct values labeled as unreliable will be returned. Default to 0.5.
#' @return vector with names of samples.
#' @export
#'
#' @examples
#' samples.to.remove <- samples_to_remove(sample.Ct.control, 0.5)
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
#' Targets with high amount of unreliable data could be considered to remove from the dataset.
#'
#' @param data object returned from control_Ct_barplot_target function.
#' @param fraction numeric from 0 to 1: a threshold fraction, targets with the higher fraction of Ct values labeled as unreliable in at least one of the compared groups will be returned. Default to 0.5.
#' @param groups character vector length of two, containing names of compared groups.
#'
#' @return vector with names of targets.
#' @export
#'
#' @examples
#' targets.to.remove <- targets_to_remove(sample.Ct.control, 0.5, groups = c("Disease", "Control"))
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
#' Filters Ct dataset according to the used filtering criteria.
#'
#' @param data object returned from read_Ct_long or read_Ct_wide function, or data frame containing column named "Sample" with sample names, column named "Target" with target names,
#' column named "Ct" with raw Ct values, column named "Group" with group names. Optionally, data frame could contain column named "Flag" with flag information (ex. "Undetermined" and "OK"),
#' which will be used for filtering.
#'
#' @param flag.Ct character of a flag used for undetermined Ct values, default to "Undetermined".
#' @param maxCt numeric, a maximum of Ct value allowed.
#' @param flag character of a flag used in Flag column for values which should be filtered out, default to "Undetermined".
#' @param remove.Target vector with names of targets which should be removed from dataset
#' @param remove.Sample vector with names of samples which should be removed from dataset
#' @param remove.Group vector with names of groups which should be removed from dataset
#'
#' @return data.frame with filtered data.
#' @export
#'
#' @examples
#'
#' data.CtF <- filter_Ct(data.Ct,
#'                       flag.Ct = "Undetermined",
#'                       maxCt = 35,
#'                       flag = "Undetermined",
#'                       remove.Target = c("Gene2","Gene6"),
#'                       remove.Sample = c("Control8","control16"))
#'
#'  data.CtF <- filter_Ct(data.Ct,
#'                        flag.Ct = "Undetermined",
#'                        maxCt = 35,
#'                        flag = "Undetermined",
#'                        remove.Target = targets.to.remove,
#'                        remove.Sample = samples.to.remove)
#'
#' @importFrom base as.numeric sum colnames
#' @importFrom dplyr filter
#' @import tidyverse
#'
filter_Ct <- function(data, flag.Ct = "Undetermined", maxCt = 35,
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






#' @title exp_Ct
#'
#' @description
#' Collapses technical replicates (if present in data) by means, counts and imputes missing data by means within groups (if so indicated),
#' and finally exponentiates Ct values by using formula 2^(-Ct).
#'
#' @param data data object returned from read_Ct_long, read_Ct_wide or filter_Ct function, or data frame containing column named "Sample" with sample names, column named "Target" with target names,
#' column named "Ct" with raw Ct values, column named "Group" with group names. Any other columns could exist, by will not be used by this function.
#' @param imput.by.mean.within.groups logical: if TRUE, missing values will be imputed by means within groups.
#' @param save.to.txt logical: if TRUE, returned dataset will be saved to .txt file.
#' @param name.txt character: name of saved .txt file.
#'
#' @return data.frame with transformed Ct values and printed information about number and percentage of missing values.
#' @export
#'
#' @examples
#' data.Ct.exp.imput <- exp_Ct(data.CtF,
#'                             imput.by.mean.within.groups = TRUE)
#'
#' @importFrom base as.data.frame as.vector sum is.na ncol nrow paste cat
#' @importFrom utils write.table
#' @importFrom dplyr mutate filter select
#' @import tidyverse
#'
exp_Ct <- function(data,
                  imput.by.mean.within.groups = FALSE,
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
#' Performs relative quantification of gene expression using 2^(-Ct) method.
#' This function calculates:
#' 1. Means (in columns with "_mean" pattern) and standard deviations (in columns with "_sd" pattern) of exponentiated Ct values of analyzed targets across compared groups.
#' 2. Normality tests (Shapiro_Wilk test) of exponentiated Ct values of analyzed targets across compared groups and returned p values in columns with "_norm_p" pattern.
#' 3. Fold Change (in "FCh" column) values together with log10 Fold change values (in "log10FCh" column).
#'    Fold change values were calculated for each target by dividing  mean of exponentiated Ct values in study group by mean of exponentiated Ct values in reference group.
#' 4. Statistical testing of differences in exponentiated Ct values between study group and reference group.
#'    Student's t test and Mann-Whitney U test were used and resulted statistics (in column with "_test_stat" pattern) and p values (in column with "_test_p" pattern) are returned.
#'
#' @param data data object returned from exp_Ct function.
#' @param group.study character: name of study group (group of interest).
#' @param group.ref character: name of reference group.
#' @param do.tests logical: if TRUE, statistical significance of differences in exponentiated Ct values between compared groups will be calculated using Student's t test and Mann-Whitney U test.
#' @param save.to.txt logical: if TRUE, returned dataset with results will be saved to .txt file.
#' @param name.txt character: name of saved .txt file.
#
#' @return data.frame with transformed Ct values and printed information about number and percentage of missing values.
#' @export
#'
#' @examples
#' RQ.exp.Ct.results <- RQ_exp_Ct(data.Ct.exp.imput,
#'                                group.study = "Disease",
#'                                group.ref = "Control",
#'                                do.tests = TRUE)
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
                     name.txt = "ddCt_RQ_results"){

  data_slim <- data %>%
    filter(Group == group.study | Group == group.ref) %>%
    pivot_longer(cols = -c(Group, Sample), names_to = "Target", values_to = "expCt")

  data_slim$Group <- as.factor(data_slim$Group)

  data_FCh <- data_slim %>%
    group_by(Group, Target) %>%
    summarise(expCt = mean(expCt, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = expCt) %>%
    mutate(FCh = .data[[group.study]]/.data[[group.ref]]) %>%
    mutate(log10FCh = log10(FCh)) %>%
    as.data.frame()

  data_FCh_sd <- data_slim %>%
    group_by(Group, Target) %>%
    summarise(expCt_sd = sd(expCt, na.rm = TRUE), .groups = "keep") %>%
    pivot_wider(names_from = Group, values_from = expCt_sd) %>%
    rename_with(~paste0(.x, "_sd", recycle0 = TRUE), all_of(c(group.study, group.ref)))

  if (do.tests == TRUE){

    data_FCh_norm <- data_slim %>%
      group_by(Group, Target) %>%
      summarise(shap_wilka_p = shapiro.test(expCt)$p.value, .groups = "keep") %>%
      pivot_wider(names_from = Group, values_from = shap_wilka_p) %>%
      rename_with(~paste0(.x, "_norm_p", recycle0 = TRUE), all_of(c(group.study, group.ref))) %>%
      full_join(data_FCh, by = c("Target")) %>%
      rename_with(~paste0(.x, "_mean", recycle0 = TRUE), all_of(c(group.study, group.ref)))

    data_FCh_tests <- data_slim %>%
      group_by(Target) %>%
      summarise(t_test_p = t.test(expCt ~ Group, alternative = "two.sided")$p.value,
                t_test_stat = t.test(expCt ~ Group, alternative = "two.sided")$statistic,
                MW_test_p = pvalue(wilcox_test(expCt ~ Group)),
                MW_test_stat = statistic(wilcox_test(expCt ~ Group)), .groups = "keep")
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
#' Calculate line plot and statistics (minimum, maximum, standard deviation, variance and colinearity coefficient VIF) that could be helpful to select the best reference gene for normalization of CT values.
#'
#' @param data object returned from read_Ct_long, read_Ct_wide or filter_Ct function, or data frame containing column named "Sample" with sample names, column named "Target" with target names,
#' column named "Ct" with raw Ct values, column named "Group" with group names. Any other columns could exist, by will not be used by this function.
#'
#' @param candidates vector of names of targets - candidates for gene reference.
#' @param colors vector of colors for targets, number of colors should be equal to number of candidate genes (elements in candidates vector).
#' @param line.width numeric: width of lines drawn in the plot.
#' @param angle integer: value of angle in which names of genes should be displayed.
#' @param x.axis.title character: title of x axis.
#' @param y.axis.title character: title of y axis.
#' @param axis.title.size integer: font size of axis titles.
#' @param axis.text.size integer: font size of axis text.
#' @param legend.title character: title of legend.
#' @param legend.title.size integer: font size of legend title.
#' @param legend.text.size integer: font size of legend text.
#' @param legend.position position of the legend, one of the following: "top", "right", "bottom", "left".
#' @param plot.title character: title of plot.
#' @param save.to.tiff logical: if TRUE, plot will be saved as .tiff file.
#' @param dpi integer: resolution of saved .tiff file.
#' @param width numeric: width (in cm) of saved .tiff file.
#' @param height integer: height (in cm) of saved .tiff file.
#' @param name.tiff character: name of saved .tiff file.
#'
#' @return data.frame in long format ready to analysis.
#' @export
#'
#' @examples
#' ref <- sel_ref(data.CtF,
#'                candidates = c("Gene4", "Gene8","Gene10","Gene16","Gene17", "Gene18"),
#'                col = c("#66c2a5", "#fc8d62","#6A6599", "#D62728", "#1F77B4", "black"),
#'                line.width = 1,
#'                angle = 35)
#'
#' @importFrom base mean print min max as.data.frame
#' @importFrom stats sd var
#' @importFrom dplyr filter
#' @importFrom car vif
#' @import ggplot2
#' @import tidyverse
#'
select_ref_gene <- function(data, candidates,
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
                    dpi = 600, width = 15, height = 15,
                    save.to.tiff = FALSE,
                    name.tiff = "Ct_reference_selection"){

  ref <- data %>%
    group_by(Group, Target, Sample) %>%
    summarise(mean = mean(Ct, na.rm = TRUE)) %>%
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
    theme(legend.text = element_text(size = legend.text.size, colour="black"))

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
              var = var(mean, na.rm = TRUE)) %>%
    as.data.frame()

  ref_vif <- data %>%
    group_by(Group, Target, Sample) %>%
    summarise(mean = base::mean(Ct, na.rm = TRUE)) %>%
    ungroup()

  ref_lm <- pivot_wider(ref_vif, names_from = Target, values_from = mean)

  ref_lm$dum <- c(1:nrow(ref_lm))
  model <- lm(dum ~ ., data = select(ref_lm, -Sample, -Group))
  vif <- vif(model)
  vif_sel <- vif[names(vif) %in% candidates]
  ref_var$VIF <- vif_sel
  ref_var$VIF.Target <- names(vif_sel)

  return(ref_var)
}






deltaCt <- function(data,
                    imput.by.mean.within.groups = FALSE,
                    ref, save.to.txt = FALSE,
                    name.txt = "dCt_results"){
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





control_dCt_boxplot_sample <- function(data, coef = 1.5,
                               col = c("#66c2a5", "#fc8d62"),
                               axis.title.size = 12,
                               axis.text.size = 12,
                               x.axis.title = "Sample", y.axis.title = "dCt",
                               legend.title = "Group",
                               legend.title.size = 12,
                               legend.text.size = 12,
                               legend.position = "right",
                               plot.title = "",
                               dpi = 600, width = 15, height = 15,
                               save.to.tiff = FALSE,
                               name.tiff = "dCt_control_boxplot"){

    data <- pivot_longer(data, !c(Sample, Group), names_to = "Target" , values_to = "dCt")
    box_control <- ggplot(data, aes(x = Sample, y = dCt, color = Group)) +
      geom_boxplot(coef = coef) +
      scale_x_discrete(limits = rev(unique(data$Sample))) +
      scale_color_manual(values = c(col)) +
      coord_flip() +
      xlab(x.axis.title) + ylab(y.axis.title) +
      labs(color = legend.title, title = plot.title) +
      theme_classic() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black"))
    print(box_control)
    if (save.to.tiff == TRUE){
      ggsave(paste(name.tiff,".tiff", sep = ""), box_control, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
    }
}





control_dCt_boxplot_target <- function(data, coef = 1.5, by.group = FALSE,
                                       col = c("#66c2a5", "#fc8d62"),
                                       axis.title.size = 12,
                                       axis.text.size = 12,
                                       x.axis.title = "Target", y.axis.title = "dCt",
                                       legend.title = "Group",
                                       legend.title.size = 12,
                                       legend.text.size = 12,
                                       legend.position = "right",
                                       plot.title = "",
                                       dpi = 600, width = 15, height = 15,
                                       save.to.tiff = FALSE,
                                       name.tiff = "dCt_control_boxplot"){

  data <- pivot_longer(data, !c(Sample, Group), names_to = "Target" , values_to = "dCt")
  if (by.group == TRUE){
    box_control <- ggplot(data, aes(x = Target, y = dCt, color = Group)) +
      geom_boxplot(coef = coef) +
      scale_x_discrete(limits = rev(unique(data$Target))) +
      scale_color_manual(values = c(col)) +
      coord_flip() +
      xlab(x.axis.title) + ylab(y.axis.title) +
      labs(color = legend.title, title = plot.title) +
      theme_classic() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black"))
  } else {
    box_control <- ggplot(data, aes(x = Target, y = dCt)) +
      geom_boxplot(coef = coef, fill = col[1]) +
      scale_x_discrete(limits = rev(unique(data$Target))) +
      coord_flip() +
      xlab(x.axis.title) + ylab(y.axis.title) +
      labs(title = plot.title) +
      theme_classic() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black"))

  }
  print(box_control)
  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), box_control, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
}






control_dCt_cluster_sample <- function(data, method.dist = "euclidean", method.clust = "complete",
                                x.axis.title = "Samples", y.axis.title = "Height",
                                plot.title = "",
                                dpi = 600, width = 15, height = 15,
                                save.to.tiff = FALSE,
                                name.tiff = "dCt_control_clust_plot"){
  cluster <- hclust(dist(data, method = method.dist), method = method.clust)
  cluster$labels <- data$Sample
  plot(cluster, xlab  = x.axis.title, ylab = y.axis.title, main = plot.title)
  if (save.to.tiff == TRUE){
    tiff(paste(name.tiff, ".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
    plot(cluster, xlab  = x.axis.title, ylab = y.axis.title, main = plot.title)
    dev.off()
  }
}





control_dCt_cluster_target <- function (data, method.dist = "euclidean", method.clust = "single",
                                        x.axis.title = "Samples", y.axis.title = "Height", plot.title = "",
                                        dpi = 600, width = 15, height = 15, save.to.tiff = FALSE,
                                        name.tiff = "dCt_control_clust_plot")
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





control_dCt_pca_sample <- function(data, point.size = 4, alpha = 0.7, label.size = 3, hjust = 0, vjust = -1,
                            col = c("#66c2a5", "#fc8d62"), point.shape = 19,
                            axis.title.size = 12,
                            axis.text.size = 10,
                            legend.text.size = 12,
                            legend.title = "Group",
                            legend.title.size = 12,
                            legend.position = "right",
                            plot.title = "",
                            dpi = 600, width = 15, height = 15,
                            save.to.tiff = FALSE,
                            name.tiff = "dCt_control_pca_plot"){
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
    scale_color_manual(values = c(col)) +
    labs(colour = legend.title, title = plot.title) +
    theme_bw() +
    labs(x = paste("PC1: ", round(var_pca1*100,2), "% variance explained", sep = ""),
         y = paste("PC2: ", round(var_pca2*100,2), "% variance explained", sep = "")) +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour="black")) +
    theme(legend.text = element_text(size = legend.text.size, colour="black")) +
    theme(legend.title = element_text(size = legend.title.size, colour="black")) +
    geom_text(aes(label = Sample), hjust = hjust, vjust = vjust, size = label.size)
  print(control_pca)
  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), control_pca, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
}




control_dCt_pca_target <- function(data, point.size = 4, alpha = 0.7, label.size = 3, hjust = 0, vjust = -1, shape = 19,
                                   col = c("#66c2a5", "#fc8d62"), point.shape = 19,
                                   axis.title.size = 12,
                                   axis.text.size = 10,
                                   legend.text.size = 12,
                                   legend.title = "Group",
                                   legend.title.size = 12,
                                   legend.position = "right",
                                   plot.title = "",
                                   dpi = 600, width = 15, height = 15,
                                   save.to.tiff = FALSE,
                                   name.tiff = "dCt_control_pca_plot"){
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
    geom_point(size = point.size, shape = point.shape, alpha = alpha, col = col[1]) +
    labs(title = plot.title) +
    theme_bw() +
    labs(x = paste("PC1: ", round(var_pca1*100,2), "% variance explained", sep = ""),
         y = paste("PC2: ", round(var_pca2*100,2), "% variance explained", sep = "")) +
    theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour="black")) +
    geom_text(aes(label = Target), hjust = hjust, vjust = vjust, size = label.size)
  print(control_pca)
  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), control_pca, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
}





control_dCt_corr_target <- function(data, add.coef = "black",
                             method = "pearson",
                             order = "hclust",
                             hclust.method = "complete",
                             size = 0.6,
                             save.to.tiff = FALSE,
                             dpi = 600, width = 15, height = 15,
                             name.tiff = "dCt_control_corr_plot",
                             save.to.txt = FALSE,
                             sort.txt = "cor",
                             name.txt = "dCt_control_corr_results",
                             p.adjust.method = "BH"){
  data <- as.data.frame(data)
  data <- select(data, -Group, -Sample)
  res_cor <- rcorr(as.matrix(data), type = method)
  if (order == "hclust"){
    corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
             order = order,
             hclust.method = hclust.method,
             addCoef.col = add.coef,
             number.cex = size)
    if (save.to.tiff == TRUE){
      tiff(paste(name.tiff,".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
      corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
               order = order,
               hclust.method = hclust.method,
               addCoef.col = "black",
               number.cex = size)
      dev.off()
    }
  } else {
    corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
             order = order,
             addCoef.col = add.coef,
             number.cex = size)
    if (save.to.tiff == TRUE){
      tiff(paste(name.tiff,".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
      corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
               order = order,
               addCoef.col = "black",
               number.cex = size)
      dev.off()
    }
  }
  if (save.to.txt == TRUE){
    corr <- upper.tri(res_cor$r)
    corr.data <- data.frame(
      row = rownames(res_cor$r)[row(res_cor$r)[corr]],
      column = rownames(res_cor$r)[col(res_cor$r)[corr]],
      cor  =(res_cor$r)[corr],
      p = res_cor$P[corr]
    )
    corr.data$p.adj <- p.adjust(corr.data$p, method = p.adjust.method)
    corr.data.sort <- arrange(corr.data, -abs(cor))
    write.table(corr.data.sort, paste(name.txt, ".txt", sep = ""))
  }
  return(corr.data.sort)
}





control_dCt_corr_sample <- function(data, add.coef = "black",
                                    method = "pearson",
                                    order = "hclust",
                                    hclust.method = "complete",
                                    size = 0.6,
                                    save.to.tiff = FALSE,
                                    dpi = 600, width = 15, height = 15,
                                    name.tiff = "dCt_control_corr_plot",
                                    save.to.txt = FALSE,
                                    sort.txt = "cor",
                                    name.txt = "dCt_control_corr_results",
                                    p.adjust.method = "BH"){
  data <- as.data.frame(data)
  data_t <- t(select(data, -Group, -Sample))
  colnames(data_t) <- data$Sample
  res_cor <- rcorr(as.matrix(data_t), type = method)
  if (order == "hclust"){
    corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
             order = order,
             hclust.method = hclust.method,
             addCoef.col = add.coef,
             number.cex = size)
    if (save.to.tiff == TRUE){
      tiff(paste(name.tiff,".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
      corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
               order = order,
               hclust.method = hclust.method,
               addCoef.col = "black",
               number.cex = size)
      dev.off()
    }
  } else {
    corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
             order = order,
             addCoef.col = add.coef,
             number.cex = size)
    if (save.to.tiff == TRUE){
      tiff(paste(name.tiff,".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
      corrplot(res_cor$r, type = "upper", tl.cex = size, tl.col = "black", cl.cex = size,
               order = order,
               addCoef.col = "black",
               number.cex = size)
      dev.off()
    }
  }
  if (save.to.txt == TRUE){
    corr <- upper.tri(res_cor$r)
    corr.data <- data.frame(
      row = rownames(res_cor$r)[row(res_cor$r)[corr]],
      column = rownames(res_cor$r)[col(res_cor$r)[corr]],
      cor  =(res_cor$r)[corr],
      p = res_cor$P[corr]
    )
    corr.data$p.adj <- p.adjust(corr.data$p, method = p.adjust.method)
    corr.data.sort <- arrange(corr.data, -abs(cor))
    write.table(corr.data.sort, paste(name.txt, ".txt", sep = ""))
  }
  return(corr.data.sort)
}








single_dCt_corr <- function(data, x, y, alpha = 0.7,
                            by.group = FALSE,
                            col = c("#66c2a5", "#fc8d62"),
                            axis.title.size = 12,
                            axis.text.size = 10,
                            legend.text.size = 12,
                            legend.title = "Group",
                            legend.title.size = 12,
                            legend.position = "right",
                            plot.title = "",
                            label.position.x = 1,
                            label.position.y = 1,
                            small.p = FALSE,
                            small.r = FALSE,
                            p.digits = 3,
                            rr.digits = 2,
                            label = c("eq", "R2", "p"),
                            dpi = 600, width = 15, height = 15,
                            save.to.tiff = FALSE,
                            name.tiff = "dCt_corr_single_plot"){
  if (by.group == TRUE){
    corr_control <- ggplot(data, aes(x = .data[[x]], y = .data[[y]], color = Group)) +
      geom_point() +
      geom_smooth(method='lm', se = FALSE, alpha = alpha) +
      scale_color_manual(values = c(col)) +
      xlab(x) + ylab(y) +
      labs(color = legend.title, title = plot.title) +
      theme_classic() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black"))
    if (labels == TRUE){
      corr_control <- corr_control +
        stat_poly_eq(use_label(label),
                     label.y = c(label.position.y),
                     label.x = c(label.position.x),
                     small.p = small.p,
                     small.r = small.r,
                     p.digits = p.digits,
                     rr.digits = rr.digits)
    }
  } else {
    corr_control <- ggplot(data, aes(x = .data[[x]], y = .data[[y]])) +
      geom_point() +
      geom_smooth(method='lm', se = FALSE, alpha = alpha) +
      stat_poly_eq(use_label(label),
                   label.y = c(label.position.y),
                   label.x = c(label.position.x)) +
      xlab(x) + ylab(y) +
      labs(title = plot.title) +
      theme_classic() +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black"))
    if (labels == TRUE){
      corr_control <- corr_control +
        stat_poly_eq(use_label(label),
                     label.y = c(label.position.y),
                     label.x = c(label.position.x),
                     small.p = small.p,
                     small.r = small.r,
                     p.digits = p.digits,
                     rr.digits = rr.digits)
      }
    }
  print(corr_control)
  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), corr_control, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
}



single_dCt_corr_target <- function(data, x, y, alpha = 0.7, labels = TRUE,
                                   col = c("#66c2a5", "#fc8d62"),
                                   axis.title.size = 12,
                                   axis.text.size = 10,
                                   legend.text.size = 12,
                                   legend.title = "Group",
                                   legend.title.size = 12,
                                   legend.position = "right",
                                   plot.title = "",
                                   label.position.x = 1,
                                   label.position.y = 1,
                                   small.p = FALSE,
                                   small.r = FALSE,
                                   p.digits = 3,
                                   rr.digits = 2,
                                   label = c("eq", "R2", "p"),
                                   dpi = 600, width = 15, height = 15,
                                   save.to.tiff = FALSE,
                                   name.tiff = "dCt_corr_single_plot"){
  data <- as.data.frame(data)
  data_t <- t(select(data, -Group, -Sample))
  colnames(data_t) <- data$Sample
  corr_control <- ggplot(as.data.frame(data_t), aes(x = .data[[x]], y = .data[[y]])) +
    geom_point() +
    geom_smooth(method='lm', se = FALSE, alpha = alpha) +
    xlab(x) + ylab(y) +
    labs(title = plot.title) +
    theme_classic() +
    theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
    theme(axis.title = element_text(size = axis.title.size, colour="black"))
  if (labels == TRUE){
    corr_control <- corr_control +
      stat_poly_eq(use_label(label),
                   label.y = c(label.position.y),
                   label.x = c(label.position.x),
                   small.p = small.p,
                   small.r = small.r,
                   p.digits = p.digits,
                   rr.digits = rr.digits)
  }
  print(corr_control)
  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff,".tiff", sep = ""), corr_control, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
}






filter_dCt <- function(data,
                       remove.Target = c(""),
                       remove.Sample = c(""),
                       remove.Group = c("")){
  data <- dplyr::filter(data,
                        !Sample %in% remove.Sample,
                        !Group %in% remove.Group)
  data <- select(data, -any_of(remove.Target))
  return(data)
}






exp_dCt <- function(data,
                    save.to.txt = FALSE,
                    name.txt = "expdCt"){

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



results_dCt_boxplot <- function(data, coef = 1.5, sel.Target = "all", by.group = TRUE,
                                faceting = TRUE, facet.row = 4, facet.col = 4,
                                angle = 0, rotate = FALSE, add.mean = FALSE,
                                add.mean.size = 2,
                                add.mean.color = "black",
                                col = c("#66c2a5", "#fc8d62"),
                                axis.title.size = 12,
                                axis.text.size = 10,
                                legend.text.size = 12,
                                x.axis.title = "Sample",
                                y.axis.title = "dCt",
                                legend.title = "Group",
                                legend.title.size = 12,
                                legend.position = "right",
                                plot.title = "",
                                dpi = 600, width = 15, height = 15,
                                save.to.tiff = FALSE,
                                name.tiff = "dCt_results_boxplot"){

  data <- pivot_longer(data, !c(Sample, Group), names_to = "Target" , values_to = "dCt")
  if (sel.Target[1] == "all"){
    data <- data
  } else {
    data <- filter(data, Target %in% sel.Target)
  }
  if (by.group == TRUE){
    box_results <- ggplot(data, aes(x = Target, y = dCt, fill = Group)) +
      geom_boxplot(coef = coef) +
      #scale_x_discrete(limits = rev(unique(data$Sample))) +
      scale_fill_manual(values = c(col)) +
      #coord_flip() +
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

  } else {
    box_results <- ggplot(data, aes(x = Target, y = dCt)) +
      geom_boxplot(coef = coef, fill = col[1]) +
      #scale_x_discrete(limits = rev(unique(data$Sample))) +
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
}








RQ_ddCt <- function(data,
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





RQ_plot_expCt_or_expdCt <- function(data, use.p = TRUE, mode, Target.sel = "all", p.threshold = 0.05,
                                    col.width = 0.8,
                                    angle = 0, rotate = FALSE,
                                    col = c("#66c2a5", "#fc8d62"),
                                    axis.title.size = 12,
                                    axis.text.size = 10,
                                    legend.text.size = 12,
                                    x.axis.title = "Target",
                                    y.axis.title = "log10(Fold change)",
                                    legend.title = "Statistically significant?",
                                    legend.title.size = 12,
                                    legend.position = "top",
                                    plot.title = "",
                                    dpi = 600, width = 15, height = 15,
                                    save.to.tiff = FALSE,
                                    name.tiff = "ddCt_RQ_plot"){
  if (Target.sel[1] != "all"){
    data <- filter(data, Target %in% Target.sel)
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
    data <- mutate(data, `Statistically significant?` = ifelse(p.used > p.threshold, yes = "No (p > 0.05)",  no = "Yes (p <= 0.05)"))
    data$`Statistically significant?` <- factor(data$`Statistically significant?`, levels = c("Yes (p <= 0.05)", "No (p > 0.05)"))
    RQ <- ggplot(data, aes(x = reorder(Target, -FCh), y = log10(FCh), fill = `Statistically significant?`)) +
      geom_col(width = col.width) +
      scale_fill_manual(values = c(col)) +
      xlab(x.axis.title) + ylab(y.axis.title) +
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
      geom_col(width = col.width, fill = col[1]) +
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
  return(data)
}








RQ_plot <- function(data, use.p = TRUE, mode, Target.sel = "all", p.threshold = 0.05,
                    col.width = 0.8,
                    angle = 0, rotate = FALSE,
                    col = c("#66c2a5", "#fc8d62"),
                    axis.title.size = 12,
                    axis.text.size = 10,
                    legend.text.size = 12,
                    x.axis.title = "Target",
                    y.axis.title = "log10(RQ)",
                    legend.title = "Statistically significant?",
                    legend.title.size = 12,
                    legend.position = "top",
                    plot.title = "",
                    dpi = 600, width = 15, height = 15,
                    save.to.tiff = FALSE,
                    name.tiff = "ddCt_RQ_plot"){
  if (Target.sel[1] != "all"){
    data <- filter(data, Target %in% Target.sel)
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
  data <- mutate(data, `Statistically significant?` = ifelse(p.used > p.threshold, yes = "No (p > 0.05)",  no = "Yes (p <= 0.05)"))
  data$`Statistically significant?` <- factor(data$`Statistically significant?`, levels = c("Yes (p <= 0.05)", "No (p > 0.05)"))
  RQ <- ggplot(data, aes(x = reorder(Target, -RQ), y = log10(RQ), fill = `Statistically significant?`)) +
    geom_col(width = col.width) +
    scale_fill_manual(values = c(col)) +
    xlab(x.axis.title) + ylab(y.axis.title) +
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
    RQ <- ggplot(data, aes(x = reorder(Target, -RQ), y = log10(RQ))) +
      geom_col(width = col.width, fill = col[1]) +
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
  return(data)
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

