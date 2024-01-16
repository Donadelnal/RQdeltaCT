
#' @title read_Ct_wide
#'
#' @description
#' Imports Ct table in wide format, imports design file and merge both files to return long table ready for analysis.
#' All parameters of this function have not default values and must be specified by user.
#'
#'
#' @param Ct.file path to wide-format table (must be in .txt format) containing target names in the first column and
#' Ct values in remaining columns named with sample names.
#' @param design.file A .txt file with two columns: column named "Sample" with names of samples
#' and column named "Group" with names of groups assigned to samples. Names of samples in this file
#' must correspond to the names of columns in file with Ct values.
#' @param sep A character of a field separator in both imported files.
#' @param dec A character used for decimal points in Ct values.
#'
#' @return A table in long format ready to analysis.
#' @export
#'
#' @examples
#' data_long <- prepare_data(Ct.file = "Ct_wide.txt",
#'                           design.file = "design.txt",
#'                           sep ="\t",
#'                           dec = ".")
#'
#' @importFrom utils read.csv
#' @importFrom tidyr pivot_longer
#'
read_Ct_wide <- function(Ct.file, design.file,  sep, dec){
  data_wide <- read.csv(Ct.file,
                        header = TRUE,
                        sep = sep,
                        dec = dec)

  data_wide_design <- read.csv(design.file,
                               header = TRUE,
                               sep = sep)

  colnames(data_wide)[1] <- "Sample"
  data_wide <- mutate(data_wide, across(everything(), as.character))
  data_slim <- pivot_longer(data_wide, -Sample, names_to = "Target", values_to = "Ct")
  data_slim[ ,"Group"] <- NA

  for (x in 1:nrow(data_wide_design)) {
    index <- which(data_slim$Sample == data_wide_design$Sample[x])
    data_slim$Group[index] <- data_wide_design$Group[x]
  }
  return(data_slim)
}





#' Title
#'
#' @param path
#' @param sep
#' @param dec
#' @param skip
#' @param col.Sample
#' @param col.Target
#' @param col.Ct
#' @param col.Group
#' @param add.col.Flag
#' @param col.Flag
#'
#' @return
#' @export
#'
#' @examples
read_Ct_long <- function(path, sep, dec, skip = 0, col.Sample, col.Target, col.Ct, col.Group, add.col.Flag = FALSE, col.Flag){
  data <- read.csv(path,
    header = TRUE,
    sep = sep, dec = dec, skip = skip)
  if (add.col.Flag == FALSE){
  data <- data[ ,c(col.Sample, col.Target, col.Ct, col.Group)]
  colnames(data) <- c("Sample", "Target", "Ct", "Group")
  }
  if (add.col.Flag == TRUE){
    data <- data[ ,c(col.Sample, col.Target, col.Ct, col.Group, col.Flag)]
    colnames(data) <- c("Sample", "Target", "Ct", "Group", "Flag")
  }
  return(data)
}





control_Ct_barplot <- function(data, mode, flag.Ct = "Undetermined", maxCt = 35, flag = "Undetermined",
                           col = c("#66c2a5", "#fc8d62"), axis.title.size = 12, axis.text.size = 10,
                           x.axis.title = "", y.axis.title = "Number",
                           legend.title = "Retained for analysis?",
                           plot.title = "", legend.title.size = 12,
                           legend.text.size = 12,
                           legend.position = "top",
                           dpi = 600, width = 15, height = 15,
                           save.to.tiff = FALSE,
                           name.tiff = "Ct_control_barplot"){

  data$Ct[data$Ct == flag.Ct] <- 40
  data$Ct <- as.numeric(data$Ct)
  if(sum(colnames(data) %in% "Flag") > 0){
    data <- mutate(data, Detected = ifelse(Ct > maxCt | Flag == flag, yes = "No",  no = "Yes"))
  } else {
  data <- mutate(data, Detected = ifelse(Ct > maxCt , yes = "No",  no = "Yes"))
  }
  if (mode == "Sample"){
  bar <- as.data.frame(table(data$Detected, data$Sample))
  order <- arrange(filter(bar, Var1 == "Yes"), Freq)$Var2
  cat("Returned object contains number of targets retained (Yes) or not retained (No) in samples after filtering.")
  barr <- ggplot(bar, aes(x = reorder(Var2, desc(Freq)), y = Freq, fill = Var1)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(breaks = c("Yes", "No"), values = c("Yes" = col[1],"No" = col[2])) +
    xlab(x.axis.title) + ylab(y.axis.title) +
    labs(fill = legend.title, title = plot.title) +
    theme_classic() + theme(legend.position = legend.position) +
    theme(axis.text = element_text(size = axis.text.size, color = 'black')) +
    theme(axis.title = element_text(size = axis.title.size, color = 'black')) +
    theme(legend.title = element_text(size = legend.title.size, colour="black")) +
    theme(legend.text = element_text(size = legend.text.size, colour="black")) +
    scale_x_discrete(limits = order)

  print(barr)

    if (save.to.tiff == TRUE){
      ggsave(paste(name.tiff, ".tiff", sep = ""), barr, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
    } else {}

    return(table(data$Detected, data$Sample))

  } else {

    bar <- as.data.frame(table(data$Detected, data$Target, data$Group))
    order <- arrange(filter(bar, Var1 == "Yes"), Freq)$Var2
    cat("Returned object contains number of samples retained (Yes) or not retained (No) for each target after filtering.")
    barr <- ggplot(bar, aes(x = reorder(Var2, desc(Freq)), y = Freq, fill = Var1)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(breaks = c("Yes", "No"), values = c("Yes" = col[1],"No" = col[2])) +
      xlab(x.axis.title) + ylab(y.axis.title) +
      labs(fill = legend.title, title = plot.title) +
      theme_classic() + theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, color = 'black')) +
      theme(axis.title = element_text(size = axis.title.size, color = 'black')) +
      theme(legend.title = element_text(size = legend.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
      scale_x_discrete(limits = rev(unique(bar$Var2))) +
      facet_wrap(vars(Var3))
    print(barr)

  if (save.to.tiff == TRUE){
    ggsave(paste(name.tiff, ".tiff", sep = ""), barr, dpi = dpi, width = width, height = height, units = "cm", compression = "lzw")
  }
    return(table(data$Detected, data$Target, data$Group))
    }
}




samples_to_remove <- function(data, fraction = 0.5){
  data <- as.data.frame(data)
  n <- sum(data$Freq[c(1:2)])
  data <- dplyr::filter(data, Var1 == "No" & Freq > n*fraction)
  return(as.vector(data$Var2))
}




targets_to_remove <- function(data, fraction = 0.5, groups){
  data <- as.data.frame(data)
  data1 <- dplyr::filter(data, Var3 == groups[1])
  data2 <- dplyr::filter(data, Var3 == groups[2])
  n1 <- sum(data1$Freq[c(3:4)])
  n2 <- sum(data2$Freq[c(3:4)])
  data1F <- dplyr::filter(data1, Var1 == "No" & Freq > n1*fraction)
  data2F <- dplyr::filter(data2, Var1 == "No" & Freq > n2*fraction)
  return(unique(c(as.vector(data1F$Var2), as.vector(data2F$Var2))))
}





filter_Ct <- function(data, flag.Ct = "Undetermined", maxCt, flag = c(""),
                      remove.Target = c(""), remove.Sample = c(""), remove.Group = c("")){
  data <- dplyr::filter(data, Ct != flag.Ct)
  data$Ct <- as.numeric(data$Ct)
  data <- dplyr::filter(data, Ct <= maxCt,
                        !Target %in% remove.Target,
                        !Sample %in% remove.Sample,
                        !Group %in% remove.Group)
  if(sum(colnames(data) %in% "Flag") > 0){
    data <- dplyr::filter(data,
                          !Flag %in% flag)
  }
  return(data)
}






sel_ref <- function(data, candidates, line.width = 1,
                    col, angle = 0,
                    axis.title.size = 12,
                    axis.text.size = 10,
                    x.axis.title = "",
                    y.axis.title = "Ct",
                    legend.title = "",
                    plot.title = "",
                    legend.title.size = 12,
                    legend.text.size = 12,
                    legend.position = "top",
                    dpi = 600, width = 15, height = 15,
                    save.to.tiff = FALSE,
                    name.tiff = "Ct_reference_selection"){

  ref <- data %>%
    group_by(Group, Target, Sample) %>%
    summarise(mean = base::mean(Ct, na.rm = TRUE)) %>%
    dplyr::filter(Target %in% candidates)
  ref_plot <- ggplot(ref, aes(x = Sample, y = mean, color = Target, group = Target)) +
    geom_line(linewidth = line.width) +
    scale_color_manual(values = c(col)) +
    theme_bw() +
    xlab(x.axis.title) + ylab(y.axis.title) +
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
  ref_var <- group_by(ref, Target)
  ref_var <- summarise(ref_var,
                       min = min(mean),
                       max = max(mean),
                       sd = stats::sd(mean, na.rm = TRUE),
                       var = stats::var(mean, na.rm = TRUE))
  return(as.data.frame(ref_var))
}






deltaCt <- function(data,
                    imput.by.mean.within.groups = FALSE,
                    ref, save.to.txt = FALSE,
                    name.txt = "dCt_results"){
  data <- data %>%
    group_by(Group, Target, Sample) %>%
    summarise(mean = base::mean(Ct, na.rm = TRUE)) %>%
    as.data.frame()
  if (imput.by.mean.within.groups == TRUE){
  data_wide <- data %>%
    select(Group, Sample, Target, mean) %>%
    pivot_wider(names_from = Target, values_from = mean)
  data_wide_imp <- data_wide %>%
      group_by(Group) %>%
      mutate_if(is.numeric, ~ replace(., is.na(.), mean(., na.rm = TRUE)))
  nas <- sum(is.na(data_wide))
  percentage <- sum(is.na(data_wide))/((ncol(data_wide)-2)*nrow(data_wide))
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

    data_wide <- data %>%
      select(Group, Sample, Target, mean) %>%
      pivot_wider(names_from = Target, values_from = mean)
    nas <- sum(is.na(data_wide))
    percentage <- sum(is.na(data_wide))/((ncol(data_wide)-2)*nrow(data_wide))
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






control_dCt_cluster <- function(data, method.dist = "euclidean", method.clust = "complete",
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





control_dCt_pca <- function(data, point.size = 4, alpha = 0.7, label.size = 3, hjust = 0, vjust = -1,
                            col = c("#66c2a5", "#fc8d62"),
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
    geom_point(size = point.size, alpha = alpha) +
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




control_dCt_corr <- function(data, add.coef = "black",
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
                            label = c("eq", "R2", "p"),
                            dpi = 600, width = 15, height = 15,
                            save.to.tiff = FALSE,
                            name.tiff = "dCt_corr_single_plot"){
  if (by.group == TRUE){
    corr_control <- ggplot(data, aes(x = .data[[x]], y = .data[[y]], color = Group)) +
      geom_point() +
      geom_smooth(method='lm', se = FALSE, alpha = alpha) +
      stat_poly_eq(use_label(c("eq", "R2", "p")),
                   label.y = c(label.position.y),
                   label.x = c(label.position.x)) +
      scale_color_manual(values = c(col)) +
      xlab(x) + ylab(y) +
      labs(color = legend.title, title = plot.title) +
      theme_classic() +
      theme(legend.position = legend.position) +
      theme(axis.text = element_text(size = axis.text.size, colour = "black")) +
      theme(axis.title = element_text(size = axis.title.size, colour="black")) +
      theme(legend.text = element_text(size = legend.text.size, colour="black")) +
      theme(legend.title = element_text(size = legend.title.size, colour="black"))
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





results_dCt_boxplot <- function(data, coef = 1.5, sel.Target = "all", by.group = TRUE,
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

  }
  if (angle != 0){
    box_results <- box_results +
      guides(x =  guide_axis(angle = angle))
  }
  if (rotate == TRUE){
    box_results <- box_results +
      coord_flip()
  }
  if (add.mean == TRUE){
    box_results <- box_results +
      stat_summary(aes(group = Group),
                   fun = mean,
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









ddCt <- function(data,
                 group.study,
                 group.ref,
                 do.tests = TRUE,
                 save.to.txt = FALSE,
                 name.txt = "ddCt_RQ_results"){

  data_slim <- data %>%
    filter(Group == group.study | Group == group.ref) %>%
    pivot_longer(cols = -c(Group, Sample), names_to = "Target", values_to = "dCt")
  data_ddCt <- data_slim %>%
    group_by(Group, Target) %>%
    summarise(ddCt = mean(dCt, na.rm = TRUE)) %>%
    pivot_wider(names_from = Group, values_from = ddCt) %>%
    mutate(ddCt = .data[[group.study]] - .data[[group.ref]]) %>%
    mutate(RQ = 2^-ddCt) %>%
    rename_with(~paste0(.x, "_mean", recycle0 = TRUE), all_of(c(group.study, group.ref)))
  if (do.tests == TRUE){
  data_ddCt_norm <- data_slim %>%
    group_by(Group, Target) %>%
    summarise(shap_wilka_p = shapiro.test(dCt)$p.value) %>%
    pivot_wider(names_from = Group, values_from = shap_wilka_p) %>%
    full_join(data_ddCt, by = c("Target")) %>%
    mutate(test.for.comparison = ifelse(.data[[group.ref]] >= 0.05 & .data[[group.study]] >= 0.05, "t.student's.test", "Mann-Whitney.test")) %>%
    rename_with(~paste0(.x, "_norm_p", recycle0 = TRUE), all_of(c(group.study, group.ref)))
  data_ddCt_tests <- data_slim %>%
    group_by(Target) %>%
    summarise(t.test.p = t.test(dCt ~ Group, alternative = "two.sided")$p.value,
              t.test.stat = t.test(dCt ~ Group, alternative = "two.sided")$statistic,
              MW.test.p = wilcox.test(dCt ~ Group, alternative = "two.sided")$p.value,
              MW.test.stat = wilcox.test(dCt ~ Group, alternative = "two.sided")$statistic)
  data_ddCt_norm_tests <- full_join(data_ddCt_norm, data_ddCt_tests, by = c("Target"))
  return(data_ddCt_norm_tests)
  } else{
    return(data_ddCt)
  }
  if (save.to.txt == TRUE){
    write.table(data_ddCt_norm_tests, paste(name.txt,".txt", sep = ""))
  }
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
    data <- mutate(data, p.used = ifelse(test.for.comparison == "t.student's.test", yes = t.test.p,  no = MW.test.p))
  }
  if (mode == "user"){
    colnames(user) <- c("Target","test")
    data <- full_join(data, user, by = c("Target"))
    data <- mutate(data, p.used = ifelse(test == "t.student's.test", yes = t.test.p,  no = MW.test.p))
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
                  smooth = FALSE, auc = TRUE, plot=FALSE, ci=TRUE, of = "auc")
    parameters <- coords(myproc, "best", ret = c("threshold", "specificity", "sensitivity","accuracy", "ppv", "npv", "youden"))
    roc_param[x,2:8] <- parameters
    roc_param[x,9] <- myproc$auc
    roc_param[x,10] <- colnames(data)[x+2]
  }
  if (save.to.tiff == TRUE){
    tiff(paste(name.tiff,".tiff", sep = ""), res = dpi, width = width, height = height, units = "cm", compression = "lzw")
    par(mfrow = c(panels.row, panels.col))
    for (x in 1:nrow(roc_param)){
      myproc <- roc(response = data$Group, predictor = as.data.frame(data)[ ,x+2], levels = c(Groups),
                    smooth = FALSE, auc = TRUE, plot=FALSE, ci=TRUE, of = "auc")
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

