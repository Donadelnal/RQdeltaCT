#' @title Gene expression dataset for RQdeltaCT package. Siutable for pairwise analysis.
#' @description Dataset that contain table with gene expression data obtained for 18 genes and 42 samples using qPCR experiments with TaqMan assays.
#' Samples are divided into two groups: Before (21 samples) and After (21 samples). This dataset is created based on the real data used in article:
#' Zalewski, D.; Chmiel, P.; Kołodziej, P.; Borowski, G.; Feldo, M.; Kocki, J.; Bogucka-Kocka, A. Dysregulations of Key Regulators of Angiogenesis and Inflammation in Abdominal Aortic Aneurysm. Int. J. Mol. Sci. 2023, 24, 12087. https://doi.org/10.3390/ijms241512087
#' @format A data frame with 756 rows and 4 variables:
#' \describe{
#'   \item{\code{Sample}}{character COLUMN_DESCRIPTION}
#'   \item{\code{Gene}}{character COLUMN_DESCRIPTION}
#'   \item{\code{Ct}}{character COLUMN_DESCRIPTION}
#'   \item{\code{Group}}{character COLUMN_DESCRIPTION}
#'}
#' @source Dataset is attached to the `RQdeltaCT` package and can be loaded using `data(data.Ct.pairwise)` command.
"data.Ct.pairwise"
