
### RQdeltaCT - an R package for relative quantification of gene expression using delta Ct methods 

`RQdeltaCT` is an R package developed to perform relative quantification of gene expression using delta Ct family methods (encompassing 2^-dCt, and 2^-ddCt method), originally proposed by Kenneth J. Livak and Thomas D. Schmittgen in [Article1](https://www.sciencedirect.com/science/article/pii/S1046202301912629?via%3Dihub) and [Article2](https://www.nature.com/articles/nprot.2008.73).  

These methods have been designed to analyse gene expression data (Ct values) obtained from real-time PCR experiments. The main idea is to:

* normalise gene expression values using endogenous control gene,
* present gene expression levels in linear form by using the 2^-(value) transformation,
* calculate differences in gene expression levels between groups of samples (or technical replicates of a single sample).  

The `RQdeltaCT` package offers functions that encompass all of these steps, together with:
1. importing qPCR datasets, 
2. performing multi-step quality control of data,
3. enabling numerous data visualisations,
4. enrichment of standard workflow with additional useful methods including correlation analysis, Receiver Operating Characteristic analysis, and logistic regression),
5. a convenient export of obtained results in table and image forms.  

To install and load the `RQdeltaCT` package, simply run in R or RStudio:  

`remotes::install_github("Donadelnal/RQdeltaCT")`

`library(RQdeltaCT)`

**The package has been designed to be friendly to non-experts in R programming users. No additional, extensive coding steps are necessary in the standard workflow. Detailed demonstration of the package functionalities using test data can be found in the prepared [vignette](https://donadelnal.r-universe.dev/articles/RQdeltaCT/my-vignette.html).**  
