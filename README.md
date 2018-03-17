
### Installation

**VIPER** relies on the following R packages: **Rcpp**, **RcppArmadillo**, **quadprog**, **glmnet**. All packagess are hosted on CRAN. 
  ```R
  install.packages("Rcpp")
  install.packages("RcppArmadillo")
  install.packages("quadprog")
  install.packages("glmnet")
  ```

**VIPER** can be installed from github directly as follows:

  ```R
  install.packages("devtools")
  library(devtools)
  install_github("ChenMengjie/VIPER")
  ```
  

### Quick Start

**VIPER()** is the main funtion that implements the imputation methods. Details about the steps can be found in later sections.
To perform imputation, simply run:
  ```R
  library(VIPER)
  data(grun)
  res <- VIPER(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE, alpha = 1, 
  report = FALSE, outdir = NULL, prefix = NULL)
  ```

Our method is not restricted to the units of measurement and is applicable to all normalized measurements such as RPM (reads per million reads), TPM (transcripts per kilobase per millions reads) or RPKM (reads per kilobase per millions reads). 

Paramter **num** is the number of random sampled genes used to fit the penalized regression model to identify the set of candidate cells. The default value is 5000. To reduce the influence of missing values in the weight estimation, the nonnegative regression model is fitted using genes with a zero rate less than a certain threshold. The threshold is set using **percentage.cutoff**. 

VIPER calls **cv.glmnet()** in **glmnet** to perform fitting cross validation. Two penalty levels are available for selection of penalty level **lambda**: **lambda.min**, the value of **lambda** that gives minimum mean cross-validated error, and **lambda.1se**, the value of **lambda** that gives the most regularized model such that error is within one standard error of the minimum. The default is **lambda.1se**, i.e., **minbool = FALSE**. 

Paramter **alpha** sets the elastic net mixing parameter. The default value is 1, which is equivalent to a lasso model. We used **alpha = 0.5** for elastic net. 

Paramter **report** controls whether to save imputed data matrix in CSV files. The default value is **FALSE**.
When **report** is **TRUE**, paramter **outdir** provide the directory to save the result files.
Paramter **prefix** sets the prefix of the result files.

**VIPER()** returns the following list of values:

**imputed_log**: A **p** by **n** matrix of log transformed gene expression levels after imputation.

**imputed**: A **p** by **n** matrix of gene expression counts converted from imputed log transformed values.

**sample_weights**: A **n** by **n** matrix of estimated imputation weights. Each row represents a cell.

**outliers**: The indexes of cells that have no selected candidate neighbors according to the penalized regression model. The zero values in these cells are not imputed.




### Contact

**Mengjie Chen** (UChicago) mengjiechen@uchicago.edu
