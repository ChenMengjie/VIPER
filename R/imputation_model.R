#' @title The main function that performs imputation for single cell RNA sequencing data.
#' @description VIPER applies a weighted penalized regression model to actively select a sparse set of local neighborhood cells
#' that are most predictive of the cell of interest. The selection of this sparse set of cells is done in a progressive manner.
#' First, for each given cell, a penalized regression model is applied using a random sampled set of genes to identify a set of
#' candidate cells that are predictive of the expression of the cell of interest. Then we use a nonnegative regression model to
#' further refine the list of neighborhood cells and estimate their imputation weights. In addition, VIPER explicitly accounts
#' for expression measurement uncertainty of the zero values in scRNAseq by modeling the dropout probability in a cell-specific
#' and gene-specific fashion.
#' @param gene.expression: A p by n matrix of gene expression levels for p genes from n samples.
#' Our method is not restricted to the units of measurement and is applicable to all normalized measurements
#' such as RPM (reads per million reads), TPM (transcripts per kilobase per millions reads) or RPKM (reads per kilobase per millions reads).
#' @param num: The number of random sampled genes used to fit the penalized regression model to identify the set of candidate cells.
#' The default value is 5000. If gene number p in the dataset is less than specified \code{num}, \code{num}will be set as 0.8*p.
#' @param percentage.cutoff: To reduce the influence of missing values in the weight estimation, the nonnegative regression model is fitted using genes
#' with a zero rate less than a certain threshold. The default value is 0.1 (10 percent).
#' @param minbool: The criteria used to select the penalty level\code{lambda} in the penalized regression model. VIPER calls \code{cv.glmnet()} in \bold{glmnet} to perform fitting cross validation.
#' Two penalty levels are available for selection: \code{lambda.min}, the value of \code{lambda} that gives minimum mean cross-validated error, and \code{lambda.1se}, the value of \code{lambda}
#' that gives the most regularized model such that error is within one standard error of the minimum. The default is \code{lambda.1se}, i.e., \code{minbool = FALSE}.
#' @param alpha: The elastic net mixing parameter. The default value is 1, which is equivalent to a lasso model.
#' @param report: Whether to save imputed data matrix in CSV files. The default value is FALSE.
#' @param outdir: The directory to save the output.
#' @param prefix: prefix of the result files.
#' @return A list of imputed data matrices and summary.
#' \item{imputed_log}{A p by n matrix of log transformed gene expression levels after imputation.}
#' \item{imputed}{A p by n matrix of gene expression levels after imputation converted from log transformed values.}
#' \item{sample_weights}{ A n by n matrix of estimated imputation weights. Each row represents a cell.}
#' \item{outliers}{The indexes of cells that have no selected candidate neighbors according to the penalized regression model.}
#' The zero values in these cells are not imputed.

VIPER <- function(gene.expression, num = 5000, percentage.cutoff = 0.1, minbool = FALSE,
                  alpha = 0.5, report = FALSE, outdir = NULL, prefix = NULL){

  xx <- gene.expression # p*n
  p <- nrow(xx)
  n <- ncol(xx)
  if(p < num) num <- round(0.8*p)
  zero.rate <- apply(xx, 1, function(x){length(x[x == 0])})/n
  flag <-  zero.rate <= percentage.cutoff

  logxx <- apply(xx, 2, function(y){log(y + 0.1)})
  data <- logxx[round(runif(num)*p), ]
  zero.matrix <- xx != 0
  zero.matrix <- apply(zero.matrix, 2, as.numeric)
  selected_logxx <- logxx[flag, ]

#  if(expectation == TRUE){
#  res_imp <- imputation_by_samples_expectation(data, selected_logxx, logxx, zero.matrix, n, p, minbool, alpha)
#  } else {
    res_imp <- imputation_by_samples(data, selected_logxx, logxx, zero.matrix, n, p, minbool, alpha)
#  }

  outlier_flag <- apply(res_imp$sample_weights, 1, function(x){any(x==-1)})
  outliers <- c(1:n)[outlier_flag]

  nopredict <- logxx
  nopredict[gene.expression==0] <- res_imp$imputed[gene.expression==0]

  imputed_counts <- round(exp(nopredict)-0.1)

  colnames(imputed_counts) <- colnames(nopredict) <- colnames(xx)
  rownames(imputed_counts) <- rownames(nopredict) <- rownames(xx)

  res <- list(imputed = imputed_counts, imputed_log = nopredict,
              sample_weights = res_imp$sample_weights, outliers = outliers)
  write.table(imputed_counts, file = paste0(outdir, "/", prefix, "imputed_counts.csv"))
  write.table(res$imputed_log, file = paste0(outdir, "/", prefix, "imputed_logvalue.csv"))
  write.table(res$sample_weights, file = paste0(outdir, "/", prefix, "sample_weights.csv"))
  write.table(res$outliers, file = paste0(outdir, "/", prefix, "outliers.csv"))

  return(res)
}

