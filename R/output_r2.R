fitting_lasso_return_r2 <- function(y, X, alpha = 1){

  require(glmnet)
  cv.lasso <- cv.glmnet(X, y, intercept = FALSE, alpha = alpha)
  r2 <- cv.lasso$glmnet.fit$dev.ratio[which(cv.lasso$glmnet.fit$lambda == cv.lasso$lambda.min)]
  coeff <- as.vector(coef(cv.lasso, s = cv.lasso$lambda.min))
  selected <- which(coeff!=0)
  res <- list(coeff = coeff[selected], selected = selected - 1, r2 = r2)
  return(res)
}

PredictCell <- function(GeneExpression, GeneNum = 5000, CellNum = NULL, ZeroRate = 0.1, SaveSelection = FALSE){

  xx <- GeneExpression # p*n
  p <- nrow(xx)
  n <- ncol(xx)

  logxx <- apply(xx, 2, function(y){log(y + 0.1)})
  data <- logxx[round(runif(GeneNum)*p), ]

  zero.rate <- apply(xx, 1, function(x){length(x[x == 0])})/n
  flag <- zero.rate <= ZeroRate

  outlier.list <- NULL
  r_square_table <- list(NULL)

  if(is.null(CellNum)){

    cell.index <- 1:n

  } else {

    cell.index <- sample(1:n, CellNum)

  }

  for(j in cell.index){

    remain <- data[, -j]
    res <- fitting_lasso_return_r2(data[, j], remain)
    selected <- res$selected

    if(length(selected) < 3){

      outlier.list <- c(outlier.list, j)

    } else {

      kk <- summary(lm(logxx[flag, j]~logxx[flag, -j][, selected]))

      if(SaveSelection == TRUE){

        r_square_table[[j]] <- list(r2lasso = res$r2, r2 = kk$r.squared, r2adj = kk$adj.r.squared, num = length(selected), name = colnames(xx)[j], selected = selected)

      } else {

        r_square_table[[j]] <- list(r2lasso = res$r2, r2 = kk$r.squared, r2adj = kk$adj.r.squared, num = length(selected), name = colnames(xx)[j])

      }

    }

  }

  res <- list(outlierList = outlier.list, R2table = r_square_table)
  return(res)
}


PredictGene <- function(GeneExpression, GeneNum = 1000, CellNum = NULL, ZeroRate = 0.1, SaveSelection = FALSE){

  xx <- GeneExpression # p*n
  p <- nrow(xx)
  n <- ncol(xx)

  logxx <- apply(xx, 2, function(y){log(y + 0.1)})

  zero.rate <- apply(xx, 1, function(x){length(x[x == 0])})/n
  flag <- zero.rate <= ZeroRate
  logxx <- logxx[flag, ]
  remainID <- c(1:p)[flag]

  if(is.null(CellNum)){

    data <- t(logxx)

  } else {

    CellNum <- min(n, CellNum)
    data <- t(logxx)[sample(1:n, CellNum), ]

  }

  outlier.list <- NULL
  r_square_table <- list(NULL)

  newNum <- length(remainID)
  if(GeneNum < newNum){

    selected.nums <- sample(1:newNum, GeneNum)

  } else {

    selected.nums <- 1:newNum

  }

  finalGeneNum <- length(selected.nums)

  for(j in 1:finalGeneNum){

    id <- selected.nums[j]
    remain <- data[, -j]
    res <- try(fitting_lasso_return_r2(data[, j], remain))

    if(class(res) != "try-error"){
      coeff <- res$coeff
      selected <- res$selected

      if(length(selected) >= 1){

        kk <- summary(lm(data[, j]~remain[, selected]))

        if(SaveSelection == TRUE){

          r_square_table[[j]] <- list(r2lasso = res$r2, r2 = kk$r.squared, r2adj = kk$adj.r.squared, num = length(selected), name = rownames(xx)[remainID[id]], selected = remainID[selected])

        } else {

          r_square_table[[j]] <- list(r2lasso = res$r2, r2 = kk$r.squared, r2adj = kk$adj.r.squared, num = length(selected), name = rownames(xx)[remainID[id]])

        }

      } else {

        outlier.list <- c(outlier.list, remainID[id])

      }

    } else {

      outlier.list <- c(outlier.list, remainID[j])

    }

  }

  res <- list(outlierList = outlier.list, R2table = r_square_table)
  return(res)
}

PlotR2 <- function(r2CellSummary, r2GeneSummary = NULL, DatasetName = NULL, r2type = c("r2adj", "r2", "r2lasso")){

  r2type <- match.arg(r2type)

  require(ggplot2)
  require(easyGgplot2)

  x1 <- r2CellSummary$R2table
  tt1 <- paste0("r2_1 <- sapply(x1, function(z){z$", r2type, "})")
  eval(parse(text=tt1))
  r2_1 <- unlist(r2_1)

  if(!is.null(r2GeneSummary)){

    x2 <- r2GeneSummary$R2table
    tt2 <- paste0("r2_2 <- sapply(x2, function(z){z$", r2type, "})")
    eval(parse(text=tt2))
    r2_2 <- unlist(r2_2)
    df1 <- data.frame(Type = c(rep("By Cell", length(r2_1)), rep("By Gene", length(r2_2))), R2 = c(r2_1, r2_2))

  } else {

    df1 <- data.frame(Type = rep("By Cell", length(r2_1)), R2 = r2_1)

  }

   ff <- ggplot2.density(data=df1, xName='R2', groupName='Type', legendPosition="top", addMeanLine=FALSE,
   size = 3, backgroundColor="white",  groupColors = c('#66B2FF', '#FFD4AA'),
   ytitle="Density", xtitle=" R squared ", mainTitle = DatasetName, alpha=0.5, fillGroupDensity=TRUE, removePanelGrid=TRUE,
   removePanelBorder=TRUE, showLegend=TRUE, legendTitle = element_blank(), legendTitleFont = c(15, "bold", "black"),
   legendTextFont = c(15, "bold", "black"), xtitleFont = c(15, "bold", "black"),  ytitleFont=c(15, "bold", "black"),
   xTickLabelFont=c(15, "bold", "black"), yTickLabelFont=c(15, "bold", "black"))

   ff
}
