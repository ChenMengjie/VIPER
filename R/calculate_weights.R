calculate_weights <- function(z, X){

  suppressMessages(require(quadprog))

	xx <- X[, 1]
	remain <- as.matrix(X[, -1])
	k <- ncol(remain)
	n <- length(xx)

	w <- z - xx
	Y <- apply(remain, 2, function(ll){ ll - xx })


	D <- matrix(0, ncol = k, nrow = k)
	for(i in 1:n){
		D <- D + Y[i, ]%*%t(Y[i, ])
	}

	d <- rep(0, k)
	for(i in 1:n){
		d <- d + w[i]*t(Y[i, ])
	}

	A <- rbind(diag(k), - diag(k))
	b <- c(rep(1, k), rep(0, k))

	res <- try(solve.QP(D, d, -t(A), -b)$solution, silent = TRUE)
  if(class(res) == "try-error"){
    res <- try(solve.QP(nearPD(D, ensureSymmetry = TRUE)$mat, d, -t(A), -b)$solution, silent = TRUE)
    if(class(res) == "try-error"){
      coef <- rep(1/(k+1), k+1)
    } else {
      coef <- c(1-sum(res), res)
    }
  } else {
    coef <- c(1-sum(res), res)
  }

	return(coef)
}


fitting_lasso <- function(y, X, type = "min", alpha = 1){

  suppressMessages(require(glmnet))
	cv.lasso <- cv.glmnet(X, y, intercept = FALSE, alpha = alpha)
	if(type == "min"){
	  coeff <- as.vector(coef(cv.lasso, s = cv.lasso$lambda.min))
	} else {
    coeff <- as.vector(coef(cv.lasso, s = cv.lasso$lambda.1se))
  }
  selected <- which(coeff!=0)
  res <- list(coeff = coeff[selected], selected = selected - 1)
	return(res)
}
