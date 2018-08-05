
predict.l0araxx <- function(obj, newx, type=c("link", "response", "coefficients", "class"), offset=NULL, ...) {
  np = dim(newx)
  type <- match.arg(type)

  # check offset
  if (is.null(offset)) {
    offset <- rep(0, np[1])
  } else if (length(offset) != np[1]) {
    stop("length of offset not equal to the length of newx")
  }

  # getting coefficients
  beta <- coef.l0araxx(obj)
  if (type=="coefficients") return(beta)

  # getting newx
  if(missing(newx)) newx <- obj$x
  if (obj$standardize) newx <- scale(newx)
  newx <- cbind(rep(1, np[1]), newx)

  # calculates link value
  eta <- newx %*% beta + offset
  if (type=="link") return(drop(eta))

  # calculates response var
  response <- switch(obj$family,
                     gaussian = eta,
                     logit = exp(eta)/(1+exp(eta)),
                     poisson = exp(eta),
                     gamma = 1/eta,
                     inv.gaussian = 1/sqrt(eta))
  if(obj$family == "gaussian" & obj$standardize) response <- response + mean(obj$y)
  if (type=="response") return(drop(response))

  # calculates binary class
  if (type=="class") {
    if (object$family=="logit") {
      return(drop(1*(eta>0)))
    } else {
      stop("type='class' can only be used with family='logit'")
    }
  }
}

coef.l0araxx <- function(obj, ...) {
  coefs <- as.vector(obj$beta)
  p <- length(coefs)
  names(coefs)[1] <- "Intercept"
  names(coefs)[2:p] <- paste0("X",1:(length(coefs)-1))
  return(coefs)
}

print.l0ara <- function(x, ...) {
  family <- switch(x$family, gaussian = "Linear regression", logit = "Logistic regression", poisson = "Poisson regression", inv.gaussian = "Inverse gaussian regression", gamma = "Gamma regression")
  cat("Lambda used : ", x$lambda, "\n")
  cat("Model : ", family, "\n")
  cat("Iterations : ", x$iter, "\n")
  cat("Degree of freedom : ", x$df,"\n")
}
