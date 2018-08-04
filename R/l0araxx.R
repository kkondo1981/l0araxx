
l0araxx <- function(x, y, weights = NULL, offset = NULL,
                    family = c("gaussian", "logit", "gamma", "poisson", "inv.gaussian"),
                    lambda, standardize = TRUE, maxit = 10^3, eps = 1e-04)
{
  x0 <- x
  y0 <- y

  # error checks
  if (class(x) != "matrix") {
    tmp <- try(x <- model.matrix(~0 + ., data = x), silent = TRUE)
    if (class(tmp)[1] == "try-error")
      stop("x must be a matrix or able to be coerced to a matrix")
  }
  np = dim(x)
  if (is.null(np) | (np[2] <= 1)){
    stop("x should be a matrix with 2 or more columns")
  }

  if (class(y) != "numeric") {
    tmp <- try(y <- as.numeric(y), silent = TRUE)
    if (class(tmp)[1] == "try-error")
      stop("y must numeric or able to be coerced to numeric")
  }
  y <- drop(y)
  if (length(y) != np[1]){
    stop("length of y not equal to the number of rows of x")
  }

  if(missing(lambda)){
    stop("'lambda' was not found(must be a postiive number)")
  }
  if(length(lambda) > 1) {
    stop("Require lenght >=1 for 'lambda'")
  }
  if(lambda < 0) {
    warning("lambda < 0; set to 0")
    lambda <- 0
  }

  if (is.null(weights)) {
    weights <- rep(1, np[1])
  } else if (length(weights) != np[1]) {
    stop("length of weights not equal to the length of y")
  }

  if (is.null(offset)) {
    offset <- rep(0, np[1])
  } else if (length(offset) != np[1]) {
    stop("length of offset not equal to the length of y")
  }

  family <- match.arg(family)

  # standardize the design matrix 'x' if needed, and enhance it
  if (standardize) x <- scale(x)
  x <- cbind(rep(1, np[1]), x)

  # In case of gaussian, centering the outcome vector 'y' 
  if(family == "gaussian" & standardize) y <- y - mean(y)

  # do estimation and return
  out <- l0araxxC(x, y, weights, offset, family, lambda, maxit, eps)

  # create obj and return it
  obj <- list(beta = drop(out$beta),
              beta_hist = out$beta_hist,
              df = sum(out$beta != 0),
              lambda = lambda,
              iter = out$iter,
              family = family,
              standardize = standardize,
              x = x0,
              y = y0)
  class(obj) <- "l0araxx"
  return(obj)
}