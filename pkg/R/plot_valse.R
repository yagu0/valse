#' Plot
#'
#' It is a function which plots relevant parameters
#'
#' @param X matrix of covariates (of size n*p)
#' @param Y matrix of responses (of size n*m)
#' @param model the model constructed by valse procedure
#' @param comp TRUE to enable pairwise clusters comparison
#' @param k1 index of the first cluster to be compared
#' @param k2 index of the second cluster to be compared
#'
#' @importFrom ggplot2 ggplot aes ggtitle geom_tile geom_line geom_point scale_fill_gradient2 geom_boxplot theme
#' @importFrom cowplot background_grid
#' @importFrom reshape2 melt
#'
#' @export
plot_valse <- function(X, Y, model, comp = FALSE, k1 = NA, k2 = NA)
{
  n <- nrow(X)
  K <- length(model$pi)
  ## regression matrices
  gReg <- list()
  for (r in 1:K)
  {
    Melt <- melt(t((model$phi[, , r])))
    gReg[[r]] <- ggplot(data = Melt, aes(x = "Var1", y = "Var2", fill = "value")) +
      geom_tile() + scale_fill_gradient2(low = "blue", high = "red", mid = "white",
      midpoint = 0, space = "Lab") + ggtitle(paste("Regression matrices in cluster", r))
  }
  print(gReg)

  ## Differences between two clusters
  if (comp)
  {
    if (is.na(k1) || is.na(k2))
      print("k1 and k2 must be integers, representing the clusters you want to compare")
    Melt <- melt(t(model$phi[, , k1] - model$phi[, , k2]))
    gDiff <- ggplot(data = Melt, aes(x = "Var1", y = "Var2", fill = "value"))
      + geom_tile()
      + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
        space = "Lab")
      + ggtitle(paste("Difference between regression matrices in cluster",
        k1, "and", k2))
    print(gDiff)
  }

  ### Covariance matrices
  matCov <- matrix(NA, nrow = dim(model$rho[, , 1])[1], ncol = K)
  for (r in 1:K)
    matCov[, r] <- diag(model$rho[, , r])
  MeltCov <- melt(matCov)
  gCov <- ggplot(data = MeltCov, aes(x = "Var1", y = "Var2", fill = "value")) + geom_tile()
    + scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
      space = "Lab")
    + ggtitle("Covariance matrices")
  print(gCov)

  ### Proportions
  gam2 <- matrix(NA, ncol = K, nrow = n)
  for (i in 1:n)
    gam2[i, ] <- c(model$proba[i, model$affec[i]], model$affec[i])

  bp <- ggplot(data.frame(gam2), aes(x = "X2", y = "X1", color = "X2", group = "X2"))
    + geom_boxplot()
    + theme(legend.position = "none")
    + background_grid(major = "xy", minor = "none")
  print(bp)

  ### Mean in each cluster
  XY <- cbind(X, Y)
  XY_class <- list()
  meanPerClass <- matrix(0, ncol = K, nrow = dim(XY)[2])
  for (r in 1:K)
  {
    XY_class[[r]] <- XY[model$affec == r, ]
    if (sum(model$affec == r) == 1) {
      meanPerClass[, r] <- XY_class[[r]]
    } else {
      meanPerClass[, r] <- apply(XY_class[[r]], 2, mean)
    }
  }
  data <- data.frame(mean = as.vector(meanPerClass),
    cluster = as.character(rep(1:K, each = dim(XY)[2])), time = rep(1:dim(XY)[2], K))
  g <- ggplot(data, aes(x = "time", y = "mean", group = "cluster", color = "cluster"))
  print(g + geom_line(aes(linetype = "cluster", color = "cluster"))
    + geom_point(aes(color = "cluster")) + ggtitle("Mean per cluster"))
}
