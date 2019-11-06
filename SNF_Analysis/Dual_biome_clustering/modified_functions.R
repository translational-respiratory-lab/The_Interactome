estimateNumberOfClustersGivenGraph<-function (W, NUMC = 2:5) 
{
  if (min(NUMC) == 1) {
    warning("Note that we always assume there are more than one cluster.")
    NUMC <- NUMC[NUMC > 1]
  }
  W <- (W + t(W))/2
  diag(W) <- 0
  if (length(NUMC) <= 0) {
    warning(paste("Invalid NUMC provided, must be an integer vector", 
                  "with atleast one other number than 1.", "Using default NUMC=c(2,3,4,5)", 
                  sep = ""))
    NUMC <- 2:5
  }
  degs <- rowSums(W)
  degs[degs == 0] <- .Machine$double.eps
  D <- diag(degs)
  L <- D - W
  Di <- diag(1/sqrt(degs))
  L <- Di %*% L %*% Di
  eigs <- eigen(L)
  eigs_order <- sort(eigs$values, index.return = T)$ix
  eigs$values <- eigs$values[eigs_order]
  eigs$vectors <- eigs$vectors[, eigs_order]
  eigengap <- abs(diff(eigs$values))
  quality <- list()
  for (c_index in 1:length(NUMC)) {
    ck <- NUMC[c_index]
    UU <- eigs$vectors[, 1:ck]
    EigenvectorsDiscrete <- discretisation(UU)[[1]]
    EigenVectors <- EigenvectorsDiscrete^2
    temp1 <- EigenVectors[do.call(order, lapply(1:ncol(EigenVectors), 
                                                function(i) EigenVectors[, i])), ]
    temp1 <- t(apply(temp1, 1, sort, TRUE))
    quality[[c_index]] <- (1 - eigs$values[ck + 1])/(1 - 
                                                       eigs$values[ck]) * sum(sum(diag(1/(temp1[, 1] + .Machine$double.eps)) %*% 
                                                                                    temp1[, 1:max(2, ck - 1)]))
  }
  t1 <- sort(eigengap[NUMC], decreasing = TRUE, index.return = T)$ix
  K1 <- NUMC[t1[1]]
  K12 <- NUMC[t1[2]]
  t2 <- sort(unlist(quality), index.return = TRUE)$ix
  K2 <- NUMC[t2[1]]
  K22 <- NUMC[t2[2]]
  output <- list(`Eigen-gap best` = K1, `Eigen-gap 2nd best` = K12, 
                 `Rotation cost best` = K2, `Rotation cost 2nd best` = K22)
  return(output)
}


discretisation<-function (eigenVectors) 
{
  normalize <- function(x) x/sqrt(sum(x^2))
  eigenVectors = t(apply(eigenVectors, 1, normalize))
  eigenVectors[is.na(eigenVectors)]<-0
  eigenVectors[is.infinite(eigenVectors)]<-0
  n = nrow(eigenVectors)
  k = ncol(eigenVectors)
  R = matrix(0, k, k)
  R[, 1] = t(eigenVectors[round(n/2), ])
  mini <- function(x) {
    i = which(x == min(x))
    return(i[1])
  }
  c = matrix(0, n, 1)
  for (j in 2:k) {
    c = c + abs(eigenVectors %*% matrix(R[, j - 1], k, 1))
    i = mini(c)
    R[, j] = t(eigenVectors[i, ])
  }
  lastObjectiveValue = 0
  for (i in 1:20) {
    eigenDiscrete = SNFtool:::.discretisationEigenVectorData(eigenVectors %*% 
                                                     R)
    svde = svd(t(eigenDiscrete) %*% eigenVectors)
    U = svde[["u"]]
    V = svde[["v"]]
    S = svde[["d"]]
    NcutValue = 2 * (n - sum(S))
    if (abs(NcutValue - lastObjectiveValue) < .Machine$double.eps) 
      break
    lastObjectiveValue = NcutValue
    R = V %*% t(U)
  }
  return(list(discrete = eigenDiscrete, continuous = eigenVectors))
}