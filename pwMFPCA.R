#################################
# Multivariate Functional PCA  #
# Berrendero, Justel, and Svarc #
# Functional Data Analysis      #
#################################

# This function returns a list whose elements are:

# 1: A roahd object with the new curves represented by ncomp components
# 2: The cumulative variability by components according to criterion pi_1
# 3: The cumulative variability by components according to criterion pi_3
# 4: A list containing the eigenvalues and eigenvectors of the Sigma(t) matrix

library(roahd)

# It is assumed that the discretized points are equidistant

# MFPCA Functional Principal ----

pwMFPCA = function(mfData, ncomp) {
  
  if (class(mfData) != "mfData") 
    stop("Parameter 'mfData' must be passed as a roahd::mfData object.")
  
  n = mfData$N
  p = mfData$L
  Taus = mfData$P
  Domain = seq(mfData$t0, mfData$tP, length.out = Taus)
  
  # Norm function ----
  
  norm = function(x) sqrt(sum(x^2))
  
  # Sign of eigenvectors ----
  
  Sign = function(Sigmat) {
    for (i in 1:p) {
      for (j in 2:Taus) {
        norms1 = NULL
        norms2 = NULL
        for (k in max(1, j - 5):(j - 1)) {
          norms1 = c(norms1, norm(Sigmat[[k]]$vectors[, i] - Sigmat[[j]]$vectors[, i]))
          norms2 = c(norms2, norm(Sigmat[[k]]$vectors[, i] + Sigmat[[j]]$vectors[, i]))
        }
        if (mean(norms1) > mean(norms2)) {
          Sigmat[[j]]$vectors[, i] = -Sigmat[[j]]$vectors[, i]
        }
      }
    }
    return(Sigmat)
  }
  
  # Eigenvalues and eigenvectors ----
  
  Xt = list()
  Sigmat = list()
  
  lambda = matrix(0, nrow = Taus, ncol = p)
  for (i in 1:Taus) {
    mat = matrix(0, nrow = n, ncol = p)
    for (j in 1:p) {
      mat[, j] = mfData$fDList[[j]]$values[, i]
    }
    Xt[[i]] = mat
    Sigmat[[i]] = eigen(cov(mat))
    lambda[i, ] = Sigmat[[i]]$values
  }
  
  # Explained variability ----
  
  vt = rowSums(lambda)
  
  pi1r = NULL
  pi2r = NULL
  
  for (i in 1:p) {
    pi1r = c(pi1r, (1 / Taus) * sum(lambda[, i] / vt))
    pi2r = c(pi2r, (sum(lambda[, i]) / sum(vt)))
  }
  
  pi1 = rbind(pi1r, cumsum(pi1r))
  c("Crit1", "Cumulative") -> rownames(pi1)
  pi2 = rbind(pi2r, cumsum(pi2r))
  c("Crit2", "Cumulative") -> rownames(pi2)
  
  # Components calculation ----
  
  Sigmat = Sign(Sigmat) # Choose signs
  z = list()
  mat1 = matrix(0, nrow = n, ncol = Taus)
  
  for (i in 1:ncomp) {
    for (j in 1:Taus) {
      for (k in 1:n) {
        mat1[k, j] = Sigmat[[j]]$vectors[, i] %*% Xt[[j]][k, ]
      }
    }
    z[[i]] = mat1
  }
  
  roahd::mfData(Domain, z) -> mfpca.mfD
  
  object = list(functions = mfpca.mfD,
                pi_1 = round(pi1, 3),
                pi_2 = round(pi2, 3),
                eigen = Sigmat)
  
  class(object) = "pwMFPCA"
  return(object)
}
