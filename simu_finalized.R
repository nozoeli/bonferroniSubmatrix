#$ -S /usr/local/bin/Rscript

#------------Approximation Net------------
bitIndex <- function(x, dig=floor(log(log(x, base=2), base=2))){
  # Returns the index set to be examed, by preassigned binary truncation
  x.bit <- which(intToBits(x) != 00)
  x.bitmax <- tail(x.bit, n=1)
  index.list <- list()
  count <- 1
  for (n.dig in 2:x.bitmax){
    if (n.dig == x.bitmax){
      if ((x.bitmax - 1) %in% x.bit){
        index.list[[count]] <- 2^(x.bitmax - 1) + 2^(x.bitmax - 2)
        count <- count + 1
      }
      index.list[[count]] <- 2^(x.bitmax - 1)
      count <- count + 1
    } else {
      index.list[[count]] <- 2^(n.dig - 1) 
      count <- count + 1
      index.list[[count]] <- 2^(n.dig - 1) + 2^(n.dig - 2)
      count <- count + 1
    }
  }
  return(unlist(index.list))
}

#------------Generate Random Normal Matrix------------
normalMatrix = function(M, N, m, n, mu){
  # M, N : dimensions of the matrix
  # m, n : dimensions of the submatrix
  # mu : positive mean
  # Creates a random matrix where the entries are independent normal with unit variance and mean 0, except for the top-left submatrix (of given size) where the mean is mu
  A = matrix(rnorm(M*N), M, N)
  A[1:m, 1:n] = A[1:m, 1:n] + mu
  return(A)
}

#------------Generate Random Poisson Matrix------------
poissonMatrix <- function(M, N, m, n, theta){
  # Generalize centered poisson distribution according to the standard exponential model
  # M, N : dimension of the matrix
  # m, n : dimension of the submatrix
  # mu stands for the mean of the elevated matrix. mu >= 1 required.
  # The base distribution is Poi(1) - 1
  mu <- exp(theta)
  mx1 <- matrix(1, M, N)
  if (m * n != 0){
    mx1[1:m, 1:n] <- matrix(mu, m, n)}
  mx <- matrix(rpois(M * N, as.vector(mx1)), M, N) - 1
  return(mx)
}

#------------Adaptive Testing------------
shabalin = function(A, m, n, cInd=NULL) {
  # A : real-valued matrix
  # m, n : dimensions of the submatrix
  # cInd : starting column indices (random by default)
  M = nrow(A)
  N = ncol(A)
  if (is.null(cInd)) {
    cSum = colSums(A)
    cInd = sort(cSum, d=TRUE, ind=TRUE)$ix[1:n] 
  }
  repeat {
    rSum = A[,cInd] %*% rep(1,n)
    rInd = sort(rSum, d=TRUE, ind=TRUE)$ix[1:m]
    cSum = rep(1,m) %*% A[rInd,]
    cIndNew = sort(cSum, d=TRUE, ind=TRUE)$ix[1:n]
    if (identical(cIndNew, cInd)) break 
    cInd = cIndNew
  }
  return(list('Rows' = rInd,'Cols' = cInd))
}

largestSbmxShabalin <- function(A, m, n, iter = m + n){
  # The function finds a m*n submatrix with largest entry sum
  # The method is Shabalin(2009) with random initial start to prevent local maximal
  # The function returns the sum of the submatrix
  M <- nrow(A)
  N <- ncol(A)
  base <- 0
  for (i in 1 : iter){
    r_temp <- shabalin(A, m, n, cInd = sample(x = 1 : N, n))  # random initial column index, preventing local maximal
    sum_temp <- sum(A[r_temp$Rows, r_temp$Cols])
    base <- base + (sum_temp - base) * (sum_temp > base)  # find the m*n matrix with largest entry sum
  }
  return(base)
}

#------------Unidimensional Permutation Test------------
uniPermuTest <- function(X, m, n, alpha = 0.05, iter = 500){
  # The function returns the test result of scan test calibrated by unidimensional permutation
  # The test is calibrated by permutation, with 'iter' iterations, default 500
  # X : the matrix being tested
  # m, n : the size of the anomaly in the alternative hypothesis
  # alpha : the test level
  M <- nrow(X)
  N <- ncol(X)
  originSum <- largestSbmxShabalin(X, m, n)
  pvalue <- 0
  for (j in 1 : iter){
    permSum <- largestSbmxShabalin(t(apply(X,1,sample)), m, n)
    pvalue <- pvalue + (permSum >= originSum) / (iter+1)
  }
  pvalue <- pvalue + 1 / (iter+1)
  permTest = (pvalue <= alpha) + 0
  return(list('pvalue' = pvalue, 'test' = permTest))
}

#------------Bidimensional Permutation Test------------
biPermuTest <- function(X, m, n, alpha = 0.05, iter = 500){
  # The function returns a result of the test of an elevated m*n submatrix
  # The test is calibrated by permutation, with 'iter' iterations, default 500
  # A : the matrix being tested
  # m, n : the size of the anomaly in the alternative hypothesis
  # alpha : the test level
  M <- nrow(X)
  N <- ncol(X)
  originSum <- largestSbmxShabalin(X, m, n)
  vectorA <- as.vector(X)
  pvalue <- 0
  for (j in 1 : iter){
    permSum <- largestSbmxShabalin(matrix(sample(vectorA), M, N), m, n)
    pvalue <- pvalue + (permSum >= originSum) / (iter+1)
  }
  pvalue <- pvalue + 1 / (iter+1)
  permTest = (pvalue <= alpha) + 0
  return(list('pvalue' = pvalue, 'test' = permTest))
}

fastTestMonteCarlo <- function(A, rowset, colset, dist, alpha=0.05){
  bof <- length(rowset) * length(colset)
  iniPvalue <- 1
  for (i in rowset){
    for (j in colset){
      renew <- uniPermuTest(A, i, j, alpha=0.05, iter=500)$pvalue
      renew <- biPermuTest(A, i, j, alpha=0.05, iter=500)$pvalue
      iniPvalue <- iniPvalue + (renew < iniPvalue) * (renew - iniPvalue)
    }
  }
  print(iniPvalue)
  print(bof)
  return(list('pvalue' = min(bof * iniPvalue, 1), 'test' = (bof * iniPvalue < alpha) + 0))
}

#------------Main------------
# Parameters
args = commandArgs(TRUE)
dist <- 'Normal'
M <- 200
N <- 100
m <- as.numeric(args[[1]])
n <- as.numeric(args[[2]])
j <- as.numeric(args[[3]]) # iterationfin j, ranging from 1 to 100
mu.len <- 8
iter <- 100

# Critial value and approximation net
crit <- sqrt(2*(m*log(M/m) + n*log(N/n))/(m*n))
mu <- (seq(mu.len)*0.125 + 0.5)*crit
rows <- bitIndex(M)
cols <- bitIndex(N)

# Set seed for reproducibility
seed.num <- 228
set.seed(seed.num)
seeds <- sample(1:10^6, iter)

# fastTestMonteCarlo for iteration j
fast <- list()
for (i in 1:mu.len){
  set.seed(seeds[j])
  if (dist == 'Normal')
    data <- normalMatrix(M, N, m, n, mu[i])
  else if (dist == 'Poisson')
    data <- poissonMatrix(M, N, m, n, mu[i])
  else
    stop("distribution must be either Normal or Poisson")
  fast[[i]] <- fastTestMonteCarlo(data, rows, cols, dist)$pvalue
}

save.image(file=paste(dist, m, n, j, '.RData', sep='_'))
