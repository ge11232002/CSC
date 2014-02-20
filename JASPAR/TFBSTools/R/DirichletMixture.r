


dirichletMixtureEMEstimation <- function(inputMatrix, K, 
                                         alpha0=NULL, pmix=NULL){
## inputMatrix: the samples summarized as counts for each of A letters, N x A
  ## A <- 4 for DNA. This matrix should be concatenated by all matrices 
  ## from Jaspar.
## K: the number of sought component.
## alpha0: the estimated Dirichlet parameters A x K
## pmix: mixing proportions 1 x K

  N <- nrow(inputMatrix)
  A <- ncol(inputMatrix)
  rowSumsInputMatrix <- rowSums(inputMatrix)
  oN <- rep(1, N)
  oK <- rep(1, K)
  oA <- rep(1, A)

  gam <- matrix(0, nrow=N, ncol=K)
  contrib_n <- matrix(0, nrow=A, ncol=K)
  contrib_d <- rep(0, K)

  ## random initialization
  if(is.null(alpha0)){
    epsilon <- 0.1
    alpha0 <- epsilon * (2 * matrix(runif(A * K), nrow=A, ncol=K) - 1) + 2
  }
  Alpha0 <- colSums(alpha0)
  if(is.null(pmix)){
    pmix <- rep(1, K) / K
  }

  ## minimum alpha0 value
  alpha0_min <- 1e-8
  ftol <- 1e-10 
  iteouter_max <- 10000 
  iteinner_max <- 1
  dll_min <- 1e-6 

  ite <- 0; 
  dll <- Inf;
  ll = numeric(iteouter_max)
  while(ite < iteouter_max && (dll > dll_min || ite == 1)){
    ite <- ite + 1;

    ## E-step
    contrib <- log(pmix) + lgamma(Alpha0) - colSums(lgamma(alpha0))
    for(k in 1:K){
      alpha0_k <- alpha0[ , k]
      gam[ ,k] <- rowSums(lgamma(inputMatrix + 
                                 matrix(rep(alpha0_k, N), nrow=N, ncol=A, 
                                        byrow=TRUE))) -
                     lgamma(rowSumsInputMatrix + rep(Alpha0[k], N))
    }
    gam <- gam + matrix(rep(contrib, N), nrow=N, ncol=K, byrow=TRUE)
    maxgam <- apply(gam, 1, max)
    gam <- exp(gam - matrix(rep(maxgam, K), nrow=N, ncol=K, byrow=FALSE))

    ll[ite] <- sum(maxgam + log(rowSums(gam)))

    dll <- ll[ite] - ll[max(ite-1,1)]
    sumgam <- 1 / rowSums(gam)
    gam <- gam * matrix(rep(sumgam, K), nrow=N, ncol=K, byrow=FALSE)

    ## M-step
    pmix <- colSums(gam) / N
    dalpha0 <- Inf  
    iteinner <- 0 
    while(sum(abs(dalpha0)) > ftol && iteinner < iteinner_max){
      print(iteinner)
      iteinner = iteinner + 1
      for(k in 1:K){
        alpha0_k = alpha0[ , k]
        contrib_n[ , k] = psigamma(t(inputMatrix) + 
                             matrix(rep(alpha0_k, N), nrow=A, 
                                    ncol=N, byrow=FALSE)) %*% gam[ , k]
        contrib_d[k] = sum(psigamma(rowSumsInputMatrix + Alpha0[k]) * 
          gam[ ,k]) - pmix[k] * N * psigamma(Alpha0[k])
      }
      contrib_n = contrib_n - N * psigamma(alpha0) * 
        matrix(rep(pmix, A), nrow=A, ncol=length(pmix), byrow=TRUE)
      dalpha0 = pmax(alpha0 * (contrib_n / 
                              matrix(rep(contrib_d, A), nrow=A, 
                                     ncol=length(contrib_d), byrow=TRUE)), 
                    alpha0_min) - alpha0
      alpha0 = alpha0 + dalpha0 
      Alpha0 = colSums(alpha0)
    }
    if(! ite %% iteouter_max /10){
      print(dll)
      print(sum(abs(dalpha0)))
      print(pmix)
      print(alpha0)
    }
  }
  return(list(alpha0=alpha0, pmix=pmix, ll=ll))
}



