

library('CIBERSORT')


#' Main functions
#'
#' The Main function of CIBERSORT
#' @param sig_matrix  sig_matrix file path to gene expression from isolated cells, or a matrix of expression profile of cells.
#'
#' @param mixture_file mixture_file file path to heterogenous mixed expression file, or a matrix of heterogenous mixed expression
#'
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#' @import utils
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom stats sd
#' @export
#' @examples
#' \dontrun{
#'   ## example 1
#'   sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
#'   mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
#'   results <- cibersort(sig_matrix, mixture_file)
#'   ## example 2
#'   data(LM22)
#'   data(mixed_expr)
#'   results <- cibersort(sig_matrix = LM22, mixture_file = mixed_expr)
#' }

cibersort_run <- function(sig_matrix, mixture_file, perm = 0, QN = TRUE){


  #read in data
  if (is.character(sig_matrix)) {
    X <- read.delim(sig_matrix, header=T, sep="\t", row.names=1, check.names = F)
    X <- data.matrix(X)
  } else {
    X <- sig_matrix
    tryCatch({       rm('sig_matrix')        }, error=function(cond){       111         })
    # X <- data.matrix(X)
    X = as.matrix(X)
  }

  if (is.character(mixture_file)) {
    Y <- read.delim(mixture_file, header=T, sep="\t", row.names=1, check.names = F)
    Y <- data.matrix(Y)
  } else {
    Y <- mixture_file
    tryCatch({       rm('mixture_file')        }, error=function(cond){       111         })
    # Y <- data.matrix(Y)
    gc()
    Y = as.matrix( Y )
  }


  #order
  X <- X[order(rownames(X)),]
  Y <- Y[order(rownames(Y)),]

  P <- perm #number of permutations

  #anti-log if max < 50 in mixture file
  if(max(Y) < 50) {Y <- 2^Y}

  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    library("preprocessCore")
    Y = as.matrix( Y )     ##### 可以不要
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]

  #standardize sig matrix
  X = as.matrix( X )    ##### 可以不要
  X <- (X - mean(X)) / sd(as.vector(X))

  #empirical null distribution of correlation coefficients
  if (P > 0) {
    print('     [sorting]...，下面可能会出错【could not find function "future_map"】')
    nulldist <- sort(  doPerm(P, X, Y)$dist  )
    print('     [sorting completed]')
  }

  #print(nulldist)

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  #iterate through mixtures
  print('     [iterate through mixtures]...    ')
  while (itor <= mixtures) {

    y <- Y[,itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    result <- CoreAlg(X, y)

    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    #calculate p-value
    if (P > 0) {
      pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
    }

    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }

    Pett_ProgressBar(itor, ncol(Y))
    itor <- itor + 1
  }

  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  return(obj)
}






#' do permutations
#'
#' Do the permutations analysis
#' @param perm Number of permutations
#' @param X cell-specific gene expression
#' @param Y mixed expression per sample
#' @importFrom purrr reduce map
#' @importFrom stats sd
#' @export

doPerm_recursion_method <- function(perm, X, Y){
  ### 运用了递归的算法
  
  # library(purrr)
  
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  itorect <- function(Ylist, X) {
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)
    
    mix_r <- result$mix_r
    
    return(mix_r)
  }
  
  if (perm == 1) {
    dist <- itorect(Ylist = Ylist, X = X)
  } else {
    ### Error in future_map(1:svn_itor, res) : 
    ### could not find function "future_map"
    dist <- purrr::map(1:perm, ~ itorect(Ylist = Ylist, X = X)) %>%
      purrr::reduce(rbind)
  }
  
  newList <- list("dist" = dist)
}



doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  
  while(itor <= perm){
    #print(itor)
    
    setTxtProgressBar( txtProgressBar( style=3 ), itor/perm )
    
    
    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])
    
    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)
    
    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)
    
    mix_r <- result$mix_r
    
    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}
    
    itor <- itor + 1
  }
  tryCatch({       rm('Ylist','yr','result','mix_r')        }, error=function(cond){       111         })
  gc()
  newList <- list("dist" = dist)
}










#' Core algorithm
#'
#' The core algorithm of CIBERSORT which is used svm
#' @param X cell-specific gene expression
#' @param y mixed expression per sample
#' @importFrom furrr future_map
#' @importFrom future availableCores
#'@importFrom stats cor
#' @import e1071
#' @export
CoreAlg_you_wen_ti__could_not_find_function_future_map <- function(X, y){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  enableParallel()
  
  if (Sys.info()['sysname'] == 'Windows') {
    out <- future_map(1:svn_itor, res)
  } else {
    if (svn_itor <= availableCores() - 2) {
      enableParallel(nThreads = svn_itor)
    } else {
      enableParallel()
    }
    out <- future_map(1:svn_itor, res)
  }
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w <- weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}



CoreAlg <- function(X, y){
  
  #try different values of nu
  svn_itor <- 3
  
  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-e1071::svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }
  
  if(Sys.info()['sysname'] == 'Windows') out <- parallel::mclapply(1:svn_itor, res, mc.cores=1) else
    out <- parallel::mclapply(1:svn_itor, res, mc.cores=svn_itor)
  
  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)
  
  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }
  
  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  
  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))
  
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]
  
  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)
  
}







#' open multiple threads process
#' @param nThreads The number of threads
#' @param maxSize The max memory size in global, default 500MB, the unit is MB
#' @param verbose output other useful information
#' @import future
#' @importFrom future availableCores plan
#' @export
#' @examples
#' \dontrun{
#'   enableParallel()
#' }
enableParallel <- function(nThreads = NULL, maxSize = 500, verbose = FALSE) {
  nCores <- future::availableCores()
  options(future.globals.maxSize = maxSize*1024^2)
  if (verbose) message("The maxSize is ", maxSize, " Mb.")
  if (is.null(nThreads)) {
    if (nCores < 4) {
      nThreads <- nCores
    } else {
      nThreads <- nCores - 2
    }
  }
  if (!is.numeric(nThreads) || nThreads < 2)
    stop("nThreads must be numeric and at least 2.")
  if (nThreads > nCores) {
    if (verbose) {
      message("Requested number of threads is higher than number of available processors (or cores)")
      message(paste("The max working processes is", nCores))
    }
  }
  if (verbose) message(paste("Allowing parallel execution with up to", nThreads, "working processes."))
  future::plan("multicore", workers = nThreads)
  invisible(nThreads)
}

















