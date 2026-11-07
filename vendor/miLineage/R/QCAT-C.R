library(abind)
library(CompQuadForm)
########################################
#                                      #
#   QCAT-C test for clustered data     #
#             add in v3.0              #
#                                      #
#                                      #
########################################
.diag2 <- function(x){
  
  if(length(x)>1){
    return(diag(x))
  }else{
    return(as.matrix(x))
    
  }
  
}

.F.test <- function(x){
  
  x.stat = -2 * sum(log(x))
  return( 1 - pchisq(x.stat, df = 2 * length(x)) )
}

.simes.test <- function(x){
  
  return( min(length(x) * x/rank(x)) )
  
}


pchisqsum2.davies <- function(Q, lambda){
  
  delta = rep(0, length(lambda))
  acc = 1e-07
  delta <- delta[lambda > 0]
  lambda <- lambda[lambda > 0]
  
  tmp <- CompQuadForm::davies(q = Q, lambda = lambda, 
                              delta = delta, acc = acc)
  if (tmp$ifault > 0) {
    lambda <- zapsmall(lambda, digits = 2)
    delta <- delta[lambda > 0]
    lambda <- lambda[lambda > 0]
    tmp <- CompQuadForm::farebrother(q = Q, lambda = lambda, 
                                     delta = delta)
  }
  Qq <- if ("Qq" %in% names(tmp)) 
    tmp$Qq
  else tmp$res
  return(list(p = Qq, errflag = 0))
  
  
}


QCAT.Cluster <- function(ID, OTU, OTU.base=NULL, X, X.index, Tax=NULL, min.depth=0, perm.type=NULL, n.perm=NULL, fdr.alpha=0.05, test="chisq"){
  
  id = ID
  n.resample = n.perm
  if(is.null(perm.type)){
    perm.type="Simple"
  }
  
  if(!is.matrix(OTU)){
    warning("OTU table is not a matrix")
    OTU = as.matrix(OTU)  
  }
  
  
  if(!is.matrix(X)){
    warning("Covariate table is not a matrix")
    X = as.matrix(X)  
  }  
  
  if(missing(X.index)){
    X.index = 1:ncol(X)
  }
  
  if(nrow(OTU)!=nrow(X)){
    stop("Samples in the OTU table and the covariate table should be the same.")  
  }
  
  remove.subject = which(rowSums(OTU)<min.depth)
  if(length(remove.subject)>0){
    
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth ))
    X = X[-remove.subject, ,drop=FALSE]
    OTU = OTU[-remove.subject, ,drop=FALSE]
    if( !is.null(OTU.base) ){
      OTU.base = OTU.base[-remove.subject, ,drop=FALSE]
    }
    
  }
  
  keep = which(colSums(OTU)>0)
  count = OTU[,keep, drop=FALSE]
  if( !is.null(OTU.base) ){
    count.base = OTU.base[,keep,drop=FALSE]
  }else{
    count.base = NULL
  }
  
  X = cbind(1, X) # add the intercept term
  X.index = X.index + 1
  
  if(is.null(Tax)){ # perform one test using all OTUs
    
    if(is.null(n.resample)){ # asymptotic test only
      
      pval = as.matrix( .Score.test.cluster(id, perm.type, count, Y.base=count.base, X, X.index, resample=FALSE, n.replicates=NULL, test=test)$score.pvalue )
      colnames(pval) = "Asymptotic"
      
    }else{ # resampling test + asymptotic test 
      
      tmp = .Score.test.cluster(id, perm.type, count, Y.base=count.base, X, X.index, resample=TRUE, n.replicates=n.resample, test=test)
      pval = c(tmp$score.pvalue, tmp$score.Rpvalue)
      names(pval) = c("Asymptotic", "Resampling")
      
    }
    
    return( list(pval=pval) )
    
  }else{ # perform tests for lineages
    
    if(!is.matrix(Tax)){
      warning("Tax table is not a matrix")
      Tax = as.matrix(Tax)  
    }
    
    
    tax = Tax[keep, ,drop=FALSE]
    
    if( sum(colnames(count)!=rownames(tax))>0 ){
      
      stop("Error: OTU IDs in OTU table are not consistent with OTU IDs in Tax table")
    }
    
    W.data = data.table(data.frame(tax, t(count)))
    n.rank = ncol(tax)
    otucols = names(W.data)[-(1:n.rank)]
    if(!is.null(OTU.base) ){
      W.data.base = data.table(data.frame(tax, t(count.base)))
      otucols.base = names(W.data.base)[-(1:n.rank)]
    }
    
    n.level = n.rank-1
    
    subtree = NULL
    pval = NULL
    
    for(k in 1:n.level){
      
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)    
      
      tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
      W.count = data.matrix(tt[, otucols, with=FALSE])      
      if(!is.null(OTU.base) ){
        tt.base = W.data.base[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols.base, by=list( get(Rank.low), get(Rank.high) )]
        W.count.base = data.matrix(tt.base[, otucols.base, with=FALSE])         
      }
      
      for(j in 1:m.level){
        
        Y = t(W.count[which(W.tax == level.uni[j]), , drop=FALSE])
        
        if(!is.null(OTU.base) ){
          Y.base = t(W.count.base[which(W.tax == level.uni[j]), , drop=FALSE])
        }
        #Y = t(W.count[which(W.tax == "o__Actinomycetales"), , drop=FALSE])
        #Y.base = t(W.count.base[which(W.tax == "o__Actinomycetales"), , drop=FALSE])
        remove.index = which(colSums(Y)==0)
        
        if(length(remove.index)==ncol(Y)){
          
          #print("==skip:0==");
          next
          
          
        }else{
          
          if(length(remove.index)>0){
            Y = Y[, -remove.index, drop=FALSE] 
            if(!is.null(OTU.base) ){
              Y.base = Y.base[, -remove.index, drop=FALSE]
            }
          }
          if(is.null(OTU.base) ){
            Y.base = NULL
          }
          
          if(ncol(Y)==1){
            
            next
            #print("==skip:1==");
            
          }else{
            
            subtree = c(subtree, level.uni[j])
            
            if(is.null(n.resample)){ # asymptotic test only
              
              pval = cbind(pval, .Score.test.cluster(id, perm.type, Y, Y.base=Y.base, X, X.index, resample=FALSE, n.replicates=NULL, test=test)$score.pvalue)
              
              
            }else{ # resampling test + asymptotic test 
              
              tmp = .Score.test.cluster(id, perm.type, Y, Y.base=Y.base, X, X.index, resample=TRUE, n.replicates=n.resample, test=test)
              pval = cbind(pval, c(tmp$score.pvalue, tmp$score.Rpvalue) )
              
              
            }
            
          }
          
        }
        
        
      }# lineage loop
      
      
    }# level loop
    
    
    colnames(pval) = subtree
    
    if(is.null(n.resample)){
      
      rownames(pval) = "Asymptotic"
      score.tmp = pval[1,]
      
    }else{
      
      rownames(pval) = c("Asymptotic", "Resampling")
      score.tmp = pval[2,]
    }
    
    #print(pval)
    
    # identify significant lineages
    subtree.tmp = subtree
    index.na = which(is.na(score.tmp))
    if(length(index.na)>0){
      score.tmp = score.tmp[-index.na]
      subtree.tmp = subtree.tmp[-index.na]
    }
    
    #score.tmp[score.tmp==0] = 1e-4
    m.test = length(score.tmp)
    
    # Benjamini-Hochberg FDR control
    index.p = order(score.tmp)
    p.sort = sort(score.tmp)
    #fdr.alpha = 0.05
    
    # change 04/17/2016
    reject = rep(0, m.test)
    tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
    if(length(tmp)>0){
      index.reject = index.p[1:max(tmp)]
      reject[index.reject] = 1      
    }
    
    sig.lineage = subtree.tmp[reject==1]
    
    
    # perform global test
    global.score.fisher = .F.test(score.tmp)
    global.score.min =  .simes.test(score.tmp)
    global.pval = c(global.score.fisher, global.score.min)
    names(global.pval) = c("Fisher", "Simes")
    
    return( list(lineage.pval=pval, sig.lineage=sig.lineage, global.pval=global.pval) )
    
    
  }
  
}


.Ei.beta.XX <- function(m, p.XX, beta, XX.i, Y.i){
  
  Ei.out = rep(NA,m)
  
  idx1 = 0
  idx2 = 0
  
  for(j in 1:(m-1)){
    
    idx1 = idx2 + 1
    idx2 = idx2 + p.XX[j]
    
    Ei.out[j] = exp(beta[idx1:idx2] %*% XX.i[,1:p.XX[j],j])
    
    
  }
  
  
  Ei.out[m] = 1
  
  
  
  
  return (Ei.out)
}


.fun.neg.loglik.beta.XX <- function(beta, data){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  
  #n.beta = (m - 1)*p
  n.beta = sum(p.XX)
  loglik = 0
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    for(i in 1:n){
      
      E.i = .Ei.beta.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Y.pos.index = which(Y[i,]>0)
      #loglik = loglik + Y[i,Y.pos.index] %*% log(P.i[Y.pos.index])
      loglik = loglik + Y[i,] %*% log(P.i)
    }
    
  }
  
  return (-loglik)
  
}


.fun.neg.score.beta.XX <- function(beta, data){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  
  n.beta = sum(p.XX)
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta = rep(0, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Score.beta = Score.beta + kronecker( matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1))
      tmp = NULL
      for(j in 1:(m-1)){
        tmp = c(tmp, (Y[i,j] - nY[i]*P.i[j]) * XX[i,1:p.XX[j],j] )
      }
      Score.beta = Score.beta + tmp
    }
    
    return (-Score.beta)
  }
  
  
  
}


.fun.score.i.beta.XX <- function(beta, data){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  
  #n.beta = (m - 1)*p
  n.beta = sum(p.XX)
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta.i = matrix(0, n, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      # add 03/28/2016
      #       if(sum.E.i==0){
      #         P.i = rep(0,m)
      #       }
      
      #Score.beta.i[i,] =  kronecker( matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1) )
      tmp = NULL
      for(j in 1:(m-1)){
        tmp = c(tmp, (Y[i,j] - nY[i]*P.i[j]) * XX[i,1:p.XX[j],j] )
      }
      Score.beta.i[i,] = tmp
      
    }
    
    return (Score.beta.i)
  }
  
  
  
}


.fun.hessian.beta.XX <- function(beta, data, save.list=FALSE){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  #n.beta = (m-1)*p
  n.beta = sum(p.XX)
  
  #idx.end = (1:(m-1))*p 
  #idx.start = idx.end - (p-1)
  
  idx.start = c(0,cumsum(p.XX[-(m-1)]))+1
  idx.end = idx.start + p.XX - 1
  
  if(length(beta)!=n.beta){
    print("Waring: dim of beta is not the same as beta\n")
    
  }else{
    
    Hessian.beta = matrix(0, nrow=n.beta, ncol=n.beta)
    nY = rowSums(Y)
    I.beta.list = list()
    
    for(i in 1:n){
      
      E.i = .Ei.beta.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      ## tmp.beta
      #tmp.beta =  (E.i[-m] %o% E.i[-m])*nY[i]/sum.E.i^2 
      tmp.beta =  as.matrix(P.i[-m] %o% P.i[-m])
      diag(tmp.beta) = diag(tmp.beta) - P.i[-m]
      tmp.beta = nY[i] * tmp.beta
      #tmp.beta[is.na(tmp.beta)] = 0  ## add 03/28/2016
      
      #Hessian.beta = Hessian.beta + kronecker( tmp.beta, ( X[i,] %o% X[i,] ) ) 
      Hessian.tmp = matrix(NA, nrow=n.beta, ncol=n.beta)
      for(j in 1:(m-1)){
        for(k in j:(m-1)){
          Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]] = tmp.beta[j,k] * (XX[i,1:p.XX[j],j] %o% XX[i,1:p.XX[k],k])
          Hessian.tmp[idx.start[k]:idx.end[k], idx.start[j]:idx.end[j]] = t(Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]])
        }
      }
      Hessian.beta = Hessian.beta + Hessian.tmp
      
      if(save.list){
        I.beta.list[[i]] = tmp.beta
        
      }
      
    }
    
    
    if(save.list){
      
      return ( list(Hessian.beta=Hessian.beta, I.beta.list = I.beta.list) )
      
    }else{
      
      return (Hessian.beta)
    }
    
  }
  
  
}


# 03/28/2018 modify to allow for uneven #obs across clusters
.sampling.strata.uneven <- function(var, id){
  
  id.unique = unique(id)
  n.cluster = length(id.unique)
  var.sample = var
  for(i in 1:n.cluster){
    index = which(id==id.unique[i])
    var.sample[index] = var.sample[sample(index)]
  }
  return(var.sample)
  
}


########## functions defined for clustered data ##########
.Score.test.stat.cluster <- function(id, Y, XX, p.XX, X.par.index){
  
  #p = dim(XX)[2]
  
  nY = rowSums(Y)
  
  n = nrow(Y)
  m = ncol(Y)  
  n.beta = sum(p.XX)
  
  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }  
  
  if(is.null(X.par.index) || n==0){
    
    score.stat.beta = NA
    
  }else{
    
    #XX.reduce = XX[,-X.par.index, , drop=FALSE]
    
    p.XX.reduce = p.XX - length(X.par.index) 
    XX.reduce = array(NA,dim=c(n,dim(XX)[2],m-1))
    for(j in 1:(m-1)){
      XX.reduce[,1:p.XX.reduce[j],j] = XX[,(1:p.XX[j])[-X.par.index],j]
    }    
    
    
    #par.interest.index.beta =  kronecker( ((0:(m-2))*p), rep(1,length(X.par.index))) + X.par.index
    par.interest.index.beta = NULL
    tmp = 0
    for(j in 1:(m-1)){
      par.interest.index.beta = c(par.interest.index.beta, tmp + X.par.index)
      tmp = sum(p.XX[1:j])
    }
    
    n.par.interest.beta = length(par.interest.index.beta) 
    
    #beta.ini.reduce = rep(0, (p.reduce*(m-1)))    
    beta.ini.reduce = rep(0, sum(p.XX.reduce) )    
    
    
    data.reduce.beta = list(Y=Y, XX=XX.reduce, p.XX=p.XX.reduce)
    
    est.reduce.beta = rep(NA, n.beta)
    est.reduce.beta[par.interest.index.beta] = 0 
    # change 04/08/2016
    est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta.XX, gr=.fun.neg.score.beta.XX, data = data.reduce.beta, method="BFGS")$par
    
    
    data.beta = list(Y=Y, XX=XX, p.XX=p.XX)
    Score.reduce.beta = .fun.score.i.beta.XX(est.reduce.beta, data.beta)
    
    # for resampling: S.beta.list, I.beta.list
    # S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, ((1:(m-1))*p-p+1)])
    idx.start = c(0,cumsum(p.XX[-(m-1)]))+1
    S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, idx.start])
    tmp = .fun.hessian.beta.XX(est.reduce.beta, data.beta, save.list=TRUE)
    I.beta.list = tmp$I.beta.list  
    
    Hess.reduce.beta = tmp$Hessian.beta
    #Hess.reduce.beta =  .fun.hessian.beta.pos(est.reduce.beta, data.beta)
    
    # re-organized the score statistics and Hessian matrix
    Score.reduce.reorg = cbind( matrix(Score.reduce.beta[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
    Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.beta[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                              cbind( matrix(Hess.reduce.beta[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
    
    
    A = colSums(Score.reduce.reorg)[1:n.par.interest.beta]
    
    B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
    
    
    B2 =  matrix(0, n.beta, n.beta)
    
    # 03/12/2018   modify for clustered data (same #observation for each cluster for now)
    # 03/28/2018   modify to allow for uneven #obs across clusters
    id.unique = unique(id)
    n.cluster = length(id.unique)
    # TT = n/length(unique(id))  
    # n.cluster = n/TT
    # idx1 = 1
    # idx2 = TT
    for(i in 1:n.cluster){
      
      #Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[idx1:idx2,,drop=FALSE] )
      Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[which(id==id.unique[i]),,drop=FALSE] )
      B2 = B2 + Score.reduce.reorg.cluster %o%  Score.reduce.reorg.cluster
      
      #idx1 = idx2 + 1
      #idx2 = idx2 + TT
      
    }
    
    B = B1 %*% B2 %*% t(B1)
    score.stat.beta = A %*% ginv(B) %*% A
    
    
  }
  
  
  return(list(score.stat.beta=score.stat.beta, S.beta.list=S.beta.list, I.beta.list=I.beta.list, U=A, V=B )   )
  
  
}


.Score.test.cluster <- function(id, perm.type, Y, Y.base=NULL, X, X.par.index, seed=11, resample=FALSE, n.replicates=NULL, test="chisq"){
  
  n = nrow(Y)
  m = ncol(Y) 
  
  
  if(is.null(Y.base)){
    
    XX = array(NA,dim=c(n,ncol(X),m-1))
    
    
    for(j in 1:(m-1)){
      XX[,,j] = X
      
    }
    
    p.XX = rep(ncol(X), m-1)
    
  }else{
    
    XX = array(NA,dim=c(n,ncol(X)+1,m-1))
    
    p.XX = rep(0, m-1)
    
    for(j in 1:(m-1)){
      
      XX[,1:ncol(X),j] = X
      
      if(length(table(Y.base[,j]))==1){
        p.XX[j] = ncol(X)
      }else{
        p.XX[j] = ncol(X)+1
      }
      
    }
    
    XX[,(ncol(X)+1),] = Y.base[,-m,drop=FALSE]
  }
  
  #p = dim(XX)[2]
  
  nY = rowSums(Y)
  
  ## remove 03/28/2016
  #   nY0.index = which(nY==0)
  #   if(length(nY0.index)>0){
  #     Y = Y[-nY0.index, , drop=FALSE]
  #     X = X[-nY0.index, , drop=FALSE]
  #   }
  
  #n.beta = (m - 1)*p
  
  if(sum(X.par.index == 1)){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }  
  
  if(is.null(X.par.index) || n==0){
    
    score.stat.beta = NULL
    score.pvalue.beta = NA
    n.par.interest.beta = 0
    
  }else{
    
    tmp.one = try( .Score.test.stat.cluster(id, Y, XX, p.XX, X.par.index) )
    #tmp.one = ( .Score.test.stat.cluster(id, Y, XX, p.XX, X.par.index) )
    if( "try-error" %in% class(tmp.one)){
      
      score.stat.beta = NA
      score.pvalue.beta = NA
      n.par.interest.beta = NA
      
    }else{
      
      n.par.interest.beta = (m-1)*length(X.par.index)
      
      if(test=="chisq"){
        score.stat.beta = tmp.one$score.stat.beta
        score.pvalue.beta = 1 - pchisq(score.stat.beta,  n.par.interest.beta) 
      }
      if(test=="vc"){
        inv.V = ginv(tmp.one$V)
        score.stat.beta =  tmp.one$U %*% inv.V %*% inv.V %*% tmp.one$U
        lambda = eigen(inv.V)$values
        score.pvalue.beta = pchisqsum2.davies(Q=score.stat.beta, lambda=lambda)$p 
      }
      
    }
    
    
  }
  
  
  beta.results = list(score.stat = score.stat.beta, score.pvalue = score.pvalue.beta, df =n.par.interest.beta)
  #  print(score.stat.beta)
  
  if(resample){
    ########################################
    #                                      #
    #        Resampling Score Test         #
    # change 07/02/2016 for adaptive resampling
    # change 07/17/2016 for adding vc test #
    #                                      #
    ########################################  
    set.seed(seed)
    if(!is.na(score.stat.beta)){
      
      n.one = 0
      one.acc = 0
      
      start.nperm = 1;
      end.nperm = min(100,n.replicates);
      flag = 1
      ### 07/19/2018, modify while loop
      while(flag){
        
        results = .resample.work.one.cluster(id, perm.type, XX, p.XX, X.par.index, score.stat.beta, tmp.one$S.beta.list, tmp.one$I.beta.list, start.nperm, end.nperm, n.one, one.acc, test=test)
        n.one = results$n.one.new
        one.acc = results$one.acc.new
        flag = results$flag
        next.end.nperm = results$next.end.nperm
        
        
        start.nperm = end.nperm + 1;
        end.nperm = next.end.nperm;
        
        if(start.nperm >= n.replicates){
          flag = 0
        }else if(end.nperm > n.replicates){ 
          #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))
          results = .resample.work.one.cluster(id, perm.type, XX, p.XX, X.par.index, score.stat.beta, tmp.one$S.beta.list, tmp.one$I.beta.list, start.nperm, n.replicates, n.one, one.acc, test=test)
          n.one = results$n.one.new
          one.acc = results$one.acc.new
          flag = 0
        }
        
      }
      
      
      #      print(paste("Final number of resamplings: ", n.one) )
      #      if(n.one<end.nperm/2){
      #        print("Number of resamplings too small for one-part test")
      #      }
      
      tmp = (one.acc+1)/(n.one+1)
      
      #print(n.one)
      #print(one.acc)
      
    }else{
      
      tmp = NA
    }
    
    beta.results = c(beta.results, score.Rpvalue = tmp)
    
    
  }
  
  return(beta.results)
  
  
}


.resample.work.one.cluster <- function(id, perm.type, XX, p.XX, X.par.index, score.stat.beta, S.beta.list, I.beta.list, start.nperm, end.nperm, n.one, one.acc, test){
  
  n = dim(XX)[1]
  id.unique = unique(id)
  n.cluster = length(id.unique)
  XX.interest = XX[,X.par.index, ,drop=FALSE]
  
  # 03/28/2018   modify to allow for uneven #obs across clusters
  if(perm.type=="BTW"){
    
    cluster.nobs = NULL
    # if covaraites of interest are not the same within cluster, cannot perform between cluster permutation
    
    for(i in 1:n.cluster){
      tmp = XX.interest[which(id==id.unique[i]),,,drop=FALSE]
      cluster.nobs = c(cluster.nobs, dim(tmp)[1])
      if(dim(unique(tmp))[1]!=1){
        stop("Cannot perform between cluster permutation: covaraites of interest are not the same within cluster")
      }
      if(i==1){
        XX.interest.unique = unique(tmp)
      }else{
        XX.interest.unique = abind(XX.interest.unique, unique(tmp), along=1)
      }
      
      
    }
    
  }
  
  n.one.new = n.one
  one.acc.new = one.acc
  
  for(k in start.nperm:end.nperm){
    
    # 03/28/2018   modify to allow for uneven #obs across clusters
    # 03/29/2018   modify to allow taxon-specific covariate
    if(perm.type=="BTW"){
      
      #XX.interest.perm = NULL
      sample.cluster.idx = sample(1:n.cluster)
      
      for(i in 1:n.cluster){
        nobs = cluster.nobs[i]
        XX.interest.perm.i = array(NA, dim=c(nobs, dim(XX.interest.unique)[2], dim(XX.interest.unique)[3]))
        for(j in 1:dim(XX.interest.unique)[3]){
          XX.interest.perm.i[,,j] =  matrix(rep(XX.interest.unique[sample.cluster.idx[i],,j],nobs),nrow=nobs, byrow = TRUE)
        }
        if(i==1){
          XX.interest.perm = XX.interest.perm.i
        }else{
          XX.interest.perm = abind(XX.interest.perm, XX.interest.perm.i, along =1)
        }
        
      }
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.interest.perm
      
    }else if(perm.type=="WTH"){
      
      sample.idx = .sampling.strata.uneven(1:n, id)
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.perm[sample.idx,X.par.index,]
      
    }else{
      
      sample.idx = sample(1:n)
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.perm[sample.idx,X.par.index,]
    }
    
    
    
    #X.perm = X
    #X.perm[,X.par.index] = X.perm[sample.idx,X.par.index]
    
    score.stat.beta.perm = try( .Score.test.stat.4Gresampling.cluster(id, XX.perm, p.XX, X.par.index, S.beta.list, I.beta.list, test=test) )
    #score.stat.beta.perm = ( .Score.test.stat.4Gresampling.cluster(id, XX.perm, p.XX, X.par.index, S.beta.list, I.beta.list) )
    
    if(!( "try-error" %in% class(score.stat.beta.perm))){
      
      n.one.new = n.one.new + 1
      if(score.stat.beta.perm >= score.stat.beta){
        one.acc.new = one.acc.new + 1
        
      }  
    }
  }
  
  ### modify 07/19/2018
  if(one.acc.new < 10){
    next.end.nperm = end.nperm * 5;
    flag = 1;
    
  }else{
    next.end.nperm = end.nperm;
    flag = 0;  
  }
  
  return(list(n.one.new=n.one.new, one.acc.new=one.acc.new, flag=flag, next.end.nperm=next.end.nperm))
  
}


.Score.test.stat.4Gresampling.cluster <- function(id, XX.perm, p.XX, X.par.index, S.beta.list, I.beta.list, test){
  
  n = dim(XX.perm)[1]
  #p = dim(XX.perm)[2]
  
  m.beta = length(S.beta.list[[1]])
  
  #idx.end = (1:m.beta)*p 
  #idx.start = idx.end - (p-1)
  idx.start = c(0,cumsum(p.XX[-m.beta]))+1
  idx.end = idx.start + p.XX - 1
  
  #n.par.interest.beta = m.beta
  #n.beta = m.beta*p
  n.beta = sum(p.XX)
  #par.interest.index.beta =  kronecker( ((0:(m.beta-1))*p), rep(1,length(X.par.index))) + X.par.index
  par.interest.index.beta = NULL
  tmp = 0
  for(j in 1:m.beta){
    par.interest.index.beta = c(par.interest.index.beta, tmp + X.par.index)
    tmp = sum(p.XX[1:j])
  }
  
  n.par.interest.beta = length(par.interest.index.beta) 
  
  
  Score.reduce.beta.perm = matrix(0, n, n.beta )
  Hess.reduce.beta.perm = matrix(0, n.beta, n.beta )
  
  for(i in 1:n){
    
    ###################################################
    #                                                 #
    #         Beta part: resampling Score test        #
    #                                                 #
    ################################################### 
    #Score.reduce.beta.perm[i,] = Score.reduce.beta.perm[i,] + kronecker(matrix(S.beta.list[[i]], ncol=1),  matrix(X.perm[i,], ncol=1))  
    tmp = NULL
    for(j in 1:m.beta){
      tmp = c(tmp, S.beta.list[[i]][j] * XX.perm[i,1:p.XX[j],j] )
    }
    Score.reduce.beta.perm[i,] = tmp
    
    #Hess.reduce.beta.perm = Hess.reduce.beta.perm + kronecker(I.beta.list[[i]], (  X.perm[i,] %o% X.perm[i,] ) )
    
    Hessian.tmp = matrix(NA, nrow=n.beta, ncol=n.beta)
    for(j in 1:m.beta){
      for(k in j:m.beta){
        Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]] = I.beta.list[[i]][j,k] * (XX.perm[i,1:p.XX[j],j] %o% XX.perm[i,1:p.XX[k],k])
        Hessian.tmp[idx.start[k]:idx.end[k], idx.start[j]:idx.end[j]] = t(Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]])
      }
    }
    Hess.reduce.beta.perm = Hess.reduce.beta.perm + Hessian.tmp
    
    
    #     if(sum(is.na(Hess.reduce.beta.perm))>0){
    #       print(i); break;
    #       
    #     }
    
  }
  
  ###################################################
  #                                                 #
  #         Beta part: resampling Score test        #
  #                                                 #
  ################################################### 
  # re-organized the score statistics and Hessian matrix
  Score.reduce.beta.perm.reorg = cbind( matrix(Score.reduce.beta.perm[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta.perm[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
  Hess.reduce.beta.perm.reorg = rbind(cbind( matrix(Hess.reduce.beta.perm[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta.perm[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                                      cbind( matrix(Hess.reduce.beta.perm[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta.perm[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
  
  
  A = colSums(Score.reduce.beta.perm.reorg)[1:n.par.interest.beta]
  
  B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.beta.perm.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.beta.perm.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
  
  
  B2 =  matrix(0, n.beta, n.beta)
  
  # change in 04/08/2016
  # 03/28/2018   modify to allow for uneven #obs across clusters
  id.unique = unique(id)
  n.cluster = length(id.unique)
  #TT = n/length(unique(id))  ## same #observation for each cluster for now
  #n.cluster = n/TT
  #idx1 = 1
  #idx2 = TT
  for(i in 1:n.cluster){
    
    #Score.reduce.beta.perm.reorg.cluster =  colSums( Score.reduce.beta.perm.reorg[idx1:idx2,,drop=FALSE] )
    Score.reduce.beta.perm.reorg.cluster =  colSums( Score.reduce.beta.perm.reorg[which(id==id.unique[i]),,drop=FALSE] )
    B2 = B2 + Score.reduce.beta.perm.reorg.cluster %o%  Score.reduce.beta.perm.reorg.cluster
    
    #idx1 = idx2 + 1
    #idx2 = idx2 + TT
    
  }
  
  
  B = B1 %*% B2 %*% t(B1)
  
  if(test=="chisq"){
    score.stat.beta.perm = A %*% ginv(B) %*% A
  }
  
  if(test=="vc"){
    
    inv.B = ginv(B)
    score.stat.beta.perm = A %*% inv.B %*% inv.B %*% A
  }
  
  
  return(score.stat.beta.perm)
  
  
  
}




########################################
#                                      #
#     QCAT_GEE for clustered data      #
#          two-part test               #
#                                      #
#                                      #
########################################

QCAT_GEE.Cluster <- function(ID, OTU, OTU.base=NULL, X, X.index, Z, Z.index, Tax=NULL, min.depth=0, perm.type=NULL, n.perm=NULL, fdr.alpha=0.05, test="chisq"){
  
  id = ID
  n.resample = n.perm
  if(test=="mix"){
    test="vc"
  }
  if(is.null(perm.type)){
    perm.type="Simple"
  }
  
  if(!is.matrix(OTU)){
    warning("OTU table is not a matrix")
    OTU = as.matrix(OTU)  
  }
  
  
  if(!is.matrix(X)){
    warning("Covariate table is not a matrix")
    X = as.matrix(X)  
  }  
  
  if(missing(X.index)){
    X.index = 1:ncol(X)
  }
  
  if(missing(Z)){
    Z = X
  }else{
    
    if(!is.matrix(Z)){
      warning("Covariate table Z is not a matrix")
      Z = as.matrix(Z)  
    }  
  }
  
  
  
  if(missing(Z.index)){
    Z.index = 1:ncol(Z)
  }
  
  if(nrow(OTU)!=nrow(X)){
    stop("Samples in the OTU table and the covariate table should be the same.")  
  }
  
  remove.subject = which(rowSums(OTU)<min.depth)
  if(length(remove.subject)>0){
    
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth ))
    X = X[-remove.subject, ]
    OTU = OTU[-remove.subject, ,drop=FALSE]
    if( !is.null(OTU.base) ){
      OTU.base = OTU.base[-remove.subject, ,drop=FALSE]
    }
  }
  
  keep = which(colSums(OTU)>0)
  count = OTU[,keep,drop=FALSE]
  if( !is.null(OTU.base) ){
    count.base = OTU.base[,keep,drop=FALSE]
  }else{
    count.base = NULL
  }
  
  X = cbind(1, X) # add the intercept term
  X.index = X.index + 1
  
  Z = cbind(1, Z) # add the intercept term
  Z.index = Z.index + 1
  
  if(is.null(Tax)){ # perform one test using all OTUs
    
    if(is.null(n.resample)){ # asymptotic test only
      
      tmp = .Score.test2.cluster(id, perm.type, count, Y.base=count.base, X, X.index, Z, Z.index, seed=11, resample=FALSE, n.replicates=NULL, test=test)
      
      pval.zero = as.matrix( tmp$zero.results$score.pvalue )
      pval.pos = as.matrix( tmp$pos.results$score.pvalue )
      colnames(pval.zero) = "Asymptotic"
      colnames(pval.pos) = "Asymptotic"
      
      if(test=="chisq"){
        pval.comb = as.matrix( tmp$comb.results$score.pvalue )
        colnames(pval.comb) = "Asymptotic"
      }
      if(test=="vc"){
        #pval.comb = as.matrix( c(tmp$comb.results$scoreF.pvalue, tmp$comb.results$scoreS.pvalue), ncol=2 )
        #colnames(pval.comb) = c( "Asymptotic.Fisher", "Asymptotic.Simes" )   
        
        pval.comb = as.matrix( tmp$comb.results$score.pvalue )
        colnames(pval.comb) = "Asymptotic"  
      }
      
    }else{ # resampling test + asymptotic test 
      
      tmp = .Score.test2.cluster(id, perm.type, count, Y.base=count.base, X, X.index, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample, test=test)
      pval.zero = c(tmp$zero.results$score.pvalue, tmp$zero.results$score.Rpvalue)
      pval.pos = c(tmp$pos.results$score.pvalue, tmp$pos.results$score.Rpvalue)
      names(pval.zero) = c("Asymptotic", "Resampling")
      names(pval.pos) = c("Asymptotic", "Resampling")
      
      if(test=="chisq"){
        pval.comb = c(tmp$comb.results$score.pvalue, tmp$comb.results$score.Rpvalue)
        names(pval.comb) = c("Asymptotic", "Resampling")
      }
      if(test=="vc"){
        #pval.comb = c(tmp$comb.results$scoreF.pvalue, tmp$comb.results$scoreS.pvalue, tmp$comb.results$scoreF.Rpvalue, tmp$comb.results$scoreS.Rpvalue )
        #names(pval.comb) = c( "Asymptotic.Fisher", "Asymptotic.Simes", "Resampling.Fisher", "Resampling.Simes" )  
        pval.comb = c(tmp$comb.results$score.pvalue, tmp$comb.results$score.Rpvalue)
        names(pval.comb) = c("Asymptotic", "Resampling")        
        
      }
      
    }
    
    return( list(pval=pval.comb, pval.zero=pval.zero, pval.pos=pval.pos) )
    
  }else{ # perform tests for lineages
    
    if(!is.matrix(Tax)){
      warning("Tax table is not a matrix")
      Tax = as.matrix(Tax)  
    }
    
    tax = Tax[keep,,drop=FALSE]
    
    if( sum(colnames(count)!=rownames(tax))>0 ){
      
      stop("Error: OTU IDs in OTU table are not consistent with OTU IDs in Tax table")
    }
    
    W.data = data.table(data.frame(tax, t(count)))
    n.rank = ncol(tax)
    otucols = names(W.data)[-(1:n.rank)]
    if(!is.null(OTU.base) ){
      W.data.base = data.table(data.frame(tax, t(count.base)))
      otucols.base = names(W.data.base)[-(1:n.rank)]
    }
    
    n.level = n.rank-1
    
    subtree = NULL
    pval.comb = NULL
    #pval.combF = NULL
    #pval.combS = NULL
    pval.zero = NULL
    pval.pos = NULL
    
    for(k in 1:n.level){
      
      #print(k)
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)    
      
      tt = W.data[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols, by=list( get(Rank.low), get(Rank.high) )]
      setnames(tt, 1:2, c(Rank.low, Rank.high))
      W.tax = as.vector(unlist(tt[, Rank.low, with=FALSE]))
      W.count = data.matrix(tt[, otucols, with=FALSE])      
      if(!is.null(OTU.base) ){
        tt.base = W.data.base[, lapply(.SD , sum, na.rm=TRUE), .SDcols=otucols.base, by=list( get(Rank.low), get(Rank.high) )]
        W.count.base = data.matrix(tt.base[, otucols.base, with=FALSE])         
      }      
      
      for(j in 1:m.level){
        
        Y = t(W.count[which(W.tax == level.uni[j]), , drop=FALSE])
        
        if(!is.null(OTU.base) ){
          Y.base = t(W.count.base[which(W.tax == level.uni[j]), , drop=FALSE])
        }
        #Y = t(W.count[which(W.tax == "f__Veillonellaceae"), , drop=FALSE])
        
        remove.index = which(colSums(Y)==0)
        
        if(length(remove.index)==ncol(Y)){
          
          #print("==skip:0==");
          next
          
          
        }else{
          
          if(length(remove.index)>0){
            Y = Y[, -remove.index, drop=FALSE] 
            if(!is.null(OTU.base) ){
              Y.base = Y.base[, -remove.index, drop=FALSE]
            }
          }
          if(is.null(OTU.base) ){
            Y.base = NULL
          }
          
          
          if(ncol(Y)==1){
            
            next
            #print("==skip:1==");
            
          }else{
            
            subtree = c(subtree, level.uni[j])
            
            if(is.null(n.resample)){ # asymptotic test only
              
              tmp = .Score.test2.cluster(id, perm.type, Y, Y.base=Y.base, X, X.index, Z, Z.index, seed=11, resample=FALSE, n.replicates=NULL, test=test)   
              
              
              pval.zero = cbind(pval.zero, tmp$zero.results$score.pvalue)
              pval.pos = cbind(pval.pos, tmp$pos.results$score.pvalue)
              if(test=="chisq"){
                pval.comb = cbind(pval.comb, tmp$comb.results$score.pvalue)
              }
              if(test=="vc"){
                pval.comb = cbind(pval.comb, tmp$comb.results$score.pvalue)
                #pval.combF = cbind(pval.combF, tmp$comb.results$scoreF.pvalue)
                #pval.combS = cbind(pval.combS, tmp$comb.results$scoreS.pvalue)
              }                
              
            }else{ # resampling test + asymptotic test 
              
              tmp = .Score.test2.cluster(id, perm.type, Y, Y.base=Y.base, X, X.index, Z, Z.index, seed=11, resample=TRUE, n.replicates=n.resample, test=test)   
              
              pval.zero = cbind(pval.zero, c(tmp$zero.results$score.pvalue, tmp$zero.results$score.Rpvalue) )
              pval.pos = cbind(pval.pos, c(tmp$pos.results$score.pvalue, tmp$pos.results$score.Rpvalue) )
              
              if(test=="chisq"){
                pval.comb = cbind(pval.comb, c(tmp$comb.results$score.pvalue, tmp$comb.results$score.Rpvalue) )
              }
              if(test=="vc"){
                pval.comb = cbind(pval.comb, c(tmp$comb.results$score.pvalue, tmp$comb.results$score.Rpvalue) )
                #pval.combF = cbind(pval.combF, c(tmp$comb.results$scoreF.pvalue, tmp$comb.results$scoreF.Rpvalue) )
                #pval.combS = cbind(pval.combS, c(tmp$comb.results$scoreS.pvalue, tmp$comb.results$scoreS.Rpvalue) )
              }
              
              
            }
            
          }
          
        }
        
        
      }# lineage loop
      
      
    }# level loop
    
    
   
    colnames(pval.zero) = subtree
    colnames(pval.pos) = subtree
    if(test=="chisq"){
      colnames(pval.comb) = subtree
    }
    if(test=="vc"){
      colnames(pval.comb) = subtree
      #colnames(pval.combF) = subtree
      #colnames(pval.combS) = subtree
    }
    
    if(is.null(n.resample)){
      
      rownames(pval.zero) = "Asymptotic"
      score.zero.tmp = pval.zero[1,]
      
      rownames(pval.pos) = "Asymptotic"
      score.pos.tmp = pval.pos[1,]
      
      if(test=="chisq"){
        rownames(pval.comb) = "Asymptotic"
        score.comb.tmp = pval.comb[1,]
      }
      
      if(test=="vc"){
        rownames(pval.comb) = "Asymptotic"
        score.comb.tmp = pval.comb[1,]
        #rownames(pval.combF) = "Asymptotic"
        #score.combF.tmp = pval.combF[1,]
        #rownames(pval.combS) = "Asymptotic"
        #score.combS.tmp = pval.combS[1,]
      }
      
    }else{

      rownames(pval.zero) = c("Asymptotic", "Resampling")
      score.zero.tmp = pval.zero[2,]
      
      rownames(pval.pos) = c("Asymptotic", "Resampling")
      score.pos.tmp = pval.pos[2,]
      
      if(test=="chisq"){
        rownames(pval.comb) = c("Asymptotic", "Resampling")
        score.comb.tmp = pval.comb[2,]
      }
      
      if(test=="vc"){
        rownames(pval.comb) = c("Asymptotic", "Resampling")
        score.comb.tmp = pval.comb[2,]
        #rownames(pval.combF) = c("Asymptotic", "Resampling")
        #score.combF.tmp = pval.combF[2,]
        #rownames(pval.combS) = c("Asymptotic", "Resampling")
        #score.combS.tmp = pval.combS[2,]
      }  
      
    }
    
    #print(pval.comb)
    
    # identify significant lineages
    
    global.pval = NULL
    
    sig.lineage <- vector("list",3)
    names(sig.lineage) <- c("Zero-Part", "Positive-Part", "Two-Part")
    
    #if(test=="chisq"){SS = 3}
    #if(test=="vc"){SS = 4}
    SS = 3
    for(i in 1:SS){
      
      if(i==1){score.tmp = score.zero.tmp}
      if(i==2){score.tmp = score.pos.tmp}
      if(i==3){score.tmp = score.comb.tmp} 
      # if(test=="chisq" & i==3){
      #   score.tmp = score.comb.tmp
      # }
      # if(test=="vc" & i==3){
      #   score.tmp = score.combF.tmp
      # }
      # if(test=="vc" & i==4){
      #   score.tmp = score.combS.tmp
      # }      
      
      subtree.tmp = subtree
      index.na = which(is.na(score.tmp))
      if(length(index.na)>0){
        score.tmp = score.tmp[-index.na]
        subtree.tmp = subtree.tmp[-index.na]
      }
      
      #score.tmp[score.tmp==0] = 1e-4
      m.test = length(score.tmp)
      
      # Benjamini-Hochberg FDR control
      index.p = order(score.tmp)
      p.sort = sort(score.tmp)
      #fdr.alpha = 0.05
      
      # change 04/17/2016
      reject = rep(0, m.test)
      tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
      if(length(tmp)>0){
        index.reject = index.p[1:max(tmp)]
        reject[index.reject] = 1      
      }
      
      sig.lineage[[i]] = subtree.tmp[reject==1]
      
      # perform global test
      global.score.min =  .simes.test(score.tmp)
      global.pval = c(global.pval, global.score.min)
      
      
    }
    
    
    if(test=="chisq"){
      names(global.pval) = c("Simes_Zero-Part", "Simes_Positive-Part", "Simes_Two-Part")
    }
    
    if(test=="vc"){
      #names(global.pval) = c("Simes_Zero-Part", "Simes_Positive-Part", "SimesF_Two-Part", "SimesS_Two-Part")
      names(global.pval) = c("Simes_Zero-Part", "Simes_Positive-Part", "Simes_Two-Part")
    }
    
    pval = list(pval.zero, pval.pos, pval.comb)
    names(pval) = c("Zero-Part", "Positive-Part", "Two-Part")
    return( list(lineage.pval=pval, sig.lineage=sig.lineage, global.pval=global.pval) )
    
    
    
  }
  
}


.Ei.beta.pos.XX <- function(m, p.XX, beta, XX.i, Y.i){
  
  Ei.out = rep(NA,m)
  
  I.i = as.numeric(Y.i>0)
  
  idx1 = 0
  idx2 = 0
  
  for(j in 1:(m-1)){
    
    idx1 = idx2 + 1
    idx2 = idx2 + p.XX[j]
    
    Ei.out[j] = I.i[j] * exp(beta[idx1:idx2] %*% XX.i[,1:p.XX[j],j])
    
  }
  
  
  Ei.out[m] = I.i[m]
  
  
  return (Ei.out)
}


.fun.neg.loglik.beta.pos.XX <- function(beta, data){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  
  #n.beta = (m - 1)*p
  n.beta = sum(p.XX)
  
  loglik = 0
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Y.pos.index = which(Y[i,]>0)
      #loglik = loglik + Y[i,Y.pos.index] %*% log(P.i[Y.pos.index])
      index = which(Y[i,]>0)
      loglik = loglik + Y[i,index] %*% log(P.i[index])
      #       if(is.na(loglik)){
      #         print(i); break;
      #       }
    }
    
  }
  
  return (-loglik)
  
}


.fun.neg.score.beta.pos.XX <- function(beta, data){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  
  #n.beta = (m - 1)*p
  n.beta = sum(p.XX)
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta = rep(0, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      #Score.beta = Score.beta + kronecker(matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1))
      tmp = NULL
      for(j in 1:(m-1)){
        tmp = c(tmp, (Y[i,j] - nY[i]*P.i[j]) * XX[i,1:p.XX[j],j] )
      }
      Score.beta = Score.beta + tmp
      
      
    }
    
    return (-Score.beta)
  }
  
  
  
}


.fun.score.i.beta.pos.XX <- function(beta, data){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  
  #n.beta = (m - 1)*p
  n.beta = sum(p.XX)
  
  if(length(beta)!=n.beta){
    
    warning("Dim of initial beta does not match the dim of covariates")
    
  }else{
    
    Score.beta.i = matrix(0, n, n.beta)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      # add 03/28/2016
      #       if(sum.E.i==0){
      #         P.i = rep(0,m)
      #       }
      
      #Score.beta.i[i,] =  kronecker(matrix(Y[i,-m] - nY[i]*P.i[-m], ncol=1), matrix(X[i,], ncol=1)) 
      tmp = NULL
      for(j in 1:(m-1)){
        tmp = c(tmp, (Y[i,j] - nY[i]*P.i[j]) * XX[i,1:p.XX[j],j] )
      }
      Score.beta.i[i,] = tmp
      
    }
    
    return (Score.beta.i)
  }
  
  
  
}


.fun.hessian.beta.pos.XX <- function(beta, data, save.list=FALSE){
  
  Y = data$Y; XX = data$XX; p.XX = data$p.XX
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(XX)[2]
  #n.beta = (m-1)*p
  n.beta = sum(p.XX)
  
  #idx.end = (1:(m-1))*p 
  #idx.start = idx.end - (p-1)
  
  idx.start = c(0,cumsum(p.XX[-(m-1)]))+1
  idx.end = idx.start + p.XX - 1
  
  if(length(beta)!=n.beta){
    print("Waring: dim of beta is not the same as beta\n")
    
  }else{
    
    Hessian.beta = matrix(0, nrow=n.beta, ncol=n.beta)
    nY = rowSums(Y)
    I.beta.list = list()
    
    for(i in 1:n){
      
      E.i = .Ei.beta.pos.XX(m, p.XX, beta, XX[i,,,drop=FALSE], Y[i,])
      sum.E.i = sum(E.i)
      P.i = E.i/sum.E.i
      
      ## tmp.beta
      #tmp.beta =  (E.i[-m] %o% E.i[-m])*nY[i]/sum.E.i^2 
      tmp.beta =  as.matrix(P.i[-m] %o% P.i[-m])
      diag(tmp.beta) = diag(tmp.beta) - P.i[-m]
      tmp.beta = nY[i] * tmp.beta
      #tmp.beta[is.na(tmp.beta)] = 0  ## add 03/28/2016
      
      #Hessian.beta = Hessian.beta + kronecker(tmp.beta, ( X[i,] %o% X[i,] ))
      Hessian.tmp = matrix(NA, nrow=n.beta, ncol=n.beta)
      for(j in 1:(m-1)){
        for(k in j:(m-1)){
          Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]] = tmp.beta[j,k] * (XX[i,1:p.XX[j],j] %o% XX[i,1:p.XX[k],k])
          Hessian.tmp[idx.start[k]:idx.end[k], idx.start[j]:idx.end[j]] = t(Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]])
        }
      }
      Hessian.beta = Hessian.beta + Hessian.tmp
      
      if(save.list){
        I.beta.list[[i]] = tmp.beta
        
      }
      
    }
    
    
    if(save.list){
      
      return ( list(Hessian.beta=Hessian.beta, I.beta.list = I.beta.list) )
      
    }else{
      
      return (Hessian.beta)
    }
    
  }
  
  
}


.Pi.alpha.ZZ<-function(m, p.ZZ, alpha, ZZ.i){
  
  Pi.out = rep(NA,m)
  
  idx1 = 0
  idx2 = 0
  
  for(j in 1:m){
    
    idx1 = idx2 + 1
    idx2 = idx2 + p.ZZ[j]
    
    tmp = exp(alpha[idx1:idx2] %*% ZZ.i[,1:p.ZZ[j],j])
    if(is.infinite(tmp)){
      Pi.out[j] = 1
    }else{
      Pi.out[j] = tmp/(tmp + 1)
    }
    
  }
  
  
  return (Pi.out)    
}


.fun.score.i.alpha.ZZ <- function(alpha, data, save.list=FALSE){
  
  Y = data$Y; ZZ = data$ZZ; p.ZZ = data$p.ZZ
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(ZZ)[2]
  
  vA.list = list()
  Vinv.list = list()
  VY.list = list()
  
  #n.alpha = m*p
  n.alpha = sum(p.ZZ)
  
  #idx.end = (1:m)*p 
  #idx.start = idx.end - (p-1)
  idx.start = c(0,cumsum(p.ZZ[-m]))+1
  idx.end = idx.start + p.ZZ - 1
  
  if(length(alpha)!=n.alpha){
    
    warning("Dim of initial alpha does not match the dim of covariates")
    
  }else{
    
    Score.alpha.i = matrix(0, n, n.alpha)
    nY = rowSums(Y)
    
    for(i in 1:n){
      
      Pi.i = .Pi.alpha.ZZ(m, p.ZZ, alpha, ZZ[i,,,drop=FALSE])
      vA.tmp = Pi.i*(1-Pi.i)
      A.i = .diag2(vA.tmp) 
      #t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      t.D.i = matrix(0, nrow=n.alpha, ncol=m)
      for(j in 1:m){
        t.D.i[idx.start[j]:idx.end[j], j] = A.i[j,j] * ZZ[i,1:p.ZZ[j],j] 
      }
      
      
      V.i = A.i # independent cor structure
      
      tmp.V.i = ginv(V.i)
      tmp.VY = tmp.V.i %*% (Y[i,] - Pi.i)
      Score.alpha.i[i,] = t.D.i %*% tmp.VY
      
      if(save.list){
        vA.list[[i]] = vA.tmp
        Vinv.list[[i]] = tmp.V.i
        VY.list[[i]] = tmp.VY
      }
    }
    
    
  }
  
  if(save.list){
    
    return ( list(Score.alpha=Score.alpha.i, vA.list = vA.list, Vinv.list = Vinv.list, VY.list=VY.list) )
    
  }else{
    
    return (Score.alpha.i)
  }  
  
}


.fun.hessian.alpha.ZZ <- function(alpha, data){
  
  Y = data$Y; ZZ = data$ZZ; p.ZZ = data$p.ZZ
  
  n = nrow(Y)
  m = ncol(Y)
  #p = dim(ZZ)[2]
  #n.alpha = m*p
  n.alpha = sum(p.ZZ)
  
  #idx.end = (1:m)*p 
  #idx.start = idx.end - (p-1)
  idx.start = c(0,cumsum(p.ZZ[-m]))+1
  idx.end = idx.start + p.ZZ - 1
  
  if(length(alpha)!=n.alpha){
    print("Waring: dim of alpha is not the same as alpha\n")
    
  }else{
    
    Hessian.alpha = matrix(0, nrow=n.alpha, ncol=n.alpha)
    nY = rowSums(Y)
    
    
    for(i in 1:n){
      
      Pi.i = .Pi.alpha.ZZ(m, p.ZZ, alpha, ZZ[i,,,drop=FALSE])
      tmp = Pi.i*(1-Pi.i)
      A.i = .diag2(tmp)
      #t.D.i = kronecker( A.i, as.matrix(Z[i,], ncol=1) )
      t.D.i = matrix(0, nrow=n.alpha, ncol=m)
      for(j in 1:m){
        t.D.i[idx.start[j]:idx.end[j], j] = A.i[j,j] * ZZ[i,1:p.ZZ[j],j] 
      }      
      V.i = A.i # independent cor structure        
      
      
      Hessian.alpha = Hessian.alpha + t.D.i %*% ginv(V.i) %*% t(t.D.i) 
      
      
    }
    
    return (Hessian.alpha)
    
    
  }
  
  
}


########## functions defined for clustered data ##########
  
.Score.test2.cluster <- function(id, perm.type, Y, Y.base=NULL, X, X.par.index, Z, Z.par.index, seed=11, resample=FALSE, n.replicates=NULL, test){
  
  n = nrow(Y)
  remove.index = which(colSums(Y)==0)
  if(length(remove.index)==ncol(Y)){
    
    pos.results = list(score.stat = NA, score.pvalue = NA, df = NA)
    zero.results = list(score.stat = NA , score.pvalue = NA, df = NA)
    if(test=="chisq"){
      comb.results = list(score.stat = NA, score.pvalue = NA, df = NA )
    }
    if(test=="vc"){
      #comb.results = list(score.stat = NA, scoreF.pvalue = NA, scoreS.pvalue = NA, df = NA )
      comb.results = list(score.stat = NA, score.pvalue = NA, df = NA )
    }    
    
    if(resample){
      pos.results = c(pos.results, score.Rpvalue = NA)
      zero.results = c(zero.results, score.Rpvalue = NA)
      if(test=="chisq"){
        comb.results = c(comb.results, score.Rpvalue = NA)
      }
      if(test=="vc"){
        #comb.results = c(comb.results, scoreF.Rpvalue = NA, scoreS.Rpvalue = NA)
        comb.results = c(comb.results, score.Rpvalue = NA)
      }  
      
    }
    
  }else{
    
    if(length(remove.index)>0){
      Y = Y[, -remove.index, drop=FALSE]
      if(!is.null(Y.base)){
        Y.base = Y.base[, -remove.index, drop=FALSE]
        
      }
    }
    
    m = ncol(Y)
    
    if(is.null(Y.base)){
      
      XX = array(NA,dim=c(n,ncol(X),m-1))
      
      
      for(j in 1:(m-1)){
        XX[,,j] = X
        
      }
      
      p.XX = rep(ncol(X), m-1)
      
    }else{
      
      XX = array(NA,dim=c(n,ncol(X)+1,m-1))
      
      p.XX = rep(0, m-1)
      
      for(j in 1:(m-1)){
        
        XX[,1:ncol(X),j] = X
        
        if(length(table(Y.base[,j]))==1){
          p.XX[j] = ncol(X)
        }else{
          p.XX[j] = ncol(X)+1
        }
        
      }
      
      XX[,(ncol(X)+1),] = Y.base[,-m,drop=FALSE]
    }
    
    ############################# Asymptotic: positive part
    index.subj.pos = which(rowSums(Y)>0)
    
    if(length(index.subj.pos)==0){
      
      score.stat.pos = NA
      score.pvalue.pos = NA
      df.pos = NA
      
    }else{
      
      Y1 = Y[index.subj.pos, , drop=FALSE] 
      XX1 = XX[index.subj.pos, , ,drop=FALSE]
      
      # add 05/03/2016 handle exception: what happend if the left subjects have the same values of covariates (e.g., all case/control)
      #index.cova = 1 + which(apply(XX1[,-1,1,drop=FALSE], 2, function(x) length(table(x)) ) > 1) # index of valid covariates
      index.cova.ls = NULL
      for(j in 1:(m-1)){
        index.cova.ls = c(index.cova.ls, list(1 + which(apply(XX1[,-1,j,drop=FALSE], 2, function(x) length(table(x)) ) > 1) ) )# index of valid covariates
      }
      X1.par.index.ava.ls = lapply(index.cova.ls, function(x) as.numeric( !is.na(match(X.par.index, x )) ) ) # use 0/1 to indicate the index is still availiable or not
      # if no covariate left; even if have covariate left, they are not covariate of interest; no taxa
      if( sum(sapply(index.cova.ls, length))<1 |   sum(sapply(X1.par.index.ava.ls, sum))<1 | m<=1 ){
        
        score.stat.pos = NA
        score.pvalue.pos = NA
        df.pos = NA 
        
      }else{
        
        index.cova.ls = lapply(index.cova.ls, function(x) c(1, x))
        #XX1 = XX1[,c(1, index.cova), ,drop=FALSE] 
        #tmp = match(X.par.index, index.cova) + 1
        #X1.par.index = tmp[!is.na(tmp)]
        #d1.pos = length(X1.par.index)
        X1.par.index.ls = lapply(index.cova.ls, function(x) match(X.par.index, x) )
        d1.pos = sum( sapply(X1.par.index.ls, length) )
        
        tmp.pos = try( .Score.test.stat.pos.cluster(id[index.subj.pos], Y1, XX1, index.cova.ls, X1.par.index.ls) )
        #tmp.pos = ( .Score.test.stat.pos.cluster(id[index.subj.pos], Y1, XX1, index.cova.ls, X1.par.index.ls) )
        if("try-error" %in% class(tmp.pos) ){
          
          score.stat.pos = NA
          score.pvalue.pos = NA
          df.pos = NA
          
        }else{
          df.pos = d1.pos 
          
          if(test=="chisq"){
            score.stat.pos = tmp.pos$score.stat.beta
            score.pvalue.pos = 1 - pchisq(score.stat.pos, d1.pos ) 
          }
          
          if(test=="vc"){
            inv.V = ginv(tmp.pos$V)
            score.stat.pos =  tmp.pos$U %*% inv.V %*% inv.V %*% tmp.pos$U
            lambda = eigen(inv.V)$values
            score.pvalue.pos = pchisqsum2.davies(Q=score.stat.pos, lambda=lambda)$p 
            
          }
          
        }
        
      }    
      
    }
    
    
    
    ############################# Asymptotic: zero part
    Y0 = Y
    Y0[Y==0] = 1  
    Y0[Y>0] = 0
    remove.index = which(colSums(Y0)==0)
    
    if(!is.null(Y.base)){
      Y0.base = Y.base
      Y0.base[Y.base==0] = 1  
      Y0.base[Y.base>0] = 0
    }
    # if all 0 in one group across across taxa, then output NA
    #if( (ncol(Y0)-length(remove.index))<=1 | sum(Y0[case==1,])==0 | sum(Y0[case==0,])==0){
    if( ncol(Y0)==length(remove.index) ){
      
      score.stat.zero = NA;
      score.pvalue.zero = NA;
      df.zero = NA    
      
    }else{
      
      if(length(remove.index)>0){
        Y0 = Y0[, -remove.index, drop=FALSE] 
        if(!is.null(Y.base)){
          Y0.base = Y0.base[, -remove.index, drop=FALSE]
          
        }
      }
      m0 = ncol(Y0)
      
      if(is.null(Y.base)){
        
        ZZ = array(NA,dim=c(n,ncol(Z),m0))
        
        
        for(j in 1:m0){
          ZZ[,,j] = Z
          
        }
        
        p.ZZ = rep(ncol(Z), m0)
        
      }else{
        
        ZZ = array(NA,dim=c(n,ncol(Z)+1,m0))
        
        p.ZZ = rep(0, m0)
        
        for(j in 1:m0){
          
          ZZ[,1:ncol(Z),j] = Z
          
          if(length(table(Y0.base[,j]))==1){
            p.ZZ[j] = ncol(Z)
          }else{
            p.ZZ[j] = ncol(Z)+1
          }
          
        }
        
        ZZ[,(ncol(Z)+1),] = Y0.base
      }
      
      
      
      tmp.zero = try( .Score.test.stat.zero.cluster(id, Y0, ZZ, p.ZZ, Z.par.index, "independence") )
      #tmp.zero = ( .Score.test.stat.zero.cluster(id, Y0, ZZ, p.ZZ, Z.par.index, "independence") )
      if("try-error" %in% class(tmp.zero)){
        
        score.stat.zero = NA;
        score.pvalue.zero = NA;
        df.zero = NA
        
      }else{
        
        df.zero = tmp.zero[[1]]
        if(test=="chisq"){
          score.stat.zero = tmp.zero[[2]]
          score.pvalue.zero = tmp.zero[[3]]
        }
        if(test=="vc"){
          ################# 07/19/2018, for the zero part, keep the chisq test
          #inv.V = ginv(tmp.zero$V)
          #score.stat.zero =  tmp.zero$U %*% inv.V %*% inv.V %*% tmp.zero$U
          #lambda = eigen(inv.V)$values
          #score.pvalue.zero = pchisqsum2.davies(Q=score.stat.zero, lambda=lambda)$p 
          score.stat.zero = tmp.zero[[2]]
          score.pvalue.zero = tmp.zero[[3]]
        }        
        
      }
    }
    
    
    ############################# Asymptotic: combined
    if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      df.comb = df.zero + df.pos
      
      if(test=="chisq"){
        score.stat.comb = score.stat.pos + score.stat.zero
        score.pvalue.comb = 1 - pchisq(score.stat.comb, df.comb ) 
      }
      if(test=="vc"){
        
        #score.pvalue.combF =  .F.test(c(score.pvalue.pos, score.pvalue.zero))
        #score.pvalue.combS =  .simes.test(c(score.pvalue.pos, score.pvalue.zero))
        score.pvalue.comb =  .F.test(c(score.pvalue.pos, score.pvalue.zero))
        
      }
      
      
    }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
      
      df.comb = df.pos
      if(test=="chisq"){
        score.stat.comb = score.stat.pos
        score.pvalue.comb = score.pvalue.pos
      }
      if(test=="vc"){
        #score.pvalue.combF = score.pvalue.pos
        #score.pvalue.combS = score.pvalue.pos
        score.pvalue.comb = score.pvalue.pos
      }
      
    }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
      
      df.comb = df.zero
      if(test=="chisq"){
        score.stat.comb = score.stat.zero
        score.pvalue.comb = score.pvalue.zero  
      }
      if(test=="vc"){
        score.pvalue.comb = score.pvalue.zero  
        #score.pvalue.combF = score.pvalue.zero
        #score.pvalue.combS = score.pvalue.zero
      }
      
    }else{
      df.comb = NA
      if(test=="chisq"){
        score.stat.comb = NA
        score.pvalue.comb =  NA
      }
      if(test=="vc"){
        score.pvalue.comb =  NA
        #score.pvalue.combF = NA
        #score.pvalue.combS = NA
      }      
    }
    
    ############################# Resampling if requested
    pos.results = list(score.stat = score.stat.pos, score.pvalue = score.pvalue.pos)
    zero.results = list(score.stat = score.stat.zero, score.pvalue = score.pvalue.zero)
    if(test=="chisq"){
      comb.results = list(score.stat = score.stat.comb, score.pvalue = score.pvalue.comb )
    }
    if(test=="vc"){
      comb.results = list(score.pvalue = score.pvalue.comb )
      #comb.results = list(scoreF.pvalue = score.pvalue.combF, scoreS.pvalue = score.pvalue.combS )
    }
    
    if(resample){
      
      #print("simulated stat:")
      set.seed(seed)
      
      if(!is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        n.zero = 0
        zero.acc = 0
        n.comb = 0
        comb.acc = 0
        
        start.nperm = 1;
        end.nperm = min(100,n.replicates);
        flag = 1
        
        ### 07/19/2018, modify while loop
        while(flag){
          
          results = .resample.work.two.cluster(id, perm.type, XX, X.par.index, X1.par.index.ls, ZZ, p.ZZ, Z.par.index, index.subj.pos, index.cova.ls, score.stat.pos, score.stat.zero, tmp.pos$S.beta.list, tmp.pos$I.beta.list, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, end.nperm, n.pos, pos.acc, n.zero, zero.acc, n.comb, comb.acc, test=test)
          
          n.pos = results$n.pos.new
          pos.acc = results$pos.acc.new
          n.zero = results$n.zero.new
          zero.acc = results$zero.acc.new
          if(test=="chisq"){
            n.comb = results$n.comb.new
            comb.acc = results$comb.acc.new
          }
          flag = results$flag
          next.end.nperm = results$next.end.nperm
          
          start.nperm = end.nperm + 1;
          end.nperm = next.end.nperm;
          
          if(start.nperm >= n.replicates){
            flag = 0
          }else if(end.nperm > n.replicates){ 
            #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))
            results = .resample.work.two.cluster(id, perm.type, XX, X.par.index, X1.par.index.ls, ZZ, p.ZZ, Z.par.index, index.subj.pos, index.cova.ls, score.stat.pos, score.stat.zero, tmp.pos$S.beta.list, tmp.pos$I.beta.list, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, n.replicates, n.pos, pos.acc, n.zero, zero.acc, n.comb, comb.acc, test=test)
            
            n.pos = results$n.pos.new
            pos.acc = results$pos.acc.new
            n.zero = results$n.zero.new
            zero.acc = results$zero.acc.new
            if(test=="chisq"){
              n.comb = results$n.comb.new
              comb.acc = results$comb.acc.new
            }
            
            flag = 0
          }
          
        }
        
        #         if(n.pos<n.replicates/2){
        #           print("#replicate too small for pos test")
        #         }
        #         if(n.zero<n.replicates/2){
        #           print("#replicate too small for zero test")
        #         }
        #         if(n.comb<n.replicates/2){
        #           print("#replicate too small for comb test")
        #         }
        
        score.Rpvalue.pos = (pos.acc+1)/(n.pos+1)
        score.Rpvalue.zero = (zero.acc+1)/(n.zero+1)
        if(test=="chisq"){
          score.Rpvalue.comb = (comb.acc+1)/(n.comb+1)
          comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.comb)
        }
        if(test=="vc"){
          score.Rpvalue.comb = .F.test(c(score.Rpvalue.pos, score.Rpvalue.zero))
          comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.comb)
          #score.Rpvalue.combF = .F.test(c(score.Rpvalue.pos, score.Rpvalue.zero))
          #score.Rpvalue.combS = .simes.test(c(score.Rpvalue.pos, score.Rpvalue.zero))
          #comb.results = c(comb.results, scoreF.Rpvalue = score.Rpvalue.combF, scoreS.Rpvalue = score.Rpvalue.combS)
        }
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        
        
        
      }else if(!is.na(score.stat.pos) & is.na(score.stat.zero)){
        
        n.pos = 0
        pos.acc = 0
        
        start.nperm = 1;
        end.nperm = min(100,n.replicates);
        flag = 1
        
        ### 07/19/2018, modify while loop
        while(flag){
          
          results = .resample.work.pos.cluster(id, perm.type, XX, X.par.index, X1.par.index.ls, index.subj.pos, index.cova.ls, score.stat.pos, tmp.pos$S.beta.list, tmp.pos$I.beta.list, start.nperm, end.nperm, n.pos, pos.acc, test)
          
          n.pos = results$n.pos.new
          pos.acc = results$pos.acc.new
          flag = results$flag
          next.end.nperm = results$next.end.nperm
          
          
          start.nperm = end.nperm + 1;
          end.nperm = next.end.nperm;
          
          if(start.nperm >= n.replicates){
            flag = 0
          }else if(end.nperm > n.replicates){ 
            #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings")) 
            results = .resample.work.pos.cluster(id, perm.type, XX, X.par.index, X1.par.index.ls, index.subj.pos, index.cova.ls, score.stat.pos, tmp.pos$S.beta.list, tmp.pos$I.beta.list, start.nperm, n.replicates, n.pos, pos.acc, test)
            
            n.pos = results$n.pos.new
            pos.acc = results$pos.acc.new
            
            flag = 0
            
          }
          
        }
        
        
        score.Rpvalue.pos = (pos.acc+1)/(n.pos+1)
        
        
        pos.results = c(pos.results, score.Rpvalue = score.Rpvalue.pos)
        zero.results = c(zero.results, score.Rpvalue = NA)
        if(test=="chisq"){
          comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.pos)
        }
        if(test=="vc"){
          #comb.results = c(comb.results, scoreF.Rpvalue = score.Rpvalue.pos, scoreS.Rpvalue = score.Rpvalue.pos)
          comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.pos)
        }        
        
        
      }else if(is.na(score.stat.pos) & !is.na(score.stat.zero)){
        
        n.zero = 0
        zero.acc = 0
        
        start.nperm = 1;
        end.nperm = min(100,n.replicates);
        flag = 1
        
        ### 07/19/2018, modify while loop
        while(flag){
          
          
          results = .resample.work.zero.cluster(id, perm.type, ZZ, p.ZZ, Z.par.index, score.stat.zero, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, end.nperm, n.zero, zero.acc, test)
          
          n.zero = results$n.zero.new
          zero.acc = results$zero.acc.new
          flag = results$flag
          next.end.nperm = results$next.end.nperm
          
          
          start.nperm = end.nperm + 1;
          end.nperm = next.end.nperm;
          
          if(start.nperm >= n.replicates){
            flag = 0
          }else if(end.nperm > n.replicates){ 
            #warning(paste( "Inaccurate pvalue with", n.replicates, "resamplings"))
            
            results = .resample.work.zero.cluster(id, perm.type, ZZ, p.ZZ, Z.par.index, score.stat.zero, tmp.zero$vA.list, tmp.zero$Vinv.list, tmp.zero$VY.list, start.nperm, n.replicates, n.zero, zero.acc, test)
            
            n.zero = results$n.zero.new
            zero.acc = results$zero.acc.new
            
            flag = 0
          }
          
        }
        
        score.Rpvalue.zero = (zero.acc+1)/(n.zero+1)
        
        zero.results = c(zero.results, score.Rpvalue = score.Rpvalue.zero)
        pos.results = c(pos.results, score.Rpvalue = NA)
        if(test=="chisq"){
          comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.zero)
        }
        if(test=="vc"){
          comb.results = c(comb.results, score.Rpvalue = score.Rpvalue.zero)
          #comb.results = c(comb.results, scoreF.Rpvalue = score.Rpvalue.zero, scoreS.Rpvalue = score.Rpvalue.zero)
        }        
        
        
      }else{
        
        pos.results = c(pos.results, score.Rpvalue = NA)
        zero.results = c(zero.results, score.Rpvalue = NA)
        if(test=="chisq"){
          comb.results = c(comb.results, score.Rpvalue = NA)
        }
        if(test=="vc"){
          comb.results = c(comb.results, score.Rpvalue = NA)
          #comb.results = c(comb.results, scoreF.Rpvalue = NA, scoreS.Rpvalue = NA)
        }
      }
      
      
      
    }
    
    
  }
  
  
  return(list(zero.results = zero.results, pos.results = pos.results, comb.results = comb.results))
  
}


.Score.test.stat.pos.cluster <- function(id1, Y1, XX1, index.cova.ls, X1.par.index.ls){
  
  n = nrow(Y1)
  m = ncol(Y1) 
  #p = dim(XX)[2]
  p.XX = sapply(index.cova.ls, length)
  
  XX = array(NA,dim=c(n,dim(XX1)[2],m-1))
  for(j in 1:(m-1)){
    XX[,1:p.XX[j],j] = XX1[,index.cova.ls[[j]],j]
  }
  
  nY = rowSums(Y1)
  
  
  #n.beta = (m - 1)*p
  n.beta = sum(p.XX)
  
  if( sum( sapply(X1.par.index.ls, function(x) sum(x==1)) )>0 ){
    stop("Error: Testing parameters for the intercept is not informative. (Beta part)")
  }  
  
  if(is.null(X1.par.index.ls) || n==0){
    
    score.stat.beta = NA
    
  }else{
    
    
    #XX.reduce = XX[,-X.par.index, , drop=FALSE]
    #p.reduce = p - length(X.par.index)    
    #par.interest.index.beta =  kronecker(((0:(m-2))*p), rep(1,length(X.par.index))) + X.par.index
    
    p.XX.reduce = p.XX - sapply(X1.par.index.ls, length)
    XX.reduce = array(NA,dim=c(n,dim(XX)[2],m-1))
    for(j in 1:(m-1)){
      XX.reduce[,1:p.XX.reduce[j],j] = XX[,(1:p.XX[j])[-X1.par.index.ls[[j]]] ,j]
    }    
    
    par.interest.index.beta = NULL
    tmp = 0
    for(j in 1:(m-1)){
      par.interest.index.beta = c(par.interest.index.beta, tmp + X1.par.index.ls[[j]])
      tmp = sum(p.XX[1:j])
    }
    
    
    n.par.interest.beta = length(par.interest.index.beta) 
    
    #beta.ini.reduce = rep(0, (p.reduce*(m-1)))    
    beta.ini.reduce = rep(0, sum(p.XX.reduce) )   
    
    data.reduce.beta = list(Y=Y1, XX=XX.reduce, p.XX=p.XX.reduce)
    
    est.reduce.beta = rep(NA, n.beta)
    est.reduce.beta[par.interest.index.beta] = 0 
    #est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta.pos, gr=.fun.neg.score.beta.pos, data = data.reduce.beta)$par
    # change 04/08/2016
    est.reduce.beta[-par.interest.index.beta] = optim(par=beta.ini.reduce, fn=.fun.neg.loglik.beta.pos.XX, gr=.fun.neg.score.beta.pos.XX, data = data.reduce.beta, method="BFGS")$par
    
    
    data.beta = list(Y=Y1, XX=XX, p.XX=p.XX)
    Score.reduce.beta = .fun.score.i.beta.pos.XX(est.reduce.beta, data.beta)
    
    # for resampling: S.beta.list, I.beta.list
    #S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, ((1:(m-1))*p-p+1)])
    idx.start = c(0,cumsum(p.XX[-(m-1)]))+1
    S.beta.list = lapply(1:n, function(j) Score.reduce.beta[j, idx.start])
    tmp = .fun.hessian.beta.pos.XX(est.reduce.beta, data.beta, save.list=TRUE)
    I.beta.list = tmp$I.beta.list  
    
    Hess.reduce.beta = tmp$Hessian.beta
    #Hess.reduce.beta =  .fun.hessian.beta.pos(est.reduce.beta, data.beta)
    
    # re-organized the score statistics and Hessian matrix
    Score.reduce.reorg = cbind( matrix(Score.reduce.beta[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
    Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.beta[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                              cbind( matrix(Hess.reduce.beta[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
    
    
    A = colSums(Score.reduce.reorg)[1:n.par.interest.beta]
    
    B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
    
    
    B2 =  matrix(0, n.beta, n.beta)
    
    # 03/12/2018   modify for clustered data
    # 03/28/2018   modify to allow for uneven #obs across clusters
    id1.unique = unique(id1)
    n.cluster = length(id1.unique)
    #TT = n/length(unique(id1))  ## same #observation for each cluster for now
    #n.cluster = n/TT
    #idx1 = 1
    #idx2 = TT
    for(i in 1:n.cluster){
      
      #Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[idx1:idx2,,drop=FALSE] )
      Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[which(id1==id1.unique[i]),,drop=FALSE] )
      B2 = B2 + Score.reduce.reorg.cluster %o%  Score.reduce.reorg.cluster
      
      #idx1 = idx2 + 1
      #idx2 = idx2 + TT
      
    }
    
    
    
    B = B1 %*% B2 %*% t(B1)
    score.stat.beta = A %*% ginv(B) %*% A
    
    
  }
  
  
  return(list(score.stat.beta=score.stat.beta, S.beta.list=S.beta.list, I.beta.list=I.beta.list, U=A, V=B )   )
  
  
}


.Score.test.stat.zero.cluster <- function(id.sam, Y0, ZZ, p.ZZ, Z.par.index, cor.stru){
  
  m = ncol(Y0)
  p.ZZ.reduce = p.ZZ - length(Z.par.index)
  ZZ.reduce = array(NA, dim=c(dim(ZZ)[1], dim(ZZ)[2],m))
  for(j in 1:m){
    ZZ.reduce[,1:p.ZZ.reduce[j],j] = ZZ[,(1:p.ZZ[j])[-Z.par.index],j]
  }
  
  n = nrow(Y0)
  
  #p = dim(ZZ)[2]
  #p.reduce = dim(ZZ.reduce)[2]
  outcome = NULL
  id = NULL
  cova.reduce = NULL
  
  for(i in 1:n){
    
    outcome = c(outcome, Y0[i,])
    
    idx1 = 0
    idx2 = 0
    
    for(j in 1:m){
      
      idx1 = idx2 + 1
      idx2 = idx2 + p.ZZ.reduce[j]
      
      tmp = rep(0, sum(p.ZZ.reduce))
      tmp[idx1:idx2] = ZZ.reduce[i,1:p.ZZ.reduce[j],j]
      cova.reduce = rbind(cova.reduce, tmp )
      
    }
    
    
    id = c(id, rep(i, m))
  }
  
  #data.full = data.frame(outcome=outcome, cova, id = id, row.names=NULL)
  data.reduce = data.frame(outcome=outcome, cova.reduce, id = id, row.names=NULL)
  #gee.full = geeglm(outcome ~ .  - id - 1, data = data.full, id = factor(id), family="binomial", corstr= "independence")
  gee.reduce = geeglm(outcome ~ . - id - 1, data = data.reduce, id = factor(id), family="binomial", corstr= "independence")
  #wald.test = anova(gee.full, gee.reduce)
  
  
  ########### perform score test
  #n.alpha = m *p
  n.alpha = sum(p.ZZ)
  
  #par.interest.index.alpha =  kronecker( ((0:(m-1))*p), rep(1,length(Z.par.index))) + Z.par.index
  par.interest.index.alpha = NULL
  tmp = 0
  for(j in 1:m){
    par.interest.index.alpha = c(par.interest.index.alpha, tmp + Z.par.index)
    tmp = sum(p.ZZ[1:j])
  }
  
  n.par.interest.alpha = length(par.interest.index.alpha) 
  est.reduce.alpha = rep(NA, n.alpha)
  est.reduce.alpha[par.interest.index.alpha] = 0 
  est.reduce.alpha[-par.interest.index.alpha] = coef(gee.reduce)
  est.reduce.scale = gee.reduce
  
  data.alpha = list(Y=Y0, ZZ=ZZ, p.ZZ=p.ZZ)
  
  tmp = .fun.score.i.alpha.ZZ(est.reduce.alpha, data.alpha, save.list=TRUE)
  Score.reduce.alpha = tmp$Score.alpha
  # for resampling test
  vA.list = tmp$vA.list
  Vinv.list = tmp$Vinv.list
  VY.list = tmp$VY.list
  
  Hess.reduce.alpha =  .fun.hessian.alpha.ZZ(est.reduce.alpha, data.alpha)
  # re-organized the score statistics and Hessian matrix
  Score.reduce.reorg = cbind( matrix(Score.reduce.alpha[,par.interest.index.alpha], ncol=n.par.interest.alpha), matrix(Score.reduce.alpha[,-par.interest.index.alpha], ncol=n.alpha - n.par.interest.alpha) )
  Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha[par.interest.index.alpha, par.interest.index.alpha], nrow=n.par.interest.alpha), matrix(Hess.reduce.alpha[par.interest.index.alpha, -par.interest.index.alpha], nrow=n.par.interest.alpha) ), 
                            cbind( matrix(Hess.reduce.alpha[-par.interest.index.alpha, par.interest.index.alpha], nrow=n.alpha - n.par.interest.alpha), matrix(Hess.reduce.alpha[-par.interest.index.alpha, -par.interest.index.alpha], nrow= n.alpha - n.par.interest.alpha)))
  
  
  A = colSums(Score.reduce.reorg)[1:n.par.interest.alpha]
  
  B1 = cbind(diag(n.par.interest.alpha), -Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha+1):n.alpha)] %*% ginv(Hess.reduce.reorg[((n.par.interest.alpha+1):n.alpha), ((n.par.interest.alpha+1):n.alpha)]) )
  
  B2 =  matrix(0, n.alpha, n.alpha)
  # 03/12/2018   modify for clustered data
  # 03/28/2018   modify to allow for uneven #obs across clusters
  id.unique = unique(id.sam)
  n.cluster = length(id.unique)
  #TT = n/length(unique(id.sam))  ## same #observation for each cluster for now
  #n.cluster = n/TT
  #idx1 = 1
  #idx2 = TT
  for(i in 1:n.cluster){
    
    #Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[idx1:idx2,,drop=FALSE] )
    Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[which(id.sam==id.unique[i]),,drop=FALSE] )
    B2 = B2 + Score.reduce.reorg.cluster %o%  Score.reduce.reorg.cluster
    
    #idx1 = idx2 + 1
    #idx2 = idx2 + TT
    
  }
  
  
  
  B = B1 %*% B2 %*% t(B1)
  score.stat.alpha = A %*% ginv(B) %*% A
  score.pvalue.alpha = 1 - pchisq(score.stat.alpha, n.par.interest.alpha) 
  
  
  
  return(list(score.df.alpha=n.par.interest.alpha, score.stat.alpha = score.stat.alpha, score.pvalue.alpha=score.pvalue.alpha, vA.list=vA.list, Vinv.list=Vinv.list, VY.list=VY.list, U=A, V=B )   )
  
}



# add 07/02/2016 for adaptive resampling
.resample.work.two.cluster <- function(id, perm.type, XX, X.par.index, X1.par.index.ls, ZZ, p.ZZ, Z.par.index, index.subj.pos, index.cova.ls, score.stat.pos, score.stat.zero, pos.S.beta.list, pos.I.beta.list, zero.vA.list, zero.Vinv.list, zero.VY.list, start.nperm, end.nperm, n.pos, pos.acc, n.zero, zero.acc, n.comb, comb.acc, test){
  
  n = dim(XX)[1]
  id.unique = unique(id)
  n.cluster = length(id.unique)
  XX.interest = XX[,X.par.index,, drop=FALSE]
  ZZ.interest = ZZ[,Z.par.index,, drop=FALSE]
  
  # 03/28/2018   modify to allow for uneven #obs across clusters
  if(perm.type=="BTW"){
    
    cluster.nobs = NULL
    
    # if covaraites of interest are not the same within cluster, cannot perform between cluster permutation
    for(i in 1:n.cluster){
      tmpX = XX.interest[which(id==id.unique[i]),,,drop=FALSE]
      tmpZ = ZZ.interest[which(id==id.unique[i]),,,drop=FALSE]
      cluster.nobs = c(cluster.nobs, dim(tmpX)[1])
      if(dim(unique(tmpX))[1]!=1 | dim(unique(tmpZ))[1]!=1){
        stop("Cannot perform between cluster permutation: covaraites of interest are not the same within cluster")
      }
      if(i==1){
        XX.interest.unique = unique(tmpX)
        ZZ.interest.unique = unique(tmpZ)
      }else{
        XX.interest.unique = abind(XX.interest.unique, unique(tmpX), along=1)
        ZZ.interest.unique = abind(ZZ.interest.unique, unique(tmpZ), along=1)
      }
      
      
    }
    
    
    
  }
  
  
  
  n.pos.new = n.pos
  pos.acc.new = pos.acc
  
  n.zero.new = n.zero
  zero.acc.new = zero.acc
  
  n.comb.new = n.comb
  comb.acc.new = comb.acc
  
  score.stat.comb = score.stat.pos + score.stat.zero
  
  
  for(k in start.nperm:end.nperm){
    
    # 03/28/2018   modify to allow for uneven #obs across clusters
    if(perm.type=="BTW"){
      
      
      sample.cluster.idx = sample(1:n.cluster)
      
      for(i in 1:n.cluster){
        nobs = cluster.nobs[i]
        XX.interest.perm.i = array(NA, dim=c(nobs, dim(XX.interest.unique)[2], dim(XX.interest.unique)[3]))
        ZZ.interest.perm.i = array(NA, dim=c(nobs, dim(ZZ.interest.unique)[2], dim(ZZ.interest.unique)[3]))
        for(j in 1:dim(XX.interest.unique)[3]){
          XX.interest.perm.i[,,j] =  matrix(rep(XX.interest.unique[sample.cluster.idx[i],,j],nobs),nrow=nobs, byrow = TRUE)
        }
        for(j in 1:dim(ZZ.interest.unique)[3]){
          ZZ.interest.perm.i[,,j] =  matrix(rep(ZZ.interest.unique[sample.cluster.idx[i],,j],nobs),nrow=nobs, byrow = TRUE)
        }
        if(i==1){
          XX.interest.perm = XX.interest.perm.i
          ZZ.interest.perm = ZZ.interest.perm.i
        }else{
          XX.interest.perm = abind(XX.interest.perm, XX.interest.perm.i, along =1)
          ZZ.interest.perm = abind(ZZ.interest.perm, ZZ.interest.perm.i, along =1)
        }
        
      }
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.interest.perm
      ZZ.perm = ZZ
      ZZ.perm[,Z.par.index,] = ZZ.interest.perm      
      
      
    }else if(perm.type=="WTH"){
      
      sample.idx = .sampling.strata.uneven(1:n, id)
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.perm[sample.idx,X.par.index,]
      ZZ.perm = ZZ
      ZZ.perm[,Z.par.index,] = ZZ.perm[sample.idx,Z.par.index,]
      
    }else{
      
      sample.idx = sample(1:n)
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.perm[sample.idx,X.par.index,]
      ZZ.perm = ZZ
      ZZ.perm[,Z.par.index,] = ZZ.perm[sample.idx,Z.par.index,]
    }
    
    #X.perm = X
    #X.perm[,X.par.index] = X.perm[sample.idx,X.par.index]
    
    XX1.perm = XX.perm[index.subj.pos, , ,drop=FALSE]
    #XX1.perm = XX1.perm[,c(1, index.cova), ,drop=FALSE] 
    p.XX = sapply(index.cova.ls, length)
    XX.tmp = array(NA,dim=c(dim(XX1.perm)[1],dim(XX1.perm)[2], dim(XX1.perm)[3] ))
    for(j in 1:length(p.XX)){
      XX.tmp[,1:p.XX[j],j] = XX1.perm[,index.cova.ls[[j]],j]
    }
    XX1.perm = XX.tmp
    
    #Z.perm = Z
    #Z.perm[,Z.par.index] = Z.perm[sample.idx,Z.par.index]
    
    score.stat.pos.perm = try( .Score.test.stat.pos.4Gresampling.cluster(id[index.subj.pos], XX1.perm, index.cova.ls, X1.par.index.ls, pos.S.beta.list,pos.I.beta.list, test) )
    score.stat.zero.perm = try( .Score.test.stat.zero.4Gresampling.cluster(id, ZZ.perm, p.ZZ, Z.par.index, zero.vA.list, zero.Vinv.list, zero.VY.list, test) )
    #score.stat.pos.perm = ( .Score.test.stat.pos.4Gresampling.cluster(id[index.subj.pos], XX1.perm, index.cova.ls, X1.par.index.ls, pos.S.beta.list,pos.I.beta.list) )
    #score.stat.zero.perm = ( .Score.test.stat.zero.4Gresampling.cluster(id, ZZ.perm, p.ZZ, Z.par.index, zero.vA.list, zero.Vinv.list, zero.VY.list) )
    
    
    if(!( "try-error" %in% class(score.stat.pos.perm))){
      
      n.pos.new = n.pos.new + 1
      if(score.stat.pos.perm >= score.stat.pos){
        pos.acc.new = pos.acc.new + 1
        
      } 
    }
    
    if(!( "try-error" %in% class(score.stat.zero.perm))){
      
      n.zero.new = n.zero.new + 1
      if(score.stat.zero.perm >= score.stat.zero){
        zero.acc.new = zero.acc.new + 1
        
      }  
    }
    
    if(test=="chisq"){
      
      if(!( "try-error" %in% class(score.stat.pos.perm)) & !( "try-error" %in% class(score.stat.zero.perm))){
        
        score.stat.comb.perm = score.stat.pos.perm + score.stat.zero.perm
        n.comb.new = n.comb.new + 1
        if(score.stat.comb.perm >= score.stat.comb){
          comb.acc.new = comb.acc.new + 1
          
        }  
      }
    }
    
  }
  
  # modify 07/19/2018
  if(pos.acc.new < 10 | zero.acc.new < 10){
    next.end.nperm = end.nperm * 5;
    flag = 1;
    
  }
  # else if(pos.acc.new<10 | zero.acc.new<10){
  #   next.end.nperm = ( end.nperm + 1) * 10 - 1;
  #   flag = 1;
  #   
  # }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = end.nperm;
    flag = 0;  
  }
  
  return(list(n.pos.new=n.pos.new, pos.acc.new=pos.acc.new, 
              n.zero.new=n.zero.new, zero.acc.new=zero.acc.new, 
              n.comb.new=n.comb.new, comb.acc.new=comb.acc.new,
              flag=flag, next.end.nperm=next.end.nperm))
  
}


# add 07/02/2016 for adaptive resampling
.resample.work.pos.cluster <- function(id, perm.type, XX, X.par.index, X1.par.index.ls, index.subj.pos, index.cova.ls, score.stat.pos, pos.S.beta.list, pos.I.beta.list, start.nperm, end.nperm, n.pos, pos.acc, test){
  
  n = dim(XX)[1]
  id.unique = unique(id)
  n.cluster = length(id.unique)
  XX.interest = XX[,X.par.index,, drop=FALSE]
  
  # 03/28/2018   modify to allow for uneven #obs across clusters
  if(perm.type=="BTW"){
    
    cluster.nobs = NULL
    
    # if covaraites of interest are not the same within cluster, cannot perform between cluster permutation
    for(i in 1:n.cluster){
      tmpX = XX.interest[which(id==id.unique[i]),,,drop=FALSE]
      cluster.nobs = c(cluster.nobs, dim(tmpX)[1])
      if(dim(unique(tmpX))[1]!=1 ){
        stop("Cannot perform between cluster permutation: covaraites of interest are not the same within cluster")
      }
      if(i==1){
        XX.interest.unique = unique(tmpX)
      }else{
        XX.interest.unique = abind(XX.interest.unique, unique(tmpX), along=1)
        
      }
      
      
    }
    
    
    
  }
  
  
  n.pos.new = n.pos
  pos.acc.new = pos.acc
  
  for(k in start.nperm:end.nperm){
    
    if(perm.type=="BTW"){
      
      sample.cluster.idx = sample(1:n.cluster)
      
      for(i in 1:n.cluster){
        nobs = cluster.nobs[i]
        XX.interest.perm.i = array(NA, dim=c(nobs, dim(XX.interest.unique)[2], dim(XX.interest.unique)[3]))
        
        for(j in 1:dim(XX.interest.unique)[3]){
          XX.interest.perm.i[,,j] =  matrix(rep(XX.interest.unique[sample.cluster.idx[i],,j],nobs),nrow=nobs, byrow = TRUE)
          
        }
        if(i==1){
          XX.interest.perm = XX.interest.perm.i
          
        }else{
          XX.interest.perm = abind(XX.interest.perm, XX.interest.perm.i, along =1)
          
        }
        
      }
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.interest.perm
      
      
      
    }else if(perm.type=="WTH"){
      
      sample.idx = .sampling.strata.uneven(1:n, id)
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.perm[sample.idx,X.par.index,]
      
    }else{
      
      sample.idx = sample(1:n)
      
      XX.perm = XX
      XX.perm[,X.par.index,] = XX.perm[sample.idx,X.par.index,]
    }
    
    
    XX1.perm = XX.perm[index.subj.pos, , ,drop=FALSE]
    #XX1.perm = XX1.perm[,c(1, index.cova), ,drop=FALSE] 
    p.XX = sapply(index.cova.ls, length)
    XX.tmp = array(NA,dim=c(dim(XX1.perm)[1],dim(XX1.perm)[2],dim(XX1.perm)[3]))
    for(j in 1:length(p.XX)){
      XX.tmp[,1:p.XX[j],j] = XX1.perm[,index.cova.ls[[j]],j]
    }
    XX1.perm = XX.tmp    
    
    
    score.stat.pos.perm = try( .Score.test.stat.pos.4Gresampling.cluster(id[index.subj.pos], XX1.perm, index.cova.ls, X1.par.index.ls, pos.S.beta.list,pos.I.beta.list, test) )
    #score.stat.pos.perm = ( .Score.test.stat.pos.4Gresampling.cluster(id[index.subj.pos], XX1.perm, index.cova.ls, X1.par.index.ls, pos.S.beta.list,pos.I.beta.list) )
    
    
    if(!( "try-error" %in% class(score.stat.pos.perm))){
      
      n.pos.new = n.pos.new + 1
      if(score.stat.pos.perm >= score.stat.pos){
        pos.acc.new = pos.acc.new + 1
        
      } 
    }
    
    
    
  }
  
  ### modify 07/19/2018
  if(pos.acc.new < 10){
    next.end.nperm = end.nperm * 5;
    flag = 1;
    
  }
  # else if(pos.acc.new<10){
  #   next.end.nperm = ( end.nperm + 1) * 10 - 1;
  #   flag = 1;
  #   
  # }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = end.nperm;
    flag = 0;  
  }
  
  return(list(n.pos.new=n.pos.new, pos.acc.new=pos.acc.new, 
              flag=flag, next.end.nperm=next.end.nperm))
  
}


.Score.test.stat.pos.4Gresampling.cluster <- function(id1, XX1.perm, index.cova.ls, X1.par.index.ls, S.beta.list, I.beta.list, test){
  
  n = dim(XX1.perm)[1]
  #p = dim(XX.perm)[2]
  p.XX = sapply(index.cova.ls, length)
  
  m.beta = length(S.beta.list[[1]])
  
  XX.perm = array(NA,dim=c(n,dim(XX1.perm)[2],m.beta))
  for(j in 1:m.beta){
    XX.perm[,1:p.XX[j],j] = XX1.perm[,index.cova.ls[[j]],j]
  }
  
  
  
  #idx.end = (1:m.beta)*p 
  #idx.start = idx.end - (p-1)
  idx.start = c(0,cumsum(p.XX[-m.beta]))+1
  idx.end = idx.start + p.XX - 1
  
  #n.par.interest.beta = m.beta
  #n.beta = m.beta*p
  n.beta = sum(p.XX)
  #par.interest.index.beta =  kronecker( ((0:(m.beta-1))*p), rep(1,length(X.par.index))) + X.par.index
  par.interest.index.beta = NULL
  tmp = 0
  for(j in 1:m.beta){
    par.interest.index.beta = c(par.interest.index.beta, tmp + X1.par.index.ls[[j]])
    tmp = sum(p.XX[1:j])
  }
  
  n.par.interest.beta = length(par.interest.index.beta) 
  
  
  Score.reduce.beta.perm = matrix(0, n, n.beta )
  Hess.reduce.beta.perm = matrix(0, n.beta, n.beta )
  
  for(i in 1:n){
    
    ###################################################
    #                                                 #
    #         Beta part: resampling Score test        #
    #                                                 #
    ################################################### 
    #Score.reduce.beta.perm[i,] = Score.reduce.beta.perm[i,] + kronecker(matrix(S.beta.list[[i]], ncol=1),  matrix(X.perm[i,], ncol=1))  
    tmp = NULL
    for(j in 1:m.beta){
      tmp = c(tmp, S.beta.list[[i]][j] * XX.perm[i,1:p.XX[j],j] )
    }
    Score.reduce.beta.perm[i,] = tmp
    
    #Hess.reduce.beta.perm = Hess.reduce.beta.perm + kronecker(I.beta.list[[i]], (  X.perm[i,] %o% X.perm[i,] ) )
    
    Hessian.tmp = matrix(NA, nrow=n.beta, ncol=n.beta)
    for(j in 1:m.beta){
      for(k in j:m.beta){
        Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]] = I.beta.list[[i]][j,k] * (XX.perm[i,1:p.XX[j],j] %o% XX.perm[i,1:p.XX[k],k])
        Hessian.tmp[idx.start[k]:idx.end[k], idx.start[j]:idx.end[j]] = t(Hessian.tmp[idx.start[j]:idx.end[j], idx.start[k]:idx.end[k]])
      }
    }
    Hess.reduce.beta.perm = Hess.reduce.beta.perm + Hessian.tmp
    
    
    #     if(sum(is.na(Hess.reduce.beta.perm))>0){
    #       print(i); break;
    #       
    #     }
    
  }
  
  ###################################################
  #                                                 #
  #         Beta part: resampling Score test        #
  #                                                 #
  ################################################### 
  # re-organized the score statistics and Hessian matrix
  Score.reduce.beta.perm.reorg = cbind( matrix(Score.reduce.beta.perm[,par.interest.index.beta], ncol=n.par.interest.beta), matrix(Score.reduce.beta.perm[,-par.interest.index.beta], ncol=n.beta - n.par.interest.beta) )
  Hess.reduce.beta.perm.reorg = rbind(cbind( matrix(Hess.reduce.beta.perm[par.interest.index.beta, par.interest.index.beta], nrow=n.par.interest.beta), matrix(Hess.reduce.beta.perm[par.interest.index.beta, -par.interest.index.beta], nrow=n.par.interest.beta) ), 
                                      cbind( matrix(Hess.reduce.beta.perm[-par.interest.index.beta, par.interest.index.beta], nrow=n.beta - n.par.interest.beta), matrix(Hess.reduce.beta.perm[-par.interest.index.beta, -par.interest.index.beta], nrow= n.beta - n.par.interest.beta)))
  
  
  A = colSums(Score.reduce.beta.perm.reorg)[1:n.par.interest.beta]
  
  B1 = cbind(diag(n.par.interest.beta), -Hess.reduce.beta.perm.reorg[(1:n.par.interest.beta), ((n.par.interest.beta+1):n.beta)] %*% ginv(Hess.reduce.beta.perm.reorg[((n.par.interest.beta+1):n.beta), ((n.par.interest.beta+1):n.beta)]) )
  
  
  B2 =  matrix(0, n.beta, n.beta)
  
  # change in 04/08/2016
  # 03/28/2018   modify to allow for uneven #obs across clusters
  id1.unique = unique(id1)
  n.cluster = length(id1.unique)
  #TT = n/length(unique(id))  ## same #observation for each cluster for now
  #n.cluster = n/TT
  #idx1 = 1
  #idx2 = TT
  for(i in 1:n.cluster){
    
    #Score.reduce.beta.perm.reorg.cluster =  colSums( Score.reduce.beta.perm.reorg[idx1:idx2,,drop=FALSE] )
    Score.reduce.beta.perm.reorg.cluster =  colSums( Score.reduce.beta.perm.reorg[which(id1==id1.unique[i]),,drop=FALSE] )
    B2 = B2 + Score.reduce.beta.perm.reorg.cluster %o%  Score.reduce.beta.perm.reorg.cluster
    
    #idx1 = idx2 + 1
    #idx2 = idx2 + TT
    
  }
  
  
  B = B1 %*% B2 %*% t(B1)
  
  if(test=="chisq"){
    score.stat.beta.perm = A %*% ginv(B) %*% A
  }
  
  if(test=="vc"){
    
    inv.B = ginv(B)
    score.stat.beta.perm = A %*% inv.B %*% inv.B %*% A
  }
  
  return(score.stat.beta.perm)
  
  
  
}


# add 07/02/2016 for adaptive resampling
.resample.work.zero.cluster <- function(id, perm.type, ZZ, p.ZZ, Z.par.index, score.stat.zero, zero.vA.list, zero.Vinv.list, zero.VY.list, start.nperm, end.nperm, n.zero, zero.acc, test){
  
  n = dim(ZZ)[1]
  id.unique = unique(id)
  n.cluster = length(id.unique)
  ZZ.interest = ZZ[,Z.par.index,, drop=FALSE]
  
  # 03/28/2018   modify to allow for uneven #obs across clusters
  if(perm.type=="BTW"){
    
    cluster.nobs = NULL
    
    # if covaraites of interest are not the same within cluster, cannot perform between cluster permutation
    for(i in 1:n.cluster){
      
      tmpZ = ZZ.interest[which(id==id.unique[i]),,,drop=FALSE]
      cluster.nobs = c(cluster.nobs, dim(tmpZ)[1])
      if(dim(unique(tmpZ))[1]!=1){
        stop("Cannot perform between cluster permutation: covaraites of interest are not the same within cluster")
      }
      if(i==1){
        
        ZZ.interest.unique = unique(tmpZ)
      }else{
        
        ZZ.interest.unique = abind(ZZ.interest.unique, unique(tmpZ), along=1)
      }
      
      
    }
    
    
  }
  
  n.zero.new = n.zero
  zero.acc.new = zero.acc
  
  for(k in start.nperm:end.nperm){
    
    if(perm.type=="BTW"){
      
      
      sample.cluster.idx = sample(1:n.cluster)
      
      for(i in 1:n.cluster){
        nobs = cluster.nobs[i]
        ZZ.interest.perm.i = array(NA, dim=c(nobs, dim(ZZ.interest.unique)[2], dim(ZZ.interest.unique)[3]))
        for(j in 1:dim(ZZ.interest.unique)[3]){
          ZZ.interest.perm.i[,,j] =  matrix(rep(ZZ.interest.unique[sample.cluster.idx[i],,j],nobs),nrow=nobs, byrow = TRUE)
        }
        if(i==1){
          ZZ.interest.perm = ZZ.interest.perm.i
        }else{
          ZZ.interest.perm = abind(ZZ.interest.perm, ZZ.interest.perm.i, along =1)
        }
        
      }
      
      ZZ.perm = ZZ
      ZZ.perm[,Z.par.index,] = ZZ.interest.perm   
      
      
    }else if(perm.type=="WTH"){
      
      sample.idx = .sampling.strata.uneven(1:n, id)
      
      ZZ.perm = ZZ
      ZZ.perm[,Z.par.index,] = ZZ.perm[sample.idx,Z.par.index,]
      
    }else{
      
      sample.idx = sample(1:n)
      
      ZZ.perm = ZZ
      ZZ.perm[,Z.par.index,] = ZZ.perm[sample.idx,Z.par.index,]
    }
    
    
    
    score.stat.zero.perm = try( .Score.test.stat.zero.4Gresampling.cluster(id, ZZ.perm, p.ZZ, Z.par.index, zero.vA.list, zero.Vinv.list, zero.VY.list, test) )
    #score.stat.zero.perm = ( .Score.test.stat.zero.4Gresampling.cluster(id, ZZ.perm, p.ZZ, Z.par.index, zero.vA.list, zero.Vinv.list, zero.VY.list) )
    
    
    
    if(!( "try-error" %in% class(score.stat.zero.perm))){
      
      n.zero.new = n.zero.new + 1
      if(score.stat.zero.perm >= score.stat.zero){
        zero.acc.new = zero.acc.new + 1
        
      }  
    }
    
    
  }
  
  ### modify 07/19/2018
  if(zero.acc.new < 10){
    next.end.nperm = end.nperm * 5;
    flag = 1;
    
  }
  # else if(zero.acc.new<10){
  #   next.end.nperm = ( end.nperm + 1) * 10 - 1;
  #   flag = 1;
  #   
  # }
  #   else if(one.acc.new<20){
  #     next.end.nperm = ( end.nperm + 1) * 5 - 1;
  #     flag = 1;
  #     
  #   }
  else{
    next.end.nperm = end.nperm;
    flag = 0;  
  }
  
  return(list(n.zero.new=n.zero.new, zero.acc.new=zero.acc.new, 
              flag=flag, next.end.nperm=next.end.nperm))
  
}


.Score.test.stat.zero.4Gresampling.cluster <- function(id, ZZ.perm, p.ZZ, Z.par.index, vA.list, Vinv.list, VY.list, test){
  
  n = dim(ZZ.perm)[1]
  #p = dim(ZZ.perm)[2]
  
  
  m.alpha = length(vA.list[[1]])
  #n.alpha = m.alpha*p
  n.alpha = sum(p.ZZ)
  
  #par.interest.index.alpha =  kronecker( ((0:(m.alpha-1))*p), rep(1,length(Z.par.index))) + Z.par.index
  par.interest.index.alpha = NULL
  tmp = 0
  for(j in 1:m.alpha){
    par.interest.index.alpha = c(par.interest.index.alpha, tmp + Z.par.index)
    tmp = sum(p.ZZ[1:j])
  }
  
  
  n.par.interest.alpha = length(par.interest.index.alpha) 
  
  Score.reduce.alpha.perm = matrix(0, n, n.alpha )
  Hess.reduce.alpha.perm = matrix(0, n.alpha, n.alpha )
  
  #idx.end = (1:m.alpha)*p 
  #idx.start = idx.end - (p-1)
  idx.start = c(0,cumsum(p.ZZ[-m.alpha]))+1
  idx.end = idx.start + p.ZZ - 1
  
  for(i in 1:n){
    
    ###################################################
    #                                                 #
    #         alpha part: resampling Score test        #
    #                                                 #
    ################################################### 
    #tD.tmp = kronecker(.diag2(vA.list[[i]]), as.matrix(Z.perm[i,], ncol=1))
    tD.tmp = matrix(0, nrow=n.alpha, ncol=m.alpha)
    for(j in 1:m.alpha){
      tD.tmp[idx.start[j]:idx.end[j], j] = vA.list[[i]][j] * ZZ.perm[i,1:p.ZZ[j],j] 
    }
    
    Score.reduce.alpha.perm[i,] = Score.reduce.alpha.perm[i,] + tD.tmp %*% VY.list[[i]]
    
    Hess.reduce.alpha.perm = Hess.reduce.alpha.perm + tD.tmp %*% Vinv.list[[i]] %*% t(tD.tmp)
    
    
  }
  
  # re-organized the score statistics and Hessian matrix
  Score.reduce.reorg = cbind( matrix(Score.reduce.alpha.perm[,par.interest.index.alpha], ncol=n.par.interest.alpha), matrix(Score.reduce.alpha.perm[,-par.interest.index.alpha], ncol=n.alpha - n.par.interest.alpha) )
  Hess.reduce.reorg = rbind(cbind( matrix(Hess.reduce.alpha.perm[par.interest.index.alpha, par.interest.index.alpha], nrow=n.par.interest.alpha), matrix(Hess.reduce.alpha.perm[par.interest.index.alpha, -par.interest.index.alpha], nrow=n.par.interest.alpha) ), 
                            cbind( matrix(Hess.reduce.alpha.perm[-par.interest.index.alpha, par.interest.index.alpha], nrow=n.alpha - n.par.interest.alpha), matrix(Hess.reduce.alpha.perm[-par.interest.index.alpha, -par.interest.index.alpha], nrow= n.alpha - n.par.interest.alpha)))
  
  
  A = colSums(Score.reduce.reorg)[1:n.par.interest.alpha]
  
  B1 = cbind(diag(n.par.interest.alpha), -Hess.reduce.reorg[(1:n.par.interest.alpha), ((n.par.interest.alpha+1):n.alpha)] %*% ginv(Hess.reduce.reorg[((n.par.interest.alpha+1):n.alpha), ((n.par.interest.alpha+1):n.alpha)]) )
  
  B2 =  matrix(0, n.alpha, n.alpha)
  # 03/12/2018   modify for clustered data
  # 03/28/2018   modify to allow for uneven #obs across clusters
  id.unique = unique(id)
  n.cluster = length(id.unique)
  #TT = n/length(unique(id))  ## same #observation for each cluster for now
  #n.cluster = n/TT
  #idx1 = 1
  #idx2 = TT
  for(i in 1:n.cluster){
    
    #Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[idx1:idx2,,drop=FALSE] )
    Score.reduce.reorg.cluster =  colSums( Score.reduce.reorg[which(id==id.unique[i]),,drop=FALSE] )
    B2 = B2 + Score.reduce.reorg.cluster %o%  Score.reduce.reorg.cluster
    
    #idx1 = idx2 + 1
    #idx2 = idx2 + TT
    
  }
  
  B = B1 %*% B2 %*% t(B1)
  
  if(test=="chisq"){
    score.stat.alpha.perm = A %*% ginv(B) %*% A
  }
  
  if(test=="vc"){
    
    ##################### 07/19/2018. for zero-part, still use chisq test
    #inv.B = ginv(B)
    #score.stat.alpha.perm = A %*% inv.B %*% inv.B %*% A
    
    score.stat.alpha.perm = A %*% ginv(B) %*% A
  }
  
  
  return(score.stat.alpha.perm)
  
}

