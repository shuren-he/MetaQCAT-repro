########################################
#                                      #
#             ZIGDM test               #
#             add in v2.1              #
#             10/17/2017               #
#                                      #
########################################

ZIGDM <- function(OTU, X4mean=NULL, X4disp=NULL, X4zero=NULL, test.type="Mean", X.index, ZI.LB=10, Tax=NULL, min.depth=0, n.perm=NULL, fdr.alpha=0.05, details=FALSE){
  
  n.resample = n.perm
  
  if(!is.matrix(OTU)){
    warning("OTU table is not a matrix")
    OTU = as.matrix(OTU)  
  }
  otu.name = colnames(OTU)
  
  zero.par.name = "Intercept"; mean.par.name = "Intercept"; disp.par.name = "Intercept"; 
  
  if(!is.null(X4zero) & !is.matrix(X4zero)){
    warning("Covariate table X4zero is not a matrix")
    X4zero = as.matrix(X4zero)  
  }  
  var.name = colnames(X4zero)
  if(!is.null(X4zero)){
    if(is.null(var.name)){
      zero.par.name = c(zero.par.name, paste("var", 1:ncol(X4zero), sep="") )
    }else{
      zero.par.name = c(zero.par.name, var.name)
    }
  }
  
  if(!is.null(X4mean) & !is.matrix(X4mean)){
    warning("Covariate table X4mean is not a matrix")
    X4mean = as.matrix(X4mean)  
  }  
  var.name = colnames(X4mean)
  if(!is.null(X4mean)){
    if(is.null(var.name)){
      mean.par.name = c(mean.par.name, paste("var", 1:ncol(X4mean), sep="") )
    }else{
      mean.par.name = c(mean.par.name, var.name)
    }
  }
  
  if(!is.null(X4disp) & !is.matrix(X4disp)){
    warning("Covariate table X4disp is not a matrix")
    X4disp = as.matrix(X4disp)  
  }  
  var.name = colnames(X4disp)
  if(!is.null(X4disp)){
    if(is.null(var.name)){
      disp.par.name = c(disp.par.name, paste("var", 1:ncol(X4disp), sep="") )
    }else{
      disp.par.name = c(disp.par.name, var.name)
    }
  }
  
  
  flag = 0
  if( !is.null(X4zero) ){
    if(nrow(OTU)!=nrow(X4zero)){
      flag = 1
    }
  }
  if( !is.null(X4mean) ){
    if(nrow(OTU)!=nrow(X4mean)){
      flag = 1
    }
  }
  if( !is.null(X4disp) ){
    if(nrow(OTU)!=nrow(X4disp)){
      flag = 1
    }
  }
  
  if(flag){
    stop("Samples in the OTU table and the covariate tables should be the same.")  
  }
  
  if( test.type == "Zero" & is.null(ZI.LB)){
    stop("Test for differential presence-absence frequency cannot be performed under GDM model (ZI.LB is NULL).")  
  }
  
  if( (test.type == "Mean" & is.null(X4mean)) | (test.type == "Disp" & is.null(X4disp)) | (test.type == "Zero" & is.null(X4zero)) | (test.type == "Omni" & (is.null(X4mean) | is.null(X4disp) | is.null(X4zero)) ) ){
    stop("The specified test cannot be performed because the required covariate matrix is not provided.")
  }
  
  if( test.type == "Omni" & (!identical(X4mean, X4disp) | !identical(X4mean, X4zero) ) ){
    stop("X4mean, X4disp, and X4zero should be the same in order to perform Omnibus test.")
  }
  
  if(missing(X.index)){
    if(test.type == "Zero"){
      X.index = 1:ncol(X4zero)
    }
    if(test.type == "Mean"){
      X.index = 1:ncol(X4mean)
    }    
    if(test.type == "Disp"){
      X.index = 1:ncol(X4disp)
    }
  }  
  
  remove.subject = which(rowSums(OTU)<min.depth)
  if(length(remove.subject)>0){
    
    print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth ))
    X4zero = X4zero[-remove.subject,  ,drop=FALSE]
    X4mean = X4mean[-remove.subject,  ,drop=FALSE]
    X4disp = X4disp[-remove.subject,  ,drop=FALSE]
    OTU = OTU[-remove.subject, ,drop=FALSE]
  }
  
  keep = which(colSums(OTU)>0)
  count = OTU[,keep,drop=FALSE]
  
  n = nrow(count)
  
  if( is.null(X4zero) ){
    X4zero = matrix(1, nrow=n, 1)
  }else{
    X4zero = cbind(1, X4zero) # add the intercept term
  }
  
  if( is.null(X4mean) ){
    X4mean = matrix(1, nrow=n, 1)
  }else{
    X4mean = cbind(1, X4mean) # add the intercept term
  }
  
  if( is.null(X4disp) ){
    X4disp = matrix(1, nrow=n, 1)
  }else{
    X4disp = cbind(1, X4disp) # add the intercept term
  }
  
  
  X.index = X.index + 1
  
  test.details = NULL
  
  if(is.null(Tax)){ # perform one test using all OTUs
    
    
    aa = order( colSums(count), decreasing = TRUE )
    Y = count[,c(aa[-1],aa[1])]
    otu.name = otu.name[c(aa[-1],aa[1])]
    
    if(is.null(ZI.LB)){
      zindex = NULL
    }else{
      zindex = which(colSums(count[,1:(ncol(count)-1),drop=FALSE]==0)>ZI.LB)
      if(length(zindex)==0){
        zindex = NULL
      } 
    }
    
    if(is.null(n.resample)){ # asymptotic test only
      
      if(test.type == "Zero"){
        tmp = .ZIGDM.freq.test.PAR2(Y, X4zero, X4mean, X4disp, X.index, zero.index=zindex, nperm=NULL)
        pval = as.matrix( tmp$score.pval.freq )
      }
      if(test.type == "Mean"){
        tmp = .ZIGDM.abun.test.PAR2(Y, X4zero, X4mean, X4disp, X.index, zero.index=zindex, nperm=NULL)
        pval = as.matrix( tmp$score.pval.abun )
      }      
      if(test.type == "Disp"){
        tmp = .ZIGDM.disp.test.PAR2(Y, X4zero, X4mean, X4disp, X.index, zero.index=zindex, nperm=NULL)
        pval = as.matrix( tmp$score.pval.disp )
      } 
      if(test.type == "Omni"){
        tmp = .ZIGDM.omni.test.PAR2(Y, X4mean, X.index, zero.index=zindex, nperm=NULL)
        pval = as.matrix( tmp$score.pval.omni )
      }
      
      colnames(pval) = "Asymptotic"
      
      
    }else{ # resampling test + asymptotic test 
      
      if(test.type == "Zero"){
        tmp = .ZIGDM.freq.test.PAR2(Y, X4zero, X4mean, X4disp, X.index, zero.index=zindex, nperm=n.resample)
        pval = c( tmp$score.pval.freq, tmp$score.pval.perm.freq )
      }
      if(test.type == "Mean"){
        tmp = .ZIGDM.abun.test.PAR2(Y, X4zero, X4mean, X4disp, X.index, zero.index=zindex, nperm=n.resample)
        pval = c( tmp$score.pval.abun, tmp$score.pval.perm.abun )
      }      
      if(test.type == "Disp"){
        tmp = .ZIGDM.disp.test.PAR2(Y, X4zero, X4mean, X4disp, X.index, zero.index=zindex, nperm=n.resample)
        pval = c( tmp$score.pval.disp, tmp$score.pval.perm.disp )
      } 
      if(test.type == "Omni"){
        tmp = .ZIGDM.omni.test.PAR2(Y, X4mean, X.index, zero.index=zindex, nperm=n.resample)
        pval = c( tmp$score.pval.omni, tmp$score.pval.perm.omni )
      }
      names(pval) = c("Asymptotic", "Resampling")
      
    }
    
    test.details = list(nuis.est = tmp$nuis.est, S1 = tmp$score.S1, S1cov = tmp$score.S1cov, test.stat = tmp$score.stat )
    rownames(test.details$nuis.est$zero.par.est) = otu.name[-length(otu.name)]
    rownames(test.details$nuis.est$mean.par.est) = otu.name[-length(otu.name)]
    rownames(test.details$nuis.est$disp.par.est) = otu.name[-length(otu.name)]
    if(test.type == "Zero"){
      colnames(test.details$nuis.est$zero.par.est) = zero.par.name[-X.index]
      colnames(test.details$nuis.est$mean.par.est) = mean.par.name
      colnames(test.details$nuis.est$disp.par.est) = disp.par.name
    }
    if(test.type == "Mean"){
      colnames(test.details$nuis.est$zero.par.est) = zero.par.name
      colnames(test.details$nuis.est$mean.par.est) = mean.par.name[-X.index]
      colnames(test.details$nuis.est$disp.par.est) = disp.par.name
    }  
    if(test.type == "Disp"){
      colnames(test.details$nuis.est$zero.par.est) = zero.par.name
      colnames(test.details$nuis.est$mean.par.est) = mean.par.name
      colnames(test.details$nuis.est$disp.par.est) = disp.par.name[-X.index] 
    }  
    
    if(details){
      return( list(pval=pval, details=test.details) )
    }else{
      return( list(pval=pval) )
    }
    
    
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
      
      
      for(j in 1:m.level){
        
        Y = t(W.count[which(W.tax == level.uni[j]), , drop=FALSE])
        Y.child = as.character( unlist(tt[which(W.tax == level.uni[j]),2]) )
        #Y = t(W.count[which(W.tax == "f__Veillonellaceae"), , drop=FALSE])
        
        remove.index = which(colSums(Y)==0)
        
        if(length(remove.index)==ncol(Y)){
          
          #print("==skip:0==");
          next
          
          
        }else{
          
          if(length(remove.index)>0){
            Y = Y[, -remove.index, drop=FALSE] 
            Y.child = Y.child[-remove.index]
          }
          
          
          if(ncol(Y)==1){
            
            next
            #print("==skip:1==");
            
          }else{
            
            aa = order( colSums(Y), decreasing = TRUE )
            Y2 = Y[,c(aa[-1],aa[1])]
            Y2.child = Y.child[c(aa[-1],aa[1])]
            
            if(is.null(ZI.LB)){
              zindex = NULL
            }else{
              zindex = which(colSums(Y2[,1:(ncol(Y2)-1),drop=FALSE]==0)>ZI.LB)
              if(length(zindex)==0){
                zindex = NULL
              } 
            }
            
            # remove subject with depth ==0
            remove.sub.index = which(rowSums(Y2)==0)
            
            if(length(remove.sub.index)==0){
              X2.freq = X4zero
              X2.mean = X4mean
              X2.disp = X4disp
              
            }else{
              
              Y2 = Y2[-remove.sub.index, ,drop=FALSE]
              X2.freq = X4zero[-remove.sub.index, ,drop=FALSE]
              X2.mean = X4mean[-remove.sub.index, ,drop=FALSE]
              X2.disp = X4disp[-remove.sub.index, ,drop=FALSE]
              
            }
            
            subtree = c(subtree, level.uni[j])
            
            if(is.null(n.resample)){ # asymptotic test only
              
              if(test.type == "Zero"){
                tmp = .ZIGDM.freq.test.PAR2(Y2, X2.freq, X2.mean, X2.disp, X.index, zero.index=zindex, nperm=NULL)
                pval = cbind(pval, tmp$score.pval.freq )
              }
              if(test.type == "Mean"){
                tmp = .ZIGDM.abun.test.PAR2(Y2, X2.freq, X2.mean, X2.disp, X.index, zero.index=zindex, nperm=NULL)
                pval = cbind(pval, tmp$score.pval.abun )
              }      
              if(test.type == "Disp"){
                tmp = .ZIGDM.disp.test.PAR2(Y2, X2.freq, X2.mean, X2.disp, X.index, zero.index=zindex, nperm=NULL)
                pval = cbind(pval, tmp$score.pval.disp )
              } 
              if(test.type == "Omni"){
                tmp = .ZIGDM.omni.test.PAR2(Y2, X2.mean, X.index, zero.index=zindex, nperm=NULL)
                pval = cbind(pval, tmp$score.pval.omni )
              }
              
              
              
            }else{ # resampling test + asymptotic test 
              
              if(test.type == "Zero"){
                tmp = .ZIGDM.freq.test.PAR2(Y2, X2.freq, X2.mean, X2.disp, X.index, zero.index=zindex, nperm=n.resample)
                pval = cbind(pval, c(tmp$score.pval.freq, tmp$score.pval.perm.freq) )
              }
              if(test.type == "Mean"){
                tmp = .ZIGDM.abun.test.PAR2(Y2, X2.freq, X2.mean, X2.disp, X.index, zero.index=zindex, nperm=n.resample)
                pval = cbind(pval, c(tmp$score.pval.abun, tmp$score.pval.perm.abun) )
              }      
              if(test.type == "Disp"){
                tmp = .ZIGDM.disp.test.PAR2(Y2, X2.freq, X2.mean, X2.disp, X.index, zero.index=zindex, nperm=n.resample)
                pval = cbind(pval, c(tmp$score.pval.disp, tmp$score.pval.perm.disp) )
              } 
              if(test.type == "Omni"){
                tmp = .ZIGDM.omni.test.PAR2(Y2, X2.mean, X.index, zero.index=zindex, nperm=n.resample)
                pval = cbind(pval, c(tmp$score.pval.omni, tmp$score.pval.perm.omni) )
              }
              
              
            }
            
            test.details.tmp = list(nuis.est = tmp$nuis.est, S1 = tmp$score.S1, S1cov = tmp$score.S1cov, test.stat = tmp$score.stat )
            rownames(test.details.tmp$nuis.est$zero.par.est) = Y2.child[-length(Y2.child)]
            rownames(test.details.tmp$nuis.est$mean.par.est) = Y2.child[-length(Y2.child)]
            rownames(test.details.tmp$nuis.est$disp.par.est) = Y2.child[-length(Y2.child)]
            if(test.type == "Zero"){
              colnames(test.details.tmp$nuis.est$zero.par.est) = zero.par.name[-X.index]
              colnames(test.details.tmp$nuis.est$mean.par.est) = mean.par.name
              colnames(test.details.tmp$nuis.est$disp.par.est) = disp.par.name
            }
            if(test.type == "Mean"){
              colnames(test.details.tmp$nuis.est$zero.par.est) = zero.par.name
              colnames(test.details.tmp$nuis.est$mean.par.est) = mean.par.name[-X.index]
              colnames(test.details.tmp$nuis.est$disp.par.est) = disp.par.name
            }  
            if(test.type == "Disp"){
              colnames(test.details.tmp$nuis.est$zero.par.est) = zero.par.name
              colnames(test.details.tmp$nuis.est$mean.par.est) = mean.par.name
              colnames(test.details.tmp$nuis.est$disp.par.est) = disp.par.name[-X.index] 
            } 
            
            test.details[[level.uni[j]]] = test.details.tmp
            
            
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
    
    if(details){
      return( list(lineage.pval=pval, sig.lineage=sig.lineage, global.pval=global.pval, details=test.details) )
    }else{
      return( list(lineage.pval=pval, sig.lineage=sig.lineage, global.pval=global.pval) )
    }
    
    
    
  }
  
}



.lbeta2 <- function(a,b){
  a[a<0] = 0
  b[b<0] = 0
  return(lbeta(a,b))
}

# parameters of posterior probability P_1...P_K | Y_1...Y_K Y_{K+1}
.rP.Y.par <- function(pv, av, bv, Y){
  
  K = length(Y)-1
  N = sum(Y)
  
  pv.post = rep(0, K)
  
  av.prim = av + Y[1:K]  
  bv.prim = bv
  for(j in 1:K){
    bv.prim[j] = bv.prim[j] + (N - sum(Y[1:j]))
    
    if(Y[j]==0){
      
      if(beta(av[j], bv[j])==0){
        
        pv.post[j] = pv[j]/(pv[j] + (1-pv[j]))
        
      }else{
        tmp = pv[j] + (1-pv[j])*beta(av.prim[j], bv.prim[j])/beta(av[j], bv[j])
        if(tmp==0){
          pv.post[j] = 1
        }else{
          pv.post[j] = pv[j]/tmp
        }
        
      }
      
      
    }
    
  }
  
  return(list(pv.post = pv.post, av.post = av.prim, bv.post = bv.prim))
}

# maximize the objective function for logistic regression
# par = Logit.par.ini
# data = Logit.data
.Logit.neg.loglik <- function(par, data){
  
  tmp = as.numeric( exp(  data$X %*% par ) )
  p = tmp/(1+tmp)
  tmp = data$Del * log(p) + (1-data$Del) * log(1-p)
  index = which(p==0 | p==1)
  if(length(index)>0){
    tmp[index] = 0
  }
  #tmp[is.na(tmp)] = 0   ## remove in 12/28/2016 will make optim function output a large estimate  
  return( -sum( tmp) )
  
}

.Logit.neg.score <- function(par, data){
  
  tmp = as.numeric( exp(  data$X %*% par ) )
  p = tmp/(1+tmp)
  return( -colSums( (data$Del - p) * data$X ) )
  
}


#.Logit.optim(Del.R[,j], X.null, c.last[j, ])
#X=X.null
#Del = Del.R[,j]
#c.ini = c.last[j,]
.Logit.optim <- function(Del, X, c.ini){
  
  Logit.par.ini = c.ini
  Logit.data = list(Del=Del, X=X)
  
  #.Logit.neg.loglik(Logit.par.ini, Logit.data)
  #.Logit.neg.score(Logit.par.ini, Logit.data)
  return( optim(par=Logit.par.ini, fn=.Logit.neg.loglik, gr=.Logit.neg.score, data = Logit.data, method="BFGS")$par)
  #return( optim(par=Logit.par.ini, fn=.Logit.neg.loglik, gr=.Logit.neg.score, data = Logit.data)$par)
  
}


# maximize the objective function for Beta regression
# par = Beta.par.ini
# data = Beta.data
.Beta.neg.loglik.PAR2 <- function(par, data){
  
  da = ncol(data$Xa)
  alpha = par[1:da]
  beta = par[-(1:da)]
  
  
  tmp = as.numeric( exp( data$Xa %*% alpha) )
  mu.tmp = as.numeric( tmp/(1+tmp) )
  
  tmp = as.numeric( exp( data$Xb %*% beta) )
  sigma.tmp = as.numeric( tmp/(1+tmp) )
  
  a = (1/sigma.tmp - 1) * mu.tmp
  b = (1/sigma.tmp - 1) * (1-mu.tmp) 
  a[a<0] = 0
  b[b<0] = 0
  
  return( -sum( (1-data$Del) * ( -.lbeta2(a,b) + data$A * (a-1) + data$B *(b-1) )  ) )
  
}


.Beta.neg.score.PAR2 <- function(par, data){
  
  da = ncol(data$Xa)
  alpha = par[1:da]
  beta = par[-(1:da)]
  
  tmp.a = as.numeric( exp( data$Xa %*% alpha) )
  mu.tmp = as.numeric( tmp.a/(1+tmp.a) )
  
  tmp.b = as.numeric( exp( data$Xb %*% beta) )
  sigma.tmp = as.numeric( tmp.b/(1+tmp.b) )
  
  a = (1/sigma.tmp - 1) * mu.tmp
  b = (1/sigma.tmp - 1) * (1-mu.tmp) 
  a[a<0] = 0
  b[b<0] = 0
  
  a.a = (1/tmp.b) * mu.tmp * (1/(1+tmp.a))
  a.b = -(1/tmp.b) * mu.tmp
  b.a = - a.a
  b.b = -(1/tmp.b) *  (1/(1+tmp.a))
  
  one = (1-data$Del)*(digamma(a+b)-digamma(a) + data$A)
  two = (1-data$Del)*(digamma(a+b)-digamma(b) + data$B)
  
  return( -c( colSums( (one*a.a + two*b.a) * data$Xa ), colSums( (one*a.b + two*b.b) * data$Xb ) ) )
}


.Beta.optim.PAR2 <- function(Del, A, B, Xa, Xb, alpha.ini, beta.ini){
  
  Beta.par.ini = c(alpha.ini, beta.ini) 
  Beta.data = list(Del=Del, A=A, B=B, Xa=Xa, Xb=Xb)
  
  #Beta.neg.loglik(Beta.par.ini, Beta.data)
  #Beta.neg.score(Beta.par.ini, Beta.data)
  return( optim(par=Beta.par.ini, fn=.Beta.neg.loglik.PAR2, gr=.Beta.neg.score.PAR2, data = Beta.data, method="BFGS")$par)
  
}


################################# Reparameterize: Estimation and Testing #################################

#make a simulated data
# EM algorithm
# Y: taxa count matrix (length: K+1 )
# X: design matrix with the first component as intercept (dim: n x d)
# c0: initial value for c (dim: K x d)
# alpha0: initial value for alpha (dim: K x d)
# beta0: initial value for beta (dim: K x d)
# tol: convergence criterion (scalar 0-1)
# max.iter: maximum number of iteration (scalar integer)
#ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)

.ZIGDM.EM.PAR2 <- function(Y, Xc, Xa, Xb, c0, alpha0, beta0, tol, max.iter){
  
  CONV = 0
  CONV.iter = max.iter
  n = nrow(Y)
  K = ncol(Y)-1
  dc = ncol(Xc)
  da = ncol(Xa)
  db = ncol(Xb)
  
  #c.all = list(c0); alpha.all = list(alpha0); beta.all = list(beta0); # save estimate from all the iterations (can be comment out in the final code)
  c.last = c0; alpha.last = alpha0; beta.last = beta0;
  c.now = c0; alpha.now = alpha0; beta.now = beta0;
  
  Del.R = matrix(NA, n, K)
  A.R = matrix(NA, n, K)
  B.R = matrix(NA, n, K)
  
  for(l in 1:max.iter){
    
    # E-step
    
    # print(paste("====== ", l, "th ======", sep=""))
    for(i in 1:n){
      
      tmp = exp(c.last %*% Xc[i,])
      tmp[is.na(tmp)] = 0    
      pv = as.numeric( tmp/(1+tmp) )
      pv[is.infinite(tmp) & tmp>0] =1  # positive inf
      
      tmp = exp(alpha.last %*% Xa[i,])
      mv = as.numeric( tmp/(1+tmp) )
      tmp = exp(beta.last %*% Xb[i,])
      sv = as.numeric( tmp/(1+tmp) )
      
      av = (1/sv - 1) * mv 
      bv = (1/sv - 1) * (1-mv)  
      
      par.post = .rP.Y.par(pv, av, bv, Y[i,])
      
      Del.R[i, ] = par.post$pv.post
      tmp = par.post$av.post + par.post$bv.post
      A.R[i, ] = digamma(par.post$av.post) - digamma(tmp)
      B.R[i, ] = digamma(par.post$bv.post) - digamma(tmp)
      
    }
    
    
    
    # M-step
    for(j in 1:K){
      
      if( !is.infinite(c.last[j,1]) ){
        c.now[j, ] = .Logit.optim(Del.R[,j], Xc, c.last[j, ])
      }
      
      tmp = .Beta.optim.PAR2(Del.R[,j], A.R[,j], B.R[,j], Xa, Xb, alpha.last[j, ], beta.last[j, ])
      alpha.now[j, ] = tmp[1:da]; beta.now[j, ] = tmp[-(1:da)]
      
    }
    
    
    diff = 0
    if( sum(!is.infinite(c.now))>0 ){
      
      diff = diff + sum(abs(c.now[!is.infinite(c.now)]-c.last[!is.infinite(c.now)])) 
      
    }
    
    diff = diff + sum(abs(alpha.now-alpha.last))
    diff = diff + sum(abs(beta.now-beta.last))
    if( diff < tol ){  
      
      CONV = 1; CONV.iter = l;
      break;
      
    }else{
      
      c.last = c.now; alpha.last = alpha.now; beta.last = beta.now;
    }
    
  }
  
  return(list(c.est = c.now, alpha.est = alpha.now, beta.est = beta.now, CONV=CONV, CONV.iter = CONV.iter))
  
}


.ZIGDM.stat.omni.perm <- function(X.perm, par.interest.index, S.list, I.list){
  
  n = nrow(X.perm)
  d = ncol(X.perm)
  K3 = length(S.list[[1]])
  Kd3 = K3 * d
  
  ES = rep(0,Kd3)
  I = matrix(0, Kd3, Kd3)
  
  for(i in 1:n){
    
    ES = ES + kronecker(S.list[[i]], X.perm[i,])
    X2 = X.perm[i,] %o% X.perm[i,]
    I = I + kronecker( I.list[[i]] , X2)     
    
  }
  
  S1 = ES[par.interest.index]
  I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
  I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
  I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
  score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
  
  return(score.stat)
  
}

#Y = Y.sim
#X = X.sim
#X.index = 2
#nperm = 1000

# default: only doing asymptotic test
# default: zero.index is the taxa with at least 1 zero observation

# ZIGDM.omni.test(Y2, X2, 2, zero.index=which(colSums(Y2[,1:(ncol(Y2)-1),drop=FALSE]==0)>0), nperm=100) 
# Y = Y2
# X = X2
# X.index = 2
# zero.index=which(colSums(Y2[,1:(ncol(Y2)-1),drop=FALSE]==0)>0)
# zero.index = NULL
#ZIGDM.omni.test(Y2, X2, 2, zero.index=NULL, nperm=100)
.ZIGDM.omni.test.PAR2 <- function(Y, X, X.index, zero.index=NULL, nperm=NULL, strata=NULL){
  
  X <- as.matrix(X[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  d = ncol(X)
  
  d.interest = length(X.index)
  X.null = X[,-X.index, drop=FALSE]
  d.null = ncol(X.null)
  
  # estimate parameters under the null model (omnibus test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, d.null) 
    
  }else{
    c0 = matrix(0, K, d.null) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, d.null) 
  beta0 = matrix(0, K, d.null) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, X.null, X.null, X.null, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  ### score function and information matrix 
  Kd3 = K * d * 3
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% X.null[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% X.null[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% X.null[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2 # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    ES = ES + kronecker(ESi, X[i,])
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      #ESS[index2, index2] = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #ESS[index3, index3] = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      #diag(ESS[index2, index2]) = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #diag(ESS[index3, index3]) = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    
    X2 = X[i,] %o% X[i,]
    # I_1 = I_1 + kronecker(EI, X2)
    # I_2 = I_2 + kronecker(ESS, X2)
    # I_3 = I_3 + kronecker(ES.2, X2)
    I.list[[i]] = (-EI - ESS + ES.2)
    I = I + kronecker( I.list[[i]] , X2)  
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  tmp = c(X.index, X.index+d, X.index+2*d)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+3*d
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  S1 = ES[par.interest.index]
  #I = - I_1 - I_2 + I_3
  I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
  I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
  I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
  score.stat.omni = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
  score.pval.omni = 1-pchisq(score.stat.omni, length(par.interest.index))
  
  
  if(!is.null(nperm)){
    
    set.seed(11)
    count.omni = 0
    
    X.perm = X
    for(l in 1:nperm){
      if(is.null(strata)){
        X.perm[,X.index] = X[sample(1:n),X.index]
      }else{
        X.perm[,X.index] = X[.sampling.strata(1:n, strata),X.index]
        
      }
      tmp.omni = .ZIGDM.stat.omni.perm(X.perm, par.interest.index, S.list, I.list)
      
      if(tmp.omni >= score.stat.omni){
        count.omni = count.omni + 1
      }      
      
    }
    
    if(count.omni<=2){
      warnings("permutation p-value for omnibus test may not be accurate.")
    }
    
    score.pval.perm.omni = (count.omni+1)/(nperm+1)
    
  }else{
    
    score.pval.perm.omni = NA
    
  }
  
  
  return(list(score.stat.omni=score.stat.omni, score.pval.omni = score.pval.omni, score.pval.perm.omni = score.pval.perm.omni, df=df))
  
}

## FOR permutation test of each component
# if target == "freq", then only test for frequency
# if target == "abun", then only test for abundance
# if target == "disp", then only test for dispersion
# ZIGDM.stat.target.perm(Xc, Xa, Xb.perm, par.interest.index, S.list, I.list)
.ZIGDM.stat.target.perm <- function(Xc, Xa, Xb, par.interest.index, S.list, I.list){
  
  n = nrow(Xc)
  dc = ncol(Xc)
  da = ncol(Xa)
  db = ncol(Xb)
  dd = da + db + dc
  K3 = length(S.list[[1]])
  K = K3/3
  Kd3 = K * dd
  
  ES = rep(0,Kd3)
  I = matrix(0, Kd3, Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  
  for(i in 1:n){
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])   
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,])   
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,])
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,])
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])
      
      
    }
    
    
  } 
  
  
  
  S1 = ES[par.interest.index]
  I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
  I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
  I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
  score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
  
  return(score.stat)
  
}


.ZIGDM.stat.freq.perm <- function(Y, Xc, Xa, Xb, Xc.index, zero.index=NULL){
  
  Xc <- as.matrix(Xc[rowSums(Y) != 0, ,drop=FALSE])
  Xa <- as.matrix(Xa[rowSums(Y) != 0, ,drop=FALSE])
  Xb <- as.matrix(Xb[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  da = ncol(Xa)
  db = ncol(Xb)
  dc = ncol(Xc)
  
  dc.interest = length(Xc.index)
  Xc.null = Xc[,-Xc.index, drop=FALSE]
  dc.null = ncol(Xc.null)
  
  # estimate parameters under the null model (absent frequency test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, dc.null) 
    
  }else{
    c0 = matrix(0, K, dc.null) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, da) 
  beta0 = matrix(0, K, db) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, Xc.null, Xa, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  ### score function and information matrix 
  dd = da + db + dc
  Kd3 = K * dd
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% Xc.null[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% Xa[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% Xb[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    #ES = ES + kronecker(ESi, X[i,])
    
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    
    I.list[[i]] = (-EI - ESS + ES.2)
    #I = I + kronecker( I.list[[i]] , X2)  
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,]) 
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,]) 
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,]) 
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])      
      
    }    
    
    
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  #tmp = c(X.index, X.index+d, X.index+2*d)
  tmp = c(Xc.index)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+dd
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  if(df>0){
    
    S1 = ES[par.interest.index]
    I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
    I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
    I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
    score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
    
  }else{
    score.stat = NA
  }
  return(score.stat)
  
}

.ZIGDM.freq.test.PAR2 <- function(Y, Xc, Xa, Xb, Xc.index, zero.index=NULL, nperm=NULL, strata=NULL){
  
  Xc <- as.matrix(Xc[rowSums(Y) != 0, ,drop=FALSE])
  Xa <- as.matrix(Xa[rowSums(Y) != 0, ,drop=FALSE])
  Xb <- as.matrix(Xb[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  S1 = NA
  S1cov = NA
  nuis.est = NA
  
  nuis = 0   # if nuisance parameters are more than intercept term, then nuis = 1
  if( ncol(Xa)>1 | ncol(Xb)>1 ){
    nuis = 1
  }
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  da = ncol(Xa)
  db = ncol(Xb)
  dc = ncol(Xc)
  
  dc.interest = length(Xc.index)
  Xc.null = Xc[,-Xc.index, drop=FALSE]
  dc.null = ncol(Xc.null)
  
  # estimate parameters under the null model (absent frequency test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, dc.null) 
    
  }else{
    c0 = matrix(0, K, dc.null) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, da) 
  beta0 = matrix(0, K, db) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, Xc.null, Xa, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  nuis.est = list(zero.par.est = c.est, mean.par.est = alpha.est, disp.par.est = beta.est)
  
  ### score function and information matrix 
  dd = da + db + dc
  Kd3 = K * dd
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% Xc.null[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% Xa[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% Xb[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    #ES = ES + kronecker(ESi, X[i,])
    
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    
    I.list[[i]] = (-EI - ESS + ES.2)
    #I = I + kronecker( I.list[[i]] , X2)  
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,]) 
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,]) 
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,]) 
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])      
      
    }    
    
    
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  #tmp = c(X.index, X.index+d, X.index+2*d)
  tmp = c(Xc.index)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+dd
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  if(df>0){
    
    S1 = ES[par.interest.index]
    I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
    I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
    I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
    S1cov = I.11 - I.12%*%ginv(I.22)%*%t(I.12)
    score.stat = S1 %*% ginv( S1cov ) %*% S1
    score.pval = 1-pchisq(score.stat, length(par.interest.index))
    
    ########################################
    #                                      #
    #        Resampling Score Test         #
    # change 04/23/2018 for adaptive resampling
    #                                      #
    ######################################## 
    
    if(!is.null(nperm) & nuis == 0){
      
      set.seed(11)
      count = 0
      Xc.perm = Xc
      
      flag = 1
      start.nperm = 1
      end.nperm = min(100,nperm);
      
      while(flag){
        
        for(l in start.nperm:end.nperm){
          
          
          if(is.null(strata)){
            Xc.perm[,Xc.index] = Xc[sample(1:n),Xc.index]
          }else{
            Xc.perm[,Xc.index] = Xc[.sampling.strata(1:n, strata),Xc.index]
            
          }
          tmp = .ZIGDM.stat.target.perm(Xc.perm, Xa, Xb, par.interest.index, S.list, I.list)
          
          if(tmp >= score.stat){
            count = count + 1
          }      
          
        }
        
        if(count>=5 | end.nperm == nperm){
          flag = 0
          
        }else{
          
          start.nperm = end.nperm + 1
          end.nperm = min(end.nperm*5, nperm)
          
        }
        
      }
      
      
      #if(count<=2){
        #warning(paste( "Inaccurate pvalue with", nperm, "resamplings"))
      #}
      
      score.pval.perm = (count+1)/(end.nperm+1)
      
      
    }
    else if(!is.null(nperm) & nuis == 1){
      
      set.seed(11)
      count = 0
      Xc.perm = Xc
      
      flag = 1
      start.nperm = 1
      end.nperm = min(100,nperm);
      
      while(flag){
        
        for(l in start.nperm:end.nperm){
          
          
          if(is.null(strata)){
            perm.idx = sample(1:n)
          }else{
            perm.idx = .sampling.strata(1:n, strata)
          }
          Xc.perm = Xc[perm.idx,]
          Xa.perm = Xa[perm.idx,]
          Xb.perm = Xb[perm.idx,]
          
          tmp = .ZIGDM.stat.freq.perm(Y, Xc.perm, Xa.perm, Xb.perm, Xc.index, zero.index=zero.index)
          
          if(tmp >= score.stat){
            count = count + 1
          }      
          
        }
        
        if(count>=5 | end.nperm == nperm){
          flag = 0
          
        }else{
          
          start.nperm = end.nperm + 1
          end.nperm = min(end.nperm*5, nperm)
          
        }
        
      }
      
      
      # if(count<=2){
      #   warning(paste( "Inaccurate pvalue with", nperm, "resamplings"))
      # }
      
      score.pval.perm = (count+1)/(end.nperm+1)
      
    }
    else{
      
      score.pval.perm = NA
      
    }
    
    
  }else{
    score.stat = NA
    score.pval = NA
    score.pval.perm = NA
  }
  
  return(list(score.pval.freq = score.pval, score.pval.perm.freq = score.pval.perm, df=df, nuis.est = nuis.est, score.S1 = S1, score.S1cov = S1cov, score.stat=as.numeric(score.stat) ))
  
}


.ZIGDM.stat.abun.perm <- function(Y, Xc, Xa, Xb, Xa.index, zero.index=NULL){
  
  
  Xc <- as.matrix(Xc[rowSums(Y) != 0, ,drop=FALSE])
  Xa <- as.matrix(Xa[rowSums(Y) != 0, ,drop=FALSE])
  Xb <- as.matrix(Xb[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  da = ncol(Xa)
  db = ncol(Xb)
  dc = ncol(Xc)
  
  da.interest = length(Xa.index)
  Xa.null = Xa[,-Xa.index, drop=FALSE]
  da.null = ncol(Xa.null)
  
  # estimate parameters under the null model (abundance test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, dc) 
    
  }else{
    c0 = matrix(0, K, dc) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, da.null) 
  beta0 = matrix(0, K, db) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, Xc, Xa.null, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  ### score function and information matrix 
  dd = da + db + dc
  Kd3 = K * dd
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% Xc[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% Xa.null[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% Xb[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    #ES = ES + kronecker(ESi, X[i,])
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      #ESS[index2, index2] = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #ESS[index3, index3] = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      #diag(ESS[index2, index2]) = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #diag(ESS[index3, index3]) = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    I.list[[i]] = (-EI - ESS + ES.2)
    #I = I + kronecker( I.list[[i]] , X2)  
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,]) 
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,]) 
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,]) 
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])      
      
    }    
    
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  #tmp = c(X.index, X.index+d, X.index+2*d)
  tmp = c(Xa.index+dc)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+dd
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  if(df>0){
    
    S1 = ES[par.interest.index]
    #I = - I_1 - I_2 + I_3
    I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
    I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
    I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
    score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
    
  }else{
    score.stat = NA
  }
  return(score.stat)
  
}

.ZIGDM.abun.test.PAR2 <- function(Y, Xc, Xa, Xb, Xa.index, zero.index=NULL, nperm=NULL, strata=NULL){
  
  Xc <- as.matrix(Xc[rowSums(Y) != 0, ,drop=FALSE])
  Xa <- as.matrix(Xa[rowSums(Y) != 0, ,drop=FALSE])
  Xb <- as.matrix(Xb[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  S1 = NA
  S1cov = NA
  nuis.est = NA
  
  nuis = 0   # if nuisance parameters are more than intercept term, then nuis = 1
  if( ncol(Xc)>1 | ncol(Xb)>1 ){
    nuis = 1
  }
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  da = ncol(Xa)
  db = ncol(Xb)
  dc = ncol(Xc)
  
  da.interest = length(Xa.index)
  Xa.null = Xa[,-Xa.index, drop=FALSE]
  da.null = ncol(Xa.null)
  
  # estimate parameters under the null model (abundance test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, dc) 
    
  }else{
    c0 = matrix(0, K, dc) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, da.null) 
  beta0 = matrix(0, K, db) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, Xc, Xa.null, Xb, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  nuis.est = list(zero.par.est = c.est, mean.par.est = alpha.est, disp.par.est = beta.est)
  
  ### score function and information matrix 
  dd = da + db + dc
  Kd3 = K * dd
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% Xc[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% Xa.null[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% Xb[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    #ES = ES + kronecker(ESi, X[i,])
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      #ESS[index2, index2] = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #ESS[index3, index3] = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      #diag(ESS[index2, index2]) = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #diag(ESS[index3, index3]) = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    I.list[[i]] = (-EI - ESS + ES.2)
    #I = I + kronecker( I.list[[i]] , X2)  
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,]) 
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,]) 
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,]) 
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])      
      
    }    
    
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  #tmp = c(X.index, X.index+d, X.index+2*d)
  tmp = c(Xa.index+dc)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+dd
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  if(df>0){
    
    S1 = ES[par.interest.index]
    #I = - I_1 - I_2 + I_3
    I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
    I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
    I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
    S1cov = I.11 - I.12%*%ginv(I.22)%*%t(I.12)
    score.stat = S1 %*% ginv( S1cov ) %*% S1
    score.pval = 1-pchisq(score.stat, length(par.interest.index))
    
    if(!is.null(nperm) & nuis == 0){
      
      set.seed(11)
      count = 0
      Xa.perm = Xa
      
      flag = 1
      start.nperm = 1
      end.nperm = min(100,nperm);
      
      while(flag){
        
        for(l in start.nperm:end.nperm){
          
          if(is.null(strata)){
            Xa.perm[,Xa.index] = Xa[sample(1:n),Xa.index]
          }else{
            Xa.perm[,Xa.index] = Xa[.sampling.strata(1:n, strata),Xa.index]
            
          }
          
          tmp = .ZIGDM.stat.target.perm(Xc, Xa.perm, Xb, par.interest.index, S.list, I.list)
          
          if(tmp >= score.stat){
            count = count + 1
          }      
          
        }
        
        if(count>=5 | end.nperm == nperm){
          flag = 0
          
        }else{
          
          start.nperm = end.nperm + 1
          end.nperm = min(end.nperm*5, nperm)
          
        }
        
      }
      
      
      # if(count<=2){
      #   warning(paste( "Inaccurate pvalue with", nperm, "resamplings"))
      # }
      
      score.pval.perm = (count+1)/(end.nperm+1)
      
      
    }
    else if(!is.null(nperm) & nuis == 1){
      
      set.seed(11)
      count = 0
      Xa.perm = Xa
      
      flag = 1
      start.nperm = 1
      end.nperm = min(100,nperm);
      
      while(flag){
        
        for(l in start.nperm:end.nperm){
          
          
          if(is.null(strata)){
            perm.idx = sample(1:n)
          }else{
            perm.idx = .sampling.strata(1:n, strata)
          }
          
          Xa.perm = Xa[perm.idx,]
          Xb.perm = Xb[perm.idx,]
          Xc.perm = Xc[perm.idx,]
          
          tmp = .ZIGDM.stat.abun.perm(Y, Xc.perm, Xa.perm, Xb.perm, Xa.index, zero.index=zero.index)
          
          if(tmp >= score.stat){
            count = count + 1
          }       
          
        }
        
        if(count>=5 | end.nperm == nperm){
          flag = 0
          
        }else{
          
          start.nperm = end.nperm + 1
          end.nperm = min(end.nperm*5, nperm)
          
        }
        
      }
      
      
      # if(count<=2){
      #   warning(paste( "Inaccurate pvalue with", nperm, "resamplings"))
      # }
      
      score.pval.perm = (count+1)/(end.nperm+1)
      
      
    }
    else{
      
      score.pval.perm = NA
      
    }
    
  }else{
    score.stat = NA
    score.pval = NA
    score.pval.perm = NA
  }
  
  return(list(score.pval.abun = score.pval, score.pval.perm.abun = score.pval.perm, df=df, nuis.est = nuis.est, score.S1 = S1, score.S1cov = S1cov, score.stat=as.numeric(score.stat) ))
  
}


.ZIGDM.stat.disp.perm <- function(Y, Xc, Xa, Xb, Xb.index, zero.index=NULL){
  
  Xc <- as.matrix(Xc[rowSums(Y) != 0, ,drop=FALSE])
  Xa <- as.matrix(Xa[rowSums(Y) != 0, ,drop=FALSE])
  Xb <- as.matrix(Xb[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  da = ncol(Xa)
  db = ncol(Xb)
  dc = ncol(Xc)
  
  db.interest = length(Xb.index)
  Xb.null = Xb[,-Xb.index, drop=FALSE]
  db.null = ncol(Xb.null)
  
  # estimate parameters under the null model (abundance test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, dc) 
    
  }else{
    c0 = matrix(0, K, dc) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, da) 
  beta0 = matrix(0, K, db.null) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, Xc, Xa, Xb.null, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  ### score function and information matrix 
  dd = da + db + dc
  Kd3 = K * dd
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% Xc[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% Xa[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% Xb.null[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    #ES = ES + kronecker(ESi, X[i,])
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      #ESS[index2, index2] = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #ESS[index3, index3] = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      #diag(ESS[index2, index2]) = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #diag(ESS[index3, index3]) = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    
    
    I.list[[i]] = (-EI - ESS + ES.2)
    #I = I + kronecker( I.list[[i]] , X2) 
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,]) 
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,]) 
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,]) 
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])      
      
    }    
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  #tmp = c(X.index, X.index+d, X.index+2*d)
  tmp = c(Xb.index+dc+da)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+dd
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  if(df>0){
    
    S1 = ES[par.interest.index]
    #I = - I_1 - I_2 + I_3
    I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
    I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
    I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
    score.stat = S1 %*% ginv( I.11 - I.12%*%ginv(I.22)%*%t(I.12) ) %*% S1
    
    
  }else{
    score.stat = NA
  }
  return(score.stat)
  
}

.ZIGDM.disp.test.PAR2 <- function(Y, Xc, Xa, Xb, Xb.index, zero.index=NULL, nperm=NULL, strata=NULL){
  
  Xc <- as.matrix(Xc[rowSums(Y) != 0, ,drop=FALSE])
  Xa <- as.matrix(Xa[rowSums(Y) != 0, ,drop=FALSE])
  Xb <- as.matrix(Xb[rowSums(Y) != 0, ,drop=FALSE])
  Y <- as.matrix(Y[rowSums(Y) != 0, colSums(Y) != 0, drop=FALSE])
  
  S1 = NA
  S1cov = NA
  nuis.est = NA
  
  nuis = 0   # if nuisance parameters are more than intercept term, then nuis = 1
  if( ncol(Xc)>1 | ncol(Xa)>1 ){
    nuis = 1
  }
  
  n = nrow(Y)
  K = ncol(Y)-1
  
  tmp0 = which(colSums(Y[,1:K, drop=FALSE]==0)>0)
  zero.index = zero.index[zero.index %in%  tmp0]
  
  da = ncol(Xa)
  db = ncol(Xb)
  dc = ncol(Xc)
  
  db.interest = length(Xb.index)
  Xb.null = Xb[,-Xb.index, drop=FALSE]
  db.null = ncol(Xb.null)
  
  # estimate parameters under the null model (abundance test)
  if( is.null(zero.index) ){
    
    c0 = matrix(-Inf, K, dc) 
    
  }else{
    c0 = matrix(0, K, dc) 
    c0[-zero.index,] = -Inf
  }
  
  alpha0 = matrix(0, K, da) 
  beta0 = matrix(0, K, db.null) 
  #tol = 0.0001
  #max.iter = 1000
  est = .ZIGDM.EM.PAR2(Y, Xc, Xa, Xb.null, c0, alpha0, beta0, tol=0.0001, max.iter=1000)
  #tol = 0.00001
  #est2 = ZIGDM.EM(Y, X.null, c0, alpha0, beta0, tol, max.iter)
  c.est = est$c.est
  alpha.est = est$alpha.est
  beta.est = est$beta.est
  
  nuis.est = list(zero.par.est = c.est, mean.par.est = alpha.est, disp.par.est = beta.est)
  
  ### score function and information matrix 
  dd = da + db + dc
  Kd3 = K * dd
  K3 = K * 3
  #I_1 = matrix(0, Kd3, Kd3)
  #I_2 = matrix(0, Kd3, Kd3)
  #I_3 = matrix(0, Kd3, Kd3)
  I = matrix(0, Kd3, Kd3)
  ES = rep(0,Kd3)
  
  index1 = ( (1:K)-1 )*3 + 1
  index2 = ( (1:K)-1 )*3 + 2
  index3 = (1:K)*3
  EI = matrix(0, K3, K3)
  ESi = rep(0,K3)
  ESS =  matrix(0, K3, K3)
  
  I.list = list()
  S.list = list()
  for(i in 1:n){
    
    tmp = exp(c.est %*% Xc[i,])
    tmp[is.na(tmp)] = 0              # negative Inf
    pv = as.numeric( tmp/(1+tmp) )
    pv[is.infinite(tmp) & tmp>0] =1  # positive inf
    
    expa = exp(alpha.est %*% Xa[i,])
    mv = as.numeric( expa/(1+expa) )
    expb = exp(beta.est %*% Xb.null[i,])
    sv = as.numeric( expb/(1+expb) )
    
    av = (1/sv - 1) * mv 
    bv = (1/sv - 1) * (1-mv)  
    abv = av + bv
    
    diga.av = digamma(av)
    diga.bv = digamma(bv)
    diga.abv = digamma(abv)
    A = diga.av - diga.abv 
    B = diga.bv - diga.abv
    
    triga.av = trigamma(av)
    triga.bv = trigamma(bv)
    triga.abv = trigamma(abv)
    A2 = triga.av - triga.abv 
    B2 = triga.bv - triga.abv
    
    # stat basaed on posterious prob
    par.post = .rP.Y.par(pv, av, bv, Y[i,])
    # post.expectation of delta
    Del.post = par.post$pv.post
    av.post = par.post$av.post
    bv.post = par.post$bv.post
    abv.post = av.post + bv.post
    
    # post.expectation of logZ conditional on delta==0
    A.post = digamma(av.post) - digamma(abv.post) 
    # post.expectation of log(1-Z) conditional on delta==0
    B.post = digamma(bv.post) - digamma(abv.post) 
    # post.expectation of logZ*logZ conditional on delta==0  
    A2.post = trigamma(av.post) - trigamma(abv.post) + A.post^2
    # post.expectation of log(1-Z)*log(1-Z) conditional on delta==0  
    B2.post = trigamma(bv.post) - trigamma(abv.post) + B.post^2
    # post.expectation of logZ*log(1-Z) conditional on delta==0  
    AB.post = -trigamma(abv.post) + A.post * B.post
    
    one = (1-Del.post)*( - A + A.post )
    two = (1-Del.post)*( - B + B.post )
    
    # derivative of a, b with respect to alpha, beta
    b.tmp = 1/expb
    a0.tmp = 1/(1+expa)
    a1.tmp = expa/(1+expa)
    a2.tmp = expa/(1+expa)^2
    a3.tmp = expa/(1+expa)^3
    a.a = b.tmp * a2.tmp
    a.b = - b.tmp * a1.tmp
    b.a = - b.tmp * a2.tmp
    b.b = - b.tmp * a0.tmp
    a.a2 = b.tmp * a3.tmp * (1-expa)  # 10/12/2017 correct
    a.ab = - b.tmp * a2.tmp
    a.b2 = b.tmp * a1.tmp
    b.a2 = - b.tmp * a3.tmp * (1-expa) # 10/12/2017 correct
    b.ab = b.tmp * a2.tmp
    b.b2 = b.tmp * a0.tmp
    
    ############ EI
    if(K==1){
      
      # second derivative associated with c
      EI[index1, index1] = -pv*(1-pv)
      # second derivative associated with alpha
      EI[index2, index2] = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      EI[index3, index3] = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      EI[index2, index3] = tmp
      EI[index3, index2] = tmp   
      
    }else{
      # second derivative associated with c
      diag(EI[index1, index1]) = -pv*(1-pv)
      # second derivative associated with alpha
      diag(EI[index2, index2]) = (1-Del.post)*( - triga.av*a.a^2 - triga.bv*b.a^2) + one*a.a2 + two*b.a2  # 10/12/2017 remove 
      # second derivative associated with beta
      diag(EI[index3, index3]) = (1-Del.post)*( triga.abv*(a.b+b.b)^2 - triga.av*a.b^2 - triga.bv*b.b^2) + one*a.b2 + two*b.b2
      # second derivative associated with alpha, beta
      tmp = (1-Del.post)*( - triga.av*a.b*a.a - triga.bv*b.b*b.a) + one*a.ab + two*b.ab  # 10/12/2017 remove 
      diag(EI[index2, index3]) = tmp
      diag(EI[index3, index2]) = tmp      
      
    }
    
    
    
    ############ ES.2
    ESi[index1] = Del.post - pv
    ESi[index2] = one*a.a + two*b.a    
    ESi[index3] = one*a.b + two*b.b
    ES.2 = ESi %o% ESi
    S.list[[i]] = ESi
    #ES = ES + kronecker(ESi, X[i,])
    
    ############ ESS
    ESS = ES.2
    if(K==1){
      
      ESS[index1, index1] = Del.post - 2 * Del.post * pv + pv * pv
      #ESS[index2, index2] = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #ESS[index3, index3] = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      ESS[index2, index2] = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      ESS[index3, index3] = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      ESS[index1, index2] = tmp
      ESS[index2, index1] = tmp
      tmp = -pv * ESi[index3]
      ESS[index1, index3] = tmp
      ESS[index3, index1] = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      ESS[index2, index3] = tmp
      ESS[index3, index2] = tmp  
      
    }else{
      diag(ESS[index1, index1]) = Del.post - 2 * Del.post * pv + pv * pv
      #diag(ESS[index2, index2]) = ((1-Del.post)*av)^2 * (A^2 - 2*A*A.post + A2.post)
      #diag(ESS[index3, index3]) = ((1-Del.post)*bv)^2 * (B^2 - 2*B*B.post + B2.post)
      tmp1 = -A * a.a - B * b.a
      diag(ESS[index2, index2]) = (1-Del.post)*(tmp1*tmp1 + A2.post*a.a^2 + B2.post*b.a^2 + 2*A.post*tmp1*a.a + 2*B.post*tmp1*b.a + 2*AB.post*a.a*b.a)
      
      tmp2 = -A * a.b - B * b.b
      diag(ESS[index3, index3]) = (1-Del.post)*(tmp2*tmp2 + A2.post*a.b^2 + B2.post*b.b^2 + 2*A.post*tmp2*a.b + 2*B.post*tmp2*b.b + 2*AB.post*a.b*b.b)
      
      tmp = -pv * ESi[index2]
      diag(ESS[index1, index2]) = tmp
      diag(ESS[index2, index1]) = tmp
      tmp = -pv * ESi[index3]
      diag(ESS[index1, index3]) = tmp
      diag(ESS[index3, index1]) = tmp   
      tmp = (1-Del.post)*( tmp1*tmp2 + (tmp1*a.b + tmp2*a.a)*A.post + (tmp1*b.b + tmp2*b.a)*B.post + (a.b*b.a + a.a*b.b)*AB.post + (a.a*a.b)*A2.post + (b.b*b.a)*B2.post )
      diag(ESS[index2, index3]) = tmp
      diag(ESS[index3, index2]) = tmp  
    }
    
    
    
    
    
    I.list[[i]] = (-EI - ESS + ES.2)
    #I = I + kronecker( I.list[[i]] , X2) 
    
    b.end = 0
    for(j in 1:K){
      
      c.start = b.end + 1; c.end = c.start + dc - 1; 
      a.start = c.end + 1; a.end = a.start + da - 1; 
      b.start = a.end + 1; b.end = b.start + db - 1;
      
      ES[c.start:c.end] = ES[c.start:c.end] + S.list[[i]][index1[j]] * Xc[i,]
      ES[a.start:a.end] = ES[a.start:a.end] + S.list[[i]][index2[j]] * Xa[i,]
      ES[b.start:b.end] = ES[b.start:b.end] + S.list[[i]][index3[j]] * Xb[i,]
      
      I[c.start:c.end, c.start:c.end] = I[c.start:c.end, c.start:c.end] + I.list[[i]][index1[j], index1[j]] * (Xc[i,] %o% Xc[i,])
      I[c.start:c.end, a.start:a.end] = I[c.start:c.end, a.start:a.end] + I.list[[i]][index1[j], index2[j]] * (Xc[i,] %o% Xa[i,])
      I[a.start:a.end, c.start:c.end] = t(I[c.start:c.end, a.start:a.end])
      I[c.start:c.end, b.start:b.end] = I[c.start:c.end, b.start:b.end] + I.list[[i]][index1[j], index3[j]] * (Xc[i,] %o% Xb[i,])
      I[b.start:b.end, c.start:c.end] = t(I[c.start:c.end, b.start:b.end])
      I[a.start:a.end, a.start:a.end] = I[a.start:a.end, a.start:a.end] + I.list[[i]][index2[j], index2[j]] * (Xa[i,] %o% Xa[i,]) 
      I[b.start:b.end, b.start:b.end] = I[b.start:b.end, b.start:b.end] + I.list[[i]][index3[j], index3[j]] * (Xb[i,] %o% Xb[i,]) 
      I[a.start:a.end, b.start:b.end] = I[a.start:a.end, b.start:b.end] + I.list[[i]][index2[j], index3[j]] * (Xa[i,] %o% Xb[i,]) 
      I[b.start:b.end, a.start:a.end] = t(I[a.start:a.end, b.start:b.end])      
      
    }    
    
  }
  
  
  ## omnibus score test
  #par.interest.index =  kronecker(X.index, 1:(3*K) )
  #tmp = c(X.index, X.index+d, X.index+2*d)
  tmp = c(Xb.index+dc+da)
  par.interest.index = tmp
  if(K>1){
    
    for(k in 1:(K-1)){
      tmp = tmp+dd
      par.interest.index = c(par.interest.index, tmp)
      
    }
    
  }
  
  par.remove.index = which(.diag2(I)[par.interest.index]==0)
  if(length(par.remove.index)>0){
    par.interest.index = par.interest.index[-par.remove.index]
  }
  
  df = length(par.interest.index)
  if(df>0){
    
    S1 = ES[par.interest.index]
    #I = - I_1 - I_2 + I_3
    I.11 = I[par.interest.index, par.interest.index, drop=FALSE]
    I.12 = I[par.interest.index, -par.interest.index, drop=FALSE]
    I.22 = I[-par.interest.index, -par.interest.index, drop=FALSE]
    S1cov = I.11 - I.12%*%ginv(I.22)%*%t(I.12)
    score.stat = S1 %*% ginv( S1cov ) %*% S1
    score.pval = 1-pchisq(score.stat, length(par.interest.index))
    
    
    if(!is.null(nperm) & nuis == 0){
      
      set.seed(11)
      count = 0
      Xb.perm = Xb
      
      flag = 1
      start.nperm = 1
      end.nperm = min(100,nperm);
      
      while(flag){
        
        for(l in start.nperm:end.nperm){
          
          if(is.null(strata)){
            Xb.perm[,Xb.index] = Xb[sample(1:n),Xb.index]
          }else{
            Xb.perm[,Xb.index] = Xb[.sampling.strata(1:n, strata),Xb.index]
            
          }
          tmp = .ZIGDM.stat.target.perm(Xc, Xa, Xb.perm, par.interest.index, S.list, I.list)
          #score.stat.perm[l] = tmp
          if(tmp >= score.stat){
            count = count + 1
          }     
          
        }
        
        if(count>=5 | end.nperm == nperm){
          flag = 0
          
        }else{
          
          start.nperm = end.nperm + 1
          end.nperm = min(end.nperm*5, nperm)
          
        }
        
      }
      
      
      # if(count<=2){
      #   warning(paste( "Inaccurate pvalue with", nperm, "resamplings"))
      # }
      
      score.pval.perm = (count+1)/(end.nperm+1)
      
      
    }
    else if(!is.null(nperm) & nuis == 1){
      
      set.seed(11)
      count = 0
      Xb.perm = Xb
      
      flag = 1
      start.nperm = 1
      end.nperm = min(100,nperm);
      
      while(flag){
        
        for(l in start.nperm:end.nperm){
          
          if(is.null(strata)){
            perm.idx = sample(1:n)
            
          }else{
            perm.idx = .sampling.strata(1:n, strata)
            
            
          }
          Xc.perm = Xc[perm.idx,]
          Xa.perm = Xa[perm.idx,]
          Xb.perm = Xb[perm.idx,]
          
          tmp = .ZIGDM.stat.disp.perm(Y, Xc.perm, Xa.perm, Xb.perm, Xb.index, zero.index = zero.index)
          #score.stat.perm[l] = tmp
          if(tmp >= score.stat){
            count = count + 1
          }      
          
        }
        
        if(count>=5 | end.nperm == nperm){
          flag = 0
          
        }else{
          
          start.nperm = end.nperm + 1
          end.nperm = min(end.nperm*5, nperm)
          
        }
        
      }
      
      
      # if(count<=2){
      #   warning(paste( "Inaccurate pvalue with", nperm, "resamplings"))
      # }
      
      score.pval.perm = (count+1)/(end.nperm+1)
      
      
    }
    else{
      
      score.pval.perm = NA
      
    }
    
  }else{
    score.stat = NA
    score.pval = NA
    score.pval.perm = NA
  }
  return(list(score.pval.disp = score.pval, score.pval.perm.disp = score.pval.perm, df=df, nuis.est = nuis.est, score.S1 = S1, score.S1cov = S1cov, score.stat=as.numeric(score.stat) ))
  
}
