# Run from the MetaQCAT-repro repository root.
if (!exists("project.dir", inherits = TRUE))
  source(file.path("R", "project_paths.R"))

library(miLineage)
library(data.table)
library(dplyr)
library(parallel)

F.test <- function(x){
  x.stat = -2 * sum(log(x))
  return( 1 - pchisq(x.stat, df = 2 * length(x)) )
}

lineage.QCAT.Cluster <- function(tax.info, W.data.list, otucols, X, X.index, ID, perm.type = NULL, n.perm=NULL, test = "chisq"){
  Rank.low = unname(tax.info['Rank.low'])
  Rank.high = unname(tax.info['Rank.high'])
  level.uni = unname(tax.info['level.uni'])
  n.OTU = length(X)
  # partition and merge OTU table according to taxonomy information
  tt = lapply(1:n.OTU, function(j) W.data.list[[j]][, lapply(.SD , sum, na.rm=TRUE), .SDcols=as.vector(unlist(otucols[j])), by=list( get(Rank.low), get(Rank.high) )])
  tt = lapply(1:n.OTU,function(j) setnames(tt[[j]], 1:2, c(Rank.low, Rank.high)))
  W.count = lapply(1:n.OTU,function(j) tt[[j]][, otucols[[j]], with=FALSE])
  W.tax = as.vector(unlist(tt[[1]][, Rank.low, with=FALSE]))
  tax.names <- as.vector(unlist(tt[[1]][, Rank.high, with=FALSE]))[which(W.tax == level.uni)]
  Y = lapply(1:n.OTU, function(i) t(W.count[[i]][which(W.tax == level.uni), , drop=FALSE]))
  if(is.null(Y)){
    return(NA)
  }
  if(ncol(Y[[1]])==1){
    return(NA)
  }
  zero.ave.lst <- lapply(Y, function(x){colMeans(x==0)})
  zero.cnt = Reduce('+', zero.ave.lst)
  zero.cnt.mat <- Reduce(rbind, zero.ave.lst)
  non.base.id <- apply(zero.cnt.mat, 2, function(x){any(x == 1)})
  if(sum(non.base.id) == ncol(zero.cnt.mat)){
    remove.stduy.cnt <- apply(zero.cnt.mat, 2, function(x){sum(x == 1)})
    reference.column.candidates = which(remove.stduy.cnt == min(remove.stduy.cnt))
    reference.column = reference.column.candidates[rev(which(zero.cnt[reference.column.candidates] == min(zero.cnt[reference.column.candidates])))[1]]
    remove.stduy.id = which(zero.cnt.mat[,reference.column] == 1)
    Y = Y[-remove.stduy.id]
    X = X[-remove.stduy.id]
    if(!is.null(ID)){
      ID = ID[-remove.stduy.id] 
    }
    # return(NA) # no proper base for positive part
  }else if(sum(non.base.id)>0){
    zero.cnt[non.base.id] <- Inf
    reference.column = rev(which(zero.cnt == min(zero.cnt)))[1]
  }else{
    reference.column = rev(which(zero.cnt == min(zero.cnt)))[1]
  }
  
  re.order <- c(c(1:ncol(Y[[1]]))[-reference.column], reference.column)
  Y =  lapply(Y, function(x){x[, re.order]})# perform one test using all OTUs
  # Y =  lapply(Y, function(x){x[, order(-zero.cnt)]})
  tax.names.reorder = tax.names[re.order]
  reference.taxon <- tax.names[reference.column]
  
  # subtree = c(subtree, level.uni[j])
  #print(level.uni[j])
  if(is.null(n.perm)){ # asymptotic test only
    # (Y.list, X.list, X.par.index, seed=11, resample=FALSE, n.perm=NULL, meta.method = "FE-MetaQCAT", Weight=NULL )
    #  run test for each lineage
    pval_lst_asy = c()
    for(i in 1:length(Y)){
      tmp = try(QCAT_GEE.Cluster(ID = ID[[i]], OTU = Y[[i]], OTU.base=NULL, X = X[[i]], X.index = X.index, Z = X[[i]], Z.index = X.index,
                                 Tax=NULL, min.depth=0, perm.type=perm.type, n.perm=n.perm, fdr.alpha=0.05, test= test))
      if(!("try-error" %in% class(tmp))){
        pval_lst_asy = append(pval_lst_asy,tmp$pval[1])
      }
    }
    if(is.null(pval_lst_asy)){
      return(NA)
    }else{
      pval = F.test(na.omit(pval_lst_asy))
      return(list(pval = pval))
    }
  }
  else{
    # if n.perm in not null, select the significant lineage according to the resampling pvalue
    #  run test for each lineage
    pval_lst_asy = c()
    pval_lst_perm = c()
    for(i in 1:length(Y)){
      tmp = try(QCAT.Cluster(ID = ID[[i]], OTU = Y[[i]], OTU.base=NULL, X = X[[i]], X.index = X.index, Tax=NULL, 
                             min.depth=0, perm.type=perm.type, n.perm=n.perm, fdr.alpha=0.05, test= test))
      if(!("try-error" %in% class(tmp))){
        pval_lst_asy = append(pval_lst_asy,tmp$pval[1])
        pval_lst_perm = append(pval_lst_perm, tmp$pval[2])
      }
    }
    if(is.null(pval_lst_asy)){
      return(NA)
    }else{
      pval = c(F.test(na.omit(pval_lst_asy)),F.test(na.omit(pval_lst_perm)))
      return(list(pval = pval))
    }
  }
  
}

QCAT_C_Meta <- function(ID, OTU, OTU.base=NULL, X, X.index, Tax=NULL, min.depth=0, 
                        perm.type=NULL, n.perm=NULL, fdr.alpha=0.05, test="chisq",nthreads = 1){
  n.OTU = length(OTU)
  n.X = length(X)
  # drop.col = NULL
  
  if(n.OTU != n.X)
  {
    stop("The study number of OTU table and Covariate should be the same")
  }
  # if(length(unique(sapply(1:n.OTU,function(j) ncol(OTU[[j]]))))!= 1){
  #   stop("The taxa in each study should be the same")
  # }
  for(i in 1:n.OTU)
  {
    if(!is.matrix(OTU[[i]])){
      OTU[[i]] = as.matrix(OTU[[i]])
    }
    
    if(!is.matrix(X[[i]])){
      X[[i]] = as.matrix(X[[i]])
    }
    
    if(nrow(OTU[[i]])!=nrow(X[[i]])){
      stop(paste0("Number of samples in the OTU table and the covariate table of study ", i,
                  " should be the same\n"))
    }
    remove.subject = which(rowSums(OTU[[i]])<min.depth)
    if(length(remove.subject)>0){
      print(paste("Remove",length(remove.subject), "samples with read depth less than", min.depth, "in OTU table", i,"\n"))
      X[[i]] = X[[i]][-remove.subject, ,drop=FALSE]
      OTU[[i]] = OTU[[i]][-remove.subject, ,drop=FALSE]
      ID[[i]] = ID[[i]][-remove.subject, ,drop=FALSE]
    }
    # drop.col = union(drop.col,which(colSums(OTU[[i]])==0))
  }
  
  if(missing(X.index)){
    X.index = 1:ncol(X[[1]])
  }
  # join all the OTU table
  OTU.comb = as.data.frame(OTU[[1]])
  for (i in 2:n.OTU){
    OTU.comb = bind_rows(OTU.comb,as.data.frame(OTU[[i]]))
  }
  OTU.comb[is.na(OTU.comb)] <- 0 # substitute NA to zero (do not matter because those columns will delete during calculation)
  OTU.comb <- as.matrix(OTU.comb)
  taxa.names = colnames(OTU.comb)
  if(!is.null(Tax)){
    # preserve those columns which have taxonomy information
    col.save = intersect(rownames(Tax),taxa.names)
    OTU.comb = OTU.comb[,taxa.names %in% col.save]
    Tax = Tax[rownames(Tax) %in% col.save,]
    Tax = Tax[match(taxa.names, rownames(Tax)),]
  }
  # divide the combined OTU table according to the number of observations for each study
  batch.cnt <- unlist(lapply(OTU, function(x) nrow(x)))
  batch.cnt <- append(1,batch.cnt)
  batch.cnt <- cumsum(batch.cnt)
  count = list()
  for (i in 1:n.OTU) {
    count[[i]] <- OTU.comb[batch.cnt[i]:(batch.cnt[i+1]-1),]
  }
  # X = lapply(1:n.OTU,function(j) cbind(1, X[[j]])) # add the intercept term
  # X.index = X.index + 1
  
  if(is.null(Tax)){ # perform one test using all OTUs
    score.pvalue.beta = NULL
    for (i in 1:n.OTU){
      tmp.q1  = miLineage::QCAT_GEE.Cluster(ID[[i]], count[[i]], X = X[[i]], X.index = X.index, Z = X[[i]], Z.index = X.index, Tax = NULL,
                                        min.depth = 0, n.perm = n.perm, test = test)
      if(! "try-error" %in% class(tmp.q1)){
        if(is.null(n.perm)){
          score.pvalue.beta = append(score.pvalue.beta,tmp.q1$pval[1])
        }else{
          score.pvalue.beta = append(score.pvalue.beta,tmp.q1$pval[2])
        }
      }
    }
    pval = as.matrix(F.test(na.omit(score.pvalue.beta)))
    colnames(pval) = paste0("Asymptotic-","Fisher")
    return( list(pval=pval) )
  }else{ # perform tests for lineages
    
    if(!is.matrix(Tax)){
      Tax = as.matrix(Tax)
    }
    
    # if(length(drop.col)>0){
    #   tax = Tax[-drop.col,,drop = FALSE]
    # }else{
    #   tax = Tax
    # }
    for(i in 1:n.OTU)
    {
      if( sum(!(colnames(count[[i]]) %in% rownames(Tax)))>0 ){
        stop(paste0("Error: OTU IDs in OTU table ",i," are not consistent with OTU IDs in Tax table"))
      }
    }
    
    # n.rank = ncol(tax)
    # # merge the tax and count together for later partition and merge OTU table according to taxonomy information
    # W.data.list = lapply(1:n.OTU,function(j) data.table(data.frame(tax, t(count[[j]]))))
    # otucols = lapply(1:n.OTU,function(j) names(W.data.list[[j]])[-(1:n.rank)])
    # n.level = n.rank-1
    # 
    # subtree = NULL
    # pval = NULL
    
    n.rank = ncol(Tax)
    # merge the tax and count together for later partition and merge OTU table according to taxonomy information
    W.data.list = NULL
    W.data.list = lapply(1:n.OTU,function(j) data.table(data.frame(Tax, t(count[[j]]))))
    otucols = lapply(1:n.OTU,function(j) names(W.data.list[[j]])[-(1:n.rank)])
    
    n.level = n.rank-1
    
    subtree = NULL
    pval = NULL
    tax.tab = NULL
    for(k in 1:n.level){
      Rank.low = paste("Rank", n.rank-k,sep="")
      Rank.high = paste("Rank", n.rank-k+1,sep="")
      
      tmp = table(Tax[,n.rank-k])
      level.uni = sort( names(tmp)[which(tmp>1)] )
      m.level = length(level.uni)
      tmp = cbind(Rank.low, Rank.high,level.uni)
      tax.tab = rbind(tax.tab, tmp)
    }
    tax.info.lst <- lapply(seq_len(nrow(tax.tab)), function(i) tax.tab[i,])
    if(nthreads > 1){
    res_lst <- mclapply(tax.info.lst, lineage.QCAT.Cluster, W.data.list, otucols, X, X.index, ID, perm.type, n.perm, test, mc.cores = nthreads)
    for(j in 1:length(tax.info.lst)){
      if(any(is.na(res_lst[[j]]))|is.null(res_lst[[j]])){
        next
      }else{
        subtree = c(subtree, unname(tax.info.lst[[j]]["level.uni"]))
        pval = cbind(pval, res_lst[[j]]$pval)
      }
    }
  }else{
    for(j in 1:length(tax.info.lst)){
      # print(tax.info.lst[[j]])
      if(j == 5){
        print('stop here')
      }
      res_tmp = lineage.QCAT.Cluster(tax.info.lst[[j]], W.data.list, otucols, X, X.index, ID, perm.type, n.perm, test)
      if(any(is.na(res_tmp))|(is.null(res_tmp))){
        next
      }else{
        subtree = c(subtree, unname(tax.info.lst[[j]]["level.uni"]))
        pval = cbind(pval, res_tmp$pval)
      }
    }
  }
  colnames(pval) = subtree
  if(is.null(n.perm)){
    rownames(pval) = paste0("Asymptotic-", "Fisher")
    score.tmp = pval[1,]
  }else{
    rownames(pval) = c(paste0("Asymptotic-", "Fisher"), paste0("Resampling-", "Fisher"))
    score.tmp = pval[2,]
  }
  
  # identify significant lineages
  subtree.tmp = subtree
  index.na = which(is.na(score.tmp))
  if(length(index.na)>0){
    # drop those lineages which have NA values
    score.tmp = score.tmp[-index.na]
    subtree.tmp = subtree.tmp[-index.na]
  }
  
  #score.tmp[score.tmp==0] = 1e-4
  m.test = length(score.tmp)
  
  # Benjamini-Hochberg FDR control
  index.p = order(score.tmp)
  p.sort = sort(score.tmp)
  #fdr.alpha = 0.05
  
  # change 2022/12/01
  reject = rep(0, m.test)
  tmp = which(p.sort<=(1:m.test)*fdr.alpha/m.test)
  if(length(tmp)>0){
    index.reject = index.p[1:max(tmp)]
    reject[index.reject] = 1
  }
  
  sig.lineage = subtree.tmp[reject==1]
  
  # return all the p-values as well as significant lineages
  return( list(lineage.pval=pval, sig.lineage=sig.lineage) )
  }
}
