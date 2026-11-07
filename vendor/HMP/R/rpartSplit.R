rpartSplit <-
function(y, wt, x, parms, continuous){
	# Get initial LL
	LL <- DM.MoM(y)$loglik
	
	uniqX <- sort(unique(x))
	numUni <- length(uniqX) - 1
	
	# Determine what we are comparing
	if(continuous){
		numTests <- length(x) - 1
		dir <- rep(-1, numTests)
	}else{
		numTests <- numUni
		dir <- uniqX
	}
	
	# Run through every comparison
	LRT <- rep(0, numTests)
	for(i in 1:numUni){
		if(continuous){
			id <- which(x <= uniqX[i])
			grp1 <- y[id,, drop=FALSE]
			grp2 <- y[-id,, drop=FALSE]
		}else{
			grp1 <- y[x == uniqX[i],, drop=FALSE]
			grp2 <- y[x != uniqX[i],, drop=FALSE]
		}
		# Skip any 1 subject groups
		if(nrow(grp1) == 1 || nrow(grp2) == 1)
			next
		
		LLgrp1 <- DM.MoM(grp1)$loglik
		LLgrp2 <- DM.MoM(grp2)$loglik
		
		# Skip any infinite LL comparisons (makes lrt 0)
		if(LLgrp1 == Inf || LLgrp2 == Inf)
			next
		
		if(continuous){
			LRT[id[length(id)]] <- -2*(LL-LLgrp1-LLgrp2)
		}else{
			LRT[i] <- -2*(LL-LLgrp1-LLgrp2)
		}
	}
	ret <- list(goodness=LRT, direction=dir)
	
	return(ret)
}
