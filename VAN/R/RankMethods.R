## Class for calculating Rank based measures
## on a data frame containing ranking information

## Uses 'RankProd' in Version 2.20.0

## function for calculating rank product described by Breitling
## since this may take very long, we suggest calculating it 
## using parallel computing
## data: data.frame woth rownames being features, columns studies 
## or method amd entries being ranks
## B: number of permutations
## cluster: SNOW cluster object
## Value: data.frame
calculateRankProduct <- function(data, B=1000, cluster=NULL){
	cat("starting rank product ensemble approach\n")
	if(is.null(cluster)){
		# calculate row-wise rank product
		RP <- apply(data, 1, prod)/nrow(data)
	}else{
		# calculate row-wise rank product
		RP <- parApply(cluster, data, 1, prod)/nrow(data)
	}
	#vector containing RP p-values
	RP.p <- vector("numeric", nrow(data))
	
	if(is.null(cluster)){
		for(id in 1:nrow(data)){
			counter <- sapply(1:B, .doRP, data, RP[id])  
			RP.p[id] <- sum(counter)/B
			#output showing to be still alive
			if(id %% ceiling(nrow(data)/100) - 1 == 0){
				cat("=")
				flush.console()
			}
		}
	}else{
		for(id in 1:nrow(data)){
			counter <- parSapply(cluster, 1:B, .doRP, data, RP[id])		 
			RP.p[id] <- sum(counter)/B
			#output showing to be still alive
			if(id %% ceiling(nrow(data)/100) - 1 == 0){
				cat("=")
				flush.console()
			}
		}
	}
	
	cat("\n")
	flush.console()
	
	RP.q <- qvalue(RP.p)$qvalue
	RP.rank <- rank(RP.q, ties.method="random")
	return(data.frame(RP, p.value= RP.p, q.value= RP.q, rank=RP.rank, row.names=rownames(data)))
}

## function for sampling a random rank product 
## and comparing it to the expected one
## b: permutation counter (is ignored)
## data: data from calculateRankProduct
## RP: expected RP
## Value: 1 if random RP <= expected RP, 0 else
.doRP <- function(b, data, RP,cluster=NULL){
	if (is.null(cluster)) {
		rRP <- prod(apply(data, 2, sample, size=1))/nrow(data)
	}else{
		rRP <- prod(parApply(cluster,data, 2, sample, size=1))/nrow(data)
	}
	if(rRP <= RP){
		cbind(1)
	} else {
		cbind(0)
	}
}

## function for calculating rank sum 
## since this may take very long, we suggest calculating it 
## using parallel computing
## data: data.frame woth rownames being features, columns studies 
## or method amd entries being ranks
## B: number of permutations
## cluster: SNOW cluster object
## Value: data.frame
calculateRankSum <- function(data, B= 10000, cluster= NULL){
	cat("starting rank sum ensemble approach\n")
	if(is.null(cluster)){
		#calculate featurewise rank sum
		RS <- apply(data, 1, sum)/nrow(data)
	}else{
		RS <- parApply(cluster, data, 1, sum)/nrow(data)
	}
	#container for p-values
	RS.p <- vector("numeric", nrow(data))
	
	if(is.null(cluster)){
		for(id in 1:nrow(data)){
			counter <- sapply(1:B, .doRS, data, RS[id])
			RS.p[id] <- sum(counter)/B
			#output to prove being alive still
			if(id %% ceiling(nrow(data)/100) - 1 == 0){
				cat("=")
				flush.console()
			}
		}
	}else{
		for(id in 1:nrow(data)){
			counter <- parSapply(cluster, 1:B, .doRS, data, RS[id])
			RS.p[id] <- sum(counter)/B
			#output to prove being alive still
			if(id %% ceiling(nrow(data)/100) - 1 == 0){
				cat("=")
				flush.console()
			}
		}
	}
	
	cat("\n")
	flush.console()
	
	RS.q <- qvalue(RS.p)$qvalue
	RS.rank <- rank(RS.q, ties.method="random")
	return(data.frame(RS, p.value= RS.p, q.value= RS.q, rank=RS.rank, row.names=rownames(data)))
}

## function for sampling a random rank sum
## and comparing it to the expected one
## b: permutation counter (is ignored)
## data: data from calculateRankSum
## RS: expected RS
## Value: 1 if random RS <= expected RS, 0 else
.doRS <- function(b, data, RS,cluster=NULL){
	if (is.null(cluster)) {
		rRS <- sum(apply(data, 2, sample, size=1))/nrow(data)
	}else{
		rRS <- sum(parApply(cluster,data, 2, sample, size=1))/nrow(data)
	}
	if(rRS <= RS){
		cbind(1)
	} else {
		cbind(0)
	}
}
