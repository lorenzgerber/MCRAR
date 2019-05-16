loadData <- function(file){
  data('win_1')
}


import_from_txt <- function(filename){

  con <- file(filename,"r")
  first_line <- readLines(con,n=1)
  close(con)


}

importDatForR <- function(fileName, mzs, samples, scans ){

  inputMatrix <- read.table(file=fileName)
  dataArray <- array(0, dim=c(scans, length(mzs), samples))

  for(i in 1:scans){
    for(j in 1:samples){
      for(k in mzs){
        dataArray[i,k,j] <- inputMatrix[ (j-1) * length(mzs) + k, i ]
      }
    }
  }
  return(dataArray)
}
### Format for java is a textfile
importDatForJava <- function(rileName, mzs, samples, scans){

  inputMatrix <- read.table(file=fileName)
  outputMatrix <- matrix(0, length(mzs)*scans, samples)

  for(i in 1:scans){
    for(j in 1:samples){
      for(k in mzs){
        outputMatrix[ (k-1)*scans + i , j ] <- inputMatrix[ (j-1) * length(mzs) + k, i ]
      }
    }
  }
  return(outputMatrix)
}


generateFiles <- function(rootFileName, mzs, samples, scans){

  pathFileNameData <- paste('./Model samples/', rootFileName, '.dat', sep='')
  pathFileNameBG <- paste('./Model samples_bg_corr/bg_', rootFileName, '.Rdata', sep ='')
  inputMatrix <- read.table(file=pathFileNameData)
  load(pathFileNameBG)

  outputMatrix <- matrix(0, length(mzs)*scans, samples)
  dataArray <- array(0, dim=c(scans, length(mzs), samples))

  ### Transforming the input data to match the background correction table
  for(i in 1:scans){
    for(j in 1:samples){
      for(k in mzs){
        #outputMatrix[ (k-1)*scans + i , j ] <- inputMatrix[ (j-1) * length(mzs) + k, i ]
        dataArray[i,k,j] <- inputMatrix[ (j-1) * length(mzs) + k, i ]
      }
    }
  }

  dataArray <- dataArray - BL
  dataArray[dataArray < 0] <- 0
  outputMatrix <- matrix(dataArray, scans*length(mzs), samples)

  fileNameRData <- paste(rootFileName, ".Rdata", sep="")
  fileNameTxt <- paste(rootFileName, ".txt", sep="")

  save(dataArray, file=fileNameRData)
  write.table(outputMatrix, file=fileNameTxt, row.names=FALSE, col.names=FALSE, sep="\t")

}





##' Function do_AR_all
##'
##' Function do_AR_all
##' Funciton to do Alternate Regression
##' @param x
##' @param projectpath
##' @param NL
##' @param RP
##' @param RT_LIMIT
##' @param SCAN_RANGE
##' @return C, S, INDEX, R2, noise, CP
do_AR_all <- function(x,NL,RP,RT_LIMIT,SCAN_RANGE){
  require(MASS)

  var <-  dim(x)[1]
  mz  <-  dim(x)[2]
  obs <-  dim(x)[3]
  cat("\n[TIMEPOINTS,M/Z-Channels,Observations] = [",var,",",mz,",",obs,"]\n\n")
  R					<-	OK	<-	1
  R2					<-	0
  C					<-	S	<-	numeric()
  OK_out				<-	3
  X 					<-  t(matrix(aperm(x,c(2,1,3)),nrow=mz,ncol=var*obs))
  ind					<-	w_mor2(X)
  noise				<-	which(ind<NL)
  x[,noise,]			<-	0
  ssX					<-  sum(apply(x,3,ss))
  Xmean				<-	apply(x,2,function(x) apply(x,1,sum))/dim(x)[3]
  rm(x)
  X[,noise]		<-	0

  # estimate how many components there might be be
  max_R				<-	qr(X[,setdiff(1:length(ind),noise)])$rank

  S						<-	abs(pca(X,1)$p)
  vars				<-	which(as.logical(SCAN_RANGE))
  cat("-----------------------------\n")
  gc()

  while(OK_out>0){
    iter	<-	0
    C 		<-  numeric()

    # proceed if number of components does not estimate number
    # of estimated components
    if(max_R >= R){
      # if there is more than one component, estimate the mass
      # mass profiles from Xmean using pure algo
      if(R > 1){
        p		<-	pure(Xmean[,vars],R,0.01)
        sp  <-  p$sp
        imp <-  p$imp
        S	<-	matrix(0,mz,R)

        for(i in 1:R)
          S[vars[imp[i]],i]	<-	1
      }

      R2			<-	0.001
      I			<-	numeric()
      dif			<-	1
      C			<-	X%*%S%*%ginv(t(S)%*%S,tol=.Machine$double.eps^2)

      C[C<0]		<-	0
      out			<-	unimodal2(C,R,var,obs,RT_LIMIT)
      C	  		<-  out$C
      CP			<-  out$CP

      while(dif > 1e-6){
        r2			<-	R2
        iter		<-	iter+1
        tC  		<-  t(C)
        S			<-	t(ginv(tC%*%C,tol=.Machine$double.eps^2)%*%(tC%*%X))
        S[S<0]		<-	0
        S 			<-  apply(S,2,function(S) S/sum(S))

        if(any(is.na(S))){
          break
        }

        C	<-	X%*%S%*%ginv(t(S)%*%S,tol=.Machine$double.eps^2)

        if(any(is.na(C)))
          break
        C[C<0]	<-	0
        out			<-	unimodal2(C,R,var,obs,RT_LIMIT)
        C       <-  out$C
        CP      <-  out$CP
        R2		<-	1-ss(X-C%*%t(S))/ssX
        #				cat(iter,": ",R2,"| sum(C) = ",sum(C),"\n")

        if(iter == 50)
          dif <-  0

        else
          dif <-  abs(R2-r2)
      }

      if(any(is.na(S)) | any(is.na(C))){

        cat("Warning: NAs found.\n")
        C		<-	matrix(1,nrow(C),ncol(C))
        S		<-	matrix(0,nrow(S),ncol(S))
        R2		<-  0
      }

      INDEX	<-	matrix(0,obs,R)

      for(i in 1:R){
        c	<-	matrix(C[,i],var,obs)
        for(j in 1:obs){
          if(max(c[,j])>0)
            INDEX[j,i]  <-  min(which.max(c[,j]))
          else
            INDEX[j,i]  <-  -99
        }
      }

      if(sum(INDEX == -99) > 0){
        I <-  row(INDEX)[INDEX == -99]
        J <-  col(INDEX)[INDEX == -99]

        for(i in 1:length(J)){
          index				<-	INDEX[,J[i]]
          index				<-	index[-(index == -99)]
          INDEX[I[i],J[i]]	<-	median(index)
        }
      }

      I		<-  order(colMeans(INDEX))
      Y		<-  sort(colMeans(INDEX))
      C		<-	C[,I]
      S		<-	S[,I]
      CP		<-	CP[I]
      INDEX	<-	INDEX[,I]
      cat("Rank: ", R,"\n")

      if(R>1)
        PERMUTED_ORDER	<-	sum(diff(t(INDEX)) < -RP)
      else
        PERMUTED_ORDER	<-	0

      cat("PERMUTED_ORDER = ",PERMUTED_ORDER,"\n")
      cat("Number of iterations: ", iter,"\n")
      cat("R2X = ",round(R2,5),"\n")

    }else{
      R2				<-	0
      S				<- 	C 	<-  numeric()
      PERMUTED_ORDER	<-	99
    }

    s		<-	S
    c		<-	C
    r2		<-	R2

    if(r2>0){
      if(!PERMUTED_ORDER){
        C_out		<-	c
        S_out		<-	s
        CP_out		<-	CP
        R2_out		<-	r2
        OK_out		<-	3
        OK2_out		<-	0
        cat("OK\n")

      }else{
        OK2_out	<-	1
        cat("NOT OK [1]\n")
      }

    }else{
      OK2_out	<-	1
      cat("NOT OK [2]\n")
    }

    OK_out	<-	OK_out-OK2_out
    R		<-	R+1

    cat("Timestamp: ",format(Sys.time(), "%X"),"\n")
    cat("-----------------------------\n")
  }

  if(!exists("C_out")){
    C_out		<- numeric()
    S_out		<-	numeric()
    R2_out	<-	numeric()
    CP_out	<-	numeric()
  }

  C		<-	as.matrix(C_out)
  S		<-	as.matrix(S_out)
  R2		<-	R2_out
  CP		<-	CP_out
  R		<-	ncol(C)

  cat("FINAL MODEL\n")
  cat("Rank:\t",R,"\n")
  cat("R2X:\t",R2,"\n")
  cat("-----------------------------")


  if(!is.null(R)){ #produces errors when R=1
    #if(R>1){
    INDEX	<-	matrix(0,obs,R)

    for(i in 1:R){
      c	<-	matrix(C[,i],var,obs)
      for(j in 1:obs)
        INDEX[j,i]	<-	min(which.max(c[,j]))

    }

    I		<-	order(round(colMeans(INDEX)))
    INDEX 	<-  sort(round(colMeans(INDEX)))
  }else
    INDEX	<-	numeric()

  gc()
  out <-  list(C=C,S=S,INDEX=INDEX,R2=R2,noise=noise,CP=CP)
}

##' Function pure
##'
##' Function pure
##' [sp,imp]=pure(d,nr,f)
##' sp purest row/column profiles
##' imp indexes of purest variables
##' d data matrix; nr (rank) number of pure components to search
##' if d(nspectra,nwave) imp gives purest nwave => sp are conc. profiles (nr,nspectra)
##' if d(nwave,nspectra) imp gives purest nspectra => sp are spectra profiles (nr,nwave)
##' f percent of noise allowed respect maximum of the average spectrum given in % (i.e. 1% or 0.1%))
##' @param d
##' @param nr
##' @param f
##' @return sp, imp
pure <-function(d,nr,f){

  #[nrow,ncol]=size(d);
  # calculation of the purity spectrum
  w   		<-  p	<-	s	<-	matrix(0,nrow=max(nr,1),ncol=ncol(d))
  f			<-	f/100
  s[1,]		<-	apply(d,2,sd)
  m			<-	colMeans(d)
  ll			<-	s[1,]^2+m^2
  f			<-	max(m)*f
  p[1,]		<-	s[1,]/(m+f)

  imp 		<-  numeric()
  mp 			<-  max(p[1,])
  imp[1]  <-  which.max(p[1,])

  # calculation of the correlation matrix
  l	<-	sqrt(s[1,]^2+(m+f)^2)
  for(i in 1:nrow(d))
    d[i,]	<-	d[i,]/l

  c	<-	(t(d)%*%d)/nrow(d)

  # calculation of the weights
  # first weight
  w[1,]	<-	ll/l^2
  p[1,]	<-	w[1,]*p[1,]
  s[1,]	<-	w[1,]*s[1,]

  # next weights
  if(nr > 1){
    for(i in 2:nr){
      for(j in 1:ncol(d)){
        dm			<-	wmat(c,imp,i,j)
        w[i,j]	<-	det(dm)
        p[i,j]	<-	p[1,j]*w[i,j]
        s[i,j]	<-	s[1,j]*w[i,j]
      }

      mp[i] 	<-  max(p[i,])
      imp[i]  <-  which.max(p[i,])
    }
  }

  sp	<-	normv2(t(d[,imp[1:nr]]))
  purelist  <-  list(sp=sp,imp=imp)
}

##' Function normv2
##'
##' Functino normv2
##' @param s
##' @return sn
normv2 <-function(s){

  if(nrow(s) == 1){
    sn  <-  s/sqrt(sum(s^2))

  }else
    sn  <-  apply(s,1,function(s) s/sqrt(sum(s^2)))
}

##' Function wmat
##'
##' Function wmat
##' @param tempmat
##' @param imp
##' @param irank
##' @param jvar
##' @return dm
wmat<-function(tempmat,imp,irank,jvar){
  dm								<-	matrix(0,irank,irank)
  dm[1,1]							<-	tempmat[jvar,jvar]
  dm[1,2:irank]					<-	tempmat[jvar,imp[1:(irank-1)]]
  dm[2:irank,1]					<-	tempmat[imp[1:(irank-1)],jvar]
  dm[2:irank,2:irank] <-  tempmat[imp[1:(irank-1)],imp[1:(irank-1)]]

  dm


}

##' Function unimodal2
##'
##' Function unimodal2
##' @param C
##' @param R
##' @param var
##' @param obs
##' @param RT_LIMIT
##' @return C, CP
unimodal2<-function(C,R,var,obs,RT_LIMIT){

  N 	<-  nrow(C)
  CP  <-  numeric(R)
  for(i in 1:R){
    c       <-  matrix(C[,i],var,obs)
    out		<-	unimodal3(c,obs,RT_LIMIT)
    C[,i]	<-	matrix(out$C,N,1)
    CP[i]	<-	out$cp
  }
  out <-  list(C=C,CP=CP)
}

##' Function w_mor2
##'
##' Function w_mor2
##' @param X
##' @return ind
w_mor2<-function(X){

  X								<-	sweep(X,2,colMeans(X))
  ind								<-  apply(X,2,function(X) sqrt(sum(X^2))/sqrt(sum(diff(X)^2)))
  ind[!is.finite(ind)]  <-  0
  ind
}

##' Function unimodal3
##'
##' Function unimodal3
##' @param C
##' @param obs
##' @param RT_LIMIT
##' @return C, CP
unimodal3<-function(C,obs,RT_LIMIT){

  ### cp = center of peak, Median index from max values of each chromatographic
  ### profile in C
  cp  <-	round(median(which.max(colMeans(t(C))))) # cp is always only 1 long !!!! This could be wrong in the R version !!!!
  if(!length(cp)){
    C		<-	C*0
    cp		<-	0

  }else{
    for(i in 1:obs){
      xpeak	<-	peak_pick(t(C[,i]))$xpeak

      ## get the indices of maximum peaks (better local maxima) found in peak_pick
      mp		<-	which(as.logical(xpeak)) # 0 to some (maybe 6 -7), mp are indices

      ## Check if there is a local maxima at all, else fill the current observation
      ## (C row) with zero's and process the next observation
      if(length(mp)){

        ## find the index of the local maxima the closest to cp, the center of peak
        ## reset mp to contain only this value
        if(which.min(abs(mp-cp))){
          mp	<-	mp[which.min(abs(mp-cp))]
        }

        ## Check whether the difference between cp (center of peak) and the closest
        ## local maxima is within the limit of Retention time precision (RT_LIMIT)
        ## else fill the current observation (C vector row) with zero's and process
        ## the next observation
        if(abs(mp-cp)<RT_LIMIT){

          #------Remove local maxima before the absolut maxima-------

          # calculate the differences D from the begin of the chromatographic profile until
          # the found local maxima mp.
          D			<-	diff(C[1:mp,i])

          # Are there any differences below 0 if so, put highest index
          # of such an occurence into 'poss'. Each such index represents
          # the first from the number number which result in a negative
          # difference
          if(any(D<0)){
            poss	<-	max(which(D<0))
          } else {
            poss	<-	1
          }


          # For each index j which runs from current poss down to zero/one
          # rewrite the corresponding index in C by the minimum of values
          # between j and mp in C, current sample row.
          for(j in poss:1)
            C[j,i]	<-	min(C[j:mp,i])

          #---------Remove local maxima after the  absolut maxima
          # calculate differences from the maxpeak (mp) to the end
          D			<-	diff(C[mp:nrow(C),i])

          # is there any difference in D that is larger than zero
          # if so find the smallest index (in case there are several such)
          # where D is larger than zero and store it in 'poss'
          if(length(which(D>0)))
            poss	<-	min(which(D>0))
          else
            poss	<-	FALSE

          # if increasing signal after peak top available
          # rewrite for each index j which runs from maxpeak (mp) minus 1
          # plus poss until the end by the minimum of values
          # between mp until j in C, current sample row.
          if(poss)
            for(j in (mp-1+poss):nrow(C))
              C[j,i]	<-	min(C[mp:j,i])  ## Crash in JAVA version



          # Correction of strange peaks (most likely to happen when above
          # corrections were applied -> negative diff before mp, positive diff before mp)
          # find if there is any zero diff in C for current sample
          zeroDiff <- which(diff(C[,i]) == 0)
          # write signal for current sample except indices with zeroDiff into knew (k new)
          # else write the whole signal into knew
          if(length(zeroDiff)) {
            knew <- C[-zeroDiff,i]
          } else {
            knew <- C[,i]
          }

          # from knew find the index of the max value in knew
          mpnew						<-	which.max(knew)

          # make new zero vector same length as signal vector
          k							<-	C[,i]*0

          # fit knew into k from index one plus mp minus mpnew until index length of knew plus
          # the difference between mp and mpnew
          k[1:length(knew)+mp-mpnew]	<-	knew
          # write k vector into C for current sample
          C[,i]	<-	k

        }else
          C[,i]	<-	C[,i]*0
      }else
        C[,i]	<-	C[,i]*0
    }
  }
  out <-  list(C=C,cp=cp)
}

##' Function peak_pick
##'
##' Function peak_pick
##' Function to find peaks. After applying a savitzky-golay smooth,
##' local maxima are detected. As condition, the 2 scans before
##' and after the current scan are tested for being smaller than
##' the current scan. As noise filter, the median value of the
##' all signals from the chromatographic profile/sample under
##' investigation is used. For sequential true conditions with a
##' difference of less than 3 in signal, the smaller true
##' condition is removed. 'xout' is a vector of the same length
##' as the chromatographic profile x with zero values except for
##' the found local maxima for which respective values of x are shown.
##'
##' @param x
##' @return xpeak, xout
peak_pick<-function(x){

  xout				<-	x	<-	as.matrix(x)

  xd1         <-  sgolayfilt(xout,3,11,1)
  NOISE_LIMIT	<-	median(x)
  N1					<-	which(x>NOISE_LIMIT)
  N2					<-	matrix(0,nrow(x),ncol(x))

  for (i in 3:(ncol(x)-2))
    if(xd1[i-2] > 0)
      if(xd1[i-1] > 0)
        if(xd1[i+1] < 0)
          if(xd1[i+2] < 0)
            if(sum(N1 == i) == 1)
              N2[i]	<-	TRUE

  N	<-	which(as.logical(N2))

  if(length(N) > 1){

    while(min(diff(N)) < 3){
      p1	<-	min(which(diff(N) < 3))
      p2	<-	p1+1

      if(xout[N[p1]] < xout[N[p2]])
        N	<-	N[-p1]
      else
        N	<-	N[-p2]
      if(length(N) == 1)
        break
    }
  }
  xpeak			<-	matrix(0,nrow(x),ncol(x))
  xpeak[N]	<-	xout[N]
  out <-  list(xpeak = xpeak, xout = xout)
}

##' Function sgolayfilt
##'
##' Function sgolayfilt
##' @param x
##' @param p
##' @param n
##' @param m
##' @param ts
sgolayfilt<-function(x, p = 3, n = p + 3 - p%%2, m = 0, ts = 1){

  len = length(x)
  if (class(p) == "sgolayFilter" || (!is.null(dim(p)) && dim(p) > 1)){
    F = p
    n = nrow(F)
  }else
    F = sgolay(p, n, m, ts)
  k = floor(n/2)

  #z = filter(F[k+1,n:1], 1, x)

  filt  <-  F[k+1,n:1]
  z 		<-	na.omit(stats:::filter(c(rep(0,length(filt) - 1), x), filt, sides = 1))
  c(F[1:k,] %*% x[1:n], z[n:len], F[(k+2):n,] %*% x[(len-n+1):len])
}

##' Function sgolay
##'
##' Function sgolay
##' @param p
##' @param n
##' @param ts
##' @return F
sgolay<-function(p, n, m = 0, ts = 1){

  library(MASS)
  if (n %% 2 != 1)
    stop("sgolay needs an odd filter length n")

  if (p >= n)
    stop("sgolay needs filter length n larger than polynomial order p")

  F = matrix(0., n, n)
  k = floor(n/2)
  for (row  in  1:(k+1)) {
    C = ( ((1:n)-row) %*% matrix(1, 1, p+1) ) ^ ( matrix(1, n) %*% (0:p) )
    A = ginv(C)
    F[row,] = A[1+m,]
  }

  F[(k+2):n,] = (-1)^m * F[k:1,n:1]

  if (m > 0)
    F = F * prod(1:m) / (ts^m)

  class(F) = "sgolayFilter"
  F
}

##' Function ss
##'
##' Function ss
##' @param X
##' @return SSX
ss<-function(X)
  SSX	<-	sum(X^2)

##' Function pca
##'
##' Function pca
##' @export
##' @param X
##' @param comp
##' @return vec, p
pca<-function(X,comp){

  #   Calculates Principal componets of X by calc eigvectors of X'*X or X*X'
  #   Depending on whats easiest to calculate....

  if(nrow(as.matrix(X)) < ncol(as.matrix(X))){
    tX  <-  t(X)
    p 	<-  eigen(X%*%tX)$vectors[,1:comp]

    if(comp>1){
      p <-  apply(p,2,function(p) tX%*%p)
      p <-  apply(p,2,function(p) p/sqrt(sum(p^2)))

    }else{
      p	<-	tX%*%p
      p 	<-	p/sqrt(sum(p^2))
    }

  }else
    p		<-	eigen(t(X)%*%X)$vectors[,1:comp]

  vec		<-	X%*%p
  vecp	<-	list(vec=vec,p=p)
}
