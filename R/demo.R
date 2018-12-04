#reformat the rawdata, chain the samples
rawData	<-  t(matrix(aperm(xIncl,c(2,1,3)),nrow=250,ncol=121*40))
summedChrom<-apply(rawData[1:121,],2,sum)

# example MassSpectra
barplot(summedChrom[35:80])

# calculate tic
tic <- apply(rawData,1,sum)

# composite plot
plot(tic[1:121], type='l', ylab='Signal', xlab='Scans')
for(i in seq(122,4840,121)){
  print(i)
  points(1:121,tic[i:(i+120)], type='l')
}

# run
results <- do_AR_all(xIncl, 1, 3, 45, 1:250)

# peaks for one sample
plot(results$C[1:121,1], type='l', ylim=c(0,500000), ylab='Signal', xlab='Scans')
for(i in 2:11){
  points(1:121, results$C[1:121,i], type='l')
}

# mass spectra
for(i in 1:11){
  barplot(results$S[35:80,i], xlab='mz', ylab='counts')
  readline(prompt="Press [enter] to continue")
}



