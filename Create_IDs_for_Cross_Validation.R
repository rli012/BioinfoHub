nfold<-10
n <- 153

foldidList<-lapply(1:10,function(i){
  #n<-length(y)
  sample(rep(1:nfold,ceiling(n/nfold))[1:n])
})
