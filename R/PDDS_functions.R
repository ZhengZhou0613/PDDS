GLP<-function(n,p,type="CD2"){
  if((type %in% c("CD2","WD2","MD2"))==F){
    stop("The type shoud be chose from 'CD2','WD2', or 'MD2'.")
  }
  fb<-c(3,5,8,13,21,34,55,89,144,233,377,610,987,1597)
  if(((n+1)%in%fb)&(p==2)){
    design0<-matrix(0,(n+1),p)
    H<-rep(1,2)
    H[2]<-fb[which(fb==(n+1))-1]
    for (j in 1:p) {
      for (i in 1:(n+1)) {
        design0[i,j]<-(2*i*H[j]-1)/(2*(n+1))-floor((2*i*H[j]-1)/(2*(n+1)))
      }
    }
    design0<-design0[-(n+1),]*(n+1)/n
  }else{
    if(p==1){
      design0<-matrix(0,n,p)
      for(i in 1:n){
        design0[i,1]<-(2*i-1)/(2*n)
      }
      return(design0)
    }
    h<-c()
    for(i in 2:min((n+1),200)){
      if(numbers::coprime((n+1),i)==T){
        h<-c(h,i)
      }
    }
    if(p>2){
      for (i in 1:100) {
        if(choose(p+i,i)>5000){
          addnumber<-i
          break
        }
      }
      h<-h[sample(1:length(h),min(length(h),(p+addnumber)))]
    }
    H<-utils::combn(h,p,simplify = F)
    if(length(H)>3000){
      H<-H[sample(3000)]
    }
    design0<-matrix(0,n,p)
    d0<-UniDOE::DesignEval(design0,crit=type)
    for (t in 1:length(H)) {
      design<-matrix(0,n,p)
      for (i in 1:p) {
        for (j in 1:n) {
          design[j,i]<-(j*H[[t]][i])%%(n+1)
        }
      }
      d1<-UniDOE::DesignEval(design,crit=type)
      if(d1<d0){
        d0<-d1
        design0<-design
      }
    }
    design0<-(design0*2-1)/(2*n)
  }
  return(design0)
}

DDS<-function(data,n,type="CD2",ratio=0.85){
  if((type %in% c("CD2","WD2","MD2"))==F){
    stop("The type shoud be chose from 'CD2','WD2', or 'MD2'.")
  }
  if(n>dim(data)[1]){
    stop("The subsample size must be less than the sample size of the data.")
  }
  x<-data
  y<-n
  N<-dim(data)[1]
  p0<-dim(data)[2]
  svddata<-svd(scale(x))
  lambda<-(svddata$d)^2
  sumlambda<-sum(lambda)
  for (i in 2:p0) {
    if(sum(lambda[1:i])>(ratio*sumlambda)){
      p<-i
      print(paste("Take ratio =",ratio,",","p after pca =",i))
      break
    }
  }
  rdata<-(svddata$u)[,1:p]
  design<-GLP(n,p,type)
  yita<-matrix(nrow=n,ncol=p)
  for (i in 1:p) {
    for (j in 1:n) {
      yita[j,i]<-stats::quantile(rdata[,i],design[j,i])
    }
  }
  kdtree<-RANN::nn2(data=rdata,query=yita,k=1,treetype="kd")
  subsample<-kdtree$nn.idx
  return(subsample)
}

seqdimchoice<-function(data,partperd,samplenumber=prod(partperd)*10*dim(data)[2],corenumber=2){
  spdata<-function(data,sp,partperd){
    N<-dim(data)[1]
    p<-dim(data)[2]
    for (i in 1:p) {
      small<-min(abs(data[which(abs(data[,i])!=0),i]))
      data[,i]<-data[,i]+small*0.000001*stats::rnorm(N,0,1)
    }
    data<-as.matrix(data)
    part<-prod(partperd)
    data<-cbind(c(1:N),data)
    N<-floor(N/part)*part
    data<-data[1:N,]
    data<-data[order(data[,sp[1]+1]),]
    for (i in 2:length(partperd)) {
      partnow<-prod(partperd[1:(i-1)])
      for (j in 1:prod(partperd[1:(i-1)])) {
        data1<-data[((j-1)*(N/partnow)+1):(j*(N/partnow)),]
        data1<-data1[order(data1[,sp[i]+1]),]
        data[((j-1)*(N/partnow)+1):(j*(N/partnow)),]<-data1
      }
    }
    return(data)
  }
  if(dim(data)[2]<corenumber){
    stop("The corenumber must be no larger than the number of columns of data.")
  }
  if(length(partperd)>dim(data)[2]){
    stop("length(partperd) must be no more than the dimensions of the data.")
  }
  if((samplenumber%%prod(partperd))!=0){
    stop("samplenumber%%prod(partperd) must be 0.")
  }
  if(samplenumber>dim(data)[1]){
    stop("samplenumber must be less than dim(data)[1].")
  }
  rs<-sample(1:dim(data)[1],samplenumber)
  data<-scale(data[rs,])
  N<-dim(data)[1]
  p<-dim(data)[2]
  for (i in 1:p) {
    small<-min(abs(data[which(abs(data[,i])!=0),i]))
    data[,i]<-data[,i]+small*0.000001*stats::rnorm(dim(data)[1],0,1)
  }
  svddata<-svd(scale(data[sample(1:min(N,3000)),]))
  lambda<-(svddata$d)^2
  sumlambda<-sum(lambda)
  for (i in 2:p) {
    if(sum(lambda[1:i])>(0.6*sumlambda)){
      pcadim<-i
      break
    }
  }
  choicedim<-rep(0,length(partperd))
  cl <- parallel::makeCluster(corenumber)
  parallel::clusterSetRNGStream(cl, sample(1:10000,1))
  parallel::clusterExport(cl,c("data","N","p","partperd","pcadim","corenumber","spdata"),envir = environment())
  hzlist<-parallel::parLapply(cl,1:corenumber,function(k){
    hz<-rep(Inf,p)
    if(k>p%%corenumber){
      pindex<-c(((k-1)*(p%/%corenumber)+1):(k*(p%/%corenumber)))
    }else{
      pindex<-c(((k-1)*(p%/%corenumber)+1):(k*(p%/%corenumber)),corenumber*(p%/%corenumber)+k)
    }
    for (i in 1:length(pindex)) {
      testdata<-data[order(data[,pindex[i]]),]
      testvalue<-rep(0,partperd[1])
      for(s in 1:partperd[1]){
        testvalue[s]<-MVN::mvn(svd(scale(testdata[((s-1)*N/partperd[1]+1):(s*N/partperd[1]),]))$u[sample(1:(N/partperd[1]),min(1000,N/partperd[1])),1:pcadim],mvnTest = "royston",scale=T)$multivariateNormality$H
      }
      hz[pindex[i]]<-sum(testvalue)
    }
    hz
  })
  parallel::stopCluster(cl)
  hz<-rep(Inf,p)
  for (i in 1:p) {
    for (k in 1:corenumber) {
      if(hzlist[[k]][i]!=Inf){
        hz[i]<-hzlist[[k]][i]
      }
    }
  }
  choicedim[1]<-which.min(hz)
  data<-data[order(data[,choicedim[1]]),]
  cl <- parallel::makeCluster(corenumber)
  parallel::clusterSetRNGStream(cl, sample(1:10000,1))
  parallel::clusterExport(cl,c("data","N","p","partperd","pcadim","corenumber","spdata"),envir = environment())
  hzlist<-parallel::parLapply(cl,1:corenumber,function(k){
    hz<-rep(Inf,p)
    if(k>p%%corenumber){
      pindex<-c(((k-1)*(p%/%corenumber)+1):(k*(p%/%corenumber)))
    }else{
      pindex<-c(((k-1)*(p%/%corenumber)+1):(k*(p%/%corenumber)),corenumber*(p%/%corenumber)+k)
    }
    partnow<-partperd[1]*partperd[2]
    for (i in 1:length(pindex)) {
      testdata<-data
      for (j in 1:partperd[1]) {
        testdata[((j-1)*N/partperd[1]+1):(j*N/partperd[1]),]<-testdata[((j-1)*N/partperd[1]+1):(j*N/partperd[1]),][order(data[((j-1)*N/partperd[1]+1):(j*N/partperd[1]),pindex[i]]),]
      }
      testvalue<-rep(0,partnow)
      for(s in 1:(partnow)){
        testvalue[s]<-MVN::mvn(svd(scale(testdata[((s-1)*N/partnow+1):(s*N/partnow),]))$u[sample(1:(N/partnow),min(1000,N/partnow)),1:pcadim],mvnTest = "royston",scale=T)$multivariateNormality$H
      }
      hz[pindex[i]]<-sum(testvalue)
    }
    hz
  })
  parallel::stopCluster(cl)
  hz<-rep(Inf,p)
  for (i in 1:p) {
    for (k in 1:corenumber) {
      if(hzlist[[k]][i]!=Inf){
        hz[i]<-hzlist[[k]][i]
      }
    }
  }
  hz[choicedim[1]]<-NaN
  choicedim[2]<-which.min(hz)
  if(length(partperd)>2){
    for (l in 3:length(partperd)) {
      data<-spdata(data,choicedim[1:(l-1)],partperd[1:(l-1)])[,-1]
      cl <- parallel::makeCluster(corenumber)
      parallel::clusterSetRNGStream(cl, sample(1:10000,1))
      parallel::clusterExport(cl,c("data","N","p","partperd","pcadim","corenumber","spdata"),envir = environment())
      hzlist<-parallel::parLapply(cl,1:corenumber,function(k){
        hz<-rep(Inf,p)
        if(k>p%%corenumber){
          pindex<-c(((k-1)*(p%/%corenumber)+1):(k*(p%/%corenumber)))
        }else{
          pindex<-c(((k-1)*(p%/%corenumber)+1):(k*(p%/%corenumber)),corenumber*(p%/%corenumber)+k)
        }
        partprev<-prod(partperd[1:(l-1)])
        partnow<-prod(partperd[1:l])
        for (i in 1:length(pindex)) {
          testdata<-data
          for (j in 1:partprev) {
            testdata[((j-1)*N/partprev+1):(j*N/partprev),]<-testdata[((j-1)*N/partprev+1):(j*N/partprev),][order(data[((j-1)*N/partprev+1):(j*N/partprev),pindex[i]]),]
          }
          testvalue<-rep(0,partnow)
          for(s in 1:(partnow)){
            testvalue[s]<-MVN::mvn(svd(scale(testdata[((s-1)*N/partnow+1):(s*N/partnow),]))$u[sample(1:(N/partnow),min(1000,N/partnow)),1:pcadim],mvnTest = "royston",scale=T)$multivariateNormality$H
          }
          hz[pindex[i]]<-sum(testvalue)
        }
        hz
      })
      parallel::stopCluster(cl)
      hz<-rep(Inf,p)
      for (i in 1:p) {
        for (k in 1:corenumber) {
          if(hzlist[[k]][i]!=Inf){
            hz[i]<-hzlist[[k]][i]
          }
        }
      }
      hz[choicedim[1:(l-1)]]<-NaN
      choicedim[l]<-which.min(hz)
    }
    }
  return(choicedim)
}

PDDS<-function(data,n,partperd,sp=seqdimchoice(data,partperd),move="T",corenumber=2){
  part<-prod(partperd)
  if((n%%part)!=0){
    stop("n%%prod(partperd) must be 0.")
  }
  if(part<corenumber){
    stop("corenumber must be no larger than prod(partperd).")
  }
  N<-dim(data)[1]
  p<-dim(data)[2]
  data<-as.matrix(data)
  if(length(partperd)>dim(data)[2]){
    stop("length(partperd) must be no more than the dimensions of the data.")
  }
  if(max(sp)>dim(data)[2]){
    stop("The dimensions in sp is out of the range.")
  }
  if(n>dim(data)[1]){
    stop("The subsample size must be less than the sample size of the data.")
  }
  for (i in 1:length(sp)) {
    small<-min(abs(data[which(abs(data[sample(1:N,min(1000,N)),sp[i]])!=0),sp[i]]))*0.000001
    data[,sp[i]]<-data[,sp[i]]+small*stats::rnorm(N,0,1)
  }
  svddata<-svd(scale(data[sample(1:min(N,3000)),]))
  lambda<-(svddata$d)^2
  sumlambda<-sum(lambda)
  for (i in 2:p) {
    if(sum(lambda[1:i])>(0.6*sumlambda)){
      p<-i
      break
    }
  }
  data<-cbind(c(1:N),data)
  N<-floor(N/part)*part
  data<-data[1:N,]
  data<-data[order(data[,sp[1]+1]),]
  for (i in 2:length(partperd)) {
    partnow<-prod(partperd[1:(i-1)])
    for (j in 1:prod(partperd[1:(i-1)])) {
      data1<-data[((j-1)*(N/partnow)+1):(j*(N/partnow)),]
      data1<-data1[order(data1[,sp[i]+1]),]
      data[((j-1)*(N/partnow)+1):(j*(N/partnow)),]<-data1
    }
  }
  subsample<-matrix(0,part,n/part)
  if((n/part)>(8)){
    design<-GLP((n/part),p,"CD2")
  }
  if(((n/part)<=8)&((n/part)>2)){
    design<-((UniDOE::GenUD((n/part),p,(n/part),crit="CD2")$final_design)*2-1)/(2*n/part)
  }
  if((n/part)==2){
    design<-matrix(0,2,p)
    design[1,]<-rep(0.25,p)
    design[2,]<-rep(0.75,p)
  }
  if((n/part)==1){
    design<-t(as.matrix(rep(0.5,p)))
  }
  design0<-design
  if(move=="T"){
    movedesign<-GLP(part,p,"CD2")
    movedesign<-movedesign[sample(1:part,part),]
  }else{
    movedesign<-matrix(0,part,p)
  }
  cl <- parallel::makeCluster(corenumber, type='PSOCK')
  parallel::clusterSetRNGStream(cl, sample(1:10000,1))
  parallel::clusterExport(cl,c("data","N","p","part","corenumber","design"),envir = environment())
  samplelist<-parallel::parLapply(cl,1:corenumber,function(k){
    subsample<-c()
    if(k<=part%%corenumber){
      partindex<-c(((k-1)*(part%/%corenumber)+1):(k*(part%/%corenumber)),corenumber*(part%/%corenumber)+k)
    }
    if(k>part%%corenumber){
      partindex<-c(((k-1)*(part%/%corenumber)+1):(k*(part%/%corenumber)))
    }
    for (i in 1:length(partindex)) {
      designstar<-design
      for (j in 1:p) {
        designstar[,j]<-design[,j]+movedesign[partindex[i],j]
      }
      design<-designstar%%1
      rdata<-as.matrix(data[((partindex[i]-1)*N/part+1):(partindex[i]*N/part),-1])
      tryCatch({
        svdrdata<-svd(scale(rdata))
      },error = function(e){
        for (a in 1:dim(rdata)[2]) {
          small<-which(abs(rdata[sample(1:dim(rdata)[1],min(dim(rdata)[1],100)),a])!=0)*0.0001
          rdata[,a]<-rdata[,a]+small**stats::rnorm(dim(rdata)[1],0,1)
        }
      })
      rdata<-svdrdata$u[,1:p]
      yita<-matrix(nrow=n/part,ncol=p)
      for (j in 1:p) {
        yita[,j]<-stats::quantile(rdata[,j],design[,j])
      }
      kdtree<-RANN::nn2(data=rdata,query=yita,k=1,treetype="kd")
      subsample<-c(subsample,data[(kdtree$nn.idx+(partindex[i]-1)*N/part),1])
    }
    subsample
  })
  parallel::stopCluster(cl)
  return(unlist(samplelist))
}
