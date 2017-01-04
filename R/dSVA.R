Get_Inv<-function(Mat){

  if(dim(Mat)[1] == 1){
    return(1/Mat)
  }

  out_eigen<-eigen(Mat)

  IDX1<-which(out_eigen$values > 0.0001)
  value<-rep(0,length(out_eigen$value))
  value[IDX1]<-1/out_eigen$values[IDX1]

  Mat_Inv<-out_eigen$vectors %*%  (t(out_eigen$vectors ) * value)
  return(Mat_Inv)


}

Get_Pvalue<-function(Y, X1, SV1){

	X.SV<-cbind(X1,SV1)
	
	p<-dim(X.SV)[2]
	n<-dim(X.SV)[1]
	m<-nrow(Y)
	df1 = n-p
	idx.test<-2:ncol(X1)
	
	
	X.SV.inv<-Get_Inv(t(X.SV)%*%X.SV)
	

	Bhat.all=X.SV.inv%*%(t(X.SV)%*%Y )
	Ehat=Y-X.SV%*%Bhat.all
	Var_Est<-apply(Ehat^2,2,sum) / df1
	BhatSE<-sqrt(diag(X.SV.inv) %*% t(Var_Est))
	
	T.stat<- Bhat.all[idx.test, ] / BhatSE[idx.test,]
	Pvalue<-pf(T.stat^2, df1=1, df2= df1, lower.tail =FALSE)
	
  	re<-list(Bhat = Bhat.all[idx.test,], BhatSE= BhatSE[idx.test,], Pvalue=Pvalue)
  	return(re)
  
}


#######################################
#
# Function to compute Surrogate Variable
#   1) Input
#   Y : n x m matrix of gene expression. n is the number of samples and m is the number of genes
#   X : n x p matrix of covariates. It shouldn't include intercept.
#   ncomp : number of residual PCs used to compute surrogate variables. If ncomp=0 (default), ncomp will be computed using the num.sv function of the sva package.
#


dSVA<-function(Y,X, ncomp=0){

	# ncomp=0; X<-X[,-1]
	X1<-model.matrix(~1+X)
	n<-dim(X1)[1]
	if(ncomp==0){
		mod<-X1
  		id <- num.sv(t(Y), mod, B = 20, method = "be", seed = NULL)
  		ncomp <- id
  	}


  	q = ncomp
  	Bhat.Org=solve(t(X1)%*%X1)%*% (t(X1)%*%Y)
	Ehat=Y-X1%*%Bhat.Org

	current.svd=svd(Ehat)
	u=current.svd$u
	d=diag(current.svd$d)
	v=current.svd$v
	
	if(q > 1){
		Delta =  d[1:q,1:q]%*%(t(v)[1:q,])
	} else {
		Delta =  rbind(d[1:q,1:q]*(t(v)[1:q,])   )
  	}
	
	# Proposed Method
	Delta1<-t(Delta - apply(Delta,1,mean))
	Bhat.Org1<-Bhat.Org - apply(Bhat.Org,1,mean)
	
	# Estimate 
  	rho_est1<-(solve(t(Delta1) %*% Delta1) %*% t(Delta1) %*%  t(Bhat.Org1) )
  	Bhat.Adj<-Bhat.Org -  t(rho_est1) %*% Delta 

   	MZ1<-rho_est1
   	Z.est = u[,1:q] + X1 %*% t(MZ1)

	# Get p-value
	re<-Get_Pvalue(Y, X1, Z.est)
	re$Z  = Z.est
	re$ncomp=ncomp
   	
   	return(re)
}




