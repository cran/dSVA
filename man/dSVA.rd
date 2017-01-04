 \name{dSVA}
 \alias{dSVA}
 \title{direct surrogate variable analysis}
 \description{
     Identify hidden factors in high dimensional biomedical data   
 }
 \usage{

dSVA(Y, X, ncomp=0)


 }
\arguments{
      \item{Y}{n x m data matrix of n samples and m features.}
      \item{X}{n x p matrix of covariates without intercept. }
      \item{ncomp}{ a number of surrogate variables to be estimated. 
      If ncomp=0 (default), ncomp will be estimated using the be method in the num.sv function of the sva package.  }
      
}
\value{
Bhat = Bhat.all[idx.test,], BhatSE= BhatSE[idx.test,], Pvalue=Pvalue
	\item{Bhat}{n x m matrix of the estimated effect sizes of X }
	\item{BhatSE}{n x m matrix of the estimated standard error of Bhat }	
	\item{Pvalue}{n x m matrix of the p-values of Bhat }
	\item{Z}{a matrix of the estimated surrogate variable}
  	\item{ncomp}{a number of surrogate variables.}
	
}

\author{Seunggeun Lee}


\examples{


data(Example)
attach(Example)
out<-dSVA(Y,X, ncomp=0)

}


