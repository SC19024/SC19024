\name{SC19024-package}
\alias{SC19024-package}
\alias{SC19024}
\docType{package}
\title{
SC19024
}
\description{
This package implemeted several useful functions in statistical computation. The first function RMonte help you to do monte carlo integration with reduced variance. The second function improve the speed of bootstrap. The last function can sample from an arbitrary exp-family distribution.
}
\details{
RMonte: mote carlo integration with reduced variance.

Cboot: increasing speed of bootstrap

C_poly generate sample from exp-family distribution by MCMC

}
\author{
SC19024

Maintainer: SC19024
}
\references{
statiscal computation
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
    RMonte(FUN=exp,low = 0,up = 5,N=1000,layer = 5, Method = 'control')
    Cboot(DATA = 1:100,sims = 10,Method = 'Perc',alpha = 0.05)
    C_poly(c(1),c(0,0,-0.5),1,10000,-1000,1000)
  }
}
