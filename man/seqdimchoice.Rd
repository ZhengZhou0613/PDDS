\name{seqdimchoice}
\alias{seqdimchoice}
\title{Choosing dimensions for sequentially splitting the data in PDDS.}
\usage{
  seqdimchoice(data,partperd,samplenumber,corenumber)
}
\description{
  This function uses a small portion of a large-scale dataset to determine which dimensions should be used to sequentially split the data when performing subsampling by PDDS.
}
\arguments{

  \item{data}{a float matrix R object. The original data.}

  \item{partperd}{an interger vector R object. The number of blocks into which data is sequentially divided along each dimension. Note that the length(partperd) must be no more than the dimensions of the data.}

  \item{samplenumber}{an interger R object. The number of samples used by this function. The default setting is prod(partperd)*10*dim(data)[2].
  }

  \item{corenumber}{an integer R object.
The number of CPU cores used for parallel processing. Please set it to no more than the dimensions of the original data. We recommend the length(partperd) is no more than 3.
  }
}

\value{The dimensions used to sequentially split the data in PDDS.}

\examples{
  x1 = rnorm(n = 10000, mean = 0, sd = 1)
  x2 = rnorm(n = 10000, mean = x1^2, sd = 1)
  x<-cbind(x1,x2)
  sp<-seqdimchoice(data = x,partperd = c(5,5))
}
