\name{PDDS}
\alias{PDDS}
\title{Obtain representative subsamples of a dataset by PDDS.}
\usage{
  PDDS(data,n,partperd,sp,move,corenumber)
}
\description{
  This function extracts n representative subsamples from the input dataset by Parallel Data-Driven Subsampling.
}
\arguments{

  \item{data}{a float matrix R object. The original data.}

  \item{n}{an integer R object. The number of subsamples.}

  \item{partperd}{an interger vector R object. The number of blocks into which data is sequentially divided along each dimension in sp. Note that the length(partperd) must be no more than the dimensions of the data, and the prod(partperd) must be a factor of n. We recommend the length(partperd) is no more than 3.
  }

  \item{sp}{an interger vector R object. The dimensions used to sequentially split the data, which is provided by the function seqdimchoice(data,partperd) by default. length(partperd) must equal to length(sp).
  }

  \item{move}{a boolean R object. It indicates whether to perturb the uniform designs used on different data blocks. The default setting is "T".
  }

  \item{corenumber}{an integer R object.
The number of CPU cores used for parallel processing. Please set it to no more than prod(partperd).
  }
}

\value{Indices of the subsamples.}

\examples{
  x1 = rnorm(n = 10000, mean = 0, sd = 1)
  x2 = rnorm(n = 10000, mean = x1^2, sd = 1)
  x<-cbind(x1,x2)
  subsample_index<-PDDS(data = x, n = 100, partperd = c(5,5))
}
