\name{DDS}
\alias{DDS}
\title{Obtain representative subsamples of a dataset by DDS.}
\usage{
  DDS(data,n,type,ratio)
}
\description{
  This function extracts n representative subsamples from the input dataset by Data-Driven Subsampling.
}
\arguments{

  \item{data}{a float matrix R object. The original data.}

  \item{n}{an integer R object. The number of subsamples.}

  \item{type}{a character R object. The type of criterion used to generate the uniform design:

    "CD2" -- Centered L2-Discrepancy (default) ;

    "WD2" -- Wrap-around L2-Discrepancy ;

    "MD2" -- Mixture L2-Discrepancy ;
  }
  \item{ratio}{a float R object. The proportion of the variance that can be explained by the components retained after PCA. We suggest taking it no less than 0.85 (default).
  }
}

\value{Indices of the subsamples.}

\examples{
  x1 <- rnorm(n = 10000, mean = 0, sd = 1)
  x2 <- rnorm(n = 10000, mean = 0, sd = 1)
  x<-cbind(x1,x2)
  subsample_index<-DDS(data = x, n = 100, type = "CD2", ratio = 0.85)
}
