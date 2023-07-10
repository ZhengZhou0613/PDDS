\name{GLP}
\alias{GLP}
\title{Generating A Uniform Design by The Method of Good Lattice Points}
\usage{
  GLP(n,p,type)
}
\description{
  This function takes n,p and the criterion to output a uniform design of order n*p.
}
\arguments{

  \item{n}{an integer R object. The number of rows of the design.}

  \item{p}{an integer R object. The number of columns of the design.}

  \item{type}{a character R object. Type of criterion to use:

    "CD2" -- Centered L2-Discrepancy (default) ;

    "WD2" -- Wrap-around L2-Discrepancy ;

    "MD2" -- Mixture L2-Discrepancy ;
  }
}

\value{A uniform design.}

\examples{
  ud<-GLP(n=10,p=2,"CD2")
}
