\name{makeValidParams}
\alias{makeValidParams}

\title{Non-standard Generic for Checking Validity of Parameter Objects}
\description{
  This function is not intended for end-users, but may be useful for
  developers extending the Hypergeometric testing capabilities
  provideded by the Category package.

  \code{makeValidParams} is intended to validate a parameter object
  instance (e.g. HyperGParams or subclass).  The idea is that unlike
  \code{validObject}, methods for this generic attempt to fix invalid
  instances when possible, and in this case issuing a warning, and
  only give an error if the object cannot be fixed.
}
\usage{
makeValidParams(object)
}
\arguments{
  \item{object}{A parameter object.  Consult \code{showMethods} to
  see signatures currently supported.}
}
\value{
  The value must have the same class as the \code{object} argument.
}

\author{Seth Falcon}

\keyword{manip}

