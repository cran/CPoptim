%% File CPoptim.Rd
%% Part of the CPoptim R package
%% Copyright 2013 Erick G.G. de Paz, Humberto Vaquera Huerta, Francisco J Albores Velasco, John R. Bauer Mengelberg, Juan Manuel Romero Padilla
%% Distributed under GPL 2
\name{CPoptim}
\alias{CPoptim}
\title{CPoptim}
\description{
  Minimise a given objective function in a bounded search space.
}
\usage{
CPoptim(rFunction,lower,upper,maxFE,sampleSize)
}
\arguments{
  \item{rFunction}{the function to minimise}
  \item{lower}{a vector providing the lower bounds of the search space}
  \item{upper}{a vector providing the upper bounds of the search space}
  \item{maxFE}{number of function evaluations to compute (default 5000*dim)}
  \item{sampleSize}{sample size per partition (default 1000)}
}
\details{
Convex Partition is a black-box optimisation algorithm for single objective functions with real parameters. The basic principle is to progressively estimate and exploit a regression tree similar to CART (Classification and Regression Tree) of the objective function. This model is computed by recursively partitioning the function domain into subsets and estimating local statistics for each one. The subsets most likely to contain the optimum are selected by the Bayesian inference method proposed in \emph{de Paz et al. (2024)}.

The optimisation strategy consists of iterating two phases: 1) Selecting the most promising subset to contain extreme values according to the current model, and 2) Updating the model by sampling over the most promising subsets.
}
\value{
  CPoptim returns two lists: \code{sample} and \code{subsets}. If the function call is not proper, CPoptim returns \code{NULL}.

  \code{sample} contains three matrices that summarise information about each evaluated point. The i-th row of these matrices contains:
  \item{sample$x}{the i-th point that was evaluated}
  \item{sample$y}{the function evaluation for the i-th point}
  \item{sample$subsetID}{the subset-id to which the i-th point belongs}

  \code{subsets} summarises information about each defined subset. The i-th row of each matrix contains
  \item{subsets$lower}{the lower bounds of the i-th subset}
  \item{subsets$upper}{the upper bounds of the i-th subset}
  \item{subsets$mean}{the mean value of the objective fun. in the i-th subset}
  \item{subsets$std}{the standard deviation of the objective fun. in the i-th subset}
  \item{subsets$aPriori}{the relative (length, area, volume, etc.) of the i-th subset}
  \item{subsets$aPosteriori}{the posteriori probability that the i-th subset contains the optimum}
}



\author{
  The design is inspired by the algorithm proposed in \emph{de Paz et al. (2024)}.  However, instead of the original regression tree based on simplexes, this implementation is based on hyper-rectangular subsets (a model similar to the continuous Classification and Regression Trees) \emph{Loh (2011)}.
}


\references{

  de Paz, Erick G.G., et al. (2024).
  \emph{A Regression Tree as Acquisition Function for Low-dimensional Optimisation.} Pattern Recognition. \emph{MCPR 2024. Lecture Notes in Computer Science, vol} \bold{14755}. Springer, Cham. \doi{10.1007/978-3-031-62836-8_3}

  Loh, Wei-Yin (2011).
  \emph{Classification and regression trees.} WIREs Data Mining Knowl Discov 2011 1 14-23. John Wiley & Sons, Inc. \doi{10.1002/widm.8}
}


\examples{
## An illustrative function
sphere<-function(X) sum(X^2)^.5
bounds<-rep(5,10)

## Calling the CPoptim function
obj<-CPoptim(sphere, -bounds, +bounds)

## Ploting the convergence curve
plot(obj$sample$y,t='l')

## Selecting the best X evaluated
optimum.x<-obj$sample$x[which.min(obj$sample$y),]
}
