% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basics_utilities.R
\name{plot_upgma}
\alias{plot_upgma}
\title{Plot UPGMA phylo objects}
\usage{
plot_upgma(tree, ...)
}
\arguments{
\item{tree}{a single tree or a set of trees in a multiPhylo object}

\item{...}{parameters of the tree_consensus function}
}
\description{
Users can directly plot a single tree or the consensus of a multiPhylo object. The consensus is controled by the function tree_consensus.
}
\examples{
df <- simulate_data()
S <- renkonen(df, 100)
trees <- upgma(S)
plot_upgma(trees, "mcc")
}
\seealso{
tree_consensus() upgma()
}
\author{
Wilson Frantine-Silva <wilsonfrantine@gmail.com>
}
