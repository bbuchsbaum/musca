% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/regal_align.R
\name{extract_node_identity}
\alias{extract_node_identity}
\title{Node Identity Extraction}
\usage{
extract_node_identity(A, k, attributes = NULL)
}
\arguments{
\item{A}{An adjacency matrix representing the graph.}

\item{k}{The maximum number of hops to consider.}

\item{attributes}{A matrix of node attributes (optional).}
}
\value{
A matrix of node identity features, where each row corresponds to a node in the graph.
}
\description{
This function extracts node identity features based on the degree distributions
of k-hop neighborhoods and optional node attributes.
}
