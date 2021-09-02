#### UPGMA ####
#' @name upgma.c
#' @title UPGMA calc core
#' @description this function is just a wraper for the hclut function which calculates an upgma dendrogram from a distance matrix. The function also converts a the hclust object in a phylo object.
#' @param df a squar data frame with similarity index as values
#' @param method same methods as in hclust function
#' @param ... ignored
#' @seealso stats::hclust()
#' @keywords internal
upgma.c <- function(df, method="average", ...){
  dist <- stats::as.dist(df)
  hc <- stats::hclust(dist, method = method)
  result <- ape::as.phylo(hc)
  result <- stats::reorder(result, "postorder")
  return(result)
}

#' @name upgma
#' @title UPGMA calculation
#' @export
#' @description it takes a tree or a list of trees and returns a phylo object (from ape package) with x trees as inputed.
#' @concept the function is based on stats:hclust.
#' @param x a distance matrix or a list of distance matrix as outputed by the renkonen function
#' @param method the UPGMA summation method. It uses the same as hclust function
#' @param ... ignored
#' @seealso renkonen(), stats::hclust()

upgma<- function(x, method="average", ...){
  x<-x$m
  if(!is.data.frame(x)){
    trees<-list()
    class(trees)<-"multiPhylo"
    for(i in 1:length(x)){
      result<-upgma.c(1-x[[i]])
      trees[[i]] <- result
    }
  }else{
    trees<-upgma.c(1-x)
  }
  return(trees)
}
