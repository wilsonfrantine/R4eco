### From raw to percentage ####

#' @name relative
#' @title Relative abundance calculation
#' @export
#' @description this function returns the relative proportion of each element based on a vector or dataframe.
#' @export
#' @author Wilson Frantine-Silva <wilsonfratine@@gmail.com>
#' @param df a dataframe or a vector of numeric values.
#' @details This function was designed to perform faster than df/colSums(df) to be used in a bootstrap contest into Renkonen function.
#' @usage relative(df)
#' @examples 
#' df <- simulate_data(5,10)
#' head(df)
#' df.rel <- relative(df)
#' head(df.rel)

relative <- function(df){
  if(is.null(N<-ncol(df))){
    N <- 1; results <- df/sum(df); return(results)
  }else{
    df.temp<-c()
    for(i in 1:(N)){
      df.temp<-cbind(df.temp, df[,i]/sum(df[,i]))
    }
    colnames(df.temp)<-colnames(df)
    results<-data.frame(df.temp)
    return(results)
  }
}


##Utilities ####
#' @nans.remove
#' @title Remove nans from an object
#' @description this function removes from a given object
#' @param x a vector, list or data.frame
nans.remove<-function(x){
  x[apply(x, 2, is.nan)]<-0
  return(x)
}

#Simulating data.sets###

#' @name simulate_data
#' @title Simulate a community data
#' @export
#' @description this function simulates a single community data using a pseudorandom root, uniform and poison distributions
#' @param n.sites the number of sites required or n of columns required
#' @param n.species the n of rows representing species
#' @param min.ab the expected minimum abundance
#' @param mean.ab mean abundance for a given column
#' @param seed an integer to be used as random seed.
#' @usage simulate_data(n.sites, n.species, min.ab, mean.ab, seed)
#' @examples 
#' simulate_data()
#' simulate_data(n.species=10,n.sites=20)
#' simulate_data(n.species=10,n.sites=20, 
#'     min.ab=10, mean.ab=100, seed=123)

simulate_data <- function (n.sites=5, n.species=10, min.ab=0, 
                           mean.ab=100, seed=1234) {
  sites.names = sprintf("P%s",1:n.sites)
  sp.names    = sprintf("sp%s",1:n.species)
  df          = data.frame(1:n.species)
  
  for(i in 1:n.sites){
    set.seed(i*seed)
    d = round(runif(n.species, min.ab, mean.ab) *0.25* rpois(n.species,mean.ab/runif(n.species, mean.ab*0.1, mean.ab*0.2)))
    df <- cbind(df, d)
  }
  
  df            <- df[,-1]
  names(df)     <- sites.names
  row.names(df) <- sp.names
  
  return(df)
}

############ Handling trees ##########

#' @name plot_upgma
#' @title Plot UPGMA phylo objects
#' @description Users can directly plot a single tree or the consensus of a multiPhylo object. The consensus is controled by the function tree_consensus. 
#' @export
#' @import ape
#' @seealso tree_consensus() upgma()
#' @author Wilson Frantine-Silva <wilsonfrantine@@gmail.com>
#' @param tree a single tree or a set of trees in a multiPhylo object
#' @param ... parameters of the tree_consensus function
#' @usage plot_upgma(tree, ...)
#' @examples 
#' df <- simulate_data()
#' S <- renkonen(df, 100)
#' trees <- upgma(S)
#' plot_upgma(trees, "mcc")

plot_upgma <- function(tree, ...){
  if(grepl("phylo",class(tree)[[1]])){
    plot(tree)
  }else if(grepl("multiPhylo",class(tree)[[1]])){
    tr <- tree_consensus(tree, ...)
    plot(tr)
    ape::nodelabels(round(tr$node.label*100), frame = "none", cex=0.8)
  }else if(base::is.list(tree)){
    stop("This is not a phylo object, are you trying to run a list of distance matrix instead? You might wanna run the upgma function")
  }else{
    stop(paste0("Your object is of the class: ",class(tree), " \n ",
                "You need to input a phylo or multiPhylo object"))
  }
}

#' @name tree_consensus
#' @title Compute tree consensus for multiple trees
#' @export
#' @description this function returns a consensus tree of the class phylo. I can be used to be saved or plot with diferent functions and pakcages
#' @param trees an object of the class multiPhylo
#' @param type a type of consensus. Default is "mcc" (Maximum Credibility Tree), but majority rule ("majority") or "strict" rules are also avaliable. In case of majority rule, the user also might use the param "p".
#' @param p the proportion of concordance for some specific node. Default is 0.5, but users might change that between 0 and 1. If 1 is seted the algoritms behaves as strict rule.
#' @details this function is an wrap of ape and phangorn packages
#' @usage tree_consensus(trees, type, p)
#' @examples 
#' df <- simulate_data()
#' S <- renkonen(df, 100)
#' trees <- upgma(S)
#' mcc_tree <- tree_consensus(trees, "mcc")

tree_consensus <- function(trees, type="mcc", p=0.5){
  if(grepl("multiPhylo", class(trees)[[1]])){
    if(type=="mcc"){
      tree <- phangorn::mcc(trees)
      bootvalues <- bootvalues<- ape::prop.clades(tree, trees)/length(trees)
      tree$node.label <- bootvalues
    }else if(type=="strict"){
      tree <- ape::consensus(trees)
      bootvalues <- bootvalues<- ape::prop.clades(tree, trees)/length(trees)
      tree$node.label <- bootvalues
    }else if(type=="majority"){
      tree <- ape::consensus(trees, p = p)
      bootvalues <- bootvalues<- ape::prop.clades(tree, trees)/length(trees)
      tree$node.label <- bootvalues
    }
  }else{
    stop("The object is not a multiPhylo")
  }
  return(tree)
}


