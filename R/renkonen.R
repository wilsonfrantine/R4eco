#### Basic equatition Reikonen Similarity ####
#' @name renkonen.index
#' @title Renkonen equation
#' @import stats utils
#' @export
#' @description the very core of the Renkonen index function. See more in concept and details.
#' @author Wilson Frantine-Silva <wilsonfrantine@@gmail.com>
#' @param a a vector representing coutings of elements in given group
#' @param b same as "a", but from a group to be compared
#' @concept Renkonen index is given by the summation of the paralella minimum of to vectors of the same size. In ecological context, it gives the summation of the minimal relative abundances of a species set from two sites. So, if one has sp1 and sp2, representing respectivelly 0.3 and 0.7 at site1 and 0.7 and 0.3 at site2, the Renkonen's index will be S = 0.3 + 0.3 = 0.6. Then, the similarity between site1 and site2 will be 0.6
#' @details this function calls the relative abundance to avoid users to use absolute countings as input. This function is called by other functions as reikonen() and reikonen.c().

renkonen.index <-function(a,b){
  df<-relative(cbind(a,b))
  return(sum(pmin(df[,1], df[,2])))
}

#' @name renkonen.c
#' @title A caller for calculating the Renkonen's index
#' @keyword internal
#' @description This function returns a list of matrix and lists of pairwise comparisons of several sampling sites and its respective Renkonen's similarities.
#' @details As renkonen.index is only able to handle vectors, this function parse data.frames and calls the renkonen.index function. This is the basic version for beeing called by a bootstrap handler.
#'This function was not designed to be used directly at the first place, but users are free to tweeck this.
#' @param df a dataframe which each column is one site and each row is the absolute abundance of a species
#' @param as.matrix a boolean whether it showd be returned as a matrix as well
#' @usage renkonen.c(df,as.matrix=T)


renkonen.c <- function (df, as.matrix=T ){
  
  if(!base::is.data.frame(df)){
    stop("you did not provide a data.frame as input")
  }
  
  N = ncol(df) #n of sites  
  n = nrow(df) #n of species
  S <- vector() #to hold the Renkonen's index results
  p1<- vector() #sampling site 1 vector
  p2<- vector() #sampling site 2 vector
  
  for( i in 1:(N) ){
    for( j in (i):N ){
      renkonen.S = renkonen.index(df[,i],df[,j])
      S<-rbind(S,renkonen.S)
      p1<-rbind(p1,i)
      p2<-rbind(p2,j)
    }
  }
  df.new<-data.frame(S, p1, p2, row.names = seq(1:length(S)))
  df.new<-nans.remove(df.new)
  return(df.new)
}

#' @name renkonen
#' @title Renkonen's Similarity index with bootstrap
#' @author Wilson Frantine-Silva <wilsonfrantine@@gmail.com>
#' @export
#' @description this function returns both single or bootstrap Renkonen's index for a data frame of counting input. It uses a internal bootstrap function to calculate N matrix. Users might then extract confidence intervals or standard deviation. 
#' @param df a data frame with sampling sites as columns, species as rows and abundances as values
#' @param boot an integer with bootstrap replications desired. The default is 0, but something around 1000 seems to be suitable in most cases.
#' @details This function is the main envelope of Renkonen's index calculations to the users interact with. implements a internal bootstrap replication loop designed to be as fast as possible in a high level language. Remember, bootstrap process are time consuming and the total time will depend on your data set.
#' @usage renkonen(df, boot)
#' @examples 
#' df <- simulate_data(5,10)
#' S <- renkonen(df)
#' head(S)

renkonen <- function(df, boot=0){
  n = nrow(df)
  samples=c(1:n)
  results <- list()
  res <- list()
  res.matrix <- list()
  if(boot != 0 ){
    pb<-txtProgressBar(min = 1, max = boot, style = 3, char = "=")
    randomvector<-sample(1:n,boot, replace = T)
    for(i in 1:boot){
      jacknife<-samples[-randomvector[i]]
      bootvector<-c(jacknife, sample(jacknife,1))
      dfs <- df[bootvector,]
      dfs <- relative(dfs)
      res[[i]]<-renkonen.c(dfs)
      res.matrix[[i]]<-reshape(data = res[[i]], idvar = "p2", timevar = "p1", direction = "wide")[,-1]
      setTxtProgressBar(pb, i)
    }    
  }else{
    dfs<-relative(df)
    res<-renkonen.c(dfs)
    res.matrix<-reshape(data = res, idvar = "p2", timevar = "p1", direction = "wide")[,-1]
    
  }
  results<-list("s"=res, "m"=res.matrix)
  return(results)
}