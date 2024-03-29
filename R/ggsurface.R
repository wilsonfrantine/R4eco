#' @title Plot Model Surface with ggplot2
#' @description a ggplot2 2D surface plot function for R glm or lm models
#' @author Wilson Frantine-Silva
#' @name ggsurface
#' @export
#' @param m an object glm, lm or loess object model from stats package
#' @param x.var the name or columns position of the variable to be plot at the x axis. This param has to match with the name or column position of the variable as input in the model object
#' @param y.var same as x.var but for the y axis
#' @param x.label a string with label for the x axis (optional). If NULL the function return x.var
#' @param y.label same as x.label for the y axis
#' @param legend.title a string with name for the response variable. If null the function returns 'response variable'.
#' @param low.col a hexadecimal or ggplot standard color name to represent the low fitting. Default is "grey80"
#' @param high.col a hexadecimal or ggplot standard color name to represent the high fitting. Default is "blue"
#' @param xgrid.size  An integer representing the size of the grid to create a prediction. This variable controls the number of elements to build a new x.var interval to fit the model. The default is 15.
#' @param ygrid.size  same as xgrid.size for the y variable. The default is also 15.
#' @param n.bins An integer. The number of division / breaks for the plot. More bins will result in a plot with more divisions.
#' @param round.legend the number of decimal cases to round the scales at the legend. The default is 0, which rounds to the next integer. Users may also uses -1 to round for the nearest 10th number.
#' @param scale.type the type of prediction required. The default is "link", and returns the scale of the linear predictors. The alternative "response" return the prediction at the reponse variable scale. So if a glm is ran with a binomial link function than the prediction are probabilities at logit scale, but if type = "response" the prediction will be returned at probabiliti scales. See more in stats::predict help.
#' @usage ggsurface(m, x.var, y.var, x.label, y.label, legend.title,
#' low.col, high.col, xgrid.size, ygrid.size, n.bins, round.legend, scale.type)
#' @examples
#' data(mtcars)
#' m <- glm(mpg ~ wt + hp, data=mtcars, family = "gaussian")
#' ggsurface(m, x.var = "wt", y.var = "hp",
#'    legend.title = "milles per galon", high.col = "darkred",
#'    round.legend = 0)
#' @seealso stats::predict, stats::glm, stats::lm, stats::lm
#' @details This function takes a model object and uses the function prediction to fit the model to a grid of x size. Therefore at this version the function needs one of the model objects from the stats package.
#' In case of any crash, plese contact-me at: <wilsonfratntine@@gmail.com>.

ggsurface <- function(m=NULL, x.var=NULL, y.var=NULL, x.label=NULL, y.label=NULL, legend.title=NULL, low.col="grey80", high.col="darkred", xgrid.size=10, ygrid.size=10, n.bins=11, round.legend=0, scale.type="link"){
  
  if(is.null(m)){
    stop("you have to provide a model object")
  }else if("loess" %in% class(m)){
    d <- base::data.frame(m$x)
  } else if("glm" %in% class(m)){
    if(!is.null(m$model)){
      if(is.data.frame(m$model)) d <- base::data.frame(m$model)
    }
  } else {
    stop(paste("The object of the class", class(m)[[1]], "is not supported yet. Please, report the bug at github.com/wilsonfrantine/ggmodel"))
  }
  
  if(is.null(x.var)){
    stop("Please, provide the name or position of the variable you want to use as predictor on the x axis")
  }
  if(is.null(y.var)){
    stop("Please, provide the name or position of the variable you want to use as predictor on the y axis")
  }
  if(is.null(x.label)){
    x.label <- names(d[x.var])
  }
  if(is.null(y.label)){
    y.label <- names(d[y.var])
  }
  if(is.null(legend.title)){
    legend.title <- "prediction at response scale"
  }
  
  x = d[x.var]
  y = d[y.var]
  
  xgrid <- base::seq(base::min(x), base::max(x), length = xgrid.size)
  ygrid <- base::seq(base::min(y), base::max(y), length = ygrid.size)
  
  nd <- base::expand.grid(xgrid, ygrid)
  base::names(nd) <- c(base::names(d[x.var]), base::names(d[y.var]))
  
  z <- stats::predict(m, newdata = nd, type = scale.type)
  
  df <- base::data.frame(nd,"z"=z)
  
  base::names(df) <- c("x", "y", "z")
  x <- df$x
  y <- df$y
  z <- df$z
  
  legend.labels <- base::round(base::seq(
    base::min(z), base::max(z), length=n.bins), round.legend)
  legend.labels[seq(ifelse( (n.bins %% 2) == 0, 1, 0), n.bins-1, 2)] <- ""
  
  p <- ggplot2::ggplot(data = df, mapping = ggplot2::aes(x=x,y=y,z=z))+
    ggplot2::geom_contour_filled(bins = n.bins)+
    ggplot2::labs(x=x.label, y=y.label)+
    ggplot2::guides(fill= ggplot2::guide_legend(
      title = legend.title,
      title.position = "left",
      title.theme = ggplot2::element_text(face = "italic",
                                          angle = 90,
                                          hjust = 0.5,
                                          vjust = 0.5),
      label = T,
      keyheight = ggplot2::unit(0.05,"npc"),
      keywidth  = ggplot2::unit(0.06,"npc"),
      reverse = T,
      title.hjust = 0.5,
      title.vjust = 0.5
    ))+
    ggplot2::scale_fill_manual(values=grDevices::colorRampPalette(c(low.col, high.col))(n.bins), labels=legend.labels)+
    ggplot2::scale_y_continuous(expand=c(0,0))+
    ggplot2::scale_x_continuous(expand = c(0,0))+
    ggplot2::theme(
      axis.text   = ggplot2::element_text(colour = "black" ),
      plot.margin = ggplot2::unit(x = c(5,5,5,5), "mm")
    )
  return(p)
}
