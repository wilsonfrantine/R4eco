#' @title lmerPredictionPlot
#' @description
#'   lmerPredictionPlot is a function for creating plots of predicted intervals from lme4 models using ggplot2. 
#'   It provides flexibility in specifying the variables for the x-axis, y-axis, and grouping, 
#'   and allows customization of plot appearance through additional parameters.
#'   
#' @param model A lme4 model object.
#' @param x The name of the x variable (independent variable). Defaults to NULL. If not provided, 
#'   the function attempts to guess the name by examining the formula of the model.
#' @param y The name of the y variable (dependent variable). Defaults to NULL. If not provided, 
#'   the function isolates the first variable from the second formula term.
#' @param grp A categorical variable of the model. Defaults to NULL. If not provided, 
#'   the function isolates the second variable from the second formula term.
#' @param size The size of the data. By default, the size is set to the number of rows in the model frame.
#' @param ... Additional parameters for the function merTools::predictInterval(); 
#'   see the reference documentation for the function.
#'
#' @examples
#'   d <- data.frame(
#'     Type = rep(c("Forest", "Regeneration", "Restoration"), each = 12),
#'     Landscape = rep(paste0("L", 1:12), times = 3),
#'     Mean_NDVI_SD_500 = rnorm(36, mean = 0.2, sd = 0.02),
#'     hill_q0 = sample(5:14, 36, replace = TRUE),
#'     Abundance = sample(5:50, 36, replace = TRUE)
#'   )
#'   modelX <- lme4::lmer(formula = hill_q0 ~ Mean_NDVI_SD_500 * Type+(1|Landscape), data=d)
#'   
#'   lmerPredictionPlot(model = modelX, type = "linear.prediction")
#'   lmerPredictionPlot(model = modelX) +
#'     scico::scale_color_scico_d(palette = "batlow", begin = 0.1, end = 0.7) +
#'     scico::scale_fill_scico_d(palette = "batlow", begin = 0.1, end = 0.7)
#'   lmerPredictionPlot(model = modelX) +
#'     scico::scale_color_scico_d(palette = "batlow", begin = 0.1, end = 0.7) +
#'     scico::scale_fill_scico_d(palette = "batlow", begin = 0.1, end = 0.7) +
#'     ggplot2::facet_wrap(~Type)
#'     
#'   modelY <- lme4::glmer.nb(formula = Abundance ~ Mean_NDVI_SD_500 * Type+(1|Landscape), data=d)
#'   
#'   lmerPredictionPlot(model = modelY, type = "probability")
#'   lmerPredictionPlot(model = modelY) +
#'     scico::scale_color_scico_d(palette = "batlow", begin = 0.1, end = 0.7) +
#'     scico::scale_fill_scico_d(palette = "batlow", begin = 0.1, end = 0.7)
#'   lmerPredictionPlot(model = modelY) +
#'     scico::scale_color_scico_d(palette = "batlow", begin = 0.1, end = 0.7) +
#'     scico::scale_fill_scico_d(palette = "batlow", begin = 0.1, end = 0.7) +
#'     ggplot2::facet_wrap(~Type)
#' @return A ggplot object.
#' @seealso 
#'   \code{\link{merTools::predictInterval}}, \code{\link{ggplot2}}
#'
#' @export
#' 

lmerPredictionPlot <- function(model, x= NULL, y=NULL, grp= NULL, size = NULL, ...) {
  d    <- model.frame(model)
  x    <- ifelse(is.null(x), attr(terms( formula(model) ), "term.labels")[1], x )
  y    <- ifelse(is.null(y), rlang::as_string(formula(model)[[2]]), y)
  grp  <- ifelse(is.null(grp), attr(terms( formula(model) ), "term.labels")[2], grp )
  size <- ifelse(is.null(size), nrow(d), size)
  
  stat <- "mean"
  level <- 0.95
  which <- "fixed"
  include.resid.var = FALSE
  type <- "probability"
  
  nd <- expand.grid(
    x = seq(min(d[[x]]), max(d[[x]]), length.out = 12),
    grp = levels(d[[grp]]) ) 
  
  colnames(nd) <- c(x, grp)
  
  nd <- data.frame(nd,d[, !colnames(d) %in% c(x[], grp[] ) ] )
  
  ls1 <- list(merMod=model, newdata=nd, stat=stat, level=level, which=which, include.resid.var=include.resid.var, type=type)
  ls2 <- list(...)
  params <- c(ls1[!names(ls1) %in% names(ls2)], ls2)
  
  prediction <- do.call(merTools::predictInterval, params)
  
  ggdata <- data.frame(nd, prediction)
  
  p <- ggplot2::ggplot(ggdata, ggplot2::aes(x = !!dplyr::sym(x), y = fit, color = !!dplyr::sym(grp))) +
    ggplot2::geom_point(data=d, ggplot2::aes(y = !!dplyr::sym(y), x = !!dplyr::sym(x), color= !!dplyr::sym(grp)),   size = 4, alpha = 0.7)+
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lwr, ymax = upr, fill = !!dplyr::sym(grp)), color=NA, alpha = 0.2)+
    ggplot2::geom_line(ggplot2::aes(y=fit))+
    ggplot2::labs(y=dplyr::sym(y)) 
  return(p)
}
