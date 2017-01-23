

# Define a custom precis plot function that allows you to alter x axis limits

custom_precis_plot <- 
  function (x, y, pars, cex, col.ci = "black", xlab = "Value", xlim...) {
    
    x <- x@output
    if (!missing(pars)) {
      x <- x[pars, ]
    }
    n <- nrow(x)
    mu <- x[n:1, 1]
    left <- x[[3]][n:1]
    right <- x[[4]][n:1]
    set_nice_margins()
    dotchart(mu, labels = rownames(x)[n:1], cex = cex, 
             xlab = xlab, xlim = xlim...)
    for (i in 1:length(mu)) lines(c(left[i], right[i]), c(i, i), 
                                  lwd = 2, col = col.ci)
    abline(v = 0, lty = 1, col = col.alpha("black", 0.15))
}


# Define a function that allows you to check for divergent iterations in 
# a fit Stan model object

get_divergences <- function(fit.stan.model) {
 
  # Get diagnostic parameters from a fit Stan model
  sampler.params <- get_sampler_params(fit.stan.model, inc_warmup = F)
  
  # Use "n_divergent__" as the default string to match
  match.string <- "n_divergent__"
  
  # If the model has instead a "divergent__" column, convert the match
  # string to this string
  if (sum(colnames(sampler.params[[1]]) %in% "divergent__") > 0)
    match.string <- "divergent__"
  
  # Report the column sums of divergent iteration columns for every chain
  # used in the fit model
  colSums(sapply(sampler.params, function(x) x[, match.string])) 
}
