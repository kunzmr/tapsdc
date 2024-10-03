#' Fit the Standard Diffusion Curve
#'
#' Fit the standard diffusion curve coefficients given a y (or a series of y) and time array.
#' This function uses the Levenberg-Marquardt algorithm for a non-linear least squares fit of the coefficients.
#'
#' @param x An array of time values in which the y was collected.
#' @param y A matrix or array of flux measurements.
#' @param n_max An integer for the total number of attempts for the LM-optimization.
#' @param fixed_values A list of fixed parameters. For example, assuming the y is irreversible, then desorption must be zero, i.e., fixed_values = list(kd = 0).
#' @param n_cores An integer for the total number cores to use in calculation.
#' @return The standard diffusion curve based on the kinetic coefficients given.
#' @importFrom parallel detectCores makeCluster parLapply stopCluster
#' @examples
#' x = (1:5000) * 0.001
#'
#' y = sdc(x) + rnorm(length(x), mean = .2, sd = .01)
#' fit = fit_single_sdc(x = x, y = y, fixed_values = list(desorption = 0))
#'
#' y2 = sdc(x) + rnorm(length(x), mean = .1, sd = .05)
#' y = as.data.frame(cbind(y, y2))
#' names(y) = c("f1", "f2")
#'
#' fit = fit_sdc(x, y, fixed_values = list(desorption = 0))
#' names(fit)
#' plot(fit$data$x, fit$data$y)
#' lines(fit$data$x, fit$data$fitted, col = "red")
#'
#' @export fit_sdc
fit_sdc = function(x, y, n_max = 1, fixed_values = NULL, n_cores = NULL) {


  if (is.null(dim(y))) {
    result = tapsdc::fit_single_sdc(x = x, y = y, n_max = n_max)
  } else {

    y_list = vector(mode = "list", length = ncol(y))
    for (i in 1:ncol(y)) {
      y_list[[i]] = y[, i]
    }

    if (is.null(n_cores)) {
      n_cores = parallel::detectCores() - 1
    }
    cl = parallel::makeCluster(n_cores)
    results = parallel::parLapply(cl, y_list, tapsdc::fit_single_sdc,
                                  x = x,
                                  fixed_values = fixed_values, n_max = n_max)
    parallel::stopCluster(cl)

    result = list(
      data = NULL,
      coefs = NULL,
      nls_model = vector(mode = "list", length = ncol(y)),
      ls_model = vector(mode = "list", length = ncol(y))
    )
    for (i in 1:length(results)) {
      results[[i]]$data$index = colnames(y)[i]
      results[[i]]$coefs$index = colnames(y)[i]

      if (i == 1) {
        result$data = results[[i]]$data
        result$coefs = results[[i]]$coefs
      } else {
        result$data = rbind(result$data, results[[i]]$data)
        result$coefs = rbind(result$coefs, results[[i]]$coefs)
      }
      result$nls_model[[i]] = results[[i]]$nls_model
      result$ls_model[[i]] = results[[i]]$ls_model
    }
  }

  return(result)
}
