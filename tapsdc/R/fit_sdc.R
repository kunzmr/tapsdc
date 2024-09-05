#' Fit the Standard Diffusion Curve
#'
#' Fit the standard diffusion curve coefficients given a flux (or a series of flux) and time array.
#' This function uses the Levenberg-Marquardt algorithm for a non-linear least squares fit of the coefficients.
#'
#' @param time_array An array of time values in which the flux was collected.
#' @param flux A matrix or array of flux measurements.
#' @param n_max An integer for the total number of attempts for the LM-optimization.
#' @param fixed_values A list of fixed parameters. For example, assuming the flux is irreversible, then desorption must be zero, i.e., fixed_values = list(desorption = 0).
#' @param n_cores An integer for the total number cores to use in calculation.
#' @return The standard diffusion curve based on the kinetic coefficients given.
#' @importFrom parallel detectCores makeCluster parLapply stopCluster
#' @examples
#' time_array = (1:5000) * 0.001
#'
#' flux = sdc(time_array) + rnorm(length(time_array), mean = .2, sd = .01)
#' fit = fit_single_sdc(flux, time_array = time_array, fixed_values = list(desorption = 0))
#'
#' flux2 = sdc(time_array) + rnorm(length(time_array), mean = .1, sd = .05)
#' flux = as.data.frame(cbind(flux, flux2))
#' names(flux) = c("f1", "f2")
#'
#' fit = fit_sdc(time_array, flux, fixed_values = list(desorption = 0))
#' names(fit)
#' plot(fit$data$flux)
#' lines(fit$data$fitted, col = "red")
#'
#' @export fit_sdc
fit_sdc = function(time_array, flux, n_max = 1, fixed_values = NULL,
                   n_cores = NULL) {


  if (is.null(dim(flux))) {
    result = tapsdc::fit_single_sdc(time_array = time_array, flux = flux, n_max = n_max)
  } else {

    flux_list = vector(mode = "list", length = ncol(flux))
    for (i in 1:ncol(flux)) {
      flux_list[[i]] = flux[, i]
    }

    if (is.null(n_cores)) {
      n_cores = parallel::detectCores() - 1
    }
    cl = parallel::makeCluster(n_cores)
    results = parallel::parLapply(cl, flux_list, tapsdc::fit_single_sdc,
                                  time_array = time_array,
                                  fixed_values = fixed_values, n_max = n_max)
    parallel::stopCluster(cl)

    result = list(
      data = NULL,
      coefs = NULL,
      nls_model = vector(mode = "list", length = ncol(flux)),
      ls_model = vector(mode = "list", length = ncol(flux))
    )
    for (i in 1:length(results)) {
      results[[i]]$data$index = colnames(flux)[i]
      results[[i]]$coefs$index = colnames(flux)[i]

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
