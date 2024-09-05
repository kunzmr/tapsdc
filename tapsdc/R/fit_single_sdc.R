#' Fit the Standard Diffusion Curve
#'
#' Fit the standard diffusion curve coefficients given a flux and time array.
#' This function uses the Levenberg-Marquardt algorithm for a non-linear least squares fit of the coefficients.
#'
#' @param time_array An array of time values in which the flux was collected.
#' @param flux An array of flux measurements.
#' @param n_max An integer for the total number of attempts for the LM-optimization.
#' @param fixed_values A list of fixed parameters. For example, assuming the flux is irreversible, then desorption must be zero, i.e., fixed_values = list(desorption = 0).
#' @return A list of results including the time series data as a dataframe, the optimized coefficients as a dataframe, the non-linear least squares fit of the coefficients, and the least squares fit of the sdc to the original data.
#' @importFrom pracma trapz
#' @importFrom minpack.lm nlsLM
#' @examples
#' time_array = (1:5000) * 0.001
#' flux = sdc(time_array) + rnorm(length(time_array), mean = .2, sd = .01)
#' fit = fit_single_sdc(flux, time_array = time_array, fixed_values = list(desorption = 0))
#'
#' plot(fit$data$time_array, fit$data$flux)
#' lines(fit$data$time_array, fit$data$fitted, col = "red")
#'
#' @export fit_single_sdc
fit_single_sdc = function(time_array, flux, n_max = 1, fixed_values = NULL) {
  result = list(
    data = data.frame(
      time_array = time_array,
      flux = flux,
      fitted = rep(0, length(flux)),
      sdc = rep(0, length(flux))
    ),
    coefs = data.frame(
      intercept = 0,
      beta = 0,
      diffusion = 0,
      adsorption = 0,
      desorption = 0
    ),
    nls_model = NULL,
    ls_model = NULL
  )

  # ensure the values of the flux non-negative
  flux = flux - min(flux)
  # area normalize the flux
  flux = flux / pracma::trapz(time_array, flux)
  time_step = median(diff(time_array))

  # Estimate baseline and area
  fit_ar = lm(flux[-1] ~ flux[-length(flux)])

  coefs = list(
    intercept = (median(flux) + coef(fit_ar)[1]) / 2,
    beta = max(c(time_step, coef(fit_ar)[2])),
    diffusion = 1, adsorption = 1, desorption = 0
  )

  if (!is.null(fixed_values)) {
    for (i in names(fixed_values)) {
      idx = which(names(coefs) == i)
      if (length(idx) != 0) {
        coefs = coefs[-idx]
        assign(i, fixed_values[[i]])
        result$coefs[[i]] = fixed_values[[i]]
      }
    }
  }

  fit_list = vector(mode = "list", length = n_max)
  fit_rmse = rep(1e6, n_max)
  for (i in 1:n_max) {
    fit_list[[i]] = try(
      minpack.lm::nlsLM(
        flux ~ intercept + beta * tapsdc::sdc(time_array, diffusion,
                                           adsorption, desorption),
        start = coefs
      ),
      silent = TRUE
    )
    if (class(fit_list[[i]]) == "nls") {
      fit_rmse[i] = sqrt(mean(summary(fit_list[[i]])$residuals^2))
    }
  }

  nls_model = fit_list[[which.min(fit_rmse)[1]]]

  if (class(nls_model) == "nls") {
    for (i in names(coef(nls_model))) {
      if (i %in% c("diffusion", "adsorption", "desorption")) {
        result$coefs[[i]] = abs(coef(nls_model)[[i]])
      }
    }
    result$data$sdc = sdc(result$data$time_array, result$coefs$diffusion,
                          result$coefs$adsorption, result$coefs$desorption)
    ls_model = lm(result$data$flux ~ result$data$sdc)
    result$coefs$intercept = coef(ls_model)[1]
    result$coefs$beta = coef(ls_model)[2]
    result$data$fitted = ls_model$fitted.values
    result$nls_model = nls_model
    result$ls_model = ls_model
  }

  return(result)
}

