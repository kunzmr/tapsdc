#' Fit the Standard Diffusion Curve
#'
#' Fit the standard diffusion curve coefficients given a flux and time array.
#' This function uses the Levenberg-Marquardt algorithm for a non-linear least squares fit of the coefficients.
#'
#' @param x An array of time values in which the flux was collected.
#' @param y An array of flux measurements.
#' @param n_max An integer for the total number of attempts for the LM-optimization.
#' @param fixed_values A list of fixed parameters. For example, assuming the flux is irreversible, then desorption must be zero, i.e., fixed_values = list(desorption = 0).
#' @return A list of results including the time series data as a dataframe, the optimized coefficients as a dataframe, the non-linear least squares fit of the coefficients, and the least squares fit of the sdc to the original data.
#' @importFrom pracma trapz
#' @importFrom minpack.lm nlsLM
#' @examples
#' x = (1:5000) * 0.001
#' flux = sdc(x) + rnorm(length(x), mean = .2, sd = .01)
#' fit = fit_single_sdc(y, x = x, fixed_values = list(desorption = 0))
#'
#' plot(fit$data$time_array, fit$data$flux)
#' lines(fit$data$time_array, fit$data$fitted, col = "red")
#'
#' @export fit_single_sdc
fit_single_sdc = function(x, y, n_max = 1, fixed_values = NULL) {
  result = list(
    data = data.frame(
      x = x,
      y = y,
      fitted = rep(0, length(y)),
      sdc = rep(0, length(y))
    ),
    coefs = data.frame(
      baseline = 0,
      beta = 0,
      tau = 0,
      ka = 0,
      kd = 0
    ),
    nls_model = NULL,
    ls_model = NULL
  )

  # ensure the values of the flux non-negative
  y = y - min(y)
  # area normalize the flux
  y = y / pracma::trapz(x, y)
  time_step = median(diff(x))

  # Estimate baseline and area
  #fit_ar = lm(y[-1] ~ y[-length(y)])

  #noise_term = rnorm(length(y), 0, sd(diff(y)))
  #fit_noise = lm(y ~ noise_term)

  init_area = pracma::trapz(x, y - median(y))
  init = list(
    baseline = median(y), # fit_noise$coefficients[1], #as.numeric((median(y) + coef(fit_ar)[1]) / 2),
    beta = init_area, # as.numeric(max(c(time_step, coef(fit_ar)[2]))),
    tau = 1 / (6 * x[which.max(y)]), # pracma::trapz(x, x * y),
    ka = abs(pracma::trapz(x, x * (y - median(y)))),
    kd = 0
  )

  # tau = 2 * sd(y)

  str_form = c("y ~ ", "baseline", " + ", "beta",
               " * tapsdc::sdc(x, ", "tau", ", ",
               "ka", ", ", "kd", ")")


  for (i in names(fixed_values)) {
    init_idx = which(names(init) == i)
    form_idx = which(str_form == i)
    if (length(form_idx) != 0) {
      init = init[-init_idx]
      str_form[form_idx] = paste0(str_form[form_idx], " = ", fixed_values[[i]])
      result$coefs[[i]] = fixed_values[[i]]
    }
  }
  str_form = paste0(str_form, collapse = "")

  fit_list = vector(mode = "list", length = n_max)
  fit_rmse = rep(1e6, n_max)

  for (i in 1:n_max) {
    fit_list[[i]] = try(minpack.lm::nlsLM(formula(str_form), start = init),
                        silent = TRUE
    )
    if (class(fit_list[[i]]) == "nls") {
      fit_rmse[i] = sqrt(mean(summary(fit_list[[i]])$residuals^2))
    }
  }

  nls_model = fit_list[[which.min(fit_rmse)[1]]]

  if (class(nls_model) == "nls") {
    for (i in names(coef(nls_model))) {
      if (i %in% c("tau", "ka", "kd")) {
        result$coefs[[i]] = abs(coef(nls_model)[[i]])
      }
    }
    result$data$sdc = sdc(result$data$x, result$coefs$tau,
                          result$coefs$ka, result$coefs$kd)
    ls_model = lm(result$data$y ~ result$data$sdc)
    result$coefs$baseline = coef(ls_model)[1]
    result$coefs$beta = coef(ls_model)[2]
    result$data$fitted = ls_model$fitted.values
    result$nls_model = nls_model
    result$ls_model = ls_model
  }

  return(result)
}

