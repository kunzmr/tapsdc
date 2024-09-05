#' Standard Diffusion Curve
#'
#' Generate the standard diffusion curve from Gleaves et al TAP-2.
#'
#' This function describes a theoretical 1-zone TAP flux including diffusion, adsorption, and desorption.
#'
#' @param time_array An array of time values in which the flux was collected.
#' @param diffusion A scalar value of the diffusion coefficient.
#' @param adsorption A scalar value of the adsorption coefficient.
#' @param desorption A scalar value of the desorption coefficient.
#' @param n_max An integer for total number of indices in the iterative calculation of the SDC. Default is 2000 iterations.
#' @return The standard diffusion curve based on the kinetic coefficients given.
#' @examples
#' # Recreation Figure 9 of Gleaves TAP-2
#'
#' time_array = (1:5000) * 0.001
#' diffusion = 1
#' adsorption = 0
#' desorption = 0
#'
#' flux1 = sdc(time_array, diffusion = diffusion,
#'   adsorption = adsorption, desorption = desorption
#' )
#'
#' adsorption = 20
#' desorption = 20
#' flux2 = sdc(time_array, diffusion = diffusion,
#'   adsorption = adsorption, desorption = desorption
#' )
#'
#' desorption = 5
#' flux3 = sdc(time_array, diffusion = diffusion,
#'   adsorption = adsorption, desorption = desorption
#' )
#'
#' plot(time_array, flux1, type = "l", ylim = c(-0.1, 2.4),
#' xlab="Time", ylab = "Flux", main = "Black:ka =0, kd =0; Blue: ka = 20, kd = 20; Red: ka = 20, kd=5")
#' lines(time_array, flux2, col = "blue")
#' lines(time_array, flux3, col = "red")
#' @export sdc

sdc = function(time_array, diffusion = 1, adsorption = 0, desorption = 0, n_max = 2000) {
  # ensuring all values are positive
  diffusion = abs(diffusion)
  adsorption = abs(adsorption)
  desorption = abs(desorption)
  # time cannot be zero
  time_array = time_array[time_array != 0]

  # Most quantities will be transformed to a matrix for computation speed.
  time_array = matrix(time_array, nrow = 1)
  # equivalent to 2n + 1
  pn = ((0:n_max) + 0.5) * pi
  sign_n = 2 * diffusion * c(1, rep(c(-1, 1), ceiling(n_max / 2)))[1:length(pn)] * pn
  kn = pn^2 + adsorption + desorption
  r_sqrt = sqrt(kn^2 - 4 * pn^2 * desorption)
  An = ((-kn + r_sqrt) / 2 + kn - desorption) / r_sqrt

  beta_neg = matrix(sign_n * An, nrow = 1)
  beta_pos = matrix(sign_n * (1 - An), nrow = 1)
  exp_neg = exp(diffusion * matrix((-kn - r_sqrt) / 2, ncol = 1) %*% time_array)
  exp_pos = exp(diffusion * matrix((-kn + r_sqrt) / 2, ncol = 1) %*% time_array)

  flux = abs(as.numeric(beta_neg %*% exp_neg + beta_pos %*% exp_pos))
  return(flux)
}
