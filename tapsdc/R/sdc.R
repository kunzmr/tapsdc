#' Standard tau Curve
#'
#' Generate the standard tau curve from Gleaves et al TAP-2.
#'
#' This function describes a theoretical 1-zone TAP flux including tau, ka, and kd.
#'
#' @param x An array of time values in which the flux was collected.
#' @param tau A scalar value of the tau coefficient.
#' @param ka A scalar value of the ka coefficient.
#' @param kd A scalar value of the kd coefficient.
#' @param n_max An integer for total number of indices in the iterative calculation of the SDC. Default is 2000 iterations.
#' @return The standard tau curve based on the kinetic coefficients given.
#' @examples
#' # Recreation Figure 9 of Gleaves TAP-2
#'
#' x = (1:5000) * 0.001
#' tau = 1
#' ka = 0
#' kd = 0
#'
#' flux1 = sdc(x, tau = tau,
#'   ka = ka, kd = kd
#' )
#'
#' ka = 20
#' kd = 20
#' flux2 = sdc(x, tau = tau,
#'   ka = ka, kd = kd
#' )
#'
#' kd = 5
#' flux3 = sdc(x, tau = tau,
#'   ka = ka, kd = kd
#' )
#'
#' plot(x, flux1, type = "l", ylim = c(-0.1, 2.4),
#' xlab="Time", ylab = "Flux", main = "Black:ka =0, kd =0; Blue: ka = 20, kd = 20; Red: ka = 20, kd=5")
#' lines(x, flux2, col = "blue")
#' lines(x, flux3, col = "red")
#' @export sdc

sdc = function(x, tau = 1, ka = 0, kd = 0, n_max = 2000) {
  # ensuring all values are positive
  tau = abs(tau)
  ka = abs(ka)
  kd = abs(kd)
  # time must be strictly greater than 0
  x = x[x > 0]

  # Most quantities will be transformed to a matrix for computation speed.
  x = matrix(x, nrow = 1)
  # equivalent to 2n + 1
  pn = ((0:n_max) + 0.5) * pi
  sign_n = 2 * tau * c(1, rep(c(-1, 1), ceiling(n_max / 2)))[1:length(pn)] * pn
  kn = pn^2 + ka + kd
  r_sqrt = sqrt(kn^2 - 4 * pn^2 * kd)
  An = ((-kn + r_sqrt) / 2 + kn - kd) / r_sqrt

  beta_neg = matrix(sign_n * An, nrow = 1)
  beta_pos = matrix(sign_n * (1 - An), nrow = 1)
  exp_neg = exp(tau * matrix((-kn - r_sqrt) / 2, ncol = 1) %*% x)
  exp_pos = exp(tau * matrix((-kn + r_sqrt) / 2, ncol = 1) %*% x)

  flux = abs(as.numeric(beta_neg %*% exp_neg + beta_pos %*% exp_pos))
  return(flux)
}
