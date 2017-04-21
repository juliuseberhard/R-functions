################################################################################
# Some Functions
# Author: Julius Eberhard
# Last Edit: 2017-04-21
################################################################################

# Package Check and Installation -----------------------------------------------
# reasonable for multiple packages only

CheckPack <- function(packages  # vector of packages
                      ) {

  for (i in 1:length(packages))
    if (!packages[i] %in% installed.packages()[, 1])
      install.packages(packages[[i]])

}

# Leap Year Check --------------------------------------------------------------
# (Author: Forester, quantitative-ecology.blogspot.de)

IsLeapyear <- function(year
                       ) {
  
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))

}

# Soil Water Pedotransfer Function: Woesten et al. 1999 - theta_s --------------
# (Author: Julien Moeys, https://r-forge.r-project.org/projects/soilwater/)

ptf.wosten.theta.s <- function(clay,  # clay content in %
                               bulkD,  # bulk density in kg.L-1
                               silt,  # silt content in %
                               om,  # soil organic matter in %
                               topSoil  # vector of [0=subsoil, 1=topsoil]
                               ) {
  
  # topsoil is usually defined as soil at a depth of less than 30 cm
  # (Dai et al., 2013)
  
  if (any(silt == 0)) {
    stop("Some 'silt' values are 0. The PTF does not handle 0% silt content (ln(silt)).")
  }
  theta.s <- 0.7919 + 0.001691 * clay - 0.29619 * bulkD - 1.491e-06 *
    (silt^2) + 8.21e-05 * (om^2) + 0.02427 * (clay^-1) +
    0.01113 * (silt^-1) + 0.01472 * logb(silt, base = exp(1)) -
    7.33e-05 * om * clay - 0.000619 * bulkD * clay - 0.001183 *
    bulkD * om - 0.0001664 * topSoil * silt
  return(theta.s)

}

# Soil Water Pedotransfer Function: Rawls & Brakensiek 1985 --------------------
# (Author: Julien Moeys, https://r-forge.r-project.org/projects/soilwater/)

pft.rawls <- function (soilprop,
                       parameters = NULL,
                       h = NULL
                       ) {

  ret_val = NULL
  sa = 100 - soilprop$silt - soilprop$clay
  p = 1 - soilprop$bulkD/2.65
  if ("theta" %in% parameters || !is.null(h)) {
    rawls_coefs = rbind(cbind(40, 0.7899, -0.0037, 0, 0.01, -0.1315),
                        cbind(70, 0.7135, -0.003, 0.0017, 0, -0.1693),
                        cbind(100, 0.4118, -0.003, 0.0023, 0.0317, 0),
                        cbind(200, 0.418, -0.0021, 0.0035, 0.0232, -0.0859),
                        cbind(330, 0.3486, -0.0018, 0.0039, 0.0228, -0.0738),
                        cbind(600, 0.2819, -0.0014, 0.0042, 0.0216, -0.0612),
                        cbind(1000, 0.2352, -0.0012, 0.0043, 0.0202, -0.0517),
                        cbind(2000, 0.1837, -9e-04, 0.0044, 0.0181, -0.0407),
                        cbind(4000, 0.1426, -7e-04, 0.0045, 0.016, -0.0315),
                        cbind(7000, 0.1155, -5e-04, 0.0045, 0.0143, -0.0253),
                        cbind(10000, 0.1005, -4e-04, 0.0045, 0.0133, -0.0218),
                        cbind(15000, 0.0854, -4e-04, 0.0044, 0.0122, -0.0182))
      dimnames(rawls_coefs)[[2]] = c("head_cm", letters[1:5])
      coefs = apply(rawls_coefs[, -1], MARGIN = 2, FUN = approx,
                    x = rawls_coefs[, "head_cm"], xout = h)
      coefs = unlist(coefs)[(1:5) * 2]
      names(coefs) = letters[1:5]
      theta = coefs["a"] + coefs["b"] * sa + coefs["c"] * soilprop$clay +
        coefs["d"] * soilprop$om/1.74 + coefs["e"] * soilprop$bulkD
      ret_val = cbind(ret_val, theta = theta)
      rm(theta)
    }
    if ("theta_r" %in% parameters) {
      ret_val = cbind(ret_val, theta_r = rep(NA, nrow(soilprop)))
      ret_val[, "theta_r"] = -0.0182482 + 0.00087269 * sa +
        0.00513488 * soilprop$clay + 0.02939286 * p - 0.00015395 *
        soilprop$clay^2 - 0.0010827 * sa * p - 0.00018233 *
        soilprop$clay^2 * p^2 + 0.00030703 * soilprop$clay^2 *
        p - 0.0023584 * p^2 * soilprop$clay
    }
    if ("S_f" %in% parameters) {
      ret_val = cbind(ret_val, S_f = rep(NA, nrow(soilprop)))
      ret_val[, "S_f"] = exp(6.53 - 7.326 * p + 0.00158 * soilprop$clay^2 +
                               3.809 * p^2 + 0.000344 * sa * soilprop$clay - 0.04989 *
                               sa * p + 0.0016 * sa^2 * p^2 + 0.0016 * soilprop$clay^2 *
                               p^2 - 1.36e-05 * sa^2 * soilprop$clay - 0.00348 *
                               soilprop$clay^2 * p + 0.000799 * sa^2 * p)
      ret_val[, "S_f"] = ret_val[, "S_f"] * 10
    }
    if ("h_b" %in% parameters) {
      ret_val = cbind(ret_val, h_b = rep(NA, nrow(soilprop)))
      ret_val[, "h_b"] = exp(5.3396738 + 0.1845038 * soilprop$clay -
                               2.48394546 * p - 0.00213853 * soilprop$clay^2 - 0.04356349 *
                               sa * p - 0.61745089 * soilprop$clay * p + 0.00143589 *
                               sa^2 * p^2 - 0.00855375 * soilprop$clay^2 * p^2 -
                               1.282e-05 * sa^2 * soilprop$clay + 0.00895359 * soilprop$clay^2 *
                               p - 0.0007247 * sa^2 * p + 5.4e-06 * soilprop$clay^2 *
                               sa + 0.5002806 * p^2 * soilprop$clay)
    }
    if ("lambda" %in% parameters) {
      ret_val = cbind(ret_val, lambda = rep(NA, nrow(soilprop)))
      ret_val[, "lambda"] = exp(-0.7842831 + 0.0177544 * sa -
                                  1.062498 * p - 5.304e-05 * sa^2 - 0.00273493 * soilprop$clay^2 +
                                  1.11134946 * p^2 - 0.03088295 * sa * p + 0.00026587 *
                                  sa^2 * p^2 - 0.00610522 * soilprop$clay^2 * p^2 -
                                  2.35e-06 * sa^2 * soilprop$clay + 0.00798746 * soilprop$clay^2 *
                                  p - 0.00674497 * p^2 * soilprop$clay)
    }
    return(ret_val)

}
