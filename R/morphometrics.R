#' Watershed Morphometric Analysis
#'
#' Computes various morphometric parameters for a watershed.
#'
#' @param A Numeric. Area of the watershed (in square km).
#' @param Lb Numeric. Watershed length (in km).
#' @param stream_orders Data frame. Stream orders with columns 'order', 'Nu' (number of streams), 'Lu' (total stream length in km).
#' @param P Numeric. Perimeter of the watershed (in km).
#' @param h Numeric. Minimum elevation (in m).
#' @param H Numeric. Maximum elevation (in m).
#' @param mean_elev Numeric. Mean elevation (in m).
#'
#' @return A list containing computed morphometric parameters.
#' @export
#'
#' @examples
#' stream_orders <- data.frame(
#'   order = c(1, 2, 3),
#'   Nu = c(50, 20, 5),
#'   Lu = c(10, 25, 30)
#' )
#' morphometric_analysis(A = 100, Lb = 15, stream_orders, P = 50, h = 100, H = 500, mean_elev = 300)
morphometric_analysis <- function(A, Lb, stream_orders, P, h, H, mean_elev) {

  # Validate input
  if (!all(c("order", "Nu", "Lu") %in% names(stream_orders))) {
    stop("stream_orders must contain columns 'order', 'Nu', 'Lu'.")
  }

  # Ensure the stream orders are sorted
  stream_orders <- stream_orders[order(stream_orders$order), ]

  # Calculate stream-based parameters
  sum_Nu <- sum(stream_orders$Nu)
  sum_Lu <- sum(stream_orders$Lu)
  highest_order <- max(stream_orders$order)

  F <- sum_Nu / A
  Dd <- sum_Lu / A
  Di <- sum_Nu / sum_Lu
  Ccm <- 1 / Dd
  In <- F * Dd

  n_orders <- nrow(stream_orders)
  Rb <- Rbm <- Rl <- Rlm <- NA

  # Bifurcation Ratios
  if (n_orders >= 2) {
    Rb <- sapply(1:(n_orders - 1), function(i) {
      stream_orders$Nu[i] / stream_orders$Nu[i + 1]
    })
    Rbm <- mean(Rb)
  }

  # Stream Length Ratios
  if (n_orders >= 2) {
    Rl <- sapply(1:(n_orders - 1), function(i) {
      stream_orders$Lu[i + 1] / stream_orders$Lu[i]
    })
    Rlm <- mean(Rl)
  }

  Rho <- if (n_orders >= 2) Rlm / Rbm else NA
  Lo <- 1 / (2 * Dd)
  Dt <- sum_Nu / P  # Drainage texture

  # Shape parameters
  Ff <- A / Lb^2
  Rc <- (4 * pi * A) / P^2
  Cc <- P / (2 * sqrt(pi * A))
  Re <- (2 * sqrt(A / pi)) / Lb
  K <- (Lb^2) / (4 * A)
  Sb <- 1 / Ff

  # Relief parameters
  B <- H - h /1000
  Rr <- (H * 100) / P
  Rn <- B * Dd

  # Hypsometric analysis
  E <- (mean_elev - h) / (H - h)

  # Compile results
  list(
    drainage_density = Dd,
    stream_frequency = F,
    drainage_intensity = Di,
    constant_channel_maintenance = Ccm,
    infiltration_number = In,
    bifurcation_ratios = Rb,
    mean_bifurcation_ratio = Rbm,
    stream_length_ratios = Rl,
    mean_stream_length_ratio = Rlm,
    rho_coefficient = Rho,
    length_overland_flow = Lo,
    drainage_texture = Dt,
    form_factor = Ff,
    circularity_ratio = Rc,
    compactness_coefficient = Cc,
    elongation_ratio = Re,
    lemniscate_ratio = K,
    shape_index = Sb,
    relief = B,
    relative_relief = Rr,
    ruggedness_number = Rn,
    hypsometric_index = E,
    highest_stream_order = highest_order,
    total_stream_number = sum_Nu
  )
}
