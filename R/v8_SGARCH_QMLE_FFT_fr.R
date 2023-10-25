fr_FFT <- function(cst, lnyuline_l, a) {
  # FFT-based method
  iT <- length(lnyuline_l)
  np2 <- nextn(2*iT, 2)
  sum_np2 <- fft(fft(c(a, rep(0, np2-iT))) * fft(c(lnyuline_l, rep(0, np2-iT))), inverse = T) / np2
  sum_iT <- cst + sum_np2[1:iT]
  return(Re(sum_iT))
}
