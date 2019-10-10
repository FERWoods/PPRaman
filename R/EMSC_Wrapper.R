#' EMSC wrapper - NOT TESTED USE WITH CAUTION
#'
#' @param spectra Matrix of spectral data to have EMSC applied.
#' @param background_spectra Matrix of spectra to be used as the background component for EMSC. This is then averaged and smoothed.
#' @param reference_spectra Matrix of spectra to be used as the reference component for EMSC. This is then averaged and smoothed.
#' @export
#'

EMSC.PPRAMAN <- function(spectra, background_spectra, reference_spectra){
  #EMSC
  ss_hiplan <- background_spectra ## may have to read in using read_in_no_shift function first
  # Averaging reference hi plan cancer and control
  ref_smooth <- savgol(apply(reference_spectra, 2, mean), fl = 5, forder = 3)
  ref_smooth_snip <- ref_smooth[-(1:((5 + 1)/2))]
  ref_smooth_snip <- ref_smooth_snip[-((length(ref_smooth_snip) - (1+(5 + 1)/2)):length(ref_smooth_snip))]
  # Extrapolating ends of spectra due to w/n loss from smoothing
  ref <- approxExtrap(wavenumber[-c((1:((5 + 1)/2)),((length(ref_smooth_snip) - (1+(5 + 1)/2)):length(ref_smooth_snip))),1],
                      ref_smooth_snip, wavenumber[,1], method = "linear", n = 50, rule = 2, f = 0, ties = "ordered", na.rm = FALSE)$y
  # Background contribution
  # Hplan SS
  hiplanss_bg <- apply(ss_hiplan, 2, mean)
  #nplan <- apply(hiplan_ctrl_cr232, 2, mean) - apply(nplan_ctrl_cr232, 2, mean)
  hiplanss_bg_smooth <- savgol(hiplanss_bg, fl = 5, forder = 3)
  hiplanss_bg_smooth_snip <- hiplanss_bg_smooth [-(1:((5 + 1)/2))]
  hiplanss_bg_smooth_snip <- hiplanss_bg_smooth_snip[-((length(hiplanss_bg_smooth_snip) - (1+(5 + 1)/2)):length(hiplanss_bg_smooth_snip))]
  # Extrapolating ends of spectra due to w/n loss from smoothing
  bg_hiplan_ss <- approxExtrap(wavenumber[-c((1:((5 + 1)/2)),((length(hiplanss_bg_smooth_snip) - (1+(5 + 1)/2)):length(hiplanss_bg_smooth_snip))),1],
                               hiplanss_bg_smooth_snip, wavenumber[,1], method = "linear", n = 50, rule = 2, f = 0, ties = "ordered", na.rm = FALSE)$y
  # EMSC step
  emsc_corrections <- EMSC(spectra, reference = ref, interferent = bg_hiplan_ss, degree = 3)

  return(emsc_corrections$corrected)
}
