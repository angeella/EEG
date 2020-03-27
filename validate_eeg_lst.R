validate_eeg_lst <- function(x, recursive = TRUE) {
  if (!is_eeg_lst(x)) {
    warning("Class is not eeg_lst", call. = FALSE)
  }
  
  if (recursive) {
    x$.signal <- validate_signal_tbl(x$.signal)
    x$.events <- validate_events_tbl(x$.events)
    x$.segments <- validate_segments(x$.segments)
  }
  diff_channels <- setdiff(x$.events$.channel, channel_names(x))
  if (length(diff_channels) != 0 & any(!is.na(diff_channels))) {
    warning("Unknown channel in table of events",
            call. = FALSE
    )
  }
  if (!all.equal(unique(x$.signal$.id), unique(x$.segments$.id))) {
    warning("The values of .id mismatch between tables.",
            call. = FALSE
    )
  }
  
  if (any(!dplyr::group_vars(x) %in% c(colnames(x$.signal), colnames(x$.segments)))) {
    warning("Grouping variables are missing.",
            call. = FALSE
    )
  }
  ## nulls should be caught by the recursive=TRUE
  if (!is.null(sampling_rate(x$.events)) && !is.null(sampling_rate(x$.signal)) &&
      sampling_rate(x$.events) != sampling_rate(x$.signal)) {
    warning("Sampling rates in events and signal table are different",
            call. = FALSE
    )
  }
  x
}
