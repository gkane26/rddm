#' Get pulse stimulus
#'
#' @description get full pulse stimulus
#'
#' @param up_sequence string; string timepoint by timepoint evidence to upper boundary
#' @param down_sequence string; string timepoint by timepoint evidence to lower boundary (optional)
#' @param pattern string; pattern to separate timepoints. Default = ""
#' @param dur numeric; duration of pulse, default = .01
#' @param isi numeric; inter-pulse interval, default = .1
#' @param pre_stim numeric; time before pulse within stimulus bin. default = 0
#' @param dt numeric; timestep for simulation, default = .001
#' @param as_array logical; if True, return as 3d array (2 x timepoints x trials)
#'
#' @return stimulus train
#'
#' @export
pulse_stimulus <- function(up_sequence, down_sequence=NULL, pattern="", dur=0.01, isi=0.1, pre_stim=0, dt=0.001, as_array=F){
  
  up_vec = stringr::str_split(up_sequence, pattern)
  up_vec = lapply(up_vec, as.numeric)
  up_stim_seq = lapply(up_vec, pulse_trial_stimulus, dur=dur, isi=isi, pre_stim=pre_stim, dt=dt)
  
  if (is.null(down_sequence)) {
    down_stim_seq = lapply(up_stim_seq, function(x) matrix(0, nrow=1, ncol=length(x)))
  } else {
    down_vec = stringr::str_split(down_sequence, pattern)
    down_vec = lapply(down_vec, as.numeric)
    down_stim_seq = lapply(down_vec, pulse_trial_stimulus, dur=dur, isi=isi, pre_stim=pre_stim, dt=dt)
  }
  
  if (as_array) {
    stim_lengths = sapply(up_stim_seq, length)
    res = array(NA, dim=c(2, max(stim_lengths), length(up_stim_seq)))
    for (i in 1:length(up_stim_seq)) {
      res[, 1:stim_lengths[i], i] = rbind(up_stim_seq[[i]], down_stim_seq[[i]])
    }
  } else {
    res = list()
    for (i in 1:length(up_stim_seq)) {
      res[[i]] = rbind(up_stim_seq[[i]], down_stim_seq[[i]])
    }
  }
  
  res
  
}
