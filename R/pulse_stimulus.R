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


#' Generate pulse sequence
#' 
#' @description generate new pulse sequence
#'
#' @param n integer; number of stimuli (or trials)
#' @param bins integer; number of stimulus bins
#' @param dist string; distribution (or method) used to create stimuli. Either "bernoulli" or "poisson". If "bernoulli" then p_left = 1 - p_right
#' @param p_right numeric; probability of right pulse in a given bin
#' @param p_left numeric; probability of left pulse in a given bin. if NA or dist = "bernoulli", p_left = 1 - p_right
#' @param add_1 logical; if TRUE, add pulse at the beginning of each sequence (e.g. if trial starts with pulse on both sides)
#' @param return_train logical; if TRUE, returns stimulus matrix (result from pulse_stimulus). if FALSE, returns character matrix. Default = FALSE
#' @param ... arguments passed to pulse_stimulus (only used if return_train = TRUE)
#' 
#' @md
#'
#' @return stimulus matrix or character matrix
#'
#' @export
pulse_sequence = function(n=1,
                          bins=25,
                          dist="bernoulli",
                          p_right=0.5,
                          p_left=NA,
                          add_1=FALSE,
                          return_train=FALSE,
                          ...) {
  
  ### checks ###
  
  if (!(dist %in% c("bernoulli", "poisson"))) {
      stop(paste0("dist = ", dist, " is not supported! Must use dist = \"bernoulli\" or dist = \"poisson\""))
  }
  
  if (any(p_right < 0 | p_right > 1)) {
    stop("p_right must be between 0 and 1")
  }
  
  if (length(p_right) != n) {
    if (length(p_right) != 1)
      stop("p_right must be length 1 or n")
    p_right = rep(p_right, n)
  }
  
  if ((dist == "bernoulli") | (length(p_left) == 1 & is.na(p_left))) {
    p_left = 1 - p_right
  } else {
    if (any(p_left < 0 | p_left > 1)) {
      stop("p_left must be between 0 and 1")
    }
    if (!(length(p_left) %in%  c(1,n))) {
      stop("p_left must be length 1 or n")
    }
    if (length(p_left) == 1) {
      p_left = rep(p_left, n)
    }
  }
  
  ### generate stimuli ###
  
  right_seq = sapply(p_right, function(x) sample(c(1, 0), bins, replace=T, prob=c(x, 1-x)))
  if (dist == "bernoulli") {
    left_seq = abs(1 - right_seq)
  } else {
    left_seq = sapply(p_left, function(x) sample(c(1, 0), bins, replace=T, prob=c(x, 1-x)))
  }
  
  if (add_1) {
    right_seq = rbind(rep(1, n), right_seq)
    left_seq = rbind(rep(1, n), left_seq)
  }
  
  right_str = apply(right_seq, 2, function(x) paste(x, collapse=""))
  left_str = apply(left_seq, 2, function(x) paste(x, collapse=""))
  
  if (return_train) {
    res = pulse_stimulus(right_str, left_str, ...)
  } else {
    res = matrix(c(right_str, left_str), ncol=2)
    colnames(res) = c("right", "left")
  }
  
  res
  
}
