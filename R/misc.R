get_first_passage_density = function(pars_list, dt=.01, ...){
  density = do.call(ddm_integral_fpt, c(pars_list, dt=dt, ...))
  rownames(density) = seq(dt, dt*nrow(density), dt)
  density
}

get_rt_liks = function(dat, density_list, min_p=1e-10){
  
  min_p = as.numeric(min_p)
  
  #get time bin for each response
  tvec = as.numeric(rownames(density_list[[1]]$density))
  dt = tvec[2] - tvec[1]
  dat[, rt_bin := rt / dt]
  
  condition_idx = 2:length(density_list[[1]])
  p_response = numeric()
  new_dat = data.table()
  for(i in 1:length(density_list)){
    
    #subset only this condition
    sub_dat=dat
    for(j in 2:length(density_list[[i]])){
      sub_dat=sub_dat[get(names(density_list[[i]])[j])==density_list[[i]][j]]
    }
    
    sub_dat[response == 0, p_rt := density_list[[i]]$density[rt_bin, 2]]
    sub_dat[response == 1, p_rt := density_list[[i]]$density[rt_bin, 1]]
    new_dat = rbind(new_dat, sub_dat)
  }
  
  new_dat[p_rt < min_p, p_rt := min_p]
  return(new_dat)
}