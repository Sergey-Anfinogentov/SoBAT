function mcmc_fit_ln_prior, pars, limits = limits
compile_opt idl2
  n_par = n_elements(pars)
  volume = product(limits[*,1]-limits[*,0]) 
  result = -alog(volume)
  
  for i =0, n_par - 1 do begin
    par = pars[i]
    if (par lt limits[i,0]) or (par gt limits[i,1]) then begin 
    ;    stop
    return, -!values.D_INFINITY 

    endif
  endfor
  return,result
end