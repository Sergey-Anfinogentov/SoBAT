function mcmc_fit_ln_like_multi, pars, model_funct = model_funct, seed = seed, x =x, y = y, _extra = _extra, ppd_sample = ppd_sample
  compile_opt idl2
  n_par = n_elements(pars)
  n = n_elements(x)
  n_ppd = lonarr(n)
  for i = 0, n-1 do begin
    n_ppd[i] = n_elements(x[i])
  endfor
  ppd_ind =total([0,n_ppd],/cum)
  ppd_sample = dblarr(total(n_ppd))
  sigma = pars[n_par -n : n_par -1]>1d-20
  li =0d
  for i = 0, n -1 do begin
    x_i = x[i]
    y_i = y[i]
    m = call_function(model_funct[i],x_i,pars[0:n_par-1-n])
    li += total(-0.5*(alog(2d*!dpi*sigma[i]^2) + (y_i - m)^2/(sigma[i]^2))) 
    ppd_sample[ppd_ind[i]:ppd_ind[i+1]-1] = m + sigma[i] * randomn(seed, n_ppd[i])   
  endfor
 
  ;li =  alog(1d/sqrt(2d*!dpi*sigma^2)) -(y-m)^2/(2d*sigma^2)
  ;li =  -0.5d*(alog(2d*!dpi) + alog(sigma^2)) -(y-m)^2/(2d*sigma^2)
  return, li
end