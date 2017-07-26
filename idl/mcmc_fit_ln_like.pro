function mcmc_fit_ln_like, pars, model_funct = model_funct,seed = seed, x =x, y = y, _extra = _extra, ppd_sample = ppd_sample
compile_opt idl2
  n_par = n_elements(pars)
  n_data = n_elements(y)
  m = call_function(model_funct,x,pars[0:n_par-2])
  sigma = pars[n_par -1] + 1d-20
  li = -0.5*(alog(2d*!dpi*sigma^2) + (y - m)^2/(sigma^2))
  ppd_sample = m + sigma * randomn(seed, n_data)
  return, total(li)
end