function mcmc_fit_ln_like, pars, model_funct = model_funct, x =x, y = y, _extra = _extra
compile_opt idl2
  n_par = n_elements(pars)
  m = call_function(model_funct,x,pars[0:n_par-2])
  sigma = pars[n_par -1] + 1d-9
  li = -0.5*(alog(2d*!dpi*sigma^2) + (y - m)^2/(sigma^2))
  ;li =  alog(1d/sqrt(2d*!dpi*sigma^2)) -(y-m)^2/(2d*sigma^2)
  ;li =  -0.5d*(alog(2d*!dpi) + alog(sigma^2)) -(y-m)^2/(2d*sigma^2)
  return, total(li)
end