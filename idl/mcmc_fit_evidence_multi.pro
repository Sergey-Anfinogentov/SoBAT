function mcmc_fit_evidence_multi, samples, x, y, limits, model_funct, n_iterations = n_iterations
  if not keyword_set(n_iterations) then n_iterations =10000l
  
  n_funct = n_elements(model_funct)
  n_par = n_elements(limits[*,0])
  limits_ = dblarr(n_par+n_funct,2)
  limits_[0:n_par-1,*] = limits
  
  for i = 0, n_funct -1 do begin
    limits_[n_par+i,*] = [0, max(y[i]) - min(y[i])]
  endfor


  sigma = mcmc_covariance_matrix(samples, mu = mu)

  r2 =  mcmc_int_monte_carlo_sigma('mcmc_fit_ln_prob_multi',mu,sigma,n_iterations,_extra = _extra,/log, limits = limits_,x =x ,y= y,model_funct = model_funct)

  return,r2;[-1]
end