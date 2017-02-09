
Function mcmc_evidence, ln_prob, samples, max_n,limits = limits,_extra = _extra,x =x ,y= y
  n_par = n_elements(limits[*,0]) 
  limits_ = dblarr(n_par+1,2)
  limits_[0:n_par-1,*] = limits
  limits_[n_par,*] = [0, max(y) - min(y)]

  sigma = mcmc_covariance_matrix(samples, mu = mu)

  r2 =  mcmc_int_monte_carlo_sigma(ln_prob,mu,sigma,max_n,_extra = _extra,/log, limits = limits_,x =x ,y= y)

return,r2;[-1]
  


end