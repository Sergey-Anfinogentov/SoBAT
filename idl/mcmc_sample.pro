

;+
; :Description:
;    Samples a given function (prob_fun) using Metropolis-Hastings sampler.
;
; :Params:
;    start -dblarr(n_params), starting point in the parameter space
;    prob_fun - a user supplied function to sample, must accept a parameter vector
;    n_samples - number of samples to retrieve
;
; :Keywords:
;    _extra - keywords that will be passed to the function "prob_fun"
;    sigma - dblarr(n_params)- Initial size of the random walk neighbourhood, will be tuned during sampling
;
; :Author: Sergey Anfinogentov
;-
function mcmc_sample, start, prob_fun, n_samples, _extra = _extra, sigma0 = sigma0, burn_in =  burn_in
compile_opt idl2  

  settings = mcmc_settings()
  
  n_par = n_elements(start)

  sigma = identity(n_par)
  ind =where(sigma)
  sigma[ind] = sigma0*100d
  mcmc_randomwalk_update_sigma, start, prob_fun,500, sigma = sigma, _extra = _extra
  s =mcmc_randomwalk(start, prob_fun, 100000, _extra = _extra, sigma = sigma)
  return,s
;  
  
end




