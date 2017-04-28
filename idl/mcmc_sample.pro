

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
function mcmc_sample, start, prob_fun, n_samples, _extra = _extra, sigma0 = sigma0, burn_in =  burn_in;, evidence = evidence
compile_opt idl2  

  settings = mcmc_settings()
  

  
  n_par = n_elements(start)
  if not keyword_set(sigma0) then sigma0 = replicate(1d, n_par)
  
  if not keyword_set(burn_in) then burn_in = 10000l
  
  sigma = identity(n_par)
  ind =where(sigma)
  sigma[ind] = sigma0*100d
  
  print, "-----------------------"
  print, "| IDL MCMC code " + settings.version + ' |'
  print, "-----------------------"
  
  message,"Starting Burn In sampling ("+strcompress(burn_in)+" samples requested)",/info
  
  mcmc_randomwalk_update_sigma, start, prob_fun,500, sigma = sigma, _extra = _extra
  s =mcmc_randomwalk(start, prob_fun, burn_in, _extra = _extra, sigma = sigma)
  start = s[*,burn_in-1]
  
    message,"Burn in finished",/info
    
    
  message,"Starting main sampling ("+strcompress(n_samples)+" samples requested)",/info
  sigma = mcmc_covariance_matrix(s, mu = mu)*10

  mcmc_randomwalk_update_sigma, start, prob_fun,500, sigma = sigma, _extra = _extra
  start = s[*,burn_in-1]
  s =mcmc_randomwalk(start, prob_fun, n_samples, _extra = _extra, sigma = sigma, mu = mu)
  message,"Main sampling finished",/info 
   
  ; s =mcmc_independend(mu, prob_fun, n_samples, _extra = _extra, sigma = sigma, mu = mu)
  ;sigma =  mcmc_covariance_matrix(s, mu = mu)
  ;s =mcmc_independend(mu, prob_fun, n_samples, _extra = _extra, sigma = sigma, mu = mu, evidence = evidence)
  
  return,s
;  
  
end




