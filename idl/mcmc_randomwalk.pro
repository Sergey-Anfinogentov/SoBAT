;+
    ; :Description:
    ;    Samples a given function (prob_fun) using the random walk Metropolis-Hastings sampler
    ;    using the multivariative normal distribution with a given covariance matrix as
    ;    a proposal distribution
    ;
    ; :Params:
    ;    start - dblarr(n_params), starting point in the parameter space
    ;    prob_fun - a user supplied function to sample, must accept a parameter vector
    ;    n_samples - number of samples to retrieve
    ;
    ; :Keywords:
    ;    _extra
    ;    sigma - covariance matrix of the proposal distribution
    ;    ppd_samples - dblarr(n_data, n_samples) will contain samples from the Posteriour Predictive Distribution
    ;    valuess - dblarr(n_samples) will contain values from the target function corresponding to the samples
    ;
    ; :Author: sergey
    ;-
function mcmc_randomwalk, start, prob_fun, n_samples, _extra = _extra, sigma = sigma, ppd_samples = ppd_samples, values = values
compile_opt idl2
  settings = mcmc_settings()
  min_rate = settings.min_acceptance_rate
  max_rate = settings.max_acceptance_rate
  seed = random_seed()
  n_par = n_elements(start)

  result = dblarr(n_par,n_samples)
  values = dblarr(n_samples)

  
  current = start

  accepted = lonarr(settings.acceptance_buffer_size)
  rejected = lonarr(settings.acceptance_buffer_size)
  
  ;Get information about ppd_samples
  foo = mcmc_randomwalk_step(seed, current,sigma,prob_fun, accepted = 0, rejected = 0, _extra = _extra,  ppd_sample = ppd_sample, out_prob = current_prob)
  n_data = n_elements(ppd_sample)
  if settings.no_ppd_samples then n_data =0
  if n_data gt 0 then begin
    if (settings.double_ppd_samples) then begin
      ppd_samples = dblarr(n_data, n_samples)
    endif else begin
      ppd_samples = fltarr(n_data, n_samples)
    endelse
    
  endif
  i = 0l
  while i lt n_samples do begin
    rejected_i = 0l
    accepted_i = 0l
    result[*,i] =  mcmc_randomwalk_step(seed, current,sigma,prob_fun,accepted = accepted_i, rejected = rejected_i, _extra = _extra, out_prob = value, ppd_sample = ppd_sample, current_prob = current_prob)
    current_prob = value
    values[i] = current_prob
    
    current = result[*,i]
    
    if n_data gt 0 then ppd_samples[*,i] = ppd_sample
    
    ;calculaton of the acceptance rate
    accepted = shift(accepted,1)
    rejected = shift(rejected,1)
    accepted[0] = accepted_i
    rejected[0] = rejected_i
    rate = (double(total(accepted))/double(total(accepted + rejected)))>0d
    
    ;Printng out diagnostic information
      
      mcmc_message,'Sampling: '+strcompress(i,/remove_all)+'('+string(float(i)/n_samples*100.,format = '(I2)') + '%) Acceptance rate: ' +$
        string(rate*100.,format = '(F5.1)') + '%, current log value: '+strcompress(current_prob)
    
    ;Check the efficiency
    if  (rate le min_rate) and (total(rejected + rejected) ge settings.acceptance_buffer_size) then begin
      message,"Acceptance rate is too low, tuning the proposal distribution  and restarting the chain...",/info
      sigma *= 1000d
      mcmc_randomwalk_update_sigma, current, prob_fun, sigma = sigma, _extra = _extra
      accepted *= 0
      rejected *= 0
      i = 0l
      continue  
    endif
    if  (rate gt max_rate) and (total(accepted + rejected) ge settings.acceptance_buffer_size) then begin
      message,"Acceptance rate is too high, tuning the proposal distribution and restarting the chain...",/info
      sigma *= 1000d
      mcmc_randomwalk_update_sigma, current, prob_fun, sigma = sigma, _extra = _extra
      accepted *= 0
      rejected *= 0
      i = 0l
      continue
    endif 
    i += 1l
  endwhile 
  return, result
end