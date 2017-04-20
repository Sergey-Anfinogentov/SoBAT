


function mcmc_randomwalk, start, prob_fun, n_samples, _extra = _extra, sigma = sigma
compile_opt idl2
  settings = mcmc_settings()
  n_par = n_elements(start)
  result = dblarr(n_par,n_samples)
  current = start
  seed = random_seed()
  time = systime(1)  
  accepted = lonarr(settings.acceptance_buffer_size)
  rejected = lonarr(settings.acceptance_buffer_size)
  i = 0l
  while i lt n_samples do begin
    rejected_i = 0l
    accepted_i = 0l
    result[*,i] =  mcmc_randomwalk_step(seed, current,sigma,prob_fun,accepted = accepted_i, rejected = rejected_i, _extra = _extra, value = value)
    current = result[*,i]
    
    ;calculaton of the acceptance rate
    accepted = shift(accepted,1)
    rejected = shift(rejected,1)
    accepted[0] = accepted_i
    rejected[0] = rejected_i
    rate = (double(total(accepted))/double(total(accepted + rejected)))>0d
    
    ;Printng out diagnostic information
    if systime(1) - time gt settings.printing_interval then begin
      
      message,'Sampling: '+strcompress(i,/remove_all)+'('+string(float(i)/n_samples*100.,format = '(I2)') + '%) Acceptance rate: ' +$
        string(rate*100.,format = '(F5.1)') + '%, current log value: '+strcompress(value),/info
      time = systime(1)
    endif
    
    ;Check the efficiency
    if  (rate le 0.1) and (total(rejected + rejected) ge settings.acceptance_buffer_size) then begin
      message,"Acceptance rate is too low, tuning the proposal distribution  and restarting the chain...",/info
      sigma *= 1000d
      mcmc_randomwalk_update_sigma, current, prob_fun,500, sigma = sigma, _extra = _extra
      accepted *= 0
      rejected *= 0
      i = 0l
      continue  
    endif
    if n_par eq 1 then rate *= 0.5d
    if  (rate gt 0.4) and (total(accepted + rejected) ge settings.acceptance_buffer_size) then begin
      message,"Acceptance rate is too high, tuning the proposal distribution and restarting the chain...",/info
      sigma *= 1000d
      mcmc_randomwalk_update_sigma, current, prob_fun,500, sigma = sigma, _extra = _extra
      accepted *= 0
      rejected *= 0
      i = 0l
      continue
    endif 
    i += 1l
  endwhile 
  return, result
end