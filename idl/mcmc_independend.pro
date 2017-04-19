function mcmc_independend, start, prob_fun, n_samples, _extra = _extra, sigma = sigma, mu = mu, evidence = evidence
compile_opt idl2
  settings = mcmc_settings()
  n_par = n_elements(start)
  result = dblarr(n_par,n_samples)
  current = start
  seed = random_seed()
  time = systime(1)  
  accepted = 0
  rejected = 0
  evidence = 0
  for i =0, n_samples -1 do begin
    result[*,i] =  mcmc_independend_step(seed, current,mu, sigma,prob_fun,accepted = accepted,$
       rejected = rejected, _extra = _extra, value = value, evidence = evidence)
    current = result[*,i]
    rate = (double(accepted)/double(accepted + rejected))>0d
    ; Printng out diagnostic information
    if systime(1) - time gt settings.printing_interval then begin
      
      Message,'Sampling: '+strcompress(i,/remove_all)+'('+string(float(i)/n_samples*100.,format = '(I2)')+ '%) Acceptance rate: ' +$
        string(rate*100.,format = '(F5.1)') + '%, current log value: '+strcompress(value) +  ', evidence: '+strcompress(evidence/(i+1d)),/info
      time = systime(1)
    endif
    
    ;Check the efficiency
;    if  (rate le 0.1) and (rejected gt 1000) then begin
;      message,"Accpetance rate is too low, tuning the proposal distribution...",/info
;      sigma = mcmc_covariance_matrix(result[*,0:i], mu = mu)
;      ;mcmc_randomwalk_update_sigma, current, prob_fun,500, sigma = sigma, _extra = _extra
;      accepted = 0
;      rejected = 0  
;    endif
;    if  (rate gt 0.5) and (accepted gt 200) then begin
;      message,"Accpetance rate is too high, tuning the proposal distribution...",/info
;      sigma *= 1000d
;      mcmc_randomwalk_update_sigma, current, prob_fun,500, sigma = sigma, _extra = _extra
;      accepted = 0
;      rejected = 0
;    endif
    
    
    
  endfor 
  evidence = evidence / n_samples
  return, result
end