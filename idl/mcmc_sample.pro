function mcmc_step_MH,seed, current,sigma,prob_fun,accepted = accepted, rejected = rejected, _extra = _extra
compile_opt idl2 
  n = n_elements(current)
  current_prob = call_function(prob_fun,current,_extra = _extra)
  zeroes = dblarr(n)
  for k =0, 1999 do begin
      new = mcmc_random_multyn(seed,current,sigma,1)
      new_prob = call_function(prob_fun,new,_extra = _extra)
      ratio = exp(new_prob - current_prob)
      if ratio ge 1. then begin
        current = new
        current_prob = new_prob
        accepted += 1l
        break
      endif
      rnd = randomu(seed)
      if rnd lt ratio then begin
        current = new
        current_prob = new_prob
        accepted += 1l
        break
      endif
      rejected += 1l
  endfor
  return, current

end


function mcmc_step_burn_in,seed, current,sigma,prob_fun,accepted = accepted, rejected = rejected, _extra = _extra
compile_opt idl2
  
  q = 0.99d
  p = q^(-77d/23d)
 
  n = n_elements(current)
  current_prob = call_function(prob_fun,current,_extra = _extra)
  for i = 0, n-1 do begin
    s = sigma[i]
    for k =0, 9 do begin
      step = s*randomn(seed)
      new = current
      new[i] += step
      new_prob = call_function(prob_fun,new,_extra = _extra)
      ln_ratio =(new_prob - current_prob)
      if ln_ratio lt -50d then ln_ratio = -!values.D_INFINITY 
      if ln_ratio gt 10d then ln_ratio = 10
      
      ratio = exp(ln_ratio)
      
      if ratio ge 1. then begin
        current = new
        current_prob = new_prob
        accepted[i] += 1l
        sigma[i] *= p
        break
      endif
      rnd = randomu(seed)
      if rnd lt ratio then begin
        current = new
        current_prob = new_prob
        accepted[i] += 1l
        sigma[i] *= p
        break
      endif
      rejected[i] += 1l
      sigma[i] *= q
    endfor
  endfor
  return, current

end


function mcmc_sample_burn_in, start, prob_fun, n_samples, _extra = _extra, sigma = sigma
compile_opt idl2  
  n_tune = 100l
  n_par = n_elements(start)
   time_info =2d
  result = dblarr(n_par,n_samples)
  min_steps = 1000l
  target_rate = 0.23
  
  current = start
  seed = systime(1)
  if not keyword_set(sigma) then sigma = current*0d + 1d
  accepted = lonarr(n_par)
  rejected = lonarr(n_par)
  time = systime(1)
  rate = sigma*0d
  n_good = 0l
  for i =0, n_samples -1 do begin
    result[*,i] = mcmc_step_burn_in(seed, current,sigma,prob_fun,accepted = accepted, rejected = rejected, _extra = _extra)
    current = result[*,i]
    if i mod n_tune eq (n_tune -1) then begin
      rate = double(accepted)/(accepted + rejected)
      if i ge  min_steps and product((rate gt target_rate*0.7) and (rate lt target_rate*1.3)) then n_good += n_tune
      if n_good ge 1000l then break
      for k = 0, n_par -1 do begin
        ;sigma[k]*=0.01d^(target_rate - rate[k])
        accepted *= 0l
        rejected *= 0l
      endfor    
    endif
    if systime(1) - time gt time_info then begin
    rate = double(accepted)/(accepted + rejected)
      print,'burning in: '+strcompress(i,/remove_all)+' ('+string(float(i)/n_samples*100.,format = '(I2)'), '%) Acceptance rates: ' ,$
         strcompress(round(rate*100.))+'%'
      ;print,'sigma:',sigma
       time = systime(1)
    endif 
  endfor
  message,'Burning in finished in '+strcompress(i)+' steps',/info
  return, result[*,0:i]
end


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
function mcmc_sample, start, prob_fun, n_samples, _extra = _extra, sigma = sigma, burn_in =  burn_in
compile_opt idl2  
  
  if not keyword_set(burn_in) then burn_in = 10000l
  n_tune = 1000
  n_par = n_elements(start)
  time_info =2d
  target_rate =0.23
  max_covar = 3000l
  
  
  burn_in_samples = mcmc_sample_burn_in(start, prob_fun, burn_in, _extra = _extra, sigma = sigma)
  ;stop
  ;Sampling=====================================================================================================
  n_burn = (size(burn_in_samples))[2]
  sigma = mcmc_covariance_matrix(burn_in_samples[*,n_burn - 1001:n_burn -1], mu = mu)*2.38d/sqrt(double(n_par))
  accepted = 0l
  rejected = 0l
  rate = 0d
  current = burn_in_samples[*,n_burn-1]
  result = dblarr(n_par,n_samples)
  seed = random_seed()
  time = systime(1)
  
  k_sigma = 2.38d/sqrt(double(n_par))
  
  for i =0, n_samples -1 do begin
    result[*,i] = mcmc_step_MH(seed, current,sigma,prob_fun,accepted = accepted, rejected = rejected, _extra = _extra)
    current = result[*,i]
    if i mod n_tune eq (n_tune -1) then begin
      rate = double(accepted)/double(accepted + rejected)
      k_sigma*=0.05d^(target_rate - rate)
      
        en = i
        st = (i - max_covar)>0l
      
       sigma = mcmc_covariance_matrix(result[*,st:en], mu = mu)*k_sigma 
       accepted *= 0l
       rejected *= 0l
    endif
    if systime(1) - time gt time_info then begin
    if (accepted + rejected) gt 100 then rate = (double(accepted)/double(accepted + rejected))>0d
      print,'Sampling: '+strcompress(i,/remove_all)+'('+string(float(i)/n_samples*100.,format = '(I2)'), '%) Acceptance rates: ' ,$
         string(rate*100.,format = '(F4.1)')+'%, sigma coefficient:',string(k_sigma,format = '(F4.2)')
      ;print,sigma
       time = systime(1)
    endif 
  endfor
  
  
  return, result
end





function mcmc_sample_old, start, prob_fun, n_samples, _extra = _extra, sigma = sigma
compile_opt idl2  
  n_tune = 100l
  n_par = n_elements(start)
  
  result = dblarr(n_par,n_samples)
  
  current = start
  seed = systime(1)
  if not keyword_set(sigma) then sigma = current*0d + 1d
  accepted = lonarr(n_par)
  rejected = lonarr(n_par)
  time = systime(1)
  rate = sigma*0d
  for i =0, n_samples -1 do begin
    result[*,i] = mcmc_step(seed, current,sigma,prob_fun,accepted = accepted, rejected = rejected, _extra = _extra)
    current = result[*,i]
    if i mod n_tune eq (n_tune -1) then begin
      rate = double(accepted)/(accepted + rejected)
      for k = 0, n_par -1 do begin
        sigma[k]*=0.5d^(0.6d - rate[k])
        accepted *= 0l
        rejected *= 0l
      endfor    
    endif
    if systime(1) - time gt 2d then begin
      print,'Sampling: ',strcompress(round(float(i)/n_samples*100.)), '% Acceptance rates: ' ,$
         strcompress(round(rate*100.))+'%';,i
      ;print,sigma
       time = systime(1)
    endif 
  endfor
  return, result
end