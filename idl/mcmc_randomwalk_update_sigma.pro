function mcmc_randomwalk_update_sigma_step, seed, current,sigma,prob_fun,accepted = accepted, rejected = rejected, _extra = _extra, value = value
  compile_opt idl2
  n = n_elements(current)
  current_prob = call_function(prob_fun,current,_extra = _extra)
  zeroes = dblarr(n)
  result = current
  for k =0, 99 do begin
    new = mcmc_random_multyn(seed,current,sigma,1)
    new_prob = call_function(prob_fun,new,_extra = _extra)
    ratio = exp(new_prob - current_prob)
    if ratio ge 1. then begin
      ratio = 1d/ratio
    endif
    rnd = randomu(seed)
    if rnd lt ratio then begin
      result = new
      current_prob = new_prob
      accepted += 1l
      break
    endif
    rejected += 1l
  endfor
  value = current_prob
  return, result
end

;+
    ; :Description:
    ;    Estimates the most optimal proposal distribution covariance matrix
    ;    for current location in the parameter space
    ;
    ; :Params:
    ;    start - location in the parameter space
    ;    prob_fun - probability  density function to be sampled
    ;    n_samples - number of samples used for ajustment of the proposal distribution
    ;
    ; :Keywords:
    ;    sigma - initial guess of the covariance matrix. It will contain the estimated
    ;             optimal covariance matrix
    ;    _extra
    ;
    ; :Author: sergey
    ;-
pro mcmc_randomwalk_update_sigma, start, prob_fun, sigma = sigma, _extra = _extra
compile_opt idl2
  settings = mcmc_settings()
  n_samples = settings.proposal_tune_samples
  current = start
  seed = random_seed()
  n_par = n_elements(start)
  samples = dblarr(n_par,n_samples)
  if not keyword_set(sigma) then sigma = identity(n_par)
  accepted_total =0l
  accepted = 0l
  rejected = 0l
  k = 1d
  ind = where(identity(n_par))
  for i =0, n_samples -1 do begin
    
    samples[*,i] = current - mcmc_randomwalk_update_sigma_step(seed, current,sigma, prob_fun,accepted = accepted, rejected = rejected, _extra = _extra)
    rate = double(accepted)/double(accepted + rejected)
    accepted_total += accepted
    
    if accepted_total gt n_par*3 then begin
      sigma =  mcmc_covariance_matrix(samples[*,0:i], mu = mu);stddev(samples[*,0:i]);
      
    endif
    if rate ge 0.5 then sigma*= 1.1d
    if accepted eq 0 then sigma*= 0.5d
    mcmc_message,'Estimating proposal distribution: '+strcompress(i,/remove_all)+'('+string(float(i)/n_samples*100.$
      ,format = '(I2)') + '%)';+' Acceptance rate: ' +string(rate*100.,format = '(F5.1)') + '%',0.5d
    accepted = 0
    rejected = 0l
  endfor
  sigma = mcmc_covariance_matrix(samples, mu = mu);*k



end