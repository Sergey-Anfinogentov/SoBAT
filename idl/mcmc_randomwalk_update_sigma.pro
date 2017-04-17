function mcmc_randomwalk_update_sigma_step, seed, current,sigma,prob_fun,accepted = accepted, rejected = rejected, _extra = _extra, value = value
  compile_opt idl2
  n = n_elements(current)
  current_prob = call_function(prob_fun,current,_extra = _extra)
  zeroes = dblarr(n)
  result = current
  for k =0, 199 do begin
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

pro mcmc_randomwalk_update_sigma, start, prob_fun, n_samples, sigma = sigma, _extra = _extra
compile_opt idl2
  current = start
  seed = random_seed()
  time = systime(1)
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
    ;if rate gt 0.9 then k *= 1.1
    ;if rate lt 0.05 then k *= 0.9
    ;k = k<1d

    
    if accepted_total gt n_par*3 then begin
      sigma =  mcmc_covariance_matrix(samples[*,0:i], mu = mu);stddev(samples[*,0:i]);
      
    endif
        if rate ge 0.5 then sigma*= 1.1d
        if accepted eq 0 then sigma*= 0.5d
    ;print, rate,sigma
   ; wait,0.1
    accepted = 0l
    rejected = 0l
   ; if i mod 100 eq 99 then print, rate
  endfor
  sigma = mcmc_covariance_matrix(samples, mu = mu);*k
;  ;print, rate,'rate'
;  window,0
;  h = histogram(samples[0,*],nbins = 50, locations = loc)
;  plot,loc,double(h)/max(h)
;  window,1
;  h = histogram(samples[1,*],nbins = 50, locations = loc)
;  plot,loc,double(h)/max(h)
  ;n = n_elements(loc)
  ;res = dblarr(n)
  ;for i = 0, n-1 do res[i] = call_function(prob_fun,loc[i],_extra = _extra)
  ; oplot,loc,exp(res), linest = 1
;  stop
;  print, sigma
;  print, '---'
;  print, mu
;  print, '---'


end