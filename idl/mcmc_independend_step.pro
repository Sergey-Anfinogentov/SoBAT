function mcmc_independend_step, seed, current, mu, sigma, prob_fun,accepted = accepted, rejected = rejected, _extra = _extra, value = value
  compile_opt idl2
  n = n_elements(current)
  current_prob = call_function(prob_fun,current,_extra = _extra)
  zeroes = dblarr(n)
  result = current
  for k =0, 199 do begin
    new = mcmc_random_multyn(seed,mu,sigma,1)
    new_prob = call_function(prob_fun,new,_extra = _extra)
    g_old = mcmc_multi_gauss(current, mu, sigma)
    g_new = mcmc_multi_gauss(new, mu, sigma)
    
    ratio = exp(new_prob - current_prob)*g_old/g_new
    if ratio ge 1. then begin
      result = new
      current_prob = new_prob
      accepted += 1l
      break
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