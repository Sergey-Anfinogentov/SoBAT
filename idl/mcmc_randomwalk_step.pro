function mcmc_randomwalk_step, seed, current,sigma,prob_fun,accepted = accepted, rejected = rejected, _extra = _extra, out_prob = out_prob,  ppd_sample = ppd_sample, current_prob = current_prob
  compile_opt idl2
  n = n_elements(current)
  if n_elements(current_prob) eq 0  then begin
    current_prob = call_function(prob_fun,current,_extra = _extra, ppd_sample = ppd_sample)
    ;print,'Fired'
  endif
  zeroes = dblarr(n)
  result = current
 ; for k =0, 0 do begin
    new = mcmc_random_multyn(seed,current,sigma,1)
    new_prob = call_function(prob_fun,new,_extra = _extra, ppd_sample = new_ppd_sample, seed = seed)
    ratio = exp(new_prob - current_prob)
    if ratio ge 1. then begin
      result = new
      out_prob = new_prob
      if n_elements(new_ppd_sample) gt 0 then ppd_sample = new_ppd_sample
      accepted += 1l
      return, result
    endif
    rnd = randomu(seed)
    if rnd lt ratio then begin
      result = new
      out_prob = new_prob
      accepted += 1l
      return, result
    endif
    rejected += 1l
    out_prob = current_prob
  return, result
end