


function mcmc_randomwalk, start, prob_fun, n_samples, _extra = _extra, sigma = sigma
compile_opt idl2
  settings = mcmc_settings()
  n_par = n_elements(start)
  result = dblarr(n_par,n_samples)
  current = start
  seed = random_seed()
  time = systime(1)  
  accepted = 0
  rejected = 0
  for i =0, n_samples -1 do begin
    result[*,i] =  mcmc_randomwalk_step(seed, current,sigma,prob_fun,accepted = accepted, rejected = rejected, _extra = _extra, value = value)
    current = result[*,i]
    if systime(1) - time gt settings.printing_interval then begin
      rate = (double(accepted)/double(accepted + rejected))>0d
      print,'Sampling: '+strcompress(i,/remove_all)+'('+string(float(i)/n_samples*100.,format = '(I2)'), '%) Acceptance rates: ' ,$
        string(rate*100.,format = '(F4.1)') + '%'
        print, value
      time = systime(1)
    endif
  endfor 
  print,  (double(accepted)/double(accepted + rejected)), value
  return, result
end