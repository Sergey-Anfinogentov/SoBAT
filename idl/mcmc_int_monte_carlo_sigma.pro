
function mcmc_int_monte_carlo_sigma, funct, mu, sigma, log = log, n_max,_extra=_extra
  if n_elements(limits) eq 2 then limits = reform(limits,1,2)
  seed = random_seed()
  sz = size(limits)
  n_par = sz[1] 
  values = dblarr(n_max)
 
  for i =0, n_max -1 do begin
    par = mcmc_random_multyn(seed,mu,sigma,1)
    g = mcmc_multi_gauss(par, mu, sigma)
    f = call_function(funct,par,_extra=_extra)
     if keyword_Set(log) then f = exp(f)
     if g eq 0 then stop
     mcmc_message,'calculating evidence...'  +string(float(i)/n_max*100.,format = '(I2)') + '%
    values[i]  = f/g
  endfor
  result = total(values,/cum)/(dindgen(n_max)+1d)
 return,total(values/n_max);result[n_max -1]
end