


function int_monte_carlo_sample,seed, limits
  if n_elements(limits) eq 2 then limits = reform(limits,1,2)
  sz = size(limits)
  n_par = sz[1]
  result = randomu(seed,n_par)
  result *= limits[*,1] - limits[*,0]
  result += limits[*,0]
  return,result
end

function mcmc_int_monte_carlo, funct, limits, log = log, n_max,_extra=_extra
  if n_elements(limits) eq 2 then limits = reform(limits,1,2)
  seed = random_seed()
  sz = size(limits)
  n_par = sz[1]
  volume = product(limits[*,1]-limits[*,0])  
  values = dblarr(n_max)
 
  for i =0, n_max -1 do begin
    values[i]  = call_function(funct,int_monte_carlo_sample(seed,limits),_extra=_extra)
  endfor
  if keyword_Set(log) then values = exp(values)
  result = total(values,/cum)*volume/(dindgen(n_max)+1d)
  window,2
  plot,result
 return,result[n_max -1]
end