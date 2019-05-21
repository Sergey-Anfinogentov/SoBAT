function mcmc_fit_start_value, priors
compile_opt idl2
  n = n_elements(priors)
  for i =0, n-1 do begin
    result[i] = priors[i]->get_start_value()
  endfor
  return , result
end