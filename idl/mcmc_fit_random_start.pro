function mcmc_fit_start_value, priors
compile_opt idl2
  sz = size(limits)
  seed = systime(1)
  n_par = sz[1]
  result = dblarr(n_par)
  rnd = randomu(seed,n_par)
  result = rnd*(limits[*,1] - limits[*,0]) + limits[*,0]
  return , result
end