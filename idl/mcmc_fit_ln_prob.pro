function mcmc_fit_ln_prob, pars, limits = limits, _extra=_extra, ppd_sample = ppd_sample, seed = seed, priors = priors
compile_opt idl2
   n = n_elements(pars) 
   prior = 0d
   for i=0, n-1 do begin
    prior += priors[i].get_log_value(pars[i])
   endfor
   if prior eq -!values.D_INFINITY then return, prior
   like = mcmc_fit_ln_like(pars, _extra=_extra, ppd_sample = ppd_sample, seed = seed)
  return,like + prior
end