function mcmc_fit_ln_prob, pars, limits = limits, _extra=_extra, ppd_sample = ppd_sample, seed = seed
compile_opt idl2
    
   prior =  mcmc_fit_ln_prior(pars, limits = limits)
   if prior eq -!values.D_INFINITY then return, prior
   like = mcmc_fit_ln_like(pars, _extra=_extra, ppd_sample = ppd_sample, seed = seed)
  return,like + prior
end