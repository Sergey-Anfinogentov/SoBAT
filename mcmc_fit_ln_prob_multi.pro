function mcmc_fit_ln_prob_multi, pars, limits = limits, _extra=_extra
  compile_opt idl2

  prior =  mcmc_fit_ln_prior(pars, limits = limits)
  if prior eq -!values.D_INFINITY then return, prior
  like = mcmc_fit_ln_like_multi(pars, _extra=_extra)
  return,like + prior
end