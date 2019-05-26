function mcmc_fit_evidence, samples, x, y, priors, model_funct, n_iterations = n_iterations, errors = errors, _extra = _extra
  if n_elements(model_funct) gt 1 then begin
    return, mcmc_fit_evidence_multi(samples, x, y, priors, model_funct, n_iterations = n_iterations, errors = errors, _extra = _extra)
  endif

  if not keyword_set(n_iterations) then n_iterations =10000l
  evidence = mcmc_evidence('mcmc_fit_ln_prob',samples, n_iterations, x=x,y=y,priors = priors,model_funct =  model_funct, errors = errors, _extra = _extra)
  return, evidence
end