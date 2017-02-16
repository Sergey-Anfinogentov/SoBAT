function mcmc_fit_evidence, samples, x, y, limits, model_funct, n_iterations = n_iterations
  if n_elements(model_funct) gt 1 then begin
    return, mcmc_fit_evidence_multi(samples, x, y, limits, model_funct, n_iterations = n_iterations)
  endif

  if not keyword_set(n_iterations) then n_iterations =10000l
  evidence = mcmc_evidence('mcmc_fit_ln_prob',samples, n_iterations, x=x,y=y,limits = limits,model_funct =  model_funct)
  return, evidence
end