function mcmc_fit_evidence, samples, x, y, limits, model_funct, n_iterations = n_iterations
  if not keyword_set(n_iterations) then n_iterations =10000l
  evidence = mcmc_evidence('mcmc_fit_ln_prob',samples, n_iterations, x=x,y=y,limits = limits,model_funct =  model_funct)
  return, evidence
end