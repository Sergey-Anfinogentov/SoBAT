;+
; :Description:
;    Estimates the confidence intervals for the dependend variable Y
;    for a given model and samples from the parameter space generated by the MCMC_FIT routine.
;    It is als possible to simulate observational noise to get confidence intervals for
;    the possibble values of observed values of Y
;
;
; :Params:
;    x - independent variable
;    samples - dblarr(n_params, n_samples), samples from the parameter space, generated by the MCMC_FIT routine
;    model_funct -  a fitted model function. It must accept 2 parameters, X and PARS.
;
; :Keywords:
;    confidence_level - confidence level. It must be a positiv number lower than 1.
;    sigma_samples - dblarr(n_samples), samples of observational noise standard deviation,generated by the
;                   MCMC_FIT routine. This keyword nust be provided for estimating confidence intervals
;                   for the observed value of Y
;    observation - If this keyword is set the confidence intervals for the observed (with observational noise)
;                   value of Y will be estimated. Otherwise, the confidence interval will be estimated for
;                   the undelying modelled values of Y, without observational noise.
;                   To model observational noise the sigma_samples keyword must be also given.
;
; :Author: Sergey Anfinogentov
;-
function mcmc_fit_estimate_y_limits_multi, x, samples, model_funct, confidence_level = confidence_level, sigma_samples = sigma_samples, observation = observation
  compile_opt idl2
  n_funct =  n_elements(model_funct)
  sz = size(samples)
  n_par = sz[1] - n_funct
  n_sampl = sz[2]
  samples_i = dblarr(n_par + 1,n_sampl)
  samples_i[0:-2,*] = samples[0:-1-n_funct,*]
  result = list()
  for i = 0, n_funct -1 do begin
    x_i = x[i]
    sigma_samples_i =  sigma_samples[i,*]

    samples_i[-1,*] = sigma_samples_i
    
   result.add, mcmc_fit_estimate_y_limits( x_i, samples_i, model_funct[i], confidence_level = confidence_level, sigma_samples = sigma_samples_i, observation = observation)
    
  endfor
  
  return, result
end
