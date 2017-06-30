
;+
  ; :Description:
  ;    Uses Byesian Inference and MCMC to fit a user suplied function to uder suplied points (X, Y) by adjasting a set of parameters (PARS)
  ;
  ; :Params:
  ;    x - independent variable. For multifuction fiting should be a list of arrays: list(x1, x2, x3)
  ;    y - measurments of the dependent variable. For multifuction fiting should be a list of arrays: list(y1, y2, y3)
  ;    pars - dblarr(n_params), Starting guess, will containt the fitted parameter
  ;           If the starting guess is not set a random starting point will be generated
  ;    limits - dblarr(n_params, 2) possible limitspars for the parameters,
  ;           will contain the confidence intervals for each parameter
  ;    model_funct - a model function to fit. It must accept 2 parameters, X and PARS.For multifuction fiting should be an array of strings ['funct1', 'funct2', 'funct3']
  ;
  ; :Keywords:
  ;    n_samples - numer of samples to generate using Metropolis-Hastings MCMC, default 10000l
  ;    burn_in - number of burn in samples, needed for the sampler to find the high probability region
  ;             and to tune sampling parameters
  ;    samples - dblarr(n_params, n_samples) will contain samples from the Posterior  Distribution
  ;    ppd_samples - dblarr(n_data_points, n_samples) will contain samples from the Posteriour Predictive Distribution
  ;    confidence_level - confidence level to define confidence intervals for each parameter.
  ;    sigma_samples - samples of the standart deviasion of the observational noise which is assumed to be  normally distributed.
  ;                 For multifunction fitting will contain samples of all sigmas
  ;
  ; :Author: Sergey Anfinogentov (sergey.istp@gmail.com)
  ;-
function mcmc_fit,x,y,pars, limits ,model_funct,n_samples = n_samples, sigma_samples = sigma_samples, burn_in = burn_in,$
   samples = samples, ppd_samples = ppd_samples, confidence_level = confidence_level,noise_limits = noise_limits,  _extra = _extra
compile_opt idl2
  
  if n_elements(model_funct) gt 1 then begin
    return, mcmc_fit_multi(x,y,pars, limits ,model_funct,n_samples = n_samples, sigma_samples = sigma_samples, burn_in = burn_in,$
            samples = samples, confidence_level = confidence_level,noise_limits = noise_limits,  _extra = _extra)
  endif

 
  if not keyword_set(n_samples) then n_samples = 10000l
  if not keyword_set(burn_in) then burn_in = 5000l
  if not keyword_set(confidence_level) then confidence_level = 0.95d
  if not keyword_set(pars) then pars = mcmc_fit_random_start(limits)
  n_par = n_elements(pars) 
  
  ;initial guess for sigma
  y_guess = call_function(model_funct, x, pars)
  
  if not keyword_set(noise_limits) then  noise_limits = [0, max(y) - min(y)]
  noise_guess = stddev(y-y_guess)<noise_limits[1]>noise_limits[0]
  
  
  
  
  pars_ = [pars,noise_guess]
  limits_ = dblarr(n_par+1,2)
  limits_[0:n_par-1,*] = limits
  limits_[n_par,*] = noise_limits
  
  sigma = (max(limits_,dim = 2) - min(limits_,dim = 2))/2d
  
  samples = mcmc_sample(pars_,'mcmc_fit_ln_prob',n_samples, burn_in =  burn_in, x = x, y = y,$
     model_funct = model_funct, limits = limits_, sigma = sigma, evidence = evidence, ppd_samples = ppd_samples)
  sigma_samples = samples[n_par,*]

  
  ;samples = samples[0:n_par-1,*]
  pars_ = median(samples,dimension = 2)
  pars = pars_[0:n_par - 1]
  
  dc = (1d - confidence_level)*0.5d
  for i =0, n_par -1 do begin
    limits[i,*] = cgpercentiles(samples[i,*],percentiles = [dc,1d - dc]) 
  endfor
  
  
  return, call_function(model_funct,x,pars)  
end