
;+
; :Description:
;    Uses Byesian Inference to fit a user suplied function to uder suplied points (X, Y) by adjasting a set of parameters (PARS)
;
; :Params:
;    x - independent variable
;    y - measurments of the dependent variable
;    pars - dblarr(n_params), Starting guess, will containt the fitted parameter
;           If the starting guess is not set a random starting point will be generated
;    limits - dblarr(n_params, 2) possible limitspars for the parameters,
;           will contain the confidence intervals for each parameter
;    model_funct - a model function to fit. It must accept 2 parameters, X and PARS.
;
; :Keywords:
;    n_samples - numer of samples to generate using Metropolis-Hastings MCMC, default 10000l
;    burn_in - number of burn in samples, needed for the sampler to find the high probability region
;             and to tune sampling parameters
;    samples - dblarr(n_params, n_samples) will contain samples from the posteriour PDF
;    confidence_level - confidence level to define confidence intervals for each parameter.
;    sigma_samples - samples of the standart deviasion of the observational noise which is assumed to be  normally distributed
;
; :Author: Sergey Anfinogentov
;-
function mcmc_fit_multi,x,y,pars, limits ,model_funct,n_samples = n_samples, sigma_samples = sigma_samples, burn_in = burn_in,$
  samples = samples, confidence_level = confidence_level,noise_limits = noise_limits,  _extra = _extra
  compile_opt idl2
  if not keyword_set(n_samples) then n_samples = 10000l
  if not keyword_set(burn_in) then burn_in = 5000l
  if not keyword_set(confidence_level) then confidence_level = 0.95d
  if not keyword_set(pars) then pars = mcmc_fit_random_start(limits)
  n_par = n_elements(pars)
  pars_ = [pars,replicate(1d,n_par)]
  
  n_funct = n_elements(model_funct)
  
  
  limits_ = dblarr(n_par+n_funct,2)
  limits_[0:n_par-1,*] = limits
  
  for i = 0, n_funct -1 do begin
    limits_[n_par+i,*] = [0, max(y[i]) - min(y[i])]
  endfor
  
  if keyword_set(noise_limits) then limits_[n_par:npar+n_funct-1,*] = noise_limits

  sigma = (max(limits_,dim = 2) - min(limits_,dim = 2))/6d
  ;sigma = [sigma,(max(y) - min(y))*0.01d]

  samples = mcmc_sample(pars_,'mcmc_fit_ln_prob_multi',n_samples, burn_in =  burn_in, x = x, y = y, model_funct = model_funct, limits = limits_, sigma = sigma)
  ; samples = samples[*, burn_in:*]
  sigma_samples = samples[n_par:*,*]

  ;  stop
  ; if not keyword_set(n_evidence_int) then n_evidence_int =10000l
  ; evidence = mcmc_evidence('mcmc_fit_ln_prob',samples, n_evidence_int, x=x,y=y,limits = limits_,model_funct = model_funct)



  ;samples = samples[0:n_par-1,*]
  pars_ = median(samples,dimension = 2)
  pars = pars_[0:n_par - 1]

  dc = (1d - confidence_level)*0.5d
  for i =0, n_par -1 do begin
    limits[i,*] = cgpercentiles(samples[i,*],percentiles = [dc,1d - dc])
    ;print, strcompress(limits[i,0], /remove_all)+' < parameter', strcompress(i, /remove_all)+' < '+strcompress(limits[i,1], /remove_all)
  endfor

  result =list(call_function(model_funct[0],x[0],pars))
  for  i = 1, n_elements(model_funct) -1 do begin
    result.add, call_function(model_funct[i],x[i],pars)   
  endfor
  return, result
end