
;+
; :Description:
;    Uses Byesian Inference to fit a set of user suplied functions to the suplied points (X, Y) by adjasting a set of parameters (PARS)
;
; :Params:
;    x - a list of independent variables, where each element is an array of X-s corresponding to one of the fitted functiom
;    y - a list of measurments of the dependent variable, , where each element is an array of Y-s corresponding to one of the fitted functiom
;    pars - dblarr(n_params), Starting guess, will containt the fitted parameter
;           If the starting guess is not set a random starting point will be generated
;    limits - dblarr(n_params, 2) possible limitspars for the parameters
;    model_funct - an array of the names of model functions to fit. The model function must accept 2 parameters, X and PARS.
;
; :Keywords:
;    n_samples - numer of samples to generate using Metropolis-Hastings MCMC, default 10000l
;    burn_in - number of burn in samples, needed for the sampler to find the high probability region
;             and to tune sampling parameters
;    samples - dblarr(n_params, n_samples) will contain samples from the posteriour PDF
;    ppd_samples - list(dblarr(n_data_points, n_samples),..) will contain samples from the Posteriour Predictive Distribution.
;                   For multifunction fitting, PPD samples will be returned as a list.
;    confidence_level - confidence level to define confidence intervals for each parameter.
;    credible_intervals - (output) will contain credible intervals for each paramet
;    sigma_samples - samples of the standart deviasion of the observational noise which is assumed to be  normally distributed
;
; :Author: Sergey Anfinogentov (sergey.istp@gmail.com)
;-
function mcmc_fit_multi,x,y,pars, limits ,model_funct,n_samples = n_samples, sigma_samples = sigma_samples, burn_in = burn_in,$
  samples = samples, confidence_level = confidence_level, credible_intervals=credible_intervals,$
  noise_limits = noise_limits, ppd_samples =ppd_samples,  _extra = _extra
  compile_opt idl2
  if not keyword_set(n_samples) then n_samples = 10000l
  if not keyword_set(burn_in) then burn_in = 5000l
  if not keyword_set(confidence_level) then confidence_level = 0.95d
  if not keyword_set(pars) then pars = mcmc_fit_random_start(limits)
  n_par = n_elements(pars)
  n_funct = n_elements(model_funct)
 ; pars_ = [pars,replicate(1d,n_funct)]


  
  
  limits_ = dblarr(n_par+n_funct,2)
  limits_[0:n_par-1,*] = limits
  
  ;defining limits for noises
   if not keyword_set(noise_limits) then begin
    noise_limits= dblarr(n_funct,2)
    for i = 0, n_funct -1 do begin
      noise_limits[i,*] = [0, max(y[i]) - min(y[i])]
    endfor
   endif
  limits_[n_par:*,*] = noise_limits
  
  ;calculating initial guesses for the noise
   noise_guesses = dblarr(n_funct)
  for i = 0, n_funct -1 do begin
      y_guess = call_function(model_funct[i], x[i], pars)
      noise_guesses[i] = stddev(y[i]-y_guess)<noise_limits[i,1]>noise_limits[i,0]
  endfor
  pars_ = [pars, noise_guesses]


  sigma = (max(limits_,dim = 2) - min(limits_,dim = 2))/6d
  ;sigma = [sigma,(max(y) - min(y))*0.01d]

  samples = mcmc_sample(pars_,'mcmc_fit_ln_prob_multi',n_samples, burn_in =  burn_in, x = x, y = y, model_funct = model_funct,$
            limits = limits_, sigma = sigma, ppd_samples =ppd_samples,  values = values)
  ; samples = samples[*, burn_in:*]
  sigma_samples = samples[n_par:*,*]
  
  n_ppd = lonarr(n_funct)
  for i = 0, n_funct-1 do begin
    n_ppd[i] = n_elements(x[i])
  endfor
  ppd_ind =total([0,n_ppd],/cum)
  ppd_samples_old =  temporary(ppd_samples)
  ppd_samples = list()
  for i = 0, n_funct -1 do begin
    x_i = x[i]
    y_i = y[i]
    ppd_samples.add,ppd_samples_old[ppd_ind[i]:ppd_ind[i+1]-1,*]
  endfor


  ;samples = samples[0:n_par-1,*])
  foo = max(values, ind)
  pars = samples[*,ind]
  pars = pars[0:n_par - 1]

  dc = (1d - confidence_level)*0.5d
  credible_intervals = limits
  for i =0, n_par -1 do begin
    credible_intervals [i,*] = cgpercentiles(samples[i,*],percentiles = [dc,1d - dc])
    ;print, strcompress(limits[i,0], /remove_all)+' < parameter', strcompress(i, /remove_all)+' < '+strcompress(limits[i,1], /remove_all)
  endfor

  result =list(call_function(model_funct[0],x[0],pars))
  for  i = 1, n_elements(model_funct) -1 do begin
    result.add, call_function(model_funct[i],x[i],pars)   
  endfor
  return, result
end