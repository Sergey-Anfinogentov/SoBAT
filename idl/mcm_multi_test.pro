function model_fun1, x, a

  return, a[1]*x + a[0]

end

function model_fun2, x, a

  return, -1.5d*a[1]*x + a[0]

end

pro generate_model,x,y

  x1 = dindgen(100)
  x2 = dindgen(100) + 150d
  a = [1., 0.2]
  
  sigma1 = 0.5
  sigma2 = 1.6
  
  s = random_seed()
  
  y1 = model_fun1(x1,a) + randomn(s,100) * sigma1
  y2 = model_fun2(x2,a) + randomn(s,100) * sigma2
  
  x = list(x1, x2)
  y = list(y1, y2)


end

pro mcm_multi_test
  generate_model,x,y
  window,0
  plot, x[0], y[0], psym = 1
  window,1
  plot, x[1], y[1], psym = 1
  model_fun = ['model_fun1','model_fun2']
  

  lim = dblarr(2,2)
  lim[0,*]  = [-10d,10d]
  lim[1,*] = [-10d,10d]
  
  limi = lim
  
  fit = mcmc_fit(x,y,pars, lim, model_fun, sigma_samples = sigma_samples, n_samples = 10d3, samples = samples)
  print, pars

  
  
  evidence=mcmc_fit_evidence(samples,x,y,limi, model_fun,n_iterations=10d3)
  help, evidence
  limits = mcmc_fit_estimate_y_limits(x, samples, model_fun, sigma_samples = sigma_samples, confidence_level=0.99d)
  
  wset,0
  plot, x[0], y[0], psym = 1
  oplot, x[0], fit[0], color = 155
  oplot, x[0],limits[0,*,0], color = 100
  oplot, x[0],limits[0,*,1], color = 100
  
  wset,1
  plot, x[1], y[1], psym = 1
  oplot, x[1], fit[1], color =155
  oplot, x[1],limits[1,*,1], color = 100
  oplot, x[1],limits[1,*,0], color = 100
;  stop

end