


function lin_model,x,pars,_extra = _extra
  k = pars[0]
  b = pars[1]
  return, k * x + b
end

function quad_model,x,pars,_extra=_extra
  k0 = pars[0]
  k1 = pars[1]
  k2 = pars[2]
  return, k1 * x + k0 + k2*x^2
end


function mcmc_fit_test::test_linear
  
  x = findgen(100)*0.1
  k = 0.5d
  b = 1.d
  sigma = 2.d
  seed = 100 ;random_seed()
  y = k * x + b + sigma*randomn(seed,100)


  ;Setting parameter limits
  limits = dblarr(2,2)
  limits[0,*] = [-5d,5d]
  limits[1,*] = [-5d,5d]

  limits0 = limits
  ;define the initial guess
  pars = [1d, 1d]
  ;define the number of samples
  n_samples = 100000
  ;define the number of burn in samples
  burn_in = 10000
  fit = mcmc_fit(x, y, pars, 'lin_model', limits  = limits,  burn_in =burn_in, n_samples = n_samples, samples = samples, credible_intervals=credible_intervals,/silent)
  
  assert, k ge credible_intervals[0,0] and k le credible_intervals[0,1], 'True value of k parameter lies outside credible intervals: ', k
  assert, b ge credible_intervals[1,0] and b le credible_intervals[1,1], 'True value of b parameter lies outside credible intervals: ', b
  

  return, 1
end

function mcmc_fit_test::test_linear_priors

  x = findgen(100)*0.1
  k = 0.5d
  b = 1.d
  sigma = 2.d
  seed = 100 ;random_seed()
  y = k * x + b + sigma*randomn(seed,100)
  
  

  ;Setting priors
  priors = [$
    prior_uniform(-5d, 5d), $; k
    prior_uniform(-5d, 5d)  $; b
  ]

  ;define the initial guess
  pars = [1d, 1d]
  ;define the number of samples
  n_samples = 100000
  ;define the number of burn in samples
  burn_in = 10000
  
  fit = mcmc_fit(x, y, pars,  'lin_model', priors = priors,  burn_in =burn_in, n_samples = n_samples, samples = samples, credible_intervals=credible_intervals,/silent)

  assert, k ge credible_intervals[0,0] and k le credible_intervals[0,1], 'True value of k parameter lies outside credible intervals: ', k
  assert, b ge credible_intervals[1,0] and b le credible_intervals[1,1], 'True value of b parameter lies outside credible intervals: ', b


  return, 1
end

function mcmc_fit_test::test_quad

  x = findgen(100)*0.1
  k1 = 0.5d
  k2 = 0.1
  k0 = 1.d
  sigma = 2.d
  seed = 100 ;random_seed()
  y =k2*x^2+ k1 * x + k0 + sigma*randomn(seed,100)


  ;Setting parameter limits
  limits = dblarr(3,2)
  limits[0,*] = [-5d,5d]
  limits[1,*] = [-5d,5d]
  limits[2,*] = [-5d,5d]

  limits0 = limits
  ;define the initial guess
  pars = [1d, 1d, 0d]
  ;define the number of samples
  n_samples = 100000
  ;define the number of burn in samples
  burn_in = 10000
  fit = mcmc_fit(x, y, pars,  'quad_model', limits  = limits, burn_in =burn_in, n_samples = n_samples, samples = samples, credible_intervals=credible_intervals,/silent)

  assert, k0 ge credible_intervals[0,0] and k0 le credible_intervals[0,1], 'True value of k0 parameter lies outside credible intervals: ', k0
  assert, k1 ge credible_intervals[1,0] and k1 le credible_intervals[1,1], 'True value of k1 parameter lies outside credible intervals: ', k1
  assert, k2 ge credible_intervals[2,0] and k2 le credible_intervals[2,1], 'True value of k2 parameter lies outside credible intervals: ', k2



  return, 1
end

pro mcmc_fit_test__define
  compile_opt strictarr

  define = {mcmc_fit_test, inherits MGutTestCase }
end