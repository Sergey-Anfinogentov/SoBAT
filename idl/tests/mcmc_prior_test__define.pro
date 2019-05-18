function mcmc_prior_test::test_uniform
;  stop
  res = mcmc_uniform(0.5d, 0d, 1d)
  assert, res eq 0., 'Wrong output of mcmc_uniform(0.5d, 0d, 1d): %f instead of 0.', res
  
  res = mcmc_uniform(1.5d, 0d, 1d)
  assert, res eq -!values.D_INFINITY, 'Wrong output of  mcmc_uniform(1.5d, 0d, 1d): %f instead of -infinity', res
  
  res = mcmc_uniform(-0.2d, 0d, 1d)
  assert, res eq -!values.D_INFINITY, 'Wrong output of  mcmc_uniform(-0.2d, 0d, 1d): %f instead of -infinity', res
  
  res = mcmc_uniform([-0.2d, 0.5], 0d, 1d)
  assert,  res eq -!values.D_INFINITY, 'Wrong output of mcmc_uniform([-0.2d, 0.5], 0d, 1d): %f instead of -infinity', res
  
  res = mcmc_uniform([0.2d, 0.5], 0d, 1d)
  assert,  res eq 0., 'Wrong output of mcmc_uniform([0.2d, 0.5], 0d, 1d): %f instead of 0.',  res

  return,1
end

pro mcmc_prior_test__define
  compile_opt strictarr

  define = {mcmc_prior_test, inherits MGutTestCase }
end