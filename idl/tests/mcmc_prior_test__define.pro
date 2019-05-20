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
  
  res = mcmc_uniform(0.5d, 0d, 5d)
  assert, res eq alog(0.2d), 'Wrong output of mcmc_uniform(0.5d, 0d, 1d): %f instead of %f', res, alog(0.2)
  res = mcmc_uniform([0.5d,2.,3.], 0d, 5d)
  assert, abs(res - alog(0.2d^3)) lt 1d-10, 'Wrong output of mcmc_uniform([0.5d,2.,3.], 0d, 5d): %f instead of %f', res, alog(0.2d^3)

  return,1
end

function mcmc_prior_test::test_normal
  ;  stop
  small_number = 1d-5
  
  res = mcmc_normal(0.0d, 0d, 1d)
  expected = alog(1d/sqrt(2d*!dpi))
  assert, abs(res - expected) lt small_number , 'Wrong output of mcmc_normal(0.0d, 0d, 1d): %f instead of %f', res, expected
  
  res = mcmc_normal(0.5d, 0d, 2d)
  expected = alog(0.193334d)
  assert, abs(res - expected) lt small_number , 'Wrong output of mcmc_normal(0.0d, 0d, 1d): %f instead of %f', res, expected
  
  res = mcmc_normal(0.5d, 2d, 1d)
  expected = alog(0.129518d)
  assert, abs(res - expected) lt small_number , 'Wrong output of mcmc_normal(0.0d, 0d, 1d): %f instead of %f', res, expected
  
  res = mcmc_normal(-2d, 1.5, 8d)
  expected = alog(0.0453165d)
  assert, abs(res - expected) lt small_number , 'Wrong output of mcmc_normal(0.0d, 0d, 1d): %f instead of %f', res, expected


  return,1
end


pro mcmc_prior_test__define
  compile_opt strictarr

  define = {mcmc_prior_test, inherits MGutTestCase }
end