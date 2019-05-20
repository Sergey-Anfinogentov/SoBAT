function mcmc_prior_test::test_uniform
;  stop
  u = prior_uniform(0d, 1d)
  res = u.get_log_value(0.5d)
  assert, res eq 0., 'Wrong output of a uniform prior (0d, 1d) at 0.5: %f instead of 0.', res
  
  u = prior_uniform(0d, 1d)
  res = u.get_log_value(1.5d)
  assert, res eq -!values.D_INFINITY, 'Wrong output  of a uniform prior (0d, 1d) at 1.5: %f instead of -infinity', res
  
  u = prior_uniform(0d, 1d)
  res = u.get_log_value(-0.2d)
  assert, res eq -!values.D_INFINITY, 'Wrong output  of a uniform prior (0d, 1d) at -0.2: %f instead of -infinity', res
  
  u = prior_uniform(0d, 1d)
  res = u.get_log_value([-0.2d, 0.5])
  assert, res eq -!values.D_INFINITY, 'Wrong output  of a uniform prior (0d, 1d) at [-0.2d, 0.5]: %f instead of -infinity', res
  
  u = prior_uniform(0d, 1d)
  res = u.get_log_value([0.2d, 0.5])
  assert, res eq 0d, 'Wrong output  of a uniform prior (0d, 1d) at [0.2d, 0.5]: %f instead of 0', res
  
  res = mcmc_uniform(0.5d, 0d, 5d)
  u = prior_uniform(0d, 5d)
  res = u.get_log_value(0.5)
  assert, res eq  alog(0.2d), 'Wrong output  of a uniform prior (0d, 5d) at 0.5: %f instead of %f', res, alog(0.2)
  
  u = prior_uniform(0d, 5d)
  res = u.get_log_value([0.5d,2.,3.])
  assert, abs(res - alog(0.2d^3)) lt 1d-10, 'Wrong output  of a uniform prior (0d, 5d) at [0.5d,2.,3.]: %f instead of %f', res, alog(0.2)


  return,1
end

function mcmc_prior_test::test_normal
  ;  stop
  small_number = 1d-5
  
  n = prior_normal(0d, 1.)
  res = n.get_log_value(0.0)
  expected = alog(1d/sqrt(2d*!dpi))
  assert, abs(res - expected) lt small_number , 'Wrong output of a normal prior (0d, 1d) at 0.0: %f instead of %f', res, expected
  
  n = prior_normal(0d, 2d)
  res = n.get_log_value(0.5)
  expected = alog(0.193334d)
  assert, abs(res - expected) lt small_number , 'Wrong output of a normal prior (0d, 2d) at 0.5: %f instead of %f', res, expected
  
  n = prior_normal(2d, 1d)
  res = n.get_log_value(0.5)
  expected = alog(0.129518d)
  assert, abs(res - expected) lt small_number , 'Wrong output of a normal prior (2d, 1d) at 0.5: %f instead of %f', res, expected
  
  n = prior_normal(1.5d, 8d)
  res = n.get_log_value(-2d)
  expected = alog(0.0453165d)
  assert, abs(res - expected) lt small_number , 'Wrong output of a normal prior (1.5d, 8d) at -2: %f instead of %f', res, expected


  return,1
end

function mcmc_prior_test::test_half_normal
  ;  stop
  small_number = 1d-5

  prior = prior_half_normal(1d)
  res = prior.get_log_value(0.0)
  expected = alog(2d/sqrt(2d*!dpi))
  assert, abs(res - expected) lt small_number , 'Wrong output of a half_normal prior (sigma = 1d) at 0.0: %f instead of %f', res, expected

  prior = prior_half_normal( 2d)
  res = prior.get_log_value(0.5)
  expected = alog(2d*0.193334d)
  assert, abs(res - expected) lt small_number , 'Wrong output of a half_normal prior (sigma = 2d) at 0.5: %f instead of %f', res, expected

  prior = prior_half_normal(1d)
  res = prior.get_log_value(0.5)
  expected = alog(2d*0.352065d)
  assert, abs(res - expected) lt small_number , 'Wrong output of a half_normal prior (sigma = 1d) at 0.5: %f instead of %f', res, expected

  prior = prior_half_normal(8d)
  res = prior.get_log_value(-2d)
  expected = alog(2d*0.0453165d)
  assert,res eq  -!values.D_INFINITY, 'Wrong output of a half_normal prior (sigma = 8d) at -2: %f instead of %f', res, -!values.D_INFINITY


  return,1
end



pro mcmc_prior_test__define
  compile_opt strictarr

  define = {mcmc_prior_test, inherits MGutTestCase }
end