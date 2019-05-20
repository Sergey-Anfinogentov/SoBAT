function mcmc_normal, x, expectation, sigma
  a = alog(1d/sqrt(2d * !dpi * sigma^2))
  exponen = - (x - expectation)^2/(2*sigma^2)
  return, total(exponen + a)
end