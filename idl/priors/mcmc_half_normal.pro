function mcmc_half_normal, x, sigma
  
  if total(x lt 0) gt 0 then return, -!values.D_INFINITY
  a = alog(1d/sqrt(0.5 * !dpi * sigma^2))
  exponen = -x^2/(2*sigma^2)
  return, total(exponen + a)
end