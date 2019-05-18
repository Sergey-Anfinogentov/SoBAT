
function mcmc_uniform, x, lower, upper
  n = n_elements(x)
  volume = (upper - lower)^n
  if total((x lt lower) or (x gt upper)) gt 0 then return, -!values.D_INFINITY
  return, -alog(volume)
end