function mcmc_exponential, x, lambd
  if total(x lt 0) gt 0 then return, -!values.D_INFINITY
  res = alog(lambd) - lambd * x
  return, total(res)
end