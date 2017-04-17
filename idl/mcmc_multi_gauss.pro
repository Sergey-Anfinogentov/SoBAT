function mcmc_multi_gauss, x, mu, sigma
  sz = size(x)
  n = sz[1]
  n_par = n_elements(mu)
  if n_par eq 1 then return,mu + 1d/(2d*!dpi*sigma)*exp( -(x - mu)^2/(2d*sigma))
  if sz[0] gt 1 then begin
    np = sz[2]
    result = dblarr(np)
    for i =0, np -1 do begin
      result[i] = multi_gauss(reform(x[*,i]),mu,sigma)
    endfor
  endif
  det_sigma = la_determ(sigma, zero = 0d)
  inv_sigma = la_invert(sigma)
  x_mu = x - mu
  ex = transpose(x_mu)#inv_sigma#x_mu
  ex = -0.5d*ex[0]
  const = 1d/((2d*!dpi)^(n/2d)*sqrt(abs(det_sigma)))
  return,const*exp(ex)
end