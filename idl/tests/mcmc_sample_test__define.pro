function histogram_2d, x, y, x_nbins=x_nbins, y_nbins=y_nbins, x_locations=x_locations, y_locations=y_locations, center = center
  if not keyword_set(x_nbins) then x_nbins =20
  if not keyword_set(y_nbins) then y_nbins =20
  max_x = max(x)
  min_x = min(x)
  max_y = max(y)
  min_y = min(y)
  bin_x = (max_x - min_x) / x_nbins
  bin_y = (max_y - min_y) / y_nbins
  x_locations = dindgen(x_nbins) * bin_X + min_x
  y_locations = dindgen(y_nbins) * bin_y + min_y
  if keyword_set(center) then begin
    x_locations += bin_x * 0.5
    y_locations += bin_y * 0.5
  endif
  x_norm = (x - min(x))/(max_x-min_x)*x_nbins
  y_norm = (y - min(y))/(max_y-min_y)*y_nbins
  result = lonarr(x_nbins,y_nbins)
  x_ind = floor(x_norm)<(x_nbins-1)
  y_ind = floor(y_norm)<(y_nbins-1)
  for i = 0, n_elements(x_ind)-1 do begin
    result[x_ind[i], y_ind[i]] += 1l
  endfor
  return, result
end
function uniform, par, _extra = _extra
  x = par[0]
  a = 0.5d
  b = 3d

  if (x gt a) and (x lt b) then return, alog(1d/(b -a))
  return, -!values.d_infinity
end

function exponential, par, _extra = _extra
  x = par[0]
  a = 1d


  if (x gt 0)  then return, alog(a*exp(-a*x))
  return, -!values.d_infinity
end

function triangular, par, _extra = _extra
  x = par[0]
  a = 0.5d
  b = 3d
  c = 2d

  if (x gt a) and (x lt c) then return, alog(2d*(x-a)/(b-a)/(c-a))
  if (x gt c) and (x lt b) then return, alog(2d*(b-x)/(b-a)/(b-c))
  return, -!values.d_infinity
end

function bimodal, par, _extra = _extra
  x = par[0]
  sigm1 = 2d
  sigm2 = 1d
  M1 = 0d
  M2 = 7d
  f1 = 1d/sqrt(2*sigm1^2*!dpi)*exp(-(x - m1)^2/(2d*sigm1^2))
  f2 = 1d/sqrt(2*sigm2^2*!dpi)*exp(-(x - m2)^2/(2d*sigm2^2))
  return, alog(f1*0.8d +f2*0.2d)
end

function pyramid, par, _extra = _extra
  x0 =[1d,2d,2d,3d]
  y0 =[1d,2d,3d,1.5d]
  z0 =[0d,1d,0d,0d]

  x = par[0]
  y = par[1]
  triangulate, x0, y0, tr
  res = griddata(x0,y0,z0,/linear, triangles = tr, xout = [x], yout = [y], missing = 0d)
  return,alog(res[0]>0d)
end


function uniform_2d, par, _extra =_extra
  x1 = 1d
  x2 = 3d
  y1 = 1d
  y2 = 3d

  x = par[0]
  y = par[1]
  if (x gt x1) and (x lt x2) and (y gt y1) and (y lt y2) then return, alog(1d/((y2-y1)*(x2-x1)))
  return, -!values.d_infinity

end

function exponential_2d, par, _extra = _extra
  x = par[0]
  y = par[1]
  a1 = 1d
  a2 = 2d


  if (x gt 0) and (y gt 0)  then return, alog(a1*exp(-a1*x)*a2*exp(-a2*y))
  return, -!values.d_infinity
end

function multimodal_2d, par, _extra = _extra

  sigma1 = [[1d,0.8d],[0.8d,1d]]*0.3
  sigma2 = [[1d,0d],[0d,1d]]*0.2
  sigma3 = [[1d,-0.6d],[-0.6d,1d]]*0.2
  mu1 = [2,2]
  mu2 = [3,1]
  mu3 = [0,3]

  v1 = mcmc_multi_gauss(par, mu1, sigma1)
  v2 = mcmc_multi_gauss(par, mu2, sigma2)
  v3 = mcmc_multi_gauss(par, mu3, sigma3)

  return, alog(v1 + v2 + v3*0.5)
end



pro mcmc_sample_test::check_mcmc_1d, fun_name

  samples =  mcmc_sample([1d],fun_name,100000, burn_in = 50000,/silent)
  nbins = 20
  h = double(histogram(samples, loc = loc, nbins = nbins))
  loc += (loc[1]-loc[0])*0.5d
  dloc = loc[1]-loc[0]
  norm = total(h)*dloc
  h /=norm
  v = h
  for i = 0, n_elements(h)-1 do v[i]=exp(call_function(fun_name,loc[i]))

  chi2 = max((v-h)^2/max(v)^2)

  assert, chi2 le 0.005, 'Apporoximation metric is above 0.005: %f', chi2

end

pro mcmc_sample_test::check_mcmc_2d, fun_name
  tol =1d-2
  nx = 20
  ny = 20
  n_samples = 5d5
  samples =  mcmc_sample([2d, 2d],fun_name,n_samples, burn_in = 50000,/silent)
  h = histogram_2d(samples[0,*], samples[1,*], x_loc = x_loc, y_loc = y_loc, x_nbins = nx, y_nbins = ny,/center)
  expected = dblarr(nx, ny)
  dx = x_loc[1]-x_loc[0]
  dy = y_loc[1]-y_loc[0]
  norm = total(h)*dx*dy
  h /=norm
  for ix =0, nx -1 do begin
    for iy = 0, ny - 1 do begin
      expected[ix,iy] = exp(call_function(fun_name,[x_loc[ix],y_loc[iy]]))
    endfor
  endfor
  expected /= total(expected)*dx*dy

  chi2 = max((expected-h)^2/total(expected^2))

  assert, chi2 le tol, 'Apporoximation metric is above the tollerance level: %f', chi2

end


function mcmc_sample_test::test_1D_uniform

  self->check_mcmc_1d, 'uniform'

  return, 1
end

function mcmc_sample_test::test_1D_exponential

  self->check_mcmc_1d, 'exponential'

  return, 1
end

function mcmc_sample_test::test_1D_triangular

  self->check_mcmc_1d, 'triangular'

  return, 1
end

function mcmc_sample_test::test_1D_bimodal

  self->check_mcmc_1d, 'bimodal'

  return, 1
end

function mcmc_sample_test::test_2D_pyramid

  self->check_mcmc_2d, 'pyramid'

  return, 1
end

function mcmc_sample_test::test_2D_uniform

  self->check_mcmc_2d, 'uniform_2d'

  return, 1
end

function mcmc_sample_test::test_2D_exponential

  self->check_mcmc_2d, 'exponential_2d'

  return, 1
end

function mcmc_sample_test::test_2D_multimodal

  self->check_mcmc_2d, 'multimodal_2d'

  return, 1
end

pro mcmc_sample_test__define
  compile_opt strictarr

  define = {mcmc_sample_test, inherits MGutTestCase }
end