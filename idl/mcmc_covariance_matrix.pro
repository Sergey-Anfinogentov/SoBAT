function mcmc_covariance_matrix,samples, mu = mu
  sz = size(samples)
  n = sz[1]
  sigma = dblarr(n,n)
  mu = dblarr(n)
  for i =0, n-1 do begin
    mu[i] = mean(samples[i,*])
    for j =0, n-1 do begin
      sigma[i,j] = correlate(samples[i,*],samples[j,*],/covariance)    
    endfor
  endfor
  return,sigma
end