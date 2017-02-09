function mcmc_random_multyn,seed,mu,sigma,n
   n_sigma = (size(sigma))[1]
   a = sigma
   LA_CHOLDC, a
   for i = 0,n_sigma - 2 Do a[i+1:*,i] = 0
   
   result = dblarr(n_sigma,n)
   for i =0, n-1 do begin
    z = randomn(seed,n_sigma)
    result[*,i] = mu + a##z
   endfor
   return,result
   
end