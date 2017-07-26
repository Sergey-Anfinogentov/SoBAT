   ;+
   ; :Description:
   ;    Simulates random values from the multivariative normal distribution
   ;
   ; :Params:
   ;    seed - seed for the random number generator
   ;    mu - expect value (1D array)
   ;    sigma - covariance matrix (2D array)
   ;    n - number of values to generate
   ;
   ;
   ;
   ; :Author:  Sergey Anfinogentov (sergey.istp@gmail.com)
   ;-
function mcmc_random_multyn,seed,mu,sigma,n
   n_sigma = (size(sigma))[1]
   a = sigma
   if n_elements(a) eq 1 then return, sqrt(Sigma) * randomn(seed,n) + mu
   ;if total(a ne 0) eq 0 then return, mu*0d
   LA_CHOLDC, a
   for i = 0,n_sigma - 2 Do a[i+1:*,i] = 0
   
   result = dblarr(n_sigma,n)
   for i =0, n-1 do begin
    z = randomn(seed,n_sigma)
    result[*,i] = mu + a##z
   endfor
   return,result
   
end