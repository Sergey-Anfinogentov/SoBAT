function prior_exponential::get_log_value,x
  return, mcmc_exponential(x, self.lambd)
end
function prior_exponential::get_start_value
  return, 1d/self.lambd
end
function prior_exponential::init, lambd
  self.lambd = lambd
  return,1
end

pro prior_exponential__define
  structure = {prior_exponential,lambd:0d}
end