function prior_uniform::get_log_value,x
  return, mcmc_uniform(x, self.lower, self.upper)
end
function prior_uniform::get_start_value
  return, (self.lower + self.upper)*0.5
end
function prior_uniform::init, lower, upper
  self.lower = lower
  self.upper = upper
  return,1
end

pro prior_uniform__define
  structure = {prior_uniform,lower:0d, upper:1d}
end