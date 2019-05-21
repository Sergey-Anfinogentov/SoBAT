function prior_half_normal::get_log_value,x
  return, mcmc_half_normal(x, self.sigma)
end
function prior_half_normal::get_start_value
  return, self.sigma*sqrt(2d)/sqrt(!dpi)
end
function prior_half_normal::init, sigma
  self.sigma = sigma
  return,1
end

pro prior_half_normal__define
  structure = {prior_half_normal, sigma:1d}
end