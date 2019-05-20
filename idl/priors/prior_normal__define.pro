function prior_normal::get_log_value,x
  return, mcmc_normal(x, self.expectation, self.sigma)
end
function prior_normal::init, expectation, sigma
  self.expectation = expectation
  self.sigma = sigma
  return,1
end

pro prior_normal__define
  structure = {prior_normal,expectation:0d, sigma:1d}
end