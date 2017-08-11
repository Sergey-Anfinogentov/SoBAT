;+
;This file contains settings for the IDL MCMC library
; :Author:  Sergey Anfinogentov (sergey.istp@gmail.com)
;-
function mcmc_settings
compile_opt idl2
  settings = {$
    acceptance_buffer_size :2000l,$;Size of the buffer used to calculate the acceptance rate
    printing_interval: 2d,$;interval for printing informational messages in seconds
    min_acceptance_rate: 0.1d,$; Minimal acceptance rate
    max_acceptance_rate: 0.5d,$; maximal acceptance rate
    proposal_tune_samples: 1000l,$; number of samples used to tune the proposal distribution
    version: 'v 0.2.1'$
   }
     
   return, settings
end