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
    max_acceptance_rate: 0.7d,$; maximal acceptance rate
    proposal_tune_samples: 1000l,$; number of samples used to tune the proposal distribution 
    double_ppd_samples: 0b,$; save ppd_samples in double precision (default 0b), computations are always done in double precision
    no_ppd_samples: 1b,$; do not save ppd_samples
    version: 'v 0.3.0'$
   }
     
   return, settings
end