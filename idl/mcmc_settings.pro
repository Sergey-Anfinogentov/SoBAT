;+
;This file contains settings for the IDL MCMC library
; :Author: sergey
;-
function mcmc_settings
compile_opt idl2
  settings = {$
    sigma_buffer_size :300,$;Size of the buffer used to tune the sampling parameters
    target_acceptance_rate: 0.23d,$;Target acceptanced for the sampler
    printing_interval: 2d$;interval for printing informational messages in seconds
   }
     
   return, settings
end