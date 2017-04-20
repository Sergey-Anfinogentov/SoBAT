;+
;This file contains settings for the IDL MCMC library
; :Author: sergey
;-
function mcmc_settings
compile_opt idl2
  settings = {$
    acceptance_buffer_size :1000l,$;Size of the buffer used to calculate the acceptance rate
    printing_interval: 2d$;interval for printing informational messages in seconds
   }
     
   return, settings
end