;+
; :Description:
;    Prints informational messages. Maximum 1 message per printing_interval can be printed.  
;
; :Params:
;    text - message to print
;    printing_interval - time interval between subsequent messages in seconds.
;                          If this parameter is not suplied printing interval
;                          defined in mcmc_settings.pro will be used
;
;
;
; :Author: Sergey Anfinogentov
;-
pro mcmc_message, text, printing_interval
  common mcmc_message, time
  if not keyword_set(time) then time = 1d
  if not keyword_set(printing_intervaal) then begin
    settings = mcmc_settings()
    printing_interval = settings.printing_interval
  endif
  if systime(1) - time gt printing_interval then begin
    message, text, level = -1 ,/info
    time = systime(1)
  endif

end