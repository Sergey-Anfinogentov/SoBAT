;+
; :Description:
;    Calculates 2D histogram from the posterior predictive distribution samples
;
; :Params:
;    x - X-value of the data points
;    y - Y-value of the data points
;    ppd_samples
;
; :Keywords:
;    nbins - (in) number of bins used to make a PPD histogram (default:256)
;    hist_x - (out) X-value for each bin centre
;    hist_y - (out) Y-value for each bin centre
;
; :Author: Sergey Anfinogentov (sergey.istp@gmail.com)
;-
function mcmc_ppd_histogram, x, y, ppd_samples, nbins = nbins, hist_x = hist_x, hist_y = hist_y
  if not keyword_set(nbins) then nbins = 256
  nx = n_elements(x)
  n_samples = double(n_elements(ppd_samples[0,*]))

  result = dblarr(nx, nbins)

  yrange = minmax(y)
  yrange += (yrange[1]-yrange[0])*[-0.25,0.25]

  for i =0, nx-1 do begin
    result[i,*] = histogram(ppd_samples[i,*],nbins = nbins, max =yrange[1], min = yrange[0], loc = loc)
  endfor
  hist_x = x
  hist_y = loc+ loc[1]-loc[0]
  dy = loc[1] - loc[0]
  result = result/(dy*n_samples)
  return, result
end