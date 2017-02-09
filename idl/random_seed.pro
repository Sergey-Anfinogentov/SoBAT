FUNCTION   random_seed
return,LONG((SYSTIME(1) - long(systime(1)/1d4)*1d4 )*1D4)
end