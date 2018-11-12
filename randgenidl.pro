function randgenidl, seed1, seed2, seed3, n
t0=systime(1)
;+
; NAME:
;       RANDGENIDL
; PURPOSE:
;       Generate pseudo random numbers uniformly distributed in [0,1]
; EXPLANATION:
;       This is Wichmann-Hill algorithm which has a period ~ 7*10^12  
;
; CALLING SEQUENCE:
;       Result = RANDGENIDL(seed1, seed2, seed3, n)
; 
; INPUTS:
;       seed1, seed2, seed3 = named variables containing a non-zero integer 
;       n = integer being the number of desired pseudo-random numbers
;
; OUTPUTS:
;       result = vector of n pseudo-random numbers
;
; EXAMPLE:
;       IDL> seed1=1
;       IDL> seed2=2
;       IDL> seed3=3
;       IDL> print,randgen(seed1,seed2,seed3,3)
;       0.0338188     0.777542    0.0527352
;       
; LIMITATIONS:
;       The generator has a cycle length ~ 7*10^12
;       0 should NOT be used as a seed, as it remains the same.
;       The routine can generate ~ 500,000 numbers per second, while
;       RANDOMU generates ~ 30,000,000 numbers per second
;
; MODIFICATION HISTORY:
;       Written, F. Lacasa, IAS, October 2011
 
result=fltarr(n)
for i=0L,n-1 do begin
seed1 = 171*(seed1 mod 177) - 2*(seed1/177)
if seed1 lt 0 then seed1+=30269
seed2 = 172*(seed2 mod 176) - 35*(seed2/176)
if seed2 lt 0 then seed2+=30307
seed3 = 170*(seed3 mod 178) - 63*(seed3/178)
if seed3 lt 0 then seed3+=30323
temp = seed1/30269. + seed2/30307. + seed3/30323.
result[i] = temp-floor(temp)
endfor
print,systime(1)-t0
return, result
end
