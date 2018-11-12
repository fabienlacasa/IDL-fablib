Function spherbess,n,X,hierarchy=hierarchy
;+
; NAME:
;       SPHERBESS
;
; PURPOSE:
;       Compute the spherical bessel function of order n. Or all
;       orders up to n if keyword /hierarchy is present
; 
; EXPLANATION:
;       Use the recursion relation, and cuts at low x where the
;       computation blows up because of round-off errors.
;
; CALLING SEQUENCE:
;       Result = spherbess(n,X [,/hierarchy])
;
; INPUT:
;       n = integer containing the bessel order
;       x = scalar or vector
;
; KEYWORD:
;       /hierarchy : returns an array[n_elements(x),n+1] containing
;       the spherical bessel function from order 0 to n
;
; OUTPUT:
;       result = j_n(x)
;         or
;       result = [[j_0(x)],[j_1(x)] ... [j_n(x)]]  if /hierarchy is present
;
; EXAMPLE:
;       IDL> print,spherbess(2,0.1)  ===> 0.00066619063
;
;       In this example the exact value ca be computed analytically as 
;       j2(x)=(3/x^2-1)*sin(x)/x-3*cos(x)/x^2
;       IDL returns 0.000671387 in single precision and 0.00066619063
;       in double precision
;
; RESTRICTIONS:
;       Cuts the function to zero at the very beginning of the first
;       peak where the recursion would otherwise blow up ; value of
;       the function at this place would normally be less than 10^-7.
;       Tested up to n=10000
;
; MODIFICATION HISTORY:
;       Written, F. Lacasa, IAS, June 2011
;

on_error,2

;ensure double precision
wh0=where(x eq 0,count0) ;avoid problems with x=0
xx=double(x)
if count0 ne 0 then xx(wh0)=1D

j0 = sin(xx)/xx
j1 = sin(xx)/xx^2 - cos(xx)/xx
if count0 ne 0 then begin
j0(wh0)=1D
j1(wh0)=0D
endif

spherbesstable=dblarr(n_elements(xx),n+1)
spherbesstable(*,0) = j0
if n ge 1 then spherbesstable(*,1)=j1

for i=2,n do begin
spherbesstable(*,i)=(2D*i-1D)*spherbesstable(*,i-1)/xx-spherbesstable(*,i-2)
if count0 ne 0 then spherbesstable(wh0,i)=0D
xblow=i-33.*(i/200.)^0.33 ;computation blows up before xblow
if i le 75 then xblow=0.75*i-5.
if i le 7 then xblow=0.7*(i/8.)^3
whb=where(xx le xblow,countb)
if countb ne 0 then spherbesstable(whb,i)=dblarr(n_elements(whb))
endfor

if ~keyword_set(hierarchy) then result=spherbesstable(*,n) $
else result=spherbesstable

return, result
end
