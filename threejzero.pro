function ThreeJzero, l1,l2,l3 ;computes 3J(l1,l2,l3)
;+
; NAME:
;       ThreeJzero
; PURPOSE:
;       The 3J symbol (l1,l2,l3 ; 0,0,0) is computed through the formula
;       (-1)^g*sqrt[((2g-2l1)!*(2g-2l2)!*(2g-2l3)!)/(2g+1)!]*g!/((g-l1)!*(g-l2)!*(g-l3)!)
;       (Condon&Shortley 1951 p76-77, Abramowitz&Stegun 1972 p1006-1010)
;       where J=l1+l2+l3 is even and J=2g
;
; CALLING SEQUENCE:
;       Result = ThreeJzero( l1, l2, l3 )  
; INPUTS:
;       l1,l2,l3 = multipoles. Must be integers
;       
;
; EXAMPLE:
;       IDL> print,ThreeJzero(1,1,1)   ===>  0.00000  parity condition (1+1+1 is odd)
;       IDL> print,ThreeJzero(1,1,3)   ===>  0.00000  triangle condition (1+1<3) 
;       IDL> print,ThreeJzero(2,2,2)   ===>  -0.239046 analytically it's -sqrt(2./35)=-0.239046
;
;
; MODIFICATION HISTORY:
;       Written, F. Lacasa, IAS, August 2011
;  

J=l1+l2+l3
g=floor(J)/2
if (2*g eq J) and (abs(l2-l1) le l3) and (l3 le l1+l2) then begin
term1=lngamma(J-2*l1+1)+lngamma(J-2*l2+1)+lngamma(J-2*l3+1)-lngamma(J+1+1)
term2=lngamma(g+1)-lngamma(g-l1+1)-lngamma(g-l2+1)-lngamma(g-l3+1)
result=(-1)^g*exp(term1/2.+term2)
endif else begin
result=0.
endelse
return,result
end
