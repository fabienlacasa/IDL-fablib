Function derivlsq, X, Y, Npts=Npts, errors=errors
;+
; NAME:
;	DERIVLSQ
;
; PURPOSE:
;	Perform numerical differentiation by chi-square fit of Y
;	with an affine function over npts. If no errors are given,
;	they are taken constant (least-square fit). Useful when the 
;	function has noise which may screw up the derivative.
;
; CATEGORY:
;	Numerical analysis.
;
; CALLING SEQUENCE:
;	Dy = Deriv(X,Y [,Npts=N, errors=err])
;
; INPUTS:
;	Y:  Variable to be differentiated.
;	X:  Variable to differentiate with respect to. 
;
; OPTIONAL INPUT PARAMETERS:
;	As above.
;
; OUTPUTS:
;	Returns the derivative in double precision.
;
; RESTRICTIONS:
;	
;
; HISTORY
;	written, F. Lacasa, IAS, 2011

on_error,2              ;Return to caller if an error occurs

if ~keyword_set(Npts) then Npts=3
if (Npts mod 2) ne 1 then message,'Npts must be odd'
k=(Npts-1)/2 ;Npts=2k+1

n = n_elements(x)
d=dblarr(n)
if n ne n_elements(y) then message,'Vectors must have same size'
if Npts gt n then message,'Npts must be less than vectors size'
if ~keyword_set(errors) then errors=dblarr(n)+1d

;dy/dx(i) = [sum_j(1/sigma(j)^2)*sum_j(x(j)*y(j)/sigma(j)^2) - sum_j(x(j)/sigma(j)^2)*sum_j(y(j)/sigma(j)^2)] / 
; [sum_j(1/sigma(j)^2)*sum_j(x(j)^2/sigma(j)^2) - (sum_j(x(j)/sigma(j)^2))^2] 
; where j=i-k..i+k

;ensure double precision
xx = double(x)
yy = double(y)-mean(double(y))  ;substract the mean to avoid round-off errors
err = double(errors)

xtable=fltarr(n,Npts)
ytable=fltarr(n,Npts)
errtable=fltarr(n,Npts)
for l=-k,k do begin
xtable(*,l+k)=shift(xx,l)-xx
ytable(*,l+k)=shift(yy,l)
errtable(*,l+k)=shift(err,l)
endfor
d=[total(1d/errtable^2,2)*total(xtable*ytable/errtable^2,2)-$
total(xtable/errtable^2,2)*total(ytable/errtable^2,2)]/$
[total(1d/errtable^2,2)*total(xtable^2/errtable^2,2)-total(xtable/errtable^2,2)^2]

;first point
x01=xx[0]/err[0]-xx[1]/err[1]
x02=xx[0]/err[0]-xx[2]/err[2]
x12=xx[1]/err[1]-xx[2]/err[2]
d[0] = yy[0]*(x01+x02)/(x01*x02) - yy[1]*x02/(x01*x12) + yy[2]*x01/(x02*x12)
;last point
x01=x[n-3]/err[n-3]-x[n-2]/err[n-2]
x02=x[n-3]/err[n-3]-x[n-1]/err[n-1]
x12=x[n-2]/err[n-2]-x[n-1]/err[n-1]
d[n-1] = -yy[n-3]*x12/(x01*x02) + yy[n-2]*x02/(x01*x12) - yy[n-1]*(x02+x12)/(x02*x12)
;make chi-square fit with number of available points for other
;extreme points
for i=1,k-1 do begin
kk=i
nnpts=2.*kk+1
xtable=xx[0:2*i]-xx[i]
ytable=yy[0:2*i]
errtable=err[0:2*i]
d[i]=[total(1d/errtable^2)*total(xtable*ytable/errtable^2)-total(xtable/errtable^2)*total(ytable/errtable^2)]/$
[total(1d/errtable^2)*total(xtable^2/errtable^2)-total(xtable/errtable^2)^2]
endfor
for i=n-k,n-2 do begin
kk=n-1-i
nnpts=2.*kk+1
xtable=xx[(i-kk):n-1]-xx[i]
ytable=yy[(i-kk):n-1]
d[i]=[nnpts*total(xtable*ytable)-total(xtable)*total(ytable)]/$
[nnpts*total(xtable^2)-total(xtable)^2]
endfor

return, d
end
