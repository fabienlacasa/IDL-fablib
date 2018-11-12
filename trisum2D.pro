FUNCTION TRISUM2D,X,Y,F               ;Fast tetrahedron summation in 2D
;+
; NAME:
;       TRISUM2D
; PURPOSE:
;       Tetrahedron summation of the volume under a surface.   
;
; CALLING SEQUENCE:
;       Result = TRISUM2D( x, y, f )  
; INPUTS:
;       x = 1D array containing monotonic independent variable.
;       y = 1D array containing monotonic independent variable.
;       f = 2D array containing f(x,y)
; OUTPUTS:
;       result = volume under the surface f(x,y) within the rectangle [xmin,xmax]*[ymin,ymax]
;
; EXAMPLE:
;       IDL> x = dindgen(51)/50
;       IDL> y = dindgen(101)/50
;       IDL> f = cos(x)#transpose(sin(y))
;       IDL> print, trisum2D(x,y,f), tsumfast(x,cos(x))*tsumfast(y,sin(y))
;             1.1839188       1.1915670
;       IDL> print,sin(1d)*(1.-cos(2d))
;             1.1916465
;
; PROCEDURE:
;       
;
; MODIFICATION HISTORY:
;       Written, F. Lacasa, IAS, December 2012
;  

x = double(temporary(x))
y = double(temporary(y))
f = double(temporary(f)) 

xx = x[1:*]-x
nx = n_elements(xx)
yy = y[1:*]-y
ny = n_elements(yy)
bary = (f[1:*,1:*]+2d*f[1:*,*]+2d*f[*,1:*]+f)

sum = total(bary*(xx#transpose(dblarr(ny)+1d))*((dblarr(nx)+1d)#transpose(yy)))/6d

return, sum
END
