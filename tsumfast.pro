FUNCTION TSUMFAST,X,Y               ;Fast trapezoidal summation
;+
; NAME:
;       TSUMFAST
; PURPOSE:
;       Trapezoidal summation of the area under a curve. 
; EXPLANATION:
;       Adapted from the procedure TSUM in the ASTRON procedure library.  
;
; CALLING SEQUENCE:
;       Result = TSUMFAST( x, y )  
; INPUTS:
;       x = array containing monotonic independent variable.
;       y = array containing dependent variable y = f(x)
; OUTPUTS:
;       result = area under the curve y=f(x) between xmin and xmax.
;
; EXAMPLE:
;       IDL> x = [0.0,0.1,0.14,0.3] 
;       IDL> y = sin(x)
;       IDL> print,tsumfast(x,y)    ===>  0.0445843
;       
;       In this example, the exact curve can be computed analytically as 
;       1.0 - cos(0.3) = 0.0446635     
;
; PROCEDURE:
;       The area is determined of individual trapezoids defined by x[i],
;       x[i+1], y[i] and y[i+1].
;
;       A more accurate integration can be done with the standard IDL
;       User Library int_tabulated.pro.
;
; MODIFICATION HISTORY:
;       Written, F. Lacasa, IAS, June 2011
;  

; Compute areas of trapezoids and sum result
    
  sum = total( ( y[1:*] + y )*( x[1:*] - x ) )/2. 

  return, sum

END    
