Function CHOLINVERT, A, status, check=check
;+
; NAME:
;       CHOLINVERT
; PURPOSE:
;       Computes the inverse of a symmetric positive-definite
;       matrix, by Cholesky decomposition. The inverse is also
;       symmetric (and positive-definite).
;
; CALLING SEQUENCE:
;       Inv = CHOLINVERT(A)  
; INPUTS:
;       A : square symmetric matrix positive-definite
;       
; OUTPUTS:
;       Inv = inverse, also symmetric, in double precision
;
; KEYWORDS:
;       /check : check if A is indeed a square symmetric matrix
;
; MODIFICATION HISTORY:
;       Written, F. Lacasa, IAS, February 2012
;  

error=-1
if keyword_set(check) then begin
sA=size(A)
if sA[0] ne 2 then begin
print,'A not a matrix'
return, error
endif
if sA[1] ne sA[2] then begin
print,'A rectangular'
return, error
endif
if total(A eq transpose(A)) ne n_elements(A) then begin
print,'A not symmetric'
return, error
endif
endif
bufA=A
nlines=(size(A))[1]
LA_CHOLDC,bufA,/double
for i=0,nlines-2 do bufA[i+1:*,i]=0
Linv=invert(bufA,stat,/double)
inv=Linv#transpose(Linv)
if arg_present(status) then status=stat
return, inv
END
