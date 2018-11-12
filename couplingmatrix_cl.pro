pro threej_init, lmax, logvec

j = indgen(3*lmax+2)
logvec = lngamma(2d*j+1d)-2d*lngamma(j+1d)

end

Function threejzerosq, l1, l2, l3, logvec
J=l1+l2+l3
g=J/2
truc=logvec[g-l1]+logvec[g-l2]+logvec[g-l3]-logvec[g]
result=exp(truc)/(J+1d)
return,result
END

PRO couplingmatrix_cl, mask, lmax, matrix, ordering=ord, verbose=verbose
t0=systime(1)

ianafast,mask,clmask,nlmax=lmax,simul_type=1,/silent,ordering=ord,/double

threej_init, lmax, logvec

ell=indgen(lmax+1)
matrix=dblarr(lmax+1,lmax+1)

for i=0,lmax do begin
   l1=ell[i]
   for j=i,lmax do begin
      l2=ell[j]
      ll=ell[abs(l1-l2):(l1+l2)<lmax]
      ll=ll[where(((l1+l2+ll) mod 2) eq 0)]
      nl=n_elements(ll)
      ll1=intarr(nl)+l1 & ll2=intarr(nl)+l2
      threejvec=threejzerosq(ll1,ll2,ll,logvec)
      truc=total(threejvec*(2d*ll+1d)*clmask[ll])
      matrix[i,j]=(2d*l2+1d)/(4d*!dpi)*truc
      matrix[j,i]=(2d*l1+1d)/(4d*!dpi)*truc

   endfor
   if keyword_set(verbose) and (i mod 100) eq 0 then begin
      print,format='($,A,I," ",F10.3)',string("15b),i,systime(1)-t0;"
   endif
endfor

if keyword_set(verbose) then print,' '

END
