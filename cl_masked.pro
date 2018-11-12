PRO Cl_masked, map, mask, lmax, clout, bin=bin ,weight=weight, ordering=ord,$
lmin=lmin, alreadymasked=already, mbb=mbb, clbin=cl_bin,mll=mll, P=P,Q=Q, $
clfull=cl_full, beam=beam, verbose=verbose

if ~keyword_set(bin) then begin
print,'Must give bin size'
return
endif
nbin=lmax/bin
lmoy=findgen(nbin)*bin+bin/2.
if keyword_set(ord) then begin
   mmap=map
   mmask=mask
endif else begin
   read_fits_map,map,mmap,/silent,nside=nsmap,ordering=ordmap
   read_fits_map,mask,mmask,/silent,nside=nsmask,ordering=ordmask
   if nsmap ne nsmask then begin
      print,'Map and Mask do not have the same nside'
      return
   endif
   if ordmap ne ordmask then begin
      print,'Map and Mask do not have the same ordering'
      return
   endif
   ord=ordmap
   nside=nsmap
endelse
if keyword_set(already) then begin
   masked_map=mmap
endif else begin
   masked_map=mmap*mmask
endelse
if ~keyword_set(lmin) then lmin=2 ;start from the quadrupole
if ~keyword_set(weight) then begin
   l=dindgen(lmax+1)
   weight=l*(l+1.)
endif
if ~keyword_set(beam) then beam=0.
gbeam=gaussbeam(beam,lmax)
beam_mat=(dblarr(lmax+1)+1)#gbeam^2

;compute full coupling matrix if necessary
if ~keyword_set(mll) then begin
   couplingmatrix_cl, mmask, lmax, mll, ordering=ord, verbose=verbose
endif

;bin it
P = dblarr(nbin,lmax+1)
Q = dblarr(lmax+1,nbin)

for i=0,nbin-1 do begin
lminbin = bin*i > lmin
lmaxbin = bin*(i+1)-1
binsize = double(lmaxbin-lminbin+1)
   for j=lminbin,lmaxbin do begin
      P[i,j] = weight[j]/binsize
      Q[j,i] = 1d/weight[j]
   endfor
endfor

mbb= P # ((mll*beam_mat) # Q)

;full power spectrum
ianafast,masked_map,cl_full,nlmax=lmax,simul_type=1,/silent,ordering=ord,/double
;remove first multipoles
if lmin gt 0 then cl_full[0:(lmin-1)]=0
;bin it
cl_bin= P # cl_full
;invert
clout=(invert(mbb) # cl_bin)/weight[lmoy]


END
