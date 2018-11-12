FUNCTION POOR_INPAINTER, map, mask, niter=nitermax, nested=nested, verbose=verbose
t0=systime(1)
npix=n_elements(map)
nside=npix2nside(npix)
if n_elements(mask) ne npix  then begin
   print,'POOR_INPAINTER : error, mask and map must have same number of elements'
   return,map
endif
wh=where((mask eq 1) or (mask eq 0),cwh)
if cwh ne n_elements(mask) then begin
   print,'POOR_INPAINTER : error, mask must have either 0 or 1'
   return,map
endif
if ~keyword_set(nitermax) then nitermax=100 ;maximum number of iterations

inpainted_map=map
goodpixel=mask

niter=0
whmasked=where(goodpixel eq 0,nwhmasked)

while (niter lt nitermax) AND (nwhmasked gt 0) do begin
   if keyword_set(verbose) AND ((niter mod 1) eq 0) then print,$
     'niter: '+strtrim(niter,2)+'  nwhmasked: '+strtrim(nwhmasked,2)

   for i=0,nwhmasked-1 do begin
      pixel=whmasked[i]
      neighbours_ring,nside,pixel,neighbours
      if keyword_set(nested) then neighbours_nest,nside,pixel,neighbours
      whneigh=where(goodpixel[neighbours] eq 1,ngoodneigh)
      if ngoodneigh ne 0 then begin
         neighbours=neighbours[whneigh]
         inpainted_map[pixel]=mean(inpainted_map[neighbours])
         goodpixel[pixel]=1
      endif
      if keyword_set(verbose) AND ((i mod 100000) eq 0) then print,$
        strtrim(i,2)+'/'+strtrim(nwhmasked,2),systime(1)-t0
   endfor

   whmasked=where(goodpixel eq 0,nwhmasked)
   niter+=1
endwhile

if (niter eq nitermax) AND keyword_set(verbose) then print,'Maximum number of iterations reached'

return, inpainted_map
END
