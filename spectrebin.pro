PRO SPECTREBIN, mapin,spectreout,clrw,lmax=lmax,bin=binell,op=op,$
filename=fileout,weights=weights,realw=realw,lmoy=lmoy,ordering=ord,verbose=verbose,$
help=help

;+
; NAME:
;       spectrebin
; PURPOSE:
;       Computes the binned power spectrum of a map in healpix .fits
;       format or online (then give ordering).
;
; CALLING SEQUENCE:
;       SPECTREBIN,map,spectreout,lmax=lmax,bin=binell
;
; INPUTS:
;       - map : Healpix map, either a fits or online (then precise
;         ordering)
;       - lmax : maximum multipole of the analysis
;         must be a multiple of binell
;       - binell : size of the ell bins (regular binning)
;
; OPTIONNAL KEYWORDS:
;       - help : displays the syntax of the procedure then stops
;       - verbose : displays progression of the computation (default : silent)
;       - ordering : precise the ordering for an online map (ring or nested)
;       - filename : precise a file where to save the power spectrum
;         (IDL save)
;       - weights : to multiply alm with to improve the binning
;         approximation
;       - op : if set uses the minimal resolution possible for a given
;         bin to speed up computations (i.e. 2*nside >= lmaxbin)
;       - lmoy : outputs the centers of the bins
;       
; OUTPUTS:
;       - spectreout : fltarr(nbin) containing the binned power spectrum
;
; MODIFICATION HISTORY:
;       Written, F. Lacasa, Institut d'Astrophysique Spatiale, February 2012
;


if keyword_set(help) then begin
print,'spectrebin,mapin,spectreout,lmax=lmax,'+$
'bin=binell [,filename=fileout,weights=weights,ordering=ord,/verbose,/help]'
return
endif
t0=systime(1)
if ~keyword_set(lmax) then begin
print,'Must declare lmax.'
return
endif
nsidemax=ell2nside(lmax)
nellbin=lmax/binell
if nellbin gt 1024 then begin
print,'Too many bins, too much memory needed.'
return
endif
if keyword_set(weights) then begin
if n_elements(weights) eq 1 then weights=fltarr(lmax)+weights[0]
if n_elements(weights) ne lmax then begin
print,'weights must have lmax elements'
return
endif
endif

lmoy=findgen(nellbin)*binell+binell/2.
npixmax=nside2npix(nsidemax)
omegapixmax=4*!pi/npixmax
tmpdir='/tmp/'

;map2alm
if keyword_set(verbose) then print,'anafasting input map...'
seed=systime(1)
randint=randomu(seed,/long)
filealmtmp=tmpdir+'alm_4spbin_tmp_'+strtrim(randint,2)+'.fits'
if keyword_set(ord) then begin
ianafast,mapin,alm1_out=filealmtmp,/silent,nlmax=lmax,ordering=ord
endif else begin
ianafast,mapin,alm1_out=filealmtmp,/silent,nlmax=lmax ;if no keyword
endelse
if keyword_set(verbose) then $
print,format='(A,F10.3)','computing scalemaps...',systime(1)-t0

spectreout=fltarr(nellbin)
fits2alm,index,alm,filealmtmp
index2lm,index,ll,mm
if keyword_set(weights) then begin
;remove monopole
wh=where(ll ne 0)
index=index[wh] & ll=ll[wh] & mm=mm[wh] & alm=alm[wh,*]
wweights=weights[ll-1]
alm=[[alm[*,0]*wweights],[alm[*,1]*wweights]]
endif


;scalemaps and spectrum
for i=0,nellbin-1 do begin
lminbin=binell*i+1
lmaxbin=binell*(i+1)
ell=findgen(binell)+lminbin ;lminbin..lmaxbin
prefspec=1./total(2.*ell+1.)
wh=where((ll ge lminbin) AND (ll le lmaxbin))
filealmscmap=tmpdir+'alm_4scmap_tmp_'+strtrim(randint,2)+'.fits'
alm2fits,index[wh],alm[wh,*],filealmscmap

if ~keyword_set(op) then begin
isynfast,0,scalemaptemp,alm_in=filealmscmap,nside=nsidemax,$
simul_type=1,lmax=lmaxbin+1,/silent
spectreout[i]=prefspec*total(scalemaptemp^2)*omegapixmax
endif else begin
nsideop=ell2nside(lmaxbin)
omegapixop=4*!pi/nside2npix(nsideop)
isynfast,0,scalemaptemp,alm_in=filealmscmap,nside=nsideop,$
simul_type=1,lmax=lmaxbin+1,/silent
spectreout[i]=prefspec*total(scalemaptemp^2)*omegapixop
endelse

if keyword_set(verbose) then $
print,format='($,A,I,"/",A," ",F10.3)',string("15b),i+1,strtrim(nellbin,2),systime(1)-t0;"
;print,strtrim(i+1,2)+'/'+strtrim(nellbin,2),systime(1)-t0
endfor
print,format='(a)','        '

;deweight if necessary
if keyword_set(weights) then begin
ell=findgen(lmax)+1
;wmoy=fltarr(nellbin)
;for i=0,nellbin-1 do wmoy[i]=weights[where(ell eq floor(lmoy[i]))]
wmoy=interpol(weights,ell,lmoy,/spline)
spectreout/=wmoy^2
endif

;if real weighted spectrum is asked
;if keyword_set(realw) then begin
if arg_present(clrw) then begin
if keyword_set(ord) then begin
ianafast,mapin,clraw,/silent,nlmax=lmax,ordering=ord
endif else begin
ianafast,mapin,clraw,/silent,nlmax=lmax ;if no keyword
endelse
if keyword_set(weights) then clfw=clraw[1:lmax]*weights^2 else clfw=clraw[1:lmax]
clrw=fltarr(nellbin)
for i=0,nellbin-1 do begin
lminbin=binell*i+1
lmaxbin=binell*(i+1)
if keyword_set(weights) then begin
clrw[i]=total(clfw[lminbin-1:lmaxbin-1])/total(weights[lminbin-1:lmaxbin-1]^2)
endif else begin
clrw[i]=total(clfw[lminbin-1:lmaxbin-1])/binell
endelse
endfor
endif

;save data if asked
if keyword_set(fileout) then begin
save,lmoy,spectreout,filename=fileout
endif

;clean temporary files
spawn,'rm -f '+filealmtmp+' '+filealmscmap

if keyword_set(verbose) then print,format='(A,F10.3)','Finished',systime(1)-t0
END
