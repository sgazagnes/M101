cspeed = 2.9979E5
!PATH = !PATH +':'+Expand_path('+/Users/johnchisholm/Documents/old_mac/admin/Desktop/idllib/') + Expand_path('+/Users/johnchisholm/Documents/old_mac/admin/Desktop/idllib/idlutils/pro/kcorrect/pro/utils/pro')
gal = 'REG4'
     
readcol, '/Users/johnchisholm/Documents/M101/'+gal+'/'+gal+'.dat',$
     owave1, iflux1, ierr1
z = 0.00052173913

wave = owave1 / (1.+z)

nflux1 = iflux1/median(iflux1[where(wave gt 1267 and wave lt 1276)])
nerr1 = ierr1/median(iflux1[where(wave gt 1267 and wave lt 1276)])

super_mods= mrdfits('/Users/johnchisholm/megasaura/super_mods.fits', 1)
;bp_mod= mrdfits('/Users/johnchisholm/megasaura/BPASS/super_mods_bp.fits', 1)

;nflux = iflux
;nerr = ierr
;
;mod001 = mrdfits('/Volumes/Transcend/HST/metal/newsb/sb99001.fits', 1)
;mod004 = mrdfits('/Volumes/Transcend/HST/metal/newsb/sb99004.fits', 1)
;mod008 = mrdfits('/Volumes/Transcend/HST/metal/newsb/sb99008.fits', 1)
;mod02 = mrdfits('/Volumes/Transcend/HST/metal/newsb/sb9902.fits', 1)
;;mod04 = mrdfits('/Volumes/Transcend/HST/metal/newsb/sb9904.fits', 1)
;mod04 = mrdfits('/Users/johnchisholm/megasaura/new_super_solar.fits', 1)
;mod001 = mrdfits('/Users/johnchisholm/megasaura/w_nebcont/sb99001_wnebcont_logU-2.5.fits', 1)
;mod004 = mrdfits('/Users/johnchisholm/megasaura/w_nebcont/sb99004_wnebcont_logU-2.5.fits', 1)
;mod008 = mrdfits('/Users/johnchisholm/megasaura/w_nebcont/sb99008_wnebcont_logU-2.5.fits', 1)
;mod02 = mrdfits('/Users/johnchisholm/megasaura/w_nebcont/sb9902_wnebcont_logU-2.5.fits', 1)
;;mod04 = mrdfits('/Volumes/Transcend/HST/metal/newsb/sb9904.fits', 1)
;mod04 = mrdfits('/Users/johnchisholm/megasaura/w_nebcont/sb9904_wnebcont_logU-2.5.fits', 1)
;modo6 = mrdfits('/Volumes/Transcend/megasaura/modelsO6.fits', 1)
;modcon = mrdfits('/Users/johnchisholm/megasaura/chuck/004CSP.fits', 1)
;nebular = mrdfits('/Users/admin/Documents/megasaura/nebular_cont_model.fits', 1)
;
;mode001 = mod001
;mode004 = mod004
;mode008 = mod008
;mode02 = mod02
;mode04 = mod04
;
;
;
;mode001 = {wave: mode001.wave, flux: mod001.FLAM_STELLAR_NEBCONT,age: mode001.age,sfr: mode001.sfr, id: mode001.id, norm: mode001.norm}
;mode004 = {wave: mode004.wave, flux: mod004.FLAM_STELLAR_NEBCONT,age: mode004.age,sfr: mode004.sfr, id: mode004.id, norm: mode004.norm}
;mode008 = {wave: mode008.wave, flux: mod008.FLAM_STELLAR_NEBCONT,age: mode008.age,sfr: mode008.sfr, id: mode008.id, norm: mode008.norm}
;mode02 = {wave: mode02.wave, flux: mod02.FLAM_STELLAR_NEBCONT,age: mode02.age,sfr: mode02.sfr, id: mode02.id, norm: mode02.norm}
;mode04 = {wave: mode04.wave, flux: mod04.FLAM_STELLAR_NEBCONT,age: mode04.age,sfr: mode04.sfr, id: mode04.id, norm: mode04.norm}
;
;
;for modnum = 0, n_elements(mode001.age)-1 do begin & $
; modnormreg = where(mode001.wave gt 1267 and mode001.wave lt 1276) & $
;  mode001.flux[modnum,*] = mode001.flux[modnum, *]/median(mode001.flux[modnum, modnormreg]) & $
;end
;
;for modnum = 0, n_elements(mode004.age)-1 do begin & $
;  modnormreg = where(mode004.wave gt 1267 and mode004.wave lt 1276) & $
;  mode004.flux[modnum,*] = mode004.flux[modnum, *]/median(mode004.flux[modnum, modnormreg]) & $
;end
;
;for modnum = 0, n_elements(mode008.age)-1 do begin & $
;  modnormreg = where(mode008.wave gt 1267 and mode008.wave lt 1276) & $
;  mode008.flux[modnum,*] = mode008.flux[modnum, *]/median(mode008.flux[modnum, modnormreg]) & $
;end
;
;for modnum = 0, n_elements(mode02.age)-1 do begin & $
;  modnormreg = where(mode02.wave gt 1267 and mode02.wave lt 1276) & $
;  mode02.flux[modnum,*] = mode02.flux[modnum, *]/median(mode02.flux[modnum, modnormreg]) & $
;end
;
;for modnum = 0, n_elements(mode04.age)-1 do begin & $
;  modnormreg = where(mode04.wave gt 1267 and mode04.wave lt 1276) & $
;  mode04.flux[modnum,*] = mode04.flux[modnum, *]/median(mode04.flux[modnum, modnormreg]) & $
;end
;
;for modnum = 0, n_elements(mode04.age)-1 do begin & $
;  modnormreg = where(mode04.wave gt 1267 and mode04.wave lt 1276) & $
;  mode04.flux[modnum,*] = mode04.flux[modnum, *]/median(mode04.flux[modnum, modnormreg]) & $
;end
;
;for modnum = 0, n_elements(modecon.age)-1 do begin & $
;  modnormreg = where(modecon.wave gt 1267 and modecon.wave lt 1276) & $
;  modecon.flux[modnum,*] = modecon.flux[modnum, *]/median(modcon.flux[modnum, modnormreg]) & $
;end
;
;for modnum = 0, n_elements(modeo6.age)-1 do begin & $
;  modnormreg = where(modeo6.wave gt 1045 and modeo6.wave lt 1055) & $
;  modeo6.flux[*,modnum] = modeo6.flux[*, modnum]/median(modeo6.flux[ modnormreg, modnum]) & $
;end

s1 = {wave: wave, flux: nflux1, err: nerr1}
s = create_struct(s1, 'name', gal, 'z', z)
readcol, '/Users/johnchisholm/megasaura/stack-A/stacked.linelist', $
  mwlines1, id, dumbywave, fval, dumby, color, zoff, type, format = 'f, a, f, f,f, a, f, a'
dv = 700.
mask = fltarr(n_elements(wave))+1.
mwlines = mwlines1[where(strmatch(type, 'ISM') eq 1 or strmatch(type, 'INTERVE') eq 1 or strmatch(type, 'EMISSION') eq 1)]
;mwlines = mwlines1[where(strmatch(type, 'INTERVE') eq 1 or strmatch(type, 'ISM') eq 1 or strmatch(type, 'FINESTR') eq 1 or strmatch(type, 'PHOTOSPHERE') eq 1)]
;mwlines = mwlines1[where(strmatch(type, 'ISM') eq 1)]
if strmatch(gal, 'REG4') eq 1 then mwlines = [mwlines, 1303, 1305, 1307, 1310, 1643., 1900.]




nmwlines = n_elements(mwlines)

for maskn = 0, nmwlines-1 do begin & $
    tempmask = where(wave gt -dv/cspeed*mwlines[maskn]+mwlines[maskn] and wave lt dv/cspeed*mwlines[maskn]+mwlines[maskn])   & $
    if mwlines[maskn] eq 1550.77 then tempmask = where(wave gt -dv/cspeed*mwlines[maskn]+mwlines[maskn] and wave lt 50./cspeed*mwlines[maskn]+mwlines[maskn])   & $
    ;if mwlines[maskn] eq 1550.77 then tempmask = where(wave gt -2700./cspeed*mwlines[maskn]+mwlines[maskn] and wave lt 50./cspeed*mwlines[maskn]+mwlines[maskn])   & $
    ;if maskn eq 41 then tempmask = where(wave gt -2700./cspeed*mwlines[maskn]+mwlines[maskn] and wave lt 50./cspeed*mwlines[maskn]+mwlines[maskn])   & $
    ;if maskn eq 41 then tempmask = where(wave gt -2700./cspeed*mwlines[maskn]+mwlines[maskn] and wave lt 50./cspeed*mwlines[maskn]+mwlines[maskn])   & $
    ;if maskn eq 56 then tempmask = where(wave gt -2700./cspeed*mwlines[maskn]+mwlines[maskn] and wave lt 50./cspeed*mwlines[maskn]+mwlines[maskn])   & $
    mask[tempmask] = 0 & $
endfor

;deflines, s.flux, interstellar = ism, stellar_photo = spho, stellar_wind = swin
;mask[ism] = 0
;mask[spho] = 0
;mask[swin] = 0

chisreg = where(mask eq 1  and s.wave gt 1300. and s.wave lt 1950 and finite(s.flux) eq 1 and s.err ne 0, comp = badc)
; and  $ ;,comp = badc); and $ ;These are galaxy specific masks
   ;(s.wave lt 1440 or s.wave gt 1455) , comp = badc)
s.err[badc] = 0

plot, s.wave, s.flux, yr = [0,2], xr = [1300, 2000]
djs_oplot, s.wave[chisreg], s.flux[chisreg], color = 'green'

vdisp = cspeed/1900.*2.355
c_super_mods = s99_continuum(super_mods, s.wave, s.flux, s.err, vdisp, yfit = cfit_super, /noplot)


age = fltarr(100)
met = fltarr(100)
ebv = fltarr(100)
for i = 0, 99 do begin & $
  temp_flux = s.flux+randomn(seed, n_elements(s.flux))*s.err & $
  c_super_mods_temp = s99_continuum(super_mods, s.wave, temp_flux, s.err, vdisp, yfit = cfit_super_temp, /noplot) & $
  ;f_rat[i] = temp_fit[where(super_mods.wave gt 894 and super_mods.wave lt 905), i]/temp_fit[where(super_mods.wave gt 1494 and super_mods.wave lt 1496), i] & $
  met[i] = total(c_super_mods_temp.light_frac[0:49]*super_mods.z[0:49])/total(c_super_mods_temp.light_Frac[0:49]) & $
  age[i] = total(c_super_mods_temp.light_frac[0:49]*super_mods.age[0:49])/total(c_super_mods_temp.light_Frac[0:49]) & $
  ebv[i] = c_super_mods_temp.ebv & $
endfor

mwrfits, age, '/Users/johnchisholm/Documents/M101/'+gal+'/age.fits
mwrfits, met, '/Users/johnchisholm/Documents/M101/'+gal+'/met.fits


newvoff = velshift(s.wave, s.flux, s.err, cfit_super, z)
newvoff_b = velshift(owave1, s.flux, s.err, bfit_super, z)

newvoff_b.voff = 0.
stdvoff = sqrt(total(newvoff.chis*(newvoff.v_off-newvoff.voff)^2))
nwave = owave1/(1.+z+newvoff.voff/cspeed)
linterp, owave1/(1.+z), cfit_super, owave1/(1.+z+newvoff.voff/cspeed), fit_super
linterp, owave1/(1.+z), bfit_super, owave1/(1.+z+newvoff_B.voff/cspeed), fit_super_b
flux = {wave: wave, flux: nflux1, err: nerr1}
flux.wave = owave1/(1.+z+newvoff.voff/cspeed)
cfit=cfit_super
print, total((s.flux[chisreg]-cfit_super[chisreg])^2/s.err[chisreg]^2)/n_elements(chisreg)
crs = {wave: wave, flux: nflux1, err: nerr1}
cr = create_struct(crs, 'name', gal, 'z', z)

plot, s.wave, s.flux,  yr = [0, 2], xr= [950, 2050]
djs_oplot, s.wave[chisreg], s.flux[chisreg], color = 'green'
djs_oplot, s.wave, cfit_super, color = 'red'

!p.multi = [0, 1, 2]
plot, s.wave, s.flux, xr = [1450, 1575], yr = [0, 2], thick = 2
djs_oplot, s.wave[chisreg], s.flux[chisreg], color = 'green'
djs_oplot, s.wave, cfit_super, color = 'light red', thick = 2

plot, s.wave, s.flux, xr = [1480, 1520], yr = [0, 2], thick = 2
djs_oplot, s.wave, cfit_super, color = 'blue', thick = 2

plot, s.wave, s.flux, xr = [1350, 1450], yr = [0, 2], thick = 2
djs_oplot, s.wave[chisreg], s.flux[chisreg], color = 'green'
djs_oplot, s.wave, cfit_super, color = 'blue', thick = 2

plot, s.wave, s.flux, xr = [1250, 1350], yr = [0, 2], thick = 2
djs_oplot, s.wave[chisreg], s.flux[chisreg], color = 'green', thick = 2
djs_oplot, s.wave, cfit_super, color = 'blue', thick = 2

plot, s.wave, s.flux, xr = [1850, 1950], yr = [0, 2]
djs_oplot, s.wave, cfit_super, color = 'blue', thick = 2

print, c_super_mods.ebv, b_super_mods.ebv
print, total(c_super_mods.light_frac*super_mods.z)/total(c_super_mods.light_Frac ),total(b_super_mods.light_frac*bp_mod.z)/total(b_super_mods.light_Frac )
exz = total(c_super_mods.light_frac *super_mods.z )/total(c_super_mods.light_Frac )
exz_b = total(b_super_mods.light_frac *bp_mod.z )/total(b_super_mods.light_Frac )

stdz = total(c_super_mods.light_frac*(super_mods.z-exz)^2)/total(c_super_mods.light_Frac)
print, total(c_super_mods.light_frac *super_mods.age )/total(c_super_mods.light_Frac )/10.^6,total(b_super_mods.light_frac *bp_mod.age )/total(b_super_mods.light_Frac )/10.^6
exage = total(c_super_mods.light_frac *super_mods.age )/total(c_super_mods.light_Frac )
exage_b = total(b_super_mods.light_frac *bp_mod.age )/total(b_super_mods.light_Frac )

  
newvoff = velshift(owave1[chisreg], s.flux[chisreg], s.err[chisreg], cfit_super[chisreg], z)
;newvoff.voff = 0.
stdvoff = sqrt(total(newvoff.chis*(newvoff.v_off-newvoff.voff)^2))
flux = {wave: wave, flux: nflux1, err: nerr1}
flux.wave = owave1/(1.+z+newvoff.voff/cspeed)
;linterp, owave1/(1.+z), c, owave1/(1.+z+newvoff.voff/cspeed), cfit


newvoff_b = velshift(owave1[chisreg], s.flux[chisreg], s.err[chisreg], bfit_super[chisreg], z)
newvoff_b.voff = 0.
stdvoff = sqrt(total(newvoff.chis*(newvoff.v_off-newvoff.voff)^2))
flux_b = {wave: wave, flux: nflux1, err: nerr1}
flux_b.wave = owave1/(1.+z+newvoff_b.voff/cspeed)



dfpsplot, '/Users/johnchisholm/Documents/M101/'+gal+'/C4-bestfit-vel.ps', /color, /create
!p.font = 0
!p.multi = [0, 1, 2]
angstrom = cgsymbol('angstrom')
plot, (flux.wave[where(s.flux gt 0)]-1548.1950)/1548.1950*cspeed, smooth(flux.flux[where(flux.flux gt 0)], 2), yr = [0, 2], xr = [-10000, 10000], psym = 10, /xs,$
  xtit = 'Velocity (km/s)', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, (flux.wave[where(flux.flux gt 0)]-1548.1950)/1548.1950*cspeed, cfit[where(flux.flux gt 0)], color = 'red', thick = 7

plot, (flux.wave[where(flux.flux gt 0)]-1548.1950)/1548.1950*cspeed, smooth(flux.flux[where(flux.flux gt 0)], 2)/smooth(cfit[where(flux.flux gt 0)],2), yr = [0, 2], xr = [-10000, 10000], psym = 10, /xs,$
  xtit = 'Velocity (km/s)', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, [-10000, 10000], [1, 1], thick = 7
dfpsclose


dfpsplot, '/Users/johnchisholm/Documents/M101/'+gal+'/S4-bestfit.ps', /color, /create
!p.font = 0
!p.multi = [0, 1, 2]
angstrom = cgsymbol('angstrom')
plot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], smooth(flux.flux[where(flux.flux gt 0 and flux.flux lt 5)], 2), yr = [0, 2], xr = [1350, 1450], psym = 10, /xs,$
  xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], cfit[where(flux.flux gt 0 and flux.flux lt 5)], color = 'red', thick = 7

plot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], smooth(flux.flux[where(flux.flux gt 0 and flux.flux lt 5)], 2)/smooth(cfit[where(flux.flux gt 0 and flux.flux lt 5)],2), yr = [0, 2], xr = [1350, 1450], psym = 10, /xs,$
  xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, [1200, 10000], [1, 1], thick = 7
dfpsclose


dfpsplot, '/Users/johnchisholm/Documents/M101/'+gal+'/photo.ps', /color, /create
!p.font = 0
!p.multi = [0, 1, 2]
angstrom = cgsymbol('angstrom')
plot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], smooth(flux.flux[where(flux.flux gt 0 and flux.flux lt 5)], 2), yr = [0, 2], xr = [1250, 1350], psym = 10, /xs,$
  xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], cfit[where(flux.flux gt 0 and flux.flux lt 5)], color = 'red', thick = 7

plot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], smooth(flux.flux[where(flux.flux gt 0 and flux.flux lt 5)], 2)/smooth(cfit[where(flux.flux gt 0 and flux.flux lt 5)],2), yr = [0, 2], xr = [1250, 1350], psym = 10, /xs,$
  xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, [1200, 10000], [1, 1], thick = 7
dfpsclose



dfpsplot,  '/Users/johnchisholm/Documents/M101/'+gal+'/full.ps', /color, /create
!p.font = 0
!p.multi = [0, 1, 2]
angstrom = cgsymbol('angstrom')
plot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], smooth(flux.flux[where(flux.flux gt 0 and flux.flux lt 5)], 2), yr = [0, 2], xr = [1200, 2500], psym = 10, /xs,$
  xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], cfit_super[where(flux.flux gt 0 and flux.flux lt 5)], color = 'red', thick = 7
plot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], smooth(flux.flux[where(flux.flux gt 0 and flux.flux lt 5)], 2)/smooth(cfit_super[where(flux.flux gt 0 and flux.flux lt 5)],2), yr = [0, 2], xr = [1200, 2500], psym = 10, /xs,$
  xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, [1200, 1700], [1, 1], thick = 7
dfpsclose


dfpsplot, '/Users/johnchisholm/Documents/M101/'+gal+'/rix.ps', /color, /create
!p.font = 0
!p.multi = [0, 1, 2]
angstrom = cgsymbol('angstrom')
plot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], smooth(flux.flux[where(flux.flux gt 0 and flux.flux lt 5)], 2), yr = [0, 1], xr = [1900, 2100], psym = 10, /xs,$
  xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], cfit[where(flux.flux gt 0 and flux.flux lt 5)], color = 'red', thick = 7
djs_oplot, [1978, 1978], [0, 2], thick =4


plot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], smooth(flux.flux[where(flux.flux gt 0 and flux.flux lt 5)], 2)/smooth(cfit[where(flux.flux gt 0 and flux.flux lt 5)],2), yr = [0, 2], xr = [1900, 2100], psym = 10, /xs,$
  xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, [1200, 10000], [1, 1], thick = 7
djs_oplot, [1978, 1978], [0, 2], thick = 4

plot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], smooth(flux.flux[where(flux.flux gt 0 and flux.flux lt 5)], 2), yr = [0, 2], xr = [1400, 1450], psym = 10, /xs,$
  xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], cfit[where(flux.flux gt 0 and flux.flux lt 5)], color = 'red', thick = 7
djs_oplot, [1978, 1978], [0, 2], thick =4


plot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], smooth(flux.flux[where(flux.flux gt 0 and flux.flux lt 5)], 2)/smooth(cfit[where(flux.flux gt 0 and flux.flux lt 5)],2), yr = [0, 2], xr = [1400, 1450], psym = 10, /xs,$
  xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, [1200, 10000], [1, 1], thick = 7
djs_oplot, [1978, 1978], [0, 2], thick = 4

dfpsclose

dfpsplot, '/Users/johnchisholm/megasaura/'+gal+'/BPASS/he2.ps', /color, /create
!p.font = 0
!p.multi = [0, 1, 2]
angstrom = cgsymbol('angstrom')
plot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], smooth(flux.flux[where(flux.flux gt 0 and flux.flux lt 5)], 2), yr = [0, 2], xr = [1620, 1680], psym = 10, /xs,$
  xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], cfit_super[where(flux.flux gt 0 and flux.flux lt 5)], color = 'red', thick = 7
plot, flux.wave[where(flux.flux gt 0 and flux.flux lt 5)], smooth(flux.flux[where(flux.flux gt 0 and flux.flux lt 5)], 2)/smooth(cfit_super[where(flux.flux gt 0 and flux.flux lt 5)],2), yr = [0, 2], xr = [1620, 1680], psym = 10, /xs,$
  xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5, thick = 10
djs_oplot, [1200, 1700], [1, 1], thick = 7
dfpsclose




dfpsplot,'/Users/johnchisholm/megasaura/'+gal+'/BPASS/v_offset.ps', /color, /landscape
!p.multi = 0
plot, newvoff.v_off, newvoff.chis, thick = 10, psym = 10, xtit = 'Velocity offset', ytit = 'Probability', charsize = 1.5, charthick = 5 

dfpsclose
get_date, date

openw, lun,   '/Users/johnchisholm/Documents/M101/'+gal+'/'+gal+'-sb99-fit.txt', /get_lun
printf, lun, '# # SB99 Continuum fits to '+gal+' with a nebular continuum model and log(U) = -2.5.' 
printf, lun, '# # Created on (yyyy-mm-dd)' + date
printf, lun, '# # The fit was created using the parameters in the '+gal+'continuum-parameters.fits file'
printf, lun, '# # The first column is the wavelength, the second and third columns are the observed flux and error (in f_lam) normalized by the median flux between 1267-1276, the last column is the linear-combination Starburst99 fit to the flux.'
printf, lun, '# # The velocity dispersion (sigma) used to create this template is ' + strtrim(string(vdisp), 2) + ' km/s'
printf, lun, '# # The stellar continuum is offset by '+ strtrim(string(newvoff.voff), 2) + ' km/s'
printf, lun, '# # The fit was only made from 1240-2000, while masking out between 1303, 1305, 1307, 1310, 1643., 1900A. Be skeptical of the fit outside of these wavelengths.'
for i =0, n_elements(flux.wave)-1 do printf, lun, flux.wave[i], flux.flux[i], flux.err[i], cfit_super[i]
free_lun, lun




;c_super_mods = mrdfits('/Users/admin/Documents/megasaura/'+gal+'/continuum-properties.fits', 1)
output = {ebv: c_super_mods.ebv, ebv_err: c_super_mods.ebv_err, weighted_z: exz, $
            weighted_age: exage/10.^6,  light_frac: c_super_mods.light_Frac/total(c_super_mods.light_frac), $
            light_frac_err: c_super_mods.light_frac_err/total(c_super_mods.light_frac), model_age: c_super_mods.model_age, $
            model_z: super_mods.z}
            

mwrfits, output, '/Users/johnchisholm/Documents/M101/'+gal+'/'+gal+'-SB99-continuum-properties.fits', /create
mwrfits, output_b, '/Users/johnchisholm/megasaura/'+gal+'/BPASS/'+gal+'-BPASS-continuum-properties.fits', /create


nvoff = {z: (z+newvoff.voff/cspeed), voff: newvoff.voff, chis: newvoff.chis, v_off: newvoff.v_off}
mwrfits, nvoff, '/Users/johnchisholm/megasaura/'+gal+'/voff.fits', /create

mwrfits, cfo6, '/Users/admin/Documents/megasaura/'+gal+'/o6/cont-o6-fits.fits', /create
openw, lun,   '/Users/admin/Documents/megasaura/'+gal+'/o6/s1226-o6.txt', /get_lun
for i =0, n_elements(o6reg)-1 do printf, lun, flux.wave[o6reg[i]], flux.flux[o6reg[i]]/median(flux.flux[where(flux.wave gt 1050 and flux.wave lt 1055)]),$
   flux.err[o6reg[i]]/median(flux.flux[where(flux.wave gt 1050 and flux.wave lt 1055)]), cfito6[i]
free_lun, lun

readcol,  '/Users/admin/Documents/megasaura/'+gal+'/s1226-o6.txt', o6w, o6f, o6e, cfito6
so6 = {wave: o6w, flux: o6f, err: o6e}


dfpsplot, '/Users/admin/Documents/megasaura/'+gal+'/o6/O6-cont1.ps', /color, /landscape
;!p.multi = [0, 1,2]
!p.font = 0
angstrom = cgsymbol('angstrom')
plot, so6.wave[o6reg], so6.flux[o6reg], xr = [1015, 1045], thick = 10, charthick = 10, charsize = 2, $
    yr = [0, 2], /xs, xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'Normalized Flux'
djs_oplot, so6.wave[o6reg], cfito6, color = 'red', thick = 10
fitlab = textoidl('O VI SB99 + Ly\beta Fit')
legend, ['Observed Flux', fitlab], textcolors = [0, djs_icolor('red')], charsize = 2.1, charthick = 10, box = 0
djs_oplot, [1031.9261, 1031.9261], [1, 2], thick = 10, linestyle = 2
djs_oplot, [1037.6167, 1037.6167], [1, 2], thick = 10, linestyle = 2
xyouts, 1027, 1.4,'O VI 1032', charthick = 10, charsize = 1.5
xyouts, 1038, 1.4,'O VI 1038', charthick = 10, charsize = 1.5
;plot, so6.wave[o6reg], so6.flux[o6reg]/cfito6, xr = [1010, 1050], thick = 10, charthick = 10, charsize = 2, $
;  yr = [0, 3.5], /xs, xtit = 'Rest Wavelength ('+angstrom+')', ytit = 'SB99 Normalized Flux', /ys, $
;  ytickv = [0, 1, 2, 3]
;  !p.multi = 0
;djs_oplot, [0, 10000], [1, 1], thick = 7
;djs_oplot, [1031.9261, 1031.9261], [1.1, 3.5], thick = 10, linestyle = 2
;djs_oplot, [1037.6167, 1037.6167], [1.1, 3.5], thick = 10, linestyle = 2
;djs_oplot, [1036.3367, 1036.3367], [1.1, 3.5], thick = 10, linestyle = 3, color = 'orange'
;djs_oplot, [1039.2394, 1039.2394], [1.1, 3.5], thick = 10, linestyle = 3, color = 'magenta'
;djs_oplot, [1048.2199, 1048.2199], [1.1, 3.5], thick = 10, linestyle = 3, color = 'red'
;djs_oplot, [1020.6989, 1020.6989], [1.1, 3.5], thick = 10, linestyle = 3, color = 'blue'
;djs_oplot, [1025.7223, 1025.7223], [1.1, 3.5], thick = 10, linestyle = 3, color = 'cyan'

dfpsclose



sregs1 = where(((flux.wave-1845.51)/1845.51*cspeed lt -400 and (flux.wave-1845.51)/1845.51*cspeed gt -1000) or ((flux.wave-1845.51)/1845.51*cspeed lt 400 and (flux.wave-1845.51)/1845.51*cspeed gt 100))
sfits1 = poly_fit(flux.wave[sregs1], flux.flux[sregs1]/cfit[sregs1], 1, yfit = reconts1)
recont_full_s1 = interpol( reconts1,flux.wave[sregs1], flux.wave) 
plot, (flux.wave-1845.51)/1845.51*cspeed, flux.flux/cfit/recont_full_s1, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
mins1 = -130


djs_oplot, [mins1, mins1], [0, 2], color = 'red'
maxs1 = 75

djs_oplot, [maxs1, maxs1], [0, 2], color = 'red'

fitregs1 = where(abs((flux.wave-1845.51)/1845.51*cspeed) lt 1000)
vs1 = mcmcv(flux.wave[fitregs1], flux.flux[fitregs1]/cfit[fitregs1]/recont_full_s1[fitregs1], flux.err[fitregs1]/cfit[fitregs1]/recont_full_s1[fitregs1], 1845.51, maxv=maxs1, minv=mins1)
djs_oplot, vs1.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vs1.v90*[1,1], [0, 2], color ='orange'
ns1_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1845.51, 0.258, mins1, maxs1)
ns1 = alog10(ns1_f.n)
ns1_e = sqrt(ns1_f.e^2/ns1_f.n^2/alog(10.)^2)
print, ns1, ns1_e


plot, (flux.wave-1302.168)/1302.168*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
mino1 = -490


djs_oplot, [mino1, mino1], [0, 2], color = 'red'
maxo1 = 100

djs_oplot, [maxo1, maxo1], [0, 2], color = 'red'

fitreg3 = where(abs((flux.wave-1302.168)/1302.168*cspeed) lt 1000)
vo1 = mcmcv(flux.wave[fitreg3], flux.flux[fitreg3]/cfit[fitreg3], flux.err[fitreg3]/cfit[fitreg3], 1302.168, maxv=maxo1, minv=mino1)
djs_oplot, vo1.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vo1.v90*[1,1], [0, 2], color ='orange'
no1_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1302.168, 0.052, mino1, maxo1)
no1 = alog10(no1_f.n)
no1_e = sqrt(no1_f.e/no1_f.n/alog(10.)^2)
print, no1, no1_e

plot, (flux.wave-1608.45708)/1608.45078*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
minfe = -375


djs_oplot, [minfe, minfe], [0, 2], color = 'red'
maxfe = 70

djs_oplot, [maxfe, maxfe], [0, 2], color = 'red'

fitregfe = where(((flux.wave-1608.45708)/1608.45708*cspeed) lt maxfe and (flux.wave-1608.45708)/1608.45708*cspeed gt minfe)
vfe = mcmcv(flux.wave[fitregfe], flux.flux[fitregfe]/cfit[fitregfe], flux.err[fitregfe]/cfit[fitregfe], 1608.45078, maxv=maxfe, minv=minfe)
djs_oplot, vfe.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vfe.v90*[1,1], [0, 2], color ='orange'
nfe_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1608.45708, 0.0591, minfe, maxfe)
nfe = alog10(nfe_f.n)
nfe_e = sqrt(nfe_f.e/nfe_f.n/alog(10.)^2)
print, nfe, nfe_e



plot, (flux.wave-1304.37)/1304.37*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
mins2_13 = -430


djs_oplot, [mins2_13, mins2_13], [0, 2], color = 'red'
maxs2_13 = 150

djs_oplot, [maxs2_13, maxs2_13], [0, 2], color = 'red'

fitregs2_13 = where(abs((flux.wave-1304.37)/1304.37*cspeed) lt 1000)
vs2_13 = mcmcv(flux.wave[fitregs2_13], flux.flux[fitregs2_13]/cfit[fitregs2_13], flux.err[fitregs2_13]/cfit[fitregs2_13], 1304.37, maxv=maxs2_13, minv=mins2_13)
djs_oplot, vs2_13.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vs2_13.v90*[1,1], [0, 2], color ='orange'
no1_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1304.37, 0.0928, mins2_13, maxs2_13)
no1 = alog10(no1_f.n)
no1_e = sqrt(no1_f.e/no1_f.n/alog(10.)^2)
print, no1, no1_e

sreg = where(((flux.wave-1808.00)/1800.00*cspeed lt -400 and (flux.wave-1808.00)/1800.00*cspeed gt -1000) or ((flux.wave-1808.00)/1800.00*cspeed lt 400 and (flux.wave-1808.00)/1800.00*cspeed gt 100))
sfit = poly_fit(flux.wave[sreg], flux.flux[sreg]/cfit[sreg], 1, yfit = recont)
recont_full = interpol( recont,flux.wave[sreg], flux.wave) 
plot, (flux.wave-1808.00)/1800.00*cspeed, flux.flux/cfit/recont_full, yr = [0 ,2], xr = [-1000, 1000] ;O I region

djs_oplot, [-1000, 1000], [1,1], color = 'blue'
mins2_18 = -270


djs_oplot, [mins2_18, mins2_18], [0, 2], color = 'red'
maxs2_18 = 70

djs_oplot, [maxs2_18, maxs2_18], [0, 2], color = 'red'

fitregs2_18 = where(abs((flux.wave-1808.00)/1808.00*cspeed) lt 500)
vs2_18 = mcmcv(flux.wave[fitregs2_18], flux.flux[fitregs2_18]/cfit[fitregs2_18]/recont_full[fitregs2_18], flux.err[fitregs2_18]/cfit[fitregs2_18]/recont_full[fitregs2_18], 1808.00, maxv=maxs2_18, minv=mins2_18)
djs_oplot, vs2_18.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vs2_18.v90*[1,1], [0, 2], color ='orange'
ns2_18_f = mcmcn(flux.wave, flux.flux/cfit/recont_full,flux.err/cfit/recont_full, 1.0, 0, 1808.00, 2.49E-3, mins2_18, maxs2_18)
ns2 = alog10(ns2_18_f.n)
ns2_e = sqrt(ns2_18_f.e^2/ns2_18_f.n^2/alog(10.)^2)
print, ns2, ns2_e


plot, (flux.wave-1264.738)/1264.738*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region

djs_oplot, [-1000, 1000], [1,1], color = 'blue'
minsfs = -220


djs_oplot, [minsfs, minsfs], [0, 2], color = 'red'
maxsfs = 300

djs_oplot, [maxsfs, maxsfs], [0, 2], color = 'red'

fitregsfs = where(abs((flux.wave-1264.738)/1264.738*cspeed) lt 500)
vsfs = mcmcv(flux.wave[fitregsfs], flux.flux[fitregsfs]/cfit[fitregsfs], flux.err[fitregsfs]/cfit[fitregsfs], 1264.738, maxv=maxsfs, minv=minsfs)
djs_oplot, vsfs.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vsfs.v90*[1,1], [0, 2], color ='orange'
nsfs_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 0.86, 0.03, 1264.738, 1.09, minsfs, maxsfs)
nsfs = alog10(nsfs_f.n)
nsfs_e = sqrt(nsfs_f.e^2/nsfs_f.n^2/alog(10.)^2)
print, nsfs, nsfs_e


plot, (flux.wave-1260.42)/1260.42*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
mins2_12 = -600


djs_oplot, [mins2_12, mins2_12], [0, 2], color = 'red'
maxs2_12 = 150

djs_oplot, [maxs2_12, maxs2_12], [0, 2], color = 'red'

fitregs2_12 = where(abs((flux.wave-1260.42)/1260.42*cspeed) lt 1000)
vs2_12 = mcmcv(flux.wave[fitregs2_12], flux.flux[fitregs2_12]/cfit[fitregs2_12], flux.err[fitregs2_12]/cfit[fitregs2_12], 1260.42, maxv=maxs2_12, minv=mins2_12)
djs_oplot, vs2_12.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vs2_12.v90*[1,1], [0, 2], color ='orange'
nsi2_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1260.42, 1.22, mins2_12, maxs2_12)
nsi2 = alog10(nsi2_f.n)
nsi2_e = sqrt(nsi2_f.e/nsi2_f.n/alog(10.)^2)
print, nsi2, nsi2_e

plot, (flux.wave-1250.578)/1250.578*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
mins2 = -310


djs_oplot, [mins2, mins2], [0, 2], color = 'red'
maxs2 = 50

djs_oplot, [maxs2, maxs2], [0, 2], color = 'red'

fitregs2 = where(abs((flux.wave-1250.578)/1250.578*cspeed) lt 1000)
vs2 = mcmcv(flux.wave[fitregs2], flux.flux[fitregs2]/cfit[fitregs2], flux.err[fitregs2]/cfit[fitregs2], 1250.578, maxv=maxs2, minv=mins2)
djs_oplot, vs2.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vs2.v90*[1,1], [0, 2], color ='orange'
ns2_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1250.578, 0.00602, mins2, maxs2)
ns2 = alog10(ns2_f.n)
ns2_e = ns2_f.e/ns2_f.n/alog(10.)
print, ns2, ns2_e


plot, (flux.wave-1670.7867)/1670.7867*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
mina2 = -550


djs_oplot, [mina2, mina2], [0, 2], color = 'red'
maxa2 = 75

djs_oplot, [maxa2, maxa2], [0, 2], color = 'red'

fitrega2 = where(abs((flux.wave-1670.7867)/1670.7867*cspeed) lt 1000)
va2 = mcmcv(flux.wave[fitrega2], flux.flux[fitrega2]/cfit[fitrega2], flux.err[fitrega2]/cfit[fitrega2], 1670.7867, maxv=maxa2, minv=mina2)
djs_oplot, va2.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, va2.v90*[1,1], [0, 2], color ='orange'
na2_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1670.7867, 1.77, mina2, maxa2)
na2 = alog10(na2_f.n)
na2_e = na2_f.e/na2_f.n/alog(10.)
print, na2, na2_e

plot, (flux.wave-1854.72)/1854.72*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
mina3_s = -400


djs_oplot, [mina3_s, mina3_s], [0, 2], color = 'red'
maxa3_s = 150

djs_oplot, [maxa3_s, maxa3_s], [0, 2], color = 'red'

fitrega3_s = where(abs((flux.wave-1854.72)/1854.72*cspeed) lt 1000)
va3_S = mcmcv(flux.wave[fitrega3_s], flux.flux[fitrega3_s]/cfit[fitrega3_s], flux.err[fitrega3_s]/cfit[fitrega3_s], 1854.72, maxv=maxa3_s, minv=mina3_s)
djs_oplot, va3_s.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, va3_s.v90*[1,1], [0, 2], color ='orange'
na3_s_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1854.72, 0.56, mina3_s, maxa3_s)
na3_s = alog10(na3_s_f.n)
na3_s_e = na3_s_f.e/na3_s_f.n/alog(10.)
print, na3_s, na3_s_e

plot, (flux.wave-1862.79)/1862.79*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
mina3_w = -400


djs_oplot, [mina3_w, mina3_w], [0, 2], color = 'red'
maxa3_w = 200

djs_oplot, [maxa3_w, maxa3_w], [0, 2], color = 'red'

fitrega3_w = where(abs((flux.wave-1862.79)/1862.79*cspeed) lt 1000)
va3_w = mcmcv(flux.wave[fitrega3_w], flux.flux[fitrega3_w]/cfit[fitrega3_w], flux.err[fitrega3_w]/cfit[fitrega3_w], 1862.79, maxv=maxa3_w, minv=mina3_w)
djs_oplot, va3_w.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, va3_w.v90*[1,1], [0, 2], color ='orange'
na3_w_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1862.79, 0.279, mina3_w, maxa3_w)
na3_w = alog10(na3_w_f.n)
na3_w_e = na3_w_f.e/na3_w_f.n/alog(10.)
print, na3_w, na3_w_e


plot, (flux.wave-1334.532)/1334.532*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
minc2 = -600


djs_oplot, [minc2, minc2], [0, 2], color = 'red'
maxc2 = 300

djs_oplot, [maxc2, maxc2], [0, 2], color = 'red'

fitregc2 = where(abs((flux.wave-1334.532)/1334.532*cspeed) lt 1000)
vc2 = mcmcv(flux.wave[fitregc2], flux.flux[fitregc2]/cfit[fitregc2], flux.err[fitregc2]/cfit[fitregc2], 1334.532, maxv=maxc2, minv=minc2)
djs_oplot, vc2.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vc2.v90*[1,1], [0, 2], color ='orange'
nc2_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1334.532, .129, minc2, maxc2)
nc2_w = alog10(nc2_f.n)
nc2_e = nc2_f.e/nc2_f.n/alog(10.)
print, nc2_w, nc2_e


lyap = s99_continuum_ly(flux.wave, flux.flux/cfit, flux.err/cfit, 130., yfit = lya, noplot= noplot)

plot, (flux.wave-1206.51)/1206.51*cspeed, flux.flux/cfit/median(flux.flux[where(flux.wave gt 1208 and flux.wave lt 1209)]), yr = [0 ,2], xr = [-2000, 2000] ;O I region
djs_oplot, [-2000, 2000], [1,1], color = 'blue'
mins3 = -600


djs_oplot, [mins3, mins3], [0, 2], color = 'red'
maxs3 = 300

djs_oplot, [maxs3, maxs3], [0, 2], color = 'red'

fitregc2 = where(abs((flux.wave-1334.532)/1334.532*cspeed) lt 1000)
vc2 = mcmcv(flux.wave[fitregc2], flux.flux[fitregc2]/cfit[fitregc2], flux.err[fitregc2]/cfit[fitregc2], 1334.532, maxv=maxc2, minv=minc2)
djs_oplot, vc2.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vc2.v90*[1,1], [0, 2], color ='orange'
nc2_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1334.532, .129, minc2, maxc2)
nc2_w = alog10(nc2_f.n)
nc2_e = nc2_f.e/nc2_f.n/alog(10.)
print, nc2_w, nc2_e

plot, (flux.wave-1393.755)/1393.755*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
mins4_s = -650


djs_oplot, [mins4_s, mins4_s], [0, 2], color = 'red'
maxs4_s = 350

djs_oplot, [maxs4_s, maxs4_s], [0, 2], color = 'red'

fitregs4_S = where(abs((flux.wave-1393.755)/1393.755*cspeed) lt 1000)
vs4_s = mcmcv(flux.wave[fitregs4_S], flux.flux[fitregs4_S]/cfit[fitregs4_s], flux.err[fitregs4_s]/cfit[fitregs4_s], 1393.755, maxv=maxs4_S, minv=mins4_S)
djs_oplot, vs4_s.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vs4_s.v90*[1,1], [0, 2], color ='orange'
ns4_S_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1393.755, .513, mins4_s, maxs4_s)
ns4_s_w = alog10(ns4_s_f.n)
ns4_s_e = ns4_s_f.e/ns4_s_f.n/alog(10.)
print, ns4_s_w, ns4_s_e

plot, (flux.wave-1402.77)/1402.77*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
mins4_w = -450


djs_oplot, [mins4_w, mins4_w], [0, 2], color = 'red'
maxs4_w = 75

djs_oplot, [maxs4_w, maxs4_w], [0, 2], color = 'red'

fitregs4_w = where(abs((flux.wave-1402.77)/1402.77*cspeed) lt 1000)
vs4_w = mcmcv(flux.wave[fitregs4_w], flux.flux[fitregs4_w]/cfit[fitregs4_w], flux.err[fitregs4_w]/cfit[fitregs4_w], 1402.77, maxv=maxs4_w, minv=mins4_w)
djs_oplot, vs4_w.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vs4_w.v90*[1,1], [0, 2], color ='orange'
ns4_w_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1402.77, .255, mins4_w, maxs4_w)
ns4_w_w = alog10(ns4_w_f.n)
ns4_w_e = ns4_w_f.e/ns4_w_f.n/alog(10.)
print, ns4_w_w, ns4_w_e


plot, (flux.wave-1548.202)/1548.202*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
minc4_s = -550


djs_oplot, [minc4_s, minc4_s], [0, 2], color = 'red'
maxc4_s = 75

djs_oplot, [maxc4_s, maxc4_s], [0, 2], color = 'red'

fitregc4_s = where(abs((flux.wave-1548.202)/1548.202*cspeed) lt 1000)
vc4_s = mcmcv(flux.wave[fitregc4_s], flux.flux[fitregc4_s]/cfit[fitregc4_s], flux.err[fitregc4_s]/cfit[fitregc4_s], 1548.202, maxv=maxc4_s, minv=minc4_s)
djs_oplot, vc4_s.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vc4_s.v90*[1,1], [0, 2], color ='orange'
nc4_s_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1548.202, .19, minc4_s, maxc4_s)
nc4_s_w = alog10(nc4_s_f.n)
nc4_s_e = nc4_s_f.e/nc4_s_f.n/alog(10.)
print, nc4_s_w, nc4_s_e


plot, (flux.wave-1550.774)/1550.774*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
minc4_w = -400


djs_oplot, [minc4_w, minc4_w], [0, 2], color = 'red'
maxc4_w = 200

djs_oplot, [maxc4_w, maxc4_w], [0, 2], color = 'red'

fitregc4_w = where(abs((flux.wave-1550.774)/1550.774*cspeed) lt 1000)
vc4_w = mcmcv(flux.wave[fitregc4_w], flux.flux[fitregc4_w]/cfit[fitregc4_w], flux.err[fitregc4_w]/cfit[fitregc4_w], 1550.774, maxv=maxc4_w, minv=minc4_w)
djs_oplot, vc4_w.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vc4_w.v90*[1,1], [0, 2], color ='orange'
nc4_w_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1550.774, .0952, minc4_w, maxc4_w)
nc4_w_w = alog10(nc4_w_f.n)
nc4_w_e = nc4_w_f.e/nc4_w_f.n/alog(10.)
print, nc4_w_w, nc4_w_e


plot, (flux.wave-1242.804)/1242.804*cspeed, flux.flux/cfit, yr = [0 ,2], xr = [-1000, 1000] ;O I region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
minn5 = -100


djs_oplot, [minn5, minn5], [0, 2], color = 'red'
maxn5 = 275

djs_oplot, [maxn5, maxn5], [0, 2], color = 'red'

fitregn5 = where(abs((flux.wave-1242.804)/1242.804*cspeed) lt 1000)
vn5 = mcmcv(flux.wave[fitregn5], flux.flux[fitregn5]/cfit[fitregn5], $
 flux.err[fitregn5]/cfit[fitregn5], 1242.8041, maxv=maxn5, minv=minn5)
djs_oplot, vn5.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vn5.v90*[1,1], [0, 2], color ='orange'
nn5_f = mcmcn(flux.wave, flux.flux/cfit,flux.err/cfit, 1.0, 0, 1242.804, .078, minn5, maxn5)
nn5_w = alog10(nn5_f.n)
nn5_e = nn5_f.e/nn5_f.n/alog(10.)
print, nn5_w, nn5_e


plot,(so6.wave-1031.912)/1031.912*cspeed, so6.flux/cfito6, yr = [0 ,2], xr = [-1000, 1000] ;O VI region
djs_oplot, [-1000, 1000], [1,1], color = 'blue'
mino6 = -630



djs_oplot, [mino6, mino6], [0, 2], color = 'red'
maxo6 = 90.

djs_oplot, [maxo6, maxo6], [0, 2], color = 'red'

nfluxo6 = so6.flux/cfito6
fitrego6 = where(abs((so6.wave-1031.912)/1031.912*cspeed) lt 1000)



vo6 = mcmcv(so6.wave[fitrego6], so6.flux[fitrego6]/cfito6[fitrego6], $
so6.err[fitrego6]/cfito6[fitrego6], 1031.912, maxv=90., minv=-630., conterr = .1)
djs_oplot, vo6.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, vo6.v90*[1,1], [0, 2], color ='orange'
no6_f = mcmcn(so6.wave[o6reg[fitrego6]], so6.flux[o6reg[fitrego6]]/cfito6[fitrego6], $
so6.err[o6reg[fitrego6]]/cfito6[fitrego6], 1.0, 0, 1031.912, .133, -630., 90.)
no6_w = alog10(no6_f.n)
no6_e = no6_f.e/no6_f.n/alog(10.)
print, no6_w, no6_e


plot, (wave-1264.738)/1264.738*2.99E5, flux/cfit, xr = [-1000, 500], xtit = 'Si II* 1264 Velocity [km/s]', ytit = 'Normalized Flux', thick = 10, charsize = 1.5, charthick = 5
djs_oplot, [-1000, 1000], [1, 1]
djs_oplot, [-230, -230], [0,2]
djs_oplot, [200, 200], [0, 2]
fitregsi2s = where(abs((wave-1264.738)/1264.738*cspeed) lt 300)
si2fsv = mcmcv(wave[fitregsi2s], flux[fitregsi2s]/cfit[fitregsi2s], err[fitregsi2s]/cfit[fitregsi2s], 1264.738, maxv=200, minv=-230)
djs_oplot, si2fsv.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, si2fsv.v90*[1,1], [0, 2], color ='orange'
nsi2fs_n = mcmcn(wave[fitregsi2s], flux[fitregsi2s]/cfit[fitregsi2s], $
err[fitregsi2s]/cfit[fitregsi2s], 1.0, 0, 1264.738, 1.09, -230., 200.)
print, nsi2fs_n.e/nsi2fs_n.n/alog(10.)


plot, (wave-1306.029)/1306.029*2.99E5, flux/cfit, xr = [-1000, 500], $
xtit = 'O I* 1264 Velocity [km/s]', ytit = 'Normalized Flux', charsize = 1.5, charthick = 5
djs_oplot, [-1000, 1000], [1, 1]
djs_oplot, [-150, -150], [0,2]
djs_oplot, [-40, -40], [0, 2]
fitrego1s = where(abs((wave-1306.029)/1306.029*cspeed) lt 300)
o1fsv = mcmcv(wave[fitrego1s], flux[fitrego1s]/cfit[fitrego1s], err[fitrego1s]/cfit[fitrego1s], 1306.029, maxv=-40, minv=-150)
djs_oplot, o1fsv.ewv*[1,1], [0, 2], color ='orange'
djs_oplot, o1fsv.v90*[1,1], [0, 2], color ='orange'
no1fs_n = mcmcn(wave[fitrego1s], flux[fitrego1s]/cfit[fitrego1s], $
err[fitregsi2s]/cfit[fitregsi2s], 1.0, 0, 1306.029,  5.19e-02, -150., -40.)
print, no1fs_n.e/no1fs_n.n/alog(10.)



fwave = flux.wave
fflux = flux.flux
fullfit = cfit

dfpsplot, '/Users/admin/Documents/megasaura/'+gal+'/o6/O6-comp1.ps', /color, /landscape
!p.font = 0
!p.multi = [0, 3, 2]
plot, (so6.wave[o6reg]-1031.912)/1031.912*cspeed, so6.flux[o6reg]/cfito6,$
  yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
  ytit = 'Normalized Flux', charthick =10, charsize = 2
oreg = where((fwave-1302.1685)/1302.1685*cspeed lt 150)
djs_oplot, (fwave[oreg]-1302.1685)/1302.1685*cspeed, fflux[oreg]/fullfit[oreg], color = 'blue', psym = 10, thick = 10
legend, ['O I 1302', '14 eV'], $
  charsize = 1.8, charthick = 10, box = 0, /left
  djs_oplot,  (so6.wave[o6reg]-1031.912)/1031.912*cspeed, flux.err[o6reg]/cfito6/median(flux.flux[where(flux.wave gt 1050 and flux.wave lt 1055)]), color = 'green', thick = 10 , psym = 10
  plot, (so6.wave[o6reg]-1031.912)/1031.912*cspeed, so6.flux[o6reg]/cfito6,$
    yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
    ytit = 'Normalized Flux', charthick =10, charsize = 2
  djs_oplot, (fwave-1260.4221)/1260.4221*cspeed, fflux/fullfit, color = 'grey', psym = 10, thick = 10
  legend, ['Si II 1260', '16 eV'], $
    charsize = 1.8, charthick = 10, box = 0, /left  
  plot, (so6.wave[o6reg]-1031.912)/1031.912*cspeed, so6.flux[o6reg]/cfito6,$
    yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
    ytit = 'Normalized Flux', charthick =10, charsize = 2
  djs_oplot, (fwave-1670.7874)/1670.7874*cspeed, fflux/fullfit, color = 'purple', psym = 10, thick = 10
  legend, ['Al II 1670', '19 eV'], $
    charsize = 1.8, charthick = 10, box = 0, /left
  plot, (so6.wave[o6reg]-1031.912)/1031.912*cspeed, so6.flux[o6reg]/cfito6,$
    yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
    ytit = 'Normalized Flux', charthick =10, charsize = 2
  djs_oplot, (fwave-1334.5323)/1334.5323*cspeed, fflux/fullfit, color = 'green', psym = 10, thick = 10
  legend, ['C II 1335', '24 eV'], $
    charsize = 1.8, charthick = 10, box = 0, /left       
plot, (so6.wave[o6reg]-1031.912)/1031.912*cspeed, so6.flux[o6reg]/cfito6,$ 
    yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
    ytit = 'Normalized Flux', charthick =10, charsize = 2
djs_oplot, (fwave-1393.755)/1393.755*cspeed, fflux/fullfit, color = 'red', psym = 10, thick = 10
legend, ['Si IV 1393', '45 eV'], $
  charsize = 1.8, charthick = 10, box = 0, /left
  plot, (so6.wave[o6reg]-1031.912)/1031.912*cspeed, so6.flux[o6reg]/cfito6,$
    yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
    ytit = 'Normalized Flux', charthick =10, charsize = 2
  creg = where((fwave-1548.1950)/1548.1950*cspeed lt 100)
  djs_oplot, (fwave[creg]-1548.195)/1548.195*cspeed, fflux[creg]/fullfit[creg], color = 'orange', psym = 10, thick = 10
  legend, ['C IV 1548', '65 eV'], $
    charsize = 1.8, charthick = 10, box = 0, /left
dfpsclose

dfpsplot, '/Users/admin/Documents/megasaura/'+gal+'/si2_1808.ps', /color, /landscape
!p.multi = [0, 1, 2]
  plot, (flux.wave-1808.00)/1808.00*cspeed, flux.flux/cfit/recont_full,$
    yr = [0, 1.5], /xs, xr = [-800, 500], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
    ytit = 'Normalized Flux', charthick =10, charsize = 2
    djs_oplot, [-1000, 1000], [1, 1], thick = 10, linestyle = 2
    djs_oplot, (flux.wave-1808.00)/1808.00*cspeed, flux.err/cfit/recont_full, color = 'green', psym = 10, thick = 10
  plot, (flux.wave-1264.738)/1264.738*cspeed, flux.flux/cfit,$
    yr = [0, 1.5], /xs, xr = [-500, 500], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
    ytit = 'Normalized Flux', charthick =10, charsize = 2
    djs_oplot, [-1000, 1000], [1, 1], thick = 10, linestyle = 2
    djs_oplot, (flux.wave-1264.738)/1264.738*cspeed, flux.err/cfit, color = 'green', psym = 10, thick = 10
!p.multi = 0
dfpsclose

;Want Si II 1334

c2reg = where(abs((flux.wave - 1334.5323)/1334.5323*cspeed) lt 1000)

linterp, (flux.wave[c2reg] - 1334.5323)/1334.5323*cspeed, flux.flux[c2reg]/cfit[c2reg], (so6.wave-1031.912)/1031.912*cspeed, c2flux

;Si II 1260 
s2reg = where(abs((flux.wave - 1260.42)/1260.42*cspeed) lt 1000)

linterp, (flux.wave[s2reg] - 1260.42)/1260.42*cspeed, flux.flux[s2reg]/cfit[s2reg], (so6.wave-1031.912)/1031.912*cspeed, s2flux

;AL II 1670
al2reg = where(abs((flux.wave - 1670.7874)/1670.7874*cspeed) lt 1000)

linterp, (flux.wave[al2reg] - 1670.7874)/1670.7874*cspeed, flux.flux[al2reg]/cfit[al2reg], (so6.wave-1031.912)/1031.912*cspeed, al2flux

;AL III 1862.7895
al3reg = where(abs((flux.wave - 1402.77)/1402.77*cspeed) lt 1000)

linterp, (flux.wave[al3reg] - 1402.77)/1402.77*cspeed, flux.flux[al3reg]/cfit[al3reg], (so6.wave-1031.912)/1031.912*cspeed, al3flux


dfpsplot, '/Users/admin/Documents/megasaura/'+gal+'/o6/O6-rat-comp.ps', /color, /landscape
!p.font = 0
!p.multi = [0, 2, 2]
s2lab = textoidl('O VI / Si II')
v9 = textoidl('v_{90}')
vc = textoidl('v_{cen}')
plot, (so6.wave-1031.912)/1031.912*cspeed, so6.flux/cfito6/s2flux,$
  yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
  ytit = s2lab, charthick =10, charsize = 1.5
  djs_oplot, -532*[1,1], [-1000, 1000], thick = 7, linestyle = 2
   xyouts, -650, 1.7, v9, charsize = 1.3
 
  c2lab = textoidl('O VI / C II')
  
  plot, (so6.wave-1031.912)/1031.912*cspeed, so6.flux/cfito6/c2flux,$
  yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
  ytit = c2lab, charthick =10, charsize = 1.5
  al2lab = textoidl('O VI / Al II')
  djs_oplot, -532*[1,1], [-1000, 1000], thick = 7, linestyle = 2
   xyouts, -650, 1.7, v9, charsize = 1.3

  plot, (so6.wave-1031.912)/1031.912*cspeed, so6.flux/cfito6/al2flux,$
  yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
  ytit = al2lab, charthick =10, charsize = 1.5
  djs_oplot, -532*[1,1], [-1000, 1000], thick = 7, linestyle = 2
  djs_oplot, -246*[1,1], [-1000, 1000], thick = 7, linestyle = 2
   xyouts, -650, 1.7, v9, charsize = 1.3
   xyouts, -375, 1.7, vc, charsize = 1.3
     al3lab = textoidl('O VI / Si IV')
   
  plot, (so6.wave-1031.912)/1031.912*cspeed, so6.flux/cfito6/al3flux,$
  yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
  ytit = al3lab, charthick =10, charsize = 1.5
  djs_oplot, -532*[1,1], [-1000, 1000], thick = 7, linestyle = 2  
  djs_oplot, -246*[1,1], [-1000, 1000], thick = 7, linestyle = 2
   xyouts, -650, 1.7, v9, charsize = 1.3
   xyouts, -375, 1.7, vc, charsize = 1.3
dfpsclose


dfpsplot, '/Users/admin/Documents/megasaura/'+gal+'/o6/O6-nv-comp.ps', /color, /landscape
  !p.multi = 0
  plot, (so6.wave[o6reg]-1031.912)/1031.912*cspeed, so6.flux[o6reg]/cfito6,$
    yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
    ytit = 'Normalized Flux', charthick =10, charsize = 2
  djs_oplot, (fwave[creg]-1242.804)/1242.804*cspeed, fflux[creg]/fullfit[creg], color = 'orange', psym = 10, thick = 10
  legend, ['N V 1243', '98 eV'], $
    charsize = 1.8, charthick = 10, box = 0, /left   
  dfpsclose
  
  
 ;AOD plot
so6.flux[where(so6.flux le 0)] = 0.001
dfpsplot, '/Users/admin/Documents/megasaura/'+gal+'/o6/O6-op-comp.ps', /color, /landscape  
!p.font = 0
ytit = textoidl('\tau_{AOD}')
plot, (so6.wave-1031.912)/1031.912*cspeed, alog(1./(so6.flux/cfito6)),$
    yr = [0, 5], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
    ytit = ytit, charthick =10, charsize = 2
djs_oplot, (so6.wave-1031.912)/1031.912*cspeed, alog(1./al2flux),$
  yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
  ytit = al2lab, charthick =10, charsize = 1.5, color = 'purple'
djs_oplot, (so6.wave-1031.912)/1031.912*cspeed, alog(1./c2flux),$
  yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
  ytit = al2lab, charthick =10, charsize = 1.5, color = 'green'
djs_oplot, (so6.wave-1031.912)/1031.912*cspeed, alog(1./s2flux),$
  yr = [0, 2], /xs, xr = [-900, 300], psym = 10, thick = 10, xtit = 'Velocity (km/s)', $
  ytit = al2lab, charthick =10, charsize = 1.5, color = 'grey'
legend, ['O VI', 'C II', 'Si II', 'Al II'], textcolor = [0, djs_icolor('green'), djs_icolor('grey'), djs_icolor('purple')], thick = 5, box = 0, charsize = 2
dfpsclose

dfpsplot,  '/Users/admin/Documents/megasaura/'+gal+'/c4_met_comp.ps', /color, /landscape
!p.multi = [0, 1, 2]
plot, mod04.wave, gconv(mod02.flux[*,9], 200./146.), xr = [1500, 1600], thick = 10, xtit = 'Wavelength', ytit = 'Normalized Flux', charthick = 10, charsize = 2
djs_oplot, mod001.wave, gconv(mod001.flux[*,9], 200./146.), xr = [1500, 1600], thick = 10, color = 'red'
highmet = textoidl('1 Z_\odot')
lowmet = textoidl('0.05 Z_\odot')
legend, [highmet, lowmet], textcolor = [0, djs_icolor('red')], charsize = 2, charthick = 10,  /bottom, /right, box = 0

plot, mod04.wave, gconv(mod02.flux[*,0], 200./146.), xr = [1500, 1600], thick = 10, xtit = 'Wavelength', yr =[0., 1.2], ytit = 'Normalized Flux', charthick = 10, charsize = 2
djs_oplot, mod001.wave, gconv(mod02.flux[*,9], 200./146.), xr = [1500, 1600], thick = 10, color = 'red'
highmet = textoidl('1 Z_\odot')
lowmet = textoidl('0.05 Z_\odot')
legend, ['1 Myr', '40 Myr'], textcolor = [0, djs_icolor('red')], charsize = 2, charthick = 10,  /top, /right, box = 0
dfpsclose
