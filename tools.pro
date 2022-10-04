function get_coord, hdr
  if ~keyword_set(hdr) then begin
     print, 'Syntax - get_coord, hdr'
     return, ''
  endif
xr = [0, sxpar(hdr, 'NAXIS1')-1]
yr = [0, sxpar(hdr, 'NAXIS2')-1]
xyad, hdr, xr, yr, lr, br 
return, [lr, br]
end





pro plt, file
if ~keyword_set(file) then begin
     print, 'Syntax - plt, file'
     return
endif
fits_read, file, data, hdr
data[where(finite(data, /nan))] = 0
cgloadct, 72, /reverse
stretch, 0, 255, 0.5
img = bytscl(data, min = min(data[where(data ne 0)]), max = max(data[where(data ne 0)]))
pos = [0.2, 0.2, 0.9, 0.9]
cgimage, img, /keep_aspect_ratio, position = pos
r = get_coord(hdr)
xr = r[0:1]
yr = r[2:3]
title = strtrim(gettok(sxpar(hdr,'CTYPE1'),'-'),2)
flag = strcmp(title, 'R', 1, /fold_case)
if flag eq 0 then begin 
    xtitle = 'Galactic Longitude (!Uo!N)'
    ytitle = 'Galactic Latitude (!Uo!N)'
endif else begin 
    xtitle = 'Right Ascension (J2000)'
    ytitle = 'Declination (J2000)'
endelse 
cgplot, 0, 0, position = pos, /noerase, xrange = xr, yrange = yr, $
         xtitle = xtitle, ytitle = ytitle;, xtickinterval = 1, ytickinterval = 1
cgcolorbar, range = [min(data[where(data ne 0)]), max(data[where(data ne 0)])], $
       /vertical, position = [0.9, 0.2, 0.92, 0.9], /right
end






pro plt_rgb, redfile, greenfile, bluefile
if n_params() lt 2  then begin
     print, 'Syntax - plt_rgb, redfile, greenfile, bluefile'
     return
endif
fits_read, redfile, red
red[where(finite(red, /nan))] = 0
fits_read, greenfile, green
green[where(finite(green, /nan))] = 0
fits_read, bluefile, blue, hdr
blue[where(finite(blue, /nan))] = 0
r = get_coord(hdr)
xr = r[0:1]
yr = r[2:3]
title = strtrim(gettok(sxpar(hdr,'CTYPE1'),'-'),2)
flag = strcmp(title, 'R', 1, /fold_case)
if flag eq 0 then begin 
    xtitle = 'Galactic Longitude (!Uo!N)'
    ytitle = 'Galactic Latitude (!Uo!N)'
endif else begin 
    xtitle = 'Right Ascension (J2000)'
    ytitle = 'Declination (J2000)'
endelse 
s = size(green, /dimension)
red = bytscl(smooth(sqrt(red), 3), min=sqrt(min(red[where(red ne 0)])), max=sqrt(10*median(red[where(red ne 0)])))
green = bytscl(smooth(sqrt(green), 3), min=sqrt(min(green[where(green ne 0)])), max=sqrt(10*median(green[where(green ne 0)])))
blue = bytscl(smooth(sqrt(blue), 3), min=sqrt(min(blue[where(blue ne 0)])), max=sqrt(10*median(blue[where(blue ne 0)])))
img = bytarr(3, s[0], s[1])
img[0, *, *]=red
img[1, *, *]=green
img[2, *, *]=blue
pos = [0.2, 0.2, 0.9, 0.9]
cgimage, img, /keep_aspect_ratio, position = pos, Stretch='Clip', Clip=0.6
cgplot, 0, 0, position = pos, /noerase, xrange = xr, yrange = yr,$
       xtitle = xtitle, ytitle = ytitle; ,xtickinterval = 1, ytickinterval = 1
cgplot, 0, 0, position = pos, /noerase, xrange = xr, yrange = yr,$
       axiscolor = 'gray', xtickformat = '(A1)', ytickformat = '(A1)'
end 





pro plt_pv, file
if ~keyword_set(file) then begin
     print, 'Syntax - plt_pv, file'
     return
endif
fits_read, file, data, hdr
data[where(finite(data, /nan))] = 0
cgloadct, 72, /reverse
stretch, 0, 255, 0.5
img = bytscl(data, min = min(data[where(data ne 0)]), max = max(data[where(data ne 0)]))
pos = [0.2, 0.2, 0.9, 0.9]
cgimage, img, /keep_aspect_ratio, position = pos
axis1 = sxpar(hdr,'CRVAL1')+(dindgen(sxpar(hdr,'NAXIS1'))-sxpar(hdr,'CRPIX1')+1)*sxpar(hdr,'CDELT1')
axis2 = sxpar(hdr,'CRVAL2')+(dindgen(sxpar(hdr,'NAXIS2'))-sxpar(hdr,'CRPIX2')+1)*sxpar(hdr,'CDELT2')
title1 = strtrim(gettok(sxpar(hdr,'CTYPE1'),'-'),2)
title2 = strtrim(gettok(sxpar(hdr,'CTYPE2'),'-'),2)
;find velocity axis
axisv = where(strcmp([title1,title2], 'v', 1, /fold_case))
title = strtrim(gettok(sxpar(hdr,'CTYPE1'),'-'),2)
flag = where(strcmp([title1,title2], 'v', 1, /fold_case))
if (abs(sxpar(hdr,'CDELT1')) gt 100) then axis1 = axis1/1000d 
if (abs(sxpar(hdr,'CDELT2')) gt 100) then axis2 = axis2/1000d
if flag eq 0 then begin 
    xtitle = 'Velocity (km s!U-1!N)'
    ytitle = 'Position (!Uo!N)'
endif else begin 
    xtitle = 'Position (!Uo!N)'
    ytitle = 'Velocity (km s!U-1!N)'
endelse 
cgplot, 0, 0, position = pos, /noerase, xrange = [axis1[0], axis1[-1]], $
       yrange = [min(axis2), max(axis2)], xtitle = xtitle, ytitle = ytitle
cgcolorbar, range = [min(data[where(data ne 0)]), max(data[where(data ne 0)])], $
      /vertical, position = [0.9, 0.2, 0.92, 0.9], /right
end 






pro im3d, file 
if ~keyword_set(file) then begin
     print, 'Syntax - im3d, file'
     return
endif
fits_read, file, data, hdr
newdata = transpose(data, [0, 2, 1])
img = bytscl(newdata, min = 1.5, max = max(data)*0.5)
xvolume, img, /INTERPOLATE, RENDERER=1
end



function stdrange, r, rmax
;standardize the range
r = r[sort(r)]
r = [round(r[0]),round(r[1])]
if (r[0] gt rmax) or (r[1] lt 0) then begin
    r=-1
    return,r
endif
r = r >0 <rmax
return,r
end


pro tex, file, vrange, rmsfile, limit
if ~keyword_set(file) then begin
     print, 'Syntax - tex, file, vrange, rmsfile, limit'
     return
endif 
fits_read, file, data, hdr
vrange = vrange * 1000d   ;km/s to m/s
    nc = sxpar(hdr,'NAXIS3')
    cv = sxpar(hdr,'CRVAL3')
    cp = sxpar(hdr,'CRPIX3')
    cd = sxpar(hdr,'CDELT3')
    v = ((findgen(nc)-cp+1)*cd+cv)
    c = (vrange-cv)/cd+cp-1
    c=stdrange(c,sxpar(hdr,'NAXIS3')-1)
 if c[0] eq -1 then begin
        print,'Error: Range is out of the given datacube.'
        return
 endif
data = data[*,*,c[0]:c[1]]
v = v[c[0]:c[1]]
s = size(data)
print,s[1],s[2]
peak = dblarr(s[1],s[2])
fits_read, rmsfile, rms
for i=0,s[1]-1 do begin
  for j=0,s[2]-1 do begin
    peak[i,j] = max(data[i,j,*])
    if (peak[i,j] gt 0) then begin 
    zz = where(data[i,j,*] eq peak[i,j])
   if ((data[i,j,zz-2] gt limit*rms[i,j]) and (data[i,j,zz-1] gt limit*rms[i,j]) and (data[i,j,zz] gt limit*rms[i,j]) and $
      (data[i,j,zz+1] gt limit*rms[i,j]) and (data[i,j,zz+2] gt limit*rms[i,j])) then begin 
    peak[i,j]=peak[i,j] 
    endif else begin
    peak[i,j]=0
    endelse
    endif 
endfor
endfor
prompt = 'at least 5 contiguous channels GT 2sigma in the velocity range [-60,-25]'
sxaddhist, prompt, hdr

T0 = 5.53213817d
JTbg = (exp(T0/2.7)-1)^(-1)
Tex = T0*(alog(1 + (peak/T0 + JTbg)^(-1)))^(-1)
Tex[where(Tex eq 2.7)] = 0
sxdelpar, hdr, ['NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA']+'3'
sxdelpar, hdr, ['NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA']+'4'
fits_write, file_basename(file,'.fits')+'_tex.fits', tex, hdr
end



pro tau13, file13CO, texfile, rmsfile13CO, vrange
if ~keyword_set(file13CO) then begin
     print, 'Syntax - tau13, file13CO, texfile, rmsfile, vrange'
     return
endif 
fits_read, file13CO, dat13, hdr13
nc13 = sxpar(hdr13, 'NAXIS3')
cv13 = sxpar(hdr13, 'CRVAL3')
cp13 = sxpar(hdr13, 'CRPIX3')
cd13 = sxpar(hdr13, 'CDELT3')
if (abs(cd13) gt 1) then begin 
    cd13 = cd13/1000d
    cv13 = cv13/1000d
endif 
v13 = ((findgen(nc13) - cp13 + 1)*cd13 + cv13)
if keyword_set(vrange) then begin
    c13 = (vrange - cv13)/cd13 + cp13 - 1
    c13 = stdrange(c13, nc13-1)  
endif else begin 
    c13 = [0, nc13-1]
endelse
if c13[0] eq -1 then begin
    print, 'Error: vrange is out of the given 13CO datacube.'
    return
endif
dat13 = dat13[*, *, c13[0]:c13[1]]
v13 = v13[c13[0]:c13[1]]
s13 = size(dat13)
peak = dblarr(s13[1], s13[2])
fits_read, rmsfile13CO, rms13
for i=0, s13[1]-1 do begin
    for j=0, s13[2]-1 do begin
        peak[i,j] = max(dat13[i,j,*])
        if (peak[i,j] gt 0) then begin 
            zz = where(dat13[i,j,*] eq peak[i,j])
;            if ((dat13[i,j,zz-1] gt 3*rms13[i,j]) and (dat13[i,j,zz] gt 3*rms13[i,j]) and $
;                (dat13[i,j,zz+1] gt 3*rms13[i,j])) then begin                
            if (dat13[i,j,zz] gt 4*rms13[i,j]) then begin 
                peak[i,j] = peak[i,j] 
            endif else begin
                peak[i,j] = 0
            endelse
        endif 
     endfor
endfor
fits_read, texfile, tex
prompt = 'peak channel above 4 sigma in ' $
        + strtrim(string(v13[0], format = '(f11.2)'), 2) $
        +' to '+ strtrim(string(v13[-1], format = '(f11.2)'), 2) + ' km/s.'
sxaddhist, prompt, hdr13 
J13 = (exp(5.29d/Tex) - 1)^(-1)
tau13 = -alog(1 - peak/(5.29*(J13 - 0.164)))
tau13[where(peak eq 0 or tau13 lt 0)] = 0
tau13[where(finite(tau13, /nan))] = 0
sxdelpar, hdr13, ['NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA']+'3'
sxdelpar, hdr13, ['NAXIS', 'CTYPE', 'CRVAL', 'CDELT', 'CRPIX', 'CROTA']+'4'
fits_write, file_basename(file13CO,'.fits')+'_tau.fits', tau13, hdr13
end 
    



pro cal_n, m0file, texfile, taufile
if ~keyword_set(m0file) then begin
     print, 'Syntax - cal_n, m0file, texfile, taufile'
     return
endif 
fits_read, m0file, m013, hdr
fits_read, texfile, tex
fits_read, taufile, tau13
N13 = 2.42d14*7e5*m013/(1-exp(-5.29d/tex))*(1 + 0.88/tex) * tau13/(1-exp(-tau13)) 
;massLTE = total(N13)*((30d/3600*!dtor)*(dis*3.08568025e18))^2*(2.8*1.6606e-24)/(1.99e33)
;print, strtrim(string(massLTE, format = '(f11.2)'), 2) + ' solar mass derived from 13CO emission.'
fits_write, file_basename(m0, 'm0.fits')+'_N.fits', N13, hdr
end
