function circle, xcenter, ycenter, radius, npoints, ratio
theta = 2 * !pi / double(npoints) * dindgen(npoints)
x = xcenter + radius * cos(theta) ; ratio is rmaj/rmin
y = ycenter + radius/ratio * sin(theta)
return, transpose([[x],[y]])
end

pro sps
;input files and parameters
file_name = 'sps.fits'
color_table = [70, 70]
distance = 2000d
nannuli = 15
psym = 15
step = 0.3 ; the interval of the sampled points on the annualli
;output files 
fftfig = file_basename(file_name, '.fits')+'_fft.eps'
spsfig = file_basename(file_name, '.fits')+'_sps.eps'

;calculate 2-d fft of the input field; plot the original map and the modulous distribution.
cgloadct, color_table[0], /reverse
p1 = [0.1, 0.1, 0.45, 0.9]
p2 = [0.55, 0.1, 0.9, 0.9]
fits_read, file_name, dat, hdr
s1 = max(size(dat, /dim))
s = size(dat, /dim)
z = dat
z = (z - mean(z[where(finite(z))]))/stdev(z[where(finite(z))])
z[where(finite(z, /nan))] = 0
;z = SMOOTH(z, 2, /edge_mirror)
ratio = float(s[0])/s[1]
cgps_open, fftfig, xsize = 18*ratio, ysize = 11, /encapsulated
!P.charsize = 1.5*ratio
cgimage, bytscl(z, min = min(z), max = max(z)*0.8), /keep_aspect_ratio, position = p1
; compute the two-dimensional fft
s = size(z, /dim)
help, z
f = fft(z, /center)
f = f;*double(s[0]*s[1])
linpower = abs(f)^2
logpower = alog10(abs(f)^2) ; log of fourier power spectrum
print, min(logpower), max(logpower)
cgplot, [0], [0], /nodata, /noerase,  xtitle = 'x (pixel)', ytitle = 'y (pixel)', $
        xrange = [0, s[0]-1], yrange = [0, s[1]-1], position = p1
cgcolorbar, position = [0.83, 0.1, 0.88, 0.45], /top, title = 'velocity [km s!u-1!n]', range = [min(z), max(z)*0.8]
cgloadct, color_table[1], /reverse
cgimage, bytscl(logpower, min = min(logpower), max = max(logpower)), /keep_aspect_ratio, position = p2, /noerase
cgplot, [0], [0], /nodata, /noerase, xrange = [0-s[0]/2d, s[0]-1-s[0]/2d], $
        yrange = [0-s[1]/2d, s[1]-1-s[1]/2d], xtitle = 'u ( k!du!n )', ytitle = 'v ( k!dv!n )', position=p2
cgcolorbar, position = [0.83, 0.55, 0.88, 0.9], /top, title = 'log!d10!n power', range = [min(logpower), max(logpower)]
;plot annuli on the 2-d power density map 
mx = max(logpower, location)
center = array_indices(logpower, location)
spatial_ex = distance/206265d * 30d * sqrt(((max(s))^2/ratio))
;print, spatial_ex, ratio
krange = [1d, max(s)/2 - 1]
lrange = spatial_ex/krange
;divide the spatial extent into evenly distributed steps in log space. 
logstep = alog10(lrange[0]) - dindgen(nannuli)*(alog10(lrange[0])-alog10(lrange[1]))/(nannuli-1)
linstep = 10d^logstep
kstep = spatial_ex/linstep
for i = 0, n_elements(kstep)-1 do begin 
;tvellipse, kstep[i], kstep[i]/ratio, 0, 0, /data, npoints = 500, color = cgcolor('black')
radius = kstep[i]
npoints = fix(2*!pi*radius/step) 
results = circle(center[0], center[1], radius, npoints, ratio)
cgplot, results[0, *] - s[0]/2, results[1, *] - s[1]/2, /overplot
endfor
cgps_close

;make interpolation according to the coordinates of the annuli in the logpower array.
de = dblarr(nannuli)
mad = dblarr(nannuli)
;length = dblarr(nannuli)
cgps_open, 'annuli_stat.eps', xsize=18, ysize=8, /encapsulated
!p.multi = [0, 5, 3]
!p.charsize = 2
for i = 0, nannuli-1 do begin 
radius = kstep[i]
npoints = fix(2*!pi*radius/step) 
results = circle(center[0], center[1], radius, npoints, ratio)
sample = interpolate(linpower, results[0,*], results[1, *], missing = 0)
good = where(sample gt 0)
de[i] = median(sample[good])
mad[i] = median(sqrt((sample[good]-de[i])^2))
cghistoplot, alog10(sample[good]), yminor = 10, /fillpolygon, polycolor = 'rose', $
            /frequency, ytickformat = '(f7.2)', xtitle = 'log!d10!n power'
cgplot, [alog10(median(sample[good])), alog10(median(sample[good]))], [0, 1], $
        color = 'forest green', /overplot, thick = 4
cgplot, [alog10(mean(sample[good])), alog10(mean(sample[good]))], [0,1], $
        color = 'dodger blue', /overplot, thick = 4
endfor
cgps_close


;----------------------------------ellipse average------------------------------
cgps_open, spsfig, xsize = 10, ysize = 7, /encapsulated
;!p.font = -1
!p.charsize = 2
;fit the spectral index
measure_errors = mad
expar = 'p[0]*x^p[1]'
start = [1.9d10, -3d]
fit = mpfitexpr(expar, kstep, de, measure_errors, start, perror=perror, $
      bestnorm = bestnorm, dof= dof, quiet = 1)
;pcerror = perror * sqrt(bestnorm / dof)
yfit = fit[0]*kstep^fit[1]
cgplot, linstep, de, /xlog, /ylog, psym = psym, color = 'red', position = [0.2, 0.2, 0.8, 0.8], $
        xrange = [max(linstep)*3, min(linstep)*0.5], yrange = [min(de)*0.2, 5*max(de)], yminor = 9,$
        xtitle = 'Linear scale (pc)', ytitle = 'Power', title = 'Power Spectrum',$
        err_yhigh = mad, err_ylow = mad, err_color = 'blu5', err_thick=5
cgplot, linstep, yfit, /overplot
legend = 'E(k)'+textoidl('\propto')+'k!u'+strtrim(string(fit[1], format = '(f7.2)'),2)+'!n'+ $
        textoidl('!u\pm')+strtrim(string(perror[1],format='(f4.2)'),2)+'!n'
al_legend, legend, /right, background_color = 'lavender'
cgps_close
end
