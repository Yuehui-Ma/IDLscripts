pro plot_slice 
;--------------------------------------------------------
;------------pv cut setting------------------------------
;--------------------------------------------------------
inputfile_12 = './datacube/cas_U.fits'
inputfile_13 = './datacube/cas_L.fits'
outputfile_12 = 'knot1_12.fits'
outputfile_13 = 'knot1_13.fits'
picname = 'knot1_slice.eps'
slice = 'a'
path = [[111.742, -2.23334], [111.742, -2.03334]] ; slice a
step = 0
width = 2.
;--------------------------------------------------------
;------------pv plot setting-----------------------------
;--------------------------------------------------------
picname = picname 
prange = 0
vrange = [-52, -32]
start_level = 0.6/sqrt(width); for 13CO
nlevel = 5. ;for 13CO
min_bg = 1.0/sqrt(width); start level for 12CO 2sigma 
;-------------------------------------------------------
cgLoadCT, 49, NColors=12, Bottom=3
pos = [0.2, 0.2, 0.8, 0.8]
;--------------------------------------------------------
cgps_open, picname, xsize=6, ysize=4, /portrait , /encapsulated
    !p.CHARSIZE=1
    !P.Thick = 1
    !P.CharThick = 1.5
    !X.Thick = 1
    !Y.Thick = 1
    !Z.Thick = 1
    !Y.Ticklen = 0.04
    !P.Font = -1
;--------------------------------------------------------
;------------extract pvslice or belt---------------------
;--------------------------------------------------------
fits_read, inputfile_12, dat12, hdr
fits_read, inputfile_13, dat13, hdr13
sxaddpar, hdr, 'CTYPE1', repstr(sxpar(hdr,'CTYPE1'), 'GLS', 'SFL')
sxaddpar, hdr, 'CTYPE2', repstr(sxpar(hdr,'CTYPE2'), 'GLS', 'SFL')
adxy, hdr, path[0,*], path[1,*], x, y
defstep = 0.5 
if (step eq 0) then step = defstep
length = sqrt( ((x-shift(x,1))[1:*])^2+((y-shift(y,1))[1:*])^2 )
nstep = ceil(total(length)/step)
step = total(length)/nstep
print, 'Use step ' + string(step) + ' * pixelsize'
xs = dindgen(nstep+1)
ys = dindgen(nstep+1)
intlen = dindgen(n_elements(length))
for i=0, n_elements(length)-1 do intlen[i]=total(length[0:i])
for i=0, nstep do begin
    way = i*step
    node = (where(way le intlen))[0]
    xs[i] = x[node+1]-(intlen[node]-way)/length[node]*(x[node+1]-x[node])
    ys[i] = y[node+1]-(intlen[node]-way)/length[node]*(y[node+1]-y[node])
endfor
if width ne 0 then begin 
    xg = fltarr(nstep+1,width)
    yg = fltarr(nstep+1,width)
    width = fix(width)
    temp = sqrt((xs[0]-xs[1])^2+(ys[0]-ys[1])^2)
    xg[0,*] = xs[0]+(findgen(width)-(width-1)/2.)/temp*(ys[1]-ys[0])
    yg[0,*] = ys[0]-(findgen(width)-(width-1)/2.)/temp*(xs[1]-xs[0])
    for i=1, nstep do begin
        temp = sqrt((xs[i]-xs[i-1])^2+(ys[i]-ys[i-1])^2)
        xg[i,*] = xs[i]+(findgen(width)-(width-1)/2.)/temp*(ys[i]-ys[i-1])
        yg[i,*] = ys[i]-(findgen(width)-(width-1)/2.)/temp*(xs[i]-xs[i-1])
    endfor
endif 
;-------------------for 12CO--------------------------------------------
nx1_12 = sxpar(hdr, 'NAXIS3')
nx2_12 = n_elements(xs)
nx1_13 = sxpar(hdr13, 'NAXIS3')
nx2_13 = n_elements(xs)
slice12 = make_array(nx1_12, nx2_12, type=size(dat12,/type))
slice13 = make_array(nx1_13, nx2_13, type=size(dat13,/type))
if width eq 0 then begin 
    for i=0, nx1_12-1 do slice12[i,*] = interpolate(dat12[*,*,i], xs, ys, missing=0)
    for j=0, nx1_13-1 do slice13[j,*] = interpolate(dat13[*,*,j], xs, ys, missing=0)
endif else begin
    for i=0, nx1_12-1 do begin
        g = interpolate(dat12[*,*,i], xg, yg, missing=0)
        slice12[i,*] = total(g,2)/width
    endfor 
    for j=0, nx1_13-1 do begin 
        g = interpolate(dat13[*,*,j], xg, yg, missing=0)
        slice13[j,*] = total(g,2)/width
    endfor
endelse 
mkhdr, pvhdr12, slice12
factor = abs(sxpar(hdr,'CDELT3') lt 100)?(1d):(1000d)
sxaddpar, pvhdr12, 'CTYPE1', 'VELOCITY'
sxaddpar, pvhdr12, 'CRPIX1', sxpar(hdr,'CRPIX3')
sxaddpar, pvhdr12, 'CRVAL1', sxpar(hdr,'CRVAL3')/factor
sxaddpar, pvhdr12, 'CDELT1', sxpar(hdr,'CDELT3')/factor
sxaddpar, pvhdr12, 'CTYPE2', 'POSITION'
sxaddpar, pvhdr12, 'CRPIX2', 1
sxaddpar, pvhdr12, 'CRVAL2', 0d
sxaddpar, pvhdr12, 'CDELT2', step*abs(sxpar(hdr,'CDELT1'))
sxaddhist, 'PV file: ' + inputfile_12, pvhdr12
sxaddhist, 'PV path:', pvhdr12
for i=0, n_elements(path[0,*])-1 do sxaddhist, string(path[0,i])+' '+string(path[1,i]), pvhdr12
sxaddhist, 'Position in Degree', pvhdr12
sxaddhist, 'Velocity in km/s', pvhdr12
fits_write, outputfile_12, slice12, pvhdr12
;-------------------------for 13CO------------------------------------
mkhdr, pvhdr13, slice13
factor = abs(sxpar(hdr13,'CDELT3') lt 100)?(1d):(1000d)
sxaddpar, pvhdr13, 'CTYPE1', 'VELOCITY'
sxaddpar, pvhdr13, 'CRPIX1', sxpar(hdr13,'CRPIX3')
sxaddpar, pvhdr13, 'CRVAL1', sxpar(hdr13,'CRVAL3')/factor
sxaddpar, pvhdr13, 'CDELT1', sxpar(hdr13,'CDELT3')/factor
sxaddpar, pvhdr13, 'CTYPE2', 'POSITION'
sxaddpar, pvhdr13, 'CRPIX2', 1
sxaddpar, pvhdr13, 'CRVAL2', 0d
sxaddpar, pvhdr13, 'CDELT2', step*abs(sxpar(hdr,'CDELT1'))
sxaddhist, 'PV file: ' + inputfile_13, pvhdr13
sxaddhist, 'PV path:', pvhdr13
for i=0, n_elements(path[0,*])-1 do sxaddhist, string(path[0,i])+' '+string(path[1,i]), pvhdr13
sxaddhist, 'Position in Degree', pvhdr13
sxaddhist, 'Velocity in km/s', pvhdr13
fits_write, outputfile_13, slice13, pvhdr13
;----------------record the path---------------------------------------
xyad, hdr, xs, ys, as, ds
posi = sxpar(pvhdr12,'CRVAL2') + (dindgen(sxpar(pvhdr12,'NAXIS2')) - sxpar(pvhdr12,'CRPIX2')+1)*sxpar(pvhdr12,'CDELT2')
openw, lun, file_basename(picname, '.eps')+'.track', /get_lun
printf, lun, string(lindgen(n_elements(as))+1)+' '+string(as)+' '+string(ds)+' '+string(posi);+' '+string(mmnt)
close,lun
free_lun,lun
;--------------------------------------------------------------------
;----------------make the plot---------------------------------------
;--------------------------------------------------------------------
fits_read, outputfile_12, slice_12, hdr12
fits_read, outputfile_13, slice_13, hdr13
nc_v12 = sxpar(hdr12,'NAXIS1')
cv_v12 = sxpar(hdr12,'CRVAL1')
cp_v12 = sxpar(hdr12,'CRPIX1')
cd_v12 = sxpar(hdr12,'CDELT1')
v12 = ((findgen(nc_v12)-cp_v12+1)*cd_v12+cv_v12)
c12 = (vrange-cv_v12)/cd_v12+cp_v12-1
;--------------------------------------------------------
nc_v13 = sxpar(hdr13,'NAXIS1')
cv_v13 = sxpar(hdr13,'CRVAL1')
cp_v13 = sxpar(hdr13,'CRPIX1')
cd_v13 = sxpar(hdr13,'CDELT1')
v13 = ((findgen(nc_v13)-cp_v13+1)*cd_v13+cv_v13)
c13 = (vrange-cv_v13)/cd_v13+cp_v13-1
;---------------------------------------------------------
nc_p12 = sxpar(hdr12,'NAXIS2')
cv_p12 = sxpar(hdr12,'CRVAL2')
cp_p12 = sxpar(hdr12,'CRPIX2')
cd_p12 = sxpar(hdr12,'CDELT2')
p12 = ((findgen(nc_p12)-cp_p12+1)*cd_p12+cv_p12)*60.
if prange eq 0 then cp12 = [0,nc_p12-1] else cp12 = (prange-cv_p12)/cd_p12+cp_p12-1
;---------------------------------------------------------
nc_p13 = sxpar(hdr13,'NAXIS2')
cv_p13 = sxpar(hdr13,'CRVAL2')
cp_p13 = sxpar(hdr13,'CRPIX2')
cd_p13 = sxpar(hdr13,'CDELT2')
p13 = ((findgen(nc_p13)-cp_p13+1)*cd_p13+cv_p13)
if prange eq 0 then cp13 = [0,nc_p13-1] else cp13 = (prange-cv_p13)/cd_p13+cp_p13-1
;----------------------------------------------------------
if prange eq 0 then yrange = [min(p12),max(p12)] else yrange = prange
img12 = bytscl(smooth(slice_12[c12[0]:c12[1],cp12[0]:cp12[1]],2), min = min_bg, max = 0.8* max(smooth(slice_12[c12[0]:c12[1],cp12[0]:cp12[1]],2)))
con = smooth(slice_13[c13[0]:c13[1],cp13[0]:cp13[1]], 2)
levels = 10
step = 0.75
userLevels = IndGen(levels) * step + 0.5
print, userlevels
SetDecomposedState, 0, CurrentState=state
cgplot, [0] , [0], position=pos, /nodata, /noerase, XSTYLE=1, YSTYLE=1, XRANGE=vrange, yrange=yrange, $
        xtitle='!17 Velocity (km/s)', ytitle='!17 Position (Arcmin)'
cgContour, smooth(slice_12[c12[0]:c12[1],cp12[0]:cp12[1]], 2), /Fill, C_Colors=Indgen(levels)+5, $
           levels=userlevels, Position = pos, /onimage, /cell_fill
max_con = max(smooth(slice_13[c13[0]:c13[1],cp13[0]:cp13[1]],2))
st = double(max_con*0.8 - start_level)/(nlevel - 1) 
level = start_level + indgen(nlevel)*st
cgcontour, con, /onimage, levels=level, label=0, position=pos, color='indian red', thick = 1.5;, c_linestyle =3
cgtext, 0.75, 0.75, slice, /normal, charsize = 1.5
cgtext, 0.56, 0.23, '!17 Contour: !U13!NCO', /normal, color = 'indian red'
cgplot, [-15, -7], [0.089, 0.089], linestyle = 1, /overplot
cgplot, [-15, -7], [0.183, 0.183], linestyle = 1, /overplot
cgplot, [0], [0], position=pos, /nodata, /noerase, XSTYLE=1, YSTYLE=1, $
        XRANGE=vrange, yrange=yrange, xtickformat = '(A1)', ytickformat = '(A1)'
cgColorBar, NColors=12, Bottom=3, Divisions=6, Range=[userlevels[0], userlevels[n_elements(userlevels)-1]], $
        Format='(f4.1)', Position = [0.2, 0.82, 0.8, 0.86], AnnotateColor='black', title = 'K', /top

center = [111.74198d, -2.1333435d]
cocenter = [111.742, -2.12917]

;slice a:
fsd = 153./3600.
rsd = 95./3600.
fshock_loc = [(center[1] - path[1,0] - fsd)* 60., (center[1] - path[1,0] + fsd)* 60.]
rshock_loc = [(center[1] - path[1,0] - rsd)* 60., (center[1] - path[1,0] + rsd)* 60.]
location = sqrt((cocenter[0]-path[0,0])^2+(cocenter[1]-path[1,0])^2) * 60. ;slice b-e
cent_loc = sqrt((center[0]-path[0,0])^2+(center[1]-path[1,0])^2) * 60.
cgplot, vrange, [location, location], /overplot, color = 'red', linestyle=2, thick=3
cgplot, vrange, [cent_loc, cent_loc], /overplot, color = 'indian red', linestyle=2, thick=3
cgplot, vrange, [rshock_loc[0], rshock_loc[0]], /overplot, color = 'purple', linestyle=2, thick=3
cgplot, vrange, [fshock_loc[0], fshock_loc[0]], /overplot, color = 'blue', linestyle=2, thick=3
cgplot, vrange, [rshock_loc[1], rshock_loc[1]], /overplot, color = 'purple', linestyle=2, thick=3
cgplot, vrange, [fshock_loc[1], fshock_loc[1]], /overplot, color = 'blue', linestyle=2, thick=3
items = ['Forward shock', 'Reverse shock', 'CO knot center', 'Expansion center']
al_legend, items, linestyle=2, colors = [ 'blue', 'purple', 'red', 'indian red'], linsize = 0.7, thick=3, /clear,charsize = 0.7;,/bottom;,/right

cgps_close 
end 