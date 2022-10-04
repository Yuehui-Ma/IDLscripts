pro channel_plot ;the output name of the plot is the input fitsname plus '_cplot.ps'
;----------------------------------------------------------------
;---------------------settings-----------------------------------
;----------------------------------------------------------------
fitsname = './mosaic/ngc7538U.fits' 
v = -60.                     ;the velocity to start  
interval = 0.5                 ;the interval of the integration  
xgrid = 5.                    ;the number of pannels in each row 
ygrid = 6.                    ;the number of pannels in each column 
trange = [2*0.5*sqrt(interval*0.16), 15]            ;the intensity scale in each panel 
center = [111.87558, 0.69368]    ;if you need to plot a sub-region instead of a whole fits, this is the center of the sub region in units of deg
psize = [0.5d, 0.5d]          ;the box size of the sub-region, you need to set this parameter to 0 if you want to plot the whole fits  
pos = [0.1, 0.1, 0.9, 0.9]    ;the positon of the entire plot frame 
textp = [0.1, 0.1]           ;the label position of the integrated velocity range in each pannel, in normal unit
textcolor = 'black'           ;the color of the label 
; you may need to modify the parameter 'xticks' in the command line 'cgcontour', to adjust the coordinates you want to show.
; by default it is set to 2 which will lead to 3 ticks in xaxis.

cgloadct, 72, /reverse, clip=[0,240]
stretch, 0, 255, 0.7 
;------------------------------------------------------------------------------------------
;------------caculate the xrange and yrange of the plot region in datacube-----------------
;------------------------------------------------------------------------------------------
fits_read, fitsname, data, hdr 
if (n_elements(psize) eq 2) then begin 
    l = [center[0] + psize[0]/2d, center[0] - psize[0]/2d]
    b = [center[1] - psize[1]/2d, center[1] + psize[1]/2d]
    adxy, hdr, l, b, cx, cy
    data = data[cx[0]:cx[1],cy[0]:cy[1],*]
    s = size(data, /dim)
    nx = s[0]
    ny = s[1]
endif else begin 
    nx = sxpar(hdr, 'NAXIS1')
    ny = sxpar(hdr, 'NAXIS2')
    nv = sxpar(hdr, 'NAXIS3')
    xr = [0,0,nx-1,nx-1]
    yr = [0,ny-1,0,ny-1]
    xyad, hdr, xr, yr, l, b
endelse
;----------------calculate the integration and make the plot-------------------
cv = sxpar(hdr,'CRVAL3')
cp = sxpar(hdr,'CRPIX3')
cd = sxpar(hdr,'CDELT3')/1000d
height = (pos[3] - pos[1])/ygrid
width = (pos[2] - pos[0])/xgrid
textpos = [max(l)-(max(l)-min(l))*textp[0], min(b)+(max(b)-min(b))*textp[1]]
cgps_open, file_basename(fitsname, '.fits')+'_cplot.eps', /portrait, xsize = xgrid, ysize = ygrid*ny/nx, /encapsulated
;-----------some color and line setting--------------------------
!p.CHARSIZE=0.5
!P.Thick = 1
!P.CharThick = 1
!X.Thick = 1
!Y.Thick = 1
!Z.Thick = 1
for i=0, xgrid-1 do begin 
    for j=0, ygrid-1 do begin 
        start_v = (v[0]+(i+j*xgrid)*interval)
        end_v = (v[0]+(i+1+j*xgrid)*interval)
        print, start_v, end_v
        start_channel = (start_v-cv)/cd+cp-1
        end_channel = (end_v-cv)/cd+cp-1
        map = total(data[*,*,start_channel:end_channel], 3)*abs(cd)
        img = bytscl(map, min=trange[0], max=trange[1])
        plot_pos = [pos[0] + i*width, pos[3] - (j+1)*height, pos[0] + (i+1)*width, pos[3] - j*height]
        cgimage, img, position = plot_pos, /keep_aspect_ratio, /noerase
        if i eq 0 and j eq ygrid-1 then begin
            cgcontour, img, level =1, position = plot_pos, /nodata, /noerase, xrange = [max(l), min(l)], yrange = [min(b), max(b)],$
                xtitle='Galactic Longitude (!Uo!N)', ytitle = 'Galactic Latitude (!Uo!N)', charsize=0.5, ytickformat='(F0.1)', $
                xticklen=0.03, yticklen=0.03, xminor=10, xticks = 2;, xtickformat = '(I)';, xtickv=[110, 111, 112, 113]
            cgtext, textpos[0], textpos[1], strtrim(string(start_v + interval/2, format = '(F8.2)'),2), /data, color=textcolor     
        endif else begin 
            cgcontour, img, level =1, position = plot_pos, /nodata, /noerase, xrange = [max(l), min(l)], yrange = [min(b), max(b)],$
                xtickformat='(A1)', ytickformat='(A1)', xticklen=0.03, yticklen=0.03, xminor=10, xticks = 2;, /overplot
            cgtext, textpos[0], textpos[1], strtrim(string(start_v + interval/2, format = '(F8.2)'),2), /data, color=textcolor
        endelse  
    endfor
endfor  
cgcolorbar, position = [plot_pos[2], plot_pos[1], plot_pos[2]+0.02, plot_pos[3]], range = [trange[0],trange[1]], /vertical, /right,$
    title = '(K km/s)'
cgps_close
end 
