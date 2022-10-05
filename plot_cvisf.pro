pro plot_cvisf
cgps_open, 'CVISF.eps', xsize = 7, ysize = 5, /encapsulated
!P.charsize = 1.3
!P.thick = 2
pos = [0.2, 0.2, 0.9, 0.9]
sym = [0, 4, 5, 6, 7, 8, 9]
cgplot, 0, 0, /nodata, xrange = [1, 1000], yrange = [-2, 4], xtitle = 'lag (pixel)', ytitle = 'lg CVISF!Dp!N (lag)',$
       /xlog, position = pos
col = ['forest green', 'pink', 'dodger blue', 'orchid', 'orange', 'chocolate']; 'lg CVISF!Dp!N (lag)'
items = ['p = 1', 'p = 2', 'p = 3', 'p = 4', 'p = 5', 'p = 6']
for i = 1, 6 do begin 
    cvi = 'cep_cvi' + strtrim(string(i, format = '(I4)'), 2) + '.fits'
    fits_read, cvi, data
    s = size(data, /dim)
    lag = 2. + 2.*dindgen(fix(sqrt(s[0]^2 + s[1]^2)/2d))
    sf = dblarr(s[2])
    err = dblarr(s[2])
    for j = 0, s[2]-1 do begin 
        slice = abs(data[*, *, j])
        good = where(slice ne 0)
        numtotal = n_elements(where(data[*, *, 0] ne 0))
        num = n_elements(good)
        ratio = num/float(numtotal)
        if good[0] ne -1 and ratio gt 0.5 then begin 
            sf[j] = mean(slice[good])
            err[j] = stddev(slice[good])
        endif
    endfor
    er = 1d/(sf*alog(10))*err
    cgplot, lag, alog10(sf), psym = -16, color = col[i-1], /overplot, symsize = 0.4;, err_yhigh = er, err_ylow = er
    result = linfit(alog10(lag[where(lag lt 35 and lag ge 8)]), alog10(sf[where(lag lt 35 and lag ge 8)]), yfit = yfit)
    cgplot, lag[where(lag lt 35 and lag ge 8)], yfit, /overplot, color = 'black', linestyle = 0
    print, result
    items[i-1] = items[i-1] + textoidl(' \alpha = ') + strtrim(string(result[1], format = '(f7.2)'), 2)
endfor 
    al_legend, items, /left, charsize = 1, box = 0, $
               textcolors = ['forest green', 'pink', 'dodger blue', 'orchid', 'orange', 'chocolate']
cgps_close
end  