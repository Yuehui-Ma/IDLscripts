function circle, xcenter, ycenter, radius, npoints
points = 2 * !pi / double(npoints) * dindgen(npoints)
x = xcenter + radius * cos(points )
y = ycenter + radius * sin(points )
return, transpose([[x],[y]])
end

pro cvi
input = 'cep_m1.fits'
step = 1.
fits_read, input, data, hdr
vmean = mean(data[where(finite(data))])
data = (data - vmean)/stdev(data[where(finite(data))])
s = size(data, /dim)
lag = 2. + 2.*dindgen(fix(sqrt(s[0]^2 + s[1]^2)/2d))
n = n_elements(lag)
dcv = dblarr(s[0], s[1], n)
for l = 1, 6 do begin 
    for i = 0, s[0]-1 do begin 
        for j = 0, s[1]-1 do begin 
            xcenter = i
            ycenter = j
            v0 = data[xcenter, ycenter]
            if finite(v0) then begin  
                for k = 0, n-1 do begin 
                    radius = lag[k]
                    npoints = fix(2*!pi*radius/step) 
                    loc = circle(xcenter, ycenter, radius, npoints)
                    vr = interpolate(data, loc[0, *], loc[1, *], missing = 0)
                    good = where(finite(vr) and vr ne 0)
                    if good[0] ne -1 then begin 
                        if l eq 1 then dcv[i, j, k] = mean((vr[good] - v0)^l) else dcv[i, j, k] = mean((abs(vr[good] - v0))^l)
                    endif 
                endfor 
            endif 
        endfor
    endfor
    fits_write, 'cep_cvi'+strtrim(string(l, format = '(I4)'), 2)+'.fits', dcv
endfor 
end