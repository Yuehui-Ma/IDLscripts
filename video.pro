PRO VIDEO
;-------------------settings------------------------------------------------
;---------------------------------------------------------------------------
;---------------------------------------------------------------------------
imfile = 'ngc7538U.fits'
center = [109.04, 0.31] ; 也可以设置某个中心，截取一定范围的 fits cube 做视频，不设置 psize 时，该参数不起作用。
psize = 0;[0.5d, 0.5d] 需要截取空间范围时设置，意义为以center 为中心，x和y方向的边长，单位是度。不需截取时设置为0.
vrange = [-60, -30];需要制作为视频的速度范围。
fps = 8
ctab = 70 ;颜色表
min = 1.5  ;颜色显示范围的最小值，默认 2 sigma。
max = 20. ;颜色显示范围的最大值。
;---------------------------------------------------------------------------
;---------------------------------------------------------------------------
;---------------------------------------------------------------------------
fits_read, imfile, data, hdr
nc = sxpar(hdr,'NAXIS3')
cv = sxpar(hdr,'CRVAL3')
cp = sxpar(hdr,'CRPIX3')
cd = sxpar(hdr,'CDELT3')
if abs(cd) gt 100 then vrange = vrange * 1000d   ;km/s to m/s
if n_elements(vrange) eq 2 then begin 
c = (vrange-cv)/cd+cp-1
c = c[sort(c)]
endif else begin 
c = [0, nc-1]
endelse 

if (n_elements(psize) eq 2) then begin 
l = [center[0] + psize[0]/2d, center[0] - psize[0]/2d]
b = [center[1] - psize[1]/2d, center[1] + psize[1]/2d]
adxy, hdr, l, b, cx, cy
data = data[cx[0]:cx[1],cy[0]:cy[1],*]
endif
s = size(data, /dim)
width = s[0]
height = s[1]  
DEVICE, GET_DECOMPOSED=old_decomposed
device, decompose = 1
cgloadct, ctab, /reverse
stretch, 0, 255, 0.5
TVLCT, R, G, B, /GET
; Each of the following lines produces a file in
; a different format.
oVid = IDLffVideoWrite(file_basename(imfile, '.fits') + '.webm')
;oVid = IDLffVideoWrite('video_example_file_format.swf')
;oVid = IDLffVideoWrite('video_example_file_format.bin', FORMAT='mp4')
; Prints out a list of supported file formats
;PRINT, "Supported file formats: ", oVid.GetFormats()
vidStream = oVid.AddVideoStream(s[0]*2, s[1]*2, fps)

v = ((findgen(nc)-cp+1)*cd+cv)/1000d
FOR i = c[0], c[1]-1 do begin
  slice = data[*, *, i]
  frame = bytscl(smooth(slice, 2), min = min, max = max)
  img = bytarr(3, s[0], s[1])
  img[0, *, *] = R[frame] 
  img[1, *, *] = G[frame]
  img[2, *, *] = B[frame]
  WINDOW, 1, XSIZE=s[0]*2, YSIZE=s[1]*2
  cgimage, img
  cgtext, 'v ='+ strtrim(string(v[i], format = '(f8.2)'), 2)+' km s!U-1!N', 0.1, 0.9, /normal, font = 1, charthick = 4, charsize = 5
  ; Create an image and display it
  ; Write a bitmap file to the temporary directory
  ; Note the use of the TRUE keywords to TVRD and WRITE_JPEG
  filename = 'movie.png'
  ;WRITE_JPEG, filename, TVRD(/TRUE), /TRUE, quality = 100
  WRITE_png, filename, TVRD(/TRUE), XRESOLUTION = 200, yresolution = 200
  ;READ_JPEG, filename, IMAGE2, /TRUE
  image2 = read_png(filename)
  File_Delete, 'movie.png'
  ; Add the high-resolution image to the video stream.
  ;void = oVid -> Put(stream, image)
  !NULL = oVid.Put(vidStream, image2)
ENDFOR
oVid = 0
DEVICE, DECOMPOSED=old_decomposed
end
