# IDLscripts
## Description

### cvi.pro 
- To calculate the increment maps of a given centroid velocity map, which are further used to calculate the velocity structure function. The velocity increment is defined as the n-th order (n=1-6 here) of the velocity difference with different lags. The output of the script is six fits files recording 1-6th order of velocity increments, respectively.

### plot_cvi.pro 
- Make a plot of velocity structure functions using the calculated centroid velocity increment maps; also, fit each order of the structure-function with a power-law function within a selected range and overlay the fitted results on the figure. The 'CVISF.eps' in the example folder is given as an example.

### sps.pro 
- Calculate the azimuthal averaged spatial power-spectrum using the FFT of a 2D image, plot the spatial power-spectrum (SPS), and fit the SPS with a power-law function. The fitted results are labeled on the output image. The figure './example/sps_fft.eps' shows the intensity image and its 2D power image converted through FFT. The figure './example/annuli_stat.eps' shows the annuali statistics of the power extracted along the concentric circles in the 2D power image. The figure './example/sps_sps.eps' shows the annuli averaged spatial power-spectrum of the input image, in which the lag is converted using physical distance.   

### video.pro
- Make a movie of an input 3D fits file. The value of each velocity channel is labeled on the movie. The media file 'ngc7538.webm' in the example folder is provided as an example.

### channel_plot.pro
- To make a velocity channel map for a given 3D data cube. Users can specify whether to plot the whole region or a subregion, spatially, of the given cube, adjust their favorite color maps, and set the corlor scale by adjusting the first few lines in the script. The figure 'ngc7538U_cplot.eps' in the example folder is provided as an example. 
 

### plot_slice.pro
- To make a p-v slice image, along any given positional tracks, from given 3D data cubes of the 12CO and 13CO data. Users can define the width of the extracted slice by providing a slice 'width'; when the width is not 0, the intensity of the extracted slice is averaged perpendicular to the track. The procedure produces two extracted slice files in fits format and makes a 'eps' figure file showing the extracted slices with 12CO and 13CO intensities as color-filled and outlined contours, respectively. The figure 'knots1_slice.eps' in the example folder is provided as an example.

### tools.pro
> General image tools for fits files.
- plt, file : display a fits image automatically.
- plt_rgb, redfile, greenfile, bluefile : make a color-coded image for three given fits files.
- plt_pv, file : display a p-v fits file (either l-v or b-v)
- im3d, file : display a 3D datacube using the xyvolumn routine.
- tex : calculate excitation temperature from a given 12CO fits cube within a given velocity interval, using only the spectra having at least five contiguous channels around their peak intensities. 
- tau13 : calculate the 13CO optical depth using spectra with S/N ratios above 4.
- cal_n : calculate the H2 column density map of a cloud under the local thermal equilibrium (LTE) assumption.  

