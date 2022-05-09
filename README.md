# Map folds
R script to produce a map fold effect. Uses an [example map](https://geoscience.nt.gov.au/gemis/ntgsjspui/handle/1/81667) and [example SRTM elevation data](https://dwtkns.com/srtm30m/).

![alt text][hermannsburg_plot_folded]

[hermannsburg_plot_folded]: https://github.com/cverdel/map_folds/blob/main/hermannsburg_folded_map.png?raw=true

R script to create a folded map effect. It will produce a 3D render in an interactive rgl window.
```
#It's important to install the latest version of rayshader from Github
#install.packages("devtools")
#devtools::install_github("tylermorganwall/rayshader")

#Load image
image_url<-"https://github.com/cverdel/rayshader_experiment/raw/main/Hermannsburg_map.tif"
temp<-tempfile()
download.file(image_url, temp, mode="wb")
rgb = raster::brick(temp)

#Load elevation data
DEM_url<-"https://github.com/cverdel/rayshader_experiment/raw/main/Hermannsburg_DEM.tif"
dem<-tempfile()
download.file(DEM_url, dem, mode="wb")
r1 = raster::raster(dem)
plot(r1)
```
**Original elevation raster**
![alt text][r1]

[r1]: https://github.com/cverdel/map_folds/blob/main/r1.png?raw=true
```
#Splits map into rgb bands. The original file has 4 bands, so there's an extra "t" band below.
names(rgb) = c("r","g","b","t")
rgb_r = rayshader::raster_to_matrix(rgb$r)
rgb_g = rayshader::raster_to_matrix(rgb$g)
rgb_b = rayshader::raster_to_matrix(rgb$b)

map_array = array(0,dim=c(nrow(rgb_r),ncol(rgb_r),3))

map_array[,,1] = rgb_r/255 #Red 
map_array[,,2] = rgb_g/255 #Blue 
map_array[,,3] = rgb_b/255 #Green 
map_array = aperm(map_array, c(2,1,3))

#Create folds
fold_raster <- r1 * 0 + 0 #Creates a new, blank raster with the same dimensions as the original DEM
rdf<-as.data.frame(fold_raster, xy = TRUE) #Creates a dataframe from raster
rdf_v<-as.data.frame(fold_raster, xy = TRUE) #Creates another dataframe from raster

#Set parameters for a 4x3 panel folded map. Can adjust steepness of folds with slope values.
vpanel=4
hpanel=3
slope=50
slope_v=40
f=100000 #Degrees to meters conversion factor

#Cutoffs
range_x<-max(rdf$x)-min(rdf$x)
range_y<-max(rdf_v$y)-min(rdf_v$y)
x0<-min(rdf$x)
x1<-x0+range_x/vpanel
x2<-x1+range_x/vpanel
x3<-x2+range_x/vpanel
x4<-max(rdf$x)

y0<-min(rdf_v$y)
y1<-y0+range_y/hpanel
y2<-y1+range_y/hpanel
y3<-max(rdf_v$y)

#Creates horizontal panels (i.e., vertical folds)
rdf$layer<-ifelse(rdf$x>=x0 & rdf$x<x1 , (f*(x1-rdf$x))*tan(slope*pi/180), 
                  ifelse(rdf$x>=x1 & rdf$x<x2 , (f*(rdf$x-x1))*tan(slope*pi/180),
                         ifelse(rdf$x>=x2 & rdf$x<x3 , (f*(x3-rdf$x))*tan(slope*pi/180), 
                                ifelse(rdf$x>=x3 & rdf$x<=x4, (f*(rdf$x-x3))*tan(slope*pi/180), NA))))

rdf$layer<-rdf$layer*1

r2<-rasterFromXYZ(rdf) #Creates raster from dataframe
r2<-subset(r2, subset=2, drop=TRUE)
plot(r2)
```
![alt text][r2]

[r2]: https://github.com/cverdel/map_folds/blob/main/r2.png?raw=true
```
#Creates vertical panels (i.e., horizontal folds)
rdf_v$layer<-ifelse(rdf_v$y>=y0 & rdf_v$y<y1 , (f*(y1-rdf_v$y))*tan(slope_v*pi/180), 
                    ifelse(rdf_v$y>=y1 & rdf_v$y<y2 , (f*(rdf_v$y-y1))*tan(slope_v*pi/180),
                           ifelse(rdf_v$y>=y2, (f*(y3-rdf_v$y))*tan(slope_v*pi/180), NA)))

rdf_v$layer<-rdf_v$layer*1

r3<-rasterFromXYZ(rdf_v) #Creates raster from dataframe
r3<-subset(r3, subset=2, drop=TRUE)
plot(r3)
```
![alt text][r3]

[r3]: https://github.com/cverdel/map_folds/blob/main/r3.png?raw=true
```
#Raster math
elevation_final<-r1+0.08*r2+0.08*r3 #Combines the 3 rasters (original elevation data, vertical folds, and horizontal folds).
plot(elevation_final)
```
Final elevation raster
![alt text][elevation_final]

[elevation_final]: https://github.com/cverdel/map_folds/blob/main/elevation_final.png?raw=true
```
#Raster to matrix conversion of elevation data
el_matrix = rayshader::raster_to_matrix(elevation_final)

#Plot maps
#Reduce the size of the elevation data to speed up processing
small_el_matrix = resize_matrix(el_matrix, scale = 0.6)

#Render
zscale1=40 #Sets vertical exaggeration for rayshader calculations
ambient_layer = ambient_shade(small_el_matrix, zscale = zscale1, multicore = TRUE, maxsearch = 200)
ray_layer = ray_shade(small_el_matrix, zscale = zscale1, sunaltitude=35, sunangle=315, multicore = TRUE)

rgl::rgl.close() #Closes the rgl window, in case it's still open

zscale2=50 #Sets vertical exaggeration for plotting

#Plot in 3D
 (map_array) %>%
  add_shadow(ray_layer,0.3) %>%
  add_shadow(ambient_layer,0.1) %>%
  plot_3d(small_el_matrix, solid=FALSE, zscale=zscale2, background='#fdfdfd')
```
