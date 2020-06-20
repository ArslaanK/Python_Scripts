# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 16:13:19 2020

@author: Arslaan Khalid
# for queries please contact @ akhalid6@gmu.edu
# following script prepare the field based plots from adcirc netcdf files
"""

# =============================================================================
# Loading Libraries
# =============================================================================

import matplotlib.pyplot as plt # Import the Matplotlib package
from mpl_toolkits.basemap import Basemap
import scipy
import scipy.interpolate
import netCDF4
import matplotlib.tri as tri
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
import matplotlib.lines as mlines

# =============================================================================
#  USER INPUTS
# =============================================================================

Input_dir = 'path to files'
Out_dir = Input_dir#'path to processed plots'
plot_type= 'WaterLevels' # choose from these [Winds&Pressure,WaterLevels] [Current and Inudnation will be added later]

stations_plot=True # True or False [this will plot the location of sensors]
Figur_title='Hurricane Maria & Jose [18Sep-29Sep 2017]'

paneltype = [3,3] # [columns, rows] [add that many timesteps in the list below, i.e. 3x3 will need 9 timesteps]
timestepstoview=[12,24,36,48,60,72,84,96,104] # add the index of time step to view on a 3x3 panel i.e. [12,24,36,48,60,72,84,96,104]


# stations info
wavegage_identifier = ['G1','G2','G3']
wavegage_lats=[37.149697,37.149598,37.149375]
wavegage_lons=[-75.940369,-75.940488,-75.94113]

wlgage_identifier=['S1T2','S1T2','S1T3','S1T4']
watergage_lat = [37.149697,37.149598,37.149375,37.148875]
watergage_lon = [-75.940369,-75.940488,-75.94113,-75.943921]


colorcodes=['green','red','blue','m']

#plotting selections
mesh_triangles_show=False # use for water level plot only
llcrnrlat0,urcrnrlat0,llcrnrlon0,urcrnrlon0= 37.1467986931,37.1523,284.05509773619997,284.062# bounding coordinates of the map box 
                                            #i.e.5,50,250,310 for US North Atlantic (EC95 domain)
                                            #i.e.37.1467986931,37.1523,284.05509773619997,284.062 for Magothy Bay
                                            
coastlineresolution = 'f'# choose from [c,l,h,f]
ESRIimageshow= True # True will show the satelittle image, False will not
ESRIimageresolution = 150 # choose from [1500,750,500,250]; smaller for oceanscale maps

contourwinds=np.arange(0, 80, 2.5) # for example (start, end, interval)
contourpressure=np.arange(970, 1040, 5)
controurwaterlevels= np.arange(-2, 2, 0.25)

grid_space_vectors = 10 # smaller values give more vectros and vise versa
vector_type= 'uniform' # choose from [magnitude, uniform]
# =============================================================================
#  Loading input files for Field plots 
# [fort.63.nc,fort.64.nc,fort.74.nc,fort.73.nc]
# =============================================================================

if plot_type =='Winds&Pressure':
    stdin1 = Input_dir + '//fort.74.nc' # wind file
    stdin2 = Input_dir + '//fort.73.nc' # pressure file
elif plot_type =='WaterLevels':
    stdin1 = Input_dir + '//fort.63.nc'
elif plot_type =='Current':
    stdin1 = Input_dir + '//fort.64.nc'
elif plot_type =='Inudnation':
    stdin1 = Input_dir + '//inundation.63.nc' 


# =============================================================================
#   MAIN CODE
# =============================================================================

m = Basemap(projection='cyl',llcrnrlat=llcrnrlat0,urcrnrlat=urcrnrlat0,llcrnrlon=llcrnrlon0,\
            urcrnrlon=urcrnrlon0,resolution=coastlineresolution, epsg = 4269)

# read the basic mesh parameters
file1 = netCDF4.Dataset(stdin1)

lat  = file1.variables['y'][:]
lon  = file1.variables['x'][:]    
depth  = file1.variables['depth'][:]  
elems = file1.variables['element'][:,:]-1  # Move to 0-indexing by subtracting 1, elements indexing starts with '1' in netcdf file         
lon[lon < 180] = lon[lon < 180] + 360 # if you want to change the longs to degrees from Easts

time_var_re = file1.variables['time'] # read the datetime
# sometimes the format isn't good so did this
mn_num = datetime.strptime(f"{time_var_re.units.split('-')[1]}",'%b').month
yr_num = f"{time_var_re.units.split(' ')[2][-4:]}"
dt_num = f"{time_var_re.units.split(' ')[2][:2]}"

Time_data_re = netCDF4.num2date(time_var_re[:],f"seconds since {yr_num}-{mn_num}-{dt_num} 00:00:00 +00:00")


# creating triangulation
triang = tri.Triangulation(lon,lat, triangles=elems)

# read the necessary file for plot
if plot_type =='Winds&Pressure': 
    file1 = netCDF4.Dataset(stdin1)
    file2 = netCDF4.Dataset(stdin2)    
else:
    file1 = netCDF4.Dataset(stdin1)

# initialize a figure

fig = plt.figure(figsize = (9,12))
fig.subplots_adjust(wspace=0.1)

# do the loop for each subplot
for y in range(0,(paneltype[0]*paneltype[1])):   
    i=timestepstoview[y]
    file_number = '%02d'%i

    plt.subplot(paneltype[0],paneltype[1],1+y)   
    # =============================================================================
    #  plotting scheme
    # =============================================================================

    if plot_type =='Winds&Pressure': 
        windX,windY = file1.variables['windx'][i,:], file1.variables['windy'][i,:]
        Pressure= file2.variables['pressure'][i,:].data*98.0638 ## converting to mbars
    
        Winds = np.sqrt(np.square(windX)+np.square(windY)).data

        plt.tricontourf(triang, Winds, levels=contourwinds,alpha=0.9,cmap='jet')
        CS=plt.tricontour(triang, Pressure, levels=contourpressure,alpha=0.9,colors='k',linestyles='--',linewidth=0.75)#,cmap='k'
        plt.clabel(CS, CS.levels, fmt='%.0f', fontsize=10, colors='w') #,fontweight='bold'
 
        m.fillcontinents(color='#B0B0B0',alpha=0.65)
        #=====Making the directionional data
        dff = pd.DataFrame()
    
        dff['lat']=lat;dff['lon']=lon;dff['depth']=depth;dff['u']=windX;dff['v']=windY
        
        lon11 = np.array(dff['lon']);lat11 = np.array(dff['lat'])
        u11 = np.array(dff['u']);v11 = np.array(dff['v'])    
        #regriding data
        grid_space=grid_space_vectors

        lon1=llcrnrlon0;lon2=urcrnrlon0;
        lat1=llcrnrlat0;lat2=urcrnrlat0
        xg = np.linspace(lon1+0.25,lon2-0.25,grid_space)
        yg = np.linspace(lat1+0.25,lat2-0.25,grid_space)
        
        xgrid,ygrid = np.meshgrid(xg,yg)
        xx2=lon;yy2=lat             
        u=windX;v=windY
        # magnitude based vector data
        ugrid = scipy.interpolate.griddata((xx2,yy2),u,(xgrid,ygrid),method='nearest')
        vgrid = scipy.interpolate.griddata((xx2,yy2),v,(xgrid,ygrid),method='nearest')
        #--normalized vector data
        u_norm = ugrid / np.sqrt(ugrid ** 2.0 + vgrid ** 2.0)
        v_norm = vgrid / np.sqrt(ugrid ** 2.0 + vgrid ** 2.0)
        
        if vector_type =='uniform':
            Q=plt.quiver(xgrid,ygrid,u_norm,v_norm, pivot='middle', scale = 20, color='w')#,headlength =20
        elif vector_type=="magnitude": 
            Q=plt.quiver(xgrid,ygrid,ugrid,vgrid, pivot='tail', color='w',scale=50)#,headlength =20,scale=25 units ='xy',
         
    elif plot_type=="WaterLevels":

        WL= file1.variables['zeta'][i,:].data
        plt.tricontourf(triang, WL, levels=controurwaterlevels,alpha=0.9,cmap='jet')
        if mesh_triangles_show == True:
            plt.triplot(triang,c='k',linewidth=0.25,alpha=0.75)    # plt the mesh triangles     

    else:
        print('capability not added yet')
          
   
    # plotting some basic info on the plot        
    m.drawcoastlines(linewidth=0.5)            
     
                     
    m.drawcountries()
    #various other basemap options here
    #url:https://basemaptutorial.readthedocs.io/en/latest/backgrounds.html

    if stations_plot==True:
        if ESRIimageshow==True:
            m.arcgisimage(service='World_Imagery', xpixels = 1500, verbose= False) #1500
            
        if len(wavegage_identifier)>0:     
            for t in range(0,len(wavegage_identifier)):
                plt.scatter(wavegage_lons[t]+360,wavegage_lats[t],marker = 'o',s=25, color=f'{colorcodes[t]}') #,facecolor='none'
                plt.text(wavegage_lons[t]+360,wavegage_lats[t]+0.0005,f'{wavegage_identifier[t]}',color=f'{colorcodes[t]}')      
        if len(wlgage_identifier)>0:        
            for t in range(0,len(wlgage_identifier)):
                plt.scatter(watergage_lon[t]+360,watergage_lat[t],marker = 's',s=25, color=f'{colorcodes[t]}') #,facecolor='none'
                plt.text(watergage_lon[t]+360,watergage_lat[t]+0.0005,f'{wlgage_identifier[t]}',color=f'{colorcodes[t]}')      
                                            

    plt.title(f'{str(Time_data_re[i])[:16]}')
        
    #===========end======       

# =============================================================================
# Adding the 1 legend
# =============================================================================


colorbar_ax = fig.add_axes([0.90, 0.1, 0.03, 0.78]) # these numbers may need to change

if plot_type=="Winds&Pressure": 
    cd2 = plt.tricontourf(triang, Winds, levels=contourwinds,alpha=0.9,cmap='jet')
    fig.colorbar(cd2, cax=colorbar_ax,label='Wind magnitude (meters per second m/s)',format = "%.1f") 
    ln00 = mlines.Line2D([], [],c='k',label='Pressure Contours',linestyle='--',mfc='none',ms=8,lw=1.5)
    fig.legend(handles=[ln00],ncol=1,fontsize=10,fancybox=False, shadow=False,frameon=False, \
               bbox_to_anchor=(0.75,0.04),loc='lower left')#bbox_to_anchor numbers may need to change if not show

    if vector_type =='uniform':    
        plt.quiverkey(Q, 0.87, 0.9,5, 'Wind Direction', coordinates='figure',color='k')
    else:
        plt.quiverkey(Q, 0.87, 0.9,40, 'Wind Direction, Magnitude = 40m/s', coordinates='figure',color='k')
  

elif plot_type=="WaterLevels":
    cd2 = plt.tricontourf(triang, WL, levels=controurwaterlevels,alpha=0.9,cmap='jet')
    fig.colorbar(cd2, cax=colorbar_ax,label='Water levels[meters above NAVD88]',format = "%.1f")


# fit subplots and save fig
#fig.tight_layout()

fig.subplots_adjust(
top=0.88,
bottom=0.11,
left=0.11,
right=0.9,
hspace=0.28,
wspace=0.0)


fig.suptitle(Figur_title)

fig.set_size_inches(w=12,h=12)

fig.savefig(r'{}\{}.png'.format(Out_dir,plot_type),  bbox_inches = 'tight',dpi = 200, pad_inches = 0.2)
print(f'Successfully created the Figure Type: {plot_type}')

