# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 17:09:02 2018

@author: Arslaan Khalid
"""

#==============================================================================
# GFS grib data downloader and reader 
#==============================================================================
 
import matplotlib.pyplot as plt # Import the Matplotlib package
from osgeo import gdal # Import the GDAL library
import sys,getopt
import pandas as pd
import numpy.core.multiarray 
import wget




# file dowloading

for i in range(0,84): 
    file_number = '%03d'%i
    print('Beginning file download with wget module: forecated hour {}'.format(file_number))
    

    url = 'http://ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs.2018050206/gfs.t06z.pgrb2.0p25.f{}'.format(file_number) # change this file name to the current advisory
    output_file = '/home/fhrl/Documents/GFS_Forced/forecast/{}.grb2'.format(i)
    
    wget.download(url,output_file )


# creating a new file 
#open('/home/fhrl/Documents/GFS_Forced/fort.22', 'w').close()


for i in range(0,85):
    # Read the GRIB file
    print(i)
    file_number = '%03d'%i
    print('Beginning file download with wget module: forecated hour {}'.format(file_number))
    

    #url = 'http://ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/gfs.2018050206/gfs.t06z.pgrb2.0p25.f{}'.format(file_number) 
   
    input_file = '/home/fhrl/Documents/GFS_Forced/forecast/{}.grb2'.format(i)

    print(input_file)
    grib = gdal.Open(input_file)         # one alternative is to download the file or else read through online, takes about the same time
    
    
    ds = grib 
    # Grid Information
    
    #lat-lon grid:(1440 x 721) units 1e-06 input WE:NS output WE:SN res 48
    #lat 90.000000 to -90.000000 by 0.250000
    #lon 0.000000 to 359.750000 by 0.250000 #points=1038240
    
    #////////////////
    
    if i == 0:
        u_value=266
        v_value=267
        p_value=352
    else:
        u_value=288
        v_value=289
        p_value=415
        
    #u_data of wind
    u_wind = grib.GetRasterBand(u_value)      # 288 if after 0 hour 266
    u_data = u_wind.ReadAsArray() 
    
    #v_data of wind
    v_wind = grib.GetRasterBand(v_value)     #289 if after 0 hour else 267
    v_data = v_wind.ReadAsArray() 
    
    
    #P_data of wind
    P = grib.GetRasterBand(p_value)         #415 if after 0 hour else 352
    P_data = P.ReadAsArray() 
    
    
    import numpy as np 
    
    wind_mag = np.sqrt(u_data**2 + v_data**2)
    
    
    ##### subsection of a area
    # Grid Information
    
    #lat-lon grid:(220 x 255)
    #lat 60 by 0.25
    #lon -110 by 0.25
    # for ADCIRC control file use the following grid information
    #255 220 60.000000 -110.000000 0.250000 0.250000 3600  ! NWLAT, NWLON, WLATMAX, WLONMIN, WLATINC, WLONINC, WTIMINC
    
    
    #========================================== u wind 10m ====================================
    domain1_u = u_data[120:375] 
    
    domain11_u = numpy.delete(domain1_u, numpy.s_[0:1000], axis=1)# right of US cut
    domain11_u = numpy.delete(domain11_u, numpy.s_[220:], axis=1)   # left of US cut
    #plt.imshow(domain11_u, cmap='coolwarm')  # to see the the Area of DATA cut
    #========================================== v wind 10m ====================================
    domain1_v = v_data[120:375]
   
    
    domain11_v = numpy.delete(domain1_v, numpy.s_[0:1000], axis=1)# right of US cut
    domain11_v = numpy.delete(domain11_v, numpy.s_[220:], axis=1)   # left of US cut
    #plt.imshow(domain11_v, cmap='coolwarm')  # to see the the Area of DATA cut
    #========================================== MSLP ====================================
    domain1_p = P_data[120:375] 
    
    domain11_p = numpy.delete(domain1_p, numpy.s_[0:1000], axis=1)# right of US cut
    domain11_p = numpy.delete(domain11_p, numpy.s_[220:], axis=1)   # left of US cut
    #plt.imshow(domain11_p, cmap='coolwarm')  # to see the the Area of DATA cut
    
    
    #==============================================================
    container_u = []
    
    for i in range(0,255):
        for j in range(0,220):
            print(round(domain11_u[i][j],1))
            key1 = round(domain11_u[i][j],1)
            container_u.append(key1)
    #============================================================    
    container_v = []
    
    for i in range(0,255):
        for j in range(0,220):
            print(round(domain11_v[i][j],1))
            key2 = round(domain11_v[i][j],1)
            container_v.append(key2)
    
    #=============================================================
    container_p = []
    
    for i in range(0,255):
        for j in range(0,220):
            print(round(domain11_p[i][j],0))
            key3 = round(domain11_p[i][j],0)
            container_p.append(key3)
     
     
    #=============================================================
    
     
    #=============================================================
    GFS_data = pd.DataFrame()
    
    GFS_data['u']=container_u
    GFS_data['v']=container_v
    GFS_data['p']=container_p
    
    # writing the data to file
    
    with open('/home/fhrl/Documents/GFS_Forced/fort.22','a') as f2:
       GFS_data.to_csv(f2, sep='\t',index=False,header = False)   





# for Data Visuals of the Entire GFS forecast


'''
    # Read the number of bands inside the GRIB file
    number_of_bands = grib.RasterCount
    print(number_of_bands)
    
    
    # Get the band name and description
    for i in range(1,418):
        print(i)
        band = grib.GetRasterBand(i)
        metadata = band.GetMetadata()
        band_name = metadata['GRIB_COMMENT']
        band_description = band.GetDescription()
        print(band_name)
        print(band_description)
        print("-----------------------------")
    
    
    
    # Show the image
    plt.imshow(u_data, cmap='coolwarm')
    
    plt.imshow(v_data, cmap='jet')
    
    plt.imshow(P_data)
    
    plt.imshow(wind_mag, cmap='coolwarm')
'''
