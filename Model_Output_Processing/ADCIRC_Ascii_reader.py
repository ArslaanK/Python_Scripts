

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 27 11:16:43 2018

Created by: Sergey @ NOAA
@author: Arslaan Khalid
"""

import os
import numpy as np
from datetime import datetime
from datetime import timedelta
import netCDF4





#==============================================================================
def readGrid ( gridFile ):
    """
    Reads ADCIRC grid file
    
    Args:
        gridFile (str): full path to fort.14 file
    Returns:
        grid (dict): field names according to ADCIRC internal variables:
    http://adcirc.org/home/documentation/users-manual-v50/
    input-file-descriptions/adcirc-grid-and-boundary-information-file-fort-14/
    """
    print ('[info]: Reading the grid from ' + gridFile)
    if not os.path.exists (gridFile):
        print ('[error]: File ' + gridFile + ' does not exist.')
        return
        
    f  = open(gridFile)
    
    myDesc     = f.readline().rstrip()
    myNE, myNP = map(int, f.readline().split())    
    print ('[info]: Grid description ' + myDesc + '.')
    print ('[info]: Grid size: NE= '   + str(myNE) + ', NP=' + str(myNP) + '.')

    myPoints   = np.zeros([myNP,3], dtype=float)
    myElements = np.zeros([myNE,3], dtype=int)
    
    print ('[info]: Reading grid points...')
    for k in range(myNP):
        line            = f.readline().split()
        myPoints[k,0] = float(line[1])
        myPoints[k,1] = float(line[2])
        myPoints[k,2] = float(line[3])

    print ('[info]: Reading grid elements...')
    for k in range(myNE):
        line              = f.readline().split()
        #myElements[k,0:2] = map(int, line[2:4])
        myElements[k,0] = int (line[2])
        myElements[k,1] = int (line[3])
        myElements[k,2] = int (line[4])

    
    myNOPE   = int(f.readline().split()[0])
    myNETA   = int(f.readline().split()[0])   
    myNVDLL  = np.zeros([myNOPE], dtype=int)
    myNBDV   = np.zeros([myNOPE, myNETA], dtype=int)
    
    print ('[info]: Reading elevation-specified boundaries...' )   
    for k in range(myNOPE):
        myNVDLL [k] = int(f.readline().split()[0])
        for j in range(myNVDLL[k]):
            myNBDV[k,j] = int(f.readline().strip())

    myNBOU = int(f.readline().split()[0])
    myNVEL = int(f.readline().split()[0])   
    myNVELL      = np.zeros([myNBOU], dtype=int)
    myIBTYPE     = np.zeros([myNBOU], dtype=int)
    myNBVV       = np.zeros([myNBOU, myNVEL], dtype=int)
    myBARLANHT   = np.zeros([myNBOU, myNVEL], dtype=float)
    myBARLANCFSP = np.zeros([myNBOU, myNVEL], dtype=float)
    myIBCONN     = np.zeros([myNBOU, myNVEL], dtype=int)
    myBARINHT    = np.zeros([myNBOU, myNVEL], dtype=float)
    myBARINCFSB  = np.zeros([myNBOU, myNVEL], dtype=float)
    myBARINCFSP  = np.zeros([myNBOU, myNVEL], dtype=float)
    myPIPEHT     = np.zeros([myNBOU, myNVEL], dtype=float)
    myPIPECOEF   = np.zeros([myNBOU, myNVEL], dtype=float)
    myPIPEDIAM   = np.zeros([myNBOU, myNVEL], dtype=float)
    
    print ('[info]: Reading normal flow-specified boundaries...')    
    for k in range(myNBOU):
        line = f.readline().split()
        myNVELL[k]  = int(line[0])
        myIBTYPE[k] = int(line[1])
        
        for j in range(myNVELL[k]):
            line = f.readline().rstrip().split()            
            if myIBTYPE[k] in   [0,1,2,10,11,12,20,21,22,30]:
                myNBVV      [k,j] = int(line[0])
            elif myIBTYPE[k] in [3,13,23]:
                myNBVV      [k,j] = int  (line[0])
                myBARLANHT  [k,j] = float(line[1])
                myBARLANCFSP[k,j] = float(line[2])
            elif myIBTYPE[k] in [4,24]:
                myNBVV      [k,j] = int  (line[0])
                myIBCONN    [k,j] = int  (line[1])
                myBARINHT   [k,j] = float(line[2])
                myBARINCFSB [k,j] = float(line[3])
                myBARINCFSP [k,j] = float(line[4])
            elif myIBTYPE[k] in [5,25]:
                myNBVV      [k,j] = int  (line[0])
                myIBCONN    [k,j] = int  (line[1])
                myBARINHT   [k,j] = float(line[2])
                myBARINCFSB [k,j] = float(line[3])
                myBARINCFSP [k,j] = float(line[4])
                myPIPEHT    [k,j] = float(line[5])
                myPIPECOEF  [k,j] = float(line[6])
                myPIPEDIAM  [k,j] = float(line[7])

    f.close()
        
    return {'GridDescription'               : myDesc, 
            'NE'                            : myNE, 
            'NP'                            : myNP, 
            'lon'                           : np.squeeze(myPoints[:,0]),
            'lat'                           : np.squeeze(myPoints[:,1]), 
            'depth'                         : np.squeeze(myPoints[:,2]), 
            'Elements'                      : np.squeeze(myElements),
            'NETA'                          : myNETA, 
            'NOPE'                          : myNOPE,
            'ElevationBoundaries'           : np.squeeze(myNBDV), 
            'NormalFlowBoundaries'          : np.squeeze(myNBVV),
            'ExternalBarrierHeights'        : np.squeeze(myBARLANHT),
            'ExternalBarrierCFSPs'          : np.squeeze(myBARLANCFSP),
            'BackFaceNodeNormalFlow'        : np.squeeze(myIBCONN),
            'InternalBarrierHeights'        : np.squeeze(myBARINHT),
            'InternallBarrierCFSPs'         : np.squeeze(myBARINCFSP),
            'InternallBarrierCFSBs'         : np.squeeze(myBARINCFSB),            
            'CrossBarrierPipeHeights'       : np.squeeze(myPIPEHT),
            'BulkPipeFrictionFactors'       : np.squeeze(myPIPECOEF),            
            'CrossBarrierPipeDiameter'      : np.squeeze(myPIPEDIAM)
            }


#==============================================================================
def readTimeSeries (ncFile, ncVar = 'zeta'):
    """
    Reads fort.61.nc-like file
    """
    print ('[info]: Reading [' + ncVar + '] from ' + ncFile)
    if not os.path.exists (ncFile):
        print ('[error]: File ' + ncFile + ' does not exist.')
        return
    
    nc = netCDF4.Dataset( ncFile )
    fld        = nc.variables[ncVar][:]
    missingVal = nc.variables[ncVar]._FillValue
    try:
        fld.unshare_mask()
    except:
        pass
    fld [fld == missingVal] = np.nan
                          
    lon  = nc.variables['x'][:]
    lat  = nc.variables['y'][:]    
    tim  = nc.variables['time'][:]    
    
    baseDate = datetime.strptime(nc.variables['time'].base_date[0:19],
                                 '%Y-%m-%d %H:%M:%S')
    realtime = np.array([baseDate + 
                         timedelta(seconds=int(tim[i])) 
                         for i in range(len(tim))])
    
    nam  = nc.variables['station_name'][:][:]
    stations = [''.join([str(x) for x in nam[y]]) for y in range(len(nam))]

    ncTitle = nc.getncattr('title')

    return  {'lat'       : lat, 
             'lon'       : lon, 
             'time'      : realtime, 
             'base_date' : baseDate, 
             'zeta'      : fld, 
             'stations'  : stations,
             'title'     : ncTitle}        
    
#==============================================================================
def readSurfaceField ( ncFile, ncVar = 'zeta_max' ):  
    """
    Reads specified variable from the ADCIRC 2D netCDF output
    and grid points along with validation time.
    Args:
        'ncFile' (str): full path to netCDF file
        'ncVar'  (str): name of netCDF field
    Returns:
        dict: 'lon', 'lat', 'time', 'base_date', 'value', 'path', 'variable'
    """
    
    print ('[info]: Reading [' + ncVar + '] from ' + ncFile)
    if not os.path.exists (ncFile):
        print ('[error]: File ' + ncFile + ' does not exist.')
        return
           
    nc   = netCDF4.Dataset (ncFile)
    lon  = nc.variables['x'][:]
    lat  = nc.variables['y'][:]
    tim  = nc.variables['time'][:]
    fld  = nc.variables[ncVar][:] 

    missingVal = nc.variables[ncVar]._FillValue
    try:
        fld.unshare_mask()
    except:
        pass
    fld [fld==missingVal] = np.nan

    baseDate = datetime.strptime(nc.variables['time'].base_date[0:19],
                                 '%Y-%m-%d %H:%M:%S')
    realtime = np.array([baseDate +
                         timedelta(seconds=int(tim[i]))
                         for i in range(len(tim))])

    return { 'lon'      : lon, 
             'lat'      : lat, 
             'time'     : realtime, 
             'base_date': baseDate,
             'value'    : fld, 
             'path'     : ncFile,              
             'variable' : ncVar}

#==============================================================================
def readSurfaceField_ascii ( asciiFile ):
    """
    Reads ADCIRC 2D output file (e.g. mmaxele)
    Args:
        'asciiFile' (str): full path to ADCIRC 2D file in ASCII format
    Returns:
        value (np.array [NP, NS]), where NP - number of grid points, 
                                     and NS - number of datasets
    """
    print ('[info]: Reading ASCII file ' + asciiFile + '.')
    f  = open(asciiFile)
    
    myDesc = f.readline().strip()
    print ('[info]: Field description [' + myDesc + '].')
    line          = f.readline().split()    
    myNDSETSE     = int(line[0])
    myNP          = int(line[1])
#    myNSPOOLGE    = int(line[3])
#    myIRTYPE      = int(line[4])
#    dtdpXnspoolge = float(line[2])   
    line          = f.readline().split()
#    myTIME        = float(line[0])
#    myIT          = float(line[1])
    value = np.zeros([myNP,myNDSETSE], dtype=float)
    for s in range(myNDSETSE):
        for n in range(myNP):
            value[n,s] = float(f.readline().split()[1])    
    value = np.squeeze(value)
    
    fill_value = -99999.0
    value[value==fill_value]=np.nan
    
    return value 

#==============================================================================
def computeMaxele (ncFile):
    nc = netCDF4.Dataset(ncFile)
    zeta  = nc.variables['zeta'][:]         
    fill_value = nc._FillValue    
    zeta = np.ma.masked_equal(zeta,fill_value)
    return {'lon'   : nc.variables['x'][:],
            'lat'   : nc.variables['y'][:],
            'value' : np.amax(zeta, axis=0)}

#==============================================================================
def readFort14 ( fort14file ):
    """
    Reads ADCIRC fort.14 file
    """
    return readGrid (fort14file)


# =============================================================================
# loaidn the file paths
# =============================================================================
gridFile = r'E:\Projects\Charleston_DesignBridge\100yr\fort.14'
asciiFile = r'E:\Projects\Charleston_DesignBridge\100yr\fort.63'


# =============================================================================
# reading the files into variables
# =============================================================================

f14 = readGrid(gridFile)
f63 = readSurfaceField_ascii(asciiFile)

