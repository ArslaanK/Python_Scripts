# -*- coding: utf-8 -*-
"""
Created on Sat Jun 29 15:52:27 2019

@author: Arslaan Khalid
"""

import pandas as pd
import numpy as np
from datetime import datetime, timezone, timedelta
import html2text
import requests
import re


# =============================================================================
#  User inputs
# =============================================================================
# find the details on the storms in the best track format at this link http://tropicalatlantic.com/models/models.cgi?basin=al&archive=2019

##----Storm (for example Irma)
name='Ophelia'
year= 2005
StormNo= 16
StormNo = '%02d'%StormNo
Region="al" # for atlantic

#specify start date of the storm

strt_date=2005090606 # make sure the date is in 06,12,18,00 UTC time
length = 11 #in days


# =============================================================================
#   code block-
# =============================================================================


tmpdate = datetime.strptime(str(strt_date),'%Y%m%d%H')+ timedelta(days=length)

end_date = int(datetime.strftime(tmpdate,'%Y%m%d%H'))

#----besttrack path




url = r'http://tropicalatlantic.com/models/archive/{}/{}/{}/atcf/nhc/Best_Track.txt'.format(Region,year,StormNo)
res = requests.get(url)
text = res.text

data = text.split('\n')
# if you have to download the file

#f = open(r'C:\Users\Arslaan Khalid\Documents\FirstStreet\Irma\cera\atcf', 'r')  # We need to re-open the file
#data = f.read().split('\n')


#----- crop the data from the user specified start and end date

for k in range(0,len(data)-1):
    print(k)
    #k=1
    tmpstrt1 = data[k].split(',')
    
    if tmpstrt1[2][1:] == str(strt_date):
        print('found the user defined start dataa at line:',k)
        break
 
#----- subset the data
data = data[k:]



for k in range(0,len(data)-1):
    print(k)
    #k=1
    linebyline = data[k].split(',')


    diff = datetime.strptime(linebyline[2][1:],'%Y%m%d%H')- datetime.strptime(str(strt_date),'%Y%m%d%H')
    duration_in_s = diff.total_seconds()
    hours = divmod(duration_in_s, 3600)[0]



    if hours < 10:
        linebyline[5] = f'   {int(hours)}'
    elif hours < 100 and hours > 10:
        linebyline[5] = f'  {int(hours)}'        
    elif hours > 100:
        linebyline[5] = f' {int(hours)}'

    makeitastring = ','.join(map(str, linebyline))
    data[k] = makeitastring


makeitastring2 = '\n'.join(map(str, data))

with open(r'C:\Users\Arslaan Khalid\Documents\FirstStreet\Ophelia2005\{}.22'.format(name),'w') as f:
        f.write(makeitastring2)

f.close()


## --- Method 2 ------ couldn't get this to work though
#
#
#original_file = pd.read_csv(url,sep='\s+',header=None, error_bad_lines=False)
#
#start_date = original_file[2][0]
#
#print('start date of the storm data is:',start_date)
#
#
##----- crop the data from the user specified start and end date
#
#
#strt1ind =np.where(original_file[2]==str(strt_date)+',')
#end1ind =np.where(original_file[2]==str(end_date)+',')
#
#
##----- subset the data
#
#
#df = original_file.drop(original_file.index[:strt1ind[0][0]])
#
#
#df = df.reset_index()
#
#df = df.drop('index', axis=1)
##------ add the hours from the start of the subset dataframe
#
#for k in range(len(df)):
#    print(k)
#    start_date_tmp = df[2][k]
#    print(start_date_tmp)
#    
#    #convert to datetime
#
#    diff = datetime.strptime(start_date_tmp[:-1],'%Y%m%d%H')- datetime.strptime(str(strt_date),'%Y%m%d%H')
#    duration_in_s = diff.total_seconds()
#    hours = divmod(duration_in_s, 3600)[0]
#    
#    df.iloc[k,df.columns.get_loc(5)] = str(int(hours))+','
#    
#    
##------------ writing the file------
#   
#    
#use_cols=[ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
#       18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,
#       36]
#
#df[1] = df[1].map('${:,.2f}'.format)
#
#with open(r'C:\Users\Arslaan Khalid\Documents\FirstStreet\Irma\cera\fort.22','w') as f:
#        df.to_string(f, index=False,header = False,col_space=1)
#

    
    
    
    
    
 
