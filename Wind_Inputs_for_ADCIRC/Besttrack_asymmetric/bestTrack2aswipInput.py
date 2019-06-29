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

year= 2017
StormNo = 11
Region="al" # for atlantic

#specify start date of the storm

strt_date=2017090506 # make sure the date is in 06,12,18,00 UTC time
length = 10 #in days


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

for k in range(0,len(data)):
    print(k)
    #k=1
    tmpstrt1 = data[k].split(',')
    
    if tmpstrt1[2][1:] == str(strt_date):
        print('found the user defined start dataa at line:',k)
        break
 
#----- subset the data
data = data[k:]



for k in range(0,len(data)):
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

with open(r'C:\Users\Arslaan Khalid\Documents\FirstStreet\Irma\cera\fort.22','w') as f:
        f.write(makeitastring2)

f.close()