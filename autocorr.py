#importing all my junk
from NordicARC import uvmultifit as uvm
import numpy as np
import time
import matplotlib.pyplot as plt
#import os
import os.path
from os import path
import sys
import numpy as np
path_to_gfunc='.'
sys.path.insert(0, path_to_gfunc)
import gfunc as gf
#os.system('import gfunc.py')

#thisant=(visdata['antenna1']==ant) | (visdata['antenna2']==ant)

#setting a dummy refant
#refant=0

#this is for Venki data
target='sgr_apr07'

#for image names
case='noshift_noflag'
auth='Venki'
refmeth='refanctfunc'
date='04032022'

#datams0='firsttwointflag_QSOB19.ms'
#os.system('rm -rf QSOB19_chanavg.ms')
#split(vis=datams0,outputvis='QSOB19_chanavg.ms',width=240,datacolumn='data')

#Initial data file name and name for channel-averaged file
datams2=target+'_flagcor.ms'
#datams2='Sagittarius_A_star_C.ms_CALIBRATED_SELFCAL'
#datams2='uvmfit_3d77_SGRA_onechan.ms'
#datams2='QSOB19_chanavg.ms'

flagfile=datams2[:-3]+'_noflags.ms'
os.system('rm -rf '+flagfile)
split(vis=datams2,outputvis=flagfile,keepflags=False,datacolumn='data')

datams1=flagfile

#Opening data and pulling necessary info
ms.open(datams1,nomodify=True)
ms.selectinit(reset=True)
visdata = ms.getdata(['antenna1','antenna2','data','axis_info','flag'],ifraxis=True)
#printing correlations
print(visdata['axis_info']['corr_axis'])

#Squeeze data then close ms
#visdata['data'] = np.where(visdata['flag']==True,np.nan,visdata['data'])

visdata['data'] = np.squeeze(visdata['data'])
print("data shape",visdata['data'].shape)
ms.close()

tchunk=visdata['axis_info']['time_axis']['MJDseconds']
smallvis=visdata['data']

allants=np.concatenate((visdata['antenna1'],visdata['antenna2']))
antlist=np.unique(allants)

nant=len(antlist)
nbl=int(nant*(nant-1)/2)
ntimes=len(visdata['axis_info']['time_axis']['MJDseconds'])


amplist=[]
a1list=[]
a2list=[]

i=0

alltimes=visdata['axis_info']['time_axis']['MJDseconds']
t1=alltimes[0]

for i in range(ntimes):
#while alltimes[i]<10000+t1:
    #amplist=[]
    #a1list=[]
    #a2list=[]
    for ant1 in np.unique(visdata['antenna1']):
        for ant2 in np.unique(visdata['antenna2']):
            if ant1 < ant2:
                thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
                if thisbase.sum()>0:
                    #potential 0 after thisbase and thistime
                    a1list.append(ant1)
                    a2list.append(ant2)
                    pt=visdata['data'][0][thisbase][0][i]
                    amp=np.absolute(pt)
                    if amp<0.1:
                        visdata['data'][0][thisbase[0][i]=np.nan
                        #plt.scatter(ant1,ant2)
                    #if amp>0.1:
                    #    plt.scatter(ant1,ant2,c='r')
                    #print(amp)
                    #amplist.append(amp)
                    #plt.scatter(ant1,ant2,)

    #plt.scatter(a1list,a2list,s=amplist)
    plt.savefig("./antplots/a1a2zeroamp_"+str(i)+".png",overwrite=True)
    plt.clf()
#sc=plt.scatter(a1list,a2list,c=amplist,vmin=np.min(amplist),vmax=np.max(amplist),cmap='viridis')
#plt.colorbar(sc)

#for ant in range(nant):
    #print("ANT",ant)
    #all baselines w/ this ant
    #thisbase=(visdata['antenna1']==ant) & (visdata['antenna2']==ant)
    #pt=visdata['data'][0][thisbase][0]
    #amp=np.absolute(pt)
    #plt.scatter(ant,ant,s=amp*2)

plt.show()
