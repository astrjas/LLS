import numpy as np
import time
import matplotlib.pyplot as plt
import os.path
from os import path
import sys
path_to_gfunc='.'
sys.path.insert(0, path_to_gfunc)
import gfunc_c5 as gf
from scipy.linalg import lstsq
import numpy.ma as ma
from math import e
import pickle
import math


attn=0
tavg='60s'
tavgint=60.
MJDdate=58232.00000000*3600.*24. #apr242018
#MJDdate=57850.00000000 #apr072017
nspw=4


#try:
#    import uvmultimodel
#except ModuleNotFoundError:
import sys
sys.path.insert(0,"/home/jasminwash/.local/lib/python2.7/site-packages/NordicARC-3.0.5-py2.7-linux-x86_64.egg")
#import uvmultimodel


from NordicARC import uvmultifit as uvm

def msflag(dfile,rmsval):
    tbdata={}
    tb.open(dfile)
    tbdata['data']=tb.getcol('DATA')
    tbdata['data']=np.squeeze(tbdata['data'])
    tbdata['scan']=tb.getcol('SCAN_NUMBER')
    tbdata['ant1']=tb.getcol('ANTENNA1')
    tbdata['ant2']=tb.getcol('ANTENNA2')
    tbdata['flag']=tb.getcol('FLAG')
    tbdata['flag']=np.squeeze(tbdata['flag'])
    tb.close()

    allants=np.concatenate((tbdata['ant1'],tbdata['ant2']))
    antlist=np.unique(allants)
    nant=len(antlist)

    flagmanager(vis=dfile,mode='save',versionname='before_pflag',merge='replace')
    #flagmanager(vis=dfile,mode='restore',versionname='before_pflag')

    allb=[]
    for a in range(nant):
        ant=antlist[a]
        bscan=[]
        #all baselines w/ this ant
        #visdata['antenna1']
        thisant=((tbdata['ant1']==ant) & (tbdata['flag'][0]==False)) | ((tbdata['ant2']==ant) & (tbdata['flag'][0]==False))
        #pull all baselines for this ant at all times w/in tfdata
        allp=tbdata['data'][0][thisant]
        pscan=tbdata['scan'][thisant]
        for s in np.unique(tbdata['scan']):
            thisscan=(pscan==s)
            if thisscan.sum()==0: continue
            sp=np.angle(allp[thisscan],deg=True)
            msp=np.mean(sp)
            stsp=np.std(sp)
            newsp=np.array([x for x in sp if (x < msp+2.*stsp)])
            newsp=np.array([x for x in newsp if (x > msp-2.*stsp)])
            ss=pscan[thisscan]
            rms=np.sqrt(np.mean(newsp**2))
            if rms>rmsval: 
                bscan.append(s)
                allb.append(s)
            plt.scatter(s,rms)
            plt.scatter(s,rmsval,color='red')
        #plotting and saving this
        plt.savefig("./antrms/ant_rms_ant%i.png"%(ant),overwrite=True)
        badscan  = str(bscan)
        badscan  = badscan.strip('[]')
        if len(badscan)>0:
            flagdata(vis=dfile,antenna=str(ant),mode='manual',scan=badscan,flagbackup=True)
        plt.clf()

    allbscan=[]
    for s in np.unique(allb):
        thisscan=(tbdata['scan']==s)
        sa1=tbdata['ant1'][thisscan]
        sa2=tbdata['ant2'][thisscan]
        sall=np.concatenate((sa1,sa2))
        nvalid=len(np.unique(sall))
        if allb.count(s)>int(nvalid/2):
            allbscan.append(s)
    if len(allbscan)>0:
        abadscan  = str(allbscan)
        abadscan  = abadscan.strip('[]')
        flagdata(vis=dfile,mode='manual',scan=abadscan,flagbackup=True)


def calflag(dfile,pval):
    os.system('rm -rf '+dfile+'.flagversions')
    tbdata={}
    tb.open(dfile)
    tbdata['data']=tb.getcol('CPARAM')
    tbdata['data']=np.squeeze(tbdata['data'])
    tbdata['scan']=tb.getcol('SCAN_NUMBER')
    tbdata['ant1']=tb.getcol('ANTENNA1')
    tbdata['ant2']=tb.getcol('ANTENNA2')
    tbdata['flag']=tb.getcol('FLAG')
    tbdata['flag']=np.squeeze(tbdata['flag'])
    tb.close()

    allants=np.concatenate((tbdata['ant1'],tbdata['ant2']))
    antlist=np.unique(allants)

    nant=len(antlist)

    flagmanager(vis=dfile,mode='save',versionname='before_pflag',merge='replace')

    allb=[]
    for a in range(nant):
        ant=antlist[a]
        bscan=[]
        #all baselines w/ this ant
        thisant=((tbdata['ant1']==ant) & (tbdata['flag'][0]==False)) | ((tbdata['ant2']==ant) & (tbdata['flag'][0]==False))
        #pull all baselines for this ant at all times w/in tfdata
        allp=tbdata['data'][0][thisant]
        pscan=tbdata['scan'][thisant]
        for s in np.unique(tbdata['scan']):
            thisscan=(pscan==s)
            if thisscan.sum()==0: continue
            sp=np.angle(allp[thisscan],deg=True)
            sa=np.abs(allp[thisscan])
            if len(sp)>1: continue
            if (np.abs(sp)>pval)==True:
                bscan.append(s)
                allb.append(s)
            if (np.abs(sa-1.)>0.3)==True: 
                bscan.append(s)
                allb.append(s)
            plt.scatter([s]*len(sp),sp)
            plt.scatter(s,pval,color='red')
            plt.scatter(s,-1.*pval,color='red')
        #plotting and saving this
        plt.savefig("./antrms/ant_gainerr_ant%i.png"%(ant),overwrite=True)
        badscan  = str(bscan)
        badscan  = badscan.strip('[]')
        if len(badscan)>0:
            flagdata(vis=dfile,antenna=str(ant),mode='manual',scan=badscan,flagbackup=True)
        plt.clf()


##################################################################
#Setting default value for TCLEAN stopcode
stop=0

#To time the code
start=time.time()

#this is for Venki data
target='sgr_apr24'

#for image names
case='allant'
auth='Jasmin'
#refmeth='ref35'
refmeth='ref5'
date='05252023_'+str(attn)
refant=5
rfname='DA46'

#fullvis='blabla'
#msflag(dfile=fullvis,rmsval=35)
#vis60=fullvis[:-3]+'_60s.ms'
#os.system('rm -rf '+vis60[:-3]+'*')
#split(vis=fullvis,
#      outputvis=vis60,
#      timebin='60s',
#      combine='state,scan',
#      datacolumn='data',
#      keepflags=False)
#msflag(dfile=vis60,rmsval=20)

#datams2=target+'_flagcor.ms'

#flagfile=datams2[:-3]+'_noflags.ms'
#os.system('rm -rf '+flagfile)
#split(vis=datams2,outputvis=flagfile,keepflags=False,datacolumn='data')

#Initial data file name and name for channel-averaged file
#datams1=target+'_flagcor_noflags.ms'
#datams1=flagfile
#datams_noavg='sgr_apr07_flag.ms'
#datams0='sgr_apr24_200scanXXYY.ms'
#sgr_apr24_XXYY
#datams_noavg='sgr_apr24_flag.ms'
#datams_noavg='sgr_apr24_mid.ms'
#datams_noavg='sgr_apr24_cal_60s.ms'
datams_noavg='sgr_apr24_rmt_60s.ms'
#datams_noavg=vis60
#datams0='sgr_apr24_final_XXYY_rdata_avg60s.ms'


dmsprefix=datams_noavg[:9]

rdavg=dmsprefix+'_noflag.ms'
os.system('rm -rf '+rdavg+' '+rdavg+'.flagversions')
split(vis=datams_noavg,
      outputvis=rdavg,
      #timebin=tavg,
      #combine='state,scan',
      #antenna='!DV24',
      #antenna='!DV17,DA52,DA59,PM02',
      datacolumn='data',keepflags=False)
clearcal(vis=rdavg,spw='0,1,2,3',addmodel=True)


#rdavg='sgr_apr24f.ms'
#dmsprefix=rdavg[:9]
#datams_nf=rdavg

#siggobble(datams_nf,tavgint)

#rdflag=datams_noavg[:-3]+'_rflag_avg'+tavg+'.ms'
#os.system('rm -rf '+rdflag+' '+rdflag+'.flagversions')
#split(vis=rdavg,
#      outputvis=rdflag,
#      #timebin=tavg,
#      #combine='state,scan',
#      datacolumn='corrected',keepflags=False)
#clearcal(vis=rdflag,spw='0,1,2,3',addmodel=True)

#datams0=rdflag
datams0=rdavg

#Getting file name
#dmsprefix=datams0[:9]

#Setting name of file for time-averaged split file
resavg=dmsprefix+'_avgresres.ms'

#Getting file name
#dmsprefix=dmsavg[:-3]
#datams1=resavg


#Splitting and channel averaging (if desired - no we need for resres)
#os.system('rm -rf '+resavg+' '+resavg+'.flagversions')
#split(vis=datams0,
#      outputvis=resavg,
#      #timebin=tavg,
#      #combine='state,scan',
#      datacolumn='data',keepflags=False)

#mstransform(vis=datams3,outputvis=dmsavg,keepflags=False,timeaverage=True,timebin='240s',timespan='scan,state',datacolumn='data')

rawdata=datams0
clearcal(vis=rawdata,spw='0,1,2,3',addmodel=True)

rawsplit=dmsprefix+'_rawsplit.ms'
os.system('rm -rf '+rawsplit+' '+rawsplit+'.flagversions')
split(vis=datams0,
      outputvis=rawsplit,
      #timebin='0s',
      #combine='scan,state',
      datacolumn='data',keepflags=False)


#Remove all previous attempts
os.system('rm -rf '+dmsprefix+'_it*_cmodel.*')
os.system('rm -rf '+dmsprefix+'_it*_umodel.*')
os.system('rm -rf mydata_'+target+'_scan*_*_it*.dat')
os.system('rm -rf '+dmsprefix+'_ext_it*.ms')
os.system('rm -rf '+dmsprefix+'_mod_it*.ms')
os.system('rm -rf '+dmsprefix+'_ext2_it*.ms')
os.system('rm -rf '+dmsprefix+'_ext3_it*.ms')


#Bustin some ghosts
os.system('rm -rf '+dmsprefix+'_selfcal*_it*_cmodel.*')
os.system('rm -rf '+dmsprefix+'_selfcal*_it*_umodel.*')
os.system('rm -rf mydata_'+target+'_selfcal*_scan*_*_it*.dat')
os.system('rm -rf '+dmsprefix+'_selfcal*_ext_it*.ms')
os.system('rm -rf '+dmsprefix+'_selfcal*_mod_it*.ms')
os.system('rm -rf '+dmsprefix+'_selfcal*_ext2_it*.ms')


clearcal(vis=datams0,spw='0,1,2,3',addmodel=True)

tbdata={}
tb.open(datams0)
tbdata['data']=tb.getcol('DATA')
tbdata['data']=np.squeeze(tbdata['data'])
tbdata['flag']=tb.getcol('FLAG')
tbdata['flag']=np.squeeze(tbdata['flag'])
tbdata['time']=tb.getcol('TIME')
tbdata['antenna1']=tb.getcol('ANTENNA1')
tbdata['antenna2']=tb.getcol('ANTENNA2')
tbdata['sigma']=tb.getcol('SIGMA')
tb.close()

ms.open(datams0,nomodify=True)
visdata = ms.getdata(['axis_info'],ifraxis=True)
#visdata['data'] = np.squeeze(visdata['data'])
ms.close()


allants=np.concatenate((tbdata['antenna1'],tbdata['antenna2']))
antlist=np.unique(allants)

print(antlist)

#Number of polarizations
npol=tbdata['data'].shape[0]
polnames=visdata['axis_info']['corr_axis']
#npol=len(polnames)
#print("NPOL",npol)
#corr=polnames[pol]


#Calculating number of antennas and baselines
nant=int(len(antlist))
nbl=int(nant*(nant-1)/2)
#ntimes=len(visdata['axis_info']['time_axis']['MJDseconds'])
ntimes=len(np.unique(tbdata['time']))

#tt=np.array(visdata['axis_info']['time_axis']['MJDseconds'],dtype=float)
tt=np.array(np.unique(tbdata['time']))


print(nant)
print(nbl)

ELO=np.unique(tbdata['time'])


#Setting initial data file for analysis
datams1=rawsplit

#Getting scan numbers from file
#tb.open(datams1,nomodify=True)
'''
tb.open(datams0,nomodify=True)
scan = tb.getcol('SCAN_NUMBER')

#Creating array to store residuals
visdata0 = tb.getcol('DATA')
resres=np.zeros_like(visdata0)
tb.close()

#tb.open(resavg,nomodify=True)
#visdata0=tb.getcol('DATA')
#resres=np.zeros_like(visdata0)
#scanavg=tb.getcol('SCAN_NUMBER')
#tb.close()


#Print scan numbers for edification
scan = np.unique(scan)
print(scan)
'''
#scanavg = np.unique(scanavg)

resresdict={}

#Calculating RMS of data for CLEAN threshold
#ms.open(datams1,nomodify=True)
#ms.open(datams0,nomodify=True)
#sigs=ms.getdata(['sigma'])
#fsigs=np.asarray(sigs['sigma']).flatten()
fsigs=tbdata['sigma']
print("Fsig shape:",fsigs.shape)
invsigs=1./(fsigs*fsigs)
sumsigs=np.sum(invsigs)
rmst=2.*np.sqrt(1./sumsigs)
print("RMS: %e" %(rmst))
#ms.close()



#Creating and initializing data file to store UVMULTIFIT pt source flux
g=open("uvmflux_"+target+".txt","w")
g.write("it scan ext mod moderr\n")

q=open("nv_"+target+".txt","w")
q.write("")



#Start with first iteration
it=0

#refant=gf.refantfinder(antlist=antlist,goodant=antinfo['goodant'])
#iref=np.where(antlist==refant)[0][0]

#refant=35 #DV20 index
#refant=5
iref=np.where(antlist==refant)[0][0]

scandict={}

#The wily (formerly) infinite while loop
#while(1):
while it<2:
    imname2=dmsprefix+'_it'+'{0:02d}'.format(it)+'_fullpresc'
    os.system('rm -rf '+imname2+'.*')
    tclean(vis=rawdata,
           imagename=imname2,
           spw='0,1,2,3',
           datacolumn='data',
           niter=0,
           pblimit=-1,
           imsize=[300,300],
           cell='0.2arcsec') 
    plt.clf()
    tb.open(datams0,nomodify=True)
    scan = tb.getcol('SCAN_NUMBER')
    fg = tb.getcol('FLAG')
    #Creating array to store residuals
    visdata0 = tb.getcol('DATA')
    resres=np.zeros_like(visdata0)
    tb.close()
    #tb.open(resavg,nomodify=True)
    #visdata0=tb.getcol('DATA')
    #resres=np.zeros_like(visdata0)
    #scanavg=tb.getcol('SCAN_NUMBER')
    #tb.close()


    #Print scan numbers for edification
    scan = np.unique(scan)
    scandict.update(key='scanit'+str(it),value=scan)
    #print(scan)
    #rdavg=dmsprefix+'_rdata_avg'+tavg+'.ms'
    #os.system('rm -rf '+rdavg+' '+rdavg+'.flagversions')
    #split(vis=rawdata,
    #      outputvis=rdavg,
    #      timebin=tavg,
    #      combine='state,scan',
    #      datacolumn='data',keepflags=False)
    #clearcal(vis=rdavg,spw='0,1,2,3',addmodel=True)

    #datams1=rdavg

    q.write("ITERATION "+str(it)+"\n")

    #ms.open(rdavg)
    #adata=ms.getdata(['data'],ifraxis=True)
    #adata['data']=np.squeeze(adata['data'])
    #ms.close()
    #tb.open(rdavg)
    #adata=tb.getcol('DATA')
    #tb.close()

    #Declaring empty arrays for mod UVM flux and error
    muvmf=np.array([])
    muvmferr=np.array([])
    ruvmf=np.array([])

    #Creating file names for residual and model data from UVM
    datams_ext=dmsprefix+'_ext_it'+'{0:02d}'.format(it)+'.ms'
    datams_mod=dmsprefix+'_mod_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files
    os.system('rm -rf '+datams_ext)
    split(vis=datams1,outputvis=datams_ext,datacolumn='data',keepflags=False)

    #cvis.append(str(datams_mod))
    os.system('rm -rf '+datams_mod)
    split(vis=datams1,outputvis=datams_mod,datacolumn='data',keepflags=False)

    #Clearing/creating corrected and model columns
    clearcal(vis=datams_ext,spw='0,1,2,3',addmodel=True)
    clearcal(vis=datams_mod,spw='0,1,2,3',addmodel=True) 

    #Cycle through each scan
    for s in scan:
    #for s in scanavg:
        #print(s)
        #if s==359:continue
        #UVM for extended structure leftover from fit
        sc=int(s)
        myuvfit_ext1 = uvm.uvmultifit(vis=datams_ext,
                    spw='0,1,2,3',
                    #scans=[[sc]],
                    scans=[sc],
                    model=['delta'],
                    var=['0,0,p[0]'],
                    p_ini=[0.],
                    OneFitPerChannel=False,
                    column='data',
                    write='residuals') # will write residuals in the 'corrected' column, but only for the subset of fitted channels !
        # to run on the continuum collapsed cube
        ruvmf=np.append(ruvmf,myuvfit_ext1.result['Parameters'][0])
        #plt.scatter(s,myuvfit_ext1.result['Parameters'][0])

        #UVM for model visibilities
        myuvfit_mod = uvm.uvmultifit(vis=datams_mod,
                    spw='0,1,2,3',
                    #scans=[[sc]], 
                    scans=[sc],
                    model=['delta'],
                    var=['0,0,p[0]'], 
                    p_ini=[0.], 
                    OneFitPerChannel=False,
                    column='data',
                    write='model') # will write best-fit model in the 'model' column, but only for the subset of fitted channels !
        #tb.open(datams_mod)
        #md=tb.getcol('MODEL_DATA')
        #tb.close()
        #Appending results of UVM to list
        muvmf=np.append(muvmf,myuvfit_mod.result['Parameters'][0])
        muvmferr=np.append(muvmferr,myuvfit_mod.result['Uncertainties'][0])

    plt.show()
    plt.clf()
    #Pulling residuals of UVM pt source fitting for analysis
    tb.open(datams_ext)
    ext1=tb.getcol('CORRECTED_DATA')
    tb.close()

    #Name of file to split these residuals into
    datams_ext2=dmsprefix+'_ext2_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files
    #The datacolumn here means nothing, it will be replaced
    os.system('rm -rf '+datams_ext2)
    split(vis=datams_ext,outputvis=datams_ext2,datacolumn='corrected',keepflags=False)

    #Making the residuals of the first UVMULTIFIT the data column of the new file
    #4/26/23 - just changed split again???
    #tb.open(datams_ext2,nomodify=False)
    #tb.putcol('DATA', ext1)
    #tb.close()

    #Clearing/creating corrected/model columns
    clearcal(vis=datams_ext2,spw='0,1,2,3',addmodel=True)

    #Fitting point source to residuals again to make sure everything is removed
    myuvfit_ext2 = uvm.uvmultifit(vis=datams_ext2,
                             spw='0,1,2,3',
                             model=['delta'],
                             var=['0,0,p[0]'],
                             p_ini=[0],
                             OneFitPerChannel=False, # to run on the continuum collapsed cube
                             column='data',
                             write='residuals') # will write residuals in the 'corrected' column, but only for the subset of fitted channels !

    #Retrieving pt source flux and error on the measurement from the fit
    pflux=myuvfit_ext2.result['Parameters'][0]
    pfluxerr=myuvfit_ext2.result['Uncertainties'][0]
    print("PFLUX",pflux)


    #Added second point source flux to previous set of point source fluxes
    tb.open(datams_mod,nomodify=False)
    md=tb.getcol('MODEL_DATA')
    mdf=md+pflux
    tb.putcol('MODEL_DATA',mdf)
    mdn=tb.getcol('MODEL_DATA')
    madmax=np.max(mdn)
    tb.close()
    muvmf1=muvmf+pflux
    muvmferr1=muvmferr+pfluxerr

    for j in range(len(scan)):
        plt.errorbar(x=scan[j],y=muvmf1[j],yerr=muvmferr1[j],ls='none',color='blue')
    plt.savefig(target+'uvmresultsit'+str(it)+'_'+date+'.png',overwrite=True)

    #Writing all the results to a file for each scan
    for i in range(len(scan)):
        g.write(str(it)+" "+str(scan[i])+" "+str(ruvmf[i])+" "+str(muvmf1[i])+" "+str(muvmferr1[i])+"\n")


    #File where final set of residuals is the data
    datams_ext3=dmsprefix+'_ext3_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files where second 
    #round of residuals is the data column
    '''
    os.system('rm -rf '+datams_ext3+' '+datams_ext3+'.flagversions')
    split(vis=datams_ext2,outputvis=datams_ext3,datacolumn='corrected')
    '''
    #tb.open(datams_ext2)
    #ext2=tb.getcol('CORRECTED_DATA')
    #tb.close()

    #ext21=np.zeros_like(tbdata['data'])

    #ms.open(datams_ext2)
    #edata=ms.getdata(['corrected_data'],ifraxis=True)
    #ext21[:,:,:]=np.squeeze(edata['corrected_data'])
    #ms.close()
    

    #Splitting and creating the files
    #The datacolumn here means nothing, it will be replaced
    #4/26/23: changed data to corrected??
    os.system('rm -rf '+datams_ext3+' '+datams_ext3+'.flagversions')
    split(vis=datams_ext2,outputvis=datams_ext3,datacolumn='corrected',keepflags=False)
    clearcal(vis=datams_ext3,spw='0,1,2,3',addmodel=True)

    #Making the residuals of the first UVMULTIFIT the data column of the new file
    #tb.open(datams_ext3,nomodify=False)
    #tb.putcol('DATA', ext2)
    #tb.close()

    #ms.open(datams_ext3,nomodify=False)
    #nv1=ms.getdata(['data'],ifraxis=True)
    #nv1['data'][:,0,:,:]=ext21
    #ms.putdata(nv1)
    #ms.close()

    #dme3avg=datams_ext3[:-3]+'avg'+tavg+'_it'+'{0:02d}'.format(it)+'.ms'
    #Splitting and channel averaging (if desired)
    #os.system('rm -rf '+dme3avg+' '+dme3avg+'.flagversions')
    #split(vis=datams_ext3,
    #      outputvis=dme3avg,
    #      #timebin=tavg,
    #      combine='state,scan',
    #      datacolumn='data',keepflags=False)

    imname_dc=dmsprefix+'_it'+'{0:02d}'.format(it)+'_dirtyext'
    os.system('rm -rf '+imname_dc+'.*')
    tclean(vis=datams_ext3,
            field="0",
            spw='0,1,2,3',
            imagename=imname_dc,
            datacolumn='data',
            niter=0,
            pblimit=-1,
            imsize=[300,300],
            cell='0.2arcsec')
    if it==0: slint='900s' 
    if it==1: slint='900s'
    flagmanager(vis=datams_ext3,mode='save',versionname='before_pcal',merge='replace')
    flagmanager(vis=rawdata,mode='save',versionname='before_pcal',merge='replace')

    #Cleaning residuals
    imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_cleanext'
    os.system('rm -rf '+imname1+'.*')
    r2c=tclean(vis=datams_ext3,
           field="0",
           imagename=imname1,
           spw='0,1,2,3',
           datacolumn='data',
           niter=5000,
           #gain=0.3,
               #uvrange='>50m',
           threshold=rmst,
           savemodel='modelcolumn',
           pblimit=-1,
           imsize=[300,300],
           interactive=0,
           cell='0.2arcsec')

    os.system("rm -rf pcal_it"+str(it)+"*")
    gaincal(vis=datams_ext3,
        caltable="pcal_it"+str(it)+".cal",
        field="0",
        solint=slint,
        calmode="ap",
        refant=rfname,
        gaintype="G",
        combine="scan,spw",
        minsnr=1.5)
        #    minblperant=6)

    calflag(dfile='pcal_it'+str(it)+'.cal',pval=20)

    plotms(vis='pcal_it'+str(it)+'.cal',
            xaxis='time',
            yaxis='phase',
            gridrows=3,
            gridcols=3,
            iteraxis='antenna',
            #uvrange='>50klambda',
            #subplot=221,
            plotrange = [0,0,-180,180],
            plotfile='pcal_p_it'+str(it)+'.png',overwrite=True)

    plotms(vis='pcal_it'+str(it)+'.cal',
            xaxis='time',
            yaxis='amp',
            gridrows=3,
            gridcols=3,
            #uvrange='>50klambda',
            iteraxis='antenna',
            #subplot=221,
            plotrange = [0,0,0,2],
            plotfile='pcal_a_it'+str(it)+'.png',overwrite=True)

    applycal(vis=rawdata,
             field='0',
             spwmap=[0,0,0,0],
             gaintable='pcal_it'+str(it)+'.cal',
             flagbackup=False)

    applycal(vis=datams_ext3,
             field='0',
             spwmap=[0,0,0,0],
             gaintable='pcal_it'+str(it)+'.cal',
             flagbackup=False)

        #rawdatag=rawdata[:-3]+'_pcal1.ms'
        #os.system('rm -rf '+rawdatap1+' '+rawdatap1+'.flagversions')
        #split(vis=rawdata,outputvis=rawdatag,datacolumn='corrected')

    if it==5:
        flagmanager(vis=rawdata,mode='save',versionname='before_pcal2',merge='replace')
        flagmanager(vis=datams_ext3,mode='save',versionname='before_pcal2',merge='replace')

        #Cleaning residuals
        imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_p1'
        os.system('rm -rf '+imname1+'.*')
        r2c=tclean(vis=datams_ext3,
               field="0",
               imagename=imname1,
               spw='0,1,2,3',
               datacolumn='data',
               niter=5000,
               threshold=rmst,
               savemodel='modelcolumn',
               pblimit=-1,
               imsize=[300,300],
               interactive=0,
               cell='0.2arcsec')

        os.system("rm -rf pcal2_it"+str(it)+".cal")
        gaincal(vis=datams_ext3,
            caltable="pcal2_it"+str(it)+".cal",
            field="0",
            solint="60s",
            combine="scan,spw",
            calmode="p",
            refant=rfname,
            gaintype="G",
            minsnr=3.0,
            minblperant=6)

        plotms(vis='pcal2_it'+str(it)+'.cal',
                xaxis='time',
                yaxis='phase',
                gridrows=3,
                gridcols=3,
                iteraxis='baseline',
                #subplot=221,
                plotrange = [0,0,-180,180],
                plotfile='pcal2_it'+str(it)+'.png',overwrite=True)

        applycal(vis=rawdata,
                 field='0',
                 spwmap=[0,0,0,0],
                 gaintable='pcal2_it'+str(it)+'.cal',
                 flagbackup=False)

        #rawdatag=rawdata[:-3]+'_pcal2.ms'
        #os.system('rm -rf '+rawdatag+' '+rawdatag+'.flagversions')
        #split(vis=rawdata,outputvis=rawdatag,datacolumn='corrected')

    if it==2:
        flagmanager(vis=rawdata,mode='save',versionname='before_pcal3',merge='replace')
        flagmanager(vis=datams_ext3,mode='save',versionname='before_pcal3',merge='replace')

        #Cleaning residuals
        imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_p2'
        os.system('rm -rf '+imname1+'.*')
        r2c=tclean(vis=datams_ext3p2,
               field="0",
               imagename=imname1,
               spw='0,1,2,3',
               datacolumn='data',
               niter=5000,
               threshold=rmst,
               savemodel='modelcolumn',
               pblimit=-1,
               imsize=[300,300],
               interactive=0,
               cell='0.2arcsec')

        os.system("rm -rf pcal3_it"+str(it)+".cal")
        gaincal(vis=datams_ext3p2,
            caltable="pcal3_it"+str(it)+".cal",
            field="0",
            solint=int,
            calmode="p",
            refant=rfname,
            gaintype="G",
            minsnr=3.0,
            minblperant=6)

        plotms(vis='pcal3_it'+str(it)+'.cal',
                xaxis='time',
                yaxis='phase',
                gridrows=3,
                gridcols=3,
                iteraxis='antenna',
                #subplot=221,
                plotrange = [0,0,-180,180],
                plotfile='pcal3_it'+str(it)+'.png',overwrite=True)

        applycal(vis=rawdata,
                 field='0',
                 gaintable='pcal3_it'+str(it)+'.cal',
                 flagbackup=False)

        #rawdataap=rawdata[:-3]+'_ap.ms'
        #os.system('rm -rf '+rawdatag+' '+rawdatag+'.flagversions')
        #split(vis=rawdata,outputvis=rawdatag,datacolumn='corrected')
    if it==3:
        flagmanager(vis=rawdataap,mode='save',versionname='before_apcal',merge='replace')
        flagmanager(vis=datams_ext3ap,mode='save',versionname='before_apcal',merge='replace')

        #Cleaning residuals
        imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_ap'
        os.system('rm -rf '+imname1+'.*')
        r2c=tclean(vis=datams_ext3ap,
               field="0",
               imagename=imname1,
               spw='0,1,2,3',
               datacolumn='data',
               niter=5000,
               threshold=rmst,
               savemodel='modelcolumn',
               pblimit=-1,
               imsize=[300,300],
               interactive=0,
               cell='0.2arcsec')

        #Self-calibrating final cleaned residuals
        os.system("rm -rf apcal_it"+str(it)+".cal")
        gaincal(vis=datams_ext3ap,
            caltable="apcal_it"+str(it)+".cal",
            field="0",
            solint="60s",
            calmode="ap",
            refant=rfname,
            minsnr=3.0,
            minblperant=6,
            gaintype="G",
            combine="scan")
            #solnorm=True)

        plotms(vis='apcal_it'+str(it)+'.cal',
                xaxis='time',
                yaxis='amp',
                gridrows=3,
                gridcols=3,
                iteraxis='antenna',
                #subplot=221,
                plotrange = [0,0,0,0],
                plotfile='apcal__a_it'+str(it)+'.png',overwrite=True)
        plotms(vis='apcal_it'+str(it)+'.cal',
                xaxis='time',
                yaxis='phase',
                gridrows=3,
                gridcols=3,
                iteraxis='antenna',
                #subplot=221,
                plotrange = [0,0,-180,180],
                plotfile='apcal_p_it'+str(it)+'.png',overwrite=True)

        #Applying the results of GAINCAL to the data
        applycal(vis=rawdataap,
                 field="0",
                 gaintable='apcal_it'+str(it)+'.cal')
                 #interp='linear')
        applycal(vis=datams_ext3ap,
                 field="0",
                 gaintable='apcal_it'+str(it)+'.cal')
                 #interp='linear')
    
        #rawdatasc=rawdata[:-3]+'_sc.ms'
        #os.system('rm -rf '+rawdatasc+' '+rawdatasc+'.flagversions')
        #split(vis=rawdataap,outputvis=rawdatasc,datacolumn='corrected')

    
    #Dirty image after selfcal
    imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_postsc'
    os.system('rm -rf '+imname1+'.*')
    tclean(vis=datams_ext3,
           field="0",
           imagename=imname1,
           spw='0,1,2,3',
           datacolumn='corrected',
           niter=0,
           threshold=rmst,
           savemodel='modelcolumn',
           pblimit=-1,
           imsize=[300,300],
           interactive=0,
           cell='0.2arcsec')
    


    '''
    #Dirty image of gain-calibrated data
    imname3=dmsprefix+'_it'+'{0:02d}'.format(it)+'_gainappl'
    os.system('rm -rf '+imname3+'.*')
    r2c1=tclean(#vis=rdavg
        vis=rawdata,
           imagename=imname3,
           spw='0,1,2,3',
           datacolumn='corrected',
           niter=0,
           threshold=rmst,
           pblimit=-1,
           imsize=[300,300],
           interactive=0,
           cell='0.2arcsec')
    '''
    '''
    #Dirty image of x gain-calibrated residuals x CLEAN model
    imname4=dmsprefix+'_it'+'{0:02d}'.format(it)+'_gainappl_ext'
    os.system('rm -rf '+imname4+'.*')
    r2c2=tclean(vis=datams_ext3,
           imagename=imname4,
           spw='0,1,2,3',
           datacolumn='model',
           niter=0,
           threshold=rmst,
           pblimit=-1,
           imsize=[300,300],
           interactive=0,
           cell='0.2arcsec')
    '''
    
    #Split again so gain-calibrated data and residuals are the data
    dms_selfcal=dmsprefix+"_selfcal_it"+str(it)+".ms"
    dms_selfcalf=dms_selfcal+'.flagversions'
    dms_selfcal_copy=dmsprefix+"_selfcal_copy_it"+str(it)+".ms"
    dms_selfcalf_copy=dms_selfcal_copy+'.flagversions'
    #dms_selfcal_ext=dmsprefix+"_selfcal_ext_it"+str(it)+".ms"
    #dms_selfcalf_ext=dms_selfcal_ext+'.flagversions'

    os.system('rm -rf '+dms_selfcal+' '+dms_selfcalf)
    split(vis=rawdata,outputvis=dms_selfcal,datacolumn='corrected',keepflags=False)
    #split(vis=rawdatasc,outputvis=dms_selfcal,datacolumn='data',keepflags=False)


    os.system('rm -rf '+dms_selfcal_copy+' '+dms_selfcalf_copy)
    split(vis=rawdata,outputvis=dms_selfcal_copy,datacolumn='corrected',keepflags=False)
    #split(vis=rawdatasc,outputvis=dms_selfcal_copy,datacolumn='data',keepflags=False)

    #os.system('rm -rf '+dms_selfcal_ext+' '+dms_selfcalf_ext)
    #split(vis=datams_ext3,outputvis=dms_selfcal_ext,datacolumn='corrected',keepflags=False)
    

    #Seeing why TCLEAN ended
    scode=r2c['stopcode']

    #Pulling (not) gain-calibrated residuals to add to previous iterations
    #tb.open(dms_selfcal_ext)
    #tb.open(datams_ext3)
    tb.open(datams_ext3)
    resids=tb.getcol('DATA')
    tb.close()
    resres=resids
    resresdict.update({'it'+str(it):resres})

    #Pulling gain-calibrated full data
    #tb.open(dms_selfcal)
    #print("Getting",it)
    #scdata=tb.getcol('DATA')
    #tb.close()

    #Subtracting residuals from full data to hopefully get point source
    #dminres=scdata-resids

    '''
    #Dummy file for plotting
    datams_ext3_sc=dmsprefix+'_scplotting_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files
    os.system('rm -rf '+datams_ext3_sc+' '+datams_ext3_sc+'.flagversions')
    split(vis=datams_ext3,outputvis=datams_ext3_sc,datacolumn='data')

    #Making the residual-subracted data the data column of this plotting file
    #tb.open(datams_ext3_sc,nomodify=False)
    #tb.putcol('DATA',dminres)
    #tb.close()

    #Declaring empty array for UVM pt src flux and its error
    muvmp=np.array([])
    muvmperr=np.array([])


    #Creating file names for residual and model data from UVM
    datams_mod=dmsprefix+'_mod_plotting_it'+'{0:02d}'.format(it)+'.ms'

    #Repeating above
    os.system('rm -rf '+datams_mod)
    split(vis=datams_ext3_sc,outputvis=datams_mod,datacolumn='data')

    #Clearing/creating corrected and model columns
    clearcal(vis=datams_mod,spw='0,1,2,3',addmodel=True) 

    #Cycle through each scan
    for s in scan:
        #UVM for residuals of fit
        sc=int(s)

        #UVM for model visibilities
        myuvfit_mod = uvm.uvmultifit(vis=datams_mod,
                    spw='0,1,2,3',
                    scans=[sc], 
                    model=['delta'],
                    var=['0,0,p[0]'], 
                    p_ini=[0.], 
                    OneFitPerChannel=False,
                    column='data',
                    write='model') # will write best-fit model in the 'model' column, but only for the subset of fitted channels !
        muvmp=np.append(muvmp,myuvfit_mod.result['Parameters'][0])
        muvmperr=np.append(muvmperr,myuvfit_mod.result['Uncertainties'][0])

    plt.errorbar(scan,muvmp,yerr=muvmperr)
    plt.savefig('sgra_selfcal_scanbyscan_'+target+'_it'+'{0:02d}'.format(it)+'.png',overwrite=True)
    '''

    
    #Imaging gain-calibrated data
    imname2=dmsprefix+'_it'+'{0:02d}'.format(it)+'_umodel'
    os.system('rm -rf '+imname2+'.*')
    tclean(vis=rawdata,
           imagename=imname2,
           spw='0,1,2,3',
           datacolumn='corrected',
           niter=0,
           pblimit=-1,
           imsize=[300,300],
           cell='0.2arcsec')    
    

    #Set new dataset for next iteration
    datams0=dms_selfcal_copy
    rawdata=dms_selfcal
    datams1=rawdata
    clearcal(vis=rawdata,spw='0,1,2,3',addmodel=True)
    #clearcal(vis=rdavg,spw='0,1,2,3',addmodel=True)

    print("MAX:",madmax)
    if stop==2: break
    if scode!=1: stop=2

    it+=1

#Print number of iterations performed
print("Performed %i iterations" %(it))


dms_selfcal1=dmsprefix+"_selfcal_copy_it"+str(it-1)+".ms"
#dms_selfcal1_ext=dmsprefix+"_selfcal_ext_it"+str(it-1)+".ms"


#tb.open(datams_ext3sc)
#resids=tb.getcol('DATA')
#tb.close()

#tb.open(dms_selfcal1)
#print("Getting",it)
#scdata=tb.getcol('DATA')
#dminres=scdata-resids
#tb.close()

#Dummy file for plotting
datams_ext3_sc=dmsprefix+'_scplotting_final.ms'

#Splitting and creating the files
os.system('rm -rf '+datams_ext3_sc+' '+datams_ext3_sc+'.flagversions')
split(vis=dms_selfcal1,outputvis=datams_ext3_sc,datacolumn='data',keepflags=False)

#tb.open(datams_ext3_sc,nomodify=False)
#tb.putcol('DATA',dminres)
#tb.close()


muvmf=np.array([])
muvmferr=np.array([])

#Creating file names for residual and model data from UVM
datams_mod=dmsprefix+'_mod_final.ms'
datams_ext=dmsprefix+'_ext_final.ms'

#cvis.append(str(datams_mod))
os.system('rm -rf '+datams_mod)
split(vis=datams_ext3_sc,outputvis=datams_mod,datacolumn='data',keepflags=False)
os.system('rm -rf '+datams_ext)
split(vis=datams_ext3_sc,outputvis=datams_ext,datacolumn='data',keepflags=False)
clearcal(vis=datams_ext3_sc,spw='0,1,2,3',addmodel=True)

#Clearing/creating corrected and model columns
clearcal(vis=datams_mod,spw='0,1,2,3',addmodel=True) 
clearcal(vis=datams_ext,spw='0,1,2,3',addmodel=True) 

tb.open(datams0,nomodify=True)
scan = tb.getcol('SCAN_NUMBER')

#Creating array to store residuals
visdata0 = tb.getcol('DATA')
#resres=np.zeros_like(visdata0)
tb.close()
scan=np.unique(scan)
#tb.open(resavg,nomodify=True)
#visdata0=tb.getcol('DATA')
#resres=np.zeros_like(visdata0)
#scanavg=tb.getcol('SCAN_NUMBER')
#tb.close()

#Cycle through each scan
for s in scan:
    sc=int(s)

    #UVM for model visibilities
    myuvfit_mod = uvm.uvmultifit(vis=datams_mod,
                spw='0,1,2,3',
                scans=[sc], 
                model=['delta'],
                var=['0,0,p[0]'], 
                p_ini=[0.], 
                OneFitPerChannel=False,
                column='data',
                write='model') # will write best-fit model in the 'model' column, but only for the subset of fitted channels !
    muvmf=np.append(muvmf,myuvfit_mod.result['Parameters'][0])
    muvmferr=np.append(muvmferr,myuvfit_mod.result['Uncertainties'][0])

    myuvfit_ext = uvm.uvmultifit(vis=datams_ext,
                spw='0,1,2,3',
                #scans=[[sc]],
                scans=[sc],
                model=['delta'],
                var=['0,0,p[0]'],
                p_ini=[0.],
                OneFitPerChannel=False,
                column='data',
                write='residuals') # will write residuals in the 'corrected' column, but only for the subset of fitted channels !
    # to run on the continuum collapsed cube
    ruvmf=np.append(ruvmf,myuvfit_ext.result['Parameters'][0])

plt.errorbar(scan,muvmf,yerr=muvmferr)
plt.savefig('sgra_selfcal_scanbyscan_'+target+'_final_'+date+'.png')

#Creating file for extended structure around Sgr A*
datams_fin=dmsprefix+'_extsum.ms'
os.system('rm -rf '+datams_fin)
split(vis=datams_ext,outputvis=datams_fin,datacolumn='corrected')

#Taking sum of residuals and putting it as the data column for this file
#tb.open(datams_fin,nomodify=False)
#print("Getting fin")
#tb.getcol('DATA')
#print("Putting fin")
#tb.putcol('DATA',resres)
#tb.close()

#Imaging the extended structure!
imname=dmsprefix+'_finalgo'
os.system('rm -rf '+imname+'.*')
tclean(vis=datams_fin,
       imagename=imname,
       spw='0,1,2,3',
       datacolumn='data',
       niter=0,
       pblimit=-1,
       imsize=[300,300],
       cell='0.2arcsec')

imnamec=dmsprefix+'_finalgo_clean'
os.system('rm -rf '+imnamec+'.*')
tclean(vis=datams_fin,
       imagename=imnamec,
       spw='0,1,2,3',
       datacolumn='data',
       niter=5000,
       savemodel='modelcolumn',
       threshold=rmst,
       pblimit=-1,
       imsize=[300,300],
       cell='0.2arcsec')

tb.open(datams_fin)
extvis=tb.getcol('MODEL_DATA')
tb.close()

tb.open(datams_ext3_sc,nomodify=False)
fullvis=tb.getcol('DATA')
tb.putcol('CORRECTED_DATA',fullvis-extvis)
myvis=tb.getcol('CORRECTED_DATA')
myvis=np.squeeze(myvis)
myscan=tb.getcol('SCAN_NUMBER')
mywt=tb.getcol('SIGMA')
mytime=tb.getcol('TIME')
tb.close()

plt.clf()

fig, ax1 = plt.subplots()

#color = 'tab:red'

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

#color = 'tab:blue'

#plt.show()

for s in scan:
    thisscan=(myscan==s)
    allbl=myvis[0][thisscan]
    blerr=mywt[0][thisscan]
    utc=(mytime[thisscan][0]-MJDdate)/3600. #apr242018
    ax1.scatter(utc, np.abs(np.mean(allbl)), color='blue')
    ax2.scatter(utc, np.sqrt(np.sum(np.abs(blerr)**2))/len(blerr), color='orange')
    #plt.errorbar(x=utc,y=np.abs(np.mean(allbl)),yerr=np.sqrt(np.sum(np.abs(blerr)**2))/len(blerr),color='blue')
ax1.set_xlabel('Hours since 0:00 UTC')
ax1.set_ylabel('Flux density [Jy]', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')
ax1.set_title('Sgr A* - April 24')
ax2.set_ylabel('Error', color='orange')  # we already handled the x-label with ax1
ax2.tick_params(axis='y', labelcolor='orange')
#fig.tight_layout()  # otherwise the right y-label is slightly clipped


#plt.title('Sgr A* - April 24')
#plt.xlabel('Hours since 0:00 UTC')
#plt.ylabel('Flux density [Jy]')
plt.savefig('sgr_apr24_lc_Jasmin_both900_'+date+'.png',overwrite=True)
plt.show()


g.close()
'''
data=np.loadtxt("uvmflux_"+target+".txt",skiprows=1)

it=np.array(data[:,0],dtype=int)
scan=data[:,1]
modf=data[:,3]
modferr=data[:,4]

marks=['.','v','s','^','x','*','D','+','>']
cs=['blue','orange','green','purple','red','cyan','olive','pink']

last3=int(len(it)/3)
chunk=len(it)-last3

#for j in range(chunk,len(it)):
    #plt.errorbar(x=scan[j]+(0.1*it[j]),y=modf[j],yerr=modferr[j],c=cs[it[j]],marker=marks[it[j]],ls='none')
for j in range(len(it)):
    plt.errorbar(x=scan[j]+(0.1*it[j]),y=modf[j],yerr=modferr[j],c=cs[it[j]],marker=marks[it[j]],ls='none')

plt.title('Sgr A time series')
plt.xlabel('scan number')
plt.ylabel('flux [Jy]')
plt.show()
plt.savefig('sgratime_scan_'+target+'_final_pswitch_'+date+'.png')
'''

#How long this struggle bus took to run
end=time.time()
print("Time to run is",(end-start)/60,"minutes",(end-start)/3600,"hours")

