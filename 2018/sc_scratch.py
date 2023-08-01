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

#outies='DA43,DA52,DA54,DA59,DV09,DV13,DV16,DV17,DV20,DV22,DV24,PM01,PM02,PM03,PM04'
exclAnts='DA43,DA52,DA54,DA59,DA65,DV09,DV10,DV20,DV24,PM01,PM02,PM03'
boundscan=[46,78,101,133,138,167,193,215,238,267,272,304,328,345,349,381,386,408,434,460,465,484,538,563,586,606,613,645,650,680]
bblocks=['46~78','101~133','138~167','193~215','238~267','272~304','328~345','349~381','386~408','434~460','465~484','538~563','586~606','613~645','650~680']

#DA43,DV17,DV09,DV22,DV24,DV20,DA63,PM01,PM02,PM03,PM04,DA52,DV04,DA48,DA65,DA54,DV16,DA59
#DA43,DA52,DA54,DA59,DV09,DV16,DV17,DV20,DV22,DV24,PM01,PM02,PM03,PM04

attn=0
tavg='60s'
tavgint=60.
#MJDdate=58233.00000000*3600.*24. #apr252018
#MJDdate=58232.00000000*3600.*24. #apr242018
#MJDdate=58230.00000000*3600.*24. #apr222018
MJDdate=58229.00000000*3600.*24. #apr212018
#MJDdate=57850.00000000 #apr072017
nspw=4


#try:
#    import uvmultimodel
#except ModuleNotFoundError:
import sys
sys.path.insert(0,'/home/jasminwash/.local/lib/python2.7/site-packages/NordicARC-3.0.5-py2.7-linux-x86_64.egg')
#import uvmultimodel


from NordicARC import uvmultifit as uvm


##################################################################
#Setting default value for TCLEAN stopcode
stop=0

#To time the code
start=time.time()

#this is for Venki data
target='sgr_apr21'

#for image names
case='innie'
auth='Jasmin'
#refmeth='ref35'
refmeth='refDV07'
date='07272023_'+str(attn)
refant=5
rfname='DA47'

fullvis='sgr_apr21_XXYY.ms'
#msflag(dfile=fullvis,rmsval=25,target=target)
#vis60=fullvis[:-3]+'_60s.ms'
#os.system('rm -rf '+vis60[:-3]+'*')
#split(vis=fullvis,
#      outputvis=vis60,
#      timebin='60s',
#      combine='state,scan',
#      datacolumn='data',
#      keepflags=False)
#msflag(dfile=vis60,rmsval=20,target=target)


#Initial data file name and name for channel-averaged file
#datams1=target+'_flagcor_noflags.ms'
#datams1=flagfile
#datams_noavg='sgr_apr07_flag.ms'
#datams0='sgr_apr24_200scanXXYY.ms'
#sgr_apr24_XXYY
#datams_noavg='sgr_apr24_flag.ms'
#datams_noavg='sgr_apr24_mid.ms'
#datams_noavg='sgr_apr24_cal_60s.ms'
#datams_noavg='sgr_apr24_rmt_60s.ms'
#datams_noavg=vis60
#datams_noavg=fullvis
#datams0='sgr_apr24_final_XXYY_rdata_avg60s.ms'


dmsprefix=fullvis[:9]

rdavg=dmsprefix+'_noflag.ms'
os.system('rm -rf '+rdavg+' '+rdavg+'.flagversions')
split(vis=fullvis,
      outputvis=rdavg,
      #timebin=tavg,
      #combine='state,scan',
      #antenna='!DV24',
      #antenna='!DV17,DA52,DA59,PM02',
      datacolumn='data',keepflags=False)
clearcal(vis=rdavg,spw='0,1,2,3',addmodel=True)

#Setting name of file for time-averaged split file
#resavg=dmsprefix+'_avgresres.ms'

#rawdata=rdavg
clearcal(vis=rdavg,spw='0,1,2,3',addmodel=True)
flagdata(vis=rdavg,scan='9~41,345,527')#9~41
flagdata(vis=rdavg,antenna='DA61,DV21')

rawsplit=dmsprefix+'_rawsplit.ms'
os.system('rm -rf '+rawsplit+' '+rawsplit+'.flagversions')
split(vis=rdavg,
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


clearcal(vis=rawsplit,spw='0,1,2,3',addmodel=True)

#tbdata={}
#tb.open(rawsplit)
#tbdata['data']=tb.getcol('DATA')
#tbdata['data']=np.squeeze(tbdata['data'])
#tbdata['flag']=tb.getcol('FLAG')
#tbdata['flag']=np.squeeze(tbdata['flag'])
#tbdata['time']=tb.getcol('TIME')
#tbdata['antenna1']=tb.getcol('ANTENNA1')
#tbdata['antenna2']=tb.getcol('ANTENNA2')
#tbdata['sigma']=tb.getcol('SIGMA')
#tb.close()

#ms.open(rawsplit,nomodify=True)
#visdata = ms.getdata(['axis_info'],ifraxis=True)
##visdata['data'] = np.squeeze(visdata['data'])
#ms.close()


#allants=np.concatenate((tbdata['antenna1'],tbdata['antenna2']))
#antlist=np.unique(allants)

#print(antlist)

#Number of polarizations
#npol=tbdata['data'].shape[0]
#polnames=visdata['axis_info']['corr_axis']


#Calculating number of antennas and baselines
#nant=int(len(antlist))
#nbl=int(nant*(nant-1)/2)
#ntimes=len(np.unique(tbdata['time']))

#tt=np.array(np.unique(tbdata['time']))


#print(nant)
#print(nbl)

#ELO=np.unique(tbdata['time'])


#Setting initial data file for analysis
datams1=rawsplit

#resresdict={}

#Calculating RMS of data for CLEAN threshold
#ms.open(datams1,nomodify=True)
#ms.open(datams0,nomodify=True)
#sigs=ms.getdata(['sigma'])
#fsigs=np.asarray(sigs['sigma']).flatten()
#fsigs=tbdata['sigma']
#print('Fsig shape:',fsigs.shape)
#invsigs=1./(fsigs*fsigs)
#sumsigs=np.sum(invsigs)
#rmst=2.*np.sqrt(1./sumsigs)
#print('RMS: %e' %(rmst))
#ms.close()



#Creating and initializing data file to store UVMULTIFIT pt source flux
#g=open('uvmflux_'+target+'.txt','w')
#g.write('it scan ext mod moderr\n')

#q=open('nv_'+target+'.txt','w')
#q.write('')



#Start with first iteration
it=0

#iref=np.where(antlist==refant)[0][0]

#scandict={}

#The wily (formerly) infinite while loop
#while(1):
#while it<1:

################
#### STEP 1 ####
################

imname2=dmsprefix+'_it'+'{0:02d}'.format(it)+'_ptdirty'
os.system('rm -rf '+imname2+'.*')
#default(tclean)
tclean(vis=datams1,
       imagename=imname2,
       spw='',
       field='0',
       uvrange='>50m',
       #datacolumn='data',
       niter=0,
       pblimit=-1,
       interactive=False,
       imsize=[300,300],
       cell='0.2arcsec') 
plt.clf()

slint='inf'
cin='scan,spw'
imax=imstat(imagename=imname2+'.image')['max'][0]
thr=imax/15.

################
#### STEP 2 ####
################

#Cleaning residuals
imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_ptclean'
os.system('rm -rf '+imname1+'*')
default(tclean)
tclean(vis=datams1,
       field='0',
       imagename=imname1,
       spw='',
       datacolumn='data',
       niter=5000,
       #gain=0.3,
       uvrange='>50m',
       threshold=thr,
       savemodel='modelcolumn',
       pblimit=-1,
       #pblimit=0.001,
       imsize=[300,300],
       cell='0.2arcsec',
       mask='ptsrc.crtf',
       interactive=False)

#flagmanager(vis=datams_ext3,mode='save',versionname='before_pcal',merge='replace')
#flagmanager(vis=rawdata,mode='save',versionname='before_pcal',merge='replace')

################
#### STEP 3 ####
################

calname='ptcal_'+target+'_it'+str(it)+'.cal'

os.system('rm -rf '+calname+'*')
gaincal(vis=datams1,
    caltable=calname,
    field='0',
    #int?
    solint=slint,
    calmode='p',
    refant=rfname,
    uvrange='>50m',
    #refantmode='strict',
    gaintype='T',
    combine=cin,
    #combine='scan',
    minsnr=1.5)
    #    minblperant=6)

#calflagT(dfile=calname,pval=50.,target=target,it=it)

plotms(vis=calname,
        xaxis='scan',
        yaxis='phase',
        gridrows=3,
        gridcols=3,
        iteraxis='antenna',
        #uvrange='>50klambda',
        #subplot=221,
        plotrange = [0,0,-180,180],
        plotfile='ptcal_p_'+target+'_it'+str(it)+'.png',overwrite=True)

################
#### STEP 4 ####
################

applycal(vis=datams1,
         field='0',
         spwmap=[0,0,0,0],
         gaintable=calname,
         applymode='calonly',
         #uvrange='>50m',
         flagbackup=True)


#Dirty image after selfcal
imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_ptpostsc'
os.system('rm -rf '+imname1+'.*')
#default(tclean)
tclean(vis=datams1,
       field='0',
       imagename=imname1,
       spw='',
       datacolumn='corrected',
       niter=0,
       #threshold=thr,
       uvrange='>50m',
       #savemodel='modelcolumn',
       pblimit=-1,
       imsize=[300,300],
       interactive=1,
       cell='0.2arcsec')
    
#Split again so gain-calibrated data and residuals are the data
dms_selfcal=dmsprefix+'_selfcal_pt1.ms'
dms_selfcalf=dms_selfcal+'.flagversions'
dms_selfcal_copy=dmsprefix+'_selfcal_copy_pt1.ms'
dms_selfcalf_copy=dms_selfcal_copy+'.flagversions'

os.system('rm -rf '+dms_selfcal+' '+dms_selfcalf)
split(vis=datams1,outputvis=dms_selfcal,datacolumn='corrected',keepflags=False)

os.system('rm -rf '+dms_selfcal_copy+' '+dms_selfcalf_copy)
split(vis=datams1,outputvis=dms_selfcal_copy,datacolumn='corrected',keepflags=False)

#Seeing why TCLEAN ended
#scode=r2c['stopcode']

#Set new dataset for next iteration
#datams0=dms_selfcal_copy
#rawdata=dms_selfcal
datams1=dms_selfcal
clearcal(vis=datams1,spw='0,1,2,3',addmodel=True)
    

################
#### STEP 5 ####
################

imname2=dmsprefix+'_it'+'{0:02d}'.format(it)+'_ptdirty1'
os.system('rm -rf '+imname2+'.*')
#default(tclean)
tclean(vis=datams1,
       imagename=imname2,
       spw='',
       field='0',
       #datacolumn='data',
       uvrange='>50m',
       niter=0,
       pblimit=-1,
       interactive=False,
       imsize=[300,300],
       cell='0.2arcsec') 
#plt.clf()

################
#### STEP 6 ####
################


imax=imstat(imagename=imname2+'.image')['max'][0]
thr=imax/15.

#Cleaning residuals
imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_ptclean1'
os.system('rm -rf '+imname1+'*')
default(tclean)
tclean(vis=datams1,
       field='0',
       imagename=imname1,
       spw='',
       datacolumn='data',
       niter=5000,
       #gain=0.3,
       uvrange='>50m',
       threshold=thr,
       savemodel='modelcolumn',
       pblimit=-1,
       #pblimit=0.001,
       imsize=[300,300],
       cell='0.2arcsec',
       mask='ptsrc.crtf',
       #usemask='auto-multithresh',
       interactive=False)

#flagmanager(vis=datams_ext3,mode='save',versionname='before_pcal',merge='replace')
#flagmanager(vis=rawdata,mode='save',versionname='before_pcal',merge='replace')

################
#### STEP 7 ####
################

calname='pt1cal_'+target+'_it'+str(it)+'.cal'
slint='inf'
cin='spw'

os.system('rm -rf '+calname+'*')
gaincal(vis=datams1,
    caltable=calname,
    field='0',
    #int?
    solint=slint,
    calmode='p',
    refant=rfname,
    uvrange='>50m',
    #refantmode='strict',
    gaintype='T',
    combine=cin,
    #combine='scan',
    minsnr=1.5)
    #    minblperant=6)

#calflagT(dfile=calname,pval=50.,target=target,it=it)

plotms(vis=calname,
        xaxis='scan',
        yaxis='phase',
        gridrows=3,
        gridcols=3,
        iteraxis='antenna',
        #uvrange='>50klambda',
        #subplot=221,
        plotrange = [0,0,-180,180],
        plotfile='ptcal1_p_'+target+'_it'+str(it)+'.png',overwrite=True)

################
#### STEP 8 ####
################

applycal(vis=datams1,
         field='0',
         spwmap=[0,0,0,0],
         gaintable=calname,
         applymode='calonly',
         #uvrange='>50m',
         flagbackup=True)


#Dirty image after selfcal
imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_ptpostsc1'
os.system('rm -rf '+imname1+'.*')
#default(tclean)
tclean(vis=datams1,
       field='0',
       imagename=imname1,
       spw='',
       datacolumn='corrected',
       niter=0,
       #threshold=thr,
       #uvrange='<50m',
       #savemodel='modelcolumn',
       pblimit=-1,
       imsize=[300,300],
       interactive=1,
       cell='0.2arcsec')
      
#Split again so gain-calibrated data and residuals are the data
dms_selfcal=dmsprefix+'_selfcal_pt2.ms'
dms_selfcalf=dms_selfcal+'.flagversions'
dms_selfcal_copy=dmsprefix+'_selfcal_copy_pt2.ms'
dms_selfcalf_copy=dms_selfcal_copy+'.flagversions'

os.system('rm -rf '+dms_selfcal+' '+dms_selfcalf)
split(vis=datams1,outputvis=dms_selfcal,datacolumn='corrected',keepflags=False)

os.system('rm -rf '+dms_selfcal_copy+' '+dms_selfcalf_copy)
split(vis=datams1,outputvis=dms_selfcal_copy,datacolumn='corrected',keepflags=False)

#Seeing why TCLEAN ended
#scode=r2c['stopcode']


#Set new dataset for next iteration
datams0=dms_selfcal_copy
#rawdata=dms_selfcal
datams1=dms_selfcal
clearcal(vis=datams1,spw='0,1,2,3',addmodel=True)
##clearcal(vis=rdavg,spw='0,1,2,3',addmodel=True)

#tb.open(datams1,nomodify=True)
#scan = tb.getcol('SCAN_NUMBER')
#fg = tb.getcol('FLAG')
#Creating array to store residuals
#visdata0 = tb.getcol('DATA')
#resres=np.zeros_like(visdata0)
#tb.close()

#Print scan numbers for edification
#scan = np.unique(scan)
##scandict.update(key='scanit'+str(it),value=scan)

#q.write('ITERATION '+str(it)+'\n')

#Declaring empty arrays for mod UVM flux and error
muvmf=np.array([])
muvmferr=np.array([])
ruvmf=np.array([])
muvmsc=np.array([])

#Creating file names for residual and model data from UVM
datams_ext=dmsprefix+'_mnsp_it'+'{0:02d}'.format(it)+'.ms'
datams_mod=dmsprefix+'_ptsrc_it'+'{0:02d}'.format(it)+'.ms'

#Splitting and creating the files
os.system('rm -rf '+datams_ext)
split(vis=datams1,outputvis=datams_ext,datacolumn='data',keepflags=False)

#cvis.append(str(datams_mod))
os.system('rm -rf '+datams_mod)
split(vis=datams1,outputvis=datams_mod,datacolumn='data',keepflags=False)

#Clearing/creating corrected and model columns
clearcal(vis=datams_ext,addmodel=True)
clearcal(vis=datams_mod,addmodel=True) 

################
#### STEP 9 ####
################

tb.open(datams1)
scans1 = tb.getcol('SCAN_NUMBER')
scan = np.unique(scans1)
field = tb.getcol('FIELD_ID')
tb.close()
cdn = field == 0 # here 0 is the field ID of the calibrator
targetScans = np.sort(np.unique(scans1[cdn]))
#targetScansIdx = np.insert( np.diff(targetScans), 0, targetScans[0] )
#cdn = targetScansIdx > 4
#endcap = np.where(cdn==True)[0]-1
#badScans = np.sort( np.append(targetScans[cdn], targetScans[endcap]) ) # the second array inside append is to select second scan from agroup
#badscan  = [b for b in badScans if b in targetScans]

#9,41

i=0
bblock=[]
while i<(len(boundscan)-1):
    if i<(len(boundscan)-2):
        ecap=boundscan[i+2]
        e=boundscan[i]
        while e<ecap:
            pscans=range(e,e+4)
            block=list([b for b in pscans if b in targetScans and b<ecap])
            if len(block)>0: bblock.append(block)
            e+=4
    else: 
        ecap=boundscan[-1]
        e=boundscan[i]
        while e<ecap:
            pscans=range(e,e+4)
            block=list([b for b in pscans if b in targetScans])
            if len(block)>0: bblock.append(block)
            e+=4
    i+=2

'''
i=0
bblock=[]
while i<(len(badscan)-1):
    if i<(len(badscan)-2):
        block=list(targetScans[np.where(targetScans==badscan[i])[0]:np.where(targetScans==badscan[i+2])[0]])
    else: block=list(targetScans[np.where(targetScans==badscan[i])[0]:])
    #print(i)
    #print(block)
    bblock.append(block)
    i+=2
'''


#################
#### STEP 10 ####
#################

#Cycle through each scan
for s in scan:
#for b in bblock:
    #print(s)
    #UVM for extended structure leftover from fit
    sc=int(s)
    #for sp in range(4):
    myuvfit_ext1 = uvm.uvmultifit(vis=datams_ext,
                spw='0,1,2,3',
                #scans=[[sc]],
                scans=[sc],
                #scans=[b],
                model=['delta'],
                var=['0,0,p[0]'],
                p_ini=[0.],
                #timewidth=4,
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
                #scans=[b],
                model=['delta'],
                var=['0,0,p[0]'], 
                p_ini=[0.], 
                OneFitPerChannel=False,
                column='data',
                write='model') # will write best-fit model in the 'model' column, but only for the subset of fitted channels !
    #Appending results of UVM to list
    #muvmf=np.append(muvmf,[myuvfit_mod.result['Parameters'][0]]*len(b))
    #muvmferr=np.append(muvmferr,[myuvfit_mod.result['Uncertainties'][0]]*len(b))
    #muvmsc=np.append(muvmsc,b)
    muvmf=np.append(muvmf,myuvfit_mod.result['Parameters'][0])
    muvmferr=np.append(muvmferr,myuvfit_mod.result['Uncertainties'][0])
    muvmsc=np.append(muvmsc,s)


#plt.show()
#plt.clf()
#Pulling residuals of UVM pt source fitting for analysis
tb.open(datams_ext)
ext1=tb.getcol('CORRECTED_DATA')
tb.close()

#################
#### STEP 11 ####
#################

#Name of file to split these residuals into
datams_ext2=dmsprefix+'_mnsp2_it'+'{0:02d}'.format(it)+'.ms'

#Splitting and creating the files
#The datacolumn here means nothing, it will be replaced
os.system('rm -rf '+datams_ext2)
split(vis=datams_ext,outputvis=datams_ext2,datacolumn='corrected',keepflags=False)


#Clearing/creating corrected/model columns
clearcal(vis=datams_ext2,spw='0,1,2,3',addmodel=True)

#################
#### STEP 12 ####
#################

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
print('PFLUX',pflux)


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

#for j in range(len(muvmsc)):
#    plt.errorbar(x=muvmsc[j],y=muvmf1[j],yerr=muvmferr1[j],ls='none',color='blue')
#plt.errorbar(x=muvmsc,y=muvmf1,yerr=muvmferr1,color='blue')
#plt.savefig(target+'uvmresultsit'+str(it)+'_'+date+'.png',overwrite=True)

#Writing all the results to a file for each scan
#for i in range(len(scan)):
#    g.write(str(it)+' '+str(scan[i])+' '+str(ruvmf[i])+' '+str(muvmf1[i])+' '+str(muvmferr1[i])+'\n')

#################
#### STEP 13 ####
#################

#File where final set of residuals is the data
datams_ext3=dmsprefix+'_mnsp3_it'+'{0:02d}'.format(it)+'.ms'


#Splitting and creating the files
#The datacolumn here means nothing, it will be replaced
#4/26/23: changed data to corrected??
os.system('rm -rf '+datams_ext3+' '+datams_ext3+'.flagversions')
split(vis=datams_ext2,outputvis=datams_ext3,datacolumn='corrected',keepflags=False)
clearcal(vis=datams_ext3,spw='0,1,2,3',addmodel=True)


#####################
#### STEPS 14-18 ####
#####################

#9~41
#ctstore=[]
for b in range(len(bblocks)):
    dataseg='mnspblock'+str(b)+'.ms'
    os.system('rm -rf '+dataseg+'*')
    split(vis=datams_ext3,outputvis=dataseg,scan=bblocks[b],keepflags=False,datacolumn='data')

    imname_dc=dmsprefix+'_it'+'{0:02d}'.format(b)+'_mnspdirty'
    os.system('rm -rf '+imname_dc+'.*')
    #default(tclean)
    tclean(vis=dataseg,
            field='0',
            spw='',
            imagename=imname_dc,
            datacolumn='data',
            niter=0,
            #scan=bblocks[b],
            uvrange='<50m',
            pblimit=-1,
            imsize=[300,300],
            interactive=False,
            cell='0.2arcsec')

    slint='inf'
    cin='spw,scan'
    imax=imstat(imagename=imname_dc+'.image')['max'][0]
    thr=imax/15.

    #Cleaning residuals
    imname1=dmsprefix+'_it'+'{0:02d}'.format(b)+'_mnspclean'
    os.system('rm -rf '+imname1+'*')
    default(tclean)
    tclean(vis=dataseg,
           field='0',
           imagename=imname1,
           spw='',
           datacolumn='data',
           niter=5000,
           #gain=0.3,
           uvrange='<50m',
           threshold=thr,
           savemodel='modelcolumn',
           pblimit=-1,
           #scan=bblocks[b],
           #pblimit=0.001,
           imsize=[300,300],
           cell='0.2arcsec',
           mask='minispiral.crtf',
           #sidelobethreshold=2.0,
           #noisethreshold=4.5,
           interactive=False)

    calname='mnspcal_'+target+'_it'+str(b)+'.cal'
    #ctstore.append(calname)

    os.system('rm -rf '+calname+'*')
    gaincal(vis=dataseg,
        caltable=calname,
        field='0',
        #int?
        solint=slint,
        calmode='a',
        refant=rfname,
        #scan=bblocks[b]
        uvrange='<50m',
        #refantmode='strict',
        gaintype='T',
        combine=cin,
        #combine='scan',
        minsnr=1.5)
        #    minblperant=6)

    #calflagT(dfile=calname,pval=50.,target=target,it=it)

    plotms(vis=calname,
            xaxis='scan',
            yaxis='amp',
            gridrows=3,
            gridcols=3,
            #uvrange='>50klambda',
            iteraxis='antenna',
            #subplot=221,
            plotrange = [0,0,0,2],
            plotfile='mnspcal_a_'+target+'_it'+str(b)+'.png',overwrite=True)

    applycal(vis=datams1,
             field='0',
             spwmap=[0,0,0,0],
             gaintable=calname,
             applymode='calonly',
             #uvrange='<50m',
             flagbackup=True)

    applycal(vis=datams_ext3,
             field='0',
             spwmap=[0,0,0,0],
             gaintable=calname,
             #uvrange='<50m',
             applymode='calonly',
             flagbackup=True)


#Dirty image after selfcal
imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_mnsppostsc'
os.system('rm -rf '+imname1+'.*')
#default(tclean)
tclean(vis=datams_ext3,
       field='0',
       imagename=imname1,
       spw='',
       datacolumn='corrected',
       niter=0,
       #threshold=thr,
       #uvrange='>50m',
       #savemodel='modelcolumn',
       pblimit=-1,
       imsize=[300,300],
       interactive=1,
       cell='0.2arcsec')

#Imaging gain-calibrated data
imname2=dmsprefix+'_it'+'{0:02d}'.format(it)+'_fullpostsc'
os.system('rm -rf '+imname2+'.*')
tclean(vis=datams1,
       spw='',
       imagename=imname2,
       #spw='0,1,2,3',
       field='0',
       datacolumn='corrected',
       niter=0,
       pblimit=-1,
       imsize=[300,300],
       cell='0.2arcsec',
       interactive=False)


#Split again so gain-calibrated data and residuals are the data
dms_selfcal=dmsprefix+'_selfcal_it'+str(it)+'.ms'
dms_selfcalf=dms_selfcal+'.flagversions'
dms_selfcal_copy=dmsprefix+'_selfcal_copy_it'+str(it)+'.ms'
dms_selfcalf_copy=dms_selfcal_copy+'.flagversions'
dm3_selfcal=dmsprefix+'_mnspselfcal_it'+str(it)+'.ms'
dm3_selfcalf=dm3_selfcal+'.flagversions'


os.system('rm -rf '+dms_selfcal+' '+dms_selfcalf)
split(vis=datams1,outputvis=dms_selfcal,datacolumn='corrected',keepflags=False)


os.system('rm -rf '+dms_selfcal_copy+' '+dms_selfcalf_copy)
split(vis=datams1,outputvis=dms_selfcal_copy,datacolumn='corrected',keepflags=False)

os.system('rm -rf '+dm3_selfcal+' '+dm3_selfcalf)
split(vis=datams_ext3,outputvis=dm3_selfcal,datacolumn='corrected',keepflags=False)



#Set new dataset for next iteration
datams0=dms_selfcal_copy
#rawdata=dms_selfcal
datams1=dms_selfcal
clearcal(vis=datams1,spw='0,1,2,3',addmodel=True)
#clearcal(vis=rdavg,spw='0,1,2,3',addmodel=True)



#Seeing why TCLEAN ended
#scode=r2c['stopcode']

#Pulling (not) gain-calibrated residuals to add to previous iterations
#tb.open(datams_ext3)
#resids=tb.getcol('DATA')
#tb.close()
#resres=resids
#resresdict.update({'it'+str(it):resres})


    

'''
applycal(vis=datams1,
         field='0',
         spwmap=[0,0,0,0],
         gaintable=calname,
         #applymode='calonly',
         #uvrange='<50m',
         flagbackup=True)

applycal(vis=datams_ext3,
         field='0',
         spwmap=[0,0,0,0],
         gaintable=calname,
         #uvrange='<50m',
         #applymode='calonly',
         flagbackup=True)

#Split again so gain-calibrated data and residuals are the data
dms_selfcal1=dmsprefix+'_selfcal1_it'+str(it)+'.ms'
dms_selfcalf1=dms_selfcal1+'.flagversions'
dms_selfcal_copy1=dmsprefix+'_selfcal1_copy_it'+str(it)+'.ms'
dms_selfcalf_copy1=dms_selfcal_copy1+'.flagversions'
dm3_selfcal1=dmsprefix+'_extselfcal1_it'+str(it)+'.ms'
dm3_selfcalf1=dm3_selfcal1+'.flagversions'


os.system('rm -rf '+dms_selfcal1+' '+dms_selfcalf1)
split(vis=datams1,outputvis=dms_selfcal1,datacolumn='corrected',keepflags=False)


os.system('rm -rf '+dms_selfcal_copy1+' '+dms_selfcalf_copy1)
split(vis=datams1,outputvis=dms_selfcal_copy1,datacolumn='corrected',keepflags=False)

os.system('rm -rf '+dm3_selfcal1+' '+dm3_selfcalf1)
split(vis=datams_ext3,outputvis=dm3_selfcal1,datacolumn='corrected',keepflags=False)

#Set new dataset for next iteration
datams01=dms_selfcal_copy1
rawdata1=dms_selfcal1
datams11=rawdata1
clearcal(vis=datams11,spw='0,1,2,3',addmodel=True)
'''


#V_s*G

print('MAX:',madmax)
#if stop==2: break
#if scode!=1: stop=2
#it+=1
    
#Print number of iterations performed
print("Performed %i iterations" %(it))

#################
#### STEP 19 ####
#################

tb.open(dm3_selfcal)
extvis=tb.getcol('DATA')
tb.close()

tb.open(datams1,nomodify=False)
fullvis=tb.getcol('DATA')
tb.putcol('CORRECTED_DATA',fullvis-extvis)
myvis=tb.getcol('CORRECTED_DATA')
myvis=np.squeeze(myvis)
#myscan1=tb.getcol('SCAN_NUMBER')
#mywt1=tb.getcol('SIGMA')
#mytime1=tb.getcol('TIME')
#myspw=tb.getcol('SPECTRAL_WINDOW_ID')
tb.close()


datams12=dmsprefix+'_vs.ms'
os.system('rm -rf '+datams12+'*')
split(vis=datams1,outputvis=datams12,datacolumn='corrected',keepflags=False)
clearcal(datams12,spw='0,1,2,3')

tb.open(datams12,nomodify=True)
#fullvis=tb.getcol('DATA')
#tb.putcol('CORRECTED_DATA',fullvis-extvis)
myvis=tb.getcol('DATA')
myvis=np.squeeze(myvis)
myscan=tb.getcol('SCAN_NUMBER')
mywt=tb.getcol('SIGMA')
mytime=tb.getcol('TIME')
#myspw=tb.getcol('SPECTRAL_WINDOW_ID')
tb.close()


#plt.clf()

fig, ax1 = plt.subplots()

#color = 'tab:red'

#ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

#color = 'tab:blue'

#plt.show()
cols=['red','blue']
cols1=['orange','green']

allbl=np.array([])
utc=np.array([])

for s in np.unique(myscan):
    thisscan=(myscan==s)
    allbl=np.append(allbl,np.mean((myvis[0][thisscan]+myvis[1][thisscan])/2.))
    utc=np.append(utc,(mytime[thisscan][0]-MJDdate)/3600.)
 
#ax2.scatter(utc, np.sqrt(np.sum(np.abs(blerr)**2))/len(blerr), color='orange')
ax1.scatter(utc, np.abs(allbl), color='blue')

    #plt.errorbar(x=utc,y=np.abs(np.mean(allbl)),yerr=np.sqrt(np.sum(np.abs(blerr)**2))/len(blerr),color='blue')
ax1.set_xlabel('Hours since 0:00 UTC')
ax1.set_ylabel('Flux density [Jy]', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')
ax1.set_title('Sgr A* - April 21')
#ax2.set_ylabel('Error', color='orange')  # we already handled the x-label with ax1
#ax2.tick_params(axis='y', labelcolor='orange')
#fig.tight_layout()  # otherwise the right y-label is slightly clipped


#plt.title('Sgr A* - April 24')
#plt.xlabel('Hours since 0:00 UTC')
#plt.ylabel('Flux density [Jy]')
plt.savefig('./curves/'+target+'_lc_Jasmin_presc_'+date+'.png',format='png',overwrite=True)
#plt.show()



#####################
#### STEPS 20-24 ####
#####################

#bblocks=['9~41','46~78','101~133','138~167','193~215','238~267','272~304','328~345','349~381','386~408','434~460','465~484','538~563','586~606','613~645','650~680']
ctstore=[]
for b in range(len(bblocks)):
    dataseg='block'+str(b)+'.ms'
    os.system('rm -rf '+dataseg+'*')
    split(vis=datams12,outputvis=dataseg,scan=bblocks[b],keepflags=False,datacolumn='data')
    imname_dc=dmsprefix+'_it'+'{0:02d}'.format(b)+'_antdirty'
    os.system('rm -rf '+imname_dc+'.*')
    #default(tclean)
    tclean(vis=dataseg,
            field='0',
            spw='',
            imagename=imname_dc,
            datacolumn='data',
            niter=0,
            #uvrange='>50m',
            antenna='DA43,DA52,DA54,DA59,DA65,DV09,DV10,DV20,DV24,PM01,PM02,PM03',
            pblimit=-1,
            imsize=[300,300],
            interactive=False,
            cell='0.2arcsec')

    bsplit=str.split(bblocks[b],'~')
    #slint=str(int((int(bsplit[1])-int(bsplit[0]))*16.12))+'s'
    slint='inf'
    cin='scan,spw'
    imax=imstat(imagename=imname_dc+'.image')['max'][0]
    thr=imax/15.

    #Cleaning residuals
    imname1=dmsprefix+'_it'+'{0:02d}'.format(b)+'_antclean'
    os.system('rm -rf '+imname1+'*')
    default(tclean)
    tclean(vis=dataseg,
           field='0',
           imagename=imname1,
           spw='',
           datacolumn='data',
           niter=5000,
           #gain=0.3,
           #uvrange='>50m',
           threshold=thr,
           antenna='DA43,DA52,DA54,DA59,DA65,DV09,DV10,DV20,DV24,PM01,PM02,PM03',
           savemodel='modelcolumn',
           pblimit=-1,
           #pblimit=0.001,
           imsize=[300,300],
           cell='0.2arcsec',
           mask='ptsrc.crtf',
           #sidelobethreshold=2.0,
           #noisethreshold=4.5,
           interactive=False)

    calname='antcal_'+target+'_it'+str(b)+'.cal'
    #ctstore.append(calname)
    os.system('rm -rf '+calname+'*')
    gaincal(vis=dataseg,
        caltable=calname,
        field='0',
        #int?
        solint=slint,
        #scan=str(bblock[b][0])+'~'+str(bblock[b][-1])
        calmode='a',
        refant=rfname,
        #uvrange='>50m',
        antenna=exclAnts,
        #refantmode='strict',
        gaintype='T',
        combine=cin,
        #combine='scan',
        minsnr=1.5)
        #    minblperant=6)

    #calflagT(dfile=calname,pval=50.,target=target,it=it)

    plotms(vis=calname,
            xaxis='scan',
            yaxis='amp',
            gridrows=3,
            gridcols=3,
            #uvrange='>50klambda',
            iteraxis='antenna',
            #subplot=221,
            plotrange = [0,0,0,2],
            plotfile='antcal_a_'+str(b)+'_'+target+'_it'+str(it)+'.png',overwrite=True)

    applycal(vis=datams12,
             field='0',
             spwmap=[0,0,0,0],
             gaintable=calname,
             antenna=exclAnts,
             applymode='calonly',
             #uvrange='<50m',
             flagbackup=True)

    #applycal(vis=datams_ext3,
    #         field='0',
    #         spwmap=[0,0,0,0],
    #         gaintable=calname,
    #         #uvrange='<50m',
    #         applymode='calonly',
    #         flagbackup=True)


#Dirty image after selfcal
imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_antpostsc'
os.system('rm -rf '+imname1+'.*')
#default(tclean)
tclean(vis=datams12,
       field='0',
       imagename=imname1,
       spw='',
       datacolumn='corrected',
       niter=0,
       #threshold=thr,
       #uvrange='>50m',
       #savemodel='modelcolumn',
       pblimit=-1,
       imsize=[300,300],
       interactive=1,
       cell='0.2arcsec')



datams13=dmsprefix+'_vs_final.ms'
os.system('rm -rf '+datams13+'*')
split(vis=datams12,outputvis=datams13,datacolumn='corrected',keepflags=False)
clearcal(datams13,spw='0,1,2,3')



#####################
#### STEPS 25-29 ####
#####################

#bblocks=['9~41','46~78','101~133','138~167','193~215','238~267','272~304','328~345','349~381','386~408','434~460','465~484','538~563','586~606','613~645','650~680']
ctstore=[]
for b in range(len(bblocks)):
    dataseg='block'+str(b)+'_1.ms'
    os.system('rm -rf '+dataseg+'*')
    split(vis=datams13,outputvis=dataseg,scan=bblocks[b],keepflags=False,datacolumn='data')
    imname_dc=dmsprefix+'_it'+'{0:02d}'.format(b)+'_antdirty1'
    os.system('rm -rf '+imname_dc+'.*')
    #default(tclean)
    tclean(vis=dataseg,
            field='0',
            spw='',
            imagename=imname_dc,
            datacolumn='data',
            niter=0,
            #uvrange='>50m',
            #antenna='DA43,DA52,DA54,DA59,DA65,DV09,DV10,DV20,DV24,PM01,PM02,PM03',
            pblimit=-1,
            imsize=[300,300],
            interactive=False,
            cell='0.2arcsec')

    bsplit=str.split(bblocks[b],'~')
    #slint=str(int((int(bsplit[1])-int(bsplit[0]))*16.))+'s'
    slint='inf'
    cin='scan,spw'
    imax=imstat(imagename=imname_dc+'.image')['max'][0]
    thr=imax/15.

    #Cleaning residuals
    imname1=dmsprefix+'_it'+'{0:02d}'.format(b)+'_antclean1'
    os.system('rm -rf '+imname1+'*')
    default(tclean)
    tclean(vis=dataseg,
           field='0',
           imagename=imname1,
           spw='',
           datacolumn='data',
           niter=5000,
           #gain=0.3,
           #uvrange='>50m',
           threshold=thr,
           #antenna='DA43,DA52,DA54,DA59,DA65,DV09,DV10,DV20,DV24,PM01,PM02,PM03',
           savemodel='modelcolumn',
           pblimit=-1,
           #pblimit=0.001,
           imsize=[300,300],
           cell='0.2arcsec',
           mask='ptsrc.crtf',
           #sidelobethreshold=2.0,
           #noisethreshold=4.5,
           interactive=False)

    calname='antcal1_'+target+'_it'+str(b)+'.cal'
    #ctstore.append(calname)
    os.system('rm -rf '+calname+'*')
    gaincal(vis=dataseg,
        caltable=calname,
        field='0',
        #int?
        solint=slint,
        #scan=str(bblock[b][0])+'~'+str(bblock[b][-1])
        calmode='a',
        refant=rfname,
        #uvrange='>50m',
        #antenna='DA43,DA52,DA54,DA59,DA65,DV09,DV10,DV20,DV24,PM01,PM02,PM03',
        #refantmode='strict',
        gaintype='T',
        combine=cin,
        #combine='scan',
        minsnr=1.5)
        #    minblperant=6)

    #calflagT(dfile=calname,pval=50.,target=target,it=it)

    plotms(vis=calname,
            xaxis='scan',
            yaxis='amp',
            gridrows=3,
            gridcols=3,
            #uvrange='>50klambda',
            iteraxis='antenna',
            #subplot=221,
            plotrange = [0,0,0,2],
            plotfile='antcal1_a_'+str(b)+'_'+target+'_it'+str(it)+'.png',overwrite=True)

    applycal(vis=datams13,
             field='0',
             spwmap=[0,0,0,0],
             gaintable=calname,
             applymode='calonly',
             #uvrange='<50m',
             flagbackup=True)

    #applycal(vis=datams_ext3,
    #         field='0',
    #         spwmap=[0,0,0,0],
    #         gaintable=calname,
    #         #uvrange='<50m',
    #         applymode='calonly',
    #         flagbackup=True)


#Dirty image after selfcal
imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_antpostsc1'
os.system('rm -rf '+imname1+'.*')
#default(tclean)
tclean(vis=datams13,
       field='0',
       imagename=imname1,
       spw='',
       datacolumn='corrected',
       niter=0,
       #threshold=thr,
       #uvrange='>50m',
       #savemodel='modelcolumn',
       pblimit=-1,
       imsize=[300,300],
       interactive=1,
       cell='0.2arcsec')

datams14=dmsprefix+'_lc.ms'
os.system('rm -rf '+datams14+'*')
split(vis=datams13,outputvis=datams14,datacolumn='corrected',keepflags=False)
clearcal(datams14,spw='0,1,2,3')



#dm313=dmsprefix+'_ext3_final.ms'
#os.system('rm -rf '+dm313+'*')
#split(vis=dm3_selfcal,outputvis=dm313,datacolumn='corrected',keepflags=False)

#tb.open(dm3_selfcal)
#extvis=tb.getcol('CORRECTED_DATA')
#tb.close()

tb.open(datams14,nomodify=False)
#fullvis=tb.getcol('DATA')
#tb.putcol('CORRECTED_DATA',fullvis-extvis)
myvis=tb.getcol('DATA')
myvis=np.squeeze(myvis)
myscan=tb.getcol('SCAN_NUMBER')
mywt=tb.getcol('SIGMA')
mytime=tb.getcol('TIME')
#myspw=tb.getcol('SPECTRAL_WINDOW_ID')
tb.close()

shortvis=dmsprefix+'_shortbl.ms'
longvis=dmsprefix+'_longbl.ms'
os.system('rm -rf '+shortvis+'*')
split(vis=datams14,outputvis=shortvis,uvrange='<50m',keepflags=False,datacolumn='data')
os.system('rm -rf '+longvis+'*')
split(vis=datams14,outputvis=longvis,uvrange='>50m',keepflags=False,datacolumn='data')

tb.open(shortvis,nomodify=False)
#fullvis=tb.getcol('DATA')
#tb.putcol('CORRECTED_DATA',fullvis-extvis)
myviss=tb.getcol('DATA')
myviss=np.squeeze(myviss)
myscans=tb.getcol('SCAN_NUMBER')
#mywt=tb.getcol('SIGMA')
mytimes=tb.getcol('TIME')
#myspw=tb.getcol('SPECTRAL_WINDOW_ID')
tb.close()

tb.open(longvis,nomodify=False)
myvisl=tb.getcol('DATA')
myvisl=np.squeeze(myvisl)
myscanl=tb.getcol('SCAN_NUMBER')
#mywt=tb.getcol('SIGMA')
mytimel=tb.getcol('TIME')
tb.close()


#################
#### STEP 30 ####
#################

rmsarr=np.array([])
for s in np.unique(myscan):
    imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_rmssnap'
    os.system('rm -rf '+imname1+'.*')
    tclean(vis=datams14,
       field='0',
       imagename=imname1,
       spw='',
       datacolumn='data',
       niter=0,
       scan=str(s),
       #threshold=thr,
       #uvrange='>50m',
       #savemodel='modelcolumn',
       pblimit=-1,
       imsize=[300,300],
       interactive=1,
       cell='0.2arcsec')

    rmsval=imstat(imagename=imname1+'.image')['rms'][0]
    rmsarr=np.append(rmsarr,rmsval)


plt.clf()

fig, ax1 = plt.subplots()

#color = 'tab:red'

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

#color = 'tab:blue'

#plt.show()
cols=['red','blue']
cols1=['orange','green']

allbl=np.array([])
allbll=np.array([])
allbls=np.array([])
utcs=np.array([])
utcl=np.array([])
utc=np.array([])

for s in np.unique(myscan):
    thisscan=(myscan==s)
    allbl=np.append(allbl,np.mean((myvis[0][thisscan]+myvis[1][thisscan])/2.))
    utc=np.append(utc,(mytime[thisscan][0]-MJDdate)/3600.)

for s in np.unique(myscans):
    thisscan=(myscans==s)
    allbls=np.append(allbls,np.mean((myviss[0][thisscan]+myviss[1][thisscan])/2.))
    utcs=np.append(utcs,(mytimes[thisscan][0]-MJDdate)/3600.)

for s in np.unique(myscanl):
    thisscan=(myscanl==s)
    allbll=np.append(allbll,np.mean((myvisl[0][thisscan]+myvisl[1][thisscan])/2.))
    utcl=np.append(utcl,(mytimel[thisscan][0]-MJDdate)/3600.)


ax2.scatter(utc, rmsarr, color='orange',facecolor='none')
ax1.scatter(utc, np.abs(allbl), color='black',label='all')
ax1.scatter(utcl, np.abs(allbll), color='blue',label='>50m',marker='^')
ax1.scatter(utcs, np.abs(allbls), color='red',label='<50m',marker='v')

#ax1.errorbar(x=utc,y=np.abs(np.mean(allbl)),yerr=rmsarr[i])
#i+=1

#plt.errorbar(x=utc,y=np.abs(np.mean(allbl)),yerr=np.sqrt(np.sum(np.abs(blerr)**2))/len(blerr),color='blue')
ax1.set_xlabel('Hours since 0:00 UTC')
ax1.set_ylabel('Flux density [Jy]', color='black')
ax1.tick_params(axis='y', labelcolor='black')
ax1.set_title('Sgr A* - April 21')
ax2.set_ylabel('Error', color='orange')  # we already handled the x-label with ax1
ax2.tick_params(axis='y', labelcolor='orange')
ax1.legend()
#fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.minorticks_on()

#plt.title('Sgr A* - April 24')
#plt.xlabel('Hours since 0:00 UTC')
#plt.ylabel('Flux density [Jy]')
plt.savefig('./curves/'+target+'_lc_Jasmin_final_'+date+'.png',format='png',overwrite=True)
plt.show()


#plt.clf()


'''
imname_dc=dmsprefix+'_it'+'{0:02d}'.format(it)+'_pt50'
os.system('rm -rf '+imname_dc+'.*')
#default(tclean)
tclean(vis=datams12,
        field='0',
        spw='',
        imagename=imname_dc,
        datacolumn='data',
        niter=0,
        uvrange='<50m',
        pblimit=-1,
        imsize=[300,300],
        interactive=False,
        cell='0.2arcsec')


imax=imstat(imagename=imname_dc+'.image')['max'][0]
thr=imax/15.


imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_cleanpt50'
os.system('rm -rf '+imname1+'*')
default(tclean)
tclean(vis=datams12,
       field='0',
       imagename=imname1,
       spw='',
       datacolumn='data',
       niter=5000,
       #gain=0.3,
       uvrange='<50m',
       threshold=thr,
       savemodel='modelcolumn',
       pblimit=-1,
       #pblimit=0.001,
       imsize=[300,300],
       cell='0.2arcsec',
       mask='ptsrc.crtf',
       #threshold=2.0,
       #noisethreshold=4.5,
       interactive=False)

calname='ptshort.cal'

os.system('rm -rf '+calname+'*')
gaincal(vis=datams12,
    caltable=calname,
    field='0',
    #int?
    solint='600s',
    calmode='a',
    #antenna='DA43,DA52,DA54,DA59,DA65,DV09,DV10,DV20,DV24,PM01,PM02,PM03',
    refant=rfname,
    uvrange='<50m',
    #refantmode='strict',
    gaintype='T',
    combine='scan,spw',
    #combine='scan',
    minsnr=1.5)
    #    minblperant=6)gaincal(

calname1='exclants.cal'

os.system('rm -rf '+calname1+'*')
gaincal(vis=datams12,
    caltable=calname1,
    field='0',
    #int?
    solint='600s',
    calmode='a',
    antenna='DA43,DA52,DA54,DA59,DA65,DV09,DV10,DV20,DV24,PM01,PM02,PM03',
    refant=rfname,
    #uvrange='<50m',
    #refantmode='strict',
    gaintype='T',
    gaintable=calname,
    combine='scan,spw',
    #combine='scan',
    minsnr=1.5)
    #    minblperant=6)gaincal(


plotms(vis=calname,
        xaxis='scan',
        yaxis='amp',
        gridrows=3,
        gridcols=3,
        #uvrange='>50klambda',
        iteraxis='antenna',
        #subplot=221,
        plotrange = [0,0,0,2],
        plotfile='ptsrccal_a.png',overwrite=True)

applycal(vis=datams1,
         field='0',
         spwmap=[0,0,0,0],
         gaintable=calname,
         #antenna='DA43,DA52,DA54,DA59,DA65,DV09,DV10,DV20,DV24,PM01,PM02,PM03',
         applymode='calonly',
         #uvrange='<50m',
         flagbackup=True)
applycal(vis=dm3_selfcal,
         field='0',
         spwmap=[0,0,0,0],
         gaintable=calname,
         #antenna='DA43,DA52,DA54,DA59,DA65,DV09,DV10,DV20,DV24,PM01,PM02,PM03',
         applymode='calonly',
         #uvrange='<50m',
         flagbackup=True)
'''




'''

#################################################################################
imname=dmsprefix+'_finalgo'
os.system('rm -rf '+imname+'.*')
tclean(vis=datams_ext3,
       imagename=imname,
       spw='0,1,2,3',
       datacolumn='data',
       niter=0,
       pblimit=-1,
       imsize=[300,300],
       cell='0.2arcsec')

imax=imstat(imagename=imname+'.image')['max'][0]
thr=imax/15.


imnamec=dmsprefix+'_finalgo_clean'
os.system('rm -rf '+imnamec+'.*')
tclean(vis=datams_ext3,
       imagename=imnamec,
       spw='0,1,2,3',
       datacolumn='data',
       niter=5000,
       savemodel='modelcolumn',
       threshold=thr,
       pblimit=-1,
       imsize=[300,300],
       cell='0.2arcsec',
       interactive=True)


tb.open(datams_ext3)
extvis=tb.getcol('MODEL_DATA')
tb.close()

tb.open(datams1,nomodify=False)
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
cols=['red','blue']
cols1=['orange','green']

for s in scan:
    thisscan=(myscan==s)
    allbl=(myvis[0][thisscan]+myvis[1][thisscan])/2
    blerr=mywt[0][thisscan]
    utc=(mytime[thisscan][0]-MJDdate)/3600. 
    ax2.scatter(utc, np.sqrt(np.sum(np.abs(blerr)**2))/len(blerr), color='orange')
    ax1.scatter(utc, np.abs(np.mean(allbl)), color='blue')

    #plt.errorbar(x=utc,y=np.abs(np.mean(allbl)),yerr=np.sqrt(np.sum(np.abs(blerr)**2))/len(blerr),color='blue')
ax1.set_xlabel('Hours since 0:00 UTC')
ax1.set_ylabel('Flux density [Jy]', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')
ax1.set_title('Sgr A* - April 21')
ax2.set_ylabel('Error', color='orange')  # we already handled the x-label with ax1
ax2.tick_params(axis='y', labelcolor='orange')
#fig.tight_layout()  # otherwise the right y-label is slightly clipped


#plt.title('Sgr A* - April 24')
#plt.xlabel('Hours since 0:00 UTC')
#plt.ylabel('Flux density [Jy]')
plt.savefig(target+'_lc_Jasmin_one600sother_'+date+'.png',overwrite=True)
plt.show()

imname=dmsprefix+'_Dantest'
os.system('rm -rf '+imname+'.*')
tclean(vis=datams1,
       imagename=imname,
       spw='0,1,2,3',
       datacolumn='corrected',
       niter=0,
       pblimit=-1,
       imsize=[300,300],
       cell='0.2arcsec')


############################################################################

dms_selfcal1=dmsprefix+"_selfcal_copy_it"+str(it-1)+".ms"

#Dummy file for plotting
datams_ext3_sc=dmsprefix+'_scplotting_final.ms'

#Splitting and creating the files
os.system('rm -rf '+datams_ext3_sc+' '+datams_ext3_sc+'.flagversions')
split(vis=dms_selfcal1,outputvis=datams_ext3_sc,datacolumn='data',keepflags=False)


muvmf=np.array([])
muvmferr=np.array([])
muvmsc=np.array([])

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
tb.close()
scan=np.unique(scan)

tb.open(datams_ext)
scans1 = tb.getcol('SCAN_NUMBER')
field = tb.getcol('FIELD_ID')
tb.close()
cdn = field == 0 # here 5 is the field ID of the calibrator
targetScans = np.sort(np.unique(scans1[cdn]))
targetScansIdx = np.insert( np.diff(targetScans), 0, targetScans[0] )
cdn = targetScansIdx > 4
endcap = np.where(cdn==True)[0]-1
badScans = np.sort( np.append(targetScans[cdn], targetScans[endcap]) ) # the second array inside append is to select second scan from agroup
badscan  = [b for b in badScans if b in targetScans]

i=0
bblock=[]
while i<(len(badscan)-1):
    if i<(len(badscan)-2):
        block=list(targetScans[np.where(targetScans==badscan[i])[0]:np.where(targetScans==badscan[i+2])[0]])
    else: block=list(targetScans[np.where(targetScans==badscan[i])[0]:])
    #print(i)
    #print(block)
    bblock.append(block)
    i+=2



#Cycle through each scan
for b in bblock:
    #sc=int(s)

    #UVM for model visibilities
    myuvfit_mod = uvm.uvmultifit(vis=datams_mod,
                spw='0,1,2,3',
                scans=[b], 
                model=['delta'],
                var=['0,0,p[0]'], 
                p_ini=[0.], 
                OneFitPerChannel=False,
                column='data',
                write='model') # will write best-fit model in the 'model' column, but only for the subset of fitted channels !
    muvmf=np.append(muvmf,myuvfit_mod.result['Parameters'][0])
    muvmferr=np.append(muvmferr,myuvfit_mod.result['Uncertainties'][0])
    muvmsc=np.append(muvmsc,b[0])

    myuvfit_ext = uvm.uvmultifit(vis=datams_ext,
                spw='0,1,2,3',
                #scans=[[sc]],
                scans=[b],
                model=['delta'],
                var=['0,0,p[0]'],
                p_ini=[0.],
                OneFitPerChannel=False,
                column='data',
                write='residuals') # will write residuals in the 'corrected' column, but only for the subset of fitted channels !
    # to run on the continuum collapsed cube
    ruvmf=np.append(ruvmf,myuvfit_ext.result['Parameters'][0])

plt.errorbar(muvmsc,muvmf,yerr=muvmferr)
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

imax=imstat(imagename=imname+'.image')['max'][0]
thr=imax/15.


imnamec=dmsprefix+'_finalgo_clean'
os.system('rm -rf '+imnamec+'.*')
tclean(vis=datams_fin,
       imagename=imnamec,
       spw='0,1,2,3',
       datacolumn='data',
       niter=5000,
       savemodel='modelcolumn',
       threshold=thr,
       pblimit=-1,
       imsize=[300,300],
       cell='0.2arcsec',
       interactive=True)

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
cols=['red','blue']
cols1=['orange','green']

for s in scan:
    thisscan=(myscan==s)
    allbl=(myvis[0][thisscan]+myvis[1][thisscan])/2
    blerr=mywt[0][thisscan]
    utc=(mytime[thisscan][0]-MJDdate)/3600. 
    ax2.scatter(utc, np.sqrt(np.sum(np.abs(blerr)**2))/len(blerr), color='orange')
    ax1.scatter(utc, np.abs(np.mean(allbl)), color='blue')

    #plt.errorbar(x=utc,y=np.abs(np.mean(allbl)),yerr=np.sqrt(np.sum(np.abs(blerr)**2))/len(blerr),color='blue')
ax1.set_xlabel('Hours since 0:00 UTC')
ax1.set_ylabel('Flux density [Jy]', color='blue')
ax1.tick_params(axis='y', labelcolor='blue')
ax1.set_title('Sgr A* - April 21')
ax2.set_ylabel('Error', color='orange')  # we already handled the x-label with ax1
ax2.tick_params(axis='y', labelcolor='orange')
#fig.tight_layout()  # otherwise the right y-label is slightly clipped


#plt.title('Sgr A* - April 24')
#plt.xlabel('Hours since 0:00 UTC')
#plt.ylabel('Flux density [Jy]')
plt.savefig(target+'_lc_Jasmin_one600s_'+date+'.png',overwrite=True)
plt.show()
'''
