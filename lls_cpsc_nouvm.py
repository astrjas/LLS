
#from NordicARC import uvmultifit as uvm
import numpy as np
import time
import matplotlib.pyplot as plt
import os.path
from os import path
import sys
import numpy as np
path_to_gfunc='.'
sys.path.insert(0, path_to_gfunc)
import gfunc as gf

#Setting default value for TCLEAN stopcode
stop=0

#To time the code
#start=time.time()

#this is for Venki data
target='sgr_apr07'

#for image names
case='noshift_noflag'
auth='Venki'
refmeth='refanctfunc'
date='04262022'

#Initial data file name and name for channel-averaged file
datams1=target+'_flagcor_noflags.ms'

#Getting file name
dmsprefix=datams1[:-3]


#Initial data file name and name for channel-averaged file
datams1=target+'_flagcor_noflags.ms'

#Getting file name
dmsprefix=datams1[:-3]


ms.open(datams1,nomodify=True)
visdata = ms.getdata(['antenna1','antenna2','data','data_desc_id','sigma','axis_info'],ifraxis=True)
visdata['data'] = np.squeeze(visdata['data'])
ms.close()

#Number of polarizations
npol=visdata['data'].shape[0]
polnames=visdata['axis_info']['corr_axis']
#corr=polnames[pol]


'''

#Setting name of file for time-averaged split file
dmsavg=dmsprefix+'_avg.ms'


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


ms.open(datams1,nomodify=True)
visdata = ms.getdata(['antenna1','antenna2','data','data_desc_id','sigma','axis_info'],ifraxis=True)
visdata['data'] = np.squeeze(visdata['data'])
ms.close()


allants=np.concatenate((visdata['antenna1'],visdata['antenna2']))
antlist=np.unique(allants)

print(antlist)


#Calculating number of antennas and baselines
nant=int(len(antlist))
nbl=int(nant*(nant-1)/2)
ntimes=len(visdata['axis_info']['time_axis']['MJDseconds'])

print(nant)
print(nbl)

ELO=np.unique(visdata['axis_info']['time_axis']['MJDseconds'])

#Creating a tf dictionary for each ant and whether time is good or bad there
antinfo=dict()
antinfo['timestamp']=visdata['axis_info']['time_axis']['MJDseconds']

allant=np.zeros((nant,ntimes),dtype=bool)
alltim=np.zeros((ntimes,nant),dtype=bool)

########Inserted from looptest#######
#Creating a tf dictionary for each ant and whether time is good or bad there
antinfo=dict()
antinfo['timestamp']=visdata['axis_info']['time_axis']['MJDseconds']
tt=np.array(visdata['axis_info']['time_axis']['MJDseconds'],dtype=float)
#####################################


print(len(visdata['axis_info']['time_axis']['MJDseconds']))
print(len(np.unique(visdata['axis_info']['time_axis']['MJDseconds'])))



#print("num intervals: %f" %((np.max(ELO)-np.min(ELO))/240))
#print(.402600*240) 
#sys.exit()


#tsize,goodtimes=t_length(visdata=visdata,ELO=ELO,nbl=nbl)
#print("TSIZE")
#print(goodtimes)

#Splitting and channel averaging (if desired)
os.system('rm -rf '+dmsavg+' '+dmsavg+'.flagversions')
split(vis=datams1,
      outputvis=dmsavg,
      datacolumn='data')

rawdata=dmsavg
rawsplit=dmsprefix+'_rawsplit.ms'
os.system('rm -rf '+rawsplit+' '+rawsplit+'.flagversions')
split(vis=dmsavg,
      outputvis=rawsplit,
      datacolumn='data')


#Setting initial data file for analysis
datams0=rawsplit

#Getting scan numbers from file
tb.open(datams0,nomodify=True)
scan = tb.getcol('SCAN_NUMBER')

#Creating array to store residuals
visdata0 = tb.getcol('DATA')
resres=np.zeros_like(visdata0)
tb.close()

#Print scan numbers for edification
scan = np.unique(scan)
print(scan)



#Calculating RMS of data for CLEAN threshold
ms.open(datams0,nomodify=True)
sigs=ms.getdata(['sigma'])
fsigs=np.asarray(sigs['sigma']).flatten()
print("Fsig shape:",fsigs.shape)
invsigs=1./(fsigs*fsigs)
sumsigs=np.sum(invsigs)
rmst=2.*np.sqrt(1./sumsigs)
print("RMS: %e" %(rmst))
ms.close()



#Creating and initializing data file to store UVMULTIFIT pt source flux
g=open("uvmflux_"+target+".txt","w")
g.write("it scan ext mod moderr\n")
'''
'''
#######from looptest########
for ant in antlist:
    #print("ANT",ant)
    #all baselines w/ this ant
    iant=np.where(antlist==ant)[0][0]
    thisant=(visdata['antenna1']==ant) | (visdata['antenna2']==ant)
    
    #pull all baselines for this ant at all times w/in tfdata
    allt1ant=lp[thisant][:]
    #printing number of baselines w/ this ant and all times to make sure right shape
    #print("allt1ant",allt1ant.shape)
    #evaluating if baseline/antenna is bad at each time
    ant_tf=np.all(allt1ant,axis=0)
    #Adding to goodant
    allant[iant,:]=ant_tf
    #print("ant_tf",ant_tf.shape)
    plt.plot(ant_tf)
    #plotting and saving this
    plt.savefig("./lpgraphs/ant_tf_ant%i_%s_it%i.png"%(ant,corr,it),overwrite=True)
    plt.clf()
    #a+=1

antinfo['goodant']=allant
############################
'''
'''
#Start with first iteration
it=0

#refant=gf.refantfinder(antlist=antlist,goodant=antinfo['goodant'])
#iref=np.where(antlist==refant)[0][0]

#refant=35 #DV20 index


#The wily (formerly) infinite while loop
#while(1):
'''
it=0
while it<=0:
    '''
    #Declaring empty arrays for mod UVM flux and error
    muvmf=np.array([])
    muvmferr=np.array([])
    ruvmf=np.array([])

    #Creating file names for residual and model data from UVM
    datams_ext=dmsprefix+'_ext_it'+'{0:02d}'.format(it)+'.ms'
    datams_mod=dmsprefix+'_mod_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files
    os.system('rm -rf '+datams_ext)
    split(vis=datams0,outputvis=datams_ext,datacolumn='data')

    #cvis.append(str(datams_mod))
    os.system('rm -rf '+datams_mod)
    split(vis=datams0,outputvis=datams_mod,datacolumn='data')

    #Clearing/creating corrected and model columns
    clearcal(vis=datams_ext,spw='0,1,2,3',addmodel=True)
    clearcal(vis=datams_mod,spw='0,1,2,3',addmodel=True) 

    #Cycle through each scan
    for s in scan:
        #UVM for residuals of fit
        sc=int(s)
        myuvfit_ext1 = uvm.uvmultifit(vis=datams_ext,
                    spw='0,1,2,3',
                    scans=[sc],
                    model=['delta'],
                    var=['0,0,p[0]'],
                    p_ini=[0.],
                    OneFitPerChannel=False,
                    column='data',
                    write='residuals') # will write residuals in the 'corrected' column, but only for the subset of fitted channels !
        # to run on the continuum collapsed cube
        ruvmf=np.append(ruvmf,myuvfit_ext1.result['Parameters'][0])

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
        tb.open(datams_mod)
        md=tb.getcol('MODEL_DATA')
        tb.close()
        #Appending results of UVM to list
        muvmf=np.append(muvmf,myuvfit_mod.result['Parameters'][0])
        muvmferr=np.append(muvmferr,myuvfit_mod.result['Uncertainties'][0])

    #Pulling residuals of UVM pt source fitting for analysis
    tb.open(datams_ext)
    ext1=tb.getcol('CORRECTED_DATA')
    tb.close()

    #Name of file to split these residuals into
    datams_ext2=dmsprefix+'_ext2_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files
    #The datacolumn here means nothing, it will be replaced
    os.system('rm -rf '+datams_ext2)
    split(vis=datams_ext,outputvis=datams_ext2,datacolumn='data')

    #Making the residuals of the first UVMULTIFIT the data column of the new file
    tb.open(datams_ext2,nomodify=False)
    tb.putcol('DATA', ext1)
    tb.close()

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

    #Writing all the results to a file for each scan
    for i in range(len(scan)):
        g.write(str(it)+" "+str(scan[i])+" "+str(ruvmf[i])+" "+str(muvmf1[i])+" "+str(muvmferr1[i])+"\n")

    '''
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

    #Splitting and creating the files
    #The datacolumn here means nothing, it will be replaced
    #os.system('rm -rf '+datams_ext3+' '+datams_ext3+'.flagversions')
    #split(vis=datams_ext2,outputvis=datams_ext3,datacolumn='data')

    #Making the residuals of the first UVMULTIFIT the data column of the new file
    #tb.open(datams_ext3,nomodify=False)
    #tb.putcol('DATA', ext2)
    #tb.close()
    
    '''
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


    #Cleaning residuals
    imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_cmodel'
    os.system('rm -rf '+imname1+'.*')
    r2c=tclean(vis=datams_ext3,
           field="0",
           imagename=imname1,
           spw='0,1,2,3',
           datacolumn='data',
           niter=50000,
           threshold=rmst,
           savemodel='modelcolumn',
           pblimit=-1,
           imsize=[300,300],
           interactive=0,
           cell='0.2arcsec')



    '''    
    '''
    #Self-calibrating final cleaned residuals
    os.system("rm -rf phaseamp"+str(it)+".cal")
    gaincal(vis=datams_ext3,
        caltable="phaseamp"+str(it)+".cal",
        field="0",
        solint="900s",
        calmode="p",
        refant="DV13",
        gaintype="T",
        combine='scan')

    #Applying the results of GAINCAL to the data
    applycal(vis=rawdata,
             field="0",
             gaintable=['phaseamp'+str(it)+'.cal'])
             #interp='linear')
    applycal(vis=datams_ext3,
             field="0",
             gaintable=['phaseamp'+str(it)+'.cal'])
             #interp='linear')
    '''
    newvis=np.zeros_like(visdata['data'])
    antinfo=dict()
    for p in range(npol):
        newvis[p],newantinfo=gf.lls(pol=p,corr=polnames[p],datams1=datams_ext3,target=target,case=case,auth=auth,refmeth=refmeth,date=date,it=it)
        antinfo.update(newantinfo)

    tb.open(datams_ext3,nomodify=False)
    tb.putcol('CORRECTED', newvis)
    tb.close()

    tb.open(rawdata,nomodify=False)
    tb.putcol('CORRECTED', newvis)
    tb.close()


    #Dirty image of gain-calibrated data
    imname3=dmsprefix+'_it'+'{0:02d}'.format(it)+'_selfcalmodel'
    os.system('rm -rf '+imname3+'.*')
    r2c1=tclean(vis=rawdata,
           imagename=imname3,
           spw='0,1,2,3',
           datacolumn='corrected',
           niter=0,
           threshold=rmst,
           pblimit=-1,
           imsize=[300,300],
           interactive=0,
           cell='0.2arcsec')

    #Dirty image of gain-calibrated residuals
    imname4=dmsprefix+'_it'+'{0:02d}'.format(it)+'_selfcalmodel_ext'
    os.system('rm -rf '+imname4+'.*')
    r2c2=tclean(vis=datams_ext3,
           imagename=imname4,
           spw='0,1,2,3',
           datacolumn='corrected',
           niter=0,
           threshold=rmst,
           pblimit=-1,
           imsize=[300,300],
           interactive=0,
           cell='0.2arcsec')


    #Split again so gain-calibrated data and residuals are the data
    dms_selfcal=dmsprefix+"_selfcal_it"+str(it)+".ms"
    dms_selfcalf=dms_selfcal+'.flagversions'
    dms_selfcal_copy=dmsprefix+"_selfcal_copy_it"+str(it)+".ms"
    dms_selfcalf_copy=dms_selfcal_copy+'.flagversions'
    dms_selfcal_ext=dmsprefix+"_selfcal_ext_it"+str(it)+".ms"
    dms_selfcalf_ext=dms_selfcal_ext+'.flagversions'

    os.system('rm -rf '+dms_selfcal+' '+dms_selfcalf)
    split(vis=rawdata,outputvis=dms_selfcal,datacolumn='corrected')

    os.system('rm -rf '+dms_selfcal_copy+' '+dms_selfcalf_copy)
    split(vis=rawdata,outputvis=dms_selfcal_copy,datacolumn='corrected')

    os.system('rm -rf '+dms_selfcal_ext+' '+dms_selfcalf_ext)
    split(vis=datams_ext3,outputvis=dms_selfcal_ext,datacolumn='corrected')


    #Seeing why TCLEAN ended
    scode=r2c['stopcode']

    #Pulling gain-calibrated residuals to add to previous iterations
    tb.open(dms_selfcal_ext)
    resids=tb.getcol('DATA')
    tb.close()
    resres+=resids

    #Pulling gain-calibrated full data
    tb.open(dms_selfcal)
    print("Getting",it)
    scdata=tb.getcol('DATA')
    tb.close()

    #Subtracting residuals from full data to hopefully get point source
    dminres=scdata-resids

    '''
    #Dummy file for plotting
    datams_ext3_sc=dmsprefix+'_scplotting_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files
    os.system('rm -rf '+datams_ext3_sc+' '+datams_ext3_sc+'.flagversions')
    split(vis=datams_ext3,outputvis=datams_ext3_sc,datacolumn='corrected')

    #Making the residual-subracted data the data column of this plotting file
    tb.open(datams_ext3_sc,nomodify=False)
    tb.putcol('DATA',dminres)
    tb.close()

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
    plt.savefig('sgra_selfcal_scanbyscan_'+target+'_it'+'{0:02d}'.format(it)+'.png')
    '''

    '''
    #Imaging gain-calibrated data
    imname2=dmsprefix+'_it'+'{0:02d}'.format(it)+'_umodel'
    os.system('rm -rf '+imname2+'.*')
    tclean(vis=dms_selfcal,
           imagename=imname2,
           spw='0,1,2,3',
           datacolumn='data',
           niter=0,
           pblimit=-1,
           imsize=[300,300],
           cell='0.2arcsec')    
    '''

    #Set new dataset for next iteration
    datams0=dms_selfcal_copy
    rawdata=dms_selfcal

    print("MAX:",madmax)
    if stop==2: break
    if scode!=1: stop=2

    it+=1

#Print number of iterations performed
print("Performed %i iterations" %(it))


dms_selfcal1=dmsprefix+"_selfcal_copy_it2.ms"
dms_selfcal1_ext=dmsprefix+"_selfcal_ext_it2.ms"


tb.open(dms_selfcal1_ext)
resids=tb.getcol('DATA')
tb.close()

tb.open(dms_selfcal1)
print("Getting",it)
scdata=tb.getcol('DATA')
dminres=scdata-resids
tb.close()

#Dummy file for plotting
datams_ext3_sc=dmsprefix+'_scplotting_final.ms'

#Splitting and creating the files
os.system('rm -rf '+datams_ext3_sc+' '+datams_ext3_sc+'.flagversions')
split(vis=dms_selfcal1,outputvis=datams_ext3_sc,datacolumn='data')

tb.open(datams_ext3_sc,nomodify=False)
tb.putcol('DATA',dminres)
tb.close()


muvmf=np.array([])
muvmferr=np.array([])

#Creating file names for residual and model data from UVM
datams_mod=dmsprefix+'_mod_final.ms'

#cvis.append(str(datams_mod))
os.system('rm -rf '+datams_mod)
split(vis=datams_ext3_sc,outputvis=datams_mod,datacolumn='data')

#Clearing/creating corrected and model columns
clearcal(vis=datams_mod,spw='0,1,2,3',addmodel=True) 

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

plt.errorbar(scan,muvmf,yerr=muvmferr)
plt.savefig('sgra_selfcal_scanbyscan_'+target+'_final.png')

#Creating file for extended structure around Sgr A*
datams_fin=dmsprefix+'_extsum.ms'
os.system('rm -rf '+datams_fin)
split(vis=dmsavg,outputvis=datams_fin,datacolumn='data')

#Taking sum of residuals and putting it as the data column for this file
tb.open(datams_fin,nomodify=False)
print("Getting fin")
tb.getcol('DATA')
print("Putting fin")
tb.putcol('DATA',resres)
tb.close()

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
       niter=50000,
       savemodel='modelcolumn',
       threshold=rmst,
       pblimit=-1,
       imsize=[300,300],
       cell='0.2arcsec')


g.close()

data=np.loadtxt("uvmflux_"+target+".txt",skiprows=1)

it=np.array(data[:,0],dtype=int)
scan=data[:,1]
modf=data[:,3]
modferr=data[:,4]

marks=['.','v','s','^','x','*','D','+','>']
cs=['blue','orange','green','purple','red','cyan','olive','pink']

last3=int(len(it)/3)
chunk=len(it)-last3

for j in range(chunk,len(it)):
    plt.errorbar(x=scan[j]+(0.1*it[j]),y=modf[j],yerr=modferr[j],c=cs[it[j]],marker=marks[it[j]])

plt.title('Sgr A time series')
plt.xlabel('scan number')
plt.ylabel('flux [Jy]')
plt.show()
plt.savefig('sgratime_scan_'+target+'_final.png')


#How long this struggle bus took to run
end=time.time()
print("Time to run is",(end-start)/60,"minutes",(end-start)/3600,"hours")

