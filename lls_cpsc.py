from NordicARC import uvmultifit as uvm
import numpy as np
import time
import matplotlib.pyplot as plt
import os.path
from os import path
import sys
import numpy as np

#Setting default value for TCLEAN stopcode
stop=0

#To time the code
start=time.time()

target='sgr_apr07'

#Initial data file name and name for channel-averaged file
datams1=target+'_flagcor.ms'

#Getting file name
dmsprefix=datams1[:-3]

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

def l_Ir(ll_r,ll_m):
    print("ll_r")
    print(ll_r)
    print("ll_m")
    print(ll_m)
    #ll_r_T=np.transpose(ll_r,(1,0))
    t0=np.matmul(ll_r.T,ll_r)
    #print(t0[np.nonzero(t0)])
    t1=np.linalg.inv(np.matmul(ll_r.T,ll_r))
    print("t1")
    print(t1)
    t2=np.matmul(t1,ll_r.T)
    print("t2")
    print(t2)
    ll_Ir=np.matmul(t2,ll_m)
    print("answer!")
    print(ll_Ir)
    return ll_Ir

def arr_length(visdata,ELO,curr_time):
    thistime=(visdata['axis_info']['time_axis']['MJDseconds']==curr_time)
    rl=0
    for ant1 in np.unique(visdata['antenna1']):
        for ant2 in np.unique(visdata['antenna2']):
            if ant1 < ant2:
                thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
                if thisbase.sum()>0:
                    #print(len(thisbase==True))
                    #print(visdata['data'][0][thisbase][0][thistime][0])
                    pt=visdata['data'][0][thisbase][0][thistime][0]
                    #print(pt)
                    amp=np.absolute(pt)
                    if amp<=0:
                        #print("ah!")
                        #nb+=1
                        continue
                    else:
                        rl+=1
    return rl
    
def t_length(visdata,ELO,nbl):
    tl=0
    gt=[]
    for time in ELO:
        itime=np.where(ELO==time)[0]
        nzl=arr_length(visdata=visdata,ELO=ELO,curr_time=time)
        if nzl==nbl:
            tl+=1
            gt.append(itime)
    gt=np.hstack(gt)
    return tl,gt


ms.open(datams1,nomodify=True)
visdata = ms.getdata(['antenna1','antenna2','data','data_desc_id','sigma','axis_info'],ifraxis=True)
visdata['data'] = np.squeeze(visdata['data'])
ms.close()

allants=np.concatenate((visdata['antenna1'],visdata['antenna2']))
antlist=np.unique(allants)
#print(antlist)

#Calculating number of antennas and baselines
nant=int(len(antlist))
nbl=int(nant*(nant-1)/2)

print(nant)
print(nbl)


print(len(visdata['axis_info']['time_axis']['MJDseconds']))
print(len(np.unique(visdata['axis_info']['time_axis']['MJDseconds'])))

ELO=np.unique(visdata['axis_info']['time_axis']['MJDseconds'])

#print("num intervals: %f" %((np.max(ELO)-np.min(ELO))/240))
#print(.402600*240) 
#sys.exit()


tsize,goodtimes=t_length(visdata=visdata,ELO=ELO,nbl=nbl)
print("TSIZE")
print(goodtimes)

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


#Start with first iteration
it=0

refant=35 #DV20 index

#The wily (formerly) infinite while loop
#while(1):
while it<=0:
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


    #File where final set of residuals is the data
    datams_ext3=dmsprefix+'_ext3_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files where second 
    #round of residuals is the data column
    '''
    os.system('rm -rf '+datams_ext3+' '+datams_ext3+'.flagversions')
    split(vis=datams_ext2,outputvis=datams_ext3,datacolumn='corrected')
    '''
    tb.open(datams_ext2)
    ext2=tb.getcol('CORRECTED_DATA')
    tb.close()

    #Splitting and creating the files
    #The datacolumn here means nothing, it will be replaced
    os.system('rm -rf '+datams_ext3+' '+datams_ext3+'.flagversions')
    split(vis=datams_ext2,outputvis=datams_ext3,datacolumn='data')

    #Making the residuals of the first UVMULTIFIT the data column of the new file
    tb.open(datams_ext3,nomodify=False)
    tb.putcol('DATA', ext2)
    tb.close()

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


    #Collect data from ms
    ms.open(datams_ext3,nomodify=True)
    visdata = ms.getdata(['antenna1','antenna2','data','model_data','data_desc_id','sigma','axis_info'],ifraxis=True)
    visdata['data'] = np.squeeze(visdata['data'])
    visdata['model_data'] = np.squeeze(visdata['model_data'])
    ms.close()


    ELO_diff=np.max(ELO)-np.min(ELO)
    nint=int(ELO_diff/240)+1
    ELO_range=np.linspace(np.min(ELO),np.max(ELO),nint)

    #Jacobian and measured phase matrix
    Theta_r=np.zeros((nbl,nant-1,nint),dtype=int)
    theta_m=np.zeros((nbl,1,nint),dtype=float)

    #Defining Jacobian and measured amp matrix
    script_L=np.zeros((nbl,nant,nint),dtype=int)
    l_m=np.zeros((nbl,1,nint),dtype=float)
    anttrack=np.full(shape=(nbl,2,nint),fill_value=-1,dtype=int)


    for t_step in range(len(ELO_range)-1):
        datachunk=[x for x in visdata['data'] if (visdata['axis_info']['time_axis']['MJDseconds']>=ELO[t_step] and visdata['axis_info']['time_axis']['MJDseconds']<=ELO[t_step+1])]
        tchunk=[x for x in goodtimes if goodtimes>=ELO[t_step] and goodtimes<=ELO[t_step+1]]
        for ant1 in np.unique(visdata['antenna1']):
            for ant2 in np.unique(visdata['antenna2']):
                if ant1 < ant2:
                    thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
                    iant1=np.where(antlist==ant1)[0]
                    iant2=np.where(antlist==ant2)[0]
                    #print(iant1)
                    #print(iant2)
                    if thisbase.sum()>0:
                        print("ph!")
                        print(datachunk['data'][0][thisbase])
                        print(np.mean(datachunk['data'][0][thisbase][0]))
                        #if just good times, can use list comprehension skjekj[][x] for x in goodtimes
                        sys.exit()
                        ph=np.angle(visdata['data'][0][thisbase][0][thistime][0],deg=True)

                        anttrack[nb,0,t_step]=ant1
                        anttrack[nb,1,t_step]=ant2

                        theta_m[nb]=ph

                        if ant1==refant: Theta_r[nb,iant2-1,t_step]=-1
                        if ant2==refant: Theta_r[nb,iant1,t_step]=1
                        if ant1!=refant and ant1>refant:
                            Theta_r[nb,iant1-1,t_step]=1
                            Theta_r[nb,iant2-1,t_step]=-1
                        if ant1!=refant and ant2<refant:
                            Theta_r[nb,iant1,t_step]=1
                            Theta_r[nb,iant2,t_step]=-1
                        if (ant1!=refant and (ant2>refant and ant1<refant)):
                            Theta_r[nb,iant1,t_step]=1
                            Theta_r[nb,iant2-1,t_step]=-1

                        nb+=1
        tb+=1

        
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