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
date='04182022'

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

#sweepfile=flagfile[:-3]+'_nansweep.ms'
#os.system('rm -rf '+sweepfile)
#split(vis=flagfile,outputvis=sweepfile,datacolumn='data')

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

allants=np.concatenate((visdata['antenna1'],visdata['antenna2']))
antlist=np.unique(allants)

nant=len(antlist)
nbl=int(nant*(nant-1)/2)
ntimes=len(visdata['axis_info']['time_axis']['MJDseconds'])

#Pulling all unique timestamps, print its length and then all times
ELO=np.unique(visdata['axis_info']['time_axis']['MJDseconds'])
print(len(ELO))
print(len(visdata['axis_info']['time_axis']['MJDseconds']))
#print([x-ELO[0] for x in ELO])

allant=np.zeros((nant,ntimes),dtype=bool)
alltim=np.zeros((ntimes,nant),dtype=bool)

lp=visdata['data'][0]<0.1
#lp=np.zeros_like(visdata['data'][0])
#lp=np.where(visdata['data']<0.1,np.nan,visdata['data'])
#xx[:,i]=np.where(thisant==True,np.nan,xx[:,i])

tt=np.array(visdata['axis_info']['time_axis']['MJDseconds'],dtype=float)

#nansweep=np.zeros_like(visdata['data'][0])
#lowpt=np.abs(visdata['data'][0]<0.1)
#nansweep=np.where(lowpt==True,np.nan,visdata['data'][0])
#datams1=sweepfile

#dtest=nansweep

#for y in range(741):
    #pp=np.array(visdata['data'][0,y],dtype=float)
    #dd1=np.array(dtest[y])
    #plt.scatter(tt,pp)
    #plt.scatter(tt,np.abs(dd1))
#plt.show()


#cycle through each ant
#a=0
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
    plt.savefig("./lpgraphs/ant_tf_ant%i.png"%(ant),overwrite=True)
    plt.clf()
    #a+=1

#antinfo['goodant']=allant

#lp1=np.zeros_like(lp)
lp1=np.where(lp==True,np.nan,visdata['data'])

'''
for y in range(741):
    #pp=np.array(visdata['data'][0,y],dtype=float)
    pp1=np.array(lp1[0,y],dtype=float)
    #plt.scatter(tt,pp)
    plt.scatter(tt,pp1)
'''

#Creating a tf dictionary for each ant and whether time is good or bad there
antinfo=dict()
antinfo['timestamp']=visdata['axis_info']['time_axis']['MJDseconds']

allant=np.zeros((nant,ntimes),dtype=bool)
alltim=np.zeros((ntimes,nant),dtype=bool)

#tfdict for each time and whether ant good or bad
#goodt=dict()


#TENTATIVE BEGINNING OF WHAT WILL BE IN MASSIVE LOOP
#XX correlation
xx=np.zeros_like(visdata['data'][0])
lowpt=np.abs(visdata['data'][0]<0.1)
#bd=xx0<=-1.0
xx=np.where(lowpt==True,np.nan,visdata['data'][0])
#autocorr=visdata['antenna1']==visdata['antenna2']
print("xxshape",xx.shape)

for y in range(nbl):
    #pp=np.array(visdata['data'][0,y],dtype=float)
    pp1=np.array(xx[y])
    #plt.scatter(tt,pp)
    plt.scatter(tt,np.abs(pp1))
plt.show()
plt.savefig('datatestplot_'+auth+refmeth+'_'+date+'.png',overwrite=True)
plt.clf()

#tf matrix for where ant is bad (True)
#tfdata=np.abs(visdata['data'])<=0
#tfdata=np.abs(xx)<=0
tfdata=np.isnan(xx)
print(tfdata)

#cycle through each ant
a=0
for ant in range(nant):
    #print("ANT",ant)
    #all baselines w/ this ant
    thisant=(visdata['antenna1']==ant) | (visdata['antenna2']==ant)
    #pull all baselines for this ant at all times w/in tfdata
    allt1ant=tfdata[thisant][:]
    #printing number of baselines w/ this ant and all times to make sure right shape
    #print("allt1ant",allt1ant.shape)
    #evaluating if baseline/antenna is bad at each time
    ant_tf=np.all(allt1ant,axis=0)
    #Adding to goodant
    allant[ant,:]=ant_tf
    #print("ant_tf",ant_tf.shape)
    plt.plot(ant_tf)
    #plotting and saving this
    plt.savefig("./tfgraphs/ant_tf_ant%i.png"%(ant),overwrite=True)
    plt.clf()
    a+=1

antinfo['goodant']=allant

   
for i in range(ntimes):
    j=0
    tfstore=np.zeros(len(antlist),dtype=bool)
    ct=visdata['axis_info']['time_axis']['MJDseconds'][i]
    #print("tstamp",ct)
    for ant in range(nant):
    #np.unique(visdata['antenna1']):
        tfstore[j]=allant[ant][i]
        #print(goodant[ant][i])
        j+=1
    #print(tfstore)
    #allant[0,i]=tfstore
    alltim[i,:]=tfstore
    #plt.plot(tfstore)
    #plotting and saving this
    #plt.savefig("./tfgraphs/time_tf_t%f.png"%(ELO[i]),overwrite=True)
    #plt.clf()

antinfo['goodt']=alltim

allgoodtime=0
allbad=0
good1=0

nb=0

nvis=[]

#Cycling through all the times
for time in antinfo['goodt']:
    if np.any(time)==False: allgoodtime+=1
    if np.any(time)==True:
        nbad=np.count_nonzero(time)
        if nbad==nant: allbad+=1
        if nbad<nant: good1+=1

#plt.savefig('ngood.png')
print("All good antennas:",allgoodtime)
print("Almost all good antennas:",good1)
print("All ants bad:",allbad)

refant=gf.refantfinder(antlist=antlist,goodant=antinfo['goodant'])
iref=np.where(antlist==refant)[0][0]
#refant=0


theta_r_dict={}
theta_m_dict={}
theta_Ir_dict={}
theta_del_dict={}
theta_rms_dict={}

script_L_dict={}
l_m_dict={}
l_r_dict={}
l_del_dict={}
l_rms_dict={}
l_Ir_conv_dict={}

smc=0

for i in range(len(antinfo['timestamp'])):
    #print(antinfo['timestamp'][i])
    #Jacobian and measured phase matrix
    Theta_r=np.zeros((nbl,nant-1),dtype=int)
    theta_m=np.zeros((nbl,1),dtype=float)

    #Defining Jacobian and measured amp matrix
    script_L=np.zeros((nbl,nant),dtype=int)
    l_m=np.zeros((nbl,1),dtype=float)


    #print(time)
    thistime=(visdata['axis_info']['time_axis']['MJDseconds']==antinfo['timestamp'][i])
    time=i
    '''
    gt=antinfo['goodt'][i]

    #Pulling antennas that are bad for this time
    badant=np.where(gt==True)
    if len(badant)!=0:
        for ant in badant[0]:
            #print("ANT",ant)
            thisant=(visdata['antenna1']==ant) | (visdata['antenna2']==ant)
            #replacing [thistime] with i, we'll see how it goes
            xx[:,i]=np.where(thisant==True,np.nan,xx[:,i])
    '''
    #Calculating number of antennas and baselines
    #-1 is to account for bad ants (ant1=-1 is not valid)    
    nb=0
    for ant1 in np.unique(visdata['antenna1']):
        for ant2 in np.unique(visdata['antenna2']):
            if ant1 < ant2:
                thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
                iant1=np.where(antlist==ant1)[0][0]
                #print(iant1)
                iant2=np.where(antlist==ant2)[0][0]
                #print(iant2)
                if thisbase.sum()>0:
                    #potential 0 after thisbase and thistime
                    pt=xx[thisbase][0][i]
                    ph=np.angle(pt,deg=True)
                    amp=np.absolute(pt)
                    
                    #if iant1==5 or iant2==5:
                        #amp/=2
                        #ph+=90
                        #print("bing")

                    l_m[nb,0]=np.log10(amp)
                    #anttrack[nb,0]=ant1
                    #anttrack[nb,1]=ant2

                    script_L[nb,iant1]=1
                    script_L[nb,iant2]=1


                    #PHASE STUFF
                    theta_m[nb,0]=ph

                    if ant1==refant: Theta_r[nb,iant2-1]=-1
                    if ant2==refant: Theta_r[nb,iant1]=1
                    if ant1!=refant and ant1>refant:
                        Theta_r[nb,iant1-1]=1
                        Theta_r[nb,iant2-1]=-1
                    if ant1!=refant and ant2<refant:
                        Theta_r[nb,iant1]=1
                        Theta_r[nb,iant2]=-1
                    if (ant1!=refant and (ant2>refant and ant1<refant)):
                        Theta_r[nb,iant1]=1
                        Theta_r[nb,iant2-1]=-1
                    nb+=1

    gf.dict_update(theta_r_dict,time,Theta_r)
    gf.dict_update(theta_m_dict,time,theta_m)


    gf.dict_update(script_L_dict,time,script_L)
    gf.dict_update(l_m_dict,time,l_m)


    script_L_nn=gf.nonan(script_L)
    l_m_nn=gf.nonan(l_m)

    l_r=gf.l_Ir(ll_r=script_L_nn,ll_m=l_m_nn)

    
    #print("Theta_r \n"+str(Theta_r))
    #print("theta_m \n"+str(theta_m))

    Theta_r_nn=gf.nonan(Theta_r)
    theta_m_nn=gf.nonan(theta_m)

    #if len(theta_m)==1: theta_Ir=theta_m/Theta_r
    theta_Ir=gf.th_Ir(Th_r=Theta_r_nn,th_m=theta_m_nn)
    #print("theta_Ir \n"+str(theta_Ir))


    #Residuals
    l_del=l_m_nn-np.matmul(script_L_nn,l_r)
    #print("l_del \n"+str(l_del))
    gf.dict_update(l_del_dict,time,l_del)
    #print(l_del.shape)

    l_Ir_res=gf.l_Ir(ll_r=script_L_nn,ll_m=l_del)

    #print("l_Ir w/ resids: \n"+str(l_Ir_res))

    theta_del=theta_m_nn-np.matmul(Theta_r_nn,theta_Ir)
    #print("theta_del \n"+str(theta_del))
    gf.dict_update(theta_del_dict,time,theta_del)

    theta_Ir_res=gf.th_Ir(Th_r=Theta_r_nn,th_m=theta_del)
    #print("theta_Ir w/ resids: \n"+str(theta_Ir_res))
    

    l_Ir_final=l_r+l_Ir_res
    #print(l_Ir_final.shape)
    gf.dict_update(l_r_dict,time,l_Ir_final)
    #print("final amps (log form) \n"+str(l_Ir_final))

    Ir_converted=np.zeros_like(l_Ir_final)

    #for t in range(tsize):
    #    Ir_converted[:,:,t]=10.0**(l_Ir_final[:,:,t])

    Ir_converted=np.vstack([10.0**x for x in l_Ir_final])
    #print(Ir_converted)
    if np.isnan(Ir_converted[0])==True: 
        print("BREAK!")
        #print(alla1)
        #print(alla2)
        print(script_L)
        print(l_m)
        print(l_r)
        print("ph")
        print(theta_m)
        print(Theta_r)
        #print(Theta_r)
        print(theta_Ir)
        break
        #continue

    gf.dict_update(l_Ir_conv_dict,time,Ir_converted)

    theta_Ir_final=theta_Ir+theta_Ir_res
    gf.dict_update(theta_Ir_dict,time,theta_Ir_final)


    #print("final amps (converted) \n"+str(Ir_converted))
    #print(Ir_converted.shape)

    #print("final phases \n"+str(theta_Ir_final))
    #print(theta_Ir_final.shape)



ELO_diff=np.max(ELO)-np.min(ELO)
#print(ELO_diff)
nint=int(ELO_diff/240)
ELO_range=np.linspace(np.min(ELO),np.max(ELO),nint)

#print(l_Ir_conv_dict[0])

fkeys=theta_Ir_dict.keys()
tkeys=antinfo['timestamp']

nonan=0
nemp=0
nbad=0

ndp=[]
nda=[]

avga=np.empty((nant,1,len(ELO_range)-1))
avgp=np.empty((nant,1,len(ELO_range)-1))

phase_arr=np.zeros((nant,1,ntimes),dtype=float)
amp_arr=np.zeros((nant,1,ntimes),dtype=float)


for t_step in range(ntimes):
    for k in range(nant):
        amp_arr[k,0,t_step]=l_Ir_conv_dict[t_step][k]
        if k==refant: 
            phase_arr[k,0,t_step]=0
        elif k<refant:
            phase_arr[k,0,t_step]=theta_Ir_dict[t_step][k]
        elif k>refant:
            phase_arr[k,0,t_step]=theta_Ir_dict[t_step][k-1]

antinfo['phase']=np.squeeze(phase_arr)
antinfo['amp']=np.squeeze(amp_arr)

for t_step in range(len(ELO_range)-1):
    nd=len([antinfo['timestamp'][x] for x in range(len(tkeys)) if (tkeys[x]>=ELO_range[t_step] and tkeys[x]<=ELO_range[t_step+1])])
    for k in range(nant):
        avga[k,0,t_step]=np.mean([l_Ir_conv_dict[t_step][k] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])
        if k==refant and nd!=0: 
            avgp[k,0,t_step]=0
        elif k<refant and nd!=0:
            avgp[k,0,t_step]=np.mean([theta_Ir_dict[t_step][k] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])
        elif k>refant and nd!=0:
            avgp[k,0,t_step]=np.mean([theta_Ir_dict[t_step][k-1] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])

 
    if np.isnan(np.mean(avga[:,0,t_step]))==True or np.isnan(np.mean(avgp[:,0,t_step]))==True:
        if nd==0:
            print("That's okay!")
            avga[:,0,t_step]=np.nan
            avgp[:,0,t_step]=np.nan
            nemp+=1
        else: 
            nbad+=1
    else: nonan+=1

antinfo['avg_phase']=np.squeeze(avgp)
antinfo['avg_amp']=np.squeeze(avga)    

#where comments ended
#ELO_range[:-1]
    
for a in range(nant):
    plt.scatter(ELO_range[:-1],avga[a,0,:],color='black',marker='x')
    plt.scatter(ELO_range[:-1],avgp[a,0,:],color='red')
plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_ap_'+auth+'_'+case+'.png')

plt.clf()

#plt.ylim(top=2,bottom=-2)
for a in range(nant):
    plt.scatter(ELO_range[:-1],avga[a,0,:],color='black',marker='x')
plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_a_'+auth+'_'+case+'.png')


plt.clf()

#plt.ylim(top=2,bottom=-2)
for a in range(nant):
    plt.scatter(ELO_range[:-1],avgp[a,0,:],color='red')
plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_p_'+auth+'_'+case+'.png')

plt.clf()

for a in range(nant):
    plt.scatter(antinfo['timestamp'],antinfo['amp'][a],color='black',marker='x')
plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_ampdi_'+auth+'_'+case+'.png')

plt.clf()

for a in range(nant):
    plt.scatter(antinfo['timestamp'],antinfo['phase'][a],color='red')
plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_phasedi_'+auth+'_'+case+'.png')

plt.clf()


print("Intervals with pts:",nonan)
print("Total intervals:",len(ELO_range)-1)
print("Bad baseline/antenna cases",nbad)
print("Empty bins",nemp)
print("Total points:",len(antinfo['timestamp']))

'''
newpt=np.zeros_like(l_m)

nb1=0
tb1=0

for time in ELO:
    thistime=(visdata['axis_info']['time_axis']['MJDseconds']==time)
    itime=np.where(ELO==time)[0]
    #print(itime)
    subset=visdata['data'][0]
    #print(subset.shape)
    #print("split")
    #print(visdata['data'][0,:][subset!=0].shape)
    rsize=arr_length(visdata=visdata,ELO=ELO,curr_time=time)
    #print("RSIZE")
    #print(rsize)
    if rsize<nbl:continue

    #script_L=np.zeros((rsize,nant),dtype=int)
    #l_m=np.zeros((rsize,1),dtype=float)

    #script_L=np.zeros((nbl,nant),dtype=int)
    #l_m=np.zeros((nbl,1),dtype=float)

    igtime=np.where(goodtimes==itime)
    nb1=0

    for ant1 in np.unique(visdata['antenna1']):
        for ant2 in np.unique(visdata['antenna2']):
            if ant1 < ant2:
                thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
                iant1=np.where(antlist==ant1)[0]
                iant2=np.where(antlist==ant2)[0]
                #print(iant1)
                #print(iant2)
                if thisbase.sum()>0:
                    #print("length!")
                    #print(visdata['data'][0][thisbase][0][thistime][0])
                    pt=visdata['data'][0][thisbase][0][thistime][0]
                    #print(type(pt))
                    newpt[nb1,0,tb1]=pt.real/(Ir_converted[iant1,0,igtime]*Ir_converted[iant2,0,igtime]*np.exp(1j*(phase1-phase2)))

                    nb1+=1
    tb1+=1


print(newpt)

for t in range(tsize):
    plt.scatter(range(nb1),newpt[:,0,t])

'''






#newpt=np.zeros((nbl,ntimes),dtype=float)
#oldpt=np.zeros((nbl,ntimes),dtype=float)
newpt=np.zeros_like(xx)

#oldpt=np.zeros_like(xx)
#oldpt=(
tb1=0

for t_step in range(len(ELO_range)-1):
    tchunk=[antinfo['timestamp'][x] for x in range(len(tkeys)) if (tkeys[x]>=ELO_range[t_step] and tkeys[x]<=ELO_range[t_step+1])]
    #nb1=0
    #tb1=0
    for i in range(len(tchunk)):
        #for k in range(len(vischunk)):
        thistime=(visdata['axis_info']['time_axis']['MJDseconds']==tchunk[i])
        time=i

        gt=antinfo['goodt'][i]
    
        #Pulling antennas that are bad for this time
        '''
        badant=np.where(gt==True)
        if len(badant)!=0:
            for ant in badant[0]:
                #print("ANT",ant)
                thisant=(visdata['antenna1']==ant) | (visdata['antenna2']==ant)
                #replacing [thistime] with i, we'll see how it goes
                xx[:,i]=np.where(thisant==True,np.nan,xx[:,i])
        '''
        #Calculating number of antennas and baselines
        #-1 is to account for bad ants (ant1=-1 is not valid)    
        nb1=0
        for ant1 in np.unique(visdata['antenna1']):
            for ant2 in np.unique(visdata['antenna2']):
                if ant1 < ant2:
                    thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
                    iant1=np.where(antlist==ant1)[0][0]
                    #print(iant1)
                    iant2=np.where(antlist==ant2)[0][0]
                    #print(iant2)
                    if thisbase.sum()>0:
                        #potential 0 after thisbase and thistime
                        pt=xx[thisbase][0][i]
                        oldpt[nb1,tb1]=pt
                        ga1=antinfo['avg_amp'][iant1][t_step]
                        ga2=antinfo['avg_amp'][iant2][t_step]
                        if iant1==iref: 
                            gp1=0.
                        else:
                            gp1=antinfo['avg_phase'][iant1][t_step]
                        if iant2==iref:
                            gp2=0.
                        else:
                            gp2=antinfo['avg_phase'][iant2][t_step]
                        newpt[nb1,tb1]=pt/(ga1*ga2*np.exp(1.0j*(gp1-gp2)))
                        nb1+=1
        tb1+=1

oldpt1 = np.where(np.isnan(oldpt)==True,-1.0,oldpt)
newpt1 = np.where(np.isnan(newpt)==True,0.0,newpt)

for t in range(nb1):
    #plt.scatter(antinfo['timestamp'],np.abs(newpt1[t,:]))
    oldp=np.array(xx[t])
    newp=np.array(newpt[t])
    #plt.scatter(antinfo['timestamp'],np.abs(op))
    plt.scatter(antinfo['timestamp'],np.abs(oldp))
    plt.scatter(antinfo['timestamp'],np.abs(newp))
    #plt.scatter(range(tb1),oldpt1[t,:]-newpt1[t,:])
plt.savefig("./dplots/visplot_"+auth+date+".png",overwrite=True)

#END OF DOCUMENTED UNIVERSE

