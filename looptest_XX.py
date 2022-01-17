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

#Initial data file name and name for channel-averaged file
#datams1=target+'_flagcor.ms'
datams1='Sagittarius_A_star_C.ms_CALIBRATED_SELFCAL'
#datams1='uvmfit_3d77_SGRA_onechan.ms'

#Opening data and pulling necessary info
ms.open(datams1,nomodify=True)
ms.selectinit(reset=True)
visdata = ms.getdata(['antenna1','antenna2','data','axis_info'],ifraxis=True)
#printing correlations
print(visdata['axis_info']['corr_axis'])

#Squeeze data then close ms
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

#Creating a tf dictionary for each ant and whether time is good or bad there
antinfo=dict()
antinfo['timestamp']=visdata['axis_info']['time_axis']['MJDseconds']

allant=np.zeros((nant,ntimes),dtype=bool)
alltim=np.zeros((ntimes,nant),dtype=bool)

#tfdict for each time and whether ant good or bad
#goodt=dict()

#XX correlation
xx=visdata['data'][0]
#tf matrix for where ant is bad (True)
#tfdata=np.abs(visdata['data'])<=0
tfdata=np.abs(xx)<=0
print(tfdata)

#cycle through each ant
a=0
for ant in antlist:
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
    for ant in np.unique(visdata['antenna1']):
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

#plt.plot(ELO,[0]*len(ELO))
#plt.plot(ELO,[nbl]*len(ELO))
#plt.show()


#plt.savefig('ngood.png')
print("All good antennas:",allgoodtime)
print("Almost all good antennas:",good1)
print("All ants bad:",allbad)

#refant=gf.refantfinder(antlist=antlist,goodant=antinfo['goodant'])
refant=0

#END OF DOCUMENTED UNIVERSE

#ant1dict,ant2dict,ant_good_v=gf.antdicts(visdata,ELO)
#tkeys=ant1dict.keys()
#print("number of ants",len(ant1dict[tkeys[0]]))
#print("number of times",len(tkeys))
#print("number of times again",len(ELO))
#thistime=(visdata['axis_info']['time_axis']['MJDseconds']==ELO[0])
#print(visdata['data'][0].shape)
#print(len(ant1dict))
#print(5672/1418)

#print(sbc)

#print(antlist)

#Calculating number of antennas and baselines
#nant=int(len(antlist))

#print(nant)
#print(nbl)

#plt.close()

#for x in range(len(ant1dict)):
#    print(ant1dict[x])
#    print(ant2dict[x])

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

    #Jacobian and measured phase matrix
    Theta_r=np.zeros((nbl,nant-1),dtype=int)
    theta_m=np.zeros((nbl,1),dtype=float)

    #Defining Jacobian and measured amp matrix
    script_L=np.zeros((nbl,nant),dtype=int)
    l_m=np.zeros((nbl,1),dtype=float)


    #print(time)
    #if (time in tkeys)==False:continue
    thistime=(visdata['axis_info']['time_axis']['MJDseconds']==antinfo['timestamp'][i])
    time=antinfo['timestamp'][i]

    #alla1=ant1dict[time]
    #alla2=ant2dict[time]

    #print(alla1)
    #print(alla2)

    #print(visdata['antenna1'])
    #print(visdata['antenna2'])

    gt=antinfo['goodt'][i]

    #Pulling antennas that are bad for this time
    badant=np.where(gt==True)
    #print(badant)
    #print(badbase)
    if len(badant)!=0:
        for ant in badant[0]:
            #print("ANT",ant)
            thisant=(visdata['antenna1']==ant) | (visdata['antenna2']==ant)
            #replacing [thistime] with i, we'll see how it goes
            xx[:,i]=np.where(thisant==True,np.nan,xx[:,i])

    #xxnonan=xx[~np.isnan(xx)]

    #allants=np.concatenate((alla1,alla2))
    #antlist=np.unique(allants)
    #print(alla1)
    #print(alla2)


    #Calculating number of antennas and baselines
    #-1 is to account for bad ants (ant1=-1 is not valid)
    #nant=int(len(antlist))
    #if len(np.unique(alla2))==1: continue
    #nbl=int(nant*(nant-1)/2)
    #print(nbl)
    
    nb=0
    for ant1 in np.unique(visdata['antenna1']):
        for ant2 in np.unique(visdata['antenna2']):
            if ant1 < ant2:
                thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
                iant1=np.where(antlist==ant1)[0]
                #print(iant1)
                iant2=np.where(antlist==ant2)[0]
                #print(iant2)
                if thisbase.sum()>0:
                    #potential 0 after thisbase and thistime
                    pt=xx[thisbase][0][i]
                    ph=np.angle(pt,deg=True)
                    amp=np.absolute(pt)
                    
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

    #here! yes here!    
    #print(script_L)
    #script_L_f=np.vstack([script_L[x] for x in range(len(script_L)) if len(np.unique(script_L[x]))>1])
    #l_m_f=np.vstack([l_m[x] for x in range(len(l_m)) if len(np.unique(script_L[x]))>1])

    #Theta_r=np.vstack([Theta_r[x] for x in range(len(Theta_r)) if len(np.unique(Theta_r[x]))>1])
    #print(Theta_r.shape)
    #print(Theta_r_f)
    #theta_m_f=np.vstack([theta_m[x] for x in range(len(theta_m)) if len(np.unique(Theta_r[x]))>1])
    #print(theta_m.shape)
    #print(theta_m_f)

    #if np.linalg.det(np.matmul(script_L.T,script_L))==0 or np.linalg.det(np.matmul(Theta_r.T,Theta_r))==0:
        #smc+=1
        #continue


    gf.dict_update(theta_r_dict,time,Theta_r)
    gf.dict_update(theta_m_dict,time,theta_m)


    gf.dict_update(script_L_dict,time,script_L)
    gf.dict_update(l_m_dict,time,l_m)

    #print("script_L \n"+str(script_L))
    #print("l_m \n"+str(l_m))

    #if len(l_m)==1: continue

    script_L_nn=gf.nonan(script_L)
    l_m_nn=gf.nonan(l_m)

    #l_r=np.dstack([gf.l_Ir(ll_r=script_L_nn,ll_m=l_m_nn[:,:,x]) for x in range(ntimes)])

    l_r=gf.l_Ir(ll_r=script_L_nn,ll_m=l_m_nn)
    #print("l_r \n"+str(l_r))

    
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
    print(Ir_converted)
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



for t_step in range(len(ELO_range)-1):
    datachunk_ph=np.mean([theta_Ir_dict[tkeys[x]] for x in range(len(tkeys)) if (tkeys[x]>=ELO_range[t_step] and tkeys[x]<=ELO_range[t_step+1])])
    ndp.append(datachunk_ph)
    datachunk_amp=np.mean([l_Ir_conv_dict[tkeys[x]] for x in range(len(tkeys)) if (tkeys[x]>=ELO_range[t_step] and tkeys[x]<=ELO_range[t_step+1])])
    nda.append(datachunk_amp)

    nd=len([antinfo['timestamp'][x] for x in range(len(tkeys)) if (tkeys[x]>=ELO_range[t_step] and tkeys[x]<=ELO_range[t_step+1])])

 
    if datachunk_ph==0 or datachunk_amp==0:
        if nd==0:
            print("That's okay!")
            nemp+=1
        else: 
            nbad+=1
    else: nonan+=1



    print(datachunk_ph)
    print(datachunk_amp)


#where comments ended
'''
for f in tkeys:
    datachunk_ph=len(theta_Ir_dict[f])
    ndp.append(datachunk_ph)
    datachunk_amp=len(l_Ir_conv_dict[f])
    nda.append(datachunk_amp)
    if datachunk_ph==0 or datachunk_amp==0: continue 
    nonan+=1
'''
#ELO_range[:-1]    
plt.scatter(ELO_range[:-1],nda,label='amp',color='black',marker='x')
plt.scatter(ELO_range[:-1],ndp,label='phase',color='red')
#plt.scatter(tkeys,nvis,label='vis',color='green',marker='*')
#plt.scatter(ELO_range,[nant]*len(ELO_range),label='vis',color='blue')
#plt.xlim(right=ELO_range[50])
plt.legend()
plt.show()

plt.savefig('ndpv_01132021_refant0.png')

print("Intervals with pts:",nonan)
print("Total intervals:",len(ELO_range)-1)
print("Bad baseline/antenna cases",nbad)
print("Empty bins",nemp)
print("Total points:",len(ELO))


'''
       if k>=ELO_range[t_step] and k<=ELO_range[t_step+1]:
            thistime=(visdata['axis_info']['time_axis']['MJDseconds']==k)
            print("tee!")
            #print(thistime)
            nv=visdata['data'][0][thistime][0]
            #['antenna1']
            print(nv)
            #print(len(nv))
            #nvs+=len(nv)
            nvs+=1
for time in ELO:
    thistime=(visdata['axis_info']['time_axis']['MJDseconds']==time)
    vtime=visdata['axis_info']['time_axis']['MJDseconds']
    itime=np.where(ELO==time)[0]
    #rsize=arr_length(visdata=visdata,ELO=ELO,curr_time=time)
    #if rsize<nbl:continue

    igtime=np.where(goodtimes==itime)

    nb=0
    for ant1 in np.unique(visdata['antenna1']):
        for ant2 in np.unique(visdata['antenna2']):
            if ant1 < ant2:
                thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
                iant1=np.where(antlist==ant1)[0]
                iant2=np.where(antlist==ant2)[0]
                if thisbase.sum()>0:
                    pt=visdata['data'][0][thisbase][0][thistime][0]


for time in ELO:
    thistime=(visdata['axis_info']['time_axis']['MJDseconds']==time)

    print(visdata['antenna1'])
    #print(visdata['antenna2'])

    allants=np.concatenate((visdata['antenna1'][thistime],visdata['antenna2'][thistime]))
    antlist=np.unique(allants)


    #Calculating number of antennas and baselines
    nant=int(len(antlist))
    nbl=int(nant*(nant-1)/2)
    #print(nbl)



tsize,goodtimes=t_length(visdata=visdata,ELO=ELO,nbl=nbl)

#print(goodtimes)

allts=visdata['axis_info']['time_axis']['MJDseconds']

ELO_diff=np.max(ELO)-np.min(ELO)
print(ELO_diff)
nint=int(ELO_diff/600)
ELO_range=np.linspace(np.min(ELO),np.max(ELO),nint)

#Jacobian and measured phase matrix
Theta_r=np.zeros((nbl,nant-1,tsize),dtype=int)
theta_m=np.zeros((nbl,1,tsize),dtype=float)

#Defining Jacobian and measured amp matrix
script_L=np.zeros((nbl,nant,tsize),dtype=int)
l_m=np.zeros((nbl,1,tsize),dtype=float)
anttrack=np.full(shape=(nbl,2,tsize),fill_value=-1,dtype=int)

print(visdata['axis_info']['time_axis']['MJDseconds'][0])

#Setting integer for cycling through baselines
nb=0
tb=0
for time in ELO:
    thistime=(visdata['axis_info']['time_axis']['MJDseconds']==time)
    vtime=visdata['axis_info']['time_axis']['MJDseconds']
    itime=np.where(ELO==time)[0]
    #rsize=arr_length(visdata=visdata,ELO=ELO,curr_time=time)
    #if rsize<nbl:continue

    igtime=np.where(goodtimes==itime)

    nb=0
    for ant1 in np.unique(visdata['antenna1']):
        for ant2 in np.unique(visdata['antenna2']):
            if ant1 < ant2:
                thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
                iant1=np.where(antlist==ant1)[0]
                iant2=np.where(antlist==ant2)[0]
                if thisbase.sum()>0:
                    pt=visdata['data'][0][thisbase][0][thistime][0]
                    ph=np.angle(pt,deg=True)
                    amp=np.absolute(pt)
                    
                    if amp<=0: continue
                    
                    l_m[nb,0,igtime]=np.log10(amp)
                    anttrack[nb,0,igtime]=ant1
                    anttrack[nb,1,igtime]=ant2

                    script_L[nb,iant1,igtime]=1
                    script_L[nb,iant2,igtime]=1

                    #PHASE STUFF
                    theta_m[nb]=ph

                    if ant1==refant: Theta_r[nb,iant2-1,igtime]=-1
                    if ant2==refant: Theta_r[nb,iant1,igtime]=1
                    if ant1!=refant and ant1>refant:
                        Theta_r[nb,iant1-1,igtime]=1
                        Theta_r[nb,iant2-1,igtime]=-1
                    if ant1!=refant and ant2<refant:
                        Theta_r[nb,iant1,igtime]=1
                        Theta_r[nb,iant2,igtime]=-1
                    if (ant1!=refant and (ant2>refant and ant1<refant)):
                        Theta_r[nb,iant1,igtime]=1
                        Theta_r[nb,iant2-1,igtime]=-1


                    nb+=1
    tb+=1

print("script_L \n"+str(script_L))
print("l_m \n"+str(l_m))

l_r=np.dstack([l_Ir(ll_r=script_L[:,:,x],ll_m=l_m[:,:,x]) for x in range(tsize)])

print("script_L \n"+str(script_L))
print("l_m \n"+str(l_m))

#print(script_L[:,:,0].shape)
#print(l_m[:,:,0].shape)

#print(range(goodtimes))

l_r=np.dstack([l_Ir(ll_r=script_L[:,:,x],ll_m=l_m[:,:,x]) for x in range(tsize)])
#t1=np.linalg.inv(np.matmul(Theta_r.T,Theta_r))
#t2=np.matmul(t1,Theta_r.T)
#theta_Ir=np.matmul(t2,theta_m)

print("l_r \n"+str(l_r))
#print(l_r[:,:,520])


#Residuals
l_del=np.dstack([l_m[:,:,x]-np.matmul(script_L[:,:,x],l_r[:,:,x]) for x in range(tsize)])

print("l_del \n"+str(l_del))
print(l_del.shape)

l_Ir_res=np.dstack([l_Ir(ll_r=script_L[:,:,x],ll_m=l_del[:,:,x]) for x in range(tsize)])

print("l_Ir w/ resids: \n"+str(l_Ir_res))

l_Ir_final=l_r+l_Ir_res
print(l_Ir_final.shape)


print("final amps (log form) \n"+str(l_Ir_final))

Ir_converted=np.zeros_like(l_Ir_final)

#for t in range(tsize):
#    Ir_converted[:,:,t]=10.0**(l_Ir_final[:,:,t])

Ir_converted=np.dstack([10.0**(l_Ir_final[:,:,x]) for x in range(tsize)])



print("Theta_r \n"+str(Theta_r))
print("theta_m \n"+str(theta_m))


#l_r=np.dstack([l_Ir(ll_r=script_L[:,:,x],ll_m=l_m[:,:,x]) for x in range(tsize)])

theta_Ir=np.dstack([th_Ir(Th_r=Theta_r[:,:,x],th_m=theta_m[:,:,x]) for x in range(tsize)])
#t1=np.linalg.inv(np.matmul(Theta_r.T,Theta_r))
#t2=np.matmul(t1,Theta_r.T)
#theta_Ir=np.matmul(t2,theta_m)

print("theta_Ir \n"+str(theta_Ir))

#Residuals
theta_del=np.dstack([theta_m[:,:,x]-np.matmul(Theta_r[:,:,x],theta_Ir[:,:,x]) for x in range(tsize)])
#theta_del=np.matmul(Theta_r,theta_Ir)-theta_m


print("theta_del \n"+str(theta_del))

theta_Ir_res=np.dstack([th_Ir(Th_r=Theta_r[:,:,x],th_m=theta_del[:,:,x]) for x in range(tsize)])

print("theta_Ir w/ resids: \n"+str(theta_Ir_res))

theta_Ir_final=theta_Ir+theta_Ir_res

print("final amps (converted) \n"+str(Ir_converted))
print(Ir_converted.shape)

print("final phases \n"+str(theta_Ir_final))
print(theta_Ir_final.shape)



#datachunk=[x for x in visdata['data'] if (visdata['axis_info']['time_axis']['MJDseconds'][y]>=ELO[t_step] and visdata['axis_info']['time_axis']['MJDseconds'][y]<=ELO[t_step+1] for y in range(len(allts)))][0]
#print(len(datachunk))
all_amp=np.zeros((nbl*len(goodtimes),2),dtype=float)
all_ph=np.zeros((nbl*len(goodtimes),2),dtype=float)
j=0
for time in ELO:
    thistime=(visdata['axis_info']['time_axis']['MJDseconds']==time)
    vtime=visdata['axis_info']['time_axis']['MJDseconds']
    itime=np.where(ELO==time)[0]
    rsize=arr_length(visdata=visdata,ELO=ELO,curr_time=time)
    if rsize<nbl:continue

    #script_L=np.zeros((rsize,nant),dtype=int)
    #l_m=np.zeros((rsize,1),dtype=float)

    #script_L=np.zeros((nbl,nant),dtype=int)
    #l_m=np.zeros((nbl,1),dtype=float)

    igtime=np.where(goodtimes==itime)

    nb=0

    for ant1 in np.unique(visdata['antenna1']):
            for ant2 in np.unique(visdata['antenna2']):
                if ant1 < ant2:
                    thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
                    iant1=np.where(antlist==ant1)[0]
                    iant2=np.where(antlist==ant2)[0]

                    ph=np.angle(visdata['data'][0][thisbase][0][thistime][0],deg=True)
                    pt=visdata['data'][0][thisbase][0][thistime][0]
                    amp=np.absolute(pt)
                    
                    if amp<=0: continue
                    
                    all_ph[j,0]=ph
                    all_ph[j,1]=time
                    all_amp[j,0]=amp
                    all_amp[j,1]=time

                    #print(iant1)
                    #print(iant2)
                    #if thisbase.sum()>0:
                        #print("ph!")
                        #print(datachunk[0][thisbase])
                        #print(datachunk[0][thisbase][0])
                    j+=1

print(np.count_nonzero(all_ph))
print(all_ph.shape)
#print(all_ph[:,0])
#print(all_ph[:,1])

print(theta_Ir_final[:,:,0])
#print(len(goodtimes))
#theta_Ir_final_new=np.expand_dims(theta_Ir_final,axis=3)
#print(theta_Ir_final_new.shape)
#for i in goodtimes:
#    theta_Ir_final_new1=np.append(theta_Ir_final_new,ELO[i],axis=3)
#print(theta_Ir_final_new1.shape)
#print(theta_Ir_final_new1[:,:,:,0])

#ampdict=dict.fromkeys(["amp","time"])
#phdict=dict.fromkeys(["phase","time"])
goodELO=[ELO[x] for x in goodtimes]
phsplit=[theta_Ir_final[:,:,x] for x in range(len(theta_Ir_final[0,0,:]))]
ampsplit=[Ir_converted[:,:,x] for x in range(len(Ir_converted[0,0,:]))]
ampdict={'amp':ampsplit,'time':goodELO}
phdict={'phase':phsplit,'time':goodELO}

#for i in range(len(goodtimes)):
    #ampdict.update({ELO[goodtimes[i]]:Ir_converted[:,:,i]})
    #phdict["phase"]=[phdict["phase"],theta_Ir_final[:,:,i]]
    #phdict["time"]=[phdict["time"],ELO[goodtimes[i]]]
    #ampdict["amp"]=[ampdict["amp"],Ir_converted[:,:,i]]
    #ampdict["time"]=[ampdict["time"],ELO[goodtimes[i]]]
    #phdict.update({"phase":phdict["phase"],theta_Ir_final[:,:,i]},"time":phdict["time"],ELO[goodtimes[i]]})
    #ampdict.update({"amp":Ir_converted[:,:,i],"time":ELO[goodtimes[i]]})
    #phdict.update({ELO[goodtimes[i]]:theta_Ir_final[:,:,i]})
    #ampdict.update({ELO[goodtimes[i]]:Ir_converted[:,:,i]})

print("eeee")
print(phdict['phase'][1])
print("peeeee")
print(phdict['phase'][0])
print(phdict['phase'][0][0])


ampkeys=ampdict.keys()
phkeys=phdict.keys()

nonan=0

for t_step in range(len(ELO_range)-1):
    datachunk_ph=np.mean([phdict['phase'][x] for x in range(len(phdict['phase'])) if (phdict['time'][x]>=ELO_range[t_step] and phdict['time'][x]<=ELO_range[t_step+1])])
    datachunk_amp=np.mean([ampdict['amp'][x] for x in range(len(ampdict['amp'])) if (ampdict['time'][x]>=ELO_range[t_step] and ampdict['time'][x]<=ELO_range[t_step+1])])
    #datachunk_ph=[x['phase'] for x in phdict if x["time"]>=ELO_range[t_step] and x["time"]<=ELO_range[t_step+1]]
    #datachunk_amp=[x['amp'] for x in ampdict if x["time"]>=ELO_range[t_step] and x["time"]<=ELO_range[t_step+1]]
    #datachunk_ph=[phdict[x] for x in range(len(phdict)-1) if (float(phkeys[x])>=ELO_range[t_step] and float(phkeys[x])<=ELO_range[t_step+1])]
    #datachunk_amp=[ampdict[x] for x in range(len(ampdict)-1) if (float(ampkeys[x])>=ELO_range[t_step] and float(ampkeys[x])<=ELO_range[t_step+1])]
    #datachunk_ph=[[:,0,x] for x in range(len(theta_Ir_final[0,0,:])) if (theta_Ir_final[:,:,y]>=ELO_range[t_step] and theta_Ir_final[:,:,y]<=ELO_range[t_step+1] for y in range(len(Ir_converted[0,0,:])))]
    #datachunk_amp=[Ir_converted[:,0,x] for x in range(len(Ir_converted[0,0,:])) if (Ir_converted[:,:,y]>=ELO_range[t_step] and Ir_converted[:,:,y]<=ELO_range[t_step+1] for y in range(len(Ir_converted[0,0,:])))]
    #datachunk_amp=[Ir_converted[x,0] for x in range(len(all_amp[:,0])) if (all_amp[x,1]>=ELO_range[t_step] and all_amp[x,1]<=ELO_range[t_step+1])]
    if np.isnan(datachunk_ph)==True or np.isnan(datachunk_amp)==True: continue
    print(ELO_range[t_step])
    print(ELO_range[t_step+1])

    print(datachunk_ph)
    print(datachunk_amp)
    nonan+=1

print(nonan)
print(len(ELO_range))
'''
