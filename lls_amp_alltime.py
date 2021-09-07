import numpy as np
import matplotlib.pyplot as plt

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


#Data file visibilities will be pulled from
datams='sgr_apr07_flagcor_tenbaseline.ms'
ms.open(datams,nomodify=True)


#Collect data from ms
visdata = ms.getdata(['antenna1','antenna2','data','data_desc_id','axis_info'],ifraxis=True)
visdata['data'] = np.squeeze(visdata['data'])
ms.close()

print(len(visdata['axis_info']['time_axis']['MJDseconds']))
print(len(np.unique(visdata['axis_info']['time_axis']['MJDseconds'])))

ELO=np.unique(visdata['axis_info']['time_axis']['MJDseconds'])
#print(ELO)

allants=np.concatenate((visdata['antenna1'],visdata['antenna2']))
antlist=np.unique(allants)
#print(antlist)

#Calculating number of antennas and baselines
nant=int(len(antlist))
nbl=int(nant*(nant-1)/2)


print(nant)
print(nbl)
#print(len(np.unique(visdata['antenna1'])))
#print(len(np.unique(allants)))

tsize,goodtimes=t_length(visdata=visdata,ELO=ELO,nbl=nbl)
print("TSIZE")
print(goodtimes)

#Defining Jacobian and measured amp matrix
script_L=np.zeros((nbl,nant,tsize),dtype=int)
#script_L=np.array([],dtype=float)
#r_size=script_L.shape

l_m=np.zeros((nbl,1,tsize),dtype=float)
#l_m=np.array([],dtype=float)

anttrack=np.full(shape=(nbl,2,tsize),fill_value=-1,dtype=int)
#anttrack=np.array([],dtype=int)
#anttrack={}

#l_r=np.zeros((nbl,1,len(ELO)),dtype=float)
#l_r=np.array([],dtype=float)

'''
for ant1 in np.unique(visdata['antenna1']):
    for ant2 in np.unique(visdata['antenna2']):
        if ant1 < ant2:
            thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
            #print(len(visdata['data'][0][thisbase]))
'''


#Setting integer for cycling through baselines
nb=0
#refant=4
tb=0
for time in ELO:
    thistime=(visdata['axis_info']['time_axis']['MJDseconds']==time)
    vtime=visdata['axis_info']['time_axis']['MJDseconds']
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

    nb=0
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
                    #print(len(pt))
                    amp=np.absolute(pt)
                    #print(amp)
                    
                    if amp<=0: continue
                        #print("ah!")
                        #nb+=1
                        #continue
                    
                    #else:
                    l_m[nb,0,igtime]=np.log10(amp)
                    anttrack[nb,0,igtime]=ant1
                    anttrack[nb,1,igtime]=ant2

                    script_L[nb,iant1,igtime]=1
                    script_L[nb,iant2,igtime]=1

                    nb+=1
    #script_L_chunk=script_L[:,:,tb]
    #l_m_chunk=l_m[:,:,tb]
    #if rsize==1: l_r=
    #l_r=l_Ir(ll_r=script_L,ll_m=l_m)
    #print(str(itime))
    tb+=1

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

print("final amps (converted) \n"+str(Ir_converted))
print(Ir_converted.shape)

bpts=range(nant)
#print(len(bpts))
#print([:,0].shape)

#print(ELO[goodtimes[0]])
#print(Ir_converted[:,0,0])

l_rms = np.dstack([np.sqrt(np.matmul(l_Ir_res[:,:,x].T,l_Ir_res[:,:,x])) for x in range(tsize)])

print(np.mean(l_rms))


for t in range(tsize):
    bits=np.zeros((nant))
    bits[:]=ELO[goodtimes[t]]
    plt.scatter(bits,Ir_converted[:,0,t])
plt.show()


