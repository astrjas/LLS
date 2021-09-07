import numpy as np
import matplotlib.pyplot as plt

#Data file visibilities will be pulled from
datams='sgr_apr07_flagcor.ms'

#Setting integer for cycling through baselines
refant=24



def th_Ir(Th_r,th_m):
    t1=np.linalg.inv(np.matmul(Th_r.T,Th_r))
    t2=np.matmul(t1,Th_r.T)
    theta_Ir=np.matmul(t2,th_m)
    return theta_Ir

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

ms.open(datams,nomodify=True)

#Collect data from ms
visdata = ms.getdata(['antenna1','antenna2','data','data_desc_id','axis_info'],ifraxis=True)
visdata['data'] = np.squeeze(visdata['data'])
ms.close()

allants=np.concatenate((visdata['antenna1'],visdata['antenna2']))
antlist=np.unique(allants)
#print(antlist)

#Calculating number of antennas and baselines
nant=int(len(antlist))
nbl=int(nant*(nant-1)/2)

ELO=np.unique(visdata['axis_info']['time_axis']['MJDseconds'])


print(nant)
print(nbl)

tsize,goodtimes=t_length(visdata=visdata,ELO=ELO,nbl=nbl)
print(goodtimes)


#Defining Jacobian and measured phase matrix
Theta_r=np.zeros((nbl,nant-1,tsize),dtype=int)
r_size=Theta_r.shape
theta_m=np.zeros((nbl,1,tsize),dtype=float)
anttrack=np.full(shape=(nbl,2,tsize),fill_value=-1,dtype=int)


nb=0
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
                    print("ph!")
                    print(visdata['data'][0][thisbase][0][thistime][0])
                    ph=np.angle(visdata['data'][0][thisbase][0][thistime][0],deg=True)

                    anttrack[nb,0,igtime]=ant1
                    anttrack[nb,1,igtime]=ant2


                    #if iant1==1: ph=ph+90
                    #if iant2==1: ph-=90
                    theta_m[nb]=ph
                    '''
                    if ant2==4: Theta_r[nb,ant1]=1
                    if ant2!=4:
                        Theta_r[nb,ant1]=1
                        Theta_r[nb,ant2]=-1
                    '''
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



                    #if ant1==0: Theta_r[nb,ant2-1]=-1
                    #if ant1!=0:
                    #    Theta_r[nb,ant1-1]=1
                    #    Theta_r[nb,ant2-1]=-1

                    nb+=1
    tb+=1

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

print("final phases \n"+str(theta_Ir_final))

bpts=range(nant-1)
print(len(bpts))
print(theta_Ir_final[:,0].shape)

theta_rms = np.dstack([np.sqrt(np.matmul(theta_Ir_res[:,:,x].T,theta_Ir_res[:,:,x])) for x in range(tsize)])

print(np.mean(theta_rms))

for t in range(tsize):
    bits=np.zeros((nant-1))
    bits[:]=ELO[goodtimes[t]]
    plt.scatter(bits,theta_Ir_final[:,0,t])
plt.show()

