import numpy as np
import matplotlib.pyplot as plt

def l_Ir(ll_r,ll_m):
    t1=np.linalg.inv(np.matmul(ll_r.T,ll_r))
    t2=np.matmul(t1,ll_r.T)
    ll_Ir=np.matmul(t2,ll_m)
    return ll_Ir

#Data file visibilities will be pulled from
datams='sgr_apr07_flagcor.ms'
ms.open(datams,nomodify=True)


#Collect data from ms
visdata = ms.getdata(['antenna1','antenna2','data','data_desc_id','axis_info'],ifraxis=True)
visdata['data'] = np.squeeze(visdata['data'])
ms.close()

print(len(visdata['axis_info']['time_axis']['MJDseconds']))
print(len(np.unique(visdata['axis_info']['time_axis']['MJDseconds'])))

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

#Defining Jacobian and measured phase matrix
script_L=np.zeros((nbl,nant-1),dtype=int)
r_size=script_L.shape
l_m=np.zeros((nbl,1),dtype=float)

'''
for ant1 in np.unique(visdata['antenna1']):
    for ant2 in np.unique(visdata['antenna2']):
        if ant1 < ant2:
            thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
            #print(len(visdata['data'][0][thisbase]))

#Setting integer for cycling through baselines
nb=0
refant=4
for ant1 in np.unique(visdata['antenna1']):
    for ant2 in np.unique(visdata['antenna2']):
        if ant1 < ant2:
            thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
            iant1=np.where(antlist==ant1)[0]
            iant2=np.where(antlist==ant2)[0]
            #print(iant1)
            #print(iant2)
            if thisbase.sum()>0:
                nb=0
                for pt in np.arange(len(visdata['data'][0][thisbase])):
                    #print(len(visdata['data'][0][thisbase]))
                    amp=np.absolute(visdata['data'][0][thisbase][pt])
                    l_m[nb]=np.log10(amp)
                    '''
'''
                    if ant2==4: Theta_r[nb,ant1]=1
                    if ant2!=4:
                        Theta_r[nb,ant1]=1
                        Theta_r[nb,ant2]=-1
                    '''
'''
                    if ant1==refant: script_L[nb,iant2-1]=1
                    if ant2==refant: script_L[nb,iant1]=1
                    if ant1!=refant and ant1>refant:
                        script_L[nb,iant1-1]=1
                        script_L[nb,iant2-1]=1
                    if ant1!=refant and ant2<refant:
                        script_L[nb,iant1]=1
                        script_L[nb,iant2]=1
                    if (ant1!=refant and (ant2>refant and ant1<refant)):
                        script_L[nb,iant1]=1
                        script_L[nb,iant2-1]=1



                    #if ant1==0: Theta_r[nb,ant2-1]=-1
                    #if ant1!=0:
                    #    Theta_r[nb,ant1-1]=1
                    #    Theta_r[nb,ant2-1]=-1

                    nb+=1

print("script_L \n"+str(script_L))
print("l_m \n"+str(l_m))

l_r=l_Ir(ll_r=script_L,ll_m=l_m)
#t1=np.linalg.inv(np.matmul(Theta_r.T,Theta_r))
#t2=np.matmul(t1,Theta_r.T)
#theta_Ir=np.matmul(t2,theta_m)

print("l_r \n"+str(l_r))

#Residuals
l_del=l_m-np.matmul(script_L,l_r)

print("l_del \n"+str(l_del))

l_Ir_res=l_Ir(ll_r=script_L,ll_m=l_del)

print("l_Ir w/ resids: \n"+str(l_Ir_res))

l_Ir_final=l_r+l_Ir_res

print("final amps (log form) \n"+str(l_Ir_final))

Ir_converted=np.vstack([10.0**x for x in l_Ir_final])

print("final amps (converted) \n"+str(Ir_converted))

bpts=range(nant-1)
#print(len(bpts))
#print([:,0].shape)


plt.scatter(bpts,Ir_converted[:,0])
plt.show()

'''
