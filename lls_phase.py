import numpy as np

def th_Ir(Th_r,th_m):
    t1=np.linalg.inv(np.matmul(Th_r.T,Th_r))
    t2=np.matmul(t1,Th_r.T)
    theta_Ir=np.matmul(t2,th_m)
    return theta_Ir

#Data file visibilities will be pulled from
datams='sgr_apr07_flagcor_tenbaseline.ms'
ms.open(datams,nomodify=True)

#Collect data from ms
visdata = ms.getdata(['antenna1','antenna2','data','data_desc_id'])
visdata['data'] = np.squeeze(visdata['data'])
ms.close()

#Calculating number of antennas and baselines
nant=int(len(np.unique(visdata['antenna1'])))+1
nbl=int(nant*(nant-1)/2)

print(nant)
print(nbl)
print(len(np.unique(visdata['antenna1'])))

#Defining Jacobian and measured phase matrix
Theta_r=np.zeros((nbl,nant-1),dtype=int)
r_size=Theta_r.shape
theta_m=np.zeros((nbl,1),dtype=float)

#Setting integer for cycling through baselines
nb=0
refant=4
for ant1 in np.unique(visdata['antenna1']):
    for ant2 in np.unique(visdata['antenna2']):
        if ant1 < ant2:
            thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
            #iant1=np.unique(visdata['antenna1']).index(ant1)
            #iant2=visdata['antenna1'].index(ant2)
            if thisbase.sum()>0:
                ph=np.angle(visdata['data'][0][thisbase][10],deg=True)
                theta_m[nb]=ph
                '''
                if ant2==4: Theta_r[nb,ant1]=1
                if ant2!=4:
                    Theta_r[nb,ant1]=1
                    Theta_r[nb,ant2]=-1
                '''
                if ant1==refant: Theta_r[nb,ant2-1]=-1
                if ant2==refant: Theta_r[nb,ant1]=1
                if ant1!=refant and ant1>refant:
                    Theta_r[nb,ant1-1]=1
                    Theta_r[nb,ant2-1]=-1
                if ant1!=refant and ant2<refant:
                    Theta_r[nb,ant1]=1
                    Theta_r[nb,ant2]=-1
                if (ant1!=refant and (ant2>refant and ant1<refant)):
                    Theta_r[nb,ant1]=1
                    Theta_r[nb,ant2-1]=-1
                    

                
                #if ant1==0: Theta_r[nb,ant2-1]=-1
                #if ant1!=0:
                #    Theta_r[nb,ant1-1]=1
                #    Theta_r[nb,ant2-1]=-1
                
                nb+=1

print("Theta_r \n"+str(Theta_r))
print("theta_m \n"+str(theta_m))

theta_Ir=th_Ir(Th_r=Theta_r,th_m=theta_m)
#t1=np.linalg.inv(np.matmul(Theta_r.T,Theta_r))
#t2=np.matmul(t1,Theta_r.T)
#theta_Ir=np.matmul(t2,theta_m)

print("theta_Ir \n"+str(theta_Ir))

#Residuals
theta_del=theta_m-np.matmul(Theta_r,theta_Ir)

print("theta_del \n"+str(theta_del))

theta_Ir_res=th_Ir(Th_r=Theta_r,th_m=theta_del)

print("theta_Ir w/ resids: \n"+str(theta_Ir_res))

theta_Ir_final=theta_Ir+theta_Ir_res

print("final phases \n"+str(theta_Ir_final))

