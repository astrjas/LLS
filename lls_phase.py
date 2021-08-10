import numpy as np

#Speed of light
c=3.0e8

#Data file visibilities will be pulled from
datams='sgr_apr07_flagcor_tenbaseline.ms'
#print(datams[:-3])
ms.open(datams,nomodify=True)

#Collect data from ms
visdata = ms.getdata(['uvw','antenna1','antenna2','data','sigma','data_desc_id'])
print(visdata['data'].shape)
visdata['data'] = np.squeeze(visdata['data'])
print(visdata['data'].shape)
ms.close()
tb.open(datams+'/SPECTRAL_WINDOW')
freqs = np.squeeze(tb.getcol('CHAN_FREQ'))
print(freqs)
print("Freq shape: "+str(freqs.shape))
tb.close()

# Get the primary beam size
tb.open(datams+'/ANTENNA')
diam = np.squeeze(tb.getcol('DISH_DIAMETER'))[0] # assumes all dishes are equal size
tb.close()
PBfwhm = 1.2*(c/np.mean(freqs))/diam * (3600*180/np.pi) # in arcsec
print("Shape [:,:]: "+str(visdata['data'].shape))
print("Shape vd[sigma]: "+str(visdata['sigma'].shape))
# Current data+sigma has two pols, average them.
np.average(visdata['data'],weights=(visdata['sigma']**-2.),axis=0)
print("bingo, ho hoh ho hoh")
#print(len(newlist))
print("New visdata shape: "+str(visdata['data'].shape))
#visdata['data'] = np.average(visdata['data'],weights=(visdata['sigma']**-2.),axis=0)
visdata['sigma'] = np.sum((visdata['sigma']**-2.),axis=0)**-0.5 # these are uncertainties, average them down

# Convert uvw coords from m to wavelengths
freqyc=freqs/c
newvis=[visdata['uvw'] * freqs[x]/c for x in range(len(freqs))]
newvis1=np.asarray(newvis)
print("uvw coords array:",newvis1.shape)

nant=int(len(np.unique(visdata['antenna2'])))+1
nbl=int(nant*(nant-1)/2)

Theta_r=np.zeros((nbl,nant-1),dtype=int)
r_size=Theta_r.shape
theta_m=np.zeros((nbl,1),dtype=float)

nb=0
#print("A1")
#print(len(visdata['data'][:,0]))
print("A2")
print(np.unique(visdata['antenna2']).shape)
for ant1 in np.unique(visdata['antenna1']):
    for ant2 in np.unique(visdata['antenna2']):
        if ant1 < ant2:
            thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
            if thisbase.sum()>0:
                ph=np.angle(visdata['data'][0][thisbase][0],deg=True)
                theta_m[nb]=ph
                if ant1==0:
                    Theta_r[nb,ant2-1]=-1
                if ant1!=0:
                    Theta_r[nb,ant1-1]=1
                    Theta_r[nb,ant2-1]=-1
                nb+=1
                '''
                for i in range(r_size[0]):
                    #theta_m[i]=np.angle(visdata['data'][thisbase])
                    for j in range(r_size[1]):
                        if j+1==ant1:
                            Theta_r[i,j]=1
                        if j+1==ant2:
                            Theta_r[i,j]=-1
                    print(Theta_r)
                '''
                '''
                reals1 = [visdata['data'][x,:].real[thisbase] for x in range(len(visdata['data'][:,0]))]
                reals=np.asarray(reals1)
                imags1 = [visdata['data'][x,:].imag[thisbase] for x in range(len(visdata['data'][:,0]))]
                imags=np.asarray(imags1)
                #print(imags.shape)
                sigs = visdata['sigma'][thisbase]
                diffrs = reals - np.roll(reals,-1); diffis = imags - np.roll(imags,-1)
                std = np.mean([diffrs.std(),diffis.std()])
                facs.append(std/sigs.mean()/np.sqrt(2)) # for calculating a single rescaling
                '''
print(Theta_r)
print(theta_m)


phases=np.angle(visdata['data'],deg=True)
#facs = np.asarray(facs)
#print(facs.shape)
#visdata['sigma'] *= facs.mean()
#print datams, facs.mean(), ((visdata['sigma']**-2).sum())**-0.5

print(phases)
print("Phase shapes:"+str(phases.shape))

