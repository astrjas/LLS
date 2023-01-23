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

attn=0

#try:
#    import uvmultimodel
#except ModuleNotFoundError:
import sys
sys.path.insert(0,"/home/jasminwash/.local/lib/python2.7/site-packages/NordicARC-3.0.5-py2.7-linux-x86_64.egg")
#import uvmultimodel


from NordicARC import uvmultifit as uvm

def antfunc_clean(dfile,it,pol,date):
    itstr=str(it)
    ms.open(dfile,nomodify=True)
    visdata = ms.getdata(['antenna1','antenna2','flag','model_data','data_desc_id','sigma','axis_info'],ifraxis=True)
    visdata['model_data'] = np.squeeze(visdata['model_data'])
    ms.close()

    allants=np.concatenate((visdata['antenna1'],visdata['antenna2']))
    antlist=np.unique(allants)

    npol=visdata['model_data'].shape[0]
    polnames=visdata['axis_info']['corr_axis']
    corr=polnames[pol]
    nant=len(antlist)
    nbl=int(nant*(nant-1)/2)
    ntimes=len(visdata['axis_info']['time_axis']['MJDseconds'])

    tt=np.array(visdata['axis_info']['time_axis']['MJDseconds'],dtype=float)

    #Creating a tf dictionary for each ant and whether time is good or bad there
    antinfo=dict()
    antinfo['timestamp']=visdata['axis_info']['time_axis']['MJDseconds']

    allant=np.zeros((nant,ntimes),dtype=bool)
    alltim=np.zeros((ntimes,nant),dtype=bool)



    #TENTATIVE BEGINNING OF WHAT WILL BE IN MASSIVE LOOP
    #XX correlation
    xx=np.zeros_like(visdata['flag'][pol])
    flagpt=(visdata['flag'][pol]==True)
    xx=np.squeeze(np.where(flagpt==True,np.nan,visdata['model_data'][pol]))
    print("xxshape",xx.shape)

    for y in range(nbl):
        pp1=np.array(xx[y])
        #plt.scatter(tt,pp)
        plt.scatter(tt,np.abs(pp1))
    plt.show()
    #plt.savefig('datatestplot_'+auth+refmeth+'_'+date+corr+str(it)+'.png',overwrite=True)
    plt.savefig('datatestplot_cmodel_'+str(it)+date+'1.png',overwrite=True)
    plt.clf()

    #tf matrix for where ant is bad (True)
    tfdata=np.isnan(np.abs(xx))
    print(tfdata)

    #cycle through each ant
    #a=0
    for a in range(nant):
        ant=antlist[a]
        #all baselines w/ this ant
        thisant=(visdata['antenna1']==ant) | (visdata['antenna2']==ant)
        #pull all baselines for this ant at all times w/in tfdata
        allt1ant=tfdata[thisant][:]
        #evaluating if baseline/antenna is bad at each time
        ant_tf=np.all(allt1ant,axis=0)
        #Adding to goodant
        allant[a,:]=ant_tf
        #print("ant_tf",ant_tf.shape)
        plt.plot(ant_tf)
        #plotting and saving this
        plt.savefig("./tfgraphs/ant_tf_ant%i_%s_it%i1.png"%(ant,corr,it),overwrite=True)
        plt.clf()
        #a+=1

    antinfo['goodant_'+corr+itstr]=allant


    for i in range(ntimes):
        j=0
        tfstore=np.zeros(len(antlist),dtype=bool)
        ct=visdata['axis_info']['time_axis']['MJDseconds'][i]
        #print("tstamp",ct)
        for ant in range(nant):
            tfstore[j]=allant[ant][i]
            j+=1
        alltim[i,:]=tfstore


    antinfo['goodt_'+corr+itstr]=alltim
    return antinfo,xx


def antfunc(dfile,it,pol,date):
    itstr=str(it)
    ms.open(dfile,nomodify=True)
    visdata = ms.getdata(['antenna1','antenna2','flag','data','data_desc_id','sigma','axis_info'],ifraxis=True)
    visdata['data'] = np.squeeze(visdata['data'])
    ms.close()

    allants=np.concatenate((visdata['antenna1'],visdata['antenna2']))
    antlist=np.unique(allants)

    npol=visdata['data'].shape[0]
    polnames=visdata['axis_info']['corr_axis']
    corr=polnames[pol]
    nant=len(antlist)
    nbl=int(nant*(nant-1)/2)
    ntimes=len(visdata['axis_info']['time_axis']['MJDseconds'])

    tt=np.array(visdata['axis_info']['time_axis']['MJDseconds'],dtype=float)

    #Creating a tf dictionary for each ant and whether time is good or bad there
    antinfo=dict()
    antinfo['timestamp']=visdata['axis_info']['time_axis']['MJDseconds']

    allant=np.zeros((nant,ntimes),dtype=bool)
    alltim=np.zeros((ntimes,nant),dtype=bool)

    xx=np.zeros_like(visdata['flag'][pol])
    flagpt=(visdata['flag'][pol]==True)
    xx=np.squeeze(np.where(flagpt==True,np.nan,visdata['data'][pol]))
    print("xxshape",xx.shape)

    for y in range(nbl):
        pp1=np.array(xx[y])
        #plt.scatter(tt,pp)
        plt.scatter(tt,np.abs(pp1))
    plt.show()
    #plt.savefig('datatestplot_'+auth+refmeth+'_'+date+corr+str(it)+'.png',overwrite=True)
    plt.savefig('datatestplot_'+str(dfile[:-3])+'_'+str(it)+date+'1.png',overwrite=True)
    plt.clf()

    #tf matrix for where amp of ant is bad (True)
    tfdata=np.isnan(np.log10(np.abs(xx)))
    #np.log10(np.abs(xx))
    print(tfdata)

    #cycle through each ant
    #a=0
    for a in range(nant):
        ant=antlist[a]
        #all baselines w/ this ant
        thisant=(visdata['antenna1']==ant) | (visdata['antenna2']==ant)
        #pull all baselines for this ant at all times w/in tfdata
        allt1ant=tfdata[thisant][:]
        #evaluating if baseline/antenna is bad at each time
        ant_tf=np.all(allt1ant,axis=0)
        #Adding to goodant
        allant[a,:]=ant_tf
        #print("ant_tf",ant_tf.shape)
        plt.plot(ant_tf)
        #plotting and saving this
        plt.savefig("./tfgraphs/ant_tf_ant%i_%s_it%i1.png"%(ant,corr,it),overwrite=True)
        plt.clf()
        #a+=1

    antinfo['goodant_'+corr+itstr]=allant


    for i in range(ntimes):
        j=0
        tfstore=np.zeros(len(antlist),dtype=bool)
        ct=visdata['axis_info']['time_axis']['MJDseconds'][i]
        #print("tstamp",ct)
        for ant in range(nant):
            tfstore[j]=allant[ant][i]
            j+=1
        alltim[i,:]=tfstore


    antinfo['goodt_'+corr+itstr]=alltim
    return antinfo,xx

####################################################################
def gsolve(nbl,nant,vdata,varr,ainfo,antlist,corr,itstr,refant):
    theta_r_dict={}
    theta_m_dict={}
    theta_Ir_dict={}
    theta_del_dict={}
    theta_rms_dict={}

    script_L_dict={}
    theta_r_nn_dict={}
    l_m_dict={}
    theta_m_nn_dict={}
    l_r_dict={}
    l_del_dict={}
    l_rms_dict={}
    l_Ir_conv_dict={}
    all_bdance={}
    #all_bdance_ph={}

    smc=0

    rfantind=np.where(antlist==refant)[0][0]

    for j in range(len(ainfo['timestamp'])):
        #Jacobian and measured phase matrix
        Theta_r=np.zeros((nbl,nant-1),dtype=int)
        theta_m=np.zeros((nbl,1),dtype=float)

        #Defining Jacobian and measured amp matrix
        script_L=np.zeros((nbl,nant),dtype=int)
        l_m=np.zeros((nbl,1),dtype=float)
        astore=np.zeros((nbl,1),dtype=float)
        bantstore=[]


        #print(time)
        thistime=(vdata['axis_info']['time_axis']['MJDseconds']==ainfo['timestamp'][j])
        time=str(j)

        #Calculating number of antennas and baselines
        #-1 is to account for bad ants (ant1=-1 is not valid)    
        nb=0
        for ant1 in np.unique(vdata['antenna1']):
            for ant2 in np.unique(vdata['antenna2']):
                if ant1 < ant2:
                    thisbase = (vdata['antenna1']==ant1) & (vdata['antenna2']==ant2)
                    iant1=np.where(antlist==ant1)[0][0]
                    #print(iant1)
                    iant2=np.where(antlist==ant2)[0][0]
                    #print(iant2)
                    if thisbase.sum()>0:
                        #potential 0 after thisbase and thistime
                        pt=varr[thisbase][0][j]
                        #print("PT",pt)
                        ph=np.angle(pt,deg=True)
                        amp=np.absolute(pt)
                        if np.isnan(amp)==True: 
                            bantstore.append(iant1)
                            #bantstore.append(iant2)
                        #print("AMP",amp)

                        #if iant1==5 or iant2==5:
                            #amp*=4
                            #ph+=90
                            #print("bing")

                        l_m[nb,0]=np.log10(amp)
                        astore[nb,0]=amp

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
        badance=np.full((nant,1),False,dtype=bool)
        for ia in range(nant):
            #print(bantstore.count(ia))
            #if ia==38: 
            #    if bantstore.count(ia)>0: badance[ia,0]=True
            #else:
            if bantstore.count(ia)==38-ia: badance[ia,0]=True
        #badance=np.unique(bantstore)
        #print(badance)
        plt.scatter(j,len(np.where(badance==True)[0]))
        gf.dict_update(all_bdance,time+corr+itstr,badance)

        gf.dict_update(theta_r_dict,time+corr+itstr,Theta_r)
        gf.dict_update(theta_m_dict,time+corr+itstr,theta_m)


        gf.dict_update(script_L_dict,time+corr+itstr,script_L)
        gf.dict_update(l_m_dict,time+corr+itstr,l_m)
        

        #scl.write("\n\n")

        #lme.write("\n\n")

        #ast.write(str(astore))
        #ast.write("\n\n")


        script_L_nn=gf.nonan(script_L)
        #gf.dict_update(script_L_nn_dict,time+corr+itstr,script_L_nn)
        #scl.write(str(script_L_nn))
        l_m_nn=gf.nonan(l_m)
        #gf.dict_update(l_m_nn_dict,time+corr+itstr,l_m_nn)
        #lme.write(str(l_m_nn))

        #l_r=gf.l_Ir(ll_r=script_L_nn,ll_m=l_m_nn)
        print(len(np.where(l_m_nn.mask==False)[0]))
        print(j)
        l_r=gf.l_Ir_mask(ll_r=script_L_nn,ll_m=l_m_nn,ainfo=ainfo,t=j,nant=nant,corr=corr,itstr=itstr,bdance=badance)


        Theta_r_nn=gf.nonan(Theta_r)
        gf.dict_update(theta_r_nn_dict,time+corr+itstr,Theta_r_nn)
        theta_m_nn=gf.nonan(theta_m)
        gf.dict_update(theta_m_nn_dict,time+corr+itstr,theta_m_nn)

        theta_Ir=gf.th_Ir_mask(Th_r=Theta_r_nn,th_m=theta_m_nn,ainfo=ainfo,t=j,nant=nant,corr=corr,itstr=itstr,rfantind=rfantind,bdance=badance)
        #theta_Ir=gf.th_Ir(Th_r=Theta_r_nn,th_m=theta_m_nn)



        #Residuals
        l_del=l_m_nn-np.matmul(script_L_nn,l_r)
        #print("l_del \n"+str(l_del))
        gf.dict_update(l_del_dict,time+corr+itstr,l_del)
        #print(l_del.shape)

        l_Ir_res=gf.l_Ir_mask(ll_r=script_L_nn,ll_m=l_del,ainfo=ainfo,t=j,nant=nant,corr=corr,itstr=itstr,bdance=badance)

        #print("l_Ir w/ resids: \n"+str(l_Ir_res))

        theta_del=theta_m_nn-np.matmul(Theta_r_nn,theta_Ir)
        #print("theta_del \n"+str(theta_del))
        gf.dict_update(theta_del_dict,time+corr+itstr,theta_del)

        theta_Ir_res=gf.th_Ir_mask(Th_r=Theta_r_nn,th_m=theta_del,ainfo=ainfo,t=j,nant=nant,corr=corr,itstr=itstr,rfantind=rfantind,bdance=badance)
        #theta_Ir_res=gf.th_Ir(Th_r=Theta_r_nn,th_m=theta_del)
        #theta_Ir_res,restd,rnktd,stdel=lstsq(a=Theta_r_nn,b=theta_del)
        #print("theta_Ir w/ resids: \n"+str(theta_Ir_res))


        l_Ir_final=l_r+l_Ir_res
        #ast.write(str(l_Ir_final))
        #l_Ir_final=l_r
        #print("l_ir_shape")
        #print(l_Ir_final.shape)
        gf.dict_update(l_r_dict,time+corr+itstr,l_Ir_final)
        #print("final amps (log form) \n"+str(l_Ir_final))

        Ir_converted=np.zeros_like(l_Ir_final)

        #Ir_converted=np.vstack([10.0**x for x in l_Ir_final])
        for x in range(l_Ir_final.shape[0]):
            lir=l_Ir_final[x,0]
            Ir_converted[x,0]=10.0**lir
            
        #print("l_conv shape")
        #print(Ir_converted.shape)

        gf.dict_update(l_Ir_conv_dict,time+corr+itstr,Ir_converted)

        theta_Ir_final=theta_Ir+theta_Ir_res
        gf.dict_update(theta_Ir_dict,time+corr+itstr,theta_Ir_final)
    #plt.show()
    plt.savefig('validants.png',overwrite=True)
    plt.show()
    plt.clf()
    return theta_Ir_dict, l_Ir_conv_dict


pol=0
corr='XX'
#datams1='sgr_apr07_flag_ext3_it00.ms'
datams1='sgr_apr07_flag_ext3_it00avg240s_it00.ms'
#datams_cm='sgr_apr07_flag_ext3_it00_cmsplit.ms'
datams_cm='sgr_apr07_flag_cmodel_avg240s_it00.ms'
target='sgr_apr07'
case='ref8'
refmeth='manual'
date='01232022'
it=0
dvis='sgr_apr07_flag.ms'
rfant=8
auth='Venki'

####################################################################

#def lls(pol,corr,datams1,target,case,auth,refmeth,date,it,dvis,rfant):
itstr=str(it)

print(datams1)

#Opening data and pulling necessary info
ms.open(datams1,nomodify=True)
ms.selectinit(reset=True)
visdata = ms.getdata(['antenna1','antenna2','data','axis_info','flag'],ifraxis=True)
#printing correlations
print(visdata['axis_info']['corr_axis'])

#Squeeze data then close ms

visdata['data'] = np.squeeze(visdata['data'])
print("data shape",visdata['data'].shape)
ms.close()

allants=np.concatenate((visdata['antenna1'],visdata['antenna2']))
antlist=np.unique(allants)

npol=visdata['data'].shape[0]
polnames=visdata['axis_info']['corr_axis']
corr=polnames[pol]
nant=len(antlist)
nbl=int(nant*(nant-1)/2)
ntimes=len(visdata['axis_info']['time_axis']['MJDseconds'])

#Pulling all unique timestamps, print its length and then all times
ELO=np.unique(visdata['axis_info']['time_axis']['MJDseconds'])
print(len(ELO))
print(len(visdata['axis_info']['time_axis']['MJDseconds']))
#print([x-ELO[0] for x in ELO])

antinfo,xx=antfunc(dfile=datams1,it=it,pol=pol,date=date)
antinfo_clean,xx_clean=antfunc(dfile=datams_cm,it=it,pol=pol,date=date)
print("XX",xx.shape)


allgoodtime=0
allbad=0
good1=0

nb=0

nvis=[]

#Cycling through all the times
for time in antinfo['goodt_'+corr+itstr]:
    if np.any(time)==False: allgoodtime+=1
    if np.any(time)==True:
        nbad=np.count_nonzero(time)
        if nbad==nant: allbad+=1
        if nbad<nant: good1+=1

print("All good antennas:",allgoodtime)
print("Almost all good antennas:",good1)
print("All ants bad:",allbad)

#refant=gf.refantfinder(antlist=antlist,goodant=antinfo['goodant_'+corr+itstr])
#iref=np.where(antlist==refant)[0]
refant=rfant
iref=np.where(antlist==refant)[0]

#scl=open("scriptL_"+target+date+corr+itstr+".txt","w")
#lme=open("lm_"+target+date+corr+itstr+".txt","w")
#ast=open("astore_"+target+date+corr+itstr+".txt","w")

tir_ext,lir_ext,l_ext,ln_ext,sl_ext,sln_ext=gsolve(nbl=nbl,nant=nant,vdata=visdata,varr=xx,ainfo=antinfo,antlist=antlist,corr=corr,itstr=itstr,refant=refant)
tir_cln,lir_cln,l_cln,ln_cln,sl_cln,sln_cln=gsolve(nbl=nbl,nant=nant,vdata=visdata,varr=xx_clean,ainfo=antinfo_clean,antlist=antlist,corr=corr,itstr=itstr,refant=refant)



'''
#def gsolve(nbl,nant,vdict,varr,
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
    astore=np.zeros((nbl,1),dtype=float)


    #print(time)
    thistime=(visdata['axis_info']['time_axis']['MJDseconds']==antinfo['timestamp'][i])
    time=str(i)

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
                    #print("PT",pt)
                    ph=np.angle(pt,deg=True)
                    amp=np.absolute(pt)
                    #print("AMP",amp)

                    #if iant1==5 or iant2==5:
                        #amp*=4
                        #ph+=90
                        #print("bing")

                    l_m[nb,0]=np.log10(amp)
                    astore[nb,0]=amp

                    script_L[nb,iant1]=1
                    script_L[nb,iant2]=1


                    #PHASE STUFF
                    theta_m[nb,0]=ph

                    #print(iant1)
                    #print(iant2)

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

    gf.dict_update(theta_r_dict,time+corr+itstr,Theta_r)
    gf.dict_update(theta_m_dict,time+corr+itstr,theta_m)


    gf.dict_update(script_L_dict,time+corr+itstr,script_L)
    gf.dict_update(l_m_dict,time+corr+itstr,l_m)


    scl.write("\n\n")

    lme.write("\n\n")

    #ast.write(str(astore))
    ast.write("\n\n")


    script_L_nn=gf.nonan(script_L)
    scl.write(str(script_L_nn))
    l_m_nn=gf.nonan(l_m)
    lme.write(str(l_m_nn))

    l_r=gf.l_Ir(ll_r=script_L_nn,ll_m=l_m_nn)

    Theta_r_nn=gf.nonan(Theta_r)
    theta_m_nn=gf.nonan(theta_m)

    theta_Ir=gf.th_Ir(Th_r=Theta_r_nn,th_m=theta_m_nn)



    #Residuals
    l_del=l_m_nn-np.matmul(script_L_nn,l_r)
    #print("l_del \n"+str(l_del))
    gf.dict_update(l_del_dict,time+corr+itstr,l_del)
    #print(l_del.shape)

    l_Ir_res=gf.l_Ir(ll_r=script_L_nn,ll_m=l_del)

    #print("l_Ir w/ resids: \n"+str(l_Ir_res))

    theta_del=theta_m_nn-np.matmul(Theta_r_nn,theta_Ir)
    #print("theta_del \n"+str(theta_del))
    gf.dict_update(theta_del_dict,time+corr+itstr,theta_del)

    theta_Ir_res=gf.th_Ir(Th_r=Theta_r_nn,th_m=theta_del)
    #theta_Ir_res,restd,rnktd,stdel=lstsq(a=Theta_r_nn,b=theta_del)
    #print("theta_Ir w/ resids: \n"+str(theta_Ir_res))


    l_Ir_final=l_r+l_Ir_res
    ast.write(str(l_Ir_final))
    #l_Ir_final=l_r
    #print("l_ir_shape")
    #print(l_Ir_final.shape)
    gf.dict_update(l_r_dict,time+corr+itstr,l_Ir_final)
    #print("final amps (log form) \n"+str(l_Ir_final))

    Ir_converted=np.zeros_like(l_Ir_final)

    #Ir_converted=np.vstack([10.0**x for x in l_Ir_final])
    for x in range(l_Ir_final.shape[0]):
        lir=l_Ir_final[x,0]
        Ir_converted[x,0]=10.0**lir

    #print("l_conv shape")
    #print(Ir_converted.shape)

    gf.dict_update(l_Ir_conv_dict,time+corr+itstr,Ir_converted)

    theta_Ir_final=theta_Ir+theta_Ir_res
    gf.dict_update(theta_Ir_dict,time+corr+itstr,theta_Ir_final)

#return theta_Ir_final, Ir_converted
'''
#scl.close()
#lme.close()

'''
plt.clf()
plt.imshow(script_L)
plt.show()
plt.savefig('beforemask.png',overwrite=True)
plt.clf()
plt.imshow(script_L_nn)
plt.show()
plt.savefig('aftermask.png',overwrite=True)
plt.clf()
'''

#print("l_r keys")
#print(l_r_dict.keys())
#print("l_del")
#print(l_del_dict['3537XX0'].shape)
#print("l_Ir_conv")
#print(l_Ir_conv_dict['3537XX0'].shape)
#print(l_r_dict['4XX0'].shape)

#ELO_diff=np.max(ELO)-np.min(ELO)
#nint=int(ELO_diff/240)
#ELO_range=np.linspace(np.min(ELO),np.max(ELO),nint)

#print("fkeys")
fkeys=tir_ext.keys()
tkeys=antinfo['timestamp']
#print(fkeys)

nonan=0
nemp=0
nbad=0

ndp=[]
nda=[]

tt=np.array(visdata['axis_info']['time_axis']['MJDseconds'],dtype=float)

phase_ext1=np.zeros((nant,1,ntimes),dtype=float)
amp_ext1=np.zeros((nant,1,ntimes),dtype=float)

for t_step in range(ntimes):
    for k in range(nant):
        amp_ext1[k,0,t_step]=lir_ext[str(t_step)+corr+itstr][k]
        if k==refant: 
            phase_ext1[k,0,t_step]=0
        elif k<refant:
            phase_ext1[k,0,t_step]=tir_ext[str(t_step)+corr+itstr][k]
        elif k>refant:
            phase_ext1[k,0,t_step]=tir_ext[str(t_step)+corr+itstr][k-1]

phase_ext=gf.nonan(phase_ext1)
amp_ext=gf.nonan(amp_ext1)

phase_cln1=np.zeros((nant,1,ntimes),dtype=float)
amp_cln1=np.zeros((nant,1,ntimes),dtype=float)

for t_step in range(ntimes):
    for k in range(nant):
        amp_cln1[k,0,t_step]=lir_cln[str(t_step)+corr+itstr][k]
        if k==refant: 
            phase_cln1[k,0,t_step]=0
        elif k<refant:
            phase_cln1[k,0,t_step]=tir_cln[str(t_step)+corr+itstr][k]
        elif k>refant:
            phase_cln1[k,0,t_step]=tir_cln[str(t_step)+corr+itstr][k-1]

phase_cln=gf.nonan(phase_cln1)
amp_cln=gf.nonan(amp_cln1)

avga=np.empty((nant,1,len(ELO_range)-1))
avgp=np.empty((nant,1,len(ELO_range)-1))

phase_arr=phase_ext-phase_cln
amp_arr=amp_ext/amp_cln
#phase_arr=phase_cln
#amp_arr=amp_cln

'''
phase_arr=np.zeros((nant,1,ntimes),dtype=float)
amp_arr=np.zeros((nant,1,ntimes),dtype=float)

for t_step in range(ntimes):
    for k in range(nant):
        amp_arr[k,0,t_step]=l_Ir_conv_dict[str(t_step)+corr+itstr][k]
        if k==refant: 
            phase_arr[k,0,t_step]=0
        elif k<refant:
            phase_arr[k,0,t_step]=theta_Ir_dict[str(t_step)+corr+itstr][k]
        elif k>refant:
            phase_arr[k,0,t_step]=theta_Ir_dict[str(t_step)+corr+itstr][k-1]
'''
antinfo['phase_'+corr+itstr]=np.squeeze(phase_arr)
antinfo['amp_'+corr+itstr]=np.squeeze(amp_arr)
plt.clf()

'''
for t_step in range(len(ELO_range)-1):
    nd=0
    for x in range(len(tkeys)):
        if tkeys[x]>=ELO_range[t_step] and tkeys[x]<ELO_range[t_step+1]: nd+=1
    #nd=len([antinfo['timestamp'][x] for x in range(len(tkeys)) if (tkeys[x]>=ELO_range[t_step] and tkeys[x]<=ELO_range[t_step+1])])
    plt.scatter(ELO_range[t_step],nd)
    #print("nd\n",nd)
    for k in range(nant):
        #print("ELO_range")
        #print(ELO_range[t_step])
        #print(ELO_range[t_step+1])
        #nda=len([l_Ir_conv_dict[str(t_step)+corr+itstr][k] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])
        #plt.scatter(ELO_range[t_step],nda)
        #print("nda",nda)
        avga[k,0,t_step]=np.nanmean([l_Ir_conv_dict[str(t_step)+corr+itstr][k] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])
        if k==refant and nd!=0: 
            avgp[k,0,t_step]=0
        elif k<refant and nd!=0:
            avgp[k,0,t_step]=np.nanmean([theta_Ir_dict[str(t_step)+corr+itstr][k] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])
        elif k>refant and nd!=0:
            avgp[k,0,t_step]=np.nanmean([theta_Ir_dict[str(t_step)+corr+itstr][k-1] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])

    if np.isnan(np.nanmean(avga[:,0,t_step]))==True or np.isnan(np.nanmean(avgp[:,0,t_step]))==True:
        if nd==0:
            #print("That's okay!")
            avga[:,0,t_step]=np.nan
            avgp[:,0,t_step]=np.nan
            nemp+=1
        else: 
            nbad+=1
    else: nonan+=1

    #plt.savefig('./dplots/ndtstep_ant'+str(k)+'_'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'.png')
    #plt.clf()
plt.savefig('./dplots/ndtstep_'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'.png')
'''
plt.clf()

'''
for k in range(nant):
    #for t_step in range(len(ELO_range)-1):
        tlist=[]
        amplist=[]
        phlist=[]
        #for x in antinfo['timestamp']:
            #if x>ELO_range[t_step] and x<=ELO_range[t_step+1]: 
                tlist.append(x)
                ind=np.where(antinfo['timestamp']==x)
                amplist.append(antinfo['amp_'+corr+itstr][k][ind])
                phlist.append(antinfo['phase_'+corr+itstr][k][ind])
            #print("ind",ind)
        #print("amplist",amplist)
        nda=len(amplist)
        plt.scatter(ELO_range[t_step],nda)
        avga[k,0,t_step]=ma.mean(amplist)
        avgp[k,0,t_step]=ma.mean(phlist)

        if np.isnan(ma.mean(amplist))==True or np.isnan(ma.mean(phlist))==True:
            if nda==0:
                print("That's okay!")
                avga[:,0,t_step]=np.nan
                avgp[:,0,t_step]=np.nan
                nemp+=1
            else: 
                nbad+=1
        else: nonan+=1


        #if k==refant and nda!=0: avgp[k,0,t_step]=0
        #elif k!=refant and nda!=0: avgp[k,0,t_step]=np.nanmean(phlist)
        #elif k<refant and nda!=0: avgp[k,0,t_step]=np.nanmean(phlist)
    plt.savefig('./dplots/ndtstep_ant'+str(k)+'_'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'1.png')
    plt.clf()
'''

antinfo['avg_phase_'+corr+itstr]=np.squeeze(phase_arr)
antinfo['avg_amp_'+corr+itstr]=np.squeeze(amp_arr)    
'''
for k in range(nant):
    for t_step in range(len(ELO_range)-1):
        #print("ELO_range")
        #print(ELO_range[t_step])
        #print(ELO_range[t_step+1])
        nda=len([antinfo['amp_'+corr+itstr][k] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])
        plt.scatter(ELO_range[t_step],nda)
        ndb=len([lir_ext[str(t_step)+corr+itstr][k] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])
        plt.scatter(ELO_range[t_step],nda)

    #plt.savefig('./dplots/ndtstep_ant'+str(k)+'_'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'1.png')
    plt.clf()
plt.clf()

#plt.clf()
for a in range(nant):
    plt.scatter(ELO_range[:-1],avga[a,0,:],color='black',marker='x')
    plt.scatter(ELO_range[:-1],avgp[a,0,:],color='red')
plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_ap_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)

plt.clf()

for a in range(nant):
    plt.scatter(ELO_range[:-1],avga[a,0,:],color='black',marker='x')
    plt.savefig('./dplots/antgains/amp/ant'+str(a)+'_gampavg_'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)
    #plt.ylim(bottom=0.0,top=3.5)
    plt.clf()
plt.clf()


for a in range(nant):
    plt.scatter(ELO_range[:-1],avga[a,0,:],color='black',marker='x')
plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_a_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)


plt.clf()

for a in range(nant):
    plt.scatter(ELO_range[:-1],avgp[a,0,:],color='red')
    plt.savefig('./dplots/antgains/phase/ant'+str(a)+'_gphavg_'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)
plt.clf()


for a in range(nant):
    plt.scatter(ELO_range[:-1],avgp[a,0,:],color='red')
plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_p_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)

plt.clf()
'''
for a in range(nant):
    plt.scatter(antinfo['timestamp'],antinfo['amp_'+corr+itstr][a],color='black',marker='x')
    plt.show()
    plt.savefig('./dplots/antgains/amp/ant'+str(a)+'_gamp'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)
    #plt.show()
    #plt.ylim(bottom=0.0,top=3.5)
    plt.clf()
plt.clf()

for a in range(nant):
    plt.scatter(antinfo['timestamp'],antinfo['amp_'+corr+itstr][a],color='black',marker='x')
plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_ampdi_'+auth+'_'+case+corr+itstr+'1.png')

plt.clf()

for a in range(nant):
    plt.scatter(antinfo['timestamp'],antinfo['phase_'+corr+itstr][a],color='red')
    plt.show()
    plt.savefig('./dplots/antgains/phase/ant'+str(a)+'_gph'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)
    #plt.show()
    plt.clf()
plt.clf()

for a in range(nant):
    plt.scatter(antinfo['timestamp'],antinfo['phase_'+corr+itstr][a],color='red')
plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_phasedi_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)

plt.clf()


print("Intervals with pts:",nonan)
#print("Total intervals:",len(ELO_range)-1)
print("Bad baseline/antenna cases",nbad)
print("Empty bins",nemp)
print("Total points:",len(antinfo['timestamp']))


#newpt=np.zeros_like(xx)
#newres=np.zeros_like(xx)
tb1=0

#Opening data and pulling necessary info
ms.open(dvis,nomodify=True)
ms.selectinit(reset=True)
visdata1 = ms.getdata(['antenna1','antenna2','data','axis_info','flag'],ifraxis=True)

#Squeeze data then close ms

visdata1['data'] = np.squeeze(visdata1['data'])
ms.close()

antinfo1,xx1=antfunc(dfile=dvis,it=it,pol=pol,date=date)

newpt=np.zeros_like(xx1)

#for t_step in range(len(ELO_range)-1):
for t_step in range(len(antinfo['timestamp'])):
    #tchunk=0
    #listt=[]
    #for t in antinfo['timestamp']:
        #if t>=ELO_range[t_step] and t<ELO_range[t_step+1]: 
            #tchunk+=1
            #listt.append(t)
    #[antinfo['timestamp'][x] for x in range(len(tkeys)) if (tkeys[x]>=ELO_range[t_step] and tkeys[x]<=ELO_range[t_step+1])]
    #for i in range(len(tchunk)):
    #for i in range(tchunk):
        #thistime=(visdata1['axis_info']['time_axis']['MJDseconds']==listt[i])
        thistime=(visdata1['axis_info']['time_axis']['MJDseconds']==antinfo['timestamp'][t_step])
        time=t_step
        #indt=np.where(visdata1['axis_info']['time_axis']['MJDseconds']==listt[i])

        gt=antinfo['goodt_'+corr+itstr][t_step]

        #Pulling antennas that are bad for this time
        #Calculating number of antennas and baselines
        #-1 is to account for bad ants (ant1=-1 is not valid) 
        #need to rework how pt is pulled
        nb1=0
        for ant1 in np.unique(visdata1['antenna1']):
            for ant2 in np.unique(visdata1['antenna2']):
                if ant1 < ant2:
                    thisbase = (visdata1['antenna1']==ant1) & (visdata1['antenna2']==ant2)
                    iant1=np.where(antlist==ant1)[0][0]
                    iant2=np.where(antlist==ant2)[0][0]
                    if thisbase.sum()>0:
                        #pt=xx[thisbase][0][tb1]
                        pt1=xx1[thisbase][0][tb1]
                        ga1=antinfo['avg_amp_'+corr+itstr][iant1][t_step]
                        ga2=antinfo['avg_amp_'+corr+itstr][iant2][t_step]
                        if iant1==iref: 
                            gp1=0.
                        else:
                            gp1=antinfo['avg_phase_'+corr+itstr][iant1][t_step]
                        if iant2==iref:
                            gp2=0.
                        else:
                            gp2=antinfo['avg_phase_'+corr+itstr][iant2][t_step]
                        #newres[nb1,tb1]=pt*(ga1*ga2*np.exp(1.0j*(gp1-gp2)))
                        newpt[nb1,tb1]=pt1/(ga1*ga2*np.exp(1.0j*(gp1-gp2)))
                        nb1+=1
        tb1+=1
#ast.write("TB1\n")
#ast.write(str(tb1))

'''
for t in range(nb1):
    oldp=np.array(xx[t])
    newp=np.array(newres[t])
    #newv=np.array(newvis[t])
    plt.plot(tt,np.abs(oldp),".",c='b')
    plt.plot(tt,np.abs(newp),".",c='r')
#plt.ylim(top=3.0,bottom=0.0)
plt.show()
plt.savefig("./dplots/resplot_"+auth+date+corr+itstr+".png",overwrite=True)
plt.clf()
'''

for t in range(nb1):
    oldp=np.array(xx1[t])
    #newp=np.array(newpt[t])
    newv=np.array(newpt[t])
    plt.plot(tt,np.abs(oldp),".",c='b')
    plt.plot(tt,np.abs(newv),".",c='r')
plt.show()
plt.savefig("./dplots/visplot_"+auth+date+corr+itstr+"1.png",overwrite=True)
plt.clf()


#gf.dict_update(antinfo,'newres_'+corr+itstr,newpt)
gf.dict_update(antinfo1,'newvis_'+corr+itstr,newpt)

#ast.close()
#return newpt,antinfo1
#newres,x,antinfo,x
