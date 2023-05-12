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
from math import e
import pickle
import pdb
from multiprocessing import Process
from contextlib import closing
import math

attn=0
MJDdate=58232.00000000*3600.*24.
#MJDdate=57850.00000000 #apr072017

##################################################################

pol=0
corr='XX'
#datams1='sgr_apr07_flag_ext3_it00.ms'
#datams1='sgr_apr07_flag_ext3_it00avg240s_it00.ms'
#datams_cm='sgr_apr07_flag_ext3_it00_cmsplit.ms'
#datams_cm='sgr_apr07_flag_cmodel_avg240s_it00.ms'
target='sgr_apr24'
case='ref6'
refmeth='manual'
date='05042023'
it=0
#dvis1='sgr_apr24_flag_XXYY.ms'
dvis1='sgr_apr24_final_XXYY.ms'
dvis='sgr_apr24_final_XXYY_rflag_avg60s.ms'
#dvis='sgr_apr24_mid_rdata_avg60s.ms'
datams1=dvis[:-3]+'_ext3_it00.ms'
datams_cm=dvis[:-3]+'_cln_it00.ms'
rfant=6
auth='Jasmin'

####################################################################



#with open('data/model_redshifts.pickle', 'wb') as f:
#    pickle.dump(z_dict, f)
#with open('data/model_redshifts.pickle', 'rb') as f:
#    redshifts = pickle.load(f)


#try:
#    import uvmultimodel
#except ModuleNotFoundError:
import sys
sys.path.insert(0,"/home/jasminwash/.local/lib/python2.7/site-packages/NordicARC-3.0.5-py2.7-linux-x86_64.egg")
#import uvmultimodel


from NordicARC import uvmultifit as uvm

def antfunc(dfile,it,pol,date,case):
    itstr=str(it)

    #with open('data/'+case+'Data.pickle', 'rb') as f:
    #    visdata = pickle.load(f)
    tb.open(dfile)
    tbData1=tb.getcol('DATA')
    tbData=np.squeeze(tbData1)
    tbFlag1=tb.getcol('FLAG')
    tbFlag=np.squeeze(tbFlag1)
    tbt=tb.getcol('TIME')
    tba1=tb.getcol('ANTENNA1')
    tba2=tb.getcol('ANTENNA2')
    tbw=tb.getcol('WEIGHT')
    tb.close()

    #ms.open(dfile,nomodify=True)
    #visdata = ms.getdata(['antenna1','antenna2','flag','data','data_desc_id','sigma','axis_info'],ifraxis=True)
    #visdata['data'] = np.squeeze(visdata['data'])
    #ms.close()

    #allants=np.concatenate((visdata['antenna1'],visdata['antenna2']))
    allants=np.concatenate((tba1,tba2))
    antlist=np.unique(allants)

    #npol=visdata['data'].shape[0]
    npol=tbData.shape[0]
    #polnames=visdata['axis_info']['corr_axis']
    polnames=['XX','YY']
    corr=polnames[pol]
    nant=len(antlist)
    nbl=int(nant*(nant-1)/2)
    #ntimes=len(visdata['axis_info']['time_axis']['MJDseconds'])
    ntimes=len(tbt)

    #tt=(np.array(visdata['axis_info']['time_axis']['MJDseconds']-MJDdate,dtype=float)-MJDdate)/3600.
    tt=(np.array(tbt-MJDdate,dtype=float))/3600.

    #Creating a tf dictionary for each ant and whether time is good or bad there
    antinfo=dict()
    #antinfo['timestamp']=visdata['axis_info']['time_axis']['MJDseconds']
    antinfo['timestamp']=tbt

    allant=np.zeros((nant,ntimes),dtype=bool)
    alltim=np.zeros((ntimes,nant),dtype=bool)

    #xx=np.zeros_like(visdata['flag'][pol])
    xx=np.zeros_like(tbData[pol])
    xxw=np.zeros_like(tbw[pol])
    #flagpt=((visdata['flag'][pol]==True)|(visdata['data'][pol]==0))
    flagpt=((tbFlag[pol]==True)|(tbData[pol]==0))
    #flagpt1=(visdata['data'][pol]==0)
    #xx0=np.squeeze(np.where(flagpt==True,np.nan,visdata['data'][pol]))
    #xx=np.squeeze(np.where(flagpt==True,np.nan,visdata['data'][pol]))
    xx=np.squeeze(np.where(flagpt==True,np.nan,tbData[pol]))
    xxw=np.squeeze(np.where(flagpt==True,np.nan,tbw[pol]))
    print("xxshape",xx.shape)

    #for y in range(nbl):
    #    pp1=np.array(xx[y])
    #    plt.scatter(tt,pp)
    plt.scatter(tt,np.abs(xx))
    plt.show()
    #plt.xlim(left=2.5,right=13.5)
    #plt.savefig('datatestplot_'+auth+refmeth+'_'+date+corr+str(it)+'.png',overwrite=True)
    plt.savefig('datatestplot_'+str(dfile[:-3])+'_'+str(it)+date+'1.png',overwrite=True)
    plt.clf()
    '''
    #tf matrix for where amp of ant is bad (True)
    tfdata=np.isnan(np.log10(np.abs(xx)))
    #np.log10(np.abs(xx))
    print(tfdata)

    #cycle through each ant
    #a=0
    for a in range(nant):
        ant=antlist[a]
        #all baselines w/ this ant
        #visdata['antenna1']
        thisant=(tba1==ant) | (tba2==ant)
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
        #ct=visdata['axis_info']['time_axis']['MJDseconds'][i]
        ct=tbt[i]
        #print("tstamp",ct)
        for ant in range(nant):
            tfstore[j]=allant[ant][i]
            j+=1
        alltim[i,:]=tfstore


    antinfo['goodt_'+corr+itstr]=alltim
    '''
    #return antinfo,xx
    return xx,xxw



####################################################################
def gsolve(nbl,nant,vdata,varr,varrcln,antlist,corr,itstr,refant,tba1,tba2,tbt,tbspw,vwt):#ainfo
    #gst0=time.time()
    plt.clf()
    #theta_r_dict={}
    #theta_m_dict={}
    #theta_Ir_dict={}
    #theta_del_dict={}
    #theta_rms_dict={}

    #script_L_dict={}
    #script_L_nn_dict={}
    #theta_r_nn_dict={}
    #l_m_dict={}
    #l_m_nn_dict={}
    #theta_m_nn_dict={}
    #l_r_dict={}
    #l_del_dict={}
    #l_rms_dict={}
    #l_Ir_conv_dict={}
    #all_bdance={}
    #all_bdance_ph={}
    vism_dict={}
    jac_dict={}
    wt_dict={}
    vism_nn_dict={}
    jac_nn_dict={}
    wt_nn_dict={}
    vis_r_dict={}


    smc=0

    rfantind=np.where(antlist==refant)[0][0]

    

    untime=np.unique(tbt)
    for s in range(nspw):
        spwstr=str(s)
        print(spwstr)
        #for j in range(len(ainfo['timestamp'])):
        for j in range(len(untime)):
            allt=(tbt==untime[j])
            #dTime=vdata[pol][allt]
            a1Time=tba1[allt]
            a2Time=tba2[allt]
            aAllTime=np.concatenate((a1Time,a2Time))
            #print(aAllTime)
            exclAnts1=np.array(list(set(antlist)-set(aAllTime)))
            #print(exclAnts1)
            goodAnts=np.sort(np.array(list(set(antlist)-set(exclAnts1))))
            #print(goodAnts)
            extraAnts=gf.bl_checker(goodAnts,varr,tba1,tba2,tbt,tbspw,untime[j],s)
            #may need to add here
            #print(extraAnts)
            exclAnts=np.unique(np.concatenate((exclAnts1,extraAnts))).astype(int)
            #print(exclAnts)
            #print(len(exclAnts))
            #pdb.set_trace()
            #nant=len(np.unique(aAllTime))
            #nbl=(nant*(nant-1))/2

            #Jacobian and measured phase matrix
            #Theta_r=np.zeros((nbl,nant-1),dtype=int)
            #bad_T=np.zeros((nbl,nant-1),dtype=int)
            #theta_m=np.zeros((nbl,1),dtype=float)

            #Defining Jacobian and measured amp matrix
            #script_L=np.zeros((nbl,nant),dtype=int)
            #bad_L=np.zeros((nbl,nant),dtype=int)
            #l_m=np.zeros((nbl,1),dtype=float)
            #astore=np.zeros((nbl,1),dtype=float)
            #bantstore=[]
            #blstore=[]

            #New Jacobian and vis chunks
            jac=np.zeros((2*nbl,nant),dtype=int)
            vism=np.zeros((2*nbl,1),dtype=np.complex128)
            wt=np.zeros((2*nbl,2*nbl),dtype=np.complex128)

            #print(time)
            #thistime=(vdata['axis_info']['time_axis']['MJDseconds']==ainfo['timestamp'][j])
            time=str(j)

            #Calculating number of antennas and baselines
            #-1 is to account for bad ants (ant1=-1 is not valid)    
            nb=0
            #for ant1 in np.unique(vdata['antenna1']):
            for ant1 in np.unique(tba1):
                for ant2 in np.unique(tba2):
                #for ant2 in np.unique(vdata['antenna2']):
                    if ant1 < ant2:
                        #thisbase = (vdata['antenna1']==ant1) & (vdata['antenna2']==ant2)
                        thispt= (tba1==ant1) & (tba2==ant2) & (tbt==untime[j]) & (tbspw==s)
                        iant1=np.where(antlist==ant1)[0][0]
                        #print(iant1)
                        iant2=np.where(antlist==ant2)[0][0]
                        #print(iant2)
                        if thispt.sum()==0 or ant1 in exclAnts or ant2 in exclAnts:
                            #l_m[nb,0]=np.nan
                            #astore[nb,0]=np.nan
                            vism[nb,0]=np.nan
                            vism[nb+nbl,0]=np.nan

                            #script_L[nb,iant1]=1#1
                            #script_L[nb,iant2]=1#1
                            jac[nb,iant1]=1
                            jac[nb,iant2]=1


                            #PHASE STUFF
                            #theta_m[nb,0]=np.nan

                            if ant1==refant: jac[nb+nbl,iant2]=-1#-1
                            if ant2==refant: jac[nb+nbl,iant1]=1#1
                            #if ant1!=refant and ant1>refant:
                            #    Theta_r[nb,iant1-1]=1#1
                            #    Theta_r[nb,iant2-1]=-1#-1
                            else:
                                jac[nb+nbl,iant1]=1#1
                                jac[nb+nbl,iant2]=-1#-1
                            #if (ant1!=refant and (ant2>refant and ant1<refant)):
                            #    Theta_r[nb,iant1]=1#1
                            #    Theta_r[nb,iant2-1]=-1
                            wt[nb,nb]=np.nan
                            wt[nb+nbl,nb+nbl]=np.nan
                            nb+=1
                            continue
                        if thispt.sum()>0:
                            #potential 0 after thisbase and thistime
                            #pt=varr[thisbase][0][j]
                            pt=varr[thispt][0]
                            ptc=varrcln[thispt][0]
                            ptw=vwt[thispt][0]
                            #print("PT",pt)
                            #print("PTC",ptc)
                            #pdb.set_trace()
                            #ph=np.angle(pt,deg=False)
                            #phc=np.angle(ptc,deg=False)
                            #amp=abs(pt)
                            #ampc=abs(ptc)
                            rl=pt.real
                            rlc=ptc.real
                            im=1.0j*pt.imag
                            imc=1.0j*ptc.imag
                            #pdb.set_trace()
                            #if np.isnan(amp)==True: 
                            #    bantstore.append(iant1)
                            #    blstore.append(nb)
                            #    #bantstore.append(iant2)
                            #    bad_L[nb,iant1]=2
                            #    bad_L[nb,iant2]=2

                            #    if ant1==refant: bad_T[nb,iant2-1]=2
                            #    if ant2==refant: bad_T[nb,iant1]=2
                            #    if ant1!=refant and ant1>refant:
                            #        bad_T[nb,iant1-1]=2
                            #        bad_T[nb,iant2-1]=2
                            #    if ant1!=refant and ant2<refant:
                            #        bad_T[nb,iant1]=2
                            #        bad_T[nb,iant2]=2
                            #    if (ant1!=refant and (ant2>refant and ant1<refant)):
                            #        bad_T[nb,iant1]=2
                            #        bad_T[nb,iant2-1]=2
                            #print("AMP",amp)

                            #if iant1==5 or iant2==5:
                                #amp*=4
                                #ph+=90
                                #print("bing")

                            vism[nb,0]=rl-rlc
                            vism[nb+nbl,0]=im-imc
                            #astore[nb,0]=amp

                            jac[nb,iant1]=1
                            jac[nb,iant2]=1

                            #PHASE STUFF
                            #theta_m[nb,0]=ph-phc

                            if ant1==refant: jac[nb+nbl,iant2]=-1
                            if ant2==refant: jac[nb+nbl,iant1]=1
                            else:
                                jac[nb+nbl,iant1]=1
                                jac[nb+nbl,iant2]=-1
                            
                            wt[nb,nb]=ptw
                            wt[nb+nbl,nb+nbl]=ptw
                                
                            nb+=1
            #print(j)
            #badance=np.full((nant,1),False,dtype=bool)
            #badbl=np.full((nbl,1),False,dtype=bool)
            #for ia in range(nant):
            #    #print(bantstore.count(ia))
            #    #if ia==38: 
            #    #    if bantstore.count(ia)>0: badance[ia,0]=True
            #    #else:
            #    #if bantstore.count(ia)==nant-1: badance[ia,0]=True
            #    #if ia in bantstore: badance[ia,0]=True
            #    #if bantstore.count(ia)==(nant-1)-ia: 
            #    #    badance[ia,0]=True
            #    #else: print("bonk!")
            #for ib in range(nbl):
            #    if ib in blstore: badbl[ib,0]=True
            #print(len(blstore))
            #pdb.set_trace()    
            #badance=np.unique(bantstore)
            #print(badance)
            #plt.scatter(j,len(np.where(badance==True)[0]))
            #gf.dict_update(all_bdance,time+corr+itstr+spwstr,badance)


            #gf.dict_update(theta_r_dict,time+corr+itstr+spwstr,Theta_r)
            #gf.dict_update(theta_m_dict,time+corr+itstr+spwstr,theta_m)
            gf.dict_update(vism_dict,time+corr+itstr+spwstr,vism)
            gf.dict_update(jac_dict,time+corr+itstr+spwstr,jac)
            gf.dict_update(wt_dict,time+corr+itstr+spwstr,wt)

            #scl.write("\n\n")

            #lme.write("\n\n")

            #ast.write(str(astore))
            #ast.write("\n\n")
            #help me mark I hope
            #pdb.set_trace()

            #script_L_nn=gf.nonan(script_L)
            #gf.dict_update(script_L_dict,time+corr+itstr+spwstr,script_L)
            #scl.write(str(script_L_nn))
            #l_m_nn=gf.nonan(l_m)
            #gf.dict_update(l_m_dict,time+corr+itstr+spwstr,l_m)
            #lme.write(str(l_m_nn))

            #gf.dict_update(script_L_nn_dict,time+corr+itstr+spwstr,script_L_nn)
            #gf.dict_update(l_m_nn_dict,time+corr+itstr+spwstr,l_m_nn)


            vism_nn=gf.nonan(vism)
            jac_nn=gf.nonan(jac)
            wt_nn=gf.nonan(wt)

            gf.dict_update(vism_nn_dict,time+corr+itstr+spwstr,vism_nn)
            gf.dict_update(jac_nn_dict,time+corr+itstr+spwstr,jac_nn)
            gf.dict_update(wt_nn_dict,time+corr+itstr+spwstr,wt_nn)

            #l_r=gf.l_Ir(ll_r=script_L_nn,ll_m=l_m_nn)
            #print("wherelmnnfalse",len(np.where(l_m_nn.mask==False)[0]))
            #print("j",j)
            vis_r=gf.realim_mask(vism=vism_nn,jac=jac_nn,nant=nant,antlist=antlist,exclAnts=exclAnts,wt=wt_nn)
            #l_r=gf.l_Ir_nan(ll_r=script_L,ll_m=l_m,ainfo=ainfo,t=j,nant=nant,corr=corr,itstr=itstr,bdance=badance,bbl=badbl)

            gf.dict_update(vis_r_dict,time+corr+itstr+spwstr,vis_r)

            #pdb.set_trace()

            #Theta_r_nn=gf.nonan(Theta_r)
            #gf.dict_update(theta_r_nn_dict,time+corr+itstr+spwstr,Theta_r_nn)
            #theta_m_nn=gf.nonan(theta_m)
            #gf.dict_update(theta_m_nn_dict,time+corr+itstr+spwstr,theta_m_nn)


            #theta_Ir=gf.th_Ir_mask(Th_r=Theta_r_nn,th_m=theta_m_nn,rfantind=rfantind,nant=nant,bbl=badbl,antlist=antlist,exclAnts=exclAnts)
            #theta_Ir=gf.th_Ir(Th_r=Theta_r_nn,th_m=theta_m_nn)



            #Residuals
            #l_del=l_m_nn-np.matmul(script_L_nn,l_r)
            #print("l_del \n"+str(l_del))
            #gf.dict_update(l_del_dict,time+corr+itstr+spwstr,l_del)

            #print(l_del.shape)

            #l_Ir_res=gf.l_Ir_mask(ll_r=script_L_nn,ll_m=l_del,nant=nant,bbl=badbl,antlist=antlist,exclAnts=exclAnts)

            #print("l_Ir w/ resids: \n"+str(l_Ir_res))

            #theta_del=theta_m_nn-ma.dot(Theta_r_nn,theta_Ir)
            #print("theta_del \n"+str(theta_del))
            #gf.dict_update(theta_del_dict,time+corr+itstr+spwstr,theta_del)

            #theta_Ir_res=gf.th_Ir_mask(Th_r=Theta_r_nn,th_m=theta_del,rfantind=rfantind,nant=nant,bbl=badbl,antlist=antlist,exclAnts=exclAnts)
            #theta_Ir_res=gf.th_Ir(Th_r=Theta_r_nn,th_m=theta_del)
            #theta_Ir_res,restd,rnktd,stdel=lstsq(a=Theta_r_nn,b=theta_del)
            #print("theta_Ir w/ resids: \n"+str(theta_Ir_res))


            #l_Ir_final=l_r+l_Ir_res
            #ast.write(str(l_Ir_final))
            #l_Ir_final=l_r
            #print("l_ir_shape")
            #print(l_Ir_final.shape)
            #gf.dict_update(l_r_dict,time+corr+itstr+spwstr,l_Ir_final)
            #print("final amps (log form) \n"+str(l_Ir_final))

            #with open('data/model_redshifts.pickle', 'rb') as f:
            #    redshifts = pickle.load(f)

            #Ir_converted=np.zeros_like(l_Ir_final)

            #Ir_converted=np.vstack([10.0**x for x in l_Ir_final])
            #for x in range(l_Ir_final.shape[0]):
            #    lir=l_Ir_final[x,0]
            #    Ir_converted[x,0]=10.0**lir

            #print("l_conv shape")
            #print(Ir_converted.shape)

            #gf.dict_update(l_Ir_conv_dict,time+corr+itstr+spwstr,Ir_converted)

            #theta_Ir_final=theta_Ir+theta_Ir_res
            #theta_Ir_final=theta_Ir
            #gf.dict_update(theta_Ir_dict,time+corr+itstr+spwstr,theta_Ir_final)


          #plt.show()
    '''
    with open('data/all_bdance_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(all_bdance, f)
    with open('data/theta_Ir_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(theta_Ir_dict, f)
    with open('data/l_Ir_conv_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(l_Ir_conv_dict, f)
    with open('data/l_r_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(l_r_dict, f)
    with open('data/theta_del_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(theta_del_dict, f)
    with open('data/theta_r_nn_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(theta_r_nn_dict, f)
    with open('data/theta_m_nn_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(theta_m_nn_dict, f)
    with open('data/l_del_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(l_del_dict, f)
    with open('data/theta_r_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(theta_r_dict, f)
    with open('data/theta_m_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(theta_m_dict, f)
    with open('data/script_L_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(script_L_dict, f)
    with open('data/l_m_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(l_m_dict, f)
    with open('data/script_L_nn_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(script_L_nn_dict, f)
    with open('data/l_m_nn_dict_'+case+'.pickle', 'wb') as f:
        pickle.dump(l_m_nn_dict, f)
    '''
    with open('data/vism_nn_dict.pickle', 'wb') as f:
        pickle.dump(vism_nn_dict, f)
    with open('data/jac_nn_dict.pickle', 'wb') as f:
        pickle.dump(jac_nn_dict, f)
    with open('data/wt_nn_dict.pickle', 'wb') as f:
        pickle.dump(wt_nn_dict, f)
    with open('data/vis_r_dict.pickle', 'wb') as f:
        pickle.dump(vis_r_dict, f)    


    plt.savefig('validants.png',overwrite=True)
    plt.show()
    plt.clf()
    #gst1=time.time()
    #print("Time to run is",(gst1-gst0)/60,"minutes",(gst1-gst0)/3600,"hours")
    return 0


#def lls(pol,corr,datams1,target,case,auth,refmeth,date,it,dvis,rfant):
itstr=str(it)
nspw=4

antinfo={}

print(datams1)

#Opening data and pulling necessary info
#ms.open(datams1,nomodify=True)
#ms.selectinit(reset=True)
#visdata = ms.getdata(['antenna1','antenna2','data','axis_info','flag'],ifraxis=True)
#printing correlations
#print(visdata['axis_info']['corr_axis'])

#Squeeze data then close ms

#visdata['data'] = np.squeeze(visdata['data'])
#print("data shape",visdata['data'].shape)
#ms.close()

#with open('data/extData.pickle', 'rb') as f:
#    visdata=pickle.load(f)

tb.open(datams1)
tbData1=tb.getcol('DATA')
tbData=np.squeeze(tbData1)
tbFlag1=tb.getcol('FLAG')
tbFlag=np.squeeze(tbFlag1)
tbt=tb.getcol('TIME')
tba1=tb.getcol('ANTENNA1')
tba2=tb.getcol('ANTENNA2')
tbspw=tb.getcol('DATA_DESC_ID')
tb.close()

utbt=np.unique(tbt)

allants=np.concatenate((tba1,tba2))
antlist=np.unique(allants)

#npol=visdata['data'].shape[0]
#polnames=visdata['axis_info']['corr_axis']
#corr=polnames[pol]
#nant=len(antlist)
#nbl=int(nant*(nant-1)/2)
#ntimes=len(visdata['axis_info']['time_axis']['MJDseconds'])

npol=tbData.shape[0]
#polnames=visdata['axis_info']['corr_axis']
polnames=['XX','YY']
corr=polnames[pol]
nant=len(antlist)
nbl=int(nant*(nant-1)/2)
#ntimes=len(visdata['axis_info']['time_axis']['MJDseconds'])
ntimes=len(tbt)


#pdb.set_trace()

#Pulling all unique timestamps, print its length and then all times
#ELO=np.unique(visdata['axis_info']['time_axis']['MJDseconds'])
ELO=np.unique(tbt)
print(len(ELO))
#print(len(visdata['axis_info']['time_axis']['MJDseconds']))
#print([x-ELO[0] for x in ELO])

#antinfo,
xx,xxw=antfunc(dfile=datams1,it=it,pol=pol,date=date,case='ext')
#antinfo_clean,
xx_clean,xxw_clean=antfunc(dfile=datams_cm,it=it,pol=pol,date=date,case='cln')
print("XX",xx.shape)


allgoodtime=0
allbad=0
good1=0

nb=0

nvis=[]

#Cycling through all the times
#for time in antinfo['goodt_'+corr+itstr+spwstr]:
#    if np.any(time)==False: allgoodtime+=1
#    if np.any(time)==True:
#        nbad=np.count_nonzero(time)
#        if nbad==nant: allbad+=1
#        if nbad<nant: good1+=1

#print("All good antennas:",allgoodtime)
#print("Almost all good antennas:",good1)
#print("All ants bad:",allbad)

#refant=gf.refantfinder(tba1=tba1,tba2=tba2,tbt=tbt,tbspw=tbspw,untime=utbt,nspw=nspw,nant=nant,antlist=antlist)
#iref=np.where(antlist==refant)[0]
refant=rfant
iref=np.where(antlist==refant)[0]
#pdb.set_trace()

#scl=open("scriptL_"+target+date+corr+itstr+".txt","w")
#lme=open("lm_"+target+date+corr+itstr+".txt","w")
#ast=open("astore_"+target+date+corr+itstr+".txt","w")

#tir_ext,lir_ext=
#gsolve(nbl=nbl,nant=nant,vdata=visdata,varr=xx,ainfo=antinfo,antlist=antlist,corr=corr,itstr=itstr,refant=refant,case='ext')

gsolve(nbl=nbl,nant=nant,vdata=tbData,varr=xx,varrcln=xx_clean,antlist=antlist,corr=corr,itstr=itstr,refant=refant,tba1=tba1,tba2=tba2,tbt=tbt,tbspw=tbspw,vwt=xxw)#ainfo=antinfo
#extargs=(nbl=nbl,nant=nant,vdata=tbData,varr=xx,ainfo=antinfo,antlist=antlist,corr=corr,itstr=itstr,refant=refant,case='ext',tba1=tba1,tba2=tba2,tbt=tbt,tbspw=tbspw),(nbl=nbl,nant=nant,vdata=tbData,varr=xx_clean,ainfo=antinfo_clean,antlist=antlist,corr=corr,itstr=itstr,refant=refant,case='clean',tba1=tba1,tba2=tba2,tbt=tbt,tbspw=tbspw,)
#clnargs=(nbl=nbl,nant=nant,vdata=tbData,varr=xx_clean,ainfo=antinfo_clean,antlist=antlist,corr=corr,itstr=itstr,refant=refant,case='clean',tba1=tba1,tba2=tba2,tbt=tbt,tbspw=tbspw,)
#dargs=(nbl,nant,tbData,xx,antinfo,antlist,corr,itstr,refant,'ext',tba1,tba2,tbt,tbspw,)
#dargs1=(nbl,nant,tbData,xx_clean,antinfo_clean,antlist,corr,itstr,refant,'cln',tba1,tba2,tbt,tbspw,)

'''
if __name__ == "__main__":
    p1 = Process(target=gsolve,args=dargs)
    p1.start()
    p2 = Process(target=gsolve,args=dargs1)
    p2.start()
    p1.join()
    p2.join()
'''

#with open('data/theta_Ir_dict_all.pickle', 'rb') as f:
#    tir_all = pickle.load(f)
#with open('data/l_Ir_conv_dict_all.pickle', 'rb') as f:
#    lir_all = pickle.load(f)
with open('data/vis_r_dict.pickle','rb') as f:
    vis_all=pickle.load(f)

#tir_cln,lir_cln=

#gsolve(nbl=nbl,nant=nant,vdata=tbData,varr=xx_clean,ainfo=antinfo_clean,antlist=antlist,corr=corr,itstr=itstr,refant=refant,case='cln',tba1=tba1,tba2=tba2,tbt=tbt,tbspw=tbspw)
#clnargs=[]
#with open('data/theta_Ir_dict_cln.pickle', 'rb') as f:
#    tir_cln = pickle.load(f)
#with open('data/l_Ir_conv_dict_cln.pickle', 'rb') as f:
#    lir_cln = pickle.load(f)




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
#fkeys=tir_all.keys()
#tkeys=antinfo['timestamp']
#print(fkeys)

#pdb.set_trace()

nonan=0
nemp=0
nbad=0

ndp=[]
nda=[]

utbt=np.unique(tbt)
#tt=(np.array(visdata['axis_info']['time_axis']['MJDseconds'],dtype=float)-MJDdate)/3600.
tt=(np.array(tbt)-MJDdate)/3600.


nutimes=len(np.unique(tbt))


phase_all1=np.zeros((nant,1,nutimes,nspw),dtype=float)
amp_all1=np.zeros((nant,1,nutimes,nspw),dtype=float)

for s in range(4):
    spwstr=str(s)
    for t_step in range(nutimes):
        for k in range(nant):
            amp_all1[k,0,t_step,s]=abs(vis_all[str(t_step)+corr+itstr+spwstr][k])
            if k==refant: 
                phase_all1[k,0,t_step,s]=0
            elif k<refant:
                phase_all1[k,0,t_step,s]=np.angle(vis_all[str(t_step)+corr+itstr+spwstr][k])
            elif k>refant:
                phase_all1[k,0,t_step,s]=np.angle(vis_all[str(t_step)+corr+itstr+spwstr][k-1])

#phase_arr=gf.nonan(phase_all1)
#amp_arr=gf.nonan(amp_all1)
'''
while np.any(np.abs(phase_ext[~np.isnan(phase_ext)])>180)==True:
    #for t in range(len(antinfo['timestamp'])):
    for t in range(nutimes):
        for a in range(len(antlist)):
            for s in range(4):
                wph=phase_ext[a,0,t,s]
                if np.isnan(wph)==False and wph>180: phase_ext[a,0,t,s]=wph-360
                if np.isnan(wph)==False and wph<-180: phase_ext[a,0,t,s]=wph+360
'''
#phase_cln1=np.zeros((nant,1,nutimes,nspw),dtype=float)
#amp_cln1=np.zeros((nant,1,nutimes,nspw),dtype=float)

'''
for s in range(nspw):
    spwstr=str(s)
    for t_step in range(nutimes):
        for k in range(nant):
            amp_cln1[k,0,t_step,s]=lir_cln[str(t_step)+corr+itstr+spwstr][k]
            if k==refant: 
                phase_cln1[k,0,t_step,s]=0
            elif k<refant:
                phase_cln1[k,0,t_step,s]=tir_cln[str(t_step)+corr+itstr+spwstr][k]
            elif k>refant:
                phase_cln1[k,0,t_step,s]=tir_cln[str(t_step)+corr+itstr+spwstr][k-1]

phase_cln=gf.nonan(phase_cln1)
amp_cln=gf.nonan(amp_cln1)
'''
'''
while np.any(np.abs(phase_cln[~np.isnan(phase_cln)])>180)==True:
    #for t in range(len(antinfo['timestamp'])):
    for t in range(nutimes):
        for a in range(len(antlist)):
            for s in range(4):
                wph=phase_cln[a,0,t,s]
                if np.isnan(wph)==False and wph>180: phase_cln[a,0,t,s]=wph-360
                if np.isnan(wph)==False and wph<-180: phase_cln[a,0,t,s]=wph+360
'''
#avga=np.empty((nant,1,len(ELO_range)-1))
#avgp=np.empty((nant,1,len(ELO_range)-1))

#phase_arr1=phase_ext-phase_cln
#amp_arr1=amp_ext*amp_cln

phase_arr=ma.masked_where(amp_all1>3,phase_all1)
amp_arr=ma.masked_where(amp_all1>3,amp_all1)

while np.any(np.abs(phase_arr[~np.isnan(phase_arr)])>math.pi)==True:
    #for t in range(len(antinfo['timestamp'])):
    for t in range(nutimes):
        for a in range(len(antlist)):
            for s in range(nspw):
                wph=phase_arr[a,0,t,s]
                if np.isnan(wph)==False and wph>math.pi: phase_arr[a,0,t,s]=wph-(2.*math.pi)#-360
                if np.isnan(wph)==False and wph<(-1.*math.pi): phase_arr[a,0,t,s]=wph+(2.*math.pi)#+360

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
    for s in range(nspw):
        spwstr=str(s)
        plt.scatter(utbt,antinfo['amp_'+corr+itstr][a,:,s],color='black',marker='x')
    plt.show()
    plt.savefig('./dplots/antgains/amp/ant'+str(a)+'_gamp'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)
    #plt.show()
    #plt.ylim(bottom=0.0,top=3.5)
    plt.clf()
plt.clf()

for s in range(nspw):
    spwstr=str(s)
    for a in range(nant):
        plt.scatter(utbt,antinfo['amp_'+corr+itstr][a,:,s],color='black',marker='x')
#plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_ampdi_'+auth+'_'+case+corr+itstr+'1.png')

plt.clf()

for a in range(nant):
    for s in range(nspw):
        spwstr=str(s)
        plt.scatter(utbt,antinfo['phase_'+corr+itstr][a,:,s],color='red')
    plt.show()
    plt.savefig('./dplots/antgains/phase/ant'+str(a)+'_gph'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)
    #plt.show()
    plt.clf()
plt.clf()

for s in range(nspw):
    spwstr=str(s)
    for a in range(nant):
        plt.scatter(utbt,(180./math.pi)*antinfo['phase_'+corr+itstr][a,:,s],color='red')
plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_phasedi_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)

plt.clf()

'''
for s in range(nspw):
    for a in range(nant):
        plt.scatter(utbt,phase_cln[a,0,:,s],color='red')
#plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_phaseclndi_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)

plt.clf()

for s in range(nspw):
    for a in range(nant):
        plt.scatter(utbt,phase_ext[a,0,:,s],color='red')
#plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_phaseextdi_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)

plt.clf()

for s in range(nspw):
    for a in range(nant):
        plt.scatter(utbt,amp_cln[a,0,:,s],color='black')
#plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_ampclndi_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)

plt.clf()

for s in range(nspw):
    for a in range(nant):
        plt.scatter(utbt,amp_ext[a,0,:,s],color='black')
#plt.legend()
plt.show()

plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_ampextdi_'+auth+'_'+case+corr+itstr+'1.png',overwrite=True)

plt.clf()
'''


#print("Intervals with pts:",nonan)
#print("Total intervals:",len(ELO_range)-1)
#print("Bad baseline/antenna cases",nbad)
#print("Empty bins",nemp)
#print("Total points:",len(antinfo['timestamp']))


#newpt=np.zeros_like(xx)
#newres=np.zeros_like(xx)
#tb1=0

#Opening data and pulling necessary info
#ms.open(dvis,nomodify=True)
#ms.selectinit(reset=True)
#visdata1 = ms.getdata(['antenna1','antenna2','data','axis_info','flag'],ifraxis=True)

#Squeeze data then close ms

#visdata1['data'] = np.squeeze(visdata1['data'])
#ms.close()

#with open('data/rawData.pickle', 'rb') as f:
#    visdata1 = pickle.load(f)

#HERE!
tb.open(dvis)
tbData1=tb.getcol('DATA')
tbDataV=np.squeeze(tbData1)
tbFlag1=tb.getcol('FLAG')
tbFlagV=np.squeeze(tbFlag1)
tbtV=tb.getcol('TIME')
tba1V=tb.getcol('ANTENNA1')
tba2V=tb.getcol('ANTENNA2')
tbVspw=tb.getcol('DATA_DESC_ID')
tb.close()

with open('data/antinfoIt'+str(it)+'.pickle', 'wb') as f:
    pickle.dump(antinfo, f)

xx1,xx1w=antfunc(dfile=dvis,it=it,pol=pol,date=date,case='raw')


newpt=np.zeros_like(xx1,dtype=np.complex128)
#newpt=np.zeros((nbl,len(utbt)),dtype=np.complex128)

#for t_step in range(len(ELO_range)-1):
for spw in range(nspw):
    spwstr=str(spw)
    for t_step in range(len(utbt)):
        #print(t_step)
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
        #thistime=(visdata1['axis_info']['time_axis']['MJDseconds']==antinfo['timestamp'][t_step])
        #thistime=(tbtV==antinfo['timestamp'][t_step])
        #thistime=(tbtV==utbt[t_step])
        time=t_step
        #indt=np.where(visdata1['axis_info']['time_axis']['MJDseconds']==listt[i])

        #gt=antinfo['goodt_'+corr+itstr+spwstr][t_step]

        #Pulling antennas that are bad for this time
        #Calculating number of antennas and baselines
        #-1 is to account for bad ants (ant1=-1 is not valid) 
        #need to rework how pt is pulled
        nb1=0
        for ant1 in np.unique(tba1V):
            for ant2 in np.unique(tba2V):
                if ant1 < ant2:
                    thispt = (tba1V==ant1) & (tba2V==ant2) & (tbtV==utbt[t_step]) & (tbVspw==spw)
                    #thistime = (tbt==t_step)
                    iant1=np.where(antlist==ant1)[0][0]
                    iant2=np.where(antlist==ant2)[0][0]
                    #if thispt.sum()==0: 
                    #    newpt[nb1,t_step]=np.nan
                    #    nb1+=1
                    if thispt.sum()>0:
                        #pt=xx[thisbase][0][tb1]
                        #pt1=xx1[thisbase][0][tb1]
                        pt1=xx1[thispt][0]
                        ga1=antinfo['avg_amp_'+corr+itstr][iant1][t_step][spw]
                        ga2=antinfo['avg_amp_'+corr+itstr][iant2][t_step][spw]
                        if iant1==iref: 
                            gp1=0.
                        else:
                            gp1=antinfo['avg_phase_'+corr+itstr][iant1][t_step][spw]
                            #gp1d=antinfo['avg_phase_'+corr+itstr][iant1][t_step][spw]
                            #gp1=radians(gp1d)
                        if iant2==iref:
                            gp2=0.
                        else:
                            gp2=antinfo['avg_phase_'+corr+itstr][iant2][t_step][spw]
                            #gp2d=antinfo['avg_phase_'+corr+itstr][iant2][t_step][spw]
                            #gp2=radians(gp2d)
                        #newres[nb1,tb1]=pt*(ga1*ga2*np.exp(1.0j*(gp1-gp2)))
                        #newpt[nb1,t_step]=pt1*(ga1*ga2*np.exp(1.0j*(gp1-gp2)))
                        #newpt[nb1,t_step]=pt1/(ga1*ga2*(e**(1.0j*(gp1-gp2))))
                        newpt[thispt]=pt1/(ga1*ga2*(e**(1.0j*(gp1-gp2))))
                        if np.abs(newpt[thispt])==0:newpt[thispt]=np.nan
                        #pdb.set_trace()
                        nb1+=1
    #tb1+=1
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
#for s in range(4):
#just need oldp then tbtV
'''
for t in range(len(tt)):
    thispt=(tbtV==utbt[t])
    oldp=np.array(xx1[thispt])
    #newp=np.array(newpt[t])
    #newv=np.array(newpt[t])
    plt.plot([tt[t]]*len(oldp),np.abs(oldp),".",c='b')
    #plt.ylim(top=4.5,bottom=0)
    #plt.plot(tt,np.abs(newv),".",c='r')
'''
#plt.xlim(left=2.5,right=13.5)
plt.plot(tt,np.abs(xx1),".",c='b')
plt.show()
plt.savefig("./dplots/oldvisplot_"+auth+date+corr+itstr+".png",overwrite=True)
plt.clf()

'''
for t in range(nbl):
    #oldp=np.array(xx1[t])
    newp=np.array(newpt[t])
    #thispt=(tbtV==utbt[t])&(tbspw==0)
    #newv=np.array(xx1[thispt])
    #plt.plot(tt,np.abs(oldp),".",c='b')
    plt.plot(tt,np.abs(newp),".",c='r')
'''
#plt.xlim(left=2.5,right=13.5)
plt.plot(tt,np.abs(newpt),".",c='r')
plt.show()
plt.savefig("./dplots/newvisplot_"+auth+date+corr+itstr+".png",overwrite=True)
plt.clf()

#gf.dict_update(antinfo,'newres_'+corr+itstr,newpt)
gf.dict_update(antinfo1,'newvis_'+corr+itstr,newpt)
with open('data/antinfo1It'+str(it)+'.pickle', 'wb') as f:
    pickle.dump(antinfo1, f)
with open('data/newvisIt'+str(it)+'.pickle', 'wb') as f:
    pickle.dump(newpt, f)
#ast.close()
#return newpt,antinfo1
#newres,x,antinfo,x
