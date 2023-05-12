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
import math


attn=0
tavg='60s'
tavgint=60.
MJDdate=58232.00000000*3600.*24. #apr242018
#MJDdate=57850.00000000 #apr072017
nspw=4


#try:
#    import uvmultimodel
#except ModuleNotFoundError:
import sys
sys.path.insert(0,"/home/jasminwash/.local/lib/python2.7/site-packages/NordicARC-3.0.5-py2.7-linux-x86_64.egg")
#import uvmultimodel


from NordicARC import uvmultifit as uvm

def siggobble(msfile,tavgflt):
    #print('31: flag phasing scans')
    tb.open(msfile)
    scans = tb.getcol('SCAN_NUMBER')
    #field = tb.getcol('FIELD_ID')
    ints=tb.getcol('INTERVAL')
    siggy = tb.getcol('SIGMA')[0]
    tbData = tb.getcol('DATA')[0]
    tb.close()
    #cdn = field == fieldID # here 1 is the field ID of the science source
    targetScans = np.sort(np.unique(scans))
    #targetScansIdx = np.insert(np.diff(targetScans),0,1)
    #cdn = ((abs(tbData)/siggy)<3)[0]|(ints<35.)
    cdn = ints<tavgflt/2.
    #badScans=[]
    badScans=scans[cdn]
    #badScans = np.sort( np.append(targetScans[cdn], targetScans[cdn]+1) ) # the second array inside append is to select second scan from agroup
    #for t in range(len(targetScansIdx)):
    #    if t+1==len(targetScansIdx):
    #        badScans.append(targetScans[t])
    #        continue
    #    if targetScansIdx[t]>3 or (targetScansIdx[t-1]<3 and targetScansIdx[t+1]>3):
    #        badScans.append(targetScans[t])
    badscan  = str([b for b in targetScans if b in badScans])
    badscan  = badscan.strip('[]')
    flagdata(vis=msfile,
             mode = 'manual',
             scan = badscan,
             flagbackup = True)

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
    tbsig=tb.getcol('SIGMA')
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
    xxsig=np.zeros_like(tbsig[pol])
    #flagpt=((visdata['flag'][pol]==True)|(visdata['data'][pol]==0))
    flagpt=((tbFlag[pol]==True)|(tbData[pol]==0))
    #flagpt1=(visdata['data'][pol]==0)
    #xx0=np.squeeze(np.where(flagpt==True,np.nan,visdata['data'][pol]))
    #xx=np.squeeze(np.where(flagpt==True,np.nan,visdata['data'][pol]))
    xx=np.squeeze(np.where(flagpt==True,np.nan,tbData[pol]))
    xxsig=np.squeeze(np.where(flagpt==True,np.nan,tbsig[pol]))
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
    return xx,xxsig



####################################################################
def gsolve(nbl,nant,vdata,varr,varrcln,antlist,corr,itstr,refant,case,tba1,tba2,tbt,tbspw,sigarr,sigarrcln):#ainfo
    #gst0=time.time()
    plt.clf()
    theta_r_dict={}
    theta_m_dict={}
    theta_Ir_dict={}
    theta_del_dict={}
    theta_rms_dict={}

    script_L_dict={}
    script_L_nn_dict={}
    theta_r_nn_dict={}
    l_m_dict={}
    l_m_nn_dict={}
    theta_m_nn_dict={}
    l_r_dict={}
    l_del_dict={}
    l_rms_dict={}
    l_Ir_conv_dict={}
    all_bdance={}
    #all_bdance_ph={}

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
            Theta_r=np.zeros((nbl,nant-1),dtype=int)
            bad_T=np.zeros((nbl,nant-1),dtype=int)
            theta_m=np.zeros((nbl,1),dtype=float)

            #Defining Jacobian and measured amp matrix
            script_L=np.zeros((nbl,nant),dtype=int)
            bad_L=np.zeros((nbl,nant),dtype=int)
            l_m=np.zeros((nbl,1),dtype=float)
            astore=np.zeros((nbl,1),dtype=float)
            bantstore=[]
            blstore=[]


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
                            l_m[nb,0]=np.nan
                            astore[nb,0]=np.nan

                            script_L[nb,iant1]=1#1
                            script_L[nb,iant2]=1#1


                            #PHASE STUFF
                            theta_m[nb,0]=np.nan

                            if ant1==refant: Theta_r[nb,iant2-1]=-1#-1
                            if ant2==refant: Theta_r[nb,iant1]=1#1
                            if ant1!=refant and ant1>refant:
                                Theta_r[nb,iant1-1]=1#1
                                Theta_r[nb,iant2-1]=-1#-1
                            if ant1!=refant and ant2<refant:
                                Theta_r[nb,iant1]=1#1
                                Theta_r[nb,iant2]=-1#-1
                            if (ant1!=refant and (ant2>refant and ant1<refant)):
                                Theta_r[nb,iant1]=1#1
                                Theta_r[nb,iant2-1]=-1
                            nb+=1
                            continue
                        if thispt.sum()>0:
                            #potential 0 after thisbase and thistime
                            #pt=varr[thisbase][0][j]
                            pt=varr[thispt][0]
                            ptc=varrcln[thispt][0]
                            #print("PT",pt)
                            #print("PTC",ptc)
                            #pdb.set_trace()
                            ph=np.angle(pt,deg=False)
                            phc=np.angle(ptc,deg=False)
                            amp=abs(pt)
                            ampc=abs(ptc)
                            if ampc/sigarrcln[thispt][0]<5:
                                l_m[nb,0]=np.nan
                                astore[nb,0]=np.nan

                                script_L[nb,iant1]=1#1
                                script_L[nb,iant2]=1#1


                                #PHASE STUFF
                                theta_m[nb,0]=np.nan

                                if ant1==refant: Theta_r[nb,iant2-1]=-1#-1
                                if ant2==refant: Theta_r[nb,iant1]=1#1
                                if ant1!=refant and ant1>refant:
                                    Theta_r[nb,iant1-1]=1#1
                                    Theta_r[nb,iant2-1]=-1#-1
                                if ant1!=refant and ant2<refant:
                                    Theta_r[nb,iant1]=1#1
                                    Theta_r[nb,iant2]=-1#-1
                                if (ant1!=refant and (ant2>refant and ant1<refant)):
                                    Theta_r[nb,iant1]=1#1
                                    Theta_r[nb,iant2-1]=-1
                                nb+=1
                            else:
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

                                l_m[nb,0]=np.log10(amp/ampc)
                                astore[nb,0]=amp

                                script_L[nb,iant1]=1
                                script_L[nb,iant2]=1


                                #PHASE STUFF
                                theta_m[nb,0]=ph-phc

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
            #print(j)
            badance=np.full((nant,1),False,dtype=bool)
            badbl=np.full((nbl,1),False,dtype=bool)
            for ia in range(nant):
                #print(bantstore.count(ia))
                #if ia==38: 
                #    if bantstore.count(ia)>0: badance[ia,0]=True
                #else:
                #if bantstore.count(ia)==nant-1: badance[ia,0]=True
                #if ia in bantstore: badance[ia,0]=True
                if bantstore.count(ia)==(nant-1)-ia: 
                    badance[ia,0]=True
                #else: print("bonk!")
            for ib in range(nbl):
                if ib in blstore: badbl[ib,0]=True
            #print(len(blstore))
            #pdb.set_trace()    
            #badance=np.unique(bantstore)
            #print(badance)
            plt.scatter(j,len(np.where(badance==True)[0]))
            gf.dict_update(all_bdance,time+corr+itstr+spwstr,badance)


            gf.dict_update(theta_r_dict,time+corr+itstr+spwstr,Theta_r)
            gf.dict_update(theta_m_dict,time+corr+itstr+spwstr,theta_m)


            #scl.write("\n\n")

            #lme.write("\n\n")

            #ast.write(str(astore))
            #ast.write("\n\n")
            #help me mark I hope
            #pdb.set_trace()

            script_L_nn=gf.nonan(script_L)
            gf.dict_update(script_L_dict,time+corr+itstr+spwstr,script_L)
            #scl.write(str(script_L_nn))
            l_m_nn=gf.nonan(l_m)
            gf.dict_update(l_m_dict,time+corr+itstr+spwstr,l_m)
            #lme.write(str(l_m_nn))

            gf.dict_update(script_L_nn_dict,time+corr+itstr+spwstr,script_L_nn)
            gf.dict_update(l_m_nn_dict,time+corr+itstr+spwstr,l_m_nn)


            #l_r=gf.l_Ir(ll_r=script_L_nn,ll_m=l_m_nn)
            #print("wherelmnnfalse",len(np.where(l_m_nn.mask==False)[0]))
            #print("j",j)
            l_r=gf.l_Ir_mask(ll_r=script_L_nn,ll_m=l_m_nn,nant=nant,bbl=badbl,antlist=antlist,exclAnts=exclAnts)
            #l_r=gf.l_Ir_nan(ll_r=script_L,ll_m=l_m,ainfo=ainfo,t=j,nant=nant,corr=corr,itstr=itstr,bdance=badance,bbl=badbl)


            Theta_r_nn=gf.nonan(Theta_r)
            gf.dict_update(theta_r_nn_dict,time+corr+itstr+spwstr,Theta_r_nn)
            theta_m_nn=gf.nonan(theta_m)
            gf.dict_update(theta_m_nn_dict,time+corr+itstr+spwstr,theta_m_nn)


            theta_Ir=gf.th_Ir_mask(Th_r=Theta_r_nn,th_m=theta_m_nn,rfantind=rfantind,nant=nant,bbl=badbl,antlist=antlist,exclAnts=exclAnts)
            #theta_Ir=gf.th_Ir(Th_r=Theta_r_nn,th_m=theta_m_nn)



            #Residuals
            l_del=l_m_nn-np.matmul(script_L_nn,l_r)
            #print("l_del \n"+str(l_del))
            gf.dict_update(l_del_dict,time+corr+itstr+spwstr,l_del)

            #print(l_del.shape)

            #l_Ir_res=gf.l_Ir_mask(ll_r=script_L_nn,ll_m=l_del,nant=nant,bbl=badbl,antlist=antlist,exclAnts=exclAnts)

            #print("l_Ir w/ resids: \n"+str(l_Ir_res))

            theta_del=theta_m_nn-ma.dot(Theta_r_nn,theta_Ir)
            #print("theta_del \n"+str(theta_del))
            gf.dict_update(theta_del_dict,time+corr+itstr+spwstr,theta_del)

            #theta_Ir_res=gf.th_Ir_mask(Th_r=Theta_r_nn,th_m=theta_del,rfantind=rfantind,nant=nant,bbl=badbl,antlist=antlist,exclAnts=exclAnts)
            #theta_Ir_res=gf.th_Ir(Th_r=Theta_r_nn,th_m=theta_del)
            #theta_Ir_res,restd,rnktd,stdel=lstsq(a=Theta_r_nn,b=theta_del)
            #print("theta_Ir w/ resids: \n"+str(theta_Ir_res))


            #l_Ir_final=l_r+l_Ir_res
            #ast.write(str(l_Ir_final))
            l_Ir_final=l_r
            #print("l_ir_shape")
            #print(l_Ir_final.shape)
            gf.dict_update(l_r_dict,time+corr+itstr+spwstr,l_Ir_final)
            #print("final amps (log form) \n"+str(l_Ir_final))

            #with open('data/model_redshifts.pickle', 'rb') as f:
            #    redshifts = pickle.load(f)

            Ir_converted=np.zeros_like(l_Ir_final)

            #Ir_converted=np.vstack([10.0**x for x in l_Ir_final])
            for x in range(l_Ir_final.shape[0]):
                lir=l_Ir_final[x,0]
                Ir_converted[x,0]=10.0**lir

            #print("l_conv shape")
            #print(Ir_converted.shape)

            gf.dict_update(l_Ir_conv_dict,time+corr+itstr+spwstr,Ir_converted)

            #theta_Ir_final=theta_Ir+theta_Ir_res
            theta_Ir_final=theta_Ir
            gf.dict_update(theta_Ir_dict,time+corr+itstr+spwstr,theta_Ir_final)


          #plt.show()
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


    plt.savefig('validants.png',overwrite=True)
    plt.show()
    plt.clf()
    #gst1=time.time()
    #print("Time to run is",(gst1-gst0)/60,"minutes",(gst1-gst0)/3600,"hours")
    return 0



####################################################################

def lls(pol,corr,datams1,target,case,auth,refmeth,date,it,dvis,rfant,datams_cm):
    itstr=str(it)
    nspw=4

    antinfo={}

    print(datams1)
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
    xx,xxsig=antfunc(dfile=datams1,it=it,pol=pol,date=date,case='ext')
    #antinfo_clean,
    xx_clean,xxsig_clean=antfunc(dfile=datams_cm,it=it,pol=pol,date=date,case='cln')
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

    gsolve(nbl=nbl,nant=nant,vdata=tbData,varr=xx,varrcln=xx_clean,antlist=antlist,corr=corr,itstr=itstr,refant=refant,case='all',tba1=tba1,tba2=tba2,tbt=tbt,tbspw=tbspw,sigarr=xxsig,sigarrcln=xxsig_clean)#ainfo=antinfo
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

    with open('data/theta_Ir_dict_all.pickle', 'rb') as f:
        tir_all = pickle.load(f)
    with open('data/l_Ir_conv_dict_all.pickle', 'rb') as f:
        lir_all = pickle.load(f)

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
    fkeys=tir_all.keys()
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
                amp_all1[k,0,t_step,s]=lir_all[str(t_step)+corr+itstr+spwstr][k]
                if k==refant: 
                    phase_all1[k,0,t_step,s]=0
                elif k<refant:
                    phase_all1[k,0,t_step,s]=tir_all[str(t_step)+corr+itstr+spwstr][k]
                elif k>refant:
                    phase_all1[k,0,t_step,s]=tir_all[str(t_step)+corr+itstr+spwstr][k-1]

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

    xx1,xx1sig=antfunc(dfile=dvis,it=it,pol=pol,date=date,case='raw')


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
    gf.dict_update(antinfo,'newvis_'+corr+itstr,newpt)
    with open('data/antinfoIt'+str(it)+'.pickle', 'wb') as f:
        pickle.dump(antinfo, f)
    with open('data/newvisIt'+str(it)+'.pickle', 'wb') as f:
        pickle.dump(newpt, f)

    #ast.close()
    return newpt,antinfo
    #newres,x,antinfo,x



##################################################################
#Setting default value for TCLEAN stopcode
stop=0

#To time the code
start=time.time()

#this is for Venki data
target='sgr_apr24'

#for image names
case='allant'
auth='Jasmin'
#refmeth='ref35'
refmeth='ref6'
date='05052023_'+str(attn)
refant=6

#datams2=target+'_flagcor.ms'

#flagfile=datams2[:-3]+'_noflags.ms'
#os.system('rm -rf '+flagfile)
#split(vis=datams2,outputvis=flagfile,keepflags=False,datacolumn='data')

#Initial data file name and name for channel-averaged file
#datams1=target+'_flagcor_noflags.ms'
#datams1=flagfile
#datams_noavg='sgr_apr07_flag.ms'
#datams0='sgr_apr24_200scanXXYY.ms'
datams_noavg='sgr_apr24_final_XXYY.ms'
#datams_noavg='sgr_apr24_mid.ms'
#datams0='sgr_apr24_final_XXYY_rdata_avg60s.ms'

rdavg=datams_noavg[:-3]+'_rdata_avg'+tavg+'.ms'
os.system('rm -rf '+rdavg+' '+rdavg+'.flagversions')
split(vis=datams_noavg,
      outputvis=rdavg,
      timebin=tavg,
      combine='state,scan',
      datacolumn='data',keepflags=False)
clearcal(vis=rdavg,spw='0,1,2,3',addmodel=True)

datams_nf=rdavg

siggobble(datams_nf,tavgint)

rdflag=datams_noavg[:-3]+'_rflag_avg'+tavg+'.ms'
os.system('rm -rf '+rdflag+' '+rdflag+'.flagversions')
split(vis=rdavg,
      outputvis=rdflag,
      #timebin=tavg,
      #combine='state,scan',
      datacolumn='corrected',keepflags=False)
clearcal(vis=rdflag,spw='0,1,2,3',addmodel=True)

datams0=rdflag

#Getting file name
dmsprefix=datams0[:-3]

#Setting name of file for time-averaged split file
resavg=dmsprefix+'_avgresres.ms'

#Getting file name
#dmsprefix=dmsavg[:-3]
datams1=resavg


#Splitting and channel averaging (if desired - no we need for resres)
os.system('rm -rf '+resavg+' '+resavg+'.flagversions')
split(vis=datams0,
      outputvis=resavg,
      #timebin=tavg,
      #combine='state,scan',
      datacolumn='data',keepflags=False)

#mstransform(vis=datams3,outputvis=dmsavg,keepflags=False,timeaverage=True,timebin='240s',timespan='scan,state',datacolumn='data')

rawdata=datams0
clearcal(vis=rawdata,spw='0,1,2,3',addmodel=True)

rawsplit=dmsprefix+'_rawsplit.ms'
os.system('rm -rf '+rawsplit+' '+rawsplit+'.flagversions')
split(vis=datams0,
      outputvis=rawsplit,
      #timebin='0s',
      #combine='scan,state',
      datacolumn='data',keepflags=False)

'''
rsavg=rawsplit[:-3]+'_avg'+tavg+'.ms'
os.system('rm -rf '+rsavg+' '+rsavg+'.flagversions')
split(vis=datams0,
      outputvis=rsavg,
      timebin=tavg,
      combine='scan,state',
      datacolumn='data',keepflags=False)
'''

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


clearcal(vis=datams0,spw='0,1,2,3',addmodel=True)

tbdata={}
tb.open(datams0)
tbdata['data']=tb.getcol('DATA')
tbdata['data']=np.squeeze(tbdata['data'])
tbdata['flag']=tb.getcol('FLAG')
tbdata['flag']=np.squeeze(tbdata['flag'])
tbdata['time']=tb.getcol('TIME')
tbdata['antenna1']=tb.getcol('ANTENNA1')
tbdata['antenna2']=tb.getcol('ANTENNA2')
tbdata['sigma']=tb.getcol('SIGMA')
tb.close()

ms.open(datams0,nomodify=True)
visdata = ms.getdata(['axis_info'],ifraxis=True)
#visdata['data'] = np.squeeze(visdata['data'])
ms.close()

#tb.open(datams0)
#a1=tb.getcol('ANTENNA1')
#a2=tb.getcol('ANTENNA2')
#tb.close()
#allants=np.concatenate((a1,a2))
#antlist=np.unique(allants)
allants=np.concatenate((tbdata['antenna1'],tbdata['antenna2']))
antlist=np.unique(allants)

print(antlist)

#Number of polarizations
npol=tbdata['data'].shape[0]
polnames=visdata['axis_info']['corr_axis']
#npol=len(polnames)
#print("NPOL",npol)
#corr=polnames[pol]


#Calculating number of antennas and baselines
nant=int(len(antlist))
nbl=int(nant*(nant-1)/2)
#ntimes=len(visdata['axis_info']['time_axis']['MJDseconds'])
ntimes=len(np.unique(tbdata['time']))

#tt=np.array(visdata['axis_info']['time_axis']['MJDseconds'],dtype=float)
tt=np.array(np.unique(tbdata['time']))

'''
fstore=np.zeros_like(visdata['data'])

for pol in range(npol):
    flnan=np.zeros_like(visdata['flag'][pol])
    fl=(visdata['flag'][pol]==True)
    flnan=np.squeeze(np.where(fl==True,np.nan,visdata['data'][pol]))
    print("flnan shape",flnan.shape)
    fstore[pol,:,:]=flnan
    #visdata['data'][pol]=flnan

ms.open(datams1,nomodify=False)
nv1=ms.getdata(['corrected_data'],ifraxis=True)
nv1['corrected_data'][:,0,:,:]=fstore
ms.putdata(nv1)
ms.close()


#ms.open(datams1,nomodify=True)
#testdata = ms.getdata(['data','flag','data_desc_id','sigma','axis_info'],ifraxis=True)
#testdata['data'] = np.squeeze(testdata['data'])
#ms.close()

#xx=testdata['data'][0]
#for y in range(nbl):
#    pp1=np.array(xx[y])
#    #plt.scatter(tt,pp)
#    plt.scatter(tt,np.abs(pp1))
#plt.show()
#plt.savefig('ftp.png')
#plt.clf()

for i in range(npol):  
    ms.open(datams1,nomodify=True)
    testdata = ms.getdata(['corrected_data','flag','data_desc_id','sigma','axis_info'],ifraxis=True)
    testdata['corrected_data'] = np.squeeze(testdata['corrected_data'])
    ms.close()

    xx=testdata['corrected_data'][i]
    for y in range(nbl):
        pp1=np.array(xx[y])
        #plt.scatter(tt,pp)
        plt.scatter(tt,np.abs(pp1))
    plt.show()
    #plt.savefig('datatestplot_'+auth+refmeth+'_'+date+corr+str(it)+'.png',overwrite=True)
    plt.savefig('flagtestplot_'+str(i)+'.png',overwrite=True)
    plt.clf()
'''

print(nant)
print(nbl)

ELO=np.unique(tbdata['time'])
'''
#Creating a tf dictionary for each ant and whether time is good or bad there
antinfo=dict()
antinfo['timestamp']=visdata['axis_info']['time_axis']['MJDseconds']

allant=np.zeros((nant,ntimes),dtype=bool)
alltim=np.zeros((ntimes,nant),dtype=bool)

########Inserted from looptest#######
#Creating a tf dictionary for each ant and whether time is good or bad there
antinfo=dict()
antinfo['timestamp']=visdata['axis_info']['time_axis']['MJDseconds']
tt=np.array(visdata['axis_info']['time_axis']['MJDseconds'],dtype=float)
#####################################
'''

#print(len(visdata['axis_info']['time_axis']['MJDseconds']))
#print(len(np.unique(visdata['axis_info']['time_axis']['MJDseconds'])))



#print("num intervals: %f" %((np.max(ELO)-np.min(ELO))/240))
#print(.402600*240) 
#sys.exit()


#tsize,goodtimes=t_length(visdata=visdata,ELO=ELO,nbl=nbl)
#print("TSIZE")
#print(goodtimes)

'''
#Splitting and channel averaging (if desired)
os.system('rm -rf '+dmsavg+' '+dmsavg+'.flagversions')
split(vis=datams1,
      outputvis=dmsavg,
      datacolumn='data',keepflags=False)

rawdata=dmsavg
clearcal(vis=rawdata,spw='0,1,2,3',addmodel=True)

rawsplit=dmsprefix+'_rawsplit.ms'
os.system('rm -rf '+rawsplit+' '+rawsplit+'.flagversions')
split(vis=dmsavg,
      outputvis=rawsplit,
      datacolumn='data',keepflags=False)
'''

#Setting initial data file for analysis
datams1=rawsplit

#Getting scan numbers from file
#tb.open(datams1,nomodify=True)
tb.open(datams0,nomodify=True)
scan = tb.getcol('SCAN_NUMBER')

#Creating array to store residuals
#visdata0 = tb.getcol('DATA')
#resres=np.zeros_like(visdata0)
tb.close()

#tb.open(resavg,nomodify=True)
#visdata0=tb.getcol('DATA')
#resres=np.zeros_like(visdata0)
#scanavg=tb.getcol('SCAN_NUMBER')
#tb.close()


#Print scan numbers for edification
scan = np.unique(scan)
print(scan)
#scanavg = np.unique(scanavg)



#Calculating RMS of data for CLEAN threshold
#ms.open(datams1,nomodify=True)
#ms.open(datams0,nomodify=True)
#sigs=ms.getdata(['sigma'])
#fsigs=np.asarray(sigs['sigma']).flatten()
fsigs=tbdata['sigma']
print("Fsig shape:",fsigs.shape)
invsigs=1./(fsigs*fsigs)
sumsigs=np.sum(invsigs)
rmst=2.*np.sqrt(1./sumsigs)
print("RMS: %e" %(rmst))
#ms.close()



#Creating and initializing data file to store UVMULTIFIT pt source flux
g=open("uvmflux_"+target+".txt","w")
g.write("it scan ext mod moderr\n")

q=open("nv_"+target+".txt","w")
q.write("")


'''
#######from looptest########
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
    plt.savefig("./lpgraphs/ant_tf_ant%i_%s_it%i.png"%(ant,corr,it),overwrite=True)
    plt.clf()
    #antinfo['goodant']=allant
############################
'''

#Start with first iteration
it=0

#refant=gf.refantfinder(antlist=antlist,goodant=antinfo['goodant'])
#iref=np.where(antlist==refant)[0][0]

#refant=35 #DV20 index
#refant=5
iref=np.where(antlist==refant)[0][0]


#The wily (formerly) infinite while loop
#while(1):
while it<=0:
    plt.clf()
    #rdavg=dmsprefix+'_rdata_avg'+tavg+'.ms'
    #os.system('rm -rf '+rdavg+' '+rdavg+'.flagversions')
    #split(vis=rawdata,
    #      outputvis=rdavg,
    #      timebin=tavg,
    #      combine='state,scan',
    #      datacolumn='data',keepflags=False)
    #clearcal(vis=rdavg,spw='0,1,2,3',addmodel=True)

    #datams1=rdavg

    q.write("ITERATION "+str(it)+"\n")

    #ms.open(rdavg)
    #adata=ms.getdata(['data'],ifraxis=True)
    #adata['data']=np.squeeze(adata['data'])
    #ms.close()
    #tb.open(rdavg)
    #adata=tb.getcol('DATA')
    #tb.close()

    #Declaring empty arrays for mod UVM flux and error
    muvmf=np.array([])
    muvmferr=np.array([])
    ruvmf=np.array([])

    #Creating file names for residual and model data from UVM
    datams_ext=dmsprefix+'_ext_it'+'{0:02d}'.format(it)+'.ms'
    datams_mod=dmsprefix+'_mod_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files
    os.system('rm -rf '+datams_ext)
    split(vis=datams1,outputvis=datams_ext,datacolumn='data',keepflags=False)

    #cvis.append(str(datams_mod))
    os.system('rm -rf '+datams_mod)
    split(vis=datams1,outputvis=datams_mod,datacolumn='data',keepflags=False)

    #Clearing/creating corrected and model columns
    clearcal(vis=datams_ext,spw='0,1,2,3',addmodel=True)
    clearcal(vis=datams_mod,spw='0,1,2,3',addmodel=True) 

    #Cycle through each scan
    for s in scan:
    #for s in scanavg:
        #print(s)
        if s==359:continue
        #UVM for extended structure leftover from fit
        sc=int(s)
        myuvfit_ext1 = uvm.uvmultifit(vis=datams_ext,
                    spw='0,1,2,3',
                    #scans=[[sc]],
                    scans=[sc],
                    model=['delta'],
                    var=['0,0,p[0]'],
                    p_ini=[0.],
                    OneFitPerChannel=False,
                    column='data',
                    write='residuals') # will write residuals in the 'corrected' column, but only for the subset of fitted channels !
        # to run on the continuum collapsed cube
        ruvmf=np.append(ruvmf,myuvfit_ext1.result['Parameters'][0])
        plt.scatter(s,myuvfit_ext1.result['Parameters'][0])

        #UVM for model visibilities
        myuvfit_mod = uvm.uvmultifit(vis=datams_mod,
                    spw='0,1,2,3',
                    #scans=[[sc]], 
                    scans=[sc],
                    model=['delta'],
                    var=['0,0,p[0]'], 
                    p_ini=[0.], 
                    OneFitPerChannel=False,
                    column='data',
                    write='model') # will write best-fit model in the 'model' column, but only for the subset of fitted channels !
        #tb.open(datams_mod)
        #md=tb.getcol('MODEL_DATA')
        #tb.close()
        #Appending results of UVM to list
        muvmf=np.append(muvmf,myuvfit_mod.result['Parameters'][0])
        muvmferr=np.append(muvmferr,myuvfit_mod.result['Uncertainties'][0])

    plt.savefig(target+'uvmresults'+date+'.png',overwrite=True)
    plt.show()
    plt.clf()
    #Pulling residuals of UVM pt source fitting for analysis
    tb.open(datams_ext)
    ext1=tb.getcol('CORRECTED_DATA')
    tb.close()

    #Name of file to split these residuals into
    datams_ext2=dmsprefix+'_ext2_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files
    #The datacolumn here means nothing, it will be replaced
    os.system('rm -rf '+datams_ext2)
    split(vis=datams_ext,outputvis=datams_ext2,datacolumn='corrected',keepflags=False)

    #Making the residuals of the first UVMULTIFIT the data column of the new file
    #4/26/23 - just changed split again???
    #tb.open(datams_ext2,nomodify=False)
    #tb.putcol('DATA', ext1)
    #tb.close()

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
    #tb.open(datams_ext2)
    #ext2=tb.getcol('CORRECTED_DATA')
    #tb.close()

    #ext21=np.zeros_like(tbdata['data'])

    #ms.open(datams_ext2)
    #edata=ms.getdata(['corrected_data'],ifraxis=True)
    #ext21[:,:,:]=np.squeeze(edata['corrected_data'])
    #ms.close()
    

    #Splitting and creating the files
    #The datacolumn here means nothing, it will be replaced
    #4/26/23: changed data to corrected??
    os.system('rm -rf '+datams_ext3+' '+datams_ext3+'.flagversions')
    split(vis=datams_ext2,outputvis=datams_ext3,datacolumn='corrected',keepflags=False)
    clearcal(vis=datams_ext3,spw='0,1,2,3',addmodel=True)

    #Making the residuals of the first UVMULTIFIT the data column of the new file
    #tb.open(datams_ext3,nomodify=False)
    #tb.putcol('DATA', ext2)
    #tb.close()

    #ms.open(datams_ext3,nomodify=False)
    #nv1=ms.getdata(['data'],ifraxis=True)
    #nv1['data'][:,0,:,:]=ext21
    #ms.putdata(nv1)
    #ms.close()

    #dme3avg=datams_ext3[:-3]+'avg'+tavg+'_it'+'{0:02d}'.format(it)+'.ms'
    #Splitting and channel averaging (if desired)
    #os.system('rm -rf '+dme3avg+' '+dme3avg+'.flagversions')
    #split(vis=datams_ext3,
    #      outputvis=dme3avg,
    #      #timebin=tavg,
    #      combine='state,scan',
    #      datacolumn='data',keepflags=False)

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
    imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_cleanext'
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

    #mscm=dmsprefix+'_cmodel_avg'+tavg+'_it'+'{0:02d}'.format(it)+'.ms'
    #os.system('rm -rf '+mscm+' '+mscm+'.flagversions')
    #split(vis=dme3avg,
    #      outputvis=mscm,
    #      #timebin=tavg,
    #      combine='state,scan',
    #      datacolumn='model',keepflags=False)

    datams_cln=dmsprefix+'_cln_it'+'{0:02d}'.format(it)+'.ms'
    #Splitting and channel averaging (if desired)
    os.system('rm -rf '+datams_cln+' '+datams_cln+'.flagversions')
    split(vis=datams_ext3,
          outputvis=datams_cln,
          #timebin=tavg,
          #combine='state,scan',
          datacolumn='model',keepflags=False)
        
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

    #ms.open(rsavg,nomodify=True)
    #avgdata=ms.getdata(['data'],ifraxis=True)
    #ms.close()

    newvis=np.zeros_like(tbdata['data'])
    #Zeros like data array for averaged data
    #newvis1=np.zeros_like(adata)
    #newvis1=np.zeros_like(avgdata['data'])
    #antinfo=dict()
    antinfo1=dict()
    plt.clf()
    for p in range(npol):
        plt.clf()
        #newvis[p],newvis1[p],newantinfo,newantinfo1
        #newvis1[p],newantinfo1=lls(pol=p,corr=polnames[p],datams1=datams_ext3,target=target,case=case,auth=auth,refmeth=refmeth,date=date,it=it,dvis=rawdata,rfant=refant,datams_cm=mscm)
        newvis[p],newantinfo1=lls(pol=p,corr=polnames[p],datams1=datams_ext3,target=target,case=case,auth=auth,refmeth=refmeth,date=date,it=it,dvis=rawdata,rfant=refant,datams_cm=datams_cln)
        plt.clf()
        #newvis1[p],newantinfo1=lls(pol=p,corr=polnames[p],datams1=datams_ext3,target=target,case=case,auth=auth,refmeth=refmeth,date=date,it=it,dvis=rawdata)
        #antinfo.update(newantinfo)
        antinfo1.update(newantinfo1)
    
    '''
    ms.open(datams_ext3,nomodify=False)
    nv=ms.getdata(['corrected_data'],ifraxis=True)
    nv['corrected_data'][:,0,:,:]=newvis
    print(nv['corrected_data'][0])
    ftnv=nv['corrected_data'][0]
    e.write(str(ftnv)+"\n")
    ms.putdata(nv)
    ms.close()
    

    ms.open(datams_ext3,nomodify=True)
    nnv=ms.getdata(['corrected_data'],ifraxis=True)
    ftnnv=nnv['corrected_data'][0]
    print(ftnnv)
    e.write(str(ftnnv)+"\n")
    ms.close()
    '''

    tb.open(rawdata,nomodify=False)
    #ms.open(rawdata,nomodify=False)
    nv1=tb.getcol('CORRECTED_DATA')
    nv1[:,0,:]=newvis1
    tb.putcol('CORRECTED_DATA',nv1)
    #ms.putdata(nv1)
    #ftnv1=nv1['corrected_data'][0]
    #q.write(str(ftnv1)+"\n")
    #print(ftnv1)
    #ms.close()
    tb.close()

    #ms.open(rawdata,nomodify=True)
    #ms.open(rawdata,nomodify=True)
    #nnv1=ms.getdata(['corrected_data'],ifraxis=True)
    #ftnnv1=nnv1['corrected_data'][0]
    #q.write(str(ftnnv1)+"\n")
    #print(ftnnv1)
    #ms.close()



    #Dirty image of gain-calibrated data
    imname3=dmsprefix+'_it'+'{0:02d}'.format(it)+'_gainappl'
    os.system('rm -rf '+imname3+'.*')
    r2c1=tclean(#vis=rdavg
        vis=rawdata,
           imagename=imname3,
           spw='0,1,2,3',
           datacolumn='corrected',
           niter=0,
           threshold=rmst,
           pblimit=-1,
           imsize=[300,300],
           interactive=0,
           cell='0.2arcsec')

    '''
    #Dirty image of x gain-calibrated residuals x CLEAN model
    imname4=dmsprefix+'_it'+'{0:02d}'.format(it)+'_gainappl_ext'
    os.system('rm -rf '+imname4+'.*')
    r2c2=tclean(vis=datams_ext3,
           imagename=imname4,
           spw='0,1,2,3',
           datacolumn='model',
           niter=0,
           threshold=rmst,
           pblimit=-1,
           imsize=[300,300],
           interactive=0,
           cell='0.2arcsec')
    '''

    #Split again so gain-calibrated data and residuals are the data
    dms_selfcal=dmsprefix+"_selfcal_it"+str(it)+".ms"
    dms_selfcalf=dms_selfcal+'.flagversions'
    dms_selfcal_copy=dmsprefix+"_selfcal_copy_it"+str(it)+".ms"
    dms_selfcalf_copy=dms_selfcal_copy+'.flagversions'
    dms_selfcal_ext=dmsprefix+"_selfcal_ext_it"+str(it)+".ms"
    dms_selfcalf_ext=dms_selfcal_ext+'.flagversions'

    os.system('rm -rf '+dms_selfcal+' '+dms_selfcalf)
    #split(vis=rawdata,outputvis=dms_selfcal,datacolumn='corrected',keepflags=False)
    split(vis=rawdata,outputvis=dms_selfcal,datacolumn='corrected',keepflags=False)


    os.system('rm -rf '+dms_selfcal_copy+' '+dms_selfcalf_copy)
    #split(vis=rawdata,outputvis=dms_selfcal_copy,datacolumn='corrected',keepflags=False)
    split(vis=rawdata,outputvis=dms_selfcal_copy,datacolumn='corrected',keepflags=False)

    #os.system('rm -rf '+dms_selfcal_ext+' '+dms_selfcalf_ext)
    #split(vis=datams_ext3,outputvis=dms_selfcal_ext,datacolumn='corrected',keepflags=False)


    #Seeing why TCLEAN ended
    scode=r2c['stopcode']

    #Pulling (not) gain-calibrated residuals to add to previous iterations
    #tb.open(dms_selfcal_ext)
    #tb.open(datams_ext3)
    tb.open(dme3avg)
    resids=tb.getcol('DATA')
    tb.close()
    resres+=resids

    #Pulling gain-calibrated full data
    tb.open(dms_selfcal)
    print("Getting",it)
    scdata=tb.getcol('DATA')
    tb.close()

    #Subtracting residuals from full data to hopefully get point source
    #dminres=scdata-resids

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
    clearcal(vis=rawdata,spw='0,1,2,3',addmodel=True)
    #clearcal(vis=rdavg,spw='0,1,2,3',addmodel=True)

    print("MAX:",madmax)
    if stop==2: break
    if scode!=1: stop=2

    it+=1

#Print number of iterations performed
print("Performed %i iterations" %(it))


dms_selfcal1=dmsprefix+"_selfcal_copy_it"+str(it-1)+".ms"
#dms_selfcal1_ext=dmsprefix+"_selfcal_ext_it"+str(it-1)+".ms"


tb.open(dme3avg)
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
split(vis=dms_selfcal1,outputvis=datams_ext3_sc,datacolumn='data',keepflags=False)

tb.open(datams_ext3_sc,nomodify=False)
tb.putcol('DATA',dminres)
tb.close()


muvmf=np.array([])
muvmferr=np.array([])

#Creating file names for residual and model data from UVM
datams_mod=dmsprefix+'_mod_final.ms'

#cvis.append(str(datams_mod))
os.system('rm -rf '+datams_mod)
split(vis=datams_ext3_sc,outputvis=datams_mod,datacolumn='data',keepflags=False)

#Clearing/creating corrected and model columns
clearcal(vis=datams_mod,spw='0,1,2,3',addmodel=True) 

#Cycle through each scan
for s in scanavg:
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
plt.savefig('sgra_selfcal_scanbyscan_'+target+'_final_'+date+'.png')

#Creating file for extended structure around Sgr A*
datams_fin=dmsprefix+'_extsum.ms'
os.system('rm -rf '+datams_fin)
split(vis=dme3avg,outputvis=datams_fin,datacolumn='data')

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
    plt.errorbar(x=scan[j]+(0.1*it[j]),y=modf[j],yerr=modferr[j],c=cs[it[j]],marker=marks[it[j]],ls='none')

plt.title('Sgr A time series')
plt.xlabel('scan number')
plt.ylabel('flux [Jy]')
plt.show()
plt.savefig('sgratime_scan_'+target+'_final_pswitch_'+date+'.png')


#How long this struggle bus took to run
end=time.time()
print("Time to run is",(end-start)/60,"minutes",(end-start)/3600,"hours")

