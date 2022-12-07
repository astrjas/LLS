import numpy as np
import time
import matplotlib.pyplot as plt
import os.path
from os import path
import sys
import numpy as np
path_to_gfunc='.'
sys.path.insert(0, path_to_gfunc)
import gfunc_c5 as gf
from scipy.linalg import lstsq
import numpy.ma as ma

attn='1'
tant=10
refant=8


#try:
#    import uvmultimodel
#except ModuleNotFoundError:
import sys
sys.path.insert(0,"/home/jasminwash/.local/lib/python2.7/site-packages/NordicARC-3.0.5-py2.7-linux-x86_64.egg")
#import uvmultimodel


from NordicARC import uvmultifit as uvm

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

    #lp=visdata['data'][pol]<0.01

    tt=np.array(visdata['axis_info']['time_axis']['MJDseconds'],dtype=float)



    #cycle through each ant
    #a=0

    #lp1=np.zeros_like(lp)
    #lp1=np.where(lp==True,np.nan,visdata['data'][pol])


    #Creating a tf dictionary for each ant and whether time is good or bad there
    antinfo=dict()
    antinfo['timestamp']=visdata['axis_info']['time_axis']['MJDseconds']

    allant=np.zeros((nant,ntimes),dtype=bool)
    alltim=np.zeros((ntimes,nant),dtype=bool)



    #TENTATIVE BEGINNING OF WHAT WILL BE IN MASSIVE LOOP
    #XX correlation
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
    plt.savefig('datatestplot_'+dfile[:-3]+str(it)+date+'.png',overwrite=True)
    plt.clf()

    #tf matrix for where ant is bad (True)
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
        plt.savefig("./tfgraphs/ant_tf_ant%i_%s_it%i.png"%(ant,corr,it),overwrite=True)
        plt.clf()
        a+=1

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

def lls(pol,corr,datams1,target,case,auth,refmeth,date,it,dvis,rfant):
    #datams1=flagfile
    #MS=casatools.ms()
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
   
    scl=open("scriptL_"+target+date+corr+itstr+".txt","w")
    lme=open("lm_"+target+date+corr+itstr+".txt","w")
    ast=open("astore_"+target+date+corr+itstr+".txt","w")


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
        
        #print("script_L")
        #print(script_L)
        #scl.write(str(script_L))
        scl.write("\n\n")

        #print("l_m")
        #print(l_m)
        #lme.write(str(l_m))
        lme.write("\n\n")

        #ast.write(str(astore))
        ast.write("\n\n")


        script_L_nn=gf.nonan(script_L)
        scl.write(str(script_L_nn))
        l_m_nn=gf.nonan(l_m)
        lme.write(str(l_m_nn))
        #scl.w

        l_r=gf.l_Ir(ll_r=script_L_nn,ll_m=l_m_nn)
        #l_r=gf.l_Ir(ll_r=script_L,ll_m=l_m)
        #l_r,resl,rnkl,sl=lstsq(a=script_L_nn,b=l_m_nn)

        Theta_r_nn=gf.nonan(Theta_r)
        theta_m_nn=gf.nonan(theta_m)

        #if len(theta_m)==1: theta_Ir=theta_m/Theta_r
        theta_Ir=gf.th_Ir(Th_r=Theta_r_nn,th_m=theta_m_nn)
        #theta_Ir,rest,rnkt,st=lstsq(a=Theta_r_nn,b=theta_m_nn)
        #print("theta_Ir \n"+str(theta_Ir))


        #Residuals
        l_del=l_m_nn-np.matmul(script_L_nn,l_r)
        #print("l_del \n"+str(l_del))
        gf.dict_update(l_del_dict,time+corr+itstr,l_del)
        #print(l_del.shape)

        l_Ir_res=gf.l_Ir(ll_r=script_L_nn,ll_m=l_del)
        #l_Ir_res,resld,rnkld,stld=lstsq(a=script_L_nn,b=l_del)

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
        '''
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
        '''
        gf.dict_update(l_Ir_conv_dict,time+corr+itstr,Ir_converted)

        theta_Ir_final=theta_Ir+theta_Ir_res
        gf.dict_update(theta_Ir_dict,time+corr+itstr,theta_Ir_final)

    scl.close()
    lme.close()

    #print("l_r keys")
    #print(l_r_dict.keys())
    #print("l_del")
    #print(l_del_dict['3537XX0'].shape)
    #print("l_Ir_conv")
    #print(l_Ir_conv_dict['3537XX0'].shape)
    #print(l_r_dict['4XX0'].shape)

    ELO_diff=np.max(ELO)-np.min(ELO)
    nint=int(ELO_diff/240)
    ELO_range=np.linspace(np.min(ELO),np.max(ELO),nint)

    #print("fkeys")
    fkeys=theta_Ir_dict.keys()
    tkeys=antinfo['timestamp']
    #print(fkeys)

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
            amp_arr[k,0,t_step]=l_Ir_conv_dict[str(t_step)+corr+itstr][k]
            if k==refant: 
                phase_arr[k,0,t_step]=0
            elif k<refant:
                phase_arr[k,0,t_step]=theta_Ir_dict[str(t_step)+corr+itstr][k]
            elif k>refant:
                phase_arr[k,0,t_step]=theta_Ir_dict[str(t_step)+corr+itstr][k-1]

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
    plt.clf()
    '''

    for k in range(nant):
        for t_step in range(len(ELO_range)-1):
            tlist=[]
            amplist=[]
            phlist=[]
            for x in antinfo['timestamp']:
                if x>ELO_range[t_step] and x<=ELO_range[t_step+1]: tlist.append(x)
                ind=np.where(antinfo['timestamp']==x)
                amplist.append(antinfo['amp_'+corr+itstr][k][ind])
                phlist.append(antinfo['phase_'+corr+itstr][k][ind])
                #print("ind",ind)
            #print("tlist",tlist)
            nda=len(amplist)
            plt.scatter(ELO_range[t_step],nda)
            avga[k,0,t_step]=np.nanmean(amplist)
            avgp[k,0,t_step]=np.nanmean(phlist)

            if np.isnan(np.mean(amplist))==True or np.isnan(np.mean(phlist))==True:
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
        plt.savefig('./dplots/ndtstep_ant'+str(k)+'_'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'.png')
        plt.clf()


    antinfo['avg_phase_'+corr+itstr]=np.squeeze(avgp)
    antinfo['avg_amp_'+corr+itstr]=np.squeeze(avga)    

    for k in range(nant):
        for t_step in range(len(ELO_range)-1):
            #print("ELO_range")
            #print(ELO_range[t_step])
            #print(ELO_range[t_step+1])
            nda=len([antinfo['amp_'+corr+itstr][k] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])
            plt.scatter(ELO_range[t_step],nda)
            ndb=len([l_Ir_conv_dict[str(t_step)+corr+itstr][k] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])
            plt.scatter(ELO_range[t_step],nda)

        plt.savefig('./dplots/ndtstep_ant'+str(k)+'_'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'.png')
        plt.clf()


    #plt.clf()
    for a in range(nant):
        plt.scatter(ELO_range[:-1],avga[a,0,:],color='black',marker='x')
        plt.scatter(ELO_range[:-1],avgp[a,0,:],color='red')
    plt.legend()
    plt.show()

    plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_ap_'+auth+'_'+case+corr+itstr+'.png')

    plt.clf()

    for a in range(nant):
        plt.scatter(ELO_range[:-1],avga[a,0,:],color='black',marker='x')
        plt.savefig('./dplots/antgains/amp/ant'+str(a)+'_gampavg_'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'.png')
        plt.ylim(bottom=0.0,top=3.5)
        plt.clf()
    plt.clf()


    for a in range(nant):
        plt.scatter(ELO_range[:-1],avga[a,0,:],color='black',marker='x')
    plt.legend()
    plt.show()

    plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_a_'+auth+'_'+case+corr+itstr+'.png')


    plt.clf()

    for a in range(nant):
        plt.scatter(ELO_range[:-1],avgp[a,0,:],color='red')
    plt.legend()
    plt.show()

    plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_p_'+auth+'_'+case+corr+itstr+'.png')

    plt.clf()

    for a in range(nant):
        plt.scatter(antinfo['timestamp'],antinfo['amp_'+corr+itstr][a],color='black',marker='x')
        plt.savefig('./dplots/antgains/amp/ant'+str(a)+'_gamp'+date+'_'+refmeth+'_'+auth+'_'+case+corr+itstr+'.png')
        plt.ylim(bottom=0.0,top=3.5)
        plt.clf()
    plt.clf()

    for a in range(nant):
        plt.scatter(antinfo['timestamp'],antinfo['amp_'+corr+itstr][a],color='black',marker='x')
    plt.legend()
    plt.show()

    plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_ampdi_'+auth+'_'+case+corr+itstr+'.png')

    plt.clf()


    for a in range(nant):
        plt.scatter(antinfo['timestamp'],antinfo['phase_'+corr+itstr][a],color='red')
    plt.legend()
    plt.show()

    plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_phasedi_'+auth+'_'+case+corr+itstr+'.png')

    plt.clf()


    print("Intervals with pts:",nonan)
    print("Total intervals:",len(ELO_range)-1)
    print("Bad baseline/antenna cases",nbad)
    print("Empty bins",nemp)
    print("Total points:",len(antinfo['timestamp']))


    newpt=np.zeros_like(xx)
    newres=np.zeros_like(xx)
    tb1=0

    #Opening data and pulling necessary info
    ms.open(dvis,nomodify=True)
    ms.selectinit(reset=True)
    visdata1 = ms.getdata(['antenna1','antenna2','data','axis_info','flag'],ifraxis=True)

    #Squeeze data then close ms

    visdata1['data'] = np.squeeze(visdata1['data'])
    ms.close()

    antinfo1,xx1=antfunc(dfile=dvis,it=it,pol=pol,date=date)

    for t_step in range(len(ELO_range)-1):
        tchunk=[antinfo['timestamp'][x] for x in range(len(tkeys)) if (tkeys[x]>=ELO_range[t_step] and tkeys[x]<=ELO_range[t_step+1])]
        for i in range(len(tchunk)):
            thistime=(visdata1['axis_info']['time_axis']['MJDseconds']==tchunk[i])
            time=i

            gt=antinfo['goodt_'+corr+itstr][i]

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
                            pt=xx[thisbase][0][tb1]
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
                            newres[nb1,tb1]=pt/(ga1*ga2*np.exp(1.0j*(gp2-gp1)))
                            newpt[nb1,tb1]=pt1/(ga1*ga2*np.exp(1.0j*(gp2-gp1)))
                            nb1+=1
            tb1+=1
    print(tb1)

    for t in range(nb1):
        #oldp=np.array(xx[t])
        newp=np.array(newres[t])
        #newv=np.array(newvis[t])
        #plt.plot(tt,np.abs(oldp),".",c='b')
        plt.plot(tt,np.abs(newp),".",c='r')
    #plt.ylim(top=3.0,bottom=0.0)
    plt.show()
    plt.savefig("./dplots/resplot_"+auth+date+corr+itstr+".png",overwrite=True)
    plt.clf()

    for t in range(nb1):
        oldp=np.array(xx1[t])
        #newp=np.array(newpt[t])
        newv=np.array(newpt[t])
        #plt.plot(tt,np.abs(oldp),".",c='b')
        plt.plot(tt,np.abs(newv),".",c='r')
    plt.show()
    plt.savefig("./dplots/visplot_"+auth+date+corr+itstr+".png",overwrite=True)
    plt.clf()
    
    
    gf.dict_update(antinfo,'newres_'+corr+itstr,newpt)
    gf.dict_update(antinfo,'newvis_'+corr+itstr,newvis)

    return newres,newpt,antinfo,antinfo1

##################################################################
#Setting default value for TCLEAN stopcode
stop=0

#To time the code
start=time.time()

#this is for Venki data
target='sgr_apr07'

#for image names
case='halfant'
auth='Venki'
refmeth='newref'
date='11092022_'+str(attn)

#datams2=target+'_flagcor.ms'

#flagfile=datams2[:-3]+'_noflags.ms'
#os.system('rm -rf '+flagfile)
#split(vis=datams2,outputvis=flagfile,keepflags=False,datacolumn='data')

#Initial data file name and name for channel-averaged file
#datams1=target+'_flagcor_noflags.ms'
#datams1=flagfile
datams1='sgr_apr07_flag.ms'
tfile='sgr_apr07_flag_amp1.ms'
#os.system('rm -rf '+tfile)
#split(vis=datams1,outputvis=tfile,keepflags=False,datacolumn='data')
datams1=tfile

#Getting file name
dmsprefix=datams1[:-3]

#Setting name of file for time-averaged split file
dmsavg=dmsprefix+'_avg.ms'


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


clearcal(vis=datams1,spw='0,1,2,3',addmodel=True)


ms.open(datams1,nomodify=True)
visdata = ms.getdata(['antenna1','antenna2','data','flag','data_desc_id','sigma','axis_info'],ifraxis=True)
visdata['data'] = np.squeeze(visdata['data'])
ms.close()


allants=np.concatenate((visdata['antenna1'],visdata['antenna2']))
antlist=np.unique(allants)

print(antlist)

#Number of polarizations
npol=visdata['data'].shape[0]
polnames=visdata['axis_info']['corr_axis']
#corr=polnames[pol]


#Calculating number of antennas and baselines
nant=int(len(antlist))
nbl=int(nant*(nant-1)/2)
ntimes=len(visdata['axis_info']['time_axis']['MJDseconds'])

tt=np.array(visdata['axis_info']['time_axis']['MJDseconds'],dtype=float)
plt.clf()
fstore=np.ones_like(visdata['data'])
#nb=0
for pol in range(npol):
    #for n in range(ntimes):
    #    thistime=(visdata['axis_info']['time_axis']['MJDseconds']==tt[n])
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
                    pt=visdata['data'][pol][thisbase][0]
                    if ant1==tant or ant2==tant:
                        fstore[pol,nb,:]=2
                        plt.scatter(tt,abs(pt))
                    nb+=1
    plt.savefig("ttestplot"+date+str(pol)+".png")
    plt.clf()
#for pol in range(npol):
#    flnan=np.zeros_like(visdata['data'][pol])
#    fl=(visdata['data'][pol]==True)
#    flnan=np.squeeze(np.where(fl==True,np.nan,visdata['data'][pol]))
#    print("flnan shape",flnan.shape)
#    fstore[pol,:,:]=flnan
#    #visdata['data'][pol]=flnan

print(nb)

ms.open(datams1,nomodify=False)
nv1=ms.getdata(['data'],ifraxis=True)
nv1['data'][:,0,:,:]=fstore
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
    testdata = ms.getdata(['data','flag','data_desc_id','sigma','axis_info'],ifraxis=True)
    testdata['data'] = np.squeeze(testdata['data'])
    ms.close()

    xx=testdata['data'][i]
    for y in range(nbl):
        pp1=np.array(xx[y])
        #plt.scatter(tt,pp)
        plt.scatter(tt,np.abs(pp1))
    plt.show()
    #plt.savefig('datatestplot_'+auth+refmeth+'_'+date+corr+str(it)+'.png',overwrite=True)
    plt.savefig('flagtestplot_'+str(i)+'.png',overwrite=True)
    plt.clf()



print(nant)
print(nbl)

ELO=np.unique(visdata['axis_info']['time_axis']['MJDseconds'])
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

print(len(visdata['axis_info']['time_axis']['MJDseconds']))
print(len(np.unique(visdata['axis_info']['time_axis']['MJDseconds'])))



#print("num intervals: %f" %((np.max(ELO)-np.min(ELO))/240))
#print(.402600*240) 
#sys.exit()


#tsize,goodtimes=t_length(visdata=visdata,ELO=ELO,nbl=nbl)
#print("TSIZE")
#print(goodtimes)

#Splitting and channel averaging (if desired)
os.system('rm -rf '+dmsavg+' '+dmsavg+'.flagversions')
split(vis=datams1,
      outputvis=dmsavg,
      datacolumn='data')

rawdata=dmsavg
clearcal(vis=rawdata,spw='0,1,2,3',addmodel=True)

rawsplit=dmsprefix+'_rawsplit.ms'
os.system('rm -rf '+rawsplit+' '+rawsplit+'.flagversions')
split(vis=dmsavg,
      outputvis=rawsplit,
      datacolumn='data')


#Setting initial data file for analysis
datams0=rawsplit

#Getting scan numbers from file
tb.open(datams0,nomodify=True)
scan = tb.getcol('SCAN_NUMBER')

#Creating array to store residuals
visdata0 = tb.getcol('DATA')
resres=np.zeros_like(visdata0)
tb.close()

#Print scan numbers for edification
scan = np.unique(scan)
print(scan)



#Calculating RMS of data for CLEAN threshold
ms.open(datams0,nomodify=True)
sigs=ms.getdata(['sigma'])
fsigs=np.asarray(sigs['sigma']).flatten()
print("Fsig shape:",fsigs.shape)
invsigs=1./(fsigs*fsigs)
sumsigs=np.sum(invsigs)
rmst=2.*np.sqrt(1./sumsigs)
print("RMS: %e" %(rmst))
ms.close()



#Creating and initializing data file to store UVMULTIFIT pt source flux
g=open("uvmflux_"+target+".txt","w")
g.write("it scan ext mod moderr\n")

e=open("nv_"+target+".txt","w")
e.write("")



#Start with first iteration
it=0

#refant=gf.refantfinder(antlist=antlist,goodant=antinfo['goodant'])
#iref=np.where(antlist==refant)[0][0]

#refant=35 #DV20 index
#refant=14
#refant=8
iref=np.where(antlist==refant)[0][0]

print(tant)
print(refant)

#The wily (formerly) infinite while loop
#while(1):
while it<=0:
    '''
    e.write("ITERATION "+str(it)+"\n")

    #Declaring empty arrays for mod UVM flux and error
    muvmf=np.array([])
    muvmferr=np.array([])
    ruvmf=np.array([])

    #Creating file names for residual and model data from UVM
    datams_ext=dmsprefix+'_ext_it'+'{0:02d}'.format(it)+'.ms'
    datams_mod=dmsprefix+'_mod_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files
    os.system('rm -rf '+datams_ext)
    split(vis=datams0,outputvis=datams_ext,datacolumn='data')

    #cvis.append(str(datams_mod))
    os.system('rm -rf '+datams_mod)
    split(vis=datams0,outputvis=datams_mod,datacolumn='data')

    #Clearing/creating corrected and model columns
    clearcal(vis=datams_ext,spw='0,1,2,3',addmodel=True)
    clearcal(vis=datams_mod,spw='0,1,2,3',addmodel=True) 
    
    #Cycle through each scan
    for s in scan:
        #UVM for residuals of fit
        sc=int(s)
        myuvfit_ext1 = uvm.uvmultifit(vis=datams_ext,
                    spw='0,1,2,3',
                    scans=[[sc]],
                    model=['delta'],
                    var=['0,0,p[0]'],
                    p_ini=[0.],
                    OneFitPerChannel=False,
                    column='data',
                    write='residuals') # will write residuals in the 'corrected' column, but only for the subset of fitted channels !
        # to run on the continuum collapsed cube
        ruvmf=np.append(ruvmf,myuvfit_ext1.result['Parameters'][0])

        #UVM for model visibilities
        myuvfit_mod = uvm.uvmultifit(vis=datams_mod,
                    spw='0,1,2,3',
                    scans=[[sc]], 
                    model=['delta'],
                    var=['0,0,p[0]'], 
                    p_ini=[0.], 
                    OneFitPerChannel=False,
                    column='data',
                    write='model') # will write best-fit model in the 'model' column, but only for the subset of fitted channels !
        tb.open(datams_mod)
        md=tb.getcol('MODEL_DATA')
        tb.close()
        #Appending results of UVM to list
        muvmf=np.append(muvmf,myuvfit_mod.result['Parameters'][0])
        muvmferr=np.append(muvmferr,myuvfit_mod.result['Uncertainties'][0])

    #Pulling residuals of UVM pt source fitting for analysis
    tb.open(datams_ext)
    ext1=tb.getcol('CORRECTED_DATA')
    tb.close()

    #Name of file to split these residuals into
    datams_ext2=dmsprefix+'_ext2_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files
    #The datacolumn here means nothing, it will be replaced
    os.system('rm -rf '+datams_ext2)
    split(vis=datams_ext,outputvis=datams_ext2,datacolumn='data')

    #Making the residuals of the first UVMULTIFIT the data column of the new file
    tb.open(datams_ext2,nomodify=False)
    tb.putcol('DATA', ext1)
    tb.close()

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
    '''
    #File where final set of residuals is the data
    datams_ext3=dmsprefix+'_ext3_it'+'{0:02d}'.format(it)+'.ms'
    
    #Splitting and creating the files where second 
    #round of residuals is the data column
    
    '''
    os.system('rm -rf '+datams_ext3+' '+datams_ext3+'.flagversions')
    split(vis=datams_ext2,outputvis=datams_ext3,datacolumn='corrected')
    '''
    '''
    tb.open(datams_ext2)
    ext2=tb.getcol('CORRECTED_DATA')
    tb.close()

    #Splitting and creating the files
    #The datacolumn here means nothing, it will be replaced
    os.system('rm -rf '+datams_ext3+' '+datams_ext3+'.flagversions')
    split(vis=datams_ext2,outputvis=datams_ext3,datacolumn='data')
    clearcal(vis=datams_ext3,spw='0,1,2,3',addmodel=True)
    
    #Making the residuals of the first UVMULTIFIT the data column of the new file
    tb.open(datams_ext3,nomodify=False)
    tb.putcol('DATA', ext2)
    tb.close()

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
    imname1=dmsprefix+'_it'+'{0:02d}'.format(it)+'_cmodel'
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
    '''
    #Here


        
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
    newvis=np.zeros_like(visdata['data'])
    newvis1=np.zeros_like(visdata['data'])
    antinfo=dict()
    antinfo1=dict()
    plt.clf()
    for p in range(npol):
        plt.clf()
        newvis[p],newvis1[p],newantinfo,newantinfo1=lls(pol=p,corr=polnames[p],datams1=rawdata,target=target,case=case,auth=auth,refmeth=refmeth,date=date,it=it,dvis=rawdata,rfant=refant)
        plt.clf()
        #newvis1[p],newantinfo1=lls(pol=p,corr=polnames[p],datams1=datams_ext3,target=target,case=case,auth=auth,refmeth=refmeth,date=date,it=it,dvis=rawdata)
        antinfo.update(newantinfo)
        antinfo1.update(newantinfo1)
    

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


    ms.open(rawdata,nomodify=False)
    nv1=ms.getdata(['corrected_data'],ifraxis=True)
    nv1['corrected_data'][:,0,:,:]=newvis1
    ms.putdata(nv1)
    ftnv1=nv1['corrected_data'][0]
    e.write(str(ftnv1)+"\n")
    print(ftnv1)
    ms.close()

    ms.open(rawdata,nomodify=True)
    nnv1=ms.getdata(['corrected_data'],ifraxis=True)
    ftnnv1=nnv1['corrected_data'][0]
    e.write(str(ftnnv1)+"\n")
    print(ftnnv1)
    ms.close()



    #Dirty image of gain-calibrated data
    imname3=dmsprefix+'_it'+'{0:02d}'.format(it)+'_selfcalmodel'
    os.system('rm -rf '+imname3+'.*')
    r2c1=tclean(vis=rawdata,
           imagename=imname3,
           spw='0,1,2,3',
           datacolumn='corrected',
           niter=0,
           threshold=rmst,
           pblimit=-1,
           imsize=[300,300],
           interactive=0,
           cell='0.2arcsec')

    #Dirty image of gain-calibrated residuals
    imname4=dmsprefix+'_it'+'{0:02d}'.format(it)+'_selfcalmodel_ext'
    os.system('rm -rf '+imname4+'.*')
    r2c2=tclean(vis=datams_ext3,
           imagename=imname4,
           spw='0,1,2,3',
           datacolumn='corrected',
           niter=0,
           threshold=rmst,
           pblimit=-1,
           imsize=[300,300],
           interactive=0,
           cell='0.2arcsec')


    #Split again so gain-calibrated data and residuals are the data
    dms_selfcal=dmsprefix+"_selfcal_it"+str(it)+".ms"
    dms_selfcalf=dms_selfcal+'.flagversions'
    dms_selfcal_copy=dmsprefix+"_selfcal_copy_it"+str(it)+".ms"
    dms_selfcalf_copy=dms_selfcal_copy+'.flagversions'
    dms_selfcal_ext=dmsprefix+"_selfcal_ext_it"+str(it)+".ms"
    dms_selfcalf_ext=dms_selfcal_ext+'.flagversions'

    os.system('rm -rf '+dms_selfcal+' '+dms_selfcalf)
    split(vis=rawdata,outputvis=dms_selfcal,datacolumn='corrected')

    os.system('rm -rf '+dms_selfcal_copy+' '+dms_selfcalf_copy)
    split(vis=rawdata,outputvis=dms_selfcal_copy,datacolumn='corrected')

    os.system('rm -rf '+dms_selfcal_ext+' '+dms_selfcalf_ext)
    split(vis=datams_ext3,outputvis=dms_selfcal_ext,datacolumn='corrected')


    #Seeing why TCLEAN ended
    scode=r2c['stopcode']

    #Pulling gain-calibrated residuals to add to previous iterations
    tb.open(dms_selfcal_ext)
    resids=tb.getcol('DATA')
    tb.close()
    resres+=resids

    #Pulling gain-calibrated full data
    tb.open(dms_selfcal)
    print("Getting",it)
    scdata=tb.getcol('DATA')
    tb.close()

    #Subtracting residuals from full data to hopefully get point source
    dminres=scdata-resids

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


    print("MAX:",madmax)
    if stop==2: break
    if scode!=1: stop=2

    it+=1

#Print number of iterations performed
print("Performed %i iterations" %(it))


dms_selfcal1=dmsprefix+"_selfcal_copy_it"+str(it-1)+".ms"
dms_selfcal1_ext=dmsprefix+"_selfcal_ext_it"+str(it-1)+".ms"


tb.open(dms_selfcal1_ext)
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
split(vis=dms_selfcal1,outputvis=datams_ext3_sc,datacolumn='data')

tb.open(datams_ext3_sc,nomodify=False)
tb.putcol('DATA',dminres)
tb.close()


muvmf=np.array([])
muvmferr=np.array([])

#Creating file names for residual and model data from UVM
datams_mod=dmsprefix+'_mod_final.ms'

#cvis.append(str(datams_mod))
os.system('rm -rf '+datams_mod)
split(vis=datams_ext3_sc,outputvis=datams_mod,datacolumn='data')

#Clearing/creating corrected and model columns
clearcal(vis=datams_mod,spw='0,1,2,3',addmodel=True) 

#Cycle through each scan
for s in scan:
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
split(vis=dmsavg,outputvis=datams_fin,datacolumn='data')

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

