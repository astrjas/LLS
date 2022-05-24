#importing all my junk
#from NordicARC import uvmultifit as uvm
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
import casatools

def nonan(arr):
    return np.where(np.isfinite(arr),arr,0)

def refantfinder(antlist,goodant):
    maxgood=0
    rant=0
    for ant in range(len(antlist)):
        ng=len(np.where(goodant[:,ant]==False)[0])
        #print(ng)
        if ng>maxgood:
            maxgood=ng
            rant=ant
    print("REFANT:%i"%(rant))
    return rant

def arr_length_new(visdata,curr_time):
    thistime=(visdata['axis_info']['time_axis']['MJDseconds']==curr_time)
    ant_goodness=[]
    #gant1=[]
    #gant2=[]
    for ant1 in np.unique(visdata['antenna1']):
        for ant2 in np.unique(visdata['antenna2']):
            if ant1 < ant2:
                thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
                if thisbase.sum()>0:
                    #print(len(thisbase==True))
                    #print(visdata['data'][0][thisbase][0][thistime][0])
                    pt=visdata['data'][0][thisbase][0][thistime][0]
                    #print(pt)
                    #print(pt.shape)
                    amp=np.absolute(pt)
                    #print(amp)
                    if amp<=0:
                        ant_goodness.append(0)
                        #continue
                        #print(pt)
                        #gant1.append(ant1)
                        #gant2.append(ant2)
                    else:
                        ant_goodness.append(1)
                        #gant1.append(ant1)
                        #gant2.append(ant2)

    return ant_goodness,gant1,gant2
    
def antdicts(visdata,ELO):
    ant1dict={}
    ant2dict={}
    ant_good_v={}
    for time in ELO:
        ant_01,a1list,a2list=arr_length_new(visdata=visdata,curr_time=time)
        
        '''
        if len(a1list)<=1 or len(a2list)<=1: 
            sbc+=1
            continue
        '''

        ant1dict.update({time:a1list})
        ant2dict.update({time:a2list})
        ant_good_v.update({time:ant_01})
    return ant1dict,ant2dict,ant_good_v

def dict_update(dic,key,val):
    dic.update({key:val})

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

def th_Ir(Th_r,th_m):
    t1=np.linalg.inv(np.matmul(Th_r.T,Th_r))
    t2=np.matmul(t1,Th_r.T)
    theta_Ir=np.matmul(t2,th_m)
    return theta_Ir


def l_Ir(ll_r,ll_m):
    #print("ll_r")
    #print(ll_r)
    #print("ll_m")
    #print(ll_m)
    #ll_r_T=np.transpose(ll_r,(1,0))
    t0=np.matmul(ll_r.T,ll_r)
    #print(t0[np.nonzero(t0)])
    t1=np.linalg.inv(np.matmul(ll_r.T,ll_r))
    #print("t1")
    #print(t1)
    t2=np.matmul(t1,ll_r.T)
    #print("t2")
    #print(t2)
    ll_Ir=np.matmul(t2,ll_m)
    #print("answer!")
    #print(ll_Ir)
    return ll_Ir
'''
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

pol=1
it=0
itstr=str(it)
#def lls(pol,corr,datams1,target,case,auth,refmeth,date,it):
#thisant=(visdata['antenna1']==ant) | (visdata['antenna2']==ant)

#setting a dummy refant
#refant=0

#this is for Venki data
target='sgr_apr07'

#for image names
case='noshift_noflag'
auth='Venki'
refmeth='refanctfunc'
date='04182022'

#datams0='firsttwointflag_QSOB19.ms'
#os.system('rm -rf QSOB19_chanavg.ms')
#split(vis=datams0,outputvis='QSOB19_chanavg.ms',width=240,datacolumn='data')

#Initial data file name and name for channel-averaged file
datams2=target+'_flagcor.ms'
#datams2='Sagittarius_A_star_C.ms_CALIBRATED_SELFCAL'
#datams2='uvmfit_3d77_SGRA_onechan.ms'
#datams2='QSOB19_chanavg.ms'

flagfile=datams2[:-3]+'_noflags.ms'
os.system('rm -rf '+flagfile)
split(vis=datams2,outputvis=flagfile,keepflags=False,datacolumn='data')
'''

###############
def lls(pol,corr,datams1,target,case,auth,refmeth,date,it):
    #datams1=flagfile
    MS=casatools.ms()


    #Opening data and pulling necessary info
    MS.open(datams1,nomodify=True)
    MS.selectinit(reset=True)
    visdata = MS.getdata(['antenna1','antenna2','data','axis_info','flag'],ifraxis=True)
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

    allant=np.zeros((nant,ntimes),dtype=bool)
    alltim=np.zeros((ntimes,nant),dtype=bool)

    lp=visdata['data'][pol]<0.1

    tt=np.array(visdata['axis_info']['time_axis']['MJDseconds'],dtype=float)



    #cycle through each ant
    #a=0
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
        #a+=1


    #lp1=np.zeros_like(lp)
    lp1=np.where(lp==True,np.nan,visdata['data'][pol])


    #Creating a tf dictionary for each ant and whether time is good or bad there
    antinfo=dict()
    antinfo['timestamp']=visdata['axis_info']['time_axis']['MJDseconds']

    allant=np.zeros((nant,ntimes),dtype=bool)
    alltim=np.zeros((ntimes,nant),dtype=bool)



    #TENTATIVE BEGINNING OF WHAT WILL BE IN MASSIVE LOOP
    #XX correlation
    xx=np.zeros_like(visdata['data'][pol])
    lowpt=np.abs(visdata['data'][pol]<0.1)
    xx=np.where(lowpt==True,np.nan,visdata['data'][pol])
    print("xxshape",xx.shape)

    for y in range(nbl):
        pp1=np.array(xx[y])
        #plt.scatter(tt,pp)
        plt.scatter(tt,np.abs(pp1))
    plt.show()
    plt.savefig('datatestplot_'+auth+refmeth+'_'+date+corr+str(it)+'.png',overwrite=True)
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

    refant=gf.refantfinder(antlist=antlist,goodant=antinfo['goodant_'+corr+itstr])
    iref=np.where(antlist==refant)[0]
    #refant=0


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
                        ph=np.angle(pt,deg=True)
                        amp=np.absolute(pt)

                        #if iant1==5 or iant2==5:
                            #amp/=2
                            #ph+=90
                            #print("bing")

                        l_m[nb,0]=np.log10(amp)

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

        gf.dict_update(theta_r_dict,time+corr+itstr,Theta_r)
        gf.dict_update(theta_m_dict,time+corr+itstr,theta_m)


        gf.dict_update(script_L_dict,time+corr+itstr,script_L)
        gf.dict_update(l_m_dict,time+corr+itstr,l_m)


        script_L_nn=gf.nonan(script_L)
        l_m_nn=gf.nonan(l_m)

        l_r=gf.l_Ir(ll_r=script_L_nn,ll_m=l_m_nn)


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
        gf.dict_update(l_del_dict,time+corr+itstr,l_del)
        #print(l_del.shape)

        l_Ir_res=gf.l_Ir(ll_r=script_L_nn,ll_m=l_del)

        #print("l_Ir w/ resids: \n"+str(l_Ir_res))

        theta_del=theta_m_nn-np.matmul(Theta_r_nn,theta_Ir)
        #print("theta_del \n"+str(theta_del))
        gf.dict_update(theta_del_dict,time+corr+itstr,theta_del)

        theta_Ir_res=gf.th_Ir(Th_r=Theta_r_nn,th_m=theta_del)
        #print("theta_Ir w/ resids: \n"+str(theta_Ir_res))


        l_Ir_final=l_r+l_Ir_res
        #print(l_Ir_final.shape)
        gf.dict_update(l_r_dict,time+corr+itstr,l_Ir_final)
        #print("final amps (log form) \n"+str(l_Ir_final))

        Ir_converted=np.zeros_like(l_Ir_final)

        Ir_converted=np.vstack([10.0**x for x in l_Ir_final])
        #print(Ir_converted)
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

        gf.dict_update(l_Ir_conv_dict,time+corr+itstr,Ir_converted)

        theta_Ir_final=theta_Ir+theta_Ir_res
        gf.dict_update(theta_Ir_dict,time+corr+itstr,theta_Ir_final)



    ELO_diff=np.max(ELO)-np.min(ELO)
    nint=int(ELO_diff/240)
    ELO_range=np.linspace(np.min(ELO),np.max(ELO),nint)


    fkeys=theta_Ir_dict.keys()
    tkeys=antinfo['timestamp']

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

    for t_step in range(len(ELO_range)-1):
        nd=len([antinfo['timestamp'][x] for x in range(len(tkeys)) if (tkeys[x]>=ELO_range[t_step] and tkeys[x]<=ELO_range[t_step+1])])
        for k in range(nant):
            avga[k,0,t_step]=np.mean([l_Ir_conv_dict[str(t_step)+corr+itstr][k] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])
            if k==refant and nd!=0: 
                avgp[k,0,t_step]=0
            elif k<refant and nd!=0:
                avgp[k,0,t_step]=np.mean([theta_Ir_dict[str(t_step)+corr+itstr][k] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])
            elif k>refant and nd!=0:
                avgp[k,0,t_step]=np.mean([theta_Ir_dict[str(t_step)+corr+itstr][k-1] for x in range(ntimes) if antinfo['timestamp'][x]>=ELO_range[t_step] and antinfo['timestamp'][x]<=ELO_range[t_step+1]])


        if np.isnan(np.mean(avga[:,0,t_step]))==True or np.isnan(np.mean(avgp[:,0,t_step]))==True:
            if nd==0:
                print("That's okay!")
                avga[:,0,t_step]=np.nan
                avgp[:,0,t_step]=np.nan
                nemp+=1
            else: 
                nbad+=1
        else: nonan+=1

    antinfo['avg_phase_'+corr+itstr]=np.squeeze(avgp)
    antinfo['avg_amp_'+corr+itstr]=np.squeeze(avga)    


    for a in range(nant):
        plt.scatter(ELO_range[:-1],avga[a,0,:],color='black',marker='x')
        plt.scatter(ELO_range[:-1],avgp[a,0,:],color='red')
    plt.legend()
    plt.show()

    plt.savefig('./dplots/ndpv_'+date+'_'+refmeth+'_ap_'+auth+'_'+case+corr+itstr+'.png')

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
    tb1=0

    for t_step in range(len(ELO_range)-1):
        tchunk=[antinfo['timestamp'][x] for x in range(len(tkeys)) if (tkeys[x]>=ELO_range[t_step] and tkeys[x]<=ELO_range[t_step+1])]
        for i in range(len(tchunk)):
            thistime=(visdata['axis_info']['time_axis']['MJDseconds']==tchunk[i])
            time=i

            gt=antinfo['goodt_'+corr+itstr][i]

            #Pulling antennas that are bad for this time
            #Calculating number of antennas and baselines
            #-1 is to account for bad ants (ant1=-1 is not valid) 
            #need to rework how pt is pulled
            nb1=0
            for ant1 in np.unique(visdata['antenna1']):
                for ant2 in np.unique(visdata['antenna2']):
                    if ant1 < ant2:
                        thisbase = (visdata['antenna1']==ant1) & (visdata['antenna2']==ant2)
                        iant1=np.where(antlist==ant1)[0][0]
                        iant2=np.where(antlist==ant2)[0][0]
                        if thisbase.sum()>0:
                            pt=xx[thisbase][0][tb1]
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
                            newpt[nb1,tb1]=pt/(ga1*ga2*np.exp(1.0j*(gp1-gp2)))
                            nb1+=1
            tb1+=1


    for t in range(nb1):
        oldp=np.array(xx[t])
        newp=np.array(newpt[t])
        plt.plot(tt,np.abs(oldp),".",c='b')
        plt.plot(tt,np.abs(newp),".",c='r')
    plt.show()
    plt.savefig("./dplots/visplot_"+auth+date+corr+itstr+".png",overwrite=True)
    plt.clf()
    
    gf.dict_update(antinfo,'newvis_'+corr+itstr,newpt)

    return newpt,antinfo

