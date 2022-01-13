import numpy as np

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
