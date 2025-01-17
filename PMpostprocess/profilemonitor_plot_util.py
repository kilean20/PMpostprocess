import sys
#import commands
import os
from os import path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
import math
import collections
from collections import OrderedDict
from datetime import datetime
import json
from epics import caput,caget
#os.envron['EPICS_CA_MAX_ARRAY_BYTES'] = '100000000'


#uA = 1.00E+6
uA = 1.00
plt.rcParams.update({'font.size':14})
bgth = 3
bgth1 = 3.

sq2 = math.sqrt(2)
deg = math.pi/180.

nave = 3
ndatBG = 30
plot_dict = {}
class profilemonitor_plot:
    #def __init__(self,debug = False, jupyter = False, bgfunc = 'pol1'):
    def __init__(self,debug = False, jupyter = False, batch = False,  bgfunc = 'const'):
        self.debug = debug
        self.bgfunc = bgfunc
        self.jupyter = jupyter
        self.batch = batch
        
        print ("profilemonitor_plot_util start")
        self.color = {"x":"blue", "y":"red", "u":"green", "v":"black"}

    def get_param(self,inp):
        par_dict = {}
        #global PM
        
        
        if "dat" in inp:
            live = False
        else:
            live = True

        script_dir = os.path.dirname(os.path.abspath(__file__))
        floc = os.path.join(script_dir,'profilemonitor_param.dat')
        if os.path.isfile(floc):
            floc = floc
        elif os.path.isfile('/media/sf_WorkonDropBox/BeamStudy/script/notebooks/ProfileMonitorAnalysis/profilemonitor_param.dat'):
            floc = '/media/sf_WorkonDropBox/BeamStudy/script/notebooks/ProfileMonitorAnalysis/profilemonitor_param.dat'
        elif os.path.isfile("/Users/maruta/Dropbox/FRIB_maruta/BeamStudy/script/ProfileMonitor/notebooks/ProfileMonitorAnalysis/profilemonitor_param.dat"):
            floc = '/Users/maruta/Dropbox/FRIB_maruta/BeamStudy/script/ProfileMonitor/notebooks/ProfileMonitorAnalysis/profilemonitor_param.dat'
        elif os.path.isfile('/files/shared/ap/phyapp_notebooks/ProfileMonitorAnalysis/profilemonitor_param.dat'):
            floc = '/files/shared/ap/phyapp_notebooks/ProfileMonitorAnalysis/profilemonitor_param.dat'
        elif os.path.isfile('/usr/lib/python3/dist-packages/ProfileMonitorAnalysis/profilemonitor_param.dat'):
            floc = '/usr/lib/python3/dist-packages/ProfileMonitorAnalysis/profilemonitor_param.dat'
        f = open(floc,"r")

        a = f.readlines()
        f.close()
        ipm = 0
        for l in a:
            if l[0] == "#":
                continue
            ipm += 1
            e = l.split()
            if (e[0].replace(":","_") in inp.replace(":","_")) or (str(ipm) == inp):
                break
        e = l.split()
        name = e[0]
        coord = e[1]
        par_dict['live'] = live
        par_dict['name'] = name
        par_dict['coord'] = coord
        par_dict['scan_startpos1'] = float(e[2])
        par_dict['scan_stoppos1'] = float(e[3])
        par_dict['scan_startpos2'] = float(e[4])
        par_dict['scan_stoppos2'] = float(e[5])
        par_dict['scan_startpos3'] = float(e[6])
        par_dict['scan_stoppos3'] = float(e[7])
        par_dict['offset1'] = [float(e[8]),float(e[9]),float(e[10])] #= offset written in param file
        par_dict['offset0'] = [0.0]*3 #= offset set online
        par_dict['ang1'] = float(e[14]) #= injection angle of arm1
        par_dict['ang2'] = float(e[15]) #= injection angle of arm2
        #par_dict['wire_diam'] = float(e[16]) #= injection angle of arm2
        par_dict['wire_diam'] = float(e[16])/np.sqrt(12) #= size resolution by wire diameter
        #return name,live,coord,offset

        return par_dict

    #===
    def get_rawdata(self,inp,par_dict):
        raw_dict = {}
        name = par_dict['name']
        if par_dict['live']: #= read EPICS record
            pass
            '''
            if par_dict['coord'][0] == "F":
                #a1p = commands.getoutput("caget -a "+name+":DRV_PPOTX")
                #a1s1= commands.getoutput("caget -a "+name+":DRV_BCX")
                #a2p = commands.getoutput("caget -a "+name+":DRV_PPOTY")
                #a2s1= commands.getoutput("caget -a "+name+":DRV_BCY")

                raw_dict['a1p'] = commands.getoutput("caget -a "+name+":DRV1_PPOT")
                raw_dict['a1p1']= commands.getoutput("caget -a "+name+":DRV1_PPOT1")
                raw_dict['a1s1']= commands.getoutput("caget -a "+name+":DRV1_BC1")
                raw_dict['a1o1']= commands.getoutput("caget -a "+name+":DRV1_OFST1")
                raw_dict['a2p'] = commands.getoutput("caget -a "+name+":DRV2_PPOT")
                raw_dict['a2p1']= commands.getoutput("caget -a "+name+":DRV2_PPOT1")
                raw_dict['a2s1']= commands.getoutput("caget -a "+name+":DRV2_BC1")
                raw_dict['a2o1']= commands.getoutput("caget -a "+name+":DRV2_OFST1")
            elif par_dict['coord'][0] == 'S':
                raw_dict['a1p'] = commands.getoutput("caget -a "+name+":DRV_PPOT")
                raw_dict['a1p1'] = commands.getoutput("caget -a "+name+":DRV_PPOT1")
                raw_dict['a1p2'] = commands.getoutput("caget -a "+name+":DRV_PPOT2")
                raw_dict['a1p3'] = commands.getoutput("caget -a "+name+":DRV_PPOT3")
                raw_dict['a1s1']= commands.getoutput("caget -a "+name+":DRV_BC1")
                raw_dict['a1s2']= commands.getoutput("caget -a "+name+":DRV_BC2")
                raw_dict['a1s3']= commands.getoutput("caget -a "+name+":DRV_BC3")
                raw_dict['a1o1']= commands.getoutput("caget -a "+name+":DRV_OFST1")
                raw_dict['a1o2']= commands.getoutput("caget -a "+name+":DRV_OFST2")
                raw_dict['a1o3']= commands.getoutput("caget -a "+name+":DRV_OFST3")
            else:
                #= 6 in fork
                raw_dict['a1p'] = commands.getoutput("caget -a "+name+":DRV1_PPOT")
                raw_dict['a1p1']= commands.getoutput("caget -a "+name+":DRV1_PPOT1")
                raw_dict['a1s1']= commands.getoutput("caget -a "+name+":DRV1_BC1")
                raw_dict['a1o1']= commands.getoutput("caget -a "+name+":DRV1_OFST1")
                #= 12 in fork
                raw_dict['a2p'] = commands.getoutput("caget -a "+name+":DRV2_PPOT")
                raw_dict['a2p1']= commands.getoutput("caget -a "+name+":DRV2_PPOT1")
                raw_dict['a2p2']= commands.getoutput("caget -a "+name+":DRV2_PPOT2")
                raw_dict['a2s1']= commands.getoutput("caget -a "+name+":DRV2_BC1")
                raw_dict['a2s2']= commands.getoutput("caget -a "+name+":DRV2_BC2")
                raw_dict['a2o1']= commands.getoutput("caget -a "+name+":DRV2_OFST1")
                raw_dict['a2o2']= commands.getoutput("caget -a "+name+":DRV2_OFST2")
            '''
        else:
            f = open(inp,"r")
            a = f.readlines()
            f.close()
            for l in a:
                # Large profile monitor
                if "DRV1_PPOT1" in l:
                    raw_dict['a1p1'] = l
                elif "DRV1_PPOT" in l:
                    raw_dict['a1p'] = l
                elif "DRV1_BC1" in l:
                    raw_dict['a1s1'] = l
                elif "DRV1_OFST1" in l:
                    raw_dict['a1o1'] = l
                elif "DRV2_PPOT1" in l:
                    raw_dict['a2p1'] = l
                elif "DRV2_PPOT2" in l:
                    raw_dict['a2p2'] = l
                elif "DRV2_PPOT" in l:
                    raw_dict['a2p'] = l
                elif "DRV2_BC1" in l:
                    raw_dict['a2s1'] = l
                elif "DRV2_BC2" in l:
                    raw_dict['a2s2'] = l
                elif "DRV2_OFST1" in l:
                    raw_dict['a2o1'] = l
                elif "DRV2_OFST2" in l:
                    raw_dict['a2o2'] = l
                #= Flapper
                elif "DRV1_PPOTX" in l:
                    raw_dict['a1p1'] = l
                elif "DRV1_PPOTY" in l:
                    raw_dict['a2p1'] = l
                elif "DRV_BCX" in l:
                    raw_dict['a1s1'] = l
                elif "DRV_BCY" in l:
                    raw_dict['a2s1'] = l
                #= Small profile monitor
                elif "DRV_PPOT1" in l:
                    raw_dict['a1p1'] = l
                elif "DRV_PPOT2" in l:
                    raw_dict['a1p2'] = l
                elif "DRV_PPOT3" in l:
                    raw_dict['a1p3'] = l
                elif "DRV_PPOT" in l:
                    raw_dict['a1p'] = l
                elif "DRV_BC1" in l:
                    raw_dict['a1s1'] = l
                elif "DRV_BC2" in l:
                    raw_dict['a1s2'] = l
                elif "DRV_BC3" in l:
                    raw_dict['a1s3'] = l
                elif "DRV_OFST1" in l:
                    raw_dict['a1o1'] = l
                elif "DRV_OFST2" in l:
                    raw_dict['a1o2'] = l
                elif "DRV_OFST3" in l:
                    raw_dict['a1o3'] = l


        #print raw_dict.keys()
        return raw_dict

    #===
    def get_middlepos(self,epos,esig1,esig2, par, off1, off2, coord):
        #print ('off1: '+str(off1)+', off2: '+str(off2))
        n = len(epos)
        dat = []
        sig1l,sig2l = [],[]
        for (post,sig1t,sig2t) in zip(epos[4:],esig1[4:],esig2[4:]):
            if post == "0":
                break
            sig1l.append(float(sig1t))
            sig2l.append(float(sig2t))

        #sigmin,sigmax = min(sig1l+sig2l),max(sig1l+sig2l)
        sigmin1,sigmax1 = min(sig1l),max(sig1l)
        sigmin2,sigmax2 = min(sig2l),max(sig2l)
        #sigth1 = 0.2*(sigmax1-sigmin1)+sigmin1
        #sigth2 = 0.2*(sigmax2-sigmin2)+sigmin2
        sigth1 = 0.2*sigmax1
        sigth2 = 0.2*sigmax2
        if self.debug:print ("sigth1: "+str(sigth1)+", sigth2: "+str(sigth2))

        if coord[0] == 'S':
            ran = 20.
        else:
            ran = 40.
            
        for (post,sig1t,sig2t) in zip(epos[4:],esig1[4:],esig2[4:]):
            if post == "0":
                break
            '''
            if (float(sig1t) > sigth1) or (float(sig2t) > sigth2):
                sigsum = float(sig1t)+float(sig2t)*par
                dat.append([float(post),float(sig1t)+float(sig2t)*par])
            '''

            if ((float(sig1t) > sigth1) and (-1*ran < float(post)+off1 < ran)) or ((float(sig2t) > sigth2) and (-1*ran < float(post)+off2 < ran)):
                dat.append([float(post),float(sig1t)+float(sig2t)*par])                
        c = self.calc_center(dat)
        if self.debug:print ("middle: "+str(c))
        return c
    #===
    def get_mean_stdev(self,dat):
        n = len(dat)
        c,stdev = 0.0, 0.0
        for i in range(n):
            c+= dat[i][1]
        c/=float(n)
        for i in range(n):
            stdev += math.pow(dat[i][1]-c,2)
        stdev = math.sqrt(stdev/float(n-1))

        return c,stdev
    #===
    def get_range(self, dat, mid1, mid2, par_dict):
        # Detemination of range for parameter calculation
        ndat = len(dat)
        #print('mid1: '+str(mid1)+'\tmid2: '+str(mid2))
        if mid1 < -900. and mid2 > 900.: # no cross-talk
            ran = [dat[0][0]+1, dat[ndat-1][0]-1]
        elif mid1 > -900. and mid2 > 900.:
            ran = [mid1, dat[ndat-1][0]-1]
        elif mid1 < -900. and mid2 < 900.:
            ran = [dat[0][0]+1, mid2]
        else:
            ran = [mid1, mid2]

        if self.debug: print ("ran: "+str(ran))

        return ran

    #===
    def calc_center(self,dat):
        ndat = len(dat)
        s,c = 0., 0.
        for i in range(ndat-1):
            pos1,sig1 = dat[i+1][0],dat[i+1][1]
            pos2,sig2 = dat[i][0],dat[i][1]

            sig = (sig1+sig2)/2.
            pos = (pos1+pos2)/2.
            dpos= math.fabs(pos2-pos1)
            s += sig*dpos
            c += sig*pos*dpos
        if s>0.:
            c/=s
        return c
    #===
    def calc_background(self, dat, ran):
        if self.debug: 'In calc_background'
        n = len(dat)
        p1s,p2s =[],[]
        s1s,s2s =[],[]
        '''
        for i in range(n):
            if (ran[0] <= dat[i][0] < ran[0]+5.):
                p1s.append(dat[i][0])
                s1s.append(dat[i][1])
            elif (ran[1] - 5. < dat[i][0] <= ran[1]):
                p2s.append(dat[i][0])
                s2s.append(dat[i][1])
        '''
        ndat = 0
        print('ran[0]: '+str(ran[0])+'\tran[1]: '+str(ran[1]))
        for i in range(n):
            #if (ran[0] <= dat[i][0] < ran[1]+10.):
            if (ran[1]-2 <= dat[i][0]) and (dat[i][0] < ran[1]+2):
                p1s.append(dat[i][0])
                s1s.append(dat[i][1])
                ndat += 1
            if ndat==ndatBG:
                break
        ndat = 0
        for i in range(n-1, -1, -1):
            #if (ran[0] - 10. < dat[i][0] <= ran[1]):
            if (ran[0]-2 < dat[i][0]) and (dat[i][0] <= ran[0]+2):
                p2s.append(dat[i][0])
                s2s.append(dat[i][1])
                ndat += 1
            if ndat==ndatBG:
                break

        n1,n2 = len(p1s),len(p2s)
        if n1>0:
            print('n1: '+str(n1)+'\t'+str(p1s[0])+'\t'+str(p1s[-1]))
        else:
            print('n1: '+str(n1))
        if n2>0:
            print('n2: '+str(n2)+'\t'+str(p2s[0])+'\t'+str(p2s[-1]))
        else:
            print('n2: '+str(n2))
            
        if len(p1s+p2s)>0:
            if 'const' in self.bgfunc: #= constant background
                cb1,stdv1 = -999.,999.
                cb2,stdv2 = -999.,999.
                if n1>0 and n2>0:
                    cb1,stdv1 = np.mean(s1s),np.std(s1s)
                    cb2,stdv2 = np.mean(s2s),np.std(s2s)
                    if abs(cb1-cb2)<bgth*min(stdv1,stdv2):
                        cb,stdv=(cb1*n1+cb2*n2)/(n1+n2),(stdv1*n1+stdv2*n2)/(n1+n2)
                    else:
                        cb,stdv=min(cb1,cb2),min(stdv1,stdv2)

                    if cb+bgth1*stdv<0.:
                        cb,stdv=max(cb1,cb2),max(stdv1,stdv2)
                elif n1>0 and n2 == 0:
                    cb1,stdv1 = np.mean(s1s),np.std(s1s)
                    cb,stdv = cb1, stdv1
                elif n2>0 and n1 == 0:
                    cb2,stdv2 = np.mean(s2s),np.std(s2s)
                    cb,stdv = cb2,stdv2

                else:
                    cb,stdv = 0.0, 0.0
                '''
                if cb+bgth1*stdv<0.:
                    cb = 0.
                    stdv = 0
                '''

                coefs = [cb]
                if self.debug:print ("cb,stdv,cb1,stdv1,cb2,stdv2: "+str(cb)+"\t"+str(stdv)+"\t"+str(cb1)+"\t"+str(stdv1)+"\t"+str(cb2)+"\t"+str(stdv2))
            elif 'pol' in self.bgfunc: #= nth order polynominal function
                norder = int(self.bgfunc[len(self.bgfunc)-1])
                coefs = np.polyfit(p1s+p2s, s1s+s2s, norder)
                sts = []
                for p,s in zip(p1s+p2s,s1s+s2s):
                    st = s
                    for j in range(norder+1):
                        st-=coefs[j]*math.pow(p,norder-j)
                    sts.append(st)
                stdv = np.std(sts)
        else:
            coefs = [0]
            stdv = 0.
                
        return coefs, stdv
        
    #===
    def calc_percent_beamsize(self, dat, sum0, cen, ratio):
        ndat = len(dat)
        dxmax1 = math.fabs(dat[0][0]-cen)
        dxmax2 = math.fabs(dat[ndat-1][0]-cen)
        dxmax = max(dxmax1, dxmax2)
        if self.debug: print ('sstart: '+str(dat[0][0])+', send: ' + str(dat[ndat-1][0])+', scen: ' + str(cen))
        rlast = 0.0
        dxlast = 0.0
        dxx = 0.
        for i in range(1,int(dxmax/0.1)):
            dx = float(i)*0.1
            sumt = 0.0
            for idat in range(ndat-1):
                sig = 0.5*(dat[idat][1]+dat[idat+1][1])*math.fabs(dat[idat][0]-dat[idat+1][0])
                pos = 0.5*(dat[idat+1][0]+dat[idat][0])

                if cen-1*dx < pos < cen+1*dx:
                    sumt += sig
            r = sumt/sum0

            #print str(dx)+'\t'+str(r)
            if rlast < ratio and r > ratio:
                #print r
                c1 = (dx-dxlast)/(r-rlast)
                c0 = dx-c1*r
                dxx = c1*ratio+c0
                if self.debug: print (str(100*ratio)+'% : '+str(dx)+'\t'+str(r)+'\t'+str(rlast)+'\t'+str(dxx))
                break
            
            rlast = r
            dxlast = dx

        return dxx
        
    #===
    def get_wire_sum_center_rms(self,epos0,esig0,mid1,mid2,icoord,par_dict,plot_dict):
        print('mid1: '+str(mid1)+'\tmid2: '+str(mid2))
        coord = par_dict['coord'][icoord+1]
        off0 = par_dict['offset0'][icoord]
        off1 = par_dict['offset1'][icoord]
        d = epos0[2].split(".")[0] #= date
        s = 0.0 #= sum
        r = 0.0 #= rms
        r90p = 0.0 #= 90% beam size
        r95p = 0.0 #= 95% beam size
        r99p = 0.0 #= 99% beam size
        c0 = 0.0 #= wire center
        c1 = 0.0 #= c0+off0
        c2 = 0.0 #= c0+off0-off1

        epos,esig = [],[]

        print("epos0:\t"+str(len(epos0)))
        for i in range(min([len(epos0[4:]),len(esig0[4:])])-nave+1):
            post,sigt = 0.0, 0.0
            for j in range(nave):
                post += float(epos0[4+i+j])
                sigt += float(esig0[4+i+j])
            epos.append(post/nave)
            esig.append(sigt/nave)
        
        #= Fill data to list and apply position offsets
        dat0,dat1 = [],[] #= dat1: background subtracted data
        sigl,pos0l,pos1l = [],[],[]
        for i,(post,sigt) in enumerate(zip(epos[4:],esig[4:])):
            if post == "0":
                continue
            pos,sig = float(post),float(sigt)
            sigl.append(sig)
            pos0l.append(pos)

            if par_dict['coord'][0] == 'L':
                if icoord == 0:
                    ang = par_dict['ang1']
                else:
                    ang = par_dict['ang2']
                    
                if coord == 'x':
                    pos1l.append((pos+off0-off1)*math.cos(ang*deg))
                elif coord == 'y':
                    pos1l.append((pos+off0-off1)*math.sin(ang*deg))
                elif coord == 'u':
                    pos1l.append((pos+off0-off1)*math.cos((ang-45.)*deg))
                elif coord == 'v':
                    pos1l.append((pos+off0-off1)*math.cos((ang-135.)*deg))
                else:
                    pos1l.append(pos+off0-off1)
            elif par_dict['coord'][0] == 'F':
                if coord == 'x':
                    pos1l.append(-1*(pos+off0-off1))
                elif coord == 'y':
                    pos1l.append(pos+off0-off1)
            elif par_dict['coord'][0] == 'S':
                if (coord == 'x'):
                    pos1l.append((pos+off0-off1)*math.cos(par_dict['ang1']*deg))
                elif (coord == 'y'):
                    pos1l.append((pos+off0-off1)*math.sin(par_dict['ang1']*deg))
                elif (coord == 'u'):
                    pos1l.append((pos+off0-off1)*math.cos((par_dict['ang1']-45.)*deg))                    
                elif (coord == 'v'):
                    pos1l.append((pos+off0-off1)*math.cos((par_dict['ang1']-135.)*deg))                    

            dat0.append([pos,sig,pos1l[i]])
                
                        

        dat0.sort()
        ndat0 = len(dat0)

        ran = self.get_range(dat0, mid1, mid2, par_dict)

        #= constant background calculation
        bgcoefs, bgstdv = self.calc_background(dat0, ran)
        ncoef = len(bgcoefs)
        #print str(bgcoefs)+'\t'+str(bgstdv)

        #= sum, center
        sigcorl, poscorl = [],[]
        for i in range(ndat0-1):
            pos1 = dat0[i+1][0]
            pos2 = dat0[i][0]
            sig1 = dat0[i+1][1]
            sig2 = dat0[i][1]
            for j in range(ncoef): 
                sig1 -= bgcoefs[j]*math.pow(pos1, ncoef-1-j)
                sig2 -= bgcoefs[j]*math.pow(pos2, ncoef-1-j)
            #if pos1 < ran[0] or pos2 < ran[0] or pos2 > ran[1] or pos1 > ran[1]: # ~ 2021/2/8
            #if (pos1 < ran[0]) or (pos2 < ran[0]) or (pos2 > ran[1]) or (pos1 > ran[1]) or (sig1 < bgth1*bgstdv): # 2021/2/8 ~ 2023/2/7
            if pos1 < ran[0] or pos2 < ran[0] or pos2 > ran[1] or pos1 > ran[1]: # ~ 2021/2/8
                continue
            poscorl.append(dat0[i][2])
            sigcorl.append(sig2)
            if sig1 < bgth1*bgstdv:
                continue
            sig = (sig1+sig2)/2.
            pos = (pos1+pos2)/2.
            dpos= math.fabs(pos2-pos1)
            s += sig*dpos
            c0 += sig*pos*dpos
            dat1.append([pos,sig])
        #print str(dat1)s
        if s>0.0:
            c0/=s
        if self.debug: print ("sum: "+str(s)+", center: "+str(c0))
        #= rms
        for i in range(ndat0-1):
            pos1 = dat0[i+1][0]
            pos2 = dat0[i][0]
            sig1 = dat0[i+1][1]
            sig2 = dat0[i][1]
            for j in range(ncoef):
                sig1 -= bgcoefs[j]*math.pow(pos1, ncoef-1-j)
                sig2 -= bgcoefs[j]*math.pow(pos2, ncoef-1-j)
            if (pos1 < ran[0]) or (pos2 < ran[0]) or (pos2 > ran[1]) or (pos1 > ran[1]) or (sig1 < bgth1*bgstdv):
            #if (pos1 < ran[0]) or (pos1 > ran[1]) or (sig1 < bgth1*bgstdv):
                continue
            sig = (sig1+sig2)/2.
            pos = (pos1+pos2)/2.
            dpos= math.fabs(pos2-pos1)
        
            r += math.pow(pos-c0,2)*sig*dpos
            #if self.debug: print str(pos1-c)+"\t"+str(sig1)+"\t"+str(r)
        if s>0.0:
            r = math.sqrt(r/s)
        else:
            r = math.sqrt(r)

        #= Subtraction of wire resolution from calculated rms
        r = np.sqrt(np.power(r,2)+np.power(par_dict['wire_diam'], 2))
            
        #= Apply center offset
        if par_dict['coord'][0] != 'S' and coord == "x":
            c1 = -1*(c0+off0)
            c = -1*(c0+off0-off1)
        else:
            c1 = c0+off0
            c = c0+off0-off1

        #= 90%,95% and 99% beam size
        try:
            r90p = self.calc_percent_beamsize(dat1, s, c0, 0.9)
        except:
            r90p = -999.
        try:
            r95p = self.calc_percent_beamsize(dat1, s, c0, 0.95)
        except:
            r95p = -999.
        try:
            r99p = self.calc_percent_beamsize(dat1, s, c0, 0.99)
        except:
            r99p = -999.
        #print str(r90p)+'\t'+str(r99p)
        
        #= fill to histgram
        plot_dict["h1a"+str(icoord)+"_x"],plot_dict["h1a"+str(icoord)+"_y"] = pos0l,sigl
        plot_dict["h1b"+str(icoord)+"_x"],plot_dict["h1b"+str(icoord)+"_y"] = pos1l,sigl

        fac = 15.
        plot_dict['h2'+str(icoord)+'_x0'],plot_dict['h2'+str(icoord)+'_y0'] = pos0l, sigl
        plot_dict['h2'+str(icoord)+'_x1'],plot_dict['h2'+str(icoord)+'_y1'] = [ran[0], ran[0]], [bgcoefs[len(bgcoefs)-1]-fac*bgstdv*0.9, bgcoefs[len(bgcoefs)-1]+fac*bgstdv*0.9]
        plot_dict['h2'+str(icoord)+'_x2'],plot_dict['h2'+str(icoord)+'_y2'] = [ran[1], ran[1]], [bgcoefs[len(bgcoefs)-1]-fac*bgstdv*0.9, bgcoefs[len(bgcoefs)-1]+fac*bgstdv*0.9]
        xs = np.arange(ran[0], ran[1])
        ys3 = []
        ys4 = []
        for xt in xs:
            yt = 0.
            for j,coef in enumerate(bgcoefs):
                yt += bgcoefs[j]*math.pow(xt, ncoef-1-j)
            ys3.append(yt)
            ys4.append(yt+bgth1*bgstdv)
        plot_dict['h2'+str(icoord)+'_x3'],plot_dict['h2'+str(icoord)+'_y3'] = xs, ys3
        plot_dict['h2'+str(icoord)+'_x4'],plot_dict['h2'+str(icoord)+'_y4'] = xs, ys4
        plot_dict['h2'+str(icoord)+'_set_ylim'] = [min(-2e-3, bgcoefs[len(bgcoefs)-1]-fac*bgstdv), max(2e-3, bgcoefs[len(bgcoefs)-1]+fac*bgstdv)]
        if self.debug:
            print('ys3: '+str(ys3[0])+'\tys4: '+str(ys4[0]))
        #= Fill results to wired
        wired = {}
        wired['d'] = d
        wired['s'] = s
        wired['r'] = r
        wired['c'] = c
        wired['c0'] = c0
        wired['c1'] = c1
        wired['r90p'] = r90p
        wired['r95p'] = r95p
        wired['r99p'] = r99p
        wired['sigl'] = sigcorl
        wired['posl'] = poscorl
        #return d,s,c,r
        return wired
    
    #===
    def get_xyuv_center(self,c6in,c12in1,c12in2,par_dict): #+ Right-hand system
        coord = par_dict['coord']
        offset = par_dict['offset1']
        u,v,x,y = -999.,-999.,-999.,-999.

        if coord[0] == 'F':
            x=c6in
            y=c12in1
            
        elif coord == 'Sbxy':
            x = (c12in1)*math.cos(par_dict['ang1']*deg)
            y = (c12in2)*math.sin(par_dict['ang1']*deg)
        #####
        if self.debug:
            f1 = open("debugout.out","a")
        if coord == "Luvx":
            u=c6in
            v=c12in1
            x1 = (u-v)/sq2
            x2 = (c12in2)/sq2
            x=(x1+x2)/2.
            y=(u+v)/sq2
            if self.debug:
                u1 = c6in
                u2 = sq2*x2+v
                print ("x1: "+str(x1)+"\tx2:"+str(x2))
                #f1.write("x1: "+str(x1)+"\tx2:"+str(x2)+"\n")
                print ("u1: "+str(u1)+"\tu2:"+str(u2))
                f1.write("u1: "+str(u1)+"\tu2:"+str(u2)+"\n")
        elif coord == "Luvy":
            u=c6in
            v=c12in1
            x=(u-v)/sq2
            y1 = (u+v)/sq2
            y2 = (c12in2)/sq2
            y = (y1+y2)/2.
            if self.debug:
                u1 = c6in
                u2 = sq2*y2-v
                #print ("y1: "+str(y1)+"\ty2: "+str(y2))
                #f1.write("y1: "+str(y1)+"\ty2: "+str(y2)+"\n")
                print ("u1: "+str(u1)+"\tu2:"+str(u2))
                f1.write("u1: "+str(u1)+"\tu2:"+str(u2)+"\n")
        elif coord == "Lyxu": # FS1 bending section
            y = c6in*math.sin(par_dict['ang1']*deg)
            x = c12in1*math.cos(par_dict['ang2']*deg)
            u1 = c12in2*math.cos((par_dict['ang2']+45.)*deg)
            u2 = (x+y)/sq2
            v = (y-x)/sq2
            u = (u1+u2)/2.
            if self.debug:
                print ("u1: "+str(u1)+"\tu2:"+str(u2))
                f1.write("u1: "+str(u1)+"\tu2:"+str(u2)+"\n")                
        elif coord == "Lyxv": # FS2 bending section and BDS dump
            y = c6in*math.sin(par_dict['ang1']*deg)
            x = c12in1*math.cos(par_dict['ang2']*deg)
            v1 = c12in2*math.cos((par_dict['ang2']-45.)*deg)
            v2 = (y-x)/sq2
            v = (v1+v2)/2.
            u = (x+y)/sq2
            if self.debug:
                print ("v1: "+str(v1)+"\tv2:"+str(v2))
                f1.write("v1: "+str(v1)+"\tv2:"+str(v2)+"\n")                
        elif coord == 'Suxy':
            x = (c12in1)*math.cos(par_dict['ang1']*deg)
            y = (c12in2)*math.sin(par_dict['ang1']*deg)
            u = c6in
            v = u-math.cos(par_dict['ang1']*deg)*x
            if self.debug:
                print ("x1: "+str(x)+"\tx2:"+str(-y*math.tan(par_dict['ang1']*deg)+u/math.cos(par_dict['ang1']*deg)))
                #f1.write("x1: "+str(x)+"\tx2:"+str(-1*y+sq2*u)+"\n")
                f1.write("x1: "+str(x)+"\tx2:"+str(-y*math.tan(par_dict['ang1']*deg)+u/math.cos(par_dict['ang1']*deg))+"\n")
        elif coord == 'Svxy':
            x = (c12in1)*math.cos(par_dict['ang1']*deg)
            y = (c12in2)*math.sin(par_dict['ang1']*deg)
            v = c6in
            u = v-math.cos(par_dict['ang1']*deg)*x
            if self.debug:
                print ("x1: "+str(x)+"\tx2:"+str(-y*math.tan(par_dict['ang1']*deg)+u/math.cos(par_dict['ang1']*deg)))
                #f1.write("x1: "+str(x)+"\tx2:"+str(-1*y+sq2*u)+"\n")
                f1.write("x1: "+str(x)+"\tx2:"+str(-y*math.tan(par_dict['ang1']*deg)+u/math.cos(par_dict['ang1']*deg))+"\n")
        elif coord == 'Svyx':
            x = (c12in2)*math.cos(par_dict['ang1']*deg)
            y = (c12in1)*math.sin(par_dict['ang1']*deg)
            v = c6in
            u = v-math.cos(par_dict['ang1']*deg)*x
            if self.debug:
                print ("x1: "+str(x)+"\tx2:"+str(-y*math.tan(par_dict['ang1']*deg)+u/math.cos(par_dict['ang1']*deg)))
                #f1.write("x1: "+str(x)+"\tx2:"+str(-1*y+sq2*u)+"\n")
                f1.write("x1: "+str(x)+"\tx2:"+str(-y*math.tan(par_dict['ang1']*deg)+u/math.cos(par_dict['ang1']*deg))+"\n")
        elif coord == 'Suyx': #= Need to check 2021/12/7
            x = (c12in2)*math.cos(par_dict['ang1']*deg)
            y = (c12in1)*math.sin(par_dict['ang1']*deg)
            u = c6in
            v = u-math.cos(par_dict['ang1']*deg)*x
            if self.debug:
                print ("x1: "+str(x)+"\tx2:"+str(-y*math.tan(par_dict['ang1']*deg)+u/math.cos(par_dict['ang1']*deg)))
                #f1.write("x1: "+str(x)+"\tx2:"+str(-1*y+sq2*u)+"\n")
                f1.write("x1: "+str(x)+"\tx2:"+str(-y*math.tan(par_dict['ang1']*deg)+u/math.cos(par_dict['ang1']*deg))+"\n")
        elif coord == 'Sxyu': #= Need to check 2021/12/7
            x = c6in*math.cos(par_dict['ang1']*deg)
            y = (c12in1)*math.sin(par_dict['ang1']*deg)
            u = (c12in2)
            v = u-math.cos(par_dict['ang1']*deg)*x
            if self.debug:
                print ("x1: "+str(x)+"\tx2:"+str(-y*math.tan(par_dict['ang1']*deg)+u/math.cos(par_dict['ang1']*deg)))
                #f1.write("x1: "+str(x)+"\tx2:"+str(-1*y+sq2*u)+"\n")
                f1.write("x1: "+str(x)+"\tx2:"+str(-y*math.tan(par_dict['ang1']*deg)+u/math.cos(par_dict['ang1']*deg))+"\n")
        elif coord == 'Sbxy':
            x = (c12in1)*math.cos(par_dict['ang1']*deg)
            y = (c12in2)*math.sin(par_dict['ang1']*deg)
        elif coord == "Fxy":
            '''
            x=-1*(c6in-offset[0])
            y=c12in1-offset[1]
            '''
            x=c6in
            y=c12in1
        
        return  x,y,u,v
        
    #===
    def get_xyuv_size(self,r6in,r12in1,r12in2,par_dict):
        coord = par_dict['coord']
        #print coord
        #print str(r6in)+'\t'+str(r12in1)+'\t'+str(r12in2)
        u,v,x,y = -999.,-999.,-999.,-999.
        if coord == "Luvx":
            u = r6in
            v = r12in1
            x = r12in2/sq2
            if u*u+v*v-x*x>0.:
                y = math.sqrt(u*u+v*v-x*x)
            else:
                y = -999.
        elif coord == "Luvy":
            u = r6in
            v = r12in1
            y = r12in2/sq2
            if u*u+v*v-y*y > 0.:
                x = math.sqrt(u*u+v*v-y*y)
            else:
                x = -999.
        elif coord == 'Lyxu':
            y = r6in*math.sin(par_dict['ang1']*deg)
            x = r12in1*math.cos(par_dict['ang2']*deg)
            u = r12in2*math.cos((par_dict['ang2']+45.)*deg)
            if x*x+y*y-u*u > 0.:
                v = math.sqrt(x*x+y*y-u*u)
            else:
                v = -999.
            
        elif coord == 'Lyxv':
            y = r6in*math.fabs(math.sin(par_dict['ang1']*deg))
            x = r12in1*math.fabs(math.cos(par_dict['ang2']*deg))
            v = r12in2*math.fabs(math.cos((par_dict['ang2']-45.)*deg))
            if x*x+y*y-v*v > 0.:
                u = math.sqrt(x*x+y*y-v*v)
            else:
                u = -999.
            
        elif coord == 'Suxy':
            u = r6in
            x = (r12in1)*math.cos(par_dict['ang1']*deg)
            y = (r12in2)*math.sin(par_dict['ang1']*deg)
            if x*x+y*y-u*u > 0.:
                v = math.sqrt(x*x+y*y-u*u)
            else:
                v = -999.
        elif coord == 'Suyx':
            u = r6in
            x = (r12in2)*math.cos(par_dict['ang1']*deg)
            y = (r12in1)*math.sin(par_dict['ang1']*deg)
            if x*x+y*y-u*u > 0.:
                v = math.sqrt(x*x+y*y-u*u)
            else:
                v = -999.
        elif coord == 'Sxyu':
            x = (r6in)*math.cos(par_dict['ang1']*deg)
            y = (r12in1)*math.sin(par_dict['ang1']*deg)
            u = r12in2
            if x*x+y*y-u*u > 0.:
                v = math.sqrt(x*x+y*y-u*u)
            else:
                v = -999.
        elif coord == 'Svxy':
            v = r6in
            x = (r12in1)*math.fabs(math.cos(par_dict['ang1']*deg))
            y = (r12in2)*math.fabs(math.sin(par_dict['ang1']*deg))
            if x*x+y*y-v*v > 0.:
                u = math.sqrt(x*x+y*y-v*v)
            else:
                u = -999.
        elif coord == 'Svyx':
            v = r6in
            x = (r12in2)*math.fabs(math.cos(par_dict['ang1']*deg))
            y = (r12in1)*math.fabs(math.sin(par_dict['ang1']*deg))
            if x*x+y*y-v*v > 0.:
                u = math.sqrt(x*x+y*y-v*v)
            else:
                u = -999.
            
        elif coord == 'Sbxy':
            x = (r12in1)*math.cos(par_dict['ang1']*deg)
            y = (r12in2)*math.sin(par_dict['ang1']*deg)
            
        elif coord == "Fxy":
            x = r6in
            y = r12in1
            u = -999.
            v = -999.
        return  x,y,u,v
    #===
    def get_cxy(self,xrms,yrms,urms,vrms):
        if min(xrms,yrms,urms,vrms)<-900.:
            cxy = -999.
        else:
            cxy = -(vrms**2-urms**2)/(2.0*xrms*yrms)
        return cxy
    #===
    def Pr_data_G(self, h1, coord="Luvy", plot_dict={}, Circleflg = 1, Saveflg=1, Showflg=0, savename = "Pr1.png", TitleData = ""):

        xc = plot_dict['xc']
        yc = plot_dict['yc']
        rms1 = plot_dict['r11']
        rms2 = plot_dict['r21']
        rms3 = plot_dict['r22']
        r90p1 = plot_dict['r90p11']
        r90p2 = plot_dict['r90p21']
        r90p3 = plot_dict['r90p22']

        urms = plot_dict['ur']
        vrms = plot_dict['vr']
        yrms = plot_dict['yr']
        xrms = plot_dict['xr']
        cxy = plot_dict['cxy']

        u90p = plot_dict['u90p']
        v90p = plot_dict['v90p']
        y90p = plot_dict['y90p']
        x90p = plot_dict['x90p']
        cxy90p = plot_dict['cxy90p']

        u99p = plot_dict['u99p']
        v99p = plot_dict['v99p']
        y99p = plot_dict['y99p']
        x99p = plot_dict['x99p']
        cxy99p = plot_dict['cxy99p']

        # urms -> rms1, vrms -> rms2, wrms -> rms3
        x = np.linspace(-1*plot_dict['h3_set_xylim'], plot_dict['h3_set_xylim'], 400)
        y = np.linspace(-1*plot_dict['h3_set_xylim'], plot_dict['h3_set_xylim'], 400)
        x, y = np.meshgrid(x, y)

        def axes():
            h1.axhline(0, alpha=.5)
            h1.axvline(0, alpha=.5)

        s22 = yrms * yrms
        s00 = xrms * xrms
        s02 = cxy*xrms*yrms
        s20 = s02
      
        s2290p = y90p * y90p
        s0090p = x90p * x90p
        s0290p = cxy90p*x90p*y90p
        s2090p = s0290p
      
      
        s2299p = y99p * y99p
        s0099p = x99p * x99p
        s0299p = cxy99p*x99p*y99p
        s2099p = s0299p
      
        if (coord == 'Sbxy'):
            h1.text(plot_dict['h3_set_xylim']-2, plot_dict['h3_set_xylim']-4*plot_dict['h3_set_xylim']/40., '('+coord[2]+','+coord[3]+') rms = ('+str(round(rms2,3))+', '+str(round(rms3,3))+') [mm]',fontsize=18)
        elif (coord == 'Fxy'):
            h1.text(plot_dict['h3_set_xylim']-2, plot_dict['h3_set_xylim']-4*plot_dict['h3_set_xylim']/40., '('+coord[1]+','+coord[2]+') rms = ('+str(round(rms1,3))+', '+str(round(rms2,3))+') [mm]',fontsize=18)
        else:
            h1.text(plot_dict['h3_set_xylim']-2, plot_dict['h3_set_xylim']-4*plot_dict['h3_set_xylim']/40., '('+coord[1]+','+coord[2]+','+coord[3]+') rms = ('+str(round(rms1,3))+', '+str(round(rms2,3))+', '+str(round(rms3,3))+') [mm]',fontsize=18)
        h1.text(plot_dict['h3_set_xylim']-2, plot_dict['h3_set_xylim']-10*plot_dict['h3_set_xylim']/40.,'(x,y) center = ('+str(round(xc,2))+', ' + str(round(yc,2))  +  ') [mm]', fontsize=18)
        h1.text(plot_dict['h3_set_xylim']-2, plot_dict['h3_set_xylim']-16*plot_dict['h3_set_xylim']/40.,'(x,y) rms, cxy = (' + str(round(xrms,2)) + ', ' + str(round(yrms,2))  + ') [mm], '+ str(round(cxy,3)), fontsize=18)

        h1.set_xlabel('X [mm]', fontsize=18)
        h1.set_ylabel('Y [mm]', fontsize=18)

        h1.axis('equal')
        h1.set_xlim(plot_dict['h3_set_xylim'], -1*plot_dict['h3_set_xylim'])
        h1.set_ylim(-1*plot_dict['h3_set_xylim'], plot_dict['h3_set_xylim'])
        h1.grid()

        if 'u' in coord:
            h1.plot([xc+urms/np.sqrt(2)-100, xc+urms/np.sqrt(2)+100],[yc+urms/np.sqrt(2)+100, yc+urms/np.sqrt(2)-100],'-',color=self.color['u'], linewidth=1)
            h1.plot([xc-urms/np.sqrt(2)+100, xc-urms/np.sqrt(2)-100],[yc-urms/np.sqrt(2)-100, yc-urms/np.sqrt(2)+100],'-',color=self.color['u'], linewidth=1)
        if 'v' in coord:
            h1.plot([xc+vrms/np.sqrt(2)-100, xc+vrms/np.sqrt(2)+100],[yc-vrms/np.sqrt(2)-100, yc-vrms/np.sqrt(2)+100],'-',color=self.color['v'], linewidth=1)
            h1.plot([xc-vrms/np.sqrt(2)-100, xc-vrms/np.sqrt(2)+100],[yc+vrms/np.sqrt(2)-100, yc+vrms/np.sqrt(2)+100],'-',color=self.color['v'], linewidth=1)
        if 'x' in coord:
            h1.axvline(xc+xrms, alpha=1,color=self.color['x'], linewidth=1)
            h1.axvline(xc-xrms, alpha=1,color=self.color['x'], linewidth=1)
        if 'y' in coord:
            h1.axhline(yc+yrms, alpha=1,color=self.color['y'], linewidth=1)
            h1.axhline(yc-yrms, alpha=1,color=self.color['y'], linewidth=1)

        h1.contour(x, y,(s22*(x-xc)**2 - 2*s02*(x-xc)*(y-yc) + s00*(y-yc)**2), [s22*s00-s02*s20], colors='black')
        h1.contour(x, y,(s2290p*(x-xc)**2 - 2*s0290p*(x-xc)*(y-yc) + s0090p*(y-yc)**2), [s2290p*s0090p-s0290p*s2090p], colors='gray')
        h1.contour(x, y,(s2299p*(x-xc)**2 - 2*s0299p*(x-xc)*(y-yc) + s0099p*(y-yc)**2), [s2299p*s0099p-s0299p*s2099p], colors='silver')

        #art = []
        #lgd=h1.legend(loc=9, bbox_to_anchor=(0.5, -0.15), ncol=3,fontsize=18)
        #art.append(lgd)
        #h1.legend(loc=9, bbox_to_anchor=(0.5, -0.15), ncol=3, fontsize=18)
        '''
        if(Circleflg == 1):
            h1.contour(x, y,((x-xc)**2 +(y-yc)**2), [rms1**2], color='gray', linewidths = 0.3)
            h1.contour(x, y,((x-xc)**2 +(y-yc)**2), [rms2**2], color='gray', linewidths = 0.3)
        '''
        return_data = [xrms, yrms, cxy]
        return return_data
  
    #==
    '''
    def Pr_data_G_old(self, h1,coord="Luvy",urms=7.3,vrms=17.1,wrms=10.6, xc=0, yc=0, xyrange = [-20, 20, -20, 20], Circleflg = 1, Saveflg=1, Showflg=0, savename = "Pr1.png", TitleData = ""):
        x = np.linspace(xyrange[0], xyrange[1], 400)
        y = np.linspace(xyrange[2], xyrange[3], 400)
        x, y = np.meshgrid(x, y)

        def axes():
            h1.axhline(0, alpha=.5)
            h1.axvline(0, alpha=.5)

        #a, b, c, d, e, f = 4, -5, 2, 4, -3, -3
        #assert b**2 - 4*a*c < 0
      
        if (coord == "Luvy"):
            xrms = np.sqrt(urms*urms+vrms*vrms-wrms*wrms)
            yrms = wrms 
            cxy = -(vrms**2-urms**2)/(2.0*xrms*yrms)
        elif (coord == "Luvx"):
            xrms = wrms
            yrms = np.sqrt(urms*urms+vrms*vrms-wrms*wrms)
            cxy = -(vrms**2-urms**2)/(2.0*xrms*yrms)
        elif (coord == "Suyx"):
            xrms = wrms
            yrms = np.sqrt(urms*urms+vrms*vrms-wrms*wrms)
            cxy = -(vrms**2-urms**2)/(2.0*xrms*yrms)
        ###############
        elif (coord == "Sxyu"):
            xrms = wrms
            yrms = np.sqrt(urms*urms+vrms*vrms-wrms*wrms)
            cxy = -(vrms**2-urms**2)/(2.0*xrms*yrms)
        ###############

        s22 = yrms * yrms
        s00 = xrms * xrms
        s02 = cxy*xrms*yrms
        s20 = s02
      
        
        h1.text(xyrange[0]+2,xyrange[3]-4,"(u,v,w)rms = ("+str(round(urms,3))+", "+str(round(vrms,3))+", "+str(round(wrms,3))+") [mm]",fontsize=18)
        h1.text(xyrange[0]+2,xyrange[3]-10,"(x,y)center = ("+str(round(xc,2))+", " + str(round(yc,2))  +  ") [mm]", fontsize=18)
        h1.text(xyrange[0]+2,xyrange[3]-16,"(x,y)rms, cxy = (" + str(round(xrms,2)) + ", " + str(round(yrms,2))  + ") [mm], "+ str(round(cxy,3)), fontsize=18)
        h1.set_xlabel('X [mm]', fontsize=18)
        h1.set_ylabel('Y [mm]', fontsize=18)

        h1.axis('equal')
        h1.axis(xyrange)
        h1.grid()
        # urms limit line
        h1.plot([xc + urms/np.sqrt(2)-100,xc +  urms/np.sqrt(2)+100],[yc + urms/np.sqrt(2)+100, yc + urms/np.sqrt(2)-100],'-',color=self.color[coord[1]], label='u-wire')
        h1.plot([xc- urms/np.sqrt(2)+100,xc - urms/np.sqrt(2)-100],[yc- urms/np.sqrt(2)-100,yc - urms/np.sqrt(2)+100],'-',color=self.color[coord[1]])
        # vrms limit line
        h1.plot([xc + vrms/np.sqrt(2)-100, xc + vrms/np.sqrt(2)+100],[yc - vrms/np.sqrt(2)-100, yc - vrms/np.sqrt(2)+100],'-',color=self.color[coord[2]], label='v-wire')
        h1.plot([xc - vrms/np.sqrt(2)-100, xc - vrms/np.sqrt(2)+100],[yc + vrms/np.sqrt(2)-100, yc + vrms/np.sqrt(2)+100],'-',color=self.color[coord[2]])
     

        if (coord == "Luvy"):
            h1.axhline(yc + wrms, alpha=1,color=self.color[coord[3]], label='y-wire')
            h1.axhline(yc - wrms, alpha=1,color=self.color[coord[3]])
        elif (coord == "Luvx"):
            # wrms limit line
            h1.axvline(xc + wrms, alpha=1,color=self.color[coord[3]], label='x-wire')
            h1.axvline(xc - wrms, alpha=1,color=self.color[coord[3]])

        h1.contour(x, y,(s22*(x-xc)**2 - 2*s02*(x-xc)*(y-yc) + s00*(y-yc)**2), [s22*s00-s02*s20], colors='k')

        art = []
        lgd=h1.legend(loc=9, bbox_to_anchor=(0.5, -0.15), ncol=3,fontsize=18)
        art.append(lgd)



        if(Circleflg == 1):
            h1.contour(x, y,((x-xc)**2 +(y-yc)**2), [urms**2], color='gray', linewidths = 0.3)
            h1.contour(x, y,((x-xc)**2 +(y-yc)**2), [vrms**2], color='gray', linewidths = 0.3)

        return_data = [xrms, yrms, cxy]
        return return_data
  
    '''
    #==
    def plot(self,par_dict, plot_dict):
        fig = plt.figure(2,facecolor="white",figsize=(15,14))
        #fig = plt.figure(facecolor="white",figsize=(15,14))
        
        h1a = fig.add_subplot(4,2,1) # distribution, raw_position
        h1b = fig.add_subplot(4,2,3) # distribution
        h2a = fig.add_subplot(6,2,2) # noise plot, 1st wire
        h2b = fig.add_subplot(6,2,4) # noise plot, 2nd wire
        h2c = fig.add_subplot(6,2,6) # noise plot, 3rd wire
        h3 = fig.add_subplot(2,2,3)  # x-y plot
        h4 = fig.add_subplot(2,2,4)  # file name
        h5 = fig.add_subplot(7,2,10) # upper table
        h6 = fig.add_subplot(8,2,14) # lower table
        '''

        #plt.rcParam['figure.figsize'] = [15,14]
        h1a = plt.subplot(4,2,1)
        h1b = plt.subplot(4,2,3)
        h2a = plt.subplot(6,2,2)
        h2b = plt.subplot(6,2,4)
        h2c = plt.subplot(6,2,6)
        h3 = plt.subplot(2,2,3)
        h4 = plt.subplot(2,2,4)
        h5 = plt.subplot(6,2,10)
        h6 = plt.subplot(6,2,12)
        '''

        pltxmin,pltxmax,pltymin,pltymax = 0.0,0.0,0.0,0.0
        h1a.cla()
        h1b.cla()
        h2a.cla()
        h2b.cla()
        h2c.cla()
        h3.cla()
        h4.cla()
        h5.cla()
        h6.cla()
        h4.axis("off")
        h5.axis("off")
        h6.axis("off")
        
        h1a.set_xlabel("Wire position [mm]")
        h1a.set_ylabel("Signal")
        h1a.grid()
        h1a.ticklabel_format(style="sci",axis="y", scilimits=(0,0))
        h1b.set_xlabel("Beam position [mm]")
        h1b.set_ylabel("Signal")
        h1b.grid()
        h1b.ticklabel_format(style="sci",axis="y", scilimits=(0,0))
        h1b.set_xlim(-1*plot_dict['h1b_set_xlim'], plot_dict['h1b_set_xlim'])
        #h2a.set_yscale("log")
        h2a.ticklabel_format(style="sci",axis="y", scilimits=(0,0))
        h2b.ticklabel_format(style="sci",axis="y", scilimits=(0,0))
        h2c.ticklabel_format(style="sci",axis="y", scilimits=(0,0))
        if par_dict['coord'] == "Fxy":
            for i in range(2):
                h1a.plot(plot_dict['h1a'+str(i)+'_x'],plot_dict['h1a'+str(i)+'_y'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][i+1]+'-wire')
                h1b.plot(plot_dict['h1b'+str(i)+'_x'],plot_dict['h1b'+str(i)+'_y'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][i+1]+'-wire')

            i = 0
            h2a.plot(plot_dict['h2'+str(i)+'_x0'],plot_dict['h2'+str(i)+'_y0'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][1]+'-wire',alpha = 0.5)
            h2a.plot(plot_dict['h2'+str(i)+'_x1'],plot_dict['h2'+str(i)+'_y1'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2a.plot(plot_dict['h2'+str(i)+'_x2'],plot_dict['h2'+str(i)+'_y2'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2a.plot(plot_dict['h2'+str(i)+'_x3'],plot_dict['h2'+str(i)+'_y3'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2a.plot(plot_dict['h2'+str(i)+'_x4'],plot_dict['h2'+str(i)+'_y4'],self.color[par_dict['coord'][i+1]],linewidth = 0.5,linestyle=':')
            h2a.set_ylim(plot_dict['h2'+str(i)+'_set_ylim'])
            i = 1
            h2b.plot(plot_dict['h2'+str(i)+'_x0'],plot_dict['h2'+str(i)+'_y0'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][1]+'-wire',alpha = 0.5)
            h2b.plot(plot_dict['h2'+str(i)+'_x1'],plot_dict['h2'+str(i)+'_y1'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2b.plot(plot_dict['h2'+str(i)+'_x2'],plot_dict['h2'+str(i)+'_y2'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2b.plot(plot_dict['h2'+str(i)+'_x3'],plot_dict['h2'+str(i)+'_y3'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2b.plot(plot_dict['h2'+str(i)+'_x4'],plot_dict['h2'+str(i)+'_y4'],self.color[par_dict['coord'][i+1]],linewidth = 0.5,linestyle=':')
            h2b.set_ylim(plot_dict['h2'+str(i)+'_set_ylim'])

            h1a.legend(ncol=2)
            h1a.set_xlim(-80,80)
            h2a.set_xlim(-80,80)
            h2b.set_xlim(-80,80)
            h2c.axis("off")
            h3.set_ylim(-40,40)
            h3.set_xlim(-40,40)
            #h3.axis("off")
        elif par_dict['coord'] == "Sbxy":
            for i in range(1,3):
                h1a.plot(plot_dict['h1a'+str(i)+'_x'],plot_dict['h1a'+str(i)+'_y'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][i+1]+'-wire')
                h1b.plot(plot_dict['h1b'+str(i)+'_x'],plot_dict['h1b'+str(i)+'_y'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][i+1]+'-wire')

            i = 1
            h2a.plot(plot_dict['h2'+str(i)+'_x0'],plot_dict['h2'+str(i)+'_y0'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][1]+'-wire',alpha = 0.5)
            h2a.plot(plot_dict['h2'+str(i)+'_x1'],plot_dict['h2'+str(i)+'_y1'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2a.plot(plot_dict['h2'+str(i)+'_x2'],plot_dict['h2'+str(i)+'_y2'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2a.plot(plot_dict['h2'+str(i)+'_x3'],plot_dict['h2'+str(i)+'_y3'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2a.plot(plot_dict['h2'+str(i)+'_x4'],plot_dict['h2'+str(i)+'_y4'],self.color[par_dict['coord'][i+1]],linewidth = 0.5,linestyle=':')
            h2a.set_ylim(plot_dict['h2'+str(i)+'_set_ylim'])
            i = 2
            h2b.plot(plot_dict['h2'+str(i)+'_x0'],plot_dict['h2'+str(i)+'_y0'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][1]+'-wire',alpha = 0.5)
            h2b.plot(plot_dict['h2'+str(i)+'_x1'],plot_dict['h2'+str(i)+'_y1'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2b.plot(plot_dict['h2'+str(i)+'_x2'],plot_dict['h2'+str(i)+'_y2'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2b.plot(plot_dict['h2'+str(i)+'_x3'],plot_dict['h2'+str(i)+'_y3'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2b.plot(plot_dict['h2'+str(i)+'_x4'],plot_dict['h2'+str(i)+'_y4'],self.color[par_dict['coord'][i+1]],linewidth = 0.5,linestyle=':')
            h2b.set_ylim(plot_dict['h2'+str(i)+'_set_ylim'])

            h1a.legend(ncol=2)
            h1a.set_xlim(-200,0)
            h2a.set_xlim(-200,0)
            h2b.set_xlim(-200,0)
            h2c.axis("off")

            h3.set_ylim(-40,40)
            h3.set_xlim(-40,40)
            #h3.axis("off")
        elif (('L' in par_dict['coord']) or ('S' in par_dict['coord'])) and len(par_dict['coord'])==4:
            for i in range(3):
                h1a.plot(plot_dict['h1a'+str(i)+'_x'],plot_dict['h1a'+str(i)+'_y'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][i+1]+'-wire')
                h1b.plot(plot_dict['h1b'+str(i)+'_x'],plot_dict['h1b'+str(i)+'_y'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][i+1]+'-wire')
            i = 0
            h2a.plot(plot_dict['h2'+str(i)+'_x0'],plot_dict['h2'+str(i)+'_y0'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][1]+'-wire',alpha = 0.5)
            h2a.plot(plot_dict['h2'+str(i)+'_x1'],plot_dict['h2'+str(i)+'_y1'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2a.plot(plot_dict['h2'+str(i)+'_x2'],plot_dict['h2'+str(i)+'_y2'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2a.plot(plot_dict['h2'+str(i)+'_x3'],plot_dict['h2'+str(i)+'_y3'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2a.plot(plot_dict['h2'+str(i)+'_x4'],plot_dict['h2'+str(i)+'_y4'],self.color[par_dict['coord'][i+1]],linewidth = 0.5,linestyle=':')
            h2a.set_ylim(plot_dict['h2'+str(i)+'_set_ylim'])
            i = 1
            h2b.plot(plot_dict['h2'+str(i)+'_x0'],plot_dict['h2'+str(i)+'_y0'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][1]+'-wire',alpha = 0.5)
            h2b.plot(plot_dict['h2'+str(i)+'_x1'],plot_dict['h2'+str(i)+'_y1'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2b.plot(plot_dict['h2'+str(i)+'_x2'],plot_dict['h2'+str(i)+'_y2'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2b.plot(plot_dict['h2'+str(i)+'_x3'],plot_dict['h2'+str(i)+'_y3'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2b.plot(plot_dict['h2'+str(i)+'_x4'],plot_dict['h2'+str(i)+'_y4'],self.color[par_dict['coord'][i+1]],linewidth = 0.5,linestyle=':')
            h2b.set_ylim(plot_dict['h2'+str(i)+'_set_ylim'])
            i = 2
            h2c.plot(plot_dict['h2'+str(i)+'_x0'],plot_dict['h2'+str(i)+'_y0'],self.color[par_dict['coord'][i+1]],label=par_dict['coord'][1]+'-wire',alpha = 0.5)
            h2c.plot(plot_dict['h2'+str(i)+'_x1'],plot_dict['h2'+str(i)+'_y1'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2c.plot(plot_dict['h2'+str(i)+'_x2'],plot_dict['h2'+str(i)+'_y2'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2c.plot(plot_dict['h2'+str(i)+'_x3'],plot_dict['h2'+str(i)+'_y3'],self.color[par_dict['coord'][i+1]],linewidth = 0.5)
            h2c.plot(plot_dict['h2'+str(i)+'_x4'],plot_dict['h2'+str(i)+'_y4'],self.color[par_dict['coord'][i+1]],linewidth = 0.5,linestyle=':')
            h2c.set_ylim(plot_dict['h2'+str(i)+'_set_ylim'])
            h3.set_ylim(-40,40)
            h3.set_xlim(-40,40)
        
        
            h1a.legend(ncol=3)#, bbox_to_anchor=(0.6, 1.05))
            h1a.set_xlim(min(plot_dict['h1a0_x'][len(plot_dict['h1a0_x'])-1],plot_dict['h1a1_x'][len(plot_dict['h1a1_x'])-1],plot_dict['h1a2_x'][len(plot_dict['h1a2_x'])-1])-10,0)
            h2a.set_xlim(min(plot_dict['h1a0_x'][len(plot_dict['h1a0_x'])-1],plot_dict['h1a1_x'][len(plot_dict['h1a1_x'])-1],plot_dict['h1a2_x'][len(plot_dict['h1a2_x'])-1])-10,0)
            h2b.set_xlim(min(plot_dict['h1a0_x'][len(plot_dict['h1a0_x'])-1],plot_dict['h1a1_x'][len(plot_dict['h1a1_x'])-1],plot_dict['h1a2_x'][len(plot_dict['h1a2_x'])-1])-10,0)
            h2c.set_xlim(min(plot_dict['h1a0_x'][len(plot_dict['h1a0_x'])-1],plot_dict['h1a1_x'][len(plot_dict['h1a1_x'])-1],plot_dict['h1a2_x'][len(plot_dict['h1a2_x'])-1])-10,0)
            
        h4.text(0.0, 0.95, par_dict['file_name'], fontsize=17, fontweight = 'bold')
        #row_label1 = ["Time","Sum","Center0","Center1","Center","RMS"]
        row_label1 = ["Time","Sum","Center","RMS"]
        row_label2 = ["Pos","RMS","cxy","R90%","cxy90%","R99%","cxy99%"]
        if (par_dict['coord']=="Fxy"):
            col_label1 = [par_dict['coord'][1] + "-wire",par_dict['coord'][2] + "-wire"]
            tab_value1 = [[plot_dict['d11'], plot_dict['d21']],
                          [plot_dict['s11'], plot_dict['s21']],
                          #[plot_dict['c011'], plot_dict['c021']],
                          #[plot_dict['c111'], plot_dict['c121']],
                          [plot_dict['c11'], plot_dict['c21']],
                          [plot_dict['r11'], plot_dict['r21']]
            ]
            col_label2 = ["x", "y"]
            tab_value2 = [[plot_dict['xc'], plot_dict['yc']],
                          [plot_dict['xr'], plot_dict['yr']],
                          ["", ""],
                          [plot_dict['x90p'], plot_dict['y90p']],
                          ["", ""],
                          [plot_dict['x99p'], plot_dict['y99p']],
                          ["", ""]
            ]
        elif (par_dict['coord']=="Sbxy"):
            col_label1 = [par_dict['coord'][2] + "-wire",par_dict['coord'][3] + "-wire"]
            tab_value1 = [[plot_dict['d21'],''],
                          [plot_dict['s21'], plot_dict['s22']],
                          #[plot_dict['c021'], plot_dict['c022']],
                          #[plot_dict['c121'], plot_dict['c122']],
                          [plot_dict['c21'], plot_dict['c22']],
                          [plot_dict['r21'], plot_dict['r22']]
            ]
            col_label2 = ["x", "y"]
            tab_value2 = [[plot_dict['xc'], plot_dict['yc']],
                          [plot_dict['xr'], plot_dict['yr']],
                          ["", ""],
                          [plot_dict['x90p'], plot_dict['y90p']],
                          ["", ""],
                          [plot_dict['x99p'], plot_dict['y99p']],
                          ["", ""]
            ]
        else:
            col_label1 = [par_dict['coord'][1] + "-wire",par_dict['coord'][2] + "-wire",par_dict['coord'][3] + "-wire"]
            tab_value1 = [[plot_dict['d11'], plot_dict['d21'], ""],
                          [plot_dict['s11'], plot_dict['s21'],plot_dict['s22']],
                          #[plot_dict['c011'], plot_dict['c021'],plot_dict['c022']],
                          #[plot_dict['c111'], plot_dict['c121'],plot_dict['c122']],
                          [plot_dict['c11'], plot_dict['c21'],plot_dict['c22']],
                          [plot_dict['r11'], plot_dict['r21'],plot_dict['r22']]
            ]
            col_label2 = ["x", "y", "u", "v"]
            tab_value2 = [[plot_dict['xc'], plot_dict['yc'], plot_dict['uc'], plot_dict['vc']],
                          [plot_dict['xr'], plot_dict['yr'], plot_dict['ur'], plot_dict['vr']],
                          [plot_dict['cxy'],"","",""],
                          [plot_dict['x90p'], plot_dict['y90p'], plot_dict['u90p'], plot_dict['v90p']],
                          [plot_dict['cxy90p'],"","",""],
                          [plot_dict['x99p'], plot_dict['y99p'], plot_dict['u99p'], plot_dict['v99p']],
                          [plot_dict['cxy99p'],"","",""]
            ]

        h5tab1 = h5.table(cellText=tab_value1,
                          colWidths=[0.27]*len(col_label1),
                          rowLabels=row_label1,
                          colLabels=col_label1,
                          loc="center right")
        h5tab1.set_fontsize(18)
        h5tab1.scale(1,1.6)
        h6tab = h6.table(cellText=tab_value2,
                         colWidths = [0.22]*len(col_label2),
                         rowLabels = row_label2,
                         colLabels = col_label2,
                         loc = "center right")
        h6tab.set_fontsize(18)
        h6tab.scale(1,1.6)

        #if len(par_dict['coord']) != 3:
            #try:
        self.Pr_data_G(h3, par_dict['coord'], plot_dict, 1, 1, 1)

            #except:
            #    print 'No 2D plot displayed'

        if self.jupyter:
            fig.canvas.draw()
            plt.show()
        else:
            #fig.canvas.draw()
            if self.batch == False:
                fig.show()
        fn = par_dict['file_name'].split('/')[-1].replace('.dat','.png')
        fig.savefig(par_dict["file_dict"]+'/'+fn, dpi=(144), bbox_inches="tight")
        fig.clear()
        
    #===
    def execute(self,inp,form = 1,file_dict = ''):
        #print (type(inp))
        if isinstance(inp, str):
            par_dict = self.get_param(inp)
            par_dict['file_name'] = inp
            raw_dict = self.get_rawdata(inp,par_dict)
            if ("D3959" in inp) or ("D4222" in inp):
                tmp = raw_dict['a1s1']
                raw_dict['a1s1'] = raw_dict['a1s3']
                raw_dict['a1s3'] = tmp
        elif isinstance(inp, dict):
            par_dict = self.get_param(inp['file_name'].replace('.dat',''))
            par_dict['file_name'] = inp['file_name']
            raw_dict = inp['raw_data']
            ## wrong cable connection
            if ("D3959" in inp['file_name']) or ("D4222" in inp['file_name']):
                tmp = raw_dict['a1s1']
                raw_dict['a1s1'] = raw_dict['a1s3']
                raw_dict['a1s3'] = tmp

                    
        if len(file_dict)==0:
            tmp = par_dict['file_name'].split('/')
            del tmp[-1]
            if len(tmp) == 0:
                dic = './'
            else:
                dic = '/'.join(tmp)
            par_dict['file_dict']=dic
        else:
            par_dict['file_dict']=file_dict
            
        
        '''        
        if par_dict['live']:
            inp = par_dict['name'].replace(":","_")+"_"+datetime.now().strftime("%Y%m%d_%H%M%S")+".dat"
            f = open(inp,"w")
            for l in raw_dict.values():
                f.write(l+'\n')
            f.close()
            print "###########"
            print "file: "+par_dict["file_name"]
        '''
            
        
        if par_dict['coord'][0] == 'L':
            e1o1 = raw_dict.get('a1o1').split()
            e2o1 = raw_dict.get('a2o1').split()
            e2o2 = raw_dict.get('a2o2').split()
            par_dict['offset0'][0]=float(e1o1[3])
            par_dict['offset0'][1]=float(e2o1[3])
            par_dict['offset0'][2]=float(e2o2[3])
            plot_dict['h1b_set_xlim'] = 40.
            plot_dict['h3_set_xylim'] = 40.
        elif par_dict['coord'][0] == 'F':
            if 'a1o1' in raw_dict.keys():
                e1o1 = raw_dict.get('a1o1').split()
                par_dict['offset0'][0]=float(e1o1[3])
            if 'a2o1' in raw_dict.keys():
                e2o1 = raw_dict.get('a2o1').split()
                par_dict['offset0'][1]=float(e2o1[3])
            plot_dict['h1b_set_xlim'] = 40.
            plot_dict['h3_set_xylim'] = 40.
        elif par_dict['coord'] == 'Sbxy':
            e1o1 = raw_dict.get('a1o1').split()
            e1o2 = raw_dict.get('a1o2').split()
            e1o3 = raw_dict.get('a1o3').split()
            par_dict['offset0'][0]=float(e1o1[3])
            par_dict['offset0'][1]=float(e1o2[3])
            par_dict['offset0'][2]=float(e1o3[3])
            plot_dict['h1b_set_xlim'] = 40.
            plot_dict['h3_set_xylim'] = 40.
        elif par_dict['coord'][0] == 'S':
            e1o1 = raw_dict.get('a1o1').split()
            e1o2 = raw_dict.get('a1o2').split()
            e1o3 = raw_dict.get('a1o3').split()
            par_dict['offset0'][0]=float(e1o1[3])
            par_dict['offset0'][1]=float(e1o2[3])
            par_dict['offset0'][2]=float(e1o3[3])
            plot_dict['h1b_set_xlim'] = 20.
            plot_dict['h3_set_xylim'] = 20.

        xr,yr,ur,vr = -999., -999., -999., -999.
        xc,yc,uc,vc = -999., -999., -999., -999.
        
        #print par_dict
        if par_dict['coord'] == "Fxy":
            wired11 = self.get_wire_sum_center_rms(raw_dict['a1p'].split(),raw_dict['a1s1'].split(),-999.,999.,0,par_dict,plot_dict) #= x-wire #= date, signal sum, center, rms
            wired21 = self.get_wire_sum_center_rms(raw_dict['a2p'].split(),raw_dict['a2s1'].split(),-999.,999.,1,par_dict,plot_dict) #= y-wire
            wired22 = {'d': 0.0, 's': 0.0, 'c': 0.0, 'r': 0.0, 'c0': 0.0, 'c1': 0.0, 'r90p': 0.0, 'r95p' : 0.0, 'r99p' : 0.0}
            xc,yc,uc,vc = self.get_xyuv_center(wired11['c'],wired21['c'],wired22['c'],par_dict)
            xr,yr,ur,vr = self.get_xyuv_size(wired11['r'],wired21['r'],wired22['r'],par_dict)
            x90p,y90p,u90p,v90p = self.get_xyuv_size(wired11['r90p'],wired21['r90p'],wired22['r90p'],par_dict)
            x95p,y95p,u95p,v95p = self.get_xyuv_size(wired11['r95p'],wired21['r95p'],wired22['r95p'],par_dict)
            x99p,y99p,u99p,v99p = self.get_xyuv_size(wired11['r99p'],wired21['r99p'],wired22['r99p'],par_dict)
            cxy = -999.
            cxy90p = -999.
            cxy95p = -999.
            cxy99p = -999.

        elif par_dict['coord'] == "Luvy" or par_dict['coord'] == "Luvx" or par_dict['coord'] == 'Lyxu' or par_dict['coord'] == 'Lyxv':
            mid = self.get_middlepos(raw_dict['a2p'].split(), raw_dict['a2s1'].split(), raw_dict['a2s2'].split(), 1/sq2, par_dict['offset0'][1], par_dict['offset0'][2], par_dict['coord'])
            wired11 = self.get_wire_sum_center_rms(raw_dict['a1p'].split(),raw_dict['a1s1'].split(),-999.,999.,0,par_dict,plot_dict) #= u-wire #= date, signal sum, center, rms
            wired21 = self.get_wire_sum_center_rms(raw_dict['a2p'].split(),raw_dict['a2s1'].split(),mid-10,999.,1,par_dict,plot_dict) #= v-wire
            wired22 = self.get_wire_sum_center_rms(raw_dict['a2p'].split(),raw_dict['a2s2'].split(),-999.,mid+10,2,par_dict,plot_dict) #= y-wire
            xc,yc,uc,vc = self.get_xyuv_center(wired11['c'],wired21['c'],wired22['c'],par_dict)
            xr,yr,ur,vr = self.get_xyuv_size(wired11['r'],wired21['r'],wired22['r'],par_dict)
            x90p,y90p,u90p,v90p = self.get_xyuv_size(wired11['r90p'],wired21['r90p'],wired22['r90p'],par_dict)
            x95p,y95p,u95p,v95p = self.get_xyuv_size(wired11['r95p'],wired21['r95p'],wired22['r95p'],par_dict)
            x99p,y99p,u99p,v99p = self.get_xyuv_size(wired11['r99p'],wired21['r99p'],wired22['r99p'],par_dict)
            try:
                cxy = self.get_cxy(xr,yr,ur,vr)
            except:
                cxy = -999.
            try:
                cxy90p = self.get_cxy(x90p,y90p,u90p,v90p)
            except:
                cxy90p = -999.
            try:
                cxy95p = self.get_cxy(x95p,y95p,u95p,v95p)
            except:
                cxy95p = -999.
            try:
                cxy99p = self.get_cxy(x99p,y99p,u99p,v99p)
            except:
                cxy99p = -999.
        elif (par_dict['coord'] == "Suxy") or (par_dict['coord'] == "Svxy") or (par_dict['coord'] == "Svyx") or (par_dict['coord'] == "Suyx") or (par_dict['coord'] == "Sxyu"): #11: u-wire, 21: y-wire and 22: x-wire
            '''
            wired11 = self.get_wire_sum_center_rms(raw_dict['a1p'].split(), raw_dict['a1s1'].split(), par_dict['scan_stoppos1']-float(raw_dict['a1o1'].split()[3])-5, par_dict['scan_startpos1']-float(raw_dict['a1o1'].split()[3])+5, 0, par_dict, plot_dict) #= u-wire #= date, signal sum, center, rms
            wired21 = self.get_wire_sum_center_rms(raw_dict['a1p'].split(), raw_dict['a1s2'].split(), par_dict['scan_stoppos2']-float(raw_dict['a1o1'].split()[3])-5, par_dict['scan_startpos2']-float(raw_dict['a1o1'].split()[3])+5, 1, par_dict, plot_dict) #= x-wire
            wired22 = self.get_wire_sum_center_rms(raw_dict['a1p'].split(), raw_dict['a1s3'].split(), par_dict['scan_stoppos3']-float(raw_dict['a1o1'].split()[3])-5, par_dict['scan_startpos3']-float(raw_dict['a1o1'].split()[3])+5, 2, par_dict, plot_dict) #= y-wire
            '''
            wired11 = self.get_wire_sum_center_rms(raw_dict['a1p'].split(), raw_dict['a1s1'].split(), par_dict['scan_stoppos1']-float(raw_dict['a1o1'].split()[3]), par_dict['scan_startpos1']-float(raw_dict['a1o1'].split()[3]), 0, par_dict, plot_dict) #= u-wire #= date, signal sum, center, rms
            wired21 = self.get_wire_sum_center_rms(raw_dict['a1p'].split(), raw_dict['a1s2'].split(), par_dict['scan_stoppos2']-float(raw_dict['a1o1'].split()[3]), par_dict['scan_startpos2']-float(raw_dict['a1o1'].split()[3]), 1, par_dict, plot_dict) #= x-wire
            wired22 = self.get_wire_sum_center_rms(raw_dict['a1p'].split(), raw_dict['a1s3'].split(), par_dict['scan_stoppos3']-float(raw_dict['a1o1'].split()[3]), par_dict['scan_startpos3']-float(raw_dict['a1o1'].split()[3]), 2, par_dict, plot_dict) #= y-wire

            xc,yc,uc,vc = self.get_xyuv_center(wired11['c'],wired21['c'],wired22['c'],par_dict)
            try:
                xr,yr,ur,vr = self.get_xyuv_size(wired11['r'],wired21['r'],wired22['r'],par_dict)
            except:
                xr,yr,ur,vr = -999.,-999.,-999.,-999.
            try:
                x90p,y90p,u90p,v90p = self.get_xyuv_size(wired11['r90p'],wired21['r90p'],wired22['r90p'],par_dict)
            except:
                x90p,y90p,u90p,v90p = -999.,-999.,-999.,-999.
            try:
                x95p,y95p,u95p,v95p = self.get_xyuv_size(wired11['r95p'],wired21['r95p'],wired22['r95p'],par_dict)
            except:
                x95p,y95p,u95p,v95p = -999.,-999.,-999.,-999.
            try:
                x99p,y99p,u99p,v99p = self.get_xyuv_size(wired11['r99p'],wired21['r99p'],wired22['r99p'],par_dict)
            except:
                x99p,y99p,u99p,v99p = -999.,-999.,-999.,-999.
            try:
                cxy = self.get_cxy(xr,yr,ur,vr)
            except:
                cxy = -999.
            try:
                cxy90p = self.get_cxy(x90p,y90p,u90p,v90p)
            except:
                cxy90p = -999.
            try:
                cxy95p = self.get_cxy(x95p,y95p,u95p,v95p)
            except:
                cxy95p = -999.
            try:
                cxy99p = self.get_cxy(x99p,y99p,u99p,v99p)
            except:
                cxy99p = -999.
        elif par_dict['coord'] == "Sbxy":
            wired11 = {'d': 0.0, 's': 0.0, 'c': 0.0, 'r': 0.0, 'c0': 0.0, 'c1': 0.0, 'r90p': 0.0, 'r95p': 0.0, 'r99p' : 0.0}
            wired21 = self.get_wire_sum_center_rms(raw_dict['a1p'].split(),raw_dict['a1s2'].split(),par_dict['scan_stoppos2']-float(raw_dict['a1o1'].split()[3]),par_dict['scan_startpos2']-float(raw_dict['a1o1'].split()[3]),1,par_dict,plot_dict) #= x-wire
            wired22 = self.get_wire_sum_center_rms(raw_dict['a1p'].split(),raw_dict['a1s3'].split(),par_dict['scan_stoppos3']-float(raw_dict['a1o1'].split()[3]),par_dict['scan_startpos3']-float(raw_dict['a1o1'].split()[3]),2,par_dict,plot_dict) #= y-wire
            xc,yc,uc,vc = self.get_xyuv_center(wired11['c'],wired21['c'],wired22['c'],par_dict)
            xr,yr,ur,vr = self.get_xyuv_size(wired11['r'],wired21['r'],wired22['r'],par_dict)
            x90p,y90p,u90p,v90p = self.get_xyuv_size(wired11['r90p'],wired21['r90p'],wired22['r90p'],par_dict)
            x95p,y95p,u95p,v95p = self.get_xyuv_size(wired11['r95p'],wired21['r95p'],wired22['r95p'],par_dict)
            x99p,y99p,u99p,v99p = self.get_xyuv_size(wired11['r99p'],wired21['r99p'],wired22['r99p'],par_dict)
            cxy = -999.
            cxy90p = -999.
            cxy95p = -999.
            cxy99p = -999.

        '''
        plot_dict['r11'],plot_dict['r21'],plot_dict['r22'] = round(r11, 3),round(r21, 3),round(r22, 3)
        plot_dict['c11'],plot_dict['c21'],plot_dict['c22'] = round(c11, 3),round(c21, 3),round(c22, 3)
        plot_dict['s11'],plot_dict['s21'],plot_dict['s22'] = round(s11*uA, 3),round(s21*uA, 3),round(s22*uA, 3)
        plot_dict['d11'],plot_dict['d21'] = d11,d21
        '''
        plot_dict['r11'],plot_dict['r21'],plot_dict['r22'] = round(wired11['r'], 3),round(wired21['r'], 3),round(wired22['r'], 3)
        plot_dict['r90p11'],plot_dict['r90p21'],plot_dict['r90p22'] = round(wired11['r90p'], 3),round(wired21['r90p'], 3),round(wired22['r90p'], 3)
        plot_dict['r95p11'],plot_dict['r95p21'],plot_dict['r95p22'] = round(wired11['r95p'], 3),round(wired21['r95p'], 3),round(wired22['r95p'], 3)
        plot_dict['r99p11'],plot_dict['r99p21'],plot_dict['r99p22'] = round(wired11['r99p'], 3),round(wired21['r99p'], 3),round(wired22['r99p'], 3)
        plot_dict['c011'],plot_dict['c021'],plot_dict['c022'] = round(wired11['c0'], 3),round(wired21['c0'], 3),round(wired22['c0'], 3)
        plot_dict['c111'],plot_dict['c121'],plot_dict['c122'] = round(wired11['c1'], 3),round(wired21['c1'], 3),round(wired22['c1'], 3)
        plot_dict['c11'],plot_dict['c21'],plot_dict['c22'] = round(wired11['c'], 3),round(wired21['c'], 3),round(wired22['c'], 3)
        plot_dict['s11'],plot_dict['s21'],plot_dict['s22'] = round(wired11['s']*uA, 3),round(wired21['s']*uA, 3),round(wired22['s']*uA, 3)
        plot_dict['d11'],plot_dict['d21'] = wired11['d'],wired21['d']
        plot_dict['xc'],plot_dict['yc'],plot_dict['uc'],plot_dict['vc'] = round(xc, 3),round(yc, 3),round(uc, 3),round(vc, 3)
        plot_dict['xr'],plot_dict['yr'],plot_dict['ur'],plot_dict['vr'] = round(xr, 3),round(yr, 3),round(ur, 3),round(vr, 3)
        plot_dict['x90p'],plot_dict['y90p'],plot_dict['u90p'],plot_dict['v90p'] = round(x90p, 3),round(y90p, 3),round(u90p, 3),round(v90p, 3)
        plot_dict['x95p'],plot_dict['y95p'],plot_dict['u95p'],plot_dict['v95p'] = round(x95p, 3),round(y95p, 3),round(u95p, 3),round(v95p, 3)
        plot_dict['x99p'],plot_dict['y99p'],plot_dict['u99p'],plot_dict['v99p'] = round(x99p, 3),round(y99p, 3),round(u99p, 3),round(v99p, 3)
        plot_dict['cxy'] = round(cxy, 3)
        plot_dict['cxy90p'] = round(cxy90p, 3)
        plot_dict['cxy95p'] = round(cxy95p, 3)
        plot_dict['cxy99p'] = round(cxy99p, 3)


        #= Data output
        l = ''
        if form == 1:
            for c in par_dict['coord'][1:]:
                l += "\t"+c+"-wire\t"
            print (l)
            print ("\t"+str(wired11['d'])+"\t\t"+str(wired21['d']))
            print ("sum:\t"+str(round(wired11['s']*uA,7))+"\t"+str(round(wired21['s']*uA,7))+"\t"+str(round(wired22['s']*uA,7)))
            print ("center:\t"+str(round(wired11['c'],7))+"\t"+str(round(wired21['c'],7))+"\t"+str(round(wired22['c'],7)))
            print ("rms:\t"+str(round(wired11['r'],7))+"\t"+str(round(wired21['r'],7))+"\t"+str(round(wired22['r'],7)))
            print ("r90p:\t"+str(round(wired11['r90p'],7))+"\t"+str(round(wired21['r90p'],7))+"\t"+str(round(wired22['r90p'],7)))
            print ("r95p:\t"+str(round(wired11['r95p'],7))+"\t"+str(round(wired21['r95p'],7))+"\t"+str(round(wired22['r95p'],7)))
            print ("r99p:\t"+str(round(wired11['r99p'],7))+"\t"+str(round(wired21['r99p'],7))+"\t"+str(round(wired22['r99p'],7)))

            print ("\tx\t\ty\t\tu\t\tv")
            print ("pos:\t"+str(round(xc,7))+"\t"+str(round(yc,7))+"\t"+str(round(uc,7))+"\t"+str(round(vc,7)))
            print ("rms:\t"+str(round(xr,7))+"\t"+str(round(yr,7))+"\t"+str(round(ur,7))+"\t"+str(round(vr,7)))
            print ("r90p:\t"+str(round(x90p,7))+"\t"+str(round(y90p,7))+"\t"+str(round(u90p,7))+"\t"+str(round(v90p,7)))
            print ("r95p:\t"+str(round(x95p,7))+"\t"+str(round(y95p,7))+"\t"+str(round(u95p,7))+"\t"+str(round(v95p,7)))
            print ("r99p:\t"+str(round(x99p,7))+"\t"+str(round(y99p,7))+"\t"+str(round(u99p,7))+"\t"+str(round(v99p,7)))
            print ("cxy:\t"+str(round(cxy,7)))
            print ("cxy90p:\t"+str(round(cxy90p,7)))
            print ("cxy95p:\t"+str(round(cxy95p,7)))
            print ("cxy99p:\t"+str(round(cxy99p,7)))
            #print "file: "+par_dict["file_name"]
            print ("###########")
        elif form == 2:
            l = par_dict["file_name"]+"\t"
            l+= str(round(wired11['s']*uA,3))+"\t"+str(round(wired21['s']*uA,3))+"\t"+str(round(wired22['s']*uA,3))+"\t"
            l+= str(round(wired11['c'],3))+"\t"+str(round(wired21['c'],3))+"\t"+str(round(wired22['c'],3))+"\t"
            l+= str(round(wired11['r'],3))+"\t"+str(round(wired21['r'],3))+"\t"+str(round(wired22['r'],3))+"\t"
            l+= str(round(xc,3))+"\t"+str(round(yc,3))+"\t"+str(round(uc,3))+"\t"+str(round(vc,3))+"\t"
            l+= str(round(xr,3))+"\t"+str(round(yr,3))+"\t"+str(round(ur,3))+"\t"+str(round(vr,3))+"\t"	
            l+= str(round(cxy,3))+"\t"
            l+= str(round(x90p,3))+"\t"+str(round(y90p,3))+"\t"+str(round(u90p,3))+"\t"+str(round(v90p,3))+"\t"
            l+= str(round(cxy90p,3))+"\t"
            l+= str(round(x95p,3))+"\t"+str(round(y95p,3))+"\t"+str(round(u95p,3))+"\t"+str(round(v95p,3))+"\t"
            l+= str(round(cxy95p,3))+"\t"
            l+= str(round(x99p,3))+"\t"+str(round(y99p,3))+"\t"+str(round(u99p,3))+"\t"+str(round(v99p,3))+"\t"
            l+= str(round(cxy99p,3))
			
            print (l)

            f = open("anaProfileMonitor.out","a")
            f.write(l+"\n")
            f.close()

        #print plot_dict

        '''
        if par_dict['live']:
            if ('LS1' not in par_dict['name']) and ('D0998' not in par_dict['name']):
                name = par_dict['name']
                caput(name+':XCEN_CSET',xc)
                caput(name+':YCEN_CSET',yc)
                caput(name+':XRMS_CSET',xr)
                caput(name+':YRMS_CSET',yr)
                caput(name+':CXY_CSET',cxy)
        '''

        
        if 'sigl' not in wired22.keys():
            r = (('name',par_dict['name']),('file',par_dict["file_name"]),('coord',par_dict['coord']),
                 ('sum1',wired11['s']),('sum2',wired21['s']),
                 ('cen1',wired11['c']),('cen2',wired21['c']),
                 ('cen01',wired11['c0']),('cen02',wired21['c0']),
                 ('rms1',wired11['r']),('rms2',wired21['r']),
                 ('r90p1',wired11['r90p']),('r90p2',wired21['r90p']),
                 ('r95p1',wired11['r95p']),('r95p2',wired21['r95p']),
                 ('xrms',xr),('yrms',yr),('urms',ur),('vrms',vr),
                 ('x90p',x90p),('y90p',y90p),('u90p',u90p),('v90p',v90p),
                 ('x95p',x95p),('y95p',y95p),('u95p',u95p),('v95p',v95p),
                 ('x99p',x99p),('y99p',y99p),('u99p',u99p),('v99p',v99p),
                 ('xcen',xc),('ycen',yc),('ucen',uc),('vcen',vc),
                 ('cxy',cxy), ('cxy90p',cxy90p), ('cxy95p',cxy95p), ('cxy99p',cxy99p),
                 ('sig1',wired11['sigl']),('pos1',wired11['posl']),
                 ('sig2',wired21['sigl']),('pos2',wired21['posl'])
            )
        elif 'sigl' not in wired11.keys():
            r = (('name',par_dict['name']),('file',par_dict["file_name"]),('coord',par_dict['coord']),
                 ('sum1',wired21['s']),('sum2',wired22['s']),
                 ('cen1',wired21['c']),('cen2',wired22['c']),
                 ('cen01',wired21['c0']),('cen02',wired22['c0']),
                 ('rms1',wired21['r']),('rms2',wired22['r']),
                 ('r90p1',wired21['r90p']),('r90p2',wired22['r90p']),
                 ('r95p1',wired21['r95p']),('r95p2',wired22['r95p']),
                 ('xrms',xr),('yrms',yr),('urms',ur),('vrms',vr),
                 ('x90p',x90p),('y90p',y90p),('u90p',u90p),('v90p',v90p),
                 ('x95p',x95p),('y95p',y95p),('u95p',u95p),('v95p',v95p),
                 ('x99p',x99p),('y99p',y99p),('u99p',u99p),('v99p',v99p),
                 ('xcen',xc),('ycen',yc),('ucen',uc),('vcen',vc),
                 ('cxy',cxy), ('cxy90p',cxy90p), ('cxy95p',cxy95p), ('cxy99p',cxy99p),
                 ('sig1',wired21['sigl']),('pos1',wired21['posl']),
                 ('sig2',wired22['sigl']),('pos2',wired22['posl'])
            )
        else:
            r = (('name',par_dict['name']),('file',par_dict["file_name"]),('coord',par_dict['coord']),
                 ('sum1',wired11['s']),('sum2',wired21['s']),('sum3',wired22['s']),
                 ('cen1',wired11['c']),('cen2',wired21['c']),('cen3',wired22['c']),
                 ('cen01',wired11['c0']),('cen02',wired21['c0']),('cen03',wired22['c0']),
                 ('rms1',wired11['r']),('rms2',wired21['r']),('rms3',wired22['r']),
                 ('r90p1',wired11['r90p']),('r90p2',wired21['r90p']),('r90p3',wired22['r90p']),
                 ('r95p1',wired11['r95p']),('r95p2',wired21['r95p']),('r95p3',wired22['r95p']),
                 ('xrms',xr),('yrms',yr),('urms',ur),('vrms',vr),
                 ('x90p',x90p),('y90p',y90p),('u90p',u90p),('v90p',v90p),
                 ('x95p',x95p),('y95p',y95p),('u95p',u95p),('v95p',v95p),
                 ('x99p',x99p),('y99p',y99p),('u99p',u99p),('v99p',v99p),
                 ('xcen',xc),('ycen',yc),('ucen',uc),('vcen',vc),
                 ('cxy',cxy), ('cxy90p',cxy90p), ('cxy95p',cxy95p), ('cxy99p',cxy99p),
                 ('sig1',wired11['sigl']),('pos1',wired11['posl']),
                 ('sig2',wired21['sigl']),('pos2',wired21['posl']),
                 ('sig3',wired22['sigl']),('pos3',wired22['posl'])
            )

             
        res = OrderedDict(r)
        #with open(par_dict['file_name'].replace('.dat','.json'), 'w') as outfile:
        #    json.dump(r, outfile, indent = 4)
        fn = par_dict['file_name'].split('/')[-1].replace('.dat','.json')
        #print(par_dict['file_name'])
        #print(par_dict['file_dict'])
        #print(fn)
        f = open(par_dict['file_dict']+'/'+fn, 'w')
        json.dump(res, f, indent=4)
            
        self.plot(par_dict, plot_dict)

        return res

    #===
    def execute_all(self):
        argc = len(sys.argv)
        if argc == 1:
            fnout = "anaProfileMonitor.out"
            f = open(fnout,"w")
            f.write("file\tsum1\tsum2\tsum3\tcenter1\tcenter2\tcenter3\trms1\trms2\trms3\tpos-x\tpos-y\tpos-u\tpos-v\trms-x\trms-y\trms-u\trms-v\tcxy\tx90p\ty90p\tu90p\tv90p\tcxy90p\tx99p\ty99p\tu99p\tv99p\tcxy99p\n")
            f.close()

        
            dir = os.path.abspath(".")

            if self.debug:
                f1 = open("debugout.out","w")
                f1.write(dir+"\n")
                f1.close()

            inpl = os.listdir(dir)
            for inp in inpl:
                if ".dat" in inp:
                    print (inp)
                    if self.debug:
                        f1 = open("debugout.out","a")
                        f1.write(inp+"\t")
                        f1.close()
                    try:
                        self.execute(inp,2)
                    except:
                        f = open("anaProfileMonitor.out","a")
                        f.write(inp+"\terror\n")

            print (dir)

            
if __name__ == '__main__':

    #pmp = profilemonitor_plot(debug = False, bgfunc = 'pol1')
    #pmp = profilemonitor_plot(debug = False, bgfunc = 'const')
    pmp = profilemonitor_plot(debug = True, batch = True, bgfunc = 'const')
    argc = len(sys.argv)
    #print sys.argv
    if argc == 1:
        pmp.execute_all()
    elif argc == 2:
        inp = sys.argv[1]
        pmp.execute(inp,1)

