# -*- coding: utf-8 -*-
'''
;Draw the airmass profile of objects
;1st coded 2004 by K.H.Lee
;2nd coded 2005/Jan/10 by H.S.Hwang
;3rd coded 2005/Jan/26 Sun rise & set time, object altitude OK
;4th coded 2005/Jan/27 Input file & Plot OK and Altitude instead of airmass OK
;5th coded 2005/Jan/27 Sun RA & DEC calculation was modified
;6th coded 2005/Feb/15 Plot of Data start time was modified
;7th coded 2005/May/03 Airmass axis was inserted & Color order was shuffled
;8th updated 2006/June/06 Debug the case of 1 object & epoch keyword inserted
;Last updated 2008/July/15 (JD + 1) and multi_targetlist file system
; ### transit altitude program ###
; 2008/07/15 by wskang
; 2017/03/03 use the data from file directly
;   
; read date
; read transit time table
; plot transit target altitude
'''
import time
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.backends.backend_tkagg as tkagg
from tkinter import *
from tslib import *

dtor = np.pi/180.0000

data = 1000
path ='./plots/'

class TransitSearch(Frame):
    
    def __init__(self, master):

       # Database of Observatory Location
        odat = np.genfromtxt('obsdb.dat', dtype=None, \
               names=['obsdb','chidb','lambdadb','gmtdb'], delimiter=',') 
        self.obsdb = np.array(odat['obsdb'], dtype='unicode')
        self.chidb = odat['chidb']
        self.lambdadb = odat['lambdadb']
        self.gmtdb = odat['gmtdb'] 
        
        self.monthdb=['Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.']
        
        # define colors =====================================================================
        self.colors = \
        [u'indigo', u'firebrick', u'indianred', u'darkolivegreen', u'olive', u'tomato',
         u'orangered', u'darkslategrey', u'dimgray', u'black', u'orange', u'darkslategray', 
         u'brown', u'dodgerblue', u'chocolate', u'crimson', u'forestgreen', u'gray', 
         u'darkturquoise', u'goldenrod', u'darkgreen', u'darkviolet',  u'saddlebrown', 
         u'grey', u'darkslateblue', u'mediumvioletred', u'red', u'deeppink', u'limegreen', 
         u'darkmagenta', u'darkgoldenrod', u'maroon', u'yellowgreen', u'navy', u'olivedrab', 
         u'blue', u'slateblue', u'darkblue', u'seagreen', u'sienna', u'mediumblue', 
         u'royalblue', u'green', u'midnightblue', u'darkcyan', u'teal', u'darkorchid', 
         u'deepskyblue', u'rebeccapurple', u'darkred', u'steelblue',  u'mediumseagreen', 
         u'cadetblue', u'purple', u'darkorange', u'blueviolet']
        
        # read transit database 
        self.params = readData2()
        
        #GUI design
        Frame.__init__(self, master)
        master.title("Transit Search 0.9")
        
        f0 = LabelFrame(master, text="Location", padx=5, pady=5)
        # option menu - observatory list 
        self.obslist = ['Location'] + list(self.obsdb)
        self.OptObs = StringVar()

        w = OptionMenu(f0, self.OptObs, *self.obslist)
        w.grid(row=1,column=1,sticky=W,padx=5,pady=5)
        self.OptObs.set(self.obslist[1])
        
        self.lon, self.lat = StringVar(), StringVar()
        Label(f0,text='Longitude').grid(row=1, column=2)
        Entry(f0,width=15, textvariable=self.lon).grid(row=1, column=3)
        Label(f0,text='Latitude').grid(row=1, column=4)
        Entry(f0,width=15, textvariable=self.lat).grid(row=1, column=5)

        ptime = time.localtime()
        self.OptDate = StringVar()
        self.OptDate.set('%4d/%02d/%02d' % (ptime[0],ptime[1],ptime[2]))
        
        Entry(f0,width=10,textvariable=self.OptDate).grid(row=1,column=6,sticky=W,padx=5,pady=5)
        Frame(f0).grid(row=1,column=7)                                    
        Button(f0, text="Draw", width=15, command=self.Draw).grid(row=1,column=8,sticky=E,padx=5,pady=5)
        
        # embed the plot figure
        f2 = LabelFrame(master, text="Plot", padx=5, pady=5) 
        self.fig = plt.figure(1,figsize=(10,7))
        self.fig.clf()
        self.canvas = tkagg.FigureCanvasTkAgg(self.fig, master=f2)
        self.canvas.get_tk_widget().pack(fill=BOTH)
        
        #implement GUI controls 
        f0.pack(fill=BOTH)
        #f1.pack(fill=BOTH)
        f2.pack(fill=BOTH)

    def Draw(self):

        # set location 
        obsid = self.obslist.index(self.OptObs.get())
        print (obsid, self.OptObs.get())
        if obsid == 0:
            self.chi= np.float64(self.lat.get()) # Latitude of Observation
            self.lam= np.float64(self.lon.get()) # Longitude of Observation
            #if (lambdad < 0) | (lambdam < 0) | (lambdas < 0): self.lam=-self.lam
            #self.obs='LAT.='+'%3d' % (chid,)+'d '+'%2d' % (chim,)+'m '+ \
            #    '%2d' % (chis,)+'s   LONG.='+'%4d' % (lambdad,)+'d '+ \
            #    '%2d' % (lambdam,)+'m '+'%2d' % (lambdas,)+'s'
            self.obs = 'LAT= %.2f LON=%.2f' % (self.chi, self.lam)
            self.timezone=int(self.lam/15)+1 # Time difference between Local Standard Time and Universal Time
        else:
            self.chi=self.chidb[obsid-1]
            self.lat.set(self.chi)
            self.lam=self.lambdadb[obsid-1]
            self.lon.set(self.lam)
            self.obs=self.obsdb[obsid-1] # Observatory Name
            self.timezone = self.gmtdb[obsid-1]
        if (self.lam < 0): self.lam=360+self.lam     
        try:
            tmp = self.OptDate.get().split('/')  
        except:
            ptime = time.localtime()
            self.OptDate.set('%4d/%02d/%02d' % (ptime[0],ptime[1],ptime[2]))
            tmp = self.OptDate.get().split('/')              
        yy = int(tmp[0])
        mm = int(tmp[1])
        dd = int(tmp[2])
                              
        #targetlist_name = '%04d%02d%02d' % (yy,mm,dd)
        month=self.monthdb[mm-1]  # Month
        
        lday=(np.arange(data, dtype=np.float64)/data)*24.+12. #; Local Standard Time
        
        # calc. Julian Date
        # for yy, mm, dd - local time / convert into UT by (-timezone)
        JD12 =JulDate(yy,mm,dd)# ,12-self.timezone,0,0)
        JD = JD12 - 0.5 # JD at UT=0 in the day 
        # make the array of JD for the day (in Local Standard Time) 
        # LST 0h in the day = (JD-timezone)
        # LST 12h in the day = (JD-timezone) + 0.5 
        JDoneday= (JD - float(self.timezone)/24.0) + 0.5 + np.arange(data, dtype=np.float64)/data 
            
        # calc. Sun RA, Dec     
        rasun, decsun = SunRADec(JD)
        rasun = np.mod(rasun, 360) / 15.0   
        
        # calc. sun set/rise time  
        sunriseUT, sunsetUT = SunRiSetUT(JD,self.lam,self.chi)
        sunriseLST = (sunriseUT + self.timezone) % 24.0
        sunsetLST = (sunsetUT + self.timezone) % 24.0
        nosun, = np.where((lday >= sunsetLST) & (lday <= sunriseLST+24.0))
        sunsetjd = min(JDoneday[nosun])
        sunrisejd = max(JDoneday[nosun])
        print ("[ %04d %02d %02d ]" % (yy, mm, dd))
        print ("   JD = ", JD, " at UT=0")
        print ("   SUNRA=", rasun, " hr, SUNDEC=", decsun, " deg")
        print ("   SUNRise=", sunriseLST, " SUNSet=", sunsetLST)
        
        # REFORM TRANSIT DATA ======================================================
        pname = self.params[0]
        pper = self.params[1]
        pperlower = self.params[2]
        pperupper = self.params[3]
        ptt = self.params[4]
        pt14 = self.params[7]
        pdepth = self.params[10]
        sRAs = self.params[11]
        sDECs = self.params[12] 
        sVs = self.params[13]
        
        # total number of planets 
        pnum = len(pname)
        
        # CHECK AVAILABLE TRANSITS ======================================================
        dJD = sunrisejd - sunsetjd 
        # check whether the transit-mid time is in night-time 
        N_revol = np.array((sunrisejd - ptt) / pper, np.int)
        JD_recent = N_revol * pper + ptt
        JD_recent1 = JD_recent - pt14/2.0
        JD_recent2 = JD_recent + pt14/2.0 
        JD_err1 = - N_revol * pperlower # err2 is negative 
        JD_err2 = N_revol * pperupper # err1 is positive
        
        oklist = np.where((JD_recent2 < sunrisejd) & (JD_recent1 > sunsetjd) & \
                          (pdepth > 0.005) & (sVs < 20))[0]

        # draw two figure frames 
        f2 = self.fig 
        ax2 = f2.add_axes([0.17,0.1,0.75,0.86])
        ax2.cla()

        # set of vertical plot position of transit timing 
        y2 = 5
        lyy = np.zeros([data])
        ytickv = []
        for pidx in oklist: # = 0, pnum-1 DO BEGIN
            # read star RA, Dec
            starRA, starDec = sRAs[pidx], sDECs[pidx]
        
            # check altitude of the star
            starAlt, starAzi = \
                StarAltAzi(JDoneday,self.lam,self.chi,starRA,starDec)
            tall = np.where((JDoneday >= JD_recent1[pidx]+JD_err1[pidx]) & \
                       (JDoneday <= JD_recent2[pidx]+JD_err2[pidx]))[0]
            tdur = np.where((JDoneday >= JD_recent1[pidx]) & \
                       (JDoneday <= JD_recent2[pidx]))[0]
            # check the timing of altitude > 30 deg 
            if np.min(starAlt[tdur]) < 25: continue
            
            # mid-transit time 
            tcen = int(tdur.mean())
            # define unique color of each planet
            pcolor = self.colors[pidx % len(self.colors)]
            # draw line and transit time and labels 
            ax2.plot(lday,lyy+y2, '-', color=pcolor, lw=2, alpha=0.5) 
            ax2.text(sunsetLST-(sunriseLST-sunsetLST+26)*0.25,y2, \
                     pname[pidx],color=pcolor, fontsize=12, )
            ax2.text(sunsetLST-(sunriseLST-sunsetLST+26)*0.27,y2+4, \
                     '%5.2f' % (sVs[pidx],),color=pcolor, fontsize=10)
            ax2.text(sunriseLST+24+(sunriseLST-sunsetLST+26)*0.08,y2, \
                     '%5.2f%%' % (pdepth[pidx]*100,), color=pcolor, fontsize=10)
            # add ttick position 
            ytickv.append(y2+3)
            # draw the transit timing including the errors of pericenter, period
            ax2.plot(lday[tall],lyy[tall]+y2, color=pcolor, lw=7, alpha=0.5)
            
            # draw the duration-based transit timing 
            oalt = np.where(starAlt > 0)[0]
            pys = np.r_[lyy[oalt]+y2+starAlt[oalt]/9,\
                        lyy[oalt[::-1]]+y2]
            pxs = np.r_[lday[oalt],lday[oalt[::-1]]]
            xypair = list(zip(pxs,pys))
            p = Polygon(xypair,closed=True,fill=True,color=pcolor,alpha=0.2)
            ax2.add_patch(p)
            # draw the mid-transit 
            ax2.plot(lday[tcen],lyy[tcen]+y2, 'o', color=pcolor, \
                     ms=12, alpha=0.5)
    
            print (pname[pidx], int(starAlt[tcen]), sVs[pidx])
            # to the next row 
            y2 = y2 + 10
            
        # check the available stars 
        if y2 == 5:
            #ax1.cla()
            ax2.text(0.5, 0.3, 'No Available Targets', \
                     fontsize=25, alpha=0.6, transform=ax2.transAxes)
            f2.canvas.draw() 
            return
          
        # draw the sunset/rise line in plot 
        ax2.plot([sunsetLST,sunsetLST],[0,y2+10], 'k-', lw=1)
        ax2.text(sunsetLST-0.27,(y2+10)*0.1,'Sunset', rotation='vertical', fontsize=12)
        ax2.plot([sunsetLST+18/15.,sunsetLST+18/15.],[0,y2+10], 'k-', lw=1) 
        ax2.text(sunsetLST-0.27+18/15.,(y2+10)*0.1, 'Evening Astronomical Twilight',\
                 rotation='vertical', fontsize=12)
        ax2.plot([sunriseLST+24-18/15.,sunriseLST+24-18/15.],[0,y2+10], 'k-', lw=1) 
        ax2.text(sunriseLST+24+0.1-18/15.,(y2+10)*0.1, 'Morning Astronomical Twilight',\
                 rotation='vertical', fontsize=12)
        ax2.plot([sunriseLST+24,sunriseLST+24],[0,y2+10], 'k-', lw=1)
        ax2.text(sunriseLST+24+0.1,(y2+10)*0.1,'Sunrise', rotation='vertical', fontsize=12)
    
        # draw the title     
        #ax2.text(0.04,0.92,'DATE: %s %02d %04d' % (month, dd, yy),fontsize=34,transform=f2.transFigure)
        #ax2.text(0.97,0.94,'SITE: %s' % (self.obs,),fontsize=30,ha='right',transform=f2.transFigure) 
        #ax2.text(0.97,0.88,'LON =%7.2f , LAT =%6.2f' % (self.lam,self.chi),ha='right',fontsize=24,transform=f2.transFigure)
        ax2.set_xlim([sunsetLST-1,sunriseLST+25])
        ax2.set_ylim([0, y2])
        ax2.set_xlabel('Local Standard Time [h]', fontsize=20)
        # read and redraw the x-labels 
        f2.canvas.draw()
        labels = [item.get_text() for item in ax2.get_xticklabels()]
        for i, x in enumerate(labels):
            try: 
                x = int(x)
                if x > 24: labels[i] = '%02d' % (x-24,)
            except: pass          
        ax2.set_xticklabels(labels)
        ax2.set_yticks(ytickv)
        ax2.set_yticklabels([''])
        ax2.grid()
        f2.canvas.draw()        

    
if __name__ == "__main__":
    root = Tk()
    ap = TransitSearch(master=root)
    ap.mainloop()


