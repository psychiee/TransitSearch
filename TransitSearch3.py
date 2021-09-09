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
; 2018/01/02 use the exoplanet.org data 
; read date
; read transit time table
; plot transit target altitude
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import time 
from tslib import *
import os

par = read_params()
PATH = par['PATH']
try:
    os.mkdir(PATH)
    print(f'{PATH} created.')
except:
    print(f'{PATH} already exists.')
DB_NAME = par['DB']
OBSNUM = int(par['OBSNUM'])
lambda_D, lambda_M, lambda_S = np.array(par['LON'].split(','), float)
chi_D, chi_M, chi_S = np.array(par['LAT'].split(','), float)
USER_TZ = int(par['TIMEZONE'])

# 2006/08/28 Kang Wonseok (multi_file mode)
with open('targetdate.txt','r') as f:
    #name,hh,rmin,rsec,deg,dmin
    iyy, imm, idd = [], [], [] 
    for line in f:
        if line.strip() == '*':
            iyy, imm, idd = [],[],[]
            ptime = time.localtime()
            jd_cur = JulDate(ptime[0],ptime[1],ptime[2])#,0,0,0)
            mjd_cur = jd_cur - 2400000.5
            for jj in np.arange(jd_cur,jd_cur+180):
                y,m,d = CalDate(jj)
                iyy.append(y)
                imm.append(m)
                idd.append(d)
        else:    
            tmp = line.split()
            if len(tmp) < 3: continue
            iyy.append(tmp[0].strip())
            imm.append(tmp[1].strip())
            idd.append(tmp[2].strip())

# Database of Observatory Location =============================================
odat = np.genfromtxt('obsdb.dat',dtype='U50, f, f, f', \
       names=['obsdb','chidb','lambdadb','gmtdb'], delimiter=',') 
obsdb = np.array(odat['obsdb'], dtype='unicode')
chidb = odat['chidb']
lambdadb = odat['lambdadb']
gmtdb = odat['gmtdb']

monthdb=['Jan.','Feb.','Mar.','Apr.','May','Jun.','Jul.','Aug.','Sep.','Oct.','Nov.','Dec.']

if OBSNUM == 0:
    chi= chi_D+chi_M/60.+chi_S/3600. # Latitude of Observation
    lam=abs(lambda_D)+abs(lambda_M)/60.+abs(lambda_S)/3600. # Longitude of Observation
    if (lambda_D < 0) | (lambda_M < 0) | (lambda_S < 0): lam = -lam
    obs = 'Other Site'
    timezone=USER_TZ
else:
    chi=chidb[OBSNUM-1]
    lam=lambdadb[OBSNUM-1]
    obs=str.strip(obsdb[OBSNUM-1]) # Observatory Name
    timezone = gmtdb[OBSNUM-1]

if (lam < 0): lam=360+lam

# define colors =====================================================================
colors = \
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

# read transit data ==================================================================
params = readData(DB_NAME)
pname = params[0]
pper = params[1]
pperlower = params[2]
pperupper = params[3]
ptt = params[4]
pt14 = params[7]
pdepth = params[10]
sRAs = params[11]
sDECs = params[12] 
sVs = params[13]

pnum = len(pname)
# date loop ==========================================================================

numdate = len(idd)
    
data=1500 #; Number of Data points
    
for yy, mm, dd in zip(iyy, imm, idd):
    yy, mm, dd = int(yy), int(mm), int(dd) 
    
    targetlist_name = f'{yy:04d}{mm:02d}{dd:02d}'
    month=monthdb[mm-1]  # Month
    
    lday=np.arange(data, dtype=np.float64)/data*24.+12. #; Local Standard Time
    
    # calc. Julian Date
    # for yy, mm, dd - local time / convert into UT by (-timezone)
    JD12 =JulDate(yy,mm,dd)# ,12-self.timezone,0,0)
    JD = JD12 - 0.5 # JD at UT=0 in the day 
    # make the array of JD for the day (in Local Standard Time) 
    # LST 0h in the day = (JD-timezone)
    # LST 12h in the day = (JD-timezone) + 0.5 
    JDoneday= (JD - float(timezone)/24.0) + 0.5 + np.arange(data, dtype=np.float64)/data 
        
    # calc. Sun RA, Dec     
    rasun, decsun = SunRADec(JD)
    rasun = np.mod(rasun, 360) / 15.0   
    
    # calc. sun set/rise time  
    sunriseUT, sunsetUT = SunRiSetUT(JD,lam,chi)
    sunriseLST = (sunriseUT + timezone) % 24.0
    sunsetLST = (sunsetUT + timezone) % 24.0
    nosun = np.where((lday >= sunsetLST) & (lday <= sunriseLST+24.0))[0]
    sunsetjd = min(JDoneday[nosun])
    sunrisejd = max(JDoneday[nosun])    
    
    print(f"[ {yy:04d} {mm:02d} {dd:02d} ]")
    print(f"   JD = {JD}, TIMEZONE = {timezone}")
    print(f"   SUNRA = {rasun}, SUNDEC = {decsun}")
    print(f"   SUNSet = {sunsetLST}, SUNRise = {sunriseLST}")

    # CHECK AVAILABLE TRANSITS ======================================================
    dJD = sunrisejd - sunsetjd 
    # check whether the transit-mid time is in night-time 
    N_revol = np.array((sunrisejd - ptt) / pper, np.int)
    
    JD_recent = N_revol * pper + ptt
    JD_recent1 = JD_recent - pt14/2.0
    JD_recent2 = JD_recent + pt14/2.0 
    JD_err1 = N_revol * pperlower
    JD_err2 = N_revol * pperupper
    
    oklist = np.where((JD_recent1 > sunsetjd) & (JD_recent2 < sunrisejd) &
                      (pdepth > 0.004) & (sVs < 20))[0]

    # draw two figure frames 
    f2 = plt.figure(2,figsize=(18,12))
    ax2 = f2.add_axes([0.2,0.1,0.70,0.75])
    # set of vertical plot position of transit timing 
    y2 = 5
    lyy = np.zeros([data])
    ytickv = []
    for pidx in oklist: # = 0, pnum-1 DO BEGIN
        # read star RA, Dec
        starRA, starDec = sRAs[pidx], sDECs[pidx]
        
        # check altitude of the star
        starAlt, starAzi = \
            StarAltAzi(JDoneday,lam,chi,starRA,starDec)
        tall = np.where((JDoneday >= JD_recent1[pidx]-JD_err1[pidx]) & \
                      (JDoneday <= JD_recent2[pidx]+JD_err2[pidx]))[0]
        tdur = np.where((JDoneday >= JD_recent1[pidx]) & \
                        (JDoneday <= JD_recent2[pidx]))[0]
        # check the timing of altitude > 30 deg
        if len(tdur) == 0: continue
        if np.min(starAlt[tdur]) < 25: continue

        # mid-transit time 
        tcen= int(np.median(tdur))
        # define unique color of each planet
        pcolor = colors[pidx % len(colors)]
        # draw line and transit time and labels 
        ax2.plot(lday,lyy+y2, '-', color=pcolor, lw=2, alpha=0.5) 
        ax2.text(sunsetLST-(sunriseLST-sunsetLST+26)*0.3,y2, \
                 pname[pidx],color=pcolor, fontsize=16, )
        ax2.text(sunsetLST-(sunriseLST-sunsetLST+26)*0.34,y2+4, \
                 f'{sVs[pidx]:5.2f}', color=pcolor, fontsize=12)
        ax2.text(sunriseLST+24+(sunriseLST-sunsetLST+26)*0.1,y2, \
                 f'{pdepth[pidx]*100:5.1f}%', color=pcolor, fontsize=14)
        # add ttick position 
        ytickv.append(y2+3)
        # draw the transit timing including the errors of pericenter, period
        ax2.plot(lday[tall],lyy[tall]+y2, color=pcolor, lw=5, alpha=0.5)
        
        ax2.plot(lday[tdur],lyy[tdur]+y2, color=pcolor, lw=12, alpha=0.5)
        # mark JD start and end time 
        ax2.text(lday[tdur[0]]-0.5,lyy[tdur[0]]+y2, \
                 f'{JD_recent1[pidx]-JD12:.2f}', \
                 rotation='vertical', color=pcolor, fontsize=12 )
        ax2.text(lday[tdur[-1]]+0.25,lyy[tdur[-1]]+y2, \
                 f'{JD_recent2[pidx]-JD12:.2f}', \
                 rotation='vertical', color=pcolor, fontsize=12 )
        
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
        
    if y2 == 5: 
        plt.close('all')
        continue
    # draw the sunset/rise line in plot 
    ax2.plot([sunsetLST,sunsetLST],[0,y2+10], 'k-', lw=1)
    ax2.text(sunsetLST-0.27,(y2+10)*0.1,'Sunset', rotation='vertical', fontsize=15)
    ax2.plot([sunsetLST+18/15.,sunsetLST+18/15.],[0,y2+10], 'k-', lw=1) 
    ax2.text(sunsetLST-0.27+18/15.,(y2+10)*0.1, 'Evening Astronomical Twilight',\
             rotation='vertical', fontsize=15)    
    ax2.plot([sunriseLST+24-18/15.,sunriseLST+24-18/15.],[0,y2+10], 'k-', lw=1) 
    ax2.text(sunriseLST+24+0.1-18/15.,(y2+10)*0.1, 'Morning Astronomical Twilight',\
             rotation='vertical', fontsize=15)    
    ax2.plot([sunriseLST+24,sunriseLST+24],[0,y2+10], 'k-', lw=1)
    ax2.text(sunriseLST+24+0.1,(y2+10)*0.1,'Sunrise', rotation='vertical', fontsize=15)

    # draw the title     
    ax2.text(0.04,0.92, f'DATE: {month} {dd:02d}, {yy:04d}', 
             fontsize=34,transform=f2.transFigure)
    ax2.text(0.97,0.94, f'SITE: {obs}',
             fontsize=30,ha='right',
             transform=f2.transFigure) 
    ax2.text(0.97,0.88, f'LON = {lam:7.2f}, LAT = {chi:+6.2f}',
             fontsize=24,ha='right',
             transform=f2.transFigure)
    ax2.set_xlim([sunsetLST-1,sunriseLST+25])
    ax2.set_ylim([0, y2])
    ax2.set_xlabel('Local Standard Time [h]', fontsize=20)
    
    
    # read and redraw the x-labels 
    f2.canvas.draw()
    labels = [item.get_text() for item in ax2.get_xticklabels()]
    for i, x in enumerate(labels):
        labels[i] = f'{int(x) % 24:02d}'
    ax2.set_xticklabels(labels, fontsize=20) 
    ax2.set_yticks(ytickv)
    ax2.set_yticklabels(['' for x in ytickv])
    ax2.grid(which='both')    
    f2.savefig(PATH+targetlist_name+'_2.png')
    plt.close('all')





