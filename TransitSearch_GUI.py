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
import sys, time
from glob import glob
from tslib import *
import pyqtgraph as pg
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import matplotlib.backends.backend_qt5agg as qtagg
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QLineEdit, \
    QWidget, QMessageBox, QVBoxLayout, QHBoxLayout, QComboBox, QGroupBox, \
    QPushButton
from PyQt5.QtGui import QFont, QColor
from PyQt5.QtCore import Qt

class TransitSearch(QMainWindow):
    
    def __init__(self):
        super().__init__()

        # Check & Read transit database =======================================
        try:
            self.params = readData2()
        except Exception as e:
            emsg = 'Error in reading transit DB file(transit-YYMMDD.dat)\n\n'
            for i, z in enumerate(e.args):
                emsg += '(%i) %s \n' % (i + 1, z)
            QMessageBox(QMessageBox.Critical, 'Error', emsg)
            raise
        # Check & Read Observatory data ========================================
        try:
            self.obsdb = np.genfromtxt('observatory.dat', dtype='U', usecols=(0,), delimiter=',')
            odat = np.genfromtxt('observatory.dat', usecols=(1,2,3), delimiter=',')
            self.chidb, self.lambdadb , self.gmtdb = odat[:,0], odat[:,1], odat[:,2]
        except Exception as e:
            emsg = 'Error in reading observatory DB file(observatory.dat) \n\n '
            for i, z in enumerate(e.args):
                emsg += '(%i) %s \n' % (i + 1, z)
            QMessageBox(QMessageBox.Critical, 'Error', emsg)
            raise

        # GUI design
        qf_small = QFont('Verdana', 10)
        qf = QFont('Verdana', 12)
        qf_bold = QFont('Verdana', 13, QFont.Bold)
        cnames = [z for z in QColor.colorNames()]
        cnames.remove('white')
        self.colors = []
        for z in cnames:
            t = QColor(z).getRgb()
            self.colors.append((t[0], t[1], t[2], 180))

        # MAIN widgets
        self.setWindowTitle('Transit Search 0.90')
        self._centralWidget = QWidget(self)
        self.setCentralWidget(self._centralWidget)
        self.mainLayout = QVBoxLayout()

        # TOP Widgets
        topLayout = QHBoxLayout()
        # - - Combobox of Observatory
        self.obsList = QComboBox(font=qf_small)
        self.obsList.addItem('User Location')
        for z in self.obsdb:
            self.obsList.addItem(z)
        self.obsList.setFixedHeight(35)
        self.obsList.currentIndexChanged.connect(self.ChangeObs)

        # - - LineEdit of LAT, LON

        self.latText = QLineEdit('', font=qf, alignment=Qt.AlignRight)
        self.latText.setFixedSize(80, 30)
        self.latText.setReadOnly(True)
        self.lonText = QLineEdit('', font=qf, alignment=Qt.AlignRight)
        self.lonText.setFixedSize(80, 30)
        self.lonText.setReadOnly(True)
        self.tzText = QLineEdit('', font=qf, alignment=Qt.AlignRight)
        self.tzText.setFixedSize(40, 30)
        self.tzText.setReadOnly(True)

        # - - LineEdit of DATE
        ptime = time.localtime()
        self.dateText = QLineEdit('%4d/%02d/%02d' % (ptime[0], ptime[1], ptime[2]), font=qf_bold)
        self.dateText.setFixedSize(140, 30)
        self.dateText.returnPressed.connect(self.DrawEnter)
        # - - Button of DRAW
        self.drawButton = QPushButton('DRAW')
        self.drawButton.setFixedSize(120, 30)
        self.drawButton.setFont(qf_bold)
        self.drawButton.clicked.connect(self.Draw)

        # TOP Layouts
        topLayout.addWidget(self.obsList)
        topLayout.addWidget(QLabel("Latitude",font=qf_small,width=30))
        topLayout.addWidget(self.latText)
        topLayout.addWidget(QLabel("Longitude", font=qf_small, width=30))
        topLayout.addWidget(self.lonText)
        topLayout.addWidget(QLabel("Timezone", font=qf_small, width=30))
        topLayout.addWidget(self.tzText)
        topLayout.addWidget(QLabel("Date", font=qf_small, width=30))
        topLayout.addWidget(self.dateText)
        topLayout.addWidget(self.drawButton)

        self.mainGraph = pg.PlotWidget()
        self.mainGraph.setBackground('w')
        self.mainLayout.addLayout(topLayout)
        self.mainLayout.addWidget(self.mainGraph)

        self.obsList.setCurrentIndex(1)
        self._centralWidget.setLayout(self.mainLayout)

    def ChangeObs(self, i):
        if i == 0:
            self.latText.setReadOnly(False)
            self.lonText.setReadOnly(False)
        else:
            self.latText.setText('%.3f' % (self.chidb[i-1]))
            self.lonText.setText('%.3f' % (self.lambdadb[i-1]))
            self.tzText.setText('%.1f' % (self.gmtdb[i-1]))
            self.latText.setReadOnly(True)
            self.lonText.setReadOnly(True)

    def DrawEnter(self):
        self.Draw()

    def Draw(self):

        # set location 
        self.chi = np.float64(self.latText.text()) # Latitude of Observation
        self.lam = np.float64(self.lonText.text()) # Longitude of Observation
        self.timezone = np.float64(self.tzText.text())
        self.obs = 'LAT= %.2f LON=%.2f' % (self.chi, self.lam)

        obsid = self.obsList.currentIndex()

        try:
            tmp = self.dateText.text().split('/')
        except:
            ptime = time.localtime()
            self.OptDate.set('%4d/%02d/%02d' % (ptime[0],ptime[1],ptime[2]))
            tmp = self.OptDate.get().split('/')

        yy = int(tmp[0])
        mm = int(tmp[1])
        dd = int(tmp[2])
                              
        #targetlist_name = '%04d%02d%02d' % (yy,mm,dd)
        #month=self.monthdb[mm-1]  # Month
        NPTS = 1000
        lday=(np.arange(NPTS, dtype=np.float64)/NPTS)*24.+12. #; Local Standard Time
        
        # calc. Julian Date
        # for yy, mm, dd - local time / convert into UT by (-timezone)
        JD12 = JulDate(yy,mm,dd)# ,12-self.timezone,0,0)
        JD = JD12 - 0.5 # JD at UT=0 in the day 
        # make the array of JD for the day (in Local Standard Time) 
        # LST 0h in the day = (JD-timezone)
        # LST 12h in the day = (JD-timezone) + 0.5 
        JDoneday= (JD - float(self.timezone)/24.0) + 0.5 + np.arange(NPTS, dtype=np.float64)/NPTS

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
        
        oklist = np.where((JD_recent2 < sunrisejd) & (JD_recent1 > sunsetjd) &
                          (pdepth > 0.005) & (sVs < 20))[0]

        # draw two figure frames 
        mg = self.mainGraph
        mg.clear()

        # set of vertical plot position of transit timing
        y2 = 5
        lyy = np.zeros([NPTS])
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
            ppen1 = pg.mkPen(pcolor, width=2, alpha=0.1)
            ppen3 = pg.mkPen(pcolor, width=15, alpha=0.2)
            pf = QFont('Verdana', 9)
            pf_small = QFont('Verdana', 7)
            # draw line and transit time and labels 
            mg.plot(lday, lyy+y2, pen=ppen1)
            # draw the transit timing including the errors of pericenter, period
            mg.plot(lday[tall], lyy[tall] + y2, pen=ppen3)

            t1 = pg.TextItem(pname[pidx], color=pcolor)
            t1.setFont(pf)
            t1.setPos(sunsetLST-(sunriseLST-sunsetLST+26)*0.27, y2+4)

            mg.addItem(t1)
            '''
            mg.text('%5.2f' % (sVs[pidx],), anchor=(sunsetLST-(sunriseLST-sunsetLST+26)*0.27, y2+4),
                     color=pcolor, font=pf_small)
            mg.text('%5.2f%%' % (pdepth[pidx]*100,), anchor=(sunriseLST+24+(sunriseLST-sunsetLST+26)*0.08, y2),
                     color=pcolor, font=pf_small)
            '''
            # add ttick position 
            ytickv.append(y2+3)
            '''
            # draw the duration-based transit timing 
            oalt = np.where(starAlt > 0)[0]
            pys = np.r_[lyy[oalt]+y2+starAlt[oalt]/9,\
                        lyy[oalt[::-1]]+y2]
            pxs = np.r_[lday[oalt],lday[oalt[::-1]]]
            xypair = list(zip(pxs,pys))
            p = Polygon(xypair,closed=True,fill=True,color=pcolor,alpha=0.2)
            ax2.add_patch(p)
            # draw the mid-transit 
            ax2.plot(lday[tcen],lyy[tcen]+y2, 'o', color=pcolor, ms=12, alpha=0.5)
            '''
            print (pname[pidx], int(starAlt[tcen]), sVs[pidx])
            # to the next row 
            y2 = y2 + 10
            
        # check the available stars 
        if y2 == 5:
            #ax1.cla()
            mg.text('No Available Targets')
            #mg.canvas.draw()
            return
        mg.plot([sunsetLST, sunsetLST], [0, y2 + 10])
        mg.plot([sunsetLST + 18 / 15., sunsetLST + 18 / 15.], [0, y2 + 10])
        mg.plot([sunriseLST+ 24 , sunriseLST+ 24], [0, y2 + 10])
        mg.plot([sunriseLST + 24 + 18 / 15., sunriseLST + 24 + 18 / 15.], [0, y2 + 10])

        '''  
        # draw the sunset/rise line in plot 
        ax2.plot([sunsetLST,sunsetLST],[0,y2+10], 'k-', lw=1)
        ax2.text(sunsetLST-0.27,(y2+10)*0.1,'Sunset', rotation='vertical', fontsize=12)
        ax2.plot([sunsetLST+18/15.,sunsetLST+18/15.],[0,y2+10], 'k-', lw=1) 
        ax2.text(sunsetLST-0.27+18/15.,(y2+10)*0.1, 'Evening Astronomical Twilight',
                 rotation='vertical', fontsize=12)
        ax2.plot([sunriseLST+24-18/15.,sunriseLST+24-18/15.],[0,y2+10], 'k-', lw=1) 
        ax2.text(sunriseLST+24+0.1-18/15.,(y2+10)*0.1, 'Morning Astronomical Twilight',
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
        '''
        xticks = [(k, '%02i' % (k % 24)) for k in range(12,36)]

        mg.getAxis('bottom').setTicks([xticks])
        mg.getPlotItem().showGrid(alpha=0.5)
        mg.setXRange(12, 36)
        mg.setYRange(0, y2)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle('Breeze')
    view = TransitSearch()
    view.show()
    sys.exit(app.exec_())

