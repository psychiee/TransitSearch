# -*- coding: utf-8 -*-
"""
Created on Tue May 25 12:46:09 2021

@author: USER
"""
import numpy as np
import xml.etree.ElementTree as ET, gzip, io
from urllib.request import urlopen, urlretrieve

url = "https://github.com/OpenExoplanetCatalogue/oec_gzip/raw/master/systems.xml.gz"
oec = ET.parse(gzip.GzipFile(fileobj=io.BytesIO(urlopen(url).read())))

# Output mass and radius of all planets 
f = open('oeclist-210525.txt', 'w')
for system in oec.findall(".//system"):
    system_name = system.findtext("name")
    print(system_name)
    try:
        RA = system.findtext("rightascension")
        r1, r2, r3 = RA.split()
        RA_DEG = 15*(float(r1) + float(r2)/60 + float(r3)/3600)
        DEC = system.findtext("declination")
        d1, d2, d3 = DEC.split()
        if d1.startswith("-"):
            DEC_DEG = float(d1) - float(d2)/60 - float(d3)/3600
        else:
            DEC_DEG = float(d1) + float(d2)/60 + float(d3)/3600
    except:
        print(system.findtext("rightascension"), system.findtext("declination"))
    for star in system.findall(".//star"):
        print(star.findtext("name"))
        try:
            RSTAR = float(star.findtext("radius"))
            VMAG = float(star.findtext("magV"))
        except:
            print(star.findtext("radius"))
            continue
        for planet in star.findall("planet"):
            try: 
                RPLANET = float(planet.findtext("radius"))
                SMA = float(planet.findtext("semimajoraxis"))
                PERIOD = float(planet.findtext("period"))
                TT = float(planet.findtext("transittime"))
            except:
                print(planet.findtext("radius"), 
                      planet.findtext("semimajoraxis"),
                      planet.findtext("period"),
                      planet.findtext("transittime"))
                continue
            NAME = planet.findtext("name")            
            for nn in planet.findall("name"):
                if nn.text.startswith(system_name):
                    NAME = nn.text
            #print(NAME)
            outs = "%20s" % NAME.replace(' ','_')
            outs += " %12.8f" % PERIOD
            outs += " %20.8f" % TT 
            T14 = 60*24*PERIOD * (2*RPLANET*0.10045 + 2*RSTAR) / (2*np.pi*SMA*214.9395)
            outs += " %10.3f" % T14
            DEPTH = (RPLANET/RSTAR*0.10045)**2
            outs += " %10.5f" % DEPTH
            outs += " %8.3f" % VMAG
            outs += " %10.6f" % RA_DEG
            outs += " %10.6f" % DEC_DEG 
            outs += " %10.5f" % RSTAR
            outs += " %10.5f" % SMA
            outs += " %10.5f" % RPLANET 
            outs += " %15s" % RA
            outs += " %15s" % DEC
            f.write(outs+'\n')
        
f.close()

    