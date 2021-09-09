# -*- coding: utf-8 -*-
"""
Created on Tue May 25 12:46:09 2021

@author: USER
"""
import numpy as np
from astroquery.exoplanet_orbit_database import ExoplanetOrbitDatabase

eod_table = ExoplanetOrbitDatabase.get_table()

vlist, = np.where((eod_table["TT"] > 0) & 
                  (eod_table["T14"].value > 0) & 
                  (eod_table["PER"].value > 0) &
                  (eod_table["DEPTH"].value > 0) &
                  (eod_table["V"].value < 13) & 
                  (eod_table["V"].value > 0))

f = open('eodlist-210526.txt', 'w')
for i in vlist:
    NAME = eod_table["NAME"][i]
    PER = eod_table["PER"][i].value
    TT = eod_table["TT"][i]
    T14 = eod_table["T14"][i].value * 24 * 60
    DEPTH = eod_table["DEPTH"][i].value
    #RSTAR = eod_table["RSTAR"][i]
    RA = eod_table["RA"][i].value * 15.0 
    DEC = eod_table["DEC"][i].value
    V = eod_table["V"][i].value
    PERUPPER = eod_table["PERUPPER"][i].value
    PERLOWER = eod_table["PERLOWER"][i].value
    
    outs = "%20s" % NAME.replace(' ','_')
    outs += " %12.8f" % PER
    outs += " %20.8f" % TT 
    #T14 = 60*24*PERIOD * (2*RPLANET*0.10045 + 2*RSTAR) / (2*np.pi*SMA*214.9395)
    outs += " %10.3f" % T14
    #DEPTH = (RPLANET/RSTAR*0.10045)**2
    outs += " %10.5f" % DEPTH
    outs += " %8.3f" % V
    outs += " %10.6f" % RA
    outs += " %10.6f" % DEC
    outs += " %12.8f" % PERUPPER
    outs += " %12.8f" % PERLOWER
    f.write(outs+'\n')
    
f.close()

# # Output mass and radius of all planets 
# f = open('oeclist-210526.txt', 'w')
# for system in oec.findall(".//system"):
#     system_name = system.findtext("name")
#     #print(system_name)
#     try:
#         RA = system.findtext("rightascension")
#         r1, r2, r3 = RA.split()
#         RA_DEG = 15*(float(r1) + float(r2)/60 + float(r3)/3600)
#         DEC = system.findtext("declination")
#         d1, d2, d3 = DEC.split()
#         if d1.startswith("-"):
#             DEC_DEG = float(d1) - float(d2)/60 - float(d3)/3600
#         else:
#             DEC_DEG = float(d1) + float(d2)/60 + float(d3)/3600
#     except:
#         print(system.findtext("rightascension"), system.findtext("declination"))
#         print('RA, Dec are not found')
#     # STAR LOOP ==============================================
#     for star in system.findall(".//star"):
#         star_name = star.findtext("name")
#         print(star_name)
#         # CHECK stellar radius
#         try:
#             RSTAR = float(star.findtext("radius"))
#         except:
#             print('RSTAR is not found:', star.findtext("radius"))
#             continue
#         # PLANET LOOP ==============================================
#         for planet in star.findall("planet"):
#             # CHECK planet parameters 
#             try: 
#                 RPLANET = float(planet.findtext("radius"))
#                 SMA = float(planet.findtext("semimajoraxis"))
#                 PERIOD = float(planet.findtext("period"))
#                 TT = float(planet.findtext("transittime"))
#             except:
#                 print(planet.findtext("radius"), 
#                       planet.findtext("semimajoraxis"),
#                       planet.findtext("period"),
#                       planet.findtext("transittime"))
#                 continue
            
#             # CHECK stellar magnitude 
#             try:
#                 MAG = float(star.findtext("magV"))
#             except:
#                 try:
#                     vot = SB.query_object(star_name)
#                     MAG = vot['FLUX_G'].data[0]
#                 except:
#                     print('MAG_V and MAG_G are not found')
#                     continue
#             # CHECK planet names 
#             NAME = planet.findtext("name")            
#             for nn in planet.findall("name"):
#                 if nn.text.startswith(system_name):
#                     NAME = nn.text
#                     break
#             #print(NAME)
#             outs = "%20s" % NAME.replace(' ','_')
#             outs += " %12.8f" % PERIOD
#             outs += " %20.8f" % TT 
#             T14 = 60*24*PERIOD * (2*RPLANET*0.10045 + 2*RSTAR) / (2*np.pi*SMA*214.9395)
#             outs += " %10.3f" % T14
#             DEPTH = (RPLANET/RSTAR*0.10045)**2
#             outs += " %10.5f" % DEPTH
#             outs += " %8.3f" % MAG
#             outs += " %10.6f" % RA_DEG
#             outs += " %10.6f" % DEC_DEG 
#             outs += " %10.5f" % RSTAR
#             outs += " %10.5f" % SMA
#             outs += " %10.5f" % RPLANET 
#             outs += " %15s" % RA
#             outs += " %15s" % DEC
#             f.write(outs+'\n')
        
# f.close()

    