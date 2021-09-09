# -*- coding: utf-8 -*-
"""
Created on Thu Mar 16 22:34:25 2017

@author: wskang
"""
from lxml import html 
from time import sleep
import requests
from urllib.request import urlopen, urlretrieve
import time 

hroot = 'http://var2.astro.cz/ETD/'

page = requests.get(hroot+'index.php')
webstring = page.content
rtree = html.fromstring(page.content)

tlist = [] 
for x in rtree.xpath('//@href'):
    print(x) 
    if not x.startswith('etd.php?STARNAME'): continue
    tlist.append(x)

tlist = sorted(set(tlist))
NUM = len(tlist)
fdat = open('etdlist-%s.txt' % time.strftime('%y%m%d'), 'w')
for i, x in enumerate(tlist[:]):
    print( '%4i/%4i, %s' % (i+1,NUM, x))
    #if i < 190: continue

    pexo = requests.get(hroot+x)
    ptree = html.fromstring(pexo.content)
    
    name = ptree.xpath('//a[@href="'+x+'"]/text()')[0]
    fname = name.replace(' ','_')
    fname = fname.replace('/','_')
    
    try:
        period = float(ptree.xpath('//input[@name="PER"]/@value')[0])    
        mjd = float(ptree.xpath('//input[@name="M"]/@value')[0])
    except:
        continue
    
    for tbl in ptree.xpath('//table'):
        dat = tbl.xpath('.//tr/td//text()')
        if dat[0] == 'RA':
            ra_str = dat[7].encode('ascii','ignore').decode()
            ra = 15*(float(ra_str[0:2])+float(ra_str[2:4])/60+float(ra_str[4:])/3600)
            dec_str = dat[8].encode('ascii','ignore').decode()
            if (dec_str[0] != '-') & (dec_str[0] != '+'):
                dec_str = '+'+dec_str
            dec = abs(float(dec_str[0:3]))+float(dec_str[3:5])/60+float(dec_str[5:])/3600
            if dec_str[0]== '-': dec = -dec
            mag = float(dat[9])
            dep = float(dat[10])
            dur = float(dat[11])
    fstr = '%20s %12.8f %16.8f %8.3f %8.3f %8.3f %12.6f %+12.6f %s %s' % \
           (fname,period,mjd,dur,dep,mag,ra,dec,ra_str,dec_str) 
    print(fstr)
    fdat.write(fstr+'\n')
    # for i, ss in enumerate(ptree.xpath('//a/@href')):
    #     if not ss.startswith('ascii-etd.php?id='): continue
    #     fout = open(fname+'.dat', 'w')
    #     dat_url = '%s%s' % (hroot,ss)
    #     dat_url = dat_url.replace(' ','')
    #     furl = urlopen(dat_url)
    #     fout.write(furl.read().decode('ascii','ignore'))
    #     fout.close()
    #     print( fname, ' completed...')
    sleep(0.5)
fdat.close()   

