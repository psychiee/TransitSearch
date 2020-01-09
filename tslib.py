import numpy as np
from glob import glob

dtor = np.pi/180.0000

def readData1():
    '''
    READ TRANSIT DATA from databases
        0 pname = params
        1,2,3 pper = params
        4,5,6 ptt = params
        7,8,9 pt14 = params
        10 pdepth = params
        11,12 sRAs = params sDECs = params
        13 sVs = params
    '''
    flist = glob('exoplanets.*.csv')
    flist.sort()
    dat = np.genfromtxt(flist[-1], names=True, dtype=None, delimiter=',')
    dat = dat[1:]  # skip header
    vv = np.where((dat['PER'] > 0) & (dat['T14'] > 0) & (dat['TT'] > 0))[0]
    dat = dat[vv]  # check valid data

    pname = dat['NAME']
    pper = np.array(dat['PER'], np.float64)
    ptt = np.array(dat['TT'], np.float64)
    pt14 = np.array(dat['T14'], np.float64)
    pperupper, pperlower = dat['PERUPPER'], dat['PERLOWER']
    pttupper, pttlower = dat['TTUPPER'], dat['TTLOWER']
    pt14upper, pt14lower = dat['T14UPPER'], dat['T14LOWER']
    sVs, pdepth = dat['DEPTH'], dat['V']
    sRAs, sDECs = 15.0 * dat['RA'], dat['DEC']

    return pname, pper, pperupper, pperlower, \
           ptt, pttupper, pttlower, pt14, pt14upper, pt14lower, \
           pdepth, sRAs, sDECs, sVs


def readData2():
    '''
	Using database of ETD sites (http://var2.astro.cz/ETD/)

    '''
    flist = glob('exolist-*.txt')
    flist.sort()
    dat1 = np.genfromtxt(flist[-1], usecols=(0,), dtype='U')
    dat2 = np.genfromtxt(flist[-1], usecols=(1, 2, 3, 4, 5, 6, 7))

    pname = dat1
    pper, ptt, pt14 = dat2[:, 0], dat2[:, 1], dat2[:, 2] / (24 * 60.)
    pdepth = dat2[:, 3]

    pperupper = np.zeros_like(pper)
    pperlower = np.zeros_like(pper)
    pttupper = np.zeros_like(pper)
    pttlower = np.zeros_like(pper)
    pt14upper = np.zeros_like(pper)
    pt14lower = np.zeros_like(pper)

    sVs = dat2[:, 4]
    sRAs = dat2[:, 5]
    sDECs = dat2[:, 6]

    return pname, pper, pperupper, pperlower, \
           ptt, pttupper, pttlower, pt14, pt14upper, pt14lower, \
           pdepth, sRAs, sDECs, sVs


def SunRADec(jd):
    #  form time in Julian centuries from 1900.0
    t = (jd - 2415020.0) / 36525.0

    #  form sun's mean longitude
    l = (279.696678 + ((36000.768925 * t) % 360.0)) * 3600.0

    #  allow for ellipticity of the orbit (equation of centre)
    #  using the Earth's mean anomoly ME
    me = 358.475844 + ((35999.049750 * t) % 360.0)
    ellcor = (6910.1 - 17.2 * t) * np.sin(me * dtor) + 72.3 * np.sin(2.0 * me * dtor)
    l = l + ellcor

    # allow for the Venus perturbations using the mean anomaly of Venus MV
    mv = 212.603219 + ((58517.803875 * t) % 360.0)
    vencorr = 4.8 * np.cos((299.1017 + mv - me) * dtor) + \
              5.5 * np.cos((148.3133 + 2.0 * mv - 2.0 * me) * dtor) + \
              2.5 * np.cos((315.9433 + 2.0 * mv - 3.0 * me) * dtor) + \
              1.6 * np.cos((345.2533 + 3.0 * mv - 4.0 * me) * dtor) + \
              1.0 * np.cos((318.15 + 3.0 * mv - 5.0 * me) * dtor)
    l = l + vencorr

    #  Allow for the Mars perturbations using the mean anomaly of Mars MM
    mm = 319.529425 + ((19139.858500 * t) % 360.0)
    marscorr = 2.0 * np.cos((343.8883 - 2.0 * mm + 2.0 * me) * dtor) + \
               1.8 * np.cos((200.4017 - 2.0 * mm + me) * dtor)
    l = l + marscorr

    # Allow for the Jupiter perturbations using the mean anomaly of
    # Jupiter MJ
    mj = 225.328328 + ((3034.6920239 * t) % 360.0)
    jupcorr = 7.2 * np.cos((179.5317 - mj + me) * dtor) + \
              2.6 * np.cos((263.2167 - mj) * dtor) + \
              2.7 * np.cos((87.1450 - 2.0 * mj + 2.0 * me) * dtor) + \
              1.6 * np.cos((109.4933 - 2.0 * mj + me) * dtor)
    l = l + jupcorr

    # Allow for the Moons perturbations using the mean elongation of
    # the Moon from the Sun D
    d = 350.7376814 + ((445267.11422 * t) % 360.0)
    mooncorr = 6.5 * np.sin(d * dtor)
    l = l + mooncorr

    # Allow for long period terms
    longterm = + 6.4 * np.sin((231.19 + 20.20 * t) * dtor)
    l = l + longterm
    l = (l + 2592000.0) % 1296000.0
    # longmed = l/3600.0

    # Allow for Aberration
    l = l - 20.5

    # Allow for Nutation using the longitude of the Moons mean node OMEGA
    omega = 259.183275 - ((1934.142008 * t) % 360.0)
    l = l - 17.2 * np.sin(omega * dtor)

    # Form the True Obliquity
    oblt = 23.452294 - 0.0130125 * t + (9.2 * np.cos(omega * dtor)) / 3600.0

    # Form Right Ascension and Declination
    l = l / 3600.0
    ra = np.arctan2(np.sin(l * dtor) * np.cos(oblt * dtor), np.cos(l * dtor))

    dec = np.arcsin(np.sin(l * dtor) * np.sin(oblt * dtor))

    ra = np.mod(ra / dtor, 360.000)
    dec = dec / dtor

    return ra, dec


def JulDate(year, month, day):  # ,hour,minute,second,timezone):
    '''
    JD at 12:00 in a day
    '''
    year, month, day = int(year), int(month), int(day)

    a = int((14 - month) / 12)
    y = int(year + 4800 - a)
    m = int(month + 12 * a - 3)
    jd = day + int((153 * m + 2) / 5) + 365 * y + int(y / 4) - int(y / 100) + int(y / 400) - 32045

    return jd


def CalDate(JD):
    """Given mjd return calendar date.

    Retrns a tuple (year,month,day,hour,minute,second). The last is a
    floating point number and others are integers. The precision in
    seconds is about 1e-4.

    To convert jd to mjd use jd - 2400000.5. In this module 2400000.5 is
    stored in MJD0.
    """
    # MJD0 = 2400000.5 # 1858 November 17, 00:00:00 hours

    a = int(JD + 0.5)
    # Julian calendar on or before 1582 October 4 and Gregorian calendar
    # afterwards.
    if a < 2299161.0:
        b = 0.0
        c = a + 1524.0
    else:
        b = np.floor((a - 1867216.25) / 36524.25)
        c = a + b - np.floor(b / 4.0) + 1525

    d = np.floor((c - 122.1) / 365.25)
    e = 365 * d + np.floor(d / 4.0)
    f = np.floor((c - e) / 30.6001)

    day = c - e - np.floor(30.6001 * f)
    month = f - 1.0 - 12.0 * np.floor(f / 14.0)
    year = d - 4715.0 - np.floor((7.0 + month) / 10.0)
    # fracofday = JD - np.floor(JD)
    # hours = fracofday * 24.0

    # sign,hour,minute,second = decimal_to_base60(hours)
    return (year, month, day)  # ,int(sign+str(hour)),minute,second)


def GMST(JD):
    # Greenwich mean sideral time from http://aa.usno.navy.mil/faq/docs/GAST.html
    # GMST = 6.697374558 + 0.06570982441908 D_0 + 1.00273790935 H + 0.000026 T2
    # It will be necessary to reduce GMST to the range 0h to 24h.
    # Setting H = 0 in the above formula yields the Greenwich mean sidereal time at 0h UT
    # // GMST= (6.697374558 + 0.06570982441908*(JD-2451545.0) + 1.00273790935*FINDGEN(data)/data*24) mod 24

    # days from 2000 <<< at UT=12 >>>
    D2000 = JD - 2451545.0 - 0.5

    iD2000 = np.array(D2000, dtype=int)
    rD2000 = D2000 - iD2000

    GMST1 = np.mod(6.697374558 + 0.06570982441908 * iD2000 + \
                   1.00273790935 * (rD2000 * 24), 24)
    # GMST2 = np.mod(18.697374558 + 24.06570982441908*D2000, 24)

    # GMST in degree
    return GMST1 * 15.0


def GAST(JD):
    # days from 2000 <<< at UT=12 >>>
    D2000 = JD - 2451545.0 - 0.5
    # calc GAST Greenwich Apparent Sidereal time
    omega = 125.04 - 0.052954 * D2000  # ;the Longitude of the ascending node of the Moon
    L = 280.47 + 0.98565 * D2000  # ;the Mean Longitude of the Sun
    delphi = -0.000319 * np.sin(omega * dtor) - 0.000024 * np.sin(
        2.0 * L * dtor)  # ; the nutation in longitude, is given in hours
    epsilon = 23.4393 - 0.0000004 * D2000  # ; the obliquity
    GAST = GMST(JD) + 15.0 * delphi * np.cos(epsilon * dtor)
    # GAST in degree
    return GAST


def SunRiSetUT(JD, lam, chi):
    RA, DEC = SunRADec(JD)
    # local hour angle (LHA) in degree
    LHAsun = (GAST(JD) - RA) + lam
    # Altitude of Sun
    ALTsun = np.arcsin(np.cos(LHAsun * dtor) * np.cos(DEC * dtor) * np.cos(chi * dtor) + \
                       np.sin(DEC * dtor) * np.sin(chi * dtor)) / dtor
    # Airmass of Sun
    # zsun = 1.0/cos(pi/2.0 - ALTsun*dtor)
    # two solution(ACOS) of LHA in degree
    LHAsun1 = np.arccos(-np.tan(DEC * dtor) * np.tan(chi * dtor)) / dtor
    LHAsun2 = 360.0 - np.arccos(-np.tan(DEC * dtor) * np.tan(chi * dtor)) / dtor

    dUT = np.mod(6.697374558 + 0.06570982441908 * (JD - 2451545.0), 24)
    sunsetGMST = 360.0 + (LHAsun1 - lam) + RA
    sunsetUT = (sunsetGMST - dUT * 15.0) / 1.00273790935
    sunsetUT = np.mod(sunsetUT / 15.0, 24)

    sunriseGMST = (LHAsun2 - lam) + RA
    sunriseUT = (sunriseGMST - dUT * 15.0) / 1.00273790935
    sunriseUT = np.mod(sunriseUT / 15.0, 24)
    return sunriseUT, sunsetUT


def StarAltAzi(JD, lam, chi, RA, DEC):
    LST = GMST(JD) + lam
    HA = np.mod(LST - RA, 360.0)
    # calculate the altitude and azimuth
    sin_ALT = np.sin(DEC * dtor) * np.sin(chi * dtor) + \
              np.cos(DEC * dtor) * np.cos(chi * dtor) * np.cos(HA * dtor)
    ALT = np.arcsin(sin_ALT) / dtor
    cos_AZI = (np.sin(DEC * dtor) - np.sin(ALT * dtor) * np.sin(chi * dtor)) / \
              (np.cos(ALT * dtor) * np.cos(chi * dtor))
    AZI = np.mod(np.arccos(cos_AZI) / dtor, 360)
    hh = np.where(np.sin(HA * dtor) > 0)[0]
    AZI[hh] = 360.0 - AZI[hh]
    return ALT, AZI
