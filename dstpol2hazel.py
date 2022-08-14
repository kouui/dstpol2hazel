#!/home/kouui/miniconda3/envs/hazel2.3.6/bin/python



##-----------------------------------------------------------------------------------------------------
## History:
## 2022.08.10    k.u.    init
## 2022.08.11    k.u.    save normalization dataset to obsdata.h5
## 2022.08.12    k.u.    calculate solar radius in arcsec using "DATE_OB2" in fits header 
##                       able to read r,p,incli from dst struct in *.sav file
##                       added ibias, ibias0, conts, clvs to obsdata.h5
## 2022.08.13    k.u.    ap['yc'][0] - int(abs(ac['ddy'][0])) to get the correct continuum position
##-----------------------------------------------------------------------------------------------------



##-----------------------------------------------------------------------------------------------------
## example :
## $ dstpol2hazel /nwork/ichimoto/DST/sp/20220420_He/ar.2.2/sav.sd /nwork/ichimoto/DST/sp/20220420_He/cal/wl_He.sav --apfile /nwork/ichimoto/DST/sp/20220420_He/cal/ap_He.sav --p 37:27 --r 12:39 --i 52:46 --outdir ./output
## $ cd output
## $ mpiexec -n 8 python script.py
##-----------------------------------------------------------------------------------------------------


##-----------------------------------------------------------------------------------------------------
## TODO:
## [1].  [2022.08.10] dst_st in *.sav
##    --> [2022.08.11] added by ichimoto. need to be implemented in this program
##    --> [2022.08.12] implemented
## (2).  [2022.08.10] normalizaiton (currently we use method 3)
##    --> [2022.08.11] method 1 : fitting with photosphere
##                     method 2 : allen * si_continuum to create boundary instead of 
##                                normalizing original profile with allen
##        [2022.08.13] method 3 : keep boundary as 1,0,0,0. and normalize iquv with 
##                                continuum(scan,slit)/clv
##     2.0 [2022.08.10] norm intensity should be provided?
##     2.1 [2022.08.10] how to apply i_bias?
##     2.2 [2022.08.10] since we normalize with disc center continuum, 
##         how to do in sunspot (slit wise?)
##     2.3 [2022.08.10] how to normalize prominence profile?
##         --> [2022.08.11] normalize with peak intensity
##     2.4 [2022.08.11] save normalziation factor to obsdata.h5
##         --> [2022.08.11] normalization dataset are added to obsdata.h5
## [3]. [2022.08.10] variable solar radius in arcsec (not fixed to 960")
##    --> [2022.08.11] reference: https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/pb0r.pro
##    --> [2022.08.12] implemented
## [4]. [2022.08.10] flattern to 1d array of pixels --> 180 ambiguity? potential field?
##    -->  [2022.08.11] with the inversion result, we solve the 180 ambiguity (potential field) by ourselves (in IDL)
## [5]. [2022.08.11] output as fits file?
##    -->  [2022.08.11] will not implement, working with .h5 file is fine
##-----------------------------------------------------------------------------------------------------

import numpy as np
from numpy import arcsin as asin
from numpy import pi as PI
from numpy import cos
from scipy.io import readsav
import h5py
import os, shutil
from glob import glob
import argparse


#----------------------------------------------------------------------
# setup logger
#----------------------------------------------------------------------
import logging as log
NAME = 'hazel2'
LOGGER_TERM = log.getLogger(NAME)
LOGHDL_TERM = log.StreamHandler()
LOGHDL_TERM.setFormatter( log.Formatter(f'[{NAME} : %(asctime)s] %(message)s',datefmt='%Y/%m/%d %H:%M:%S %p') )
LOGGER_TERM.addHandler(LOGHDL_TERM)
LOGGER_TERM.setLevel(log.INFO)
del LOGHDL_TERM, NAME, log

#----------------------------------------------------------------------
# setup script arguments
#----------------------------------------------------------------------

parser = argparse.ArgumentParser(prog="dstpol2hazel", description=f'given calibrated dstpol stokes .sav files, make *.h5 file and configuration files for hazel inversion', add_help=True)
parser.add_argument('savdir', type=str, metavar='SAVDIR',
                    help='path to *.sav directory')
parser.add_argument('wlfile', type=str, metavar='WLFILE',
                    help='path to wavelength *.sav file')
parser.add_argument('--istart', type=int, default=0, metavar='ISTART', 
                    help='start index of *.sav file for inversion, default: 0')
parser.add_argument('--iend', type=int, default=-1, metavar='IEND', 
                    help='end index of *.sav file for inversion, defulat: -1')
parser.add_argument('--outdir', type=str, default='', 
                    help='output directory, default: ./output.[datetime]')
parser.add_argument('--p', type=str, default='',
                    help='dst obs polar angle, in format of "Degree:Minute". must be specified if --usedstst is False')
parser.add_argument('--r', type=str, default='',
                    help='dst obs radius, in format of "ArcMinute:ArcSecond". must be specified if --usedstst is False')
parser.add_argument('--i', type=str, default='',
                    help='dst obs inclination, in format of "Degree:Minute". must be specified if --usedstst is False')
parser.add_argument('--usedstst', action='store_true', default=False, 
                    help='whether to use dst status from fits header, default: False')
parser.add_argument('--wstart', type=float, default=10828.0,
                    help='left bound of the wavelength range in inversion. default: 10828.0 [A] ')
parser.add_argument('--wend', type=float, default=10831.5,
                    help='right bound of the wavelength range in inversion. default: 10831.5 [A] ')
parser.add_argument('--noise', type=float, default=3E-4,
                    help='noise level of QUV profile, default: 3E-4')
parser.add_argument('--height', type=float, default=3.0,
                    help='height of the atmosphere/slab from photosphere, default: 3.0 [arcsec]')
parser.add_argument('--usage', action="store_true", default=False,
                    help="print out an example of how to use this program. e.x. $dstpol2hazel . . --usage ")
parser.add_argument('--apfile', type=str, default='',
                    help='path to alignment parameter file in DSTPOL calibration. must provide either --apfile or --wcont ')
parser.add_argument('--wcont', type=float, default=0.0,
                    help='continuum wavelength [A] to calculate normalization factor. must provide either --apfile or --wcont ')
parser.add_argument('--ibias', action='store_true', default=False, 
                    help='whether to subtract ibias before normalization, default: False')
#parser.add_argument('--normtype', type=int, default=0,
#                    help='0: stokes I --> continuum(slit,scan), clv; boundary --> 1,0,0,0 | 1: stokes I --> max(continuum(slit,scan)); boundary --> clv*continuum(slit,scan)/max(continuum(slit,scan))')
ARGS = parser.parse_args()
del parser, argparse

## customize default values
def make_parser_(args):
    from datetime import datetime
    if args.outdir == '':
        dtstr = datetime.now().strftime("%Y%m%d%H%M%S")
        args.outdir = f"./output.{dtstr}"
        #args.outdir = f"./output"
        
    args.savdir = os.path.abspath(args.savdir)
    args.wlfile = os.path.abspath(args.wlfile)
    args.outdir = os.path.abspath(args.outdir)
    return None

## print input arguments to terminal
def print_args_(args):
    bar = '-'*60
    print(bar)
    print(f"input arguments:")
    for k, v in vars(args).items():
        print(f"{k:>8s}    {v}")
    print(bar)
    
    return None
#----------------------------------------------------------------------
# given utc datetime object, calculate solar radius in arcsec
#----------------------------------------------------------------------
from datetime import datetime, timedelta
from collections import namedtuple
def dt2julian_(utc_datetime):
    if not isinstance(utc_datetime, datetime):
        raise TypeError("input argument must be a datetime.datetime object")
    seconds_oneday = 24*3600
    julian_epoch = datetime(2000, 1, 1, 12) # noon (the epoch name is unrelated)
    j2000_jd = timedelta(2451545) # julian epoch in julian dates
    oneday = timedelta(days=1)

    if utc_datetime.hour < 12:
        julian_int = (utc_datetime.replace(hour=12) - julian_epoch + j2000_jd) //oneday - 1
        julian_frac= ((utc_datetime.hour+12) * 3600 + utc_datetime.minute*60 + utc_datetime.second) / seconds_oneday
    else:
        julian_int = (utc_datetime.replace(hour=12) - julian_epoch + j2000_jd) //oneday
        julian_frac= ((utc_datetime.hour-12) * 3600 + utc_datetime.minute*60 + utc_datetime.second) / seconds_oneday

    Julian = namedtuple('Julian', ('int', 'frac'))
    
    return Julian(julian_int, julian_frac)

## reference : https://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/pb0r.pro
def solarradius_(utc_datetime):
    
    ##: calculate julian day given utc datetime object
    julian = dt2julian_(utc_datetime)
    de = julian.int - 2415020 + julian.frac

    t = de/36525.0

    ##: Form the mean anomalies of Venus(MV),Earth(ME),Mars(MM),Jupiter(MJ)
    ##  and the mean elongation of the Moon from the Sun(D).
    mv = 212.6   + ( (58517.80   * t) % 360.0 )
    me = 358.476 + ( (35999.0498 * t) % 360.0 )
    mm = 319.5   + ( (19139.86   * t) % 360.0 )
    mj = 225.3   + ( ( 3034.69   * t) % 360.0 )
    d = 350.7    + ( (445267.11  * t) % 360.0 )

    ##: Form the geocentric distance(r) [ratio to 1AU] and semi-diameter(sd) [arcmin]
    dtor = PI / 180.
    r = 1.000141 - (0.016748 - 0.0000418*t)*cos(me*dtor) \
      - 0.000140 * cos(2.0*me*dtor)                      \
      + 0.000016 * cos((58.3 + 2.0*mv - 2.0*me)*dtor)    \
      + 0.000005 * cos((209.1 + mv - me)*dtor)           \
      + 0.000005 * cos((253.8 - 2.0*mm + 2.0*me)*dtor)   \
      + 0.000016 * cos(( 89.5 - mj + me)*dtor)           \
      + 0.000009 * cos((357.1 - 2.0*mj + 2.0*me)*dtor)   \
      + 0.000031 * cos(d*dtor)
    
    wcs_rsun = 695508000.00000000   # [m]
    wcs_au   = 149597870691.00000   # [m]
    sd_const = wcs_rsun / wcs_au
    sd = asin(sd_const/r)*10800. / PI  # [arcmin]
    Spos = namedtuple('SolarPos', ('geocentric_distance', 'distance_cm', 'radius_asec'))
    
    return Spos( r, wcs_au*r*100., sd*60. )
#----------------------------------------------------------------------
# functions
#----------------------------------------------------------------------

## plot IQUV heliogram for debug
def check_plot_(sarr):
    
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1,1, dpi=100)
    _, _, nw, _ = sarr.shape
    ax.imshow(sarr[:,0,nw//2,:], cmap="gray")
    plt.show()

    return None

## template directory
TEMPDIR = "/tmp_mnt/home/kouui/miniconda3/envs/hazel2.3.6/lib/hazel/dstpol2hazel.templates/"

## read a single .sav IQUV data file
def read_savfile_(path, args, wr, verbose=False, reverse_wave=False, wl0=None, ap=None):
    usedstst=args.usedstst
    out = {}
    data = readsav(path, verbose=verbose)
    dinfo = data['dinfo']
    s     = data['s']
    i_bias= data['i_bias']
    headers = [s.decode() for s in data['h']]
    obsdt   = [s for s in headers if s.startswith("DATE_OB2")]
    if len(obsdt) == 0: raise ValueError(f"DATE_OB2 not found in header\n{headers}")
    obsdt = obsdt[0].split("'")[1]
    obsdt = obsdt.replace('T',' ')
    obsdt = datetime.strptime(obsdt, '%Y-%m-%d %H:%M:%S.%f')
    obsdt -= timedelta(hours=9)   ## JST --> UTC
    
    if usedstst:
        if 'dst' not in data.keys():
            raise ValueError("these *.sav data do not contain dst structure")
        dstst = data['dst']
        deg_incli = dstst['incli'][0] * 180./PI  # rad to deg
        deg_p     = dstst['p'][0]                # deg
        asec_r    = dstst['p'][0] * 60.          # amin to asec
    else:
        if args.p == '': raise ValueError("Please provide polar angle value p")
        if args.r == '': raise ValueError("Please provide radius value r")
        if args.i == '': raise ValueError("Please provide inclination value i")
        
        dd, mm = [float(w) for w in args.p.split(':')]
        deg_p  = dd + mm/60.
        amin, asec = [float(w) for w in args.r.split(':')]
        asec_r = amin*60. + asec
        
        #deg_incli = dinfo['incli'][0]
        dd, mm = [float(w) for w in args.i.split(':')]
        deg_incli  = dd + mm/60.
        #import pdb; pdb.set_trace()
    del data
    if verbose:
        LOGGER_TERM.info(f"p: {deg_p:.2f} [deg], r: {asec_r:.2f} [asec], incli: {deg_incli:.2f} [deg]")
    

    ##: calculate normalization factor using original s and wl
    if ap is not None:
        ic0, ic1 = ap['yc'][0]
        if ic0 > ic1: ic0, ic1 = ic1, ic0
        ic0 -= int(abs(ap['ddy'][0]))
        ic1 -= int(abs(ap['ddy'][0]))        
    elif wl0 is not None:
        wcont = args.wcont
        if wcont <= wl0.min() or wcont >= wl0.max():
            raise ValueError(f"args.wcont={wcont}, not in range of wl array[{wl0.min(),wl0.max()}]")
        ic0 = np.argmin( np.abs(wl0-wcont) ) 
        ic1 = ic0+3
        if ic1 >= wl0.size-1: ic1 = wl0.size-1
        ic0 -= 2
        if ic0 < 1: ic0 = 1 
    else:
        raise ValueError("must provide ap or wl0 to calculate normalization factor")

    conts = s[0,ic0:ic1,:].mean(axis=0)

    ##: reverse wavelength direction
    if reverse_wave: s = s[:,::-1,:]#np.flip(s, axis=-2)
    ##: select necessary wavelength range
    s = s[:,wr[0]:wr[1]+1,:]

    return {
        's' : s,
        'i_bias' : i_bias,
        'dinfo'  : dinfo,
        'deg_p'  : deg_p,
        'deg_incli': deg_incli,
        'asec_r'  : asec_r,
        'obsdt' : obsdt,
        'cont_intensity': conts
    }

## given left/right bound of wavelength, select index range
def select_wrange_(wl, wstart, wend):
    
    if wstart > wend : raise ValueError(f"wstart = {wstart: .2f} must smaller than wend={wend:.2f}")
    
    if wstart < wl[1]: wstart = wl[1]
    if wend   > wl[-2]: wend  = wl[-2]
    i1 = np.where(wl <= wstart)[0][-1]
    i2 = np.where(wl >= wend)[0][0]

    LOGGER_TERM.info(f"will use wavelength range ({wl[i1]:.2f},{wl[i2]:.2f}) [A], including {i2-i1+1} pixels in wavelength direction")

    return i1, i2

## given r,p,i, calculate LOS (theta. chi, gamma)
def rpi2los_(deg_p, deg_i, asec_r, utc_datetime):
    #sun_r = 960
    sun_r = solarradius_(utc_datetime).radius_asec  #[arcsec]
    if asec_r > sun_r: theta = 90.
    else: theta = asin(asec_r / sun_r) * 180. / PI
    chi   = 0.  ##: chi is fixed to 0.0
    gamma = deg_i - deg_p
    return theta, chi, gamma

##: pack sarr (nscan,4,nwave,nslit) and losarr (nscan,nslit,3) 
def pack_data_(data_list):

    ns,ny,nx = data_list[0]['s'].shape
    nf = len(data_list)
    sarr  = np.empty((nf,ns,ny,nx), dtype='float64')
    ibias = np.empty((nf,), dtype='float64')
    losarr= np.empty((nf,nx,3), dtype='float64')
    conts = np.empty((nf,nx), dtype='float64')


    for i, dl in enumerate(data_list):
        utc_datetime = dl['obsdt']
        sarr[i,:,:,:] = dl.pop('s')
        ibias[i]      = dl['i_bias']
        conts[i,:]    = dl['cont_intensity'][:]
        losarr[i,:,:] = rpi2los_(dl['deg_p'], dl['deg_incli'], dl['asec_r'], utc_datetime)
        print(f"[{i+1:03d}/{nf:03d}]  UTC : {utc_datetime}, (theta,chi,gamma)=({losarr[i,0,0]:.2f},{losarr[i,0,1]:.2f},{losarr[i,0,2]:.2f})[deg]")
    return sarr, losarr, ibias, conts

## make .ini file from template
def make_ini_(outfiles, los, wstart, wend, nwave, height=3.0, wave0=10830):
    tempfile = os.path.join(TEMPDIR, "conf.ini")
    with open(tempfile, 'r') as f:
        content = f.read()
    

    content = content.replace("{$wave-start}", f"{wstart:.2f}")
    content = content.replace("{$wave-end}", f"{wend:.2f}")
    content = content.replace("{$nwave}", f"{nwave}")

    content = content.replace("{$output-file}", f"{outfiles['outfile']}")
    content = content.replace("{$wavelength-file}", f"{outfiles['wavefile']}")
    content = content.replace("{$wavelength-weight-file}", f"{outfiles['wavewtfile']}")
    content = content.replace("{$observation-file}", f"{outfiles['obsfile']}")
    content = content.replace("{$wave0}", f"{wave0}")
    content = content.replace("{$height}", f"{height:.1f}")
    content = content.replace("{$theta}", f"{los[0]:.2f}")
    content = content.replace("{$chi}", f"{los[1]:.2f}")
    content = content.replace("{$gamma}", f"{los[2]:.2f}")


    with open(outfiles['ini'], 'w') as f:
        f.write(content)
    LOGGER_TERM.info(f"saved as : {outfiles['ini']}")
    return 0

## calculate limb darkening dilution factor (A,B,C value of 10830)
def limb_dark_factor_(theta_deg, A=0.97, B=0.18, C=-0.53):
    theta = theta_deg * PI / 180.
    norm_factor = A + B *  np.cos(theta) + C*(1- np.cos(theta) * np.log(1 + 1/np.cos(theta)))
    return norm_factor

#----------------------------------------------------------------------
# main function
#
# step 1. terminal input arguments
# step 2. make filenames
# step 3. read wl
# step 4. read s
# step 5. normalize and save sarr,... to .h5 file
# step 6. create .ini file
# step 7. copy model_chromosphere.1d
# step 8. save .wavelength file
# step 9. save .weights file (wavelength dependent weights)
# step 10. copy script.py
#----------------------------------------------------------------------
def main_(args):

    if args.usage:
        print(
            "$ dstpol2hazel /nwork/ichimoto/DST/sp/20220420_He/ar.2.2/sav.sd /nwork/ichimoto/DST/sp/20220420_He/cal/wl_He.sav --apfile /nwork/ichimoto/DST/sp/20220420_He/cal/ap_He.sav --p 37:27 --r 12:39 --i 52:46 --outdir ./output\n"
            "$ cd output\n"
            "$ mpiexec -n 8 python script.py\n"
        )
        return 0

    if args.wcont < 1.0 and args.apfile == '':
        raise ValueError(f"must provide either --apfile or --wcont to calculate normalization factor")
    if args.apfile != '': ap = readsav(args.apfile)['ap']
    else:                 ap = None

    ##: step 1. terminal input arguments
    make_parser_(args)
    print_args_(args)

    ##: step 2. make filenames
    os.makedirs(args.outdir, exist_ok=True)
    outfiles = {
        "ini" : os.path.join(args.outdir, 'conf.ini'),
        "outfile" : "output.h5",
        "wavefile": "10830.wavelength",
        "wavewtfile" : "10830.weights",
        "obsfile"    : "obsdata.h5"
    }

    ##: step 3. read wl
    data = readsav(args.wlfile)
    wl = data['wl']
    wl0=wl.copy()
    LOGGER_TERM.info(f"original shape of wl array : {wl.shape}")
    del data
    reverse_wave=(wl[10]<wl[0])
    if reverse_wave: wl = wl[::-1]
    i1, i2 = select_wrange_(wl, args.wstart, args.wend)
    wl = wl[i1:i2+1]
    wstart = wl[0]
    wend   = wl[-1]
    nwave  = wl.size
    
    ##: step 4. read s
    savtofind = os.path.join(args.savdir, '*.sav')
    savfiles = glob( savtofind )
    nsavfiles = len(savfiles)
    if nsavfiles == 0:
        LOGGER_TERM.error(f"0 .sav files found in globbing {savtofind}")
        return 1
    del savtofind
    LOGGER_TERM.info(f"found {nsavfiles} .sav files in {args.savdir}")

    if args.iend == -1 : args.iend = nsavfiles
    savfiles = savfiles[args.istart:args.iend]
    nsavfiles = len(savfiles)
    LOGGER_TERM.info(f"{nsavfiles} .sav files for inversion")

    data_list = []
    for i, file in enumerate(savfiles,1):
        out = read_savfile_(file, args, (i1,i2), verbose=False, reverse_wave=reverse_wave, wl0=wl0, ap=ap)
        data_list.append( out )
        if i==1: LOGGER_TERM.info(f"original shape of s array : {data_list[0]['s'].shape}")
        #print(f"\r {i:03d}/{nsavfiles:03d}", end='\r')
    #print()
    del out, i
    LOGGER_TERM.info(f"read all {nsavfiles} files")


    sarr, losarr, ibias, conts = pack_data_(data_list)
    LOGGER_TERM.info(f"(nscan,nstokes,nwavelength,nslit) = {sarr.shape}")

    ##: step 5. normalize and save sarr,... to .h5 file
    outf = os.path.join(args.outdir, outfiles['obsfile'])
    sarr = np.transpose(sarr, axes=(0,3,2,1)) # to (nscan, nslit, nwave, 4)
    nscan, nslit, _, _ = sarr.shape
    norms = np.empty((nscan,nslit),dtype='float64')
    sbtibias = args.ibias
    ibias0 = ibias.copy()
    if not sbtibias: ibias[:] = 0.0
    clvs   = np.empty((nscan,), dtype="float64")
    for i in range(nscan):
        norms[i,:] = (conts[i,:]-ibias[i])
        
        theta = losarr[i, 0, 0]
        clv = limb_dark_factor_(theta)
        clvs[i] = clv

        sarr[i,:,:,0] -= ibias[i]
        sarr[i,:,:,:] /= (norms[i,:].reshape(nslit,1,1)/clv)
        

    bounds = np.zeros((nscan*nslit,nwave,4), dtype="float64")
    bounds[:,:,0] = 1.0
    
    sarr = sarr.reshape(nscan*nslit, nwave,4)
    errs = np.empty((nscan*nslit, nwave, 4), dtype="float64")
    errs[:,:,:] = args.noise
    losarr = losarr.reshape(nscan*nslit, 3)
    

    with h5py.File(outf, 'w') as f:
        dset = f.create_dataset('stokes', data=sarr)
        dset.attrs['description'] = "sarr[i,:,:,0] -= ibias[i];sarr[i,:,:,:] /= (norms[i,:].reshape(nslit,1,1)/clv)"
        dset.attrs['shape'] = '(nscan*nslit, nwave, 4)'
        dset.attrs['nscan'] = nscan
        dset.attrs['nslit'] = nslit
        dset.attrs['nwave'] = nwave
       
        dset = f.create_dataset('sigma', data=errs)
        dset.attrs['shape'] = '(nscan*nslit, nwave, 4)'
       
        dset = f.create_dataset('LOS', data=losarr)
        dset.attrs['shape'] = '(nscan*nslit, 3)'
        
        dset = f.create_dataset('boundary', data=bounds)
        dset.attrs['shape'] = '(nscan*nslit, nwave, 4)'
        dset.attrs['description'] = "boundary all set to (1,0,0,0)"

        dset = f.create_dataset('norm', data=norms)
        dset.attrs['description'] = "norms[i,:] = (conts[i,:]-ibias[i])"
        dset.attrs['shape'] = '(nscan,nslit)'

        dset = f.create_dataset('ibias0', data=ibias0)
        dset.attrs['shape'] = '(nscan,)'
        dset.attrs['description'] = "ibias calculated in calibration"

        dset = f.create_dataset('ibias', data=ibias)
        dset.attrs['shape'] = '(nscan,)'
        dset.attrs['description'] = "if ibias not used in creating inversion data, then this are all filled with 0.0"

        dset = f.create_dataset('conts',data=conts)
        dset.attrs['shape'] = '(nscan,nslit)'
        dset.attrs['description'] = "continuum level"

        dset = f.create_dataset('clvs',data=clvs)
        dset.attrs['shape'] = '(nscan,)'
        dset.attrs['description'] = "center-to-limb variation dilution factor"

        dset = f.create_dataset('wl', data=wl)
        dset.attrs['shape'] = '(nwave,)'
        dset.attrs['description'] = "wavelength"

    LOGGER_TERM.info(f"saved as : {outf}")
    del outf, i, nscan, nslit, theta, errs, bounds, norms

    ##: step 6. create .ini file
    make_ini_(outfiles, losarr[0,:], wstart, wend, nwave, height=args.height)

    ##: step 7. copy model_chromosphere.1d
    outf = os.path.join(args.outdir, "model_chromosphere.1d")
    shutil.copy( 
        os.path.join(TEMPDIR, "model_chromosphere.1d"),
        outf,
    )
    LOGGER_TERM.info(f"saved as : {outf}")
    del outf

    ##: step 8. save .wavelength file
    outf = os.path.join(args.outdir,outfiles['wavefile'])
    np.savetxt(outf, wl, header='lambda')
    # with open(outf, 'w') as f:
    #     f.write("# lambda\n")
    #     for i in range(nwave):
    #         f.write(f"{wl[i]:+.10E}\n")
    LOGGER_TERM.info(f"saved as : {outf}")
    del outf

    ##: step 9. save .weights file (wavelength dependent weights)
    outf = os.path.join(args.outdir,outfiles['wavewtfile'])
    with open(outf, 'w') as f:
        f.write("# WeightI WeightQ WeightU WeightV\n")
        for _ in range(nwave):
            f.write('1.0    1.0   1.0   1.0\n')
    LOGGER_TERM.info(f"saved as : {outf}")
    del outf
    
    ##: step 10. copy script.py
    outf = os.path.join(args.outdir, "script.py")
    with open(os.path.join(TEMPDIR, "script.txt"), 'r') as f:
        content = f.read()
    content = content.replace("{$inifile}", os.path.basename(outfiles['ini']))
    with open(outf, 'w') as f:
        f.write(content)
    LOGGER_TERM.info(f"saved as : {outf}")
    del outf, content
    
    ##: step 11. save debug IQUV spectrohelio image
    #check_plot_(sarr)
    
    LOGGER_TERM.info("completed.")
    return 0



if __name__ == "__main__":
    main_(ARGS)
