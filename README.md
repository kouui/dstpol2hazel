# dstpol2hazel
transform calibrated IQUV spectra (IDL .sav file) observed in Hida/DST to the input .h5 data file and configuration file for hazel2(https://github.com/aasensio/hazel2) inversion

hazel2:
- github link : https://github.com/aasensio/hazel2
- documentation : https://aasensio.github.io/hazel2/


## 1. setup environment

in HIDA server, the python environment of hazel is installed in `/home/kouui/miniconda3/envs/hazel2.3.6`, what you need is to add this environment to your path. 

add the following line to your `~/.cshrc` file.
```shell
alias add-hazel-path setenv PATH /home/kouui/miniconda3/envs/hazel2.3.6/bin:${PATH}
```

then every time you log into HIDA server, just run 
```shell
$ add-hazel-path
```
to make hazel2.3.6 environment visible to your current terminal.

to check whether it works well, execute
```shell
$ which python
/home/kouui/miniconda3/envs/hazel2.3.6/bin/python
```
which should print out the python belonging to hazel2.3.6 environment. and also make sure that
```python
$ python
Python 3.6.13 |Anaconda, Inc.| (default, Jun  4 2021, 14:25:59)
[GCC 7.5.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import hazel
>>>
```
`hazel` package could be imported correctly in your python interpreter.


## 2. run hazel in synthesis mode

to run one synthesis simulation using hazel, the simpliest way is as following

```python
import numpy as np
import hazel

##: initialize a model
mod = hazel.Model(working_mode='synthesis', verbose=3)

##: define a spectral called 'spec1' 
##  - 'topology' :  a topology of one atmosphere called 'ch1' (if 'ph1 -> ch1' then is a photosphere+chromosphere, refer to hazel documentation for details)
##  - 'LOS' [theta,chi,gamma] = [0,0,0], the target is located at disc center
##        LOS   : [theta, chi, gamma], in [deg]
##        - theta : angle between LOS and Z(normal to solar surface). 0 : disc center; 90: limb
##        - chi   : keep this value to 0. then X is pointing towards disc center
##        - gamma : angle between X0(vector of X projected onto the plane normal to LOS) and +Q(slit direction in DSTPOL)    
##  - 'Wavelength' : [left bound, right bound, numer of point of wavelength]
##  - 'Boundary condition' : background stokes vector, (1,0,0,0) on disk and (0,0,0,0) off limb. (refer to hazel documentation for details) 
nwave = 150
mod.add_spectral({'Name': 'spec1', 'Wavelength': [10826, 10833, nwave], 'topology': 'ch1',
    'LOS': [0.0,0.0,0.0], 'Boundary condition': [1.0,0.0,0.0,0.0]})

##: and then add the atmosphere 'ch1' to the model
##  - 'Height': height from photosphere in unit of [arcsec] (how it is converted to [km] is still not sure)
##  - 'Line' : related to which atmoc model to use
##  - 'Wavelength' : the wavelength range for simulating this specific 'Line'
mod.add_chromosphere({'Name': 'ch1', 'Spectral region': 'spec1', 'Height': 10.0, 'Line': '10830',
                      'Wavelength': [10826, 10833]})


##: then setup the model with the given parameters
mod.setup()

##: the necessary free parameter for stokes profile generation in a chromosphere is 
##    - Bxyz : magnetic filed, [Gauss]
##    - tau  : optical depth
##    - v    : doppler velocity, [km/s]
##    - delta: turbulent velocity, [km/s]
##    - beta : description not found (amplification factor?)
##    - a    : damping constant
##    - ff   : filling factor
Bx = By = Bz = 10
tau = 0.2
v   = 0.0
delta = 8.0
beta = 1.0
a = 0.1
ff = 1.0
mod.atmospheres['ch1'].set_parameters([Bx,By,Bz,tau,v,delta,beta,a],ff)

##: then we are able to synthesize
mod.synthesize()
##: the result wavelength array is mod.spectrum['spec1'].wavelength_axis[:]
##: the result stokes profile array is mod.spectrum['spec1'].stokes[:,:]

##: to plot the profiles as an example
import matplotlib.pyplot as plt
fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10,10), sharex=True)
axs = axs.flatten()
label = ['I', 'Q', 'U', 'V']
for i in range(4):
    ax = axs[i]
    ax.plot(mod.spectrum['spec1'].wavelength_axis, mod.spectrum['spec1'].stokes[i,:])
    if i>0: ## set a symmetric ylim according to the original ylim
        ylim = np.abs(ax.get_ylim()).max() * 1.05 
        ax.set_ylim(-ylim,+ylim)
    ax.set_xlabel('Wavelength - 10830[$\AA$]')
    ax.set_ylabel('{0}/Ic'.format(label[i]))
plt.show(block=False)

```

## 3. use `dstpol2hazel` to prepare inversion data

the `dstpol2hazel` command line tool is for generating input data and configuration files of hazel inversion given the DSTPOL calibrated IDL `*.sav` files.
which is a symbolic to a python script.
```
$ which dstpol2hazel
/home/kouui/miniconda3/envs/hazel2.3.6/bin/dstpol2hazel
$ ls -l /home/kouui/miniconda3/envs/hazel2.3.6/bin/dstpol2hazel
lrwxrwxrwx. 1 kouui users 67  8月 12 14:44 /home/kouui/miniconda3/envs/hazel2.3.6/bin/dstpol2hazel -> /home/kouui/miniconda3/envs/hazel2.3.6/lib/hazel/dstpol2hazel.v2.py
```
please feel free to check the source code.

to check the meaning of command line arguments
```
$ dstpol2hazel -h
usage: dstpol2hazel [-h] [--istart ISTART] [--iend IEND] [--outdir OUTDIR]
                    [--p P] [--r R] [--i I] [--usedstst] [--wstart WSTART]
                    [--wend WEND] [--noise NOISE] [--height HEIGHT] [--usage]
                    [--apfile APFILE] [--wcont WCONT] [--ibias]
                    SAVDIR WLFILE

given calibrated dstpol stokes .sav files, make *.h5 file and configuration
files for hazel inversion

positional arguments:
  SAVDIR           path to *.sav directory
  WLFILE           path to wavelength *.sav file

optional arguments:
  -h, --help       show this help message and exit
  --istart ISTART  start index of *.sav file for inversion, default: 0
  --iend IEND      end index of *.sav file for inversion, defulat: -1
  --outdir OUTDIR  output directory, default: ./output.[datetime]
  --p P            dst obs polar angle, in format of "Degree:Minute". must be
                   specified if --usedstst is False
  --r R            dst obs radius, in format of "ArcMinute:ArcSecond". must be
                   specified if --usedstst is False
  --i I            dst obs inclination, in format of "Degree:Minute". must be
                   specified if --usedstst is False
  --usedstst       whether to use dst status from fits header, default: False
  --wstart WSTART  left bound of the wavelength range in inversion. default:
                   10828.0 [A]
  --wend WEND      right bound of the wavelength range in inversion. default:
                   10831.5 [A]
  --noise NOISE    noise level of QUV profile, default: 3E-4
  --height HEIGHT  height of the atmosphere/slab from photosphere, default:
                   3.0 [arcsec]
  --usage          print out an example of how to use this program. e.x.
                   $dstpol2hazel . . --usage
  --apfile APFILE  path to alignment parameter file in DSTPOL calibration.
                   must provide either --apfile or --wcont
  --wcont WCONT    continuum wavelength [A] to calculate normalization factor.
                   must provide either --apfile or --wcont
  --ibias          whether to subtract ibias before normalization, default:
                   False
```
where
```
SAVDIR : path to *.sav directory created by DSTPOL calibration program, either /path/to/sav.s or /path/to/sav.sd. must be provided
WLFILE : path to wavelength .sav file created by DSTPOL calibration program, must be provided
--istart ISTART --iend IEND : specify the the range of *.sav file you want to process, if not provided will process all *.sav file found in SAVDIR
--outdir OUTDIR : the output directory, if not provided, then files will be generated in ./output.[datetime]
--p P : dst obs polar angle, in format of "Degree:Minute". must be provided without --usedstst
--r R : dst obs radius, in format of "ArcMinute:ArcSecond". must be provided without --usedstst
--i I : dst obs inclination, in format of "Degree:Minute". must be provided without --usedstst
--usedstst : whether to use dst status from fits header. if not specified, the must provide --p P --i I --r R
--wstart WSTART --wend WEND : this is the wavelength range to cut in observed data that will be used in inversion, defulat is 10828.0 and 10831.5, [Angstrom]
--noise NOISE : the noise level of the stokes profile, default is 3E-4
--height HEIGHT : the height of the slab above photospere, default is 3.0, [arcsec]
--apfile APFIE : the alignment parameter .sav file created by DSTPOL calibration program. we use the 'yc' and 'ddy' parameter in this file to extract continuum intensity. one must provide either --apfile or --wcont. if both are provided, then --apfile will be used
--wcont WCONT : the wavelength to extract continuum intensity if --apfile is not provided, [Angstrom]
--ibias : if specified, then will subtract the `i_bias` value stored in *.sav files before continuum normalization
```

to print out an execution example
```
$ dstpol2hazel . . --usage
[the following lines are the output in the terminal]
$ dstpol2hazel /nwork/ichimoto/DST/sp/20220420_He/ar.2.2/sav.sd /nwork/ichimoto/DST/sp/20220420_He/cal/wl_He.sav --apfile /nwork/ichimoto/DST/sp/20220420_He/cal/ap_He.sav --p 37:27 --r 12:39 --i 52:46 --outdir ./output
$ cd output
$ mpiexec -n 8 python script.py  ## do inversion with 8 cpu cores in parallel
```

the above example has p,i,r mannually specified. to use the p,i,r values stored in `dst` struct in the *.sav file and subtract ibias
```
$ dstpol2hazel /nwork/ichimoto/DST/sp/20220420_He/ar.2.2/sav.sd /nwork/ichimoto/DST/sp/20220420_He/cal/wl_He.sav --apfile /nwork/ichimoto/DST/sp/20220420_He/cal/ap_He.sav --usedstst --ibias --outdir ./output
```

or to do inversion with the 3rd (file index of 2) .sav file only
```
$ dstpol2hazel /nwork/ichimoto/DST/sp/20220420_He/ar.2.2/sav.sd /nwork/ichimoto/DST/sp/20220420_He/cal/wl_He.sav --apfile /nwork/ichimoto/DST/sp/20220420_He/cal/ap_He.sav --usedstst --ibias --outdir ./output --istart 2 --iend 3
```

or to change the inversion wavelength range to [10827.5, 10832.0]
```
$ dstpol2hazel /nwork/ichimoto/DST/sp/20220420_He/ar.2.2/sav.sd /nwork/ichimoto/DST/sp/20220420_He/cal/wl_He.sav --apfile /nwork/ichimoto/DST/sp/20220420_He/cal/ap_He.sav --usedstst --ibias --outdir ./output --wstart 10827.5 --wend 10832.0
```
and so on.

the generated folder `./output` contains files
```
$ cd output
$ ls
10830.wavelength  10830.weights  conf.ini  model_chromosphere.1d  obsdata.h5 script.py
```
where
```
10830.wavelength : the 1d wavelength array
10830.weights    : the I/Q/U/V weights in for each wavelength point
conf.ini         : the configuration file for inversion, mainly including the fitting range, number of node Z(?) direction of each free parameter, how many cycle do perform inversion and what parameter to fit in each cycle.
model_chromosphere.1d : the initial value of the parameters in chromosphere?
script.py        : the python script to run inversion in parallel
obsdata.h5       : the necessary data to perform inversion (and some other data we used during data generation) in HDF5 format
    - necessary data
        - 'stokes' : the observed stoke profile, (nscan*nslit,nwave,4)
            - has attributes 'nscan', 'nslit', 'nwave', which could be use to restore the spatial shape in (nslit,nscan) from the 1d array (nscan*nslit)
        - 'sigma'  : fitting error, (nscan*nslit,nwave,4)
        - 'LOS'    : [theta,chi,gamma], (nscan*nslit,3)
        - 'boundary' : boundary condition, all set to (1,0,0,0), (nscan*nslit,nwave,4)
    - else
        - 'ibias0' : the original ibias data restored from *.sav files, (nscan,)
        - 'ibias'  : the actual ibias used in normalization, equals to 0 if --ibias not specified, (nscan,)
        - 'conts'  : the continuum intensity, which depends on (scan,slit) position, without ibias subtraction, (nscan,nslit)
        - 'clvs'   : the center-to-limb dilution factor used in normalization, (nscan,)
        - 'wl'     : the wavelength array, (nwave)
```

the result of inversion will be stored in `output.h5`, check hazel documentation for details.

## 4. post process

with the data in `obsdata.h5`, to restore observed stokes profile
```python
import h5py

with h5py.File('obsdata.h5','r') as f:
    s = f['stokes'][:,:,:]
    nscan = f['stokes'].attrs['nscan']
    nslit = f['stokes'].attrs['nslit']
    nwave = f['stokes'].attrs['nwave']
    ibias = f['ibias'][:]
    conts = f['conts'][:,:]
    clvs  = f['clvs'][:]


s = s.reshape(nscan,nslit,nwave,4)
s_obs = s * ( ( conts - ibias.reshape(nscan,1) )/clvs.reshape(nscan,1) ).reshape(nscan,nslit,1,1)
s_obs[:,:,:,0] += ibias.reshape(nscan,1,1)
##: then s_obs is the observed stokes 2d spectra

```

to convert inversion result Bxyz to LOS coordinate
```python
import h5py
import numpy as np

with h5py.File('obsdata.h5','r') as f:
    nscan = f['stokes'].attrs['nscan']
    nslit = f['stokes'].attrs['nslit']
    los = f['LOS'][:,:].reshape(nscan,nslit,3)  # (nscan,nslit,3)
with h5py.File('output.h5','r') as f:
    Bx = f['ch1']['Bx'][:]
    By = f['ch1']['By'][:]
    Bz = f['ch1']['Bz'][:]
Bx = Bx.reshape(nscan,nslit)  # (nscan,nslit)
By = By.reshape(nscan,nslit)
Bz = Bz.reshape(nscan,nslit)
theta = np.deg2rad(los[:,:,0])  # (nscan,nslit)
gamma = np.deg2rad(los[:,:,2])  # (nscan,nslit)

B_los = Bz * np.cos( theta ) + Bx * np.sin( theta )
B_plusQ  = By * np.sin( gamma ) + ( Bx * np.cos( theta ) - Bz * np.sin( theta ) ) * np.cos( gamma )  # +Q is in slit direction
B_minusQ = By * np.cos( gamma ) - ( Bx * np.cos( theta ) - Bz * np.sin( theta ) ) * np.sin( gamma )  # -Q is in scan direction

## then have fun with data visualization
...
```

to read .h5 data in IDL, refer to https://www.l3harrisgeospatial.com/docs/hdf5_overview.html
