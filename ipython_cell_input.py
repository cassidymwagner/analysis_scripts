import sys
sys.path.append("/Users/cassswagner/Library/CloudStorage/OneDrive-OregonStateUniversity/oceans_sf")
import oceans_sf as ocsf

from scipy.io import loadmat
import pandas as pd
import numpy as np
import h5py
# import oceans_sf as ocsf
import collections
from calculate_spectral_fluxes import SpectralFlux
from calculate_sfs import StructureFunctions
from flux_sf_figures import *
from matplotlib import pyplot as plt
import time

import warnings
warnings.filterwarnings("ignore")

import seaborn as sns
sns.set_style(style='white')
sns.set_context('notebook')

plt.rcParams['figure.figsize'] = [12,6]
plt.rcParams['figure.dpi'] = 300
# %config InlineBackend.figure_format = 'svg'

import os
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
plt.rcParams['xtick.bottom'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True

import importlib
# import sys
if 'flux_sf_figures' not in sys.modules:
    import flux_sf_figures            # import module on first run 
else:                      
    import flux_sf_figures      
    importlib.reload(flux_sf_figures)

if 'calculate_sfs' not in sys.modules:
    import calculate_sfs            # import module on first run 
else:                      
    import calculate_sfs      
    importlib.reload(calculate_sfs)

if 'oceans_sf' not in sys.modules:
    import oceans_sf as ocsf            # import module on first run 
else:                      
    import oceans_sf as ocsf     
    importlib.reload(ocsf)

# Read in the data
# filename = "/Users/cassswagner/Library/CloudStorage/OneDrive-OregonStateUniversity/oceans-research/structure-function-turbulence/data_analysis/2layer_128.jld2"
filename = "../GFD_Sims_Output/Data/SurfaceQG_extended_dt_0.00075_n_1024_visc_4.0e-21_order_8_kf_50.0_F_1.0e-5.jld2"
# filename = "/Users/cassswagner/Library/CloudStorage/OneDrive-OregonStateUniversity/oceans-research/2D-forced-dissipative/singlelayerqg_forcedbeta.jld2"

f = h5py.File(filename, "r")
grid = f["grid"]
snapshots = f["snapshots"]

# Set some parameters
Lx = grid["Lx"][()]
Ly = grid["Ly"][()]
nx = grid["nx"][()]
ny = grid["ny"][()]
x = np.asarray(grid["x"])
y = np.asarray(grid["y"])

dx = Lx / nx
dy = Ly / ny

# Sort the snapshots
u_snapshots = collections.OrderedDict(
    sorted(snapshots["u"].items(), key=lambda x: int(x[0]))
)
v_snapshots = collections.OrderedDict(
    sorted(snapshots["v"].items(), key=lambda x: int(x[0]))
)
scalar_snapshots = collections.OrderedDict(
    sorted(snapshots["b"].items(), key=lambda x: int(x[0]))
)

def PostProcess(sfs,dict=True):
    mean_SF_velocity_meridional,mean_SF_velocity_zonal,mean_SF_3rd_velocity_meridional,mean_SF_3rd_velocity_zonal = flux_sf_figures.SFMean(sfs,dict)
    xdiffs, xdiffs_3rd = flux_sf_figures.SeperationDifferences(sfs,dict)
    SF_velocity_zonals,SF_velocity_meridionals,SF_3rd_velocity_zonals,SF_3rd_velocity_meridionals = flux_sf_figures.ReformatSF(sfs,dict)
    boot_SF_vz, boot_SF_vm, boot_SF_3rd_vz, boot_SF_3rd_vm = flux_sf_figures.BootstrapSF(SF_velocity_zonals,SF_velocity_meridionals,SF_3rd_velocity_zonals,SF_3rd_velocity_meridionals)
    return(mean_SF_velocity_zonal,mean_SF_velocity_meridional,xdiffs,boot_SF_vz,boot_SF_vm)

def SwathReduction(dir='zonal',dict=True,boundary="Periodic",even=True,
                   cut=(nx/[256,128,64,32,16]).astype(int)):

    processed_data_swath = []

    # sfs_swath = ocsf.advection_velocity(u_[20:50,:], v_[20:50,:], x, y[20:50],boundary="Periodic")

    for i in cut:
        start = time.time()
        u_snapshots_sliced = u_snapshots.copy()

        if dir == 'zonal':
            for j in u_snapshots_sliced:
                u_snapshots_sliced[j] = u_snapshots[j][:i,:]

            v_snapshots_sliced = v_snapshots.copy()

            for j in v_snapshots_sliced:
                v_snapshots_sliced[j] = v_snapshots[j][:i,:]

            scalar_snapshots_sliced = scalar_snapshots.copy()

            for j in scalar_snapshots_sliced:
                scalar_snapshots_sliced[j] = scalar_snapshots[j][:i,:]

            # nx_slice = np.shape(u_snapshots_sliced['1920000'])[0]
            # ny_slice = np.shape(u_snapshots_sliced['1920000'])[1]
            x_slice = x
            y_slice = y[:i]

        elif dir == 'meridional':
            for j in u_snapshots_sliced:
                u_snapshots_sliced[j] = u_snapshots[j][:,:i]

            v_snapshots_sliced = v_snapshots.copy()

            for j in v_snapshots_sliced:
                v_snapshots_sliced[j] = v_snapshots[j][:,:i]

            scalar_snapshots_sliced = scalar_snapshots.copy()

            for j in scalar_snapshots_sliced:
                scalar_snapshots_sliced[j] = scalar_snapshots[j][:,:i]

            # nx_slice = np.shape(u_snapshots_sliced['1920000'])[0]
            # ny_slice = np.shape(u_snapshots_sliced['1920000'])[1]
            x_slice = x[:i]
            y_slice = y
        
        else:
            raise ValueError('Need to specify meridional or zonal as dir, the direction of swath reduction.')


        # end = time.time()
        # print('time for slice at %s slice: %.2f'%(i,end-start))
        # start = time.time()

        sf_sliced = calculate_sfs.StructureFunctions(u_snapshots_sliced,v_snapshots_sliced,scalar_snapshots_sliced,grid,
                                                    x_slice,y_slice,boundary,even)
        print("cut=%s" %i)
        # end = time.time()
        # print('time for SFs at %s slice: %.2f'%(i,end-start))
        # start = time.time()

        meanSFz,meanSFm,xd,SFzboot,SFmboot = PostProcess(sf_sliced,dict)
        
        # end = time.time()
        # print('time for postprocess at %s slice: %.2f'%(i,end-start))    
        
        df = pd.DataFrame(
            {
                "meanSFz": pd.Series(meanSFz),
                "meanSFm": pd.Series(meanSFm),
                "xd": pd.Series(xd),
                "SFzboot": pd.Series(SFzboot),
                "SFmboot": pd.Series(SFmboot)

            }
        ) 

        processed_data_swath.append(df)
    return(processed_data_swath)

processed_data_swath_zonal = SwathReduction(dict=False,boundary=None,even=False)
