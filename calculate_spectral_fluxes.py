import pandas as pd
import h5py
import collections
import numpy as np
import pickle
from scipy import fft

# Read in the data
# filename = "/Users/cassswagner/Library/CloudStorage/OneDrive-OregonStateUniversity/oceans-research/structure-function-turbulence/data_analysis/2layer_512.jld2"
# filename = "/Users/cassswagner/OneDrive - Oregon State University/oceans-research/GFD_Sims_Output/Data/SurfaceQG_extended_dt_0.00075_n_1024_visc_4.0e-21_order_8_kf_50.0_F_1.0e-5.jld2"
# filename = "/Users/cassswagner/Library/CloudStorage/OneDrive-OregonStateUniversity/oceans-research/2D-forced-dissipative/twodturb_forced_example.jld2"
# filename = "/Users/cassswagner/Library/CloudStorage/OneDrive-OregonStateUniversity/oceans-research/2d-forced-dissipative/Output/Data/Isotropic2D_n_256_drag_0.04_order_-2_visc_3.0e-19_order_8_kf_50.0_F_0.01.jld2"

# f = h5py.File(filename, "r")
# grid = f["grid"]
# snapshots = f["snapshots"]

def SpectralFlux(u_snapshots,v_snapshots,s_snapshots,grid):

    # Set some parameters
    Lx = grid["Lx"]
    Ly = grid["Ly"]
    nx = grid["nx"]
    ny = grid["ny"]
    x = np.asarray(grid["x"])
    y = np.asarray(grid["y"])

    dx = Lx[()] / nx[()]

# # Sort the snapshots
# u_snapshots = collections.OrderedDict(
#     sorted(snapshots["u"].items(), key=lambda x: int(x[0]))
# )
# v_snapshots = collections.OrderedDict(
#     sorted(snapshots["v"].items(), key=lambda x: int(x[0]))
# )
# scalar_snapshots = collections.OrderedDict(
#     sorted(snapshots["u"].items(), key=lambda x: int(x[0]))
# )


    # Set up some values for the spectral fluxes
    k_int = 1 / dx
    N = len(x)
    k = 2 * np.pi * (k_int / N) * np.linspace(-N / 2, N / 2 - 1, len(x))
    l = k
    kk, ll = np.meshgrid(k, l)


    # Calculate the SFs

    appended_data = []

    for i in u_snapshots:
        u = u_snapshots[i][0][()]
        v = v_snapshots[i][0][()]
        s = s_snapshots[i][0][()]
        # Calculate gradients exactly via spectral calculation
        dsdx = np.real(fft.ifft2(fft.ifftshift(1j * kk * fft.fftshift(fft.fft2(s)))))
        dsdy = np.real(fft.ifft2(fft.ifftshift(1j * ll * fft.fftshift(fft.fft2(s)))))
        dudx = np.real(fft.ifft2(fft.ifftshift(1j * kk * fft.fftshift(fft.fft2(u)))))
        dudy = np.real(fft.ifft2(fft.ifftshift(1j * ll * fft.fftshift(fft.fft2(u)))))
        dvdx = np.real(fft.ifft2(fft.ifftshift(1j * kk * fft.fftshift(fft.fft2(v)))))
        dvdy = np.real(fft.ifft2(fft.ifftshift(1j * ll * fft.fftshift(fft.fft2(v)))))

        # Diagnose Fourier transformed buoyancy advection
        J_s = fft.fft2(u * dsdx + v * dsdy) * dx**2 / (2 * np.pi)
        J_u = fft.fft2(u * dudx + v * dudy) * dx**2 / (2 * np.pi)
        J_v = fft.fft2(u * dvdx + v * dvdy) * dx**2 / (2 * np.pi)

        s_f = fft.fft2(s) * dx**2 / (2 * np.pi)
        u_f = fft.fft2(u) * dx**2 / (2 * np.pi)
        v_f = fft.fft2(v) * dx**2 / (2 * np.pi)

        # Derive the spectral flux divergences from advection and streamfunction
        s_flux_div = fft.fftshift(np.real(np.conj(s_f) * J_s))
        KE_flux_div = fft.fftshift(np.real(np.conj(u_f) * J_u + np.conj(v_f) * J_v))

        # s_spectra = fft.fftshift(np.real(np.conj(s_f) * s_f))
        KE_spectra = fft.fftshift(np.real(np.conj(u_f) * v_f))  # maybe incorrect

        # Set up the wavenumbers
        k_max_mat = max(max(-k[0 : int(len(k) / 2 + 1)]), max(l))
        dk = 2 * np.pi / Lx[()]
        dl = dk
        dkr_mat = np.sqrt(dk**2 + dl**2)
        tmp_kr = np.arange(dkr_mat / 2, k_max_mat + dkr_mat, dkr_mat)
        tmp_kr = np.conj(tmp_kr).T  # get the complex conjugate transpose

        # Calculate fluxes by summing flux divergences over all wavenumbers
        s_Flux = np.zeros(np.shape(tmp_kr))
        KE_Flux = np.zeros(np.shape(tmp_kr))

        for kkk in range(len(tmp_kr)):
            s_Flux[kkk] = (
                -2
                * sum(s_flux_div[kk**2 + ll**2 >= tmp_kr[kkk] ** 2] * dk * dl)
                / ((2 * np.pi) ** 2)
            )
            KE_Flux[kkk] = -sum(
                KE_flux_div[kk**2 + ll**2 >= tmp_kr[kkk] ** 2] * dk * dl
            ) / ((2 * np.pi) ** 2)

        df = pd.DataFrame(
            {
                "scalar_Flux": pd.Series(s_Flux),
                "KE_Flux": pd.Series(KE_Flux),
                "kr": pd.Series(tmp_kr),
                # "KE_spectra": pd.Series(KE_spectra),
            }
        )

        appended_data.append(df)

    return(appended_data)

# # Output the list of SF dataframes
# with open("%s_spectral_fluxes.pickle" % filename, "wb") as handle:
#     pickle.dump(appended_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
