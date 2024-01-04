import pandas as pd
from matplotlib import pyplot as plt
import os
import matplotlib as mpl
import seaborn as sns
import pickle
from scipy.stats import bootstrap
import numpy as np

# Set plotting paramaters
# sns.set_style(style="white")
sns.set_context("talk")
sns.set_style("ticks")

plt.rcParams["xtick.bottom"] = True
plt.rcParams["ytick.left"] = True
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True
# plt.rcParams["figure.dpi"] = 100

os.environ["PATH"] = os.environ["PATH"] + ":/Library/TeX/texbin"
mpl.rcParams["text.usetex"] = True


def ReadData(filepath, filename):
    # Read in the structure functions
    # filepath = "/Users/cassswagner/Library/CloudStorage/OneDrive-OregonStateUniversity/oceans-research/structure-function-turbulence/data_analysis/2layer_512.jld2"
    # filename = "2layer_512"
    filepath = "/Users/cassswagner/Library/CloudStorage/OneDrive-OregonStateUniversity/oceans-research/GFD_Sims_Output/Data/SurfaceQG_extended_dt_0.00075_n_1024_visc_4.0e-21_order_8_kf_50.0_F_1.0e-5.jld2"
    filename = "SurfaceQG"
    # filepath = "/Users/cassswagner/Library/CloudStorage/OneDrive-OregonStateUniversity/oceans-research/2D-forced-dissipative/singlelayerqg_forcedbeta.jld2"
    # filename = "SingleLayerQG"
    # filepath = "/Users/cassswagner/Library/CloudStorage/OneDrive-OregonStateUniversity/oceans-research/2d-forced-dissipative/Output/Data/Isotropic2D_n_256_drag_0.04_order_-2_visc_3.0e-19_order_8_kf_50.0_F_0.01.jld2"
    # filename = "Isotropic2D"

    with open("%s_SFs.pickle" % filepath, "rb") as f:
        appended_data = pickle.load(f)

    with open("%s_spectral_fluxes.pickle" % filepath, "rb") as f:
        spectral_data = pickle.load(f)

    return (appended_data, spectral_data)


def SFMean(appended_data, dict=True, trad=True, scalar=True):

    if dict == True:
        # Take the mean of the structure functions
        mean_SF_velocity_meridional = sum(
            d["SFv_adv"]["SF_meridional"] for d in appended_data
        ) / len(appended_data)
        mean_SF_velocity_zonal = sum(
            d["SFv_adv"]["SF_zonal"] for d in appended_data
        ) / len(appended_data)

        if trad == True:
            mean_SF_3rd_velocity_meridional = np.mean(
                appended_data["SFv_3rd"]["SF_meridional"]
            )
            mean_SF_3rd_velocity_zonal = np.mean(appended_data["SFv_3rd"]["SF_zonal"])

        if scalar == True:
            mean_SF_scalar_meridional = np.mean(appended_data["SFs_adv"]["SF_meridional"])
            mean_SF_scalar_zonal = np.mean(appended_data["SFs_adv"]["SF_zonal"])

        # mean_SF_3rd_velocity_meridional = sum(
        #     d["SFv_3rd"]["SF_meridional"] for d in appended_data
        # ) / len(appended_data)
        # mean_SF_3rd_velocity_zonal = sum(
        #     d["SFv_3rd"]["SF_zonal"] for d in appended_data
        # ) / len(appended_data)

        # mean_SF_scalar_meridional = sum(
        #     d["SFs_adv"]["SF_meridional"] for d in appended_data
        # ) / len(appended_data)
        # mean_SF_scalar_zonal = sum(d["SFs_adv"]["SF_zonal"] for d in appended_data) / len(
        #     appended_data
        # )

        # mean_SF_3rd_scalar_meridional = sum(
        #     d["SFs_3rd"]["SF_meridional"] for d in appended_data
        # ) / len(appended_data)
        # mean_SF_3rd_scalar_zonal = sum(d["SFs_3rd"]["SF_zonal"] for d in appended_data) / len(
        #     appended_data
        # )
    else:
        mean_SF_velocity_meridional = appended_data["SFv_adv"]["SF_meridional"].mean(
            axis=0
        )
        mean_SF_velocity_zonal = appended_data["SFv_adv"]["SF_zonal"].mean(axis=0)
        
        if trad == True:
            mean_SF_3rd_velocity_meridional = appended_data["SFv_3rd"][
                "SF_meridional"
            ].mean(axis=0)
            mean_SF_3rd_velocity_zonal = appended_data["SFv_3rd"]["SF_meridional"].mean(
                axis=0
            )

        if scalar == True:
            mean_SF_scalar_meridional = appended_data["SFs_adv"]["SF_meridional"].mean(axis=0)
            mean_SF_scalar_zonal = appended_data["SFs_adv"]["SF_zonal"].mean(axis=0)

    if (trad == True) and (scalar == True):
        return (
            mean_SF_velocity_meridional,
            mean_SF_velocity_zonal,
            mean_SF_3rd_velocity_meridional,
            mean_SF_3rd_velocity_zonal,
            mean_SF_scalar_meridional,
            mean_SF_scalar_zonal,
        )

    elif (trad == True) and (scalar != True):
        
        return (
            mean_SF_velocity_meridional,
            mean_SF_velocity_zonal,
            mean_SF_3rd_velocity_meridional,
            mean_SF_3rd_velocity_zonal,
        )

    elif (trad != True) and (scalar == True):
        return (
            mean_SF_velocity_meridional,
            mean_SF_velocity_zonal,
            mean_SF_scalar_meridional,
            mean_SF_scalar_zonal,
        )
    
    else:
        return(mean_SF_velocity_meridional,mean_SF_velocity_zonal)


def FluxMean(spectral_data):

    # Take the mean of the spectral fluxes
    mean_energy_flux = sum(d["KE_Flux"] for d in spectral_data) / len(spectral_data)
    # mean_scalar_flux = sum(d["scalar_Flux"] for d in spectral_data) / len(spectral_data)
    mean_kr = sum(d["kr"] for d in spectral_data) / len(spectral_data)

    return (mean_energy_flux, mean_kr)


def SeperationDifferences(appended_data, dict=True,trad=True):

    if dict == True:
        # Get the separation distances for structure functions
        xdiffs = appended_data[0]["SFv_adv"]["x-diffs"]
        ydiffs = appended_data[0]["SFv_adv"]["y-diffs"]

        if trad == True:

            xdiffs_3rd = appended_data[0]["SFv_3rd"]["x-diffs"]
            ydiffs_3rd = appended_data[0]["SFv_3rd"][""]
    else:
        xdiffs = appended_data["SFv_adv"]["x-diffs"].mean(axis=0)
        ydiffs = appended_data["SFv_adv"]["y-diffs"].mean(axis=0)

        if trad ==  True:
            xdiffs_3rd = appended_data["SFv_3rd"]["x-diffs"].mean(axis=0)
            ydiffs_3rd = appended_data["SFv_3rd"]["y-diffs"].mean(axis=0)

    if trad == True:
        return (xdiffs, xdiffs_3rd)
    
    else:
        return(xdiffs)


def ReformatSF(appended_data, dict=True,trad=True):

    # Reformat the data
    SF_velocity_zonals = []
    SF_velocity_meridionals = []
    SF_3rd_velocity_zonals = []
    SF_3rd_velocity_meridionals = []

    # SF_scalar_zonals = []
    # SF_scalar_meridionals = []
    # SF_3rd_scalar_zonals = []
    # SF_3rd_scalar_meridionals = []

    if dict == True:
        for d in appended_data:
            SF_velocity_zonals.append(d["SFv_adv"]["SF_zonal"])
            SF_velocity_meridionals.append(d["SFv_adv"]["SF_meridional"])

            if trad == True:

                SF_3rd_velocity_zonals.append(d["SFv_3rd"]["SF_zonal"])
                SF_3rd_velocity_meridionals.append(d["SFv_3rd"]["SF_meridional"])

            # SF_scalar_zonals.append(d["SFs_adv"]["SF_zonal"])
            # SF_scalar_meridionals.append(d["SFs_adv"]["SF_meridional"])
            # SF_3rd_scalar_zonals.append(d["SFs_3rd"]["SF_zonal"])
            # SF_3rd_scalar_meridionals.append(d["SFs_3rd"]["SF_meridional"])

    else:
        for d in range(len(appended_data["SFv_adv"]["SF_zonal"])):
            SF_velocity_zonals.append(appended_data["SFv_adv"]["SF_zonal"][d])
            SF_velocity_meridionals.append(appended_data["SFv_adv"]["SF_meridional"][d])
            
            if trad == True:
                SF_3rd_velocity_zonals.append(appended_data["SFv_3rd"]["SF_zonal"][d])
                SF_3rd_velocity_meridionals.append(
                    appended_data["SFv_3rd"]["SF_meridional"][d]
                )
    if trad == True:
        return (
            SF_velocity_zonals,
            SF_velocity_meridionals,
            SF_3rd_velocity_zonals,
            SF_3rd_velocity_meridionals,
        )

    else:
        return(SF_velocity_zonals,SF_velocity_meridionals)


def ReformatFlux(spectral_data):

    # Fluxes
    energy_fluxes = []
    # scalar_fluxes = []
    # krs = []

    # Fluxes
    for d in spectral_data:
        energy_fluxes.append(d["KE_Flux"])
        # scalar_fluxes.append(d["scalar_Flux"])
        # krs.append(d["kr"])

    return energy_fluxes


def BootstrapSF(
    SF_velocity_zonals,
    SF_velocity_meridionals,
    SF_3rd_velocity_zonals,
    SF_3rd_velocity_meridionals,trad=True
):

    # Bootstrap the snapshots to get quartiles
    boot_SF_vz = bootstrap((SF_velocity_zonals,), np.mean, confidence_level=0.9, axis=0)
    boot_SF_vm = bootstrap(
        (SF_velocity_meridionals,), np.mean, confidence_level=0.9, axis=0
    )
    
    if trad == True:
        boot_SF_3rd_vz = bootstrap(
            (SF_3rd_velocity_zonals,), np.mean, confidence_level=0.9, axis=0
        )
        boot_SF_3rd_vm = bootstrap(
            (SF_3rd_velocity_meridionals,), np.mean, confidence_level=0.9, axis=0
        )

    # boot_SF_sz = bootstrap((SF_scalar_zonals,), np.mean, confidence_level=0.9, axis=0)
    # boot_SF_sm = bootstrap((SF_scalar_meridionals,), np.mean, confidence_level=0.9, axis=0)
    # boot_SF_3rd_sz = bootstrap(
    #     (SF_3rd_scalar_zonals,), np.mean, confidence_level=0.9, axis=0
    # )
    # boot_SF_3rd_sm = bootstrap(
    #     (SF_3rd_scalar_meridionals,), np.mean, confidence_level=0.9, axis=0
    # )

    if trad == True:
        return (boot_SF_vz, boot_SF_vm, boot_SF_3rd_vz, boot_SF_3rd_vm)
    
    else:
        return(boot_SF_vz, boot_SF_vm)


def BootstrapFlux(energy_fluxes):

    # Fluxes
    boot_energy_flux = bootstrap(
        (energy_fluxes,), np.mean, confidence_level=0.9, axis=0
    )
    # boot_scalar_flux = bootstrap((scalar_fluxes,), np.mean, confidence_level=0.9, axis=0)
    # boot_kr = bootstrap((krs,), np.mean, confidence_level=0.9, axis=0)

    return boot_energy_flux


# Define a plotting function for mean SFs
def SF_mean(SFz, SFm, xdiffs, title, ax):

    ax.loglog(
        xdiffs, abs(SFm), color="k", label="Alternate separation direction", alpha=0.3
    )
    ax.loglog(
        xdiffs,
        SFz,
        "x",
        color="k",
        markerfacecolor="white",
        label="Positive structure function",
        alpha=0.7,
    )
    ax.loglog(xdiffs, -SFz, "o", color="k", label="Negative structure function")

    ax.set_ylabel(r"Structure functions (s$^{-3}$)")
    ax.set_title(title)
    ax.set_xlim(xdiffs.min(), xdiffs.max())

    ax.legend()


# Define a plotting function with bootstrapping
def SF_bootstrap_plot(
    SFz, SFm, xdiffs, ydiffs, bootz0, bootz1, bootm0, bootm1, title, label1, label2, ax
):

    diffmin = min(np.nanmin(xdiffs), np.nanmin(ydiffs))
    diffmax = max(np.nanmax(xdiffs), np.nanmax(ydiffs))

    ax.hlines(
        0,
        xmin=diffmin,
        xmax=diffmax,
        color="k",
        linestyle="dashed",
        alpha=0.5,
    )
    try:
        ax.fill_between(
            xdiffs,
            bootz0,
            bootz1,
            color="k",
            alpha=0.3,
            edgecolor=None,
            # label="90\% confidence zonal",
        )
    except:
        pass

    try:
        ax.fill_between(
            ydiffs,
            bootm0,
            bootm1,
            color="tab:blue",
            alpha=0.3,
            edgecolor=None,
            # label="90\% confidence zonal",
        )
    except:
        pass

    ax.plot(xdiffs, SFz, color="k", label=label1)
    ax.plot(ydiffs, SFm, color="tab:blue", label=label2)

    ax.set_xlabel("Separation distance [m]")

    ax.set_ylabel(r"Dissipation rate [m$^2$ s$^{-3}$]")
    # ax.set_title(title)
    ax.set_xscale("log")
    ax.set_xlim(diffmin, diffmax)

    ax.legend(loc="upper right")


# Define a spectral flux plotting function
def SpectralFluxMean(kr, flux, ylabel, ax):

    ax.semilogx(kr, flux, color="k")
    ax.hlines(
        0,
        xmin=kr.min(),
        xmax=kr.max(),
        color="k",
        linestyle="dashed",
        alpha=0.5,
    )
    ax.set_ylabel(ylabel)
    ax.set_xlabel(r"Wavenumber [$m^{-1}$]")
    ax.set_xlim(kr.min(), kr.max())


# Define a spectral flux plotting function with bootstrapping
def SpectralFluxBootstrap(kr, flux, bootflux, ylabel, ax):

    ax.hlines(
        0,
        xmin=kr.min(),
        xmax=kr.max(),
        color="k",
        linestyle="dashed",
        alpha=0.5,
    )

    ax.fill_between(
        kr,
        bootflux.confidence_interval[0],
        bootflux.confidence_interval[1],
        color="k",
        alpha=0.3,
        edgecolor=None,
        # label="90\% confidence",
    )
    ax.semilogx(kr, flux, color="k")  # , label="Mean")

    ax.set_xlabel(r"Wavenumber [m$^{-1}$]")

    ax.set_ylabel(ylabel)
    ax.set_xlim(kr.min(), kr.max())

    # ax.legend(loc="lower left")


def RunAllPlots(
    filename,
    mean_kr,
    mean_energy_flux,
    boot_energy_flux,
    mean_SF_velocity_zonal,
    mean_SF_velocity_meridional,
    boot_SF_vz,
    boot_SF_vm,
    mean_SF_3rd_velocity_zonal,
    mean_SF_3rd_velocity_meridional,
    boot_SF_3rd_vz,
    boot_SF_3rd_vm,
):

    # Spectral energy flux mean plot
    fig, (ax1) = plt.subplots(figsize=(9, 6))
    SpectralFluxMean(mean_kr, mean_energy_flux, r"Energy flux [m$^2$ s$^{-3}$]", ax1)
    plt.savefig(fname="figures/spectral_energyflux_mean_%s.png" % filename, dpi=400)

    # # Spectral scalar flux mean plot
    # fig, (ax1) = plt.subplots(figsize=(9, 6))
    # SpectralFluxMean(mean_kr, mean_scalar_flux, r"Buoyancy flux [?]", ax1)
    # plt.savefig(fname="figures/spectral_scalarflux_mean_%s.png" % filename, dpi=400)

    # Spectral energy flux with bootstrapping
    fig, (ax1) = plt.subplots(figsize=(10, 7))
    SpectralFluxBootstrap(
        mean_kr,
        mean_energy_flux,
        boot_energy_flux,
        r"Energy flux [m$^2$ s$^{-3}$]",
        ax1,
    )
    ax1.tick_params(direction="in", which="both")
    ax1.xaxis.get_ticklocs(minor=True)
    ax1.minorticks_on()
    ax1.set_ylim(-0.125, 0.025)
    plt.savefig(
        fname="figures/spectral_energyflux_bootstrap_%s.png" % filename, dpi=400
    )

    # # Spectral scalar flux with bootstrapping
    # fig, (ax1) = plt.subplots(figsize=(9, 6))
    # SpectralFluxBootstrap(mean_kr, mean_scalar_flux, boot_scalar_flux, r"Buoyancy flux [?]", ax1)
    # plt.savefig(fname="figures/spectral_scalarflux_bootstrap_%s.png" % filename, dpi=400)

    # Advective velocity SF plot with bootstrapping
    fig, (ax1) = plt.subplots(figsize=(10, 7))
    SF_bootstrap_plot(
        mean_SF_velocity_zonal,
        mean_SF_velocity_meridional,
        boot_SF_vz,
        boot_SF_vm,
        "Advective velocity structure functions",
        ax=ax1,
    )
    ax1.tick_params(direction="in", which="both")
    ax1.xaxis.get_ticklocs(minor=True)
    ax1.minorticks_on()
    ax1.set_ylim(-0.025, 0.125)
    plt.savefig(fname="figures/SF_vadv_bootstrap_%s.png" % filename, dpi=400)

    # Traditional 3rd order velocity SF plot with bootstrapping
    fig, (ax1) = plt.subplots(figsize=(9, 6))
    SF_bootstrap_plot(
        mean_SF_3rd_velocity_zonal,
        mean_SF_3rd_velocity_meridional,
        boot_SF_3rd_vz,
        boot_SF_3rd_vm,
        "Traditional 3rd order velocity structure functions",
        ax=ax1,
    )
    plt.savefig(fname="figures/SF_v3rd_bootstrap_%s.png" % filename, dpi=400)

    # Advective velocity SF mean
    fig, (ax1) = plt.subplots(figsize=(9, 6))
    SF_mean(
        mean_SF_velocity_zonal,
        mean_SF_velocity_meridional,
        "Mean advective velocity structure functions",
        ax=ax1,
    )
    plt.savefig(fname="figures/SF_vadv_mean_snapshot_%s.png" % filename, dpi=400)

    # Traditional 3rd order velocity SF mean
    fig, (ax1) = plt.subplots(figsize=(9, 6))
    SF_mean(
        mean_SF_3rd_velocity_zonal,
        mean_SF_3rd_velocity_meridional,
        "Mean traditional 3rd order velocity structure functions",
        ax=ax1,
    )
    plt.savefig(fname="figures/SF_v3rd_mean_snapshot_%s.png" % filename, dpi=400)

    # # Advective scalar SF plot with bootstrapping
    # fig, (ax1) = plt.subplots(figsize=(9, 6))
    # SF_bootstrap_plot(
    #     mean_SF_scalar_zonal,
    #     mean_SF_scalar_meridional,
    #     boot_SF_sz,
    #     boot_SF_sm,
    #     "Advective buoyancy structure functions",
    #     ax=ax1,
    # )
    # plt.savefig(fname="figures/SF_sadv_bootstrap_%s.png" % filename, dpi=400)

    # # Traditional 3rd order scalar SF plot with bootstrapping
    # fig, (ax1) = plt.subplots()
    # SF_bootstrap_plot(
    #     mean_SF_3rd_scalar_zonal,
    #     mean_SF_3rd_scalar_meridional,
    #     boot_SF_3rd_sz,
    #     boot_SF_3rd_sm,
    #     "Traditional 3rd order buoyancy structure functions",
    #     ax=ax1,
    # )
    # plt.savefig(fname="figures/SF_s3rd_bootstrap_%s.png" % filename, dpi=400)

    # # Advective scalar SF mean
    # fig, (ax1) = plt.subplots(figsize=(9, 6))
    # SF_mean(
    #     mean_SF_scalar_zonal,
    #     mean_SF_scalar_meridional,
    #     "Mean advective buoyancy structure functions",
    #     ax=ax1,
    # )
    # plt.savefig(fname="figures/SF_sadv_mean_snapshot_%s.png" % filename, dpi=400)

    # # Traditional 3rd order scalar SF mean
    # fig, (ax1) = plt.subplots()
    # SF_mean(
    #     mean_SF_3rd_scalar_zonal,
    #     mean_SF_3rd_scalar_meridional,
    #     "Mean traditional 3rd order buoyancy structure functions",
    #     ax=ax1,
    # )
    # plt.savefig(fname="figures/SF_s3rd_mean_snapshot_%s.png" % filename, dpi=400)

    # Compare velocity SF mean
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    SF_mean(
        mean_SF_velocity_zonal,
        mean_SF_velocity_meridional,
        title="Advective velocity structure functions",
        ax=ax1,
    )
    SF_mean(
        mean_SF_3rd_velocity_zonal,
        mean_SF_3rd_velocity_meridional,
        title="Traditional 3rd order velocity structure functions",
        ax=ax2,
    )
    fig.text(0.5, 0.03, "Separation distance (m)", ha="center", va="center")
    ax1.set_xlabel("")
    ax2.set_xlabel("")
    ax2.set_ylabel("")
    ax2.get_legend().remove()
    plt.savefig(fname="figures/SF_v_mean_comparison_%s.png" % filename, dpi=400)

    # # Compare scalar SF mean
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    # SF_mean(
    #     mean_SF_scalar_zonal,
    #     mean_SF_scalar_meridional,
    #     title="Advective scalar structure functions",
    #     ax=ax1,
    # )
    # SF_mean(
    #     mean_SF_3rd_scalar_zonal,
    #     mean_SF_3rd_scalar_meridional,
    #     title="Traditional 3rd order scalar structure functions",
    #     ax=ax2,
    # )
    # fig.text(0.5, 0.03, "Separation distance (m)", ha="center", va="center")
    # ax1.set_xlabel("")
    # ax2.set_xlabel("")
    # ax2.set_ylabel("")
    # ax2.get_legend().remove()
    # plt.savefig(fname="figures/SF_s_mean_comparison_%s.png" % filename, dpi=400)

    # Compare velocity SF bootstrap
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    SF_bootstrap_plot(
        mean_SF_velocity_zonal,
        mean_SF_velocity_meridional,
        boot_SF_vz,
        boot_SF_vm,
        "Advective velocity structure functions",
        ax=ax1,
    )
    SF_bootstrap_plot(
        mean_SF_3rd_velocity_zonal,
        mean_SF_3rd_velocity_meridional,
        boot_SF_3rd_vz,
        boot_SF_3rd_vm,
        "Traditional 3rd order velocity structure functions",
        ax=ax2,
    )
    fig.text(0.5, 0.03, "Separation distance (m)", ha="center", va="center")
    ax1.set_ylim(-0.1, 0.3)
    ax2.set_ylim(-0.1, 0.3)
    ax1.set_xlabel("")
    ax2.set_xlabel("")
    ax2.set_ylabel("")
    ax2.get_legend().remove()
    plt.savefig(fname="figures/SF_v_bootstrap_comparison_%s.png" % filename, dpi=400)

    # # Compare scalar SF bootstrap
    # fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    # SF_bootstrap_plot(
    #     mean_SF_scalar_zonal,
    #     mean_SF_scalar_meridional,
    #     boot_SF_sz,
    #     boot_SF_sm,
    #     "Advective scalar structure functions",
    #     ax=ax1,
    # )
    # SF_bootstrap_plot(
    #     mean_SF_3rd_scalar_zonal,
    #     mean_SF_3rd_scalar_meridional,
    #     boot_SF_3rd_sz,
    #     boot_SF_3rd_sm,
    #     "Traditional 3rd order scalar structure functions",
    #     ax=ax2,
    # )
    # fig.text(0.5, 0.03, "Separation distance (m)", ha="center", va="center")
    # ax1.set_xlabel("")
    # ax2.set_xlabel("")
    # ax2.set_ylabel("")
    # ax2.get_legend().remove()
    # plt.savefig(fname="figures/SF_s_bootstrap_comparison_%s.png" % filename, dpi=400)
