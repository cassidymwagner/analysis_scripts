import pandas as pd
import oceans_sf as ocsf
import h5py
import collections
import numpy as np
import pickle

# Read in the data
# filename = "/Users/cassswagner/Library/CloudStorage/OneDrive-OregonStateUniversity/oceans-research/structure-function-turbulence/data_analysis/2layer_512.jld2"
# filename = "/Users/cassswagner/OneDrive - Oregon State University/oceans-research/GFD_Sims_Output/Data/SurfaceQG_extended_dt_0.00075_n_1024_visc_4.0e-21_order_8_kf_50.0_F_1.0e-5.jld2"
# filename = "/Users/cassswagner/Library/CloudStorage/OneDrive-OregonStateUniversity/oceans-research/2d-forced-dissipative/Output/Data/Isotropic2D_n_256_drag_0.04_order_-2_visc_3.0e-19_order_8_kf_50.0_F_0.01.jld2"

# f = h5py.File(filename, "r")
# grid = f["grid"]
# snapshots = f["snapshots"]


def StructureFunctions(
    u_snapshots,
    v_snapshots,
    s_snapshots,
    grid=None,
    x=None,
    y=None,
    boundary="Periodic",
    even=True,
):

    # Set some parameters
    # Lx = grid["Lx"]
    # Ly = grid["Ly"]
    # nx = grid["nx"]
    # ny = grid["ny"]
    # x = np.asarray(grid["x"])
    # y = np.asarray(grid["y"])

    # dx = Lx[()] / nx

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

    # Calculate the SFs

    appended_data = []

    for i in u_snapshots:
        # Calculate the velocity SFs
        SFv_adv = ocsf.advection_velocity(
            u_snapshots[i], v_snapshots[i], x, y, boundary, even
        )
        SFv_2nd = ocsf.traditional_velocity(
            u_snapshots[i],
            v_snapshots[i],
            x,
            y,
            boundary,
            order=2,
            even=even,
        )
        SFv_3rd = ocsf.traditional_velocity(
            u_snapshots[i], v_snapshots[i], x, y, boundary, order=3, even=even
        )

        # Calculate the scalar SFs
        # SFs_adv = ocsf.advection_scalar(
        #     scalar_snapshots[i], u_snapshots[i], v_snapshots[i], x, y
        # )
        # SFs_2nd = ocsf.traditional_scalar(
        #     scalar_snapshots[i], u_snapshots[i], v_snapshots[i], x, y, order=2
        # )
        # SFs_3rd = ocsf.traditional_scalar(
        #     scalar_snapshots[i], u_snapshots[i], v_snapshots[i], x, y, order=3
        # )
        try:
            df
        except NameError:

            df = pd.DataFrame(
                {
                    "Snapshot": i,
                    "SFv_adv": pd.Series(SFv_adv),
                    "SFv_2nd": pd.Series(SFv_2nd),
                    "SFv_3rd": pd.Series(SFv_3rd),
                    # "SFs_adv": pd.Series(SFs_adv),
                    # "SFs_2nd": pd.Series(SFs_2nd),
                    # "SFs_3rd": pd.Series(SFs_3rd),
                }
            )

        df_new = pd.DataFrame(
            {
                "Snapshot": i,
                "SFv_adv": pd.Series(SFv_adv),
                "SFv_2nd": pd.Series(SFv_2nd),
                "SFv_3rd": pd.Series(SFv_3rd),
                # "SFs_adv": pd.Series(SFs_adv),
                # "SFs_2nd": pd.Series(SFs_2nd),
                # "SFs_3rd": pd.Series(SFs_3rd),
            }
        )

        df = pd.concat([df, df_new])
        # appended_data.append(df)
    return df
    # return appended_data


# # Output the list of SF dataframes
# with open("%s_SFs.pickle" % filename, "wb") as handle:
#     pickle.dump(appended_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
