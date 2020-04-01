#!/usr/bin/env python3
# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import os
import numpy as np
import pandas as pd
import scipy.interpolate as spi
import vtk


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == "__main__":

    # Parse arguments
    parser = argparse.ArgumentParser(description="A post-processing tool")
    parser.add_argument(
        "-f", "--fdir", help="Folder to post-process", type=str, required=True
    )
    args = parser.parse_args()

    # Setup
    vname = os.path.join(args.fdir, "flow.vtu")
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vname)
    reader.Update()
    output = reader.GetOutput()
    print()
    coord = np.array(output.GetPoints().GetData())
    vel = np.array(output.GetPointData().GetArray("Velocity"))
    df = pd.DataFrame(
        {
            "p": np.array(output.GetPointData().GetArray("Pressure")),
            "u": vel[:, 0],
            "v": vel[:, 1],
            "w": vel[:, 2],
            "x": coord[:, 0],
            "y": coord[:, 1],
            "z": coord[:, 2],
        }
    )
    df.drop_duplicates(subset=["x", "z"], inplace=True)

    xs = [0.65, 0.8, 1.0, 1.1, 1.2, 1.3, -2.14]
    zlo = [0.116101, 0.0245493, 0.00476345, 0, 0, 0, 0]
    zhi = 0.9
    ninterp = 1001

    lst = []
    for k, x in enumerate(xs):
        xline = np.array([x])
        zline = np.linspace(zlo[k], zhi, ninterp)
        u = spi.griddata(
            (df.x, df.z), df.u, (xline[:, None], zline[None, :]), method="linear"
        ).flatten()
        v = spi.griddata(
            (df.x, df.z), df.v, (xline[:, None], zline[None, :]), method="linear"
        ).flatten()
        w = spi.griddata(
            (df.x, df.z), df.w, (xline[:, None], zline[None, :]), method="linear"
        ).flatten()
        sdf = pd.DataFrame(
            {"u": u, "v": v, "w": w, "x": x * np.ones(ninterp), "z": zline}
        )
        lst.append(sdf)

    odf = pd.concat(lst)
    odf.dropna(inplace=True)
    odf.to_csv(os.path.join(args.fdir, "profiles.csv"), index=False)
