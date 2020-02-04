#!/usr/bin/env python3

# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import utilities


# ========================================================================
#
# Some defaults variables
#
# ========================================================================
plt.rc("text", usetex=True)
cmap_med = [
    "#F15A60",
    "#7AC36A",
    "#5A9BD4",
    "#FAA75B",
    "#9E67AB",
    "#CE7058",
    "#D77FB4",
    "#737373",
]
cmap = [
    "#EE2E2F",
    "#008C48",
    "#185AA9",
    "#F47D23",
    "#662C91",
    "#A21D21",
    "#B43894",
    "#010202",
]
dashseq = [
    (None, None),
    [10, 5],
    [10, 4, 3, 4],
    [3, 3],
    [10, 4, 3, 4, 3, 4],
    [3, 3],
    [3, 3],
]
markertype = ["s", "d", "o", "p", "h"]


# ========================================================================
#
# Function definitions
#
# ========================================================================
def momentum_thickness(u, y):
    return np.trapz(u * (1 - u), y)


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == "__main__":

    # ========================================================================
    # Parse arguments
    parser = argparse.ArgumentParser(description="A simple plot tool")
    parser.add_argument("-s", "--show", help="Show the plots", action="store_true")
    parser.add_argument(
        "-f",
        "--fdirs",
        help="Folders to post-process",
        type=str,
        required=True,
        nargs="+",
    )
    args = parser.parse_args()

    # ======================================================================
    # Setup stuff
    nasa_cp_shift = 0.0  # -0.065
    xs = [0.65, 0.8, 1.0, 1.1, 1.2, 1.3, -2.14]

    # ======================================================================
    # NASA experiment
    nasadir = os.path.abspath("nasa_data")
    fname = os.path.join(nasadir, "profiles_exp.dat")
    df = pd.read_csv(fname, comment="#")

    cnt = 0
    for k, x in enumerate(xs):

        subdf = df[df["x"] == x]
        plt.figure(cnt)
        cnt += 1
        p = plt.plot(
            subdf["u"],
            subdf["y"],
            ls="None",
            lw=1,
            color=cmap[-1],
            marker=markertype[0],
            mec=cmap[-1],
            mfc=cmap[-1],
            ms=6,
            label="Exp.",
        )

        # plt.figure(cnt)
        # cnt += 1
        # p = plt.plot(
        #     subdf["uv"],
        #     subdf["y"],
        #     ls="None",
        #     lw=1,
        #     color=cmap[-1],
        #     marker=markertype[0],
        #     mec=cmap[-1],
        #     mfc=cmap[-1],
        #     ms=6,
        #     label="Exp.",
        # )

        # Momentum thickness at inflow
        if x == -2.14:
            theta = momentum_thickness(subdf["u"], subdf["y"])
            print(f"Experimental momentum thickness at x={x} is {theta}")

    # Cp and Cf
    fname = os.path.join(nasadir, "cp_exp.dat")
    df = pd.read_csv(fname, comment="#")

    plt.figure("cp")
    p = plt.plot(
        df["x"],
        -df["cp"],
        ls="None",
        marker=markertype[0],
        mec=cmap[-1],
        mfc=cmap[-1],
        ms=6,
        zorder=0,
        label="Exp.",
    )

    fname = os.path.join(nasadir, "cf_exp.dat")
    df = pd.read_csv(fname, comment="#")

    plt.figure("cf")
    p = plt.errorbar(
        df["x"],
        np.fabs(df["cf"]),
        yerr=df["error"],
        fmt="o",
        capsize=3,
        color=cmap[-1],
        marker=markertype[0],
        mec=cmap[-1],
        mfc=cmap[-1],
        ms=6,
        zorder=0,
        label="Exp.",
    )

    # ======================================================================
    # NASA CFL3D
    fname = os.path.join(nasadir, "profiles_cfl3d.dat")
    df = pd.read_csv(fname, comment="#")

    cnt = 0
    for k, x in enumerate(xs):
        subdf = df[np.fabs(df["x"] - x) < 1e-3]

        plt.figure(cnt)
        cnt += 1
        p = plt.plot(subdf["u"], subdf["y"], ls="-", lw=2, color=cmap[0], label="CFL3D")

        # plt.figure(cnt)
        # cnt += 1
        # p = plt.plot(
        #     subdf["uv"], subdf["y"], ls="-", lw=2, color=cmap[0], label="CFL3D"
        # )

        if x == -2.14:
            theta = momentum_thickness(subdf["u"], subdf["y"])
            print(f"CFL3D momentum thickness at x={x} is {theta}")

    # Cp and Cf
    fname = os.path.join(nasadir, "cp_cfl3d.dat")
    df = pd.read_csv(fname, comment="#")

    plt.figure("cp")
    p = plt.plot(df["x"], -df["cp"], color=cmap[0], lw=2, label="CFL3D")

    fname = os.path.join(nasadir, "cf_cfl3d.dat")
    df = pd.read_csv(fname, comment="#")

    plt.figure("cf")
    p = plt.plot(df["x"], np.fabs(df["cf"]), color=cmap[0], lw=2, label="CFL3D")

    # ======================================================================
    # Nalu output
    for i, fdir in enumerate(args.fdirs):
        rdir = os.path.abspath(os.path.join(fdir, "results"))
        name = os.path.splitext(
            os.path.basename(glob.glob(os.path.join(fdir, "*.yaml"))[0])
        )[0]
        yname = os.path.join(fdir, name + ".yaml")
        u0, rho0, mu = utilities.parse_ic(yname)

        fname = os.path.join(rdir, "profiles.dat")
        try:
            df = pd.read_csv(fname)
            cnt = 0
            for k, x in enumerate(xs):

                ssdf = df[np.fabs(df["x"] - x) < 1e-3].sort_values(by=["z"])

                plt.figure(cnt)
                cnt += 1
                p = plt.plot(
                    ssdf["u"], ssdf["z"], ls="-", lw=2, color=cmap[i + 1], label="Nalu"
                )
                p[0].set_dashes(dashseq[i + 1])

                # plt.figure(cnt)
                # cnt += 1
                # p = plt.plot(ssdf["uv"], ssdf["y"], ls="-", lw=2, color=cmap[1], label="Nalu")
                # p[0].set_dashes(dashseq[1])

                # Momentum thickness at inflow
                if x == -2.14:
                    theta = momentum_thickness(ssdf["u"], ssdf["z"])
                    print(f"Nalu momentum thickness at x={x} is {theta}")

        except pd.errors.EmptyDataError:
            pass

        # Cp and Cf
        fname = os.path.join(rdir, "wall_vars.dat")
        df = pd.read_csv(fname)

        plt.figure("cp")
        p = plt.plot(
            df["x"], -df["cp"] - nasa_cp_shift, color=cmap[i + 1], lw=2, label="Nalu"
        )
        p[0].set_dashes(dashseq[i + 1])

        plt.figure("cf")
        p = plt.plot(df["x"], df["cf"], color=cmap[i + 1], lw=2, label="Nalu")
        p[0].set_dashes(dashseq[i + 1])

    # ======================================================================
    # Format the plots
    with PdfPages("profiles.pdf") as pdf:
        xs = [0.65, 0.8, 1.0, 1.1, 1.2, 1.3, -2.14]
        cnt = 0
        for k, x in enumerate(xs):
            plt.figure(cnt)
            cnt += 1
            ax = plt.gca()
            plt.xlabel(r"$u / u_r$", fontsize=22, fontweight="bold")
            plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
            plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
            plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
            ax.set_xlim([-0.5, 1.5])
            ax.set_ylim([0, 0.2])
            legend = ax.legend(loc="best")
            plt.title("x={0:.2f}".format(x))
            plt.tight_layout()
            pdf.savefig(dpi=300)

            # plt.figure(cnt)
            # cnt += 1
            # ax = plt.gca()
            # plt.xlabel(r"$u'v' / u_r$", fontsize=22, fontweight="bold")
            # plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
            # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
            # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
            # ax.set_xlim([-0.005, 0])
            # ax.set_ylim([0, 0.2])
            # legend = ax.legend(loc="best")
            # plt.title("x={0:.2f}".format(x))
            # plt.tight_layout()
            # pdf.savefig(dpi=300)

        plt.figure("cp")
        ax = plt.gca()
        plt.xlabel(r"$x$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$C_p$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        # ax.set_xlim([-0.5, 2])
        ax.set_ylim([-0.4, 1])
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)

        plt.figure("cf")
        ax = plt.gca()
        plt.xlabel(r"$x$", fontsize=22, fontweight="bold")
        plt.ylabel(r"$|C_f|$", fontsize=22, fontweight="bold")
        plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
        plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
        # ax.set_xlim([-0.5, 2])
        ax.set_ylim([-0.004, 0.008])
        legend = ax.legend(loc="best")
        plt.tight_layout()
        pdf.savefig(dpi=300)
