#!/usr/bin/env python3

# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import sys
import os
import glob
import numpy as np
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import yaml

# ========================================================================
#
# Some defaults variables
#
# ========================================================================
plt.rc("text", usetex=True)
plt.rc("font", family="serif", serif="Times")
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
def parse_ic(fname):
    """Parse the Nalu yaml input file for the initial conditions"""
    with open(fname, "r") as stream:
        try:
            dat = yaml.load(stream)
            u0 = float(
                dat["realms"][0]["initial_conditions"][0]["value"]["velocity"][0]
            )
            rho0 = float(
                dat["realms"][0]["material_properties"]["specifications"][0]["value"]
            )
            mu = float(
                dat["realms"][0]["material_properties"]["specifications"][1]["value"]
            )

            return u0, rho0, mu

        except yaml.YAMLError as exc:
            print(exc)


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
    args = parser.parse_args()

    # ======================================================================
    # Nalu data stuff
    porder = 2
    # fdir = os.path.abspath('103x28')
    # fdir = os.path.abspath('205x55')
    # fdir = os.path.abspath('409x109')
    fdir = os.path.abspath("817x217")
    # fdir = os.path.abspath('1633x433')
    rdir = os.path.abspath(os.path.join(fdir, "results_p{0:d}".format(porder)))
    yname = os.path.join(fdir, "wallHump_p{0:d}.i".format(porder))
    u0, rho0, mu = parse_ic(yname)
    chord = 420
    nasa_cp_shift = -0.015

    # ======================================================================
    # NASA experiment
    fname = os.path.join(os.path.abspath("../nasa_data"), "profiles_exp.dat")
    df = pd.read_csv(fname, comment="#")

    xs = [0.65, 0.8, 1.0, 1.1, 1.2, 1.3, -2.14]
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
            label="Exp. x={0:.2f}".format(x),
        )

        plt.figure(cnt)
        cnt += 1
        p = plt.plot(
            subdf["uv"],
            subdf["y"],
            ls="None",
            lw=1,
            color=cmap[-1],
            marker=markertype[0],
            mec=cmap[-1],
            mfc=cmap[-1],
            ms=6,
            label="Exp. x={0:.2f}".format(x),
        )

    # ======================================================================
    # NASA CFL3D
    fname = os.path.join(os.path.abspath("../nasa_data"), "profiles_cfl3d.dat")
    df = pd.read_csv(fname, comment="#")

    xs = [0.65, 0.8, 1.0, 1.1, 1.2, 1.3, -2.14]
    cnt = 0
    for k, x in enumerate(xs):
        subdf = df[np.fabs(df["x"] - x) < 1e-3]
        plt.figure(cnt)
        cnt += 1
        p = plt.plot(
            subdf["u"],
            subdf["y"],
            ls="-",
            lw=2,
            color=cmap[0],
            label="CFL3D x={0:.2f}".format(x),
        )

        plt.figure(cnt)
        cnt += 1
        p = plt.plot(
            subdf["uv"],
            subdf["y"],
            ls="-",
            lw=2,
            color=cmap[0],
            label="CFL3D x={0:.2f}".format(x),
        )

    # ======================================================================
    # Nalu output
    fname = "profiles.dat"
    df = pd.read_csv(fname)

    # First line
    subdf = df[(df.nx == 205) & (df.order == 1)]

    xs = [0.65, 0.8, 1.0, 1.1, 1.2, 1.3, -2.14]
    cnt = 0
    for k, x in enumerate(xs):

        ssdf = subdf[np.fabs(subdf["x"] - x) < 1e-3]

        plt.figure(cnt)
        cnt += 1
        p = plt.plot(
            ssdf["u"],
            ssdf["y"],
            ls="-",
            lw=2,
            color=cmap[1],
            label="Nalu x={0:.2f}".format(x),
        )
        p[0].set_dashes(dashseq[1])

        plt.figure(cnt)
        cnt += 1
        p = plt.plot(
            ssdf["uv"],
            ssdf["y"],
            ls="-",
            lw=2,
            color=cmap[1],
            label="Nalu x={0:.2f}".format(x),
        )
        p[0].set_dashes(dashseq[1])

    # First line
    subdf = df[(df.nx == 205) & (df.order == 2)]

    xs = [0.65, 0.8, 1.0, 1.1, 1.2, 1.3, -2.14]
    cnt = 0
    for k, x in enumerate(xs):

        ssdf = subdf[np.fabs(subdf["x"] - x) < 1e-3]

        plt.figure(cnt)
        cnt += 1
        p = plt.plot(
            ssdf["u"],
            ssdf["y"],
            ls="-",
            lw=2,
            color=cmap[2],
            label="Nalu x={0:.2f}".format(x),
        )
        p[0].set_dashes(dashseq[2])

        plt.figure(cnt)
        cnt += 1
        p = plt.plot(
            ssdf["uv"],
            ssdf["y"],
            ls="-",
            lw=2,
            color=cmap[2],
            label="Nalu x={0:.2f}".format(x),
        )
        p[0].set_dashes(dashseq[2])

    # ======================================================================
    # Format the plots
    plt.figure(0)
    ax = plt.gca()
    plt.xlabel(r"$u / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.5, 1.5])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('u_0.pdf', format='pdf')
    plt.savefig("u_0.png", format="png")

    plt.figure(1)
    ax = plt.gca()
    plt.xlabel(r"$u'v' / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.005, 0])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('uv_0.pdf', format='pdf')
    plt.savefig("uv_0.png", format="png")

    plt.figure(2)
    ax = plt.gca()
    plt.xlabel(r"$u / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.5, 1.5])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('u_1.pdf', format='pdf')
    plt.savefig("u_1.png", format="png")

    plt.figure(3)
    ax = plt.gca()
    plt.xlabel(r"$u'v' / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.04, 0])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('uv_1.pdf', format='pdf')
    plt.savefig("uv_1.png", format="png")

    plt.figure(4)
    ax = plt.gca()
    plt.xlabel(r"$u / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.5, 1.5])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('u_2.pdf', format='pdf')
    plt.savefig("u_2.png", format="png")

    plt.figure(5)
    ax = plt.gca()
    plt.xlabel(r"$u'v' / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.04, 0])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('uv_2.pdf', format='pdf')
    plt.savefig("uv_2.png", format="png")

    plt.figure(6)
    ax = plt.gca()
    plt.xlabel(r"$u / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.5, 1.5])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('u_3.pdf', format='pdf')
    plt.savefig("u_3.png", format="png")

    plt.figure(7)
    ax = plt.gca()
    plt.xlabel(r"$u'v' / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.04, 0])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('uv_3.pdf', format='pdf')
    plt.savefig("uv_3.png", format="png")

    plt.figure(8)
    ax = plt.gca()
    plt.xlabel(r"$u / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.5, 1.5])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('u_4.pdf', format='pdf')
    plt.savefig("u_4.png", format="png")

    plt.figure(9)
    ax = plt.gca()
    plt.xlabel(r"$u'v' / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.04, 0])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('uv_4.pdf', format='pdf')
    plt.savefig("uv_4.png", format="png")

    plt.figure(10)
    ax = plt.gca()
    plt.xlabel(r"$u / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.5, 1.5])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('u_5.pdf', format='pdf')
    plt.savefig("u_5.png", format="png")

    plt.figure(11)
    ax = plt.gca()
    plt.xlabel(r"$u'v' / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.04, 0])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('uv_5.pdf', format='pdf')
    plt.savefig("uv_5.png", format="png")

    plt.figure(12)
    ax = plt.gca()
    plt.xlabel(r"$u / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.5, 1.5])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('u_6.pdf', format='pdf')
    plt.savefig("u_6.png", format="png")

    plt.figure(13)
    ax = plt.gca()
    plt.xlabel(r"$u'v' / u_r$", fontsize=22, fontweight="bold")
    plt.ylabel(r"$y$", fontsize=22, fontweight="bold")
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight="bold")
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight="bold")
    ax.set_xlim([-0.04, 0])
    ax.set_ylim([0, 0.2])
    legend = ax.legend(loc="best")
    plt.tight_layout()
    # plt.savefig('uv_6.pdf', format='pdf')
    plt.savefig("uv_6.png", format="png")
