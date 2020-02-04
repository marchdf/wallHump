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
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import yaml

# ========================================================================
#
# Some defaults variables
#
# ========================================================================
plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Times')
cmap_med = ['#F15A60', '#7AC36A', '#5A9BD4', '#FAA75B',
            '#9E67AB', '#CE7058', '#D77FB4', '#737373']
cmap = ['#EE2E2F', '#008C48', '#185AA9', '#F47D23',
        '#662C91', '#A21D21', '#B43894', '#010202']
dashseq = [(None, None), [10, 5], [10, 4, 3, 4], [
    3, 3], [10, 4, 3, 4, 3, 4], [3, 3], [3, 3]]
markertype = ['s', 'd', 'o', 'p', 'h']


# ========================================================================
#
# Function definitions
#
# ========================================================================


# ========================================================================
#
# Main
#
# ========================================================================
if __name__ == '__main__':

    # ========================================================================
    # Parse arguments
    parser = argparse.ArgumentParser(
        description='A simple plot tool')
    parser.add_argument(
        '-s', '--show', help='Show the plots', action='store_true')
    args = parser.parse_args()

    # ======================================================================
    # Nalu data stuff
    nasa_cp_shift = -0.065

    # ======================================================================
    # NASA CFL3D

    # Cp
    fname = os.path.join(os.path.abspath('../nasa_data'), 'cp_cfl3d.dat')
    df = pd.read_csv(fname, comment='#')

    plt.figure(0)
    p = plt.plot(df['x'],
                 -df['cp'],
                 color=cmap[0],
                 lw=2,
                 label='CFL3D')

    # Cf
    fname = os.path.join(os.path.abspath('../nasa_data'), 'cf_cfl3d.dat')
    df = pd.read_csv(fname, comment='#')

    plt.figure(1)
    p = plt.plot(df['x'],
                 np.fabs(df['cf']),
                 color=cmap[0],
                 lw=2,
                 label='CFL3D')

    # ======================================================================
    # Nalu output
    fname = 'wall_vars.dat'
    df = pd.read_csv(fname)

    # First line
    subdf = df[(df.nx == 103) & (df.order == 1)]
    plt.figure(0)
    p = plt.plot(subdf['x'],
                 -subdf['cp'] - nasa_cp_shift,
                 color=cmap[1],
                 lw=2,
                 label='Nalu p=1')
    p[0].set_dashes(dashseq[1])

    plt.figure(1)
    p = plt.plot(subdf['x'],
                 subdf['cf'],
                 color=cmap[1],
                 lw=2,
                 label='Nalu p=1')
    p[0].set_dashes(dashseq[1])

    # Second line
    subdf = df[(df.nx == 103) & (df.order == 2)]
    subdf = subdf.reset_index()
    subdf = subdf.iloc[::2, :]

    plt.figure(0)
    p = plt.plot(subdf['x'],
                 -subdf['cp'] - nasa_cp_shift,
                 color=cmap[2],
                 lw=2,
                 label='Nalu p=2')
    p[0].set_dashes(dashseq[2])

    plt.figure(1)
    p = plt.plot(subdf['x'],
                 subdf['cf'],
                 color=cmap[2],
                 lw=2,
                 label='Nalu p=2')
    p[0].set_dashes(dashseq[2])

    # ======================================================================
    # NASA experiment

    # Cp
    fname = os.path.join(os.path.abspath('../nasa_data'), 'cp_exp.dat')
    df = pd.read_csv(fname, comment='#')

    plt.figure(0)
    p = plt.plot(df['x'],
                 -df['cp'],
                 ls='None',
                 marker=markertype[0],
                 mec=cmap[-1],
                 mfc=cmap[-1],
                 ms=6,
                 zorder=0,
                 label='Exp.')

    # Cf
    fname = os.path.join(os.path.abspath('../nasa_data'), 'cf_exp.dat')
    df = pd.read_csv(fname, comment='#')

    plt.figure(1)
    p = plt.errorbar(df['x'],
                     np.fabs(df['cf']),
                     yerr=df['error'],
                     fmt='o',
                     capsize=3,
                     color=cmap[-1],
                     marker=markertype[0],
                     mec=cmap[-1],
                     mfc=cmap[-1],
                     ms=6,
                     zorder=0,
                     label='Exp.')

    # ======================================================================
    # Format the plots
    plt.figure(0)
    ax = plt.gca()
    plt.xlabel(r"$x$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$C_p$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    ax.set_xlim([-0.5, 2])
    ax.set_ylim([-0.4, 1])
    legend = ax.legend(loc='best')
    plt.tight_layout()
    # plt.savefig('cp.pdf', format='pdf')
    plt.savefig('cp.png', format='png')

    plt.figure(1)
    ax = plt.gca()
    plt.xlabel(r"$x$", fontsize=22, fontweight='bold')
    plt.ylabel(r"$|C_f|$", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    ax.set_xlim([-0.5, 2])
    ax.set_ylim([-0.004, 0.008])
    legend = ax.legend(loc='best')
    plt.tight_layout()
    # plt.savefig('cf.pdf', format='pdf')
    plt.savefig('cf.png', format='png')
