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
from netCDF4 import Dataset, chartostring
import re

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
    # Nalu data
    sname = 'points.dat'
    separation_exp = 0.665
    reattachement_exp = 1.1
    points_df = pd.read_csv(sname)

    print(points_df)

    # # Plot Nalu data
    # plt.figure(0)
    # p = plt.plot(points_df['nx'],
    #              points_df['separation'],
    #              lw=2,
    #              color=cmap[0],
    #              marker=markertype[0],
    #              mec=cmap[0],
    #              mfc=cmap[0],
    #              ms=10,
    #              label='Nalu')

    # plt.figure(1)
    # p = plt.plot(points_df['nx'],
    #              points_df['reattachement'],
    #              lw=2,
    #              color=cmap[0],
    #              marker=markertype[0],
    #              mec=cmap[0],
    #              mfc=cmap[0],
    #              ms=10,
    #              label='Nalu')

    # # Plot experimental data
    # nmin = np.min(points_df.nx)
    # nmax = np.max(points_df.nx)
    # plt.figure(0)
    # p = plt.plot([nmin, nmax],
    #              [separation_exp, separation_exp],
    #              color=cmap[-1],
    #              label='Exp')
    # plt.figure(1)
    # p = plt.plot([nmin, nmax],
    #              [reattachement_exp, reattachement_exp],
    #              color=cmap[-1],
    #              label='Exp')

    # # Format plots
    # plt.figure(0)
    # ax = plt.gca()
    # plt.xlabel(r"$N$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"Separation", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # # plt.savefig('separation.pdf', format='pdf')
    # plt.savefig('separation.png', format='png')

    # plt.figure(1)
    # ax = plt.gca()
    # plt.xlabel(r"$N$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"Reattachement", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # # plt.savefig('reattachement.pdf', format='pdf')
    # plt.savefig('reattachement.png', format='png')
