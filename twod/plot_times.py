#!/usr/bin/env python3

# ========================================================================
#
# Imports
#
# ========================================================================
import argparse
import sys
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

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
    points_df = pd.read_csv(sname)
    points_df['norm_time'] = points_df['walltime'] * points_df['procs'] / \
        (points_df['iterations'] * points_df['dofs'])
    points_df['norm_time'] /= points_df['norm_time'][0]
    print(points_df)
    print(points_df['walltime'] / points_df['iterations'])

    # Normalized time
    plt.figure(0)
    subdf = points_df[points_df.order == 1]
    p = plt.semilogx(subdf['dofs'],
                     subdf['norm_time'],
                     lw=2,
                     color=cmap[0],
                     marker=markertype[0],
                     mec=cmap[0],
                     mfc=cmap[0],
                     ms=10,
                     label='Nalu p=1')

    subdf = points_df[points_df.order == 2]
    p = plt.semilogx(subdf['dofs'],
                     subdf['norm_time'],
                     lw=2,
                     color=cmap[1],
                     marker=markertype[1],
                     mec=cmap[1],
                     mfc=cmap[1],
                     ms=10,
                     label='Nalu p=2')

    # Format plots
    plt.figure(0)
    ax = plt.gca()
    plt.xlabel(r"DoF", fontsize=22, fontweight='bold')
    plt.ylabel(r"$t$ per iter. per DoF", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.tight_layout()
    plt.savefig('time.pdf', format='pdf')
    plt.savefig('time.png', format='png')
