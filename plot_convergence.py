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
# Function definitions
#
# ========================================================================
def parse_ic(fname):
    """Parse the Nalu yaml input file for the initial conditions"""
    with open(fname, 'r') as stream:
        try:
            dat = yaml.load(stream)
            u0 = float(dat['realms'][0]['initial_conditions']
                       [0]['value']['velocity'][0])
            rho0 = float(dat['realms'][0]['material_properties']
                         ['specifications'][0]['value'])
            mu = float(dat['realms'][0]['material_properties']
                       ['specifications'][1]['value'])

            return u0, rho0, mu

        except yaml.YAMLError as exc:
            print(exc)


def parse_surface_probe(fname, yname):
    """Read a Nalu exodus file to get surface data."""

    wall_bdy = 'bottomwall'
    tauwall_field = 'tau_wall'
    pressure_field = 'pressure'
    epsilon = -1e-2

    # Indices of wall face in element
    ss_node_ids = np.array([0, 1, 2, 3])

    # Read in the Exodus II mesh
    msh = Dataset(fname, 'r')
    ss_names = get_name_list(msh, "ss_names")
    field_names = get_name_list(msh, "name_nod_var")
    wall_idx = ss_names.index(wall_bdy)
    tau_idx = field_names.index(tauwall_field)
    pressure_idx = field_names.index(pressure_field)

    # Get the coordinates and time
    x = msh.variables['coordx'][:]
    y = msh.variables['coordy'][:]
    z = msh.variables['coordz'][:]
    time = msh.variables['time_whole'][1:]

    # Element mapping and wall node ids
    nids = msh.variables['connect1'][:]
    wall_elems = msh.variables['elem_ss%d' % (wall_idx + 1)][:] - 1
    wall_nids_all = np.unique(
        nids[np.ix_(wall_elems, ss_node_ids)].flatten()) - 1
    wall_nids = wall_nids_all[y[wall_nids_all] > epsilon]

    # Get tau_wall and pressure on the wall
    tau_wall_all = msh.variables['vals_nod_var%d' %
                                 (tau_idx + 1)][:][1:, wall_nids]
    pressure_all = msh.variables['vals_nod_var%d' %
                                 (pressure_idx + 1)][:][1:, wall_nids]

    # Keep only the last time step in a dataframe
    df = pd.DataFrame()
    df['tau_wall'] = tau_wall_all[-1, :]
    df['pressure'] = pressure_all[-1, :]
    df['time'] = time[-1]
    df['x'] = x[wall_nids]

    # Calculate coefficients
    u0, rho0, mu = parse_ic(yname)
    dynPres = rho0 * 0.5 * u0 * u0
    df['cf'] = df['tau_wall'] / dynPres
    df['cp'] = df['pressure'] / dynPres

    return df


# ========================================================================
def get_name_list(msh, varname):
    """Return a list of python strings converted from Exodus char format

    Args:
        msh (Dataset): A netCDF4 dataset instance
        varname (str): A string instance
    """
    return [str(chartostring(v)) for v in msh.variables[varname]]


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
    fdirs = ['103x28', '205x55', '409x109', '817x217', '1633x433']
    chord = 420
    nasa_cp_shift = -0.015
    separation_exp = 0.665
    reattachement_exp = 1.1

    points_df = pd.DataFrame()
    for k, fdir in enumerate(fdirs):

        resolution = [int(i) for i in re.findall(r'\d+', fdir)]
        rdir = os.path.abspath(os.path.join(os.path.abspath(fdir), 'results'))
        yname = os.path.join(os.path.abspath(fdir), 'wallHump.i')
        u0, rho0, mu = parse_ic(yname)

        # Cp and Cf
        fname = os.path.join(rdir, 'wallHump.e')
        df = parse_surface_probe(fname, yname)
        df['abs_cp_grad'] = np.fabs(np.gradient(df['cf'], df['x']))

        # Search for separation point near the experimental one
        interval = 0.1
        ubnd = (1 + interval) * separation_exp
        lbnd = (1 - interval) * separation_exp
        subdf = df[(df.x >= lbnd) & (df.x <= ubnd)]
        idx = subdf['abs_cp_grad'].idxmax()
        separation = subdf.loc[idx].x

        # Search for reattachement point near the experimental one
        interval = 0.2
        ubnd = (1 + interval) * reattachement_exp
        lbnd = (1 - interval) * reattachement_exp
        subdf = df[(df.x >= lbnd) & (df.x <= ubnd)]
        idx = subdf['abs_cp_grad'].idxmax()
        reattachement = subdf.loc[idx].x

        points_df = points_df.append({'res': fdir,
                                      'nx': resolution[0],
                                      'ny': resolution[1],
                                      'separation': separation,
                                      'reattachement': reattachement},
                                     ignore_index=True)

    # Plot Nalu data
    plt.figure(0)
    p = plt.plot(points_df['nx'],
                 points_df['separation'],
                 lw=2,
                 color=cmap[0],
                 marker=markertype[0],
                 mec=cmap[0],
                 mfc=cmap[0],
                 ms=10,
                 label='Nalu')

    plt.figure(1)
    p = plt.plot(points_df['nx'],
                 points_df['reattachement'],
                 lw=2,
                 color=cmap[0],
                 marker=markertype[0],
                 mec=cmap[0],
                 mfc=cmap[0],
                 ms=10,
                 label='Nalu')

    # Plot experimental data
    nmin = np.min(points_df.nx)
    nmax = np.max(points_df.nx)
    plt.figure(0)
    p = plt.plot([nmin, nmax],
                 [separation_exp, separation_exp],
                 color=cmap[-1],
                 label='Exp')
    plt.figure(1)
    p = plt.plot([nmin, nmax],
                 [reattachement_exp, reattachement_exp],
                 color=cmap[-1],
                 label='Exp')

    # Format plots
    plt.figure(0)
    ax = plt.gca()
    plt.xlabel(r"$N$", fontsize=22, fontweight='bold')
    plt.ylabel(r"Separation", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.tight_layout()
    # plt.savefig('separation.pdf', format='pdf')
    plt.savefig('separation.png', format='png')

    plt.figure(1)
    ax = plt.gca()
    plt.xlabel(r"$N$", fontsize=22, fontweight='bold')
    plt.ylabel(r"Reattachement", fontsize=22, fontweight='bold')
    plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    legend = ax.legend(loc='best')
    plt.tight_layout()
    # plt.savefig('reattachement.pdf', format='pdf')
    plt.savefig('reattachement.png', format='png')
