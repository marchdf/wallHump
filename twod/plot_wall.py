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
import seaborn as sns
import yaml
from netCDF4 import Dataset, chartostring

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


# ========================================================================
def parse_velocity_probe(rdir, uref):
    """Read a Nalu velocity probe file."""

    fnames = glob.glob(os.path.join(rdir, 'probe_profile*.dat'))

    lst = []
    for fname in fnames:

        df = pd.read_csv(fname, delim_whitespace=True)
        df.rename(columns={'Time': 'time',
                           'coordinates[0]': 'x',
                           'coordinates[1]': 'y',
                           'velocity_probe[0]': 'u',
                           'velocity_probe[1]': 'v',
                           'reynolds_stress_probe[0]': 'uu',
                           'reynolds_stress_probe[1]': 'uv',
                           'reynolds_stress_probe[2]': 'vv'},
                  inplace=True)
        lst.append(df)
    df = pd.concat(lst)

    # Keep only the last time step
    time = np.unique(df['time'])
    df = df.loc[df['time'] == time[-1]]

    # Non-dimensionalize some quantities
    df['u'] = df['u'] / uref
    df['v'] = df['v'] / uref
    df['uu'] = df['uu'] / uref**2
    df['uv'] = df['uv'] / uref**2
    df['vv'] = df['vv'] / uref**2
    return df.reset_index(drop=True)


# ========================================================================
def parse_flat_wall_probe(fname, yname):
    """Read a Nalu wall probe file."""

    df = pd.read_csv(fname, delim_whitespace=True)
    df.rename(columns={'Time': 'time',
                       'coordinates[0]': 'x',
                       'coordinates[1]': 'y',
                       'tau_wall_probe[0]': 'tau_wall',
                       'pressure_probe[0]': 'pressure'},
              inplace=True)

    # Keep only the last time step
    time = np.unique(df['time'])
    df = df.loc[df['time'] == time[-1]]

    # Calculate coefficients
    u0, rho0, mu = parse_ic(yname)
    dynPres = rho0 * 0.5 * u0 * u0
    df['cf'] = df['tau_wall'] / dynPres
    df['cp'] = df['pressure'] / dynPres

    return df.reset_index(drop=True)


# ========================================================================
def parse_surface_probe(fname, yname):
    """Read a Nalu exodus file to get surface data."""

    wall_bdy = 'bottomwall'
    tauwall_field = 'tau_wall'
    pressure_field = 'pressure'

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
    time = msh.variables['time_whole'][1:]

    # Element mapping and wall node ids
    nids = msh.variables['connect1'][:]
    wall_elems = msh.variables['elem_ss%d' % (wall_idx + 1)][:] - 1
    wall_nids_all = np.unique(
        nids[np.ix_(wall_elems, ss_node_ids)].flatten()) - 1

    # Get tau_wall and pressure on the wall
    tau_wall_all = msh.variables['vals_nod_var%d' %
                                 (tau_idx + 1)][:][1:, wall_nids_all]
    pressure_all = msh.variables['vals_nod_var%d' %
                                 (pressure_idx + 1)][:][1:, wall_nids_all]

    # Keep only the last time step in a dataframe
    df = pd.DataFrame()
    df['tau_wall'] = tau_wall_all[-1, :]
    df['pressure'] = pressure_all[-1, :]
    df['time'] = time[-1]
    df['x'] = x[wall_nids_all]

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
    porder = 2
    fdir = os.path.abspath('103x28')
    # fdir = os.path.abspath('205x55')
    # fdir = os.path.abspath('409x109')
    # fdir = os.path.abspath('817x217')
    # fdir = os.path.abspath('1633x433')
    rdir = os.path.abspath(os.path.join(fdir, 'results_p{0:d}'.format(porder)))
    yname = os.path.join(fdir, 'wallHump_p{0:d}.i'.format(porder))
    u0, rho0, mu = parse_ic(yname)
    chord = 420
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

    # # ======================================================================
    # # Nalu output
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

    plt.figure(1)
    p = plt.plot(subdf['x'],
                 subdf['cf'],
                 color=cmap[1],
                 lw=2,
                 label='Nalu p=1')

    subdf = df[(df.nx == 103) & (df.order == 2)]
    subdf = subdf.reset_index()
    subdf = subdf.iloc[::2, :]

    plt.figure(0)
    p = plt.plot(subdf['x'],
                 -subdf['cp'] - nasa_cp_shift,
                 color=cmap[2],
                 lw=2,
                 label='Nalu p=2')

    plt.figure(1)
    p = plt.plot(subdf['x'],
                 subdf['cf'],
                 color=cmap[2],
                 lw=2,
                 label='Nalu p=2')

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
