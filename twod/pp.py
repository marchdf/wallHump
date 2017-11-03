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
import re
import numpy as np
import pandas as pd
import yaml
from netCDF4 import Dataset, chartostring


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
    df['y'] = y[wall_nids_all]
    print(x[wall_nids_all])
    print(y[wall_nids_all])

    # Calculate coefficients
    u0, rho0, mu = parse_ic(yname)
    dynPres = rho0 * 0.5 * u0 * u0
    df['cf'] = df['tau_wall'] / dynPres
    df['cp'] = df['pressure'] / dynPres

    return df


# ========================================================================
def tauwall_hack(fname, yname):
    """
    Hack to find all the wall quantities by looking at tau_wall > 0.

    Identical to parse_surface_probe but without relying on sidesets.
    """

    tauwall_field = 'tau_wall'
    pressure_field = 'pressure'

    # Read in the Exodus II mesh
    msh = Dataset(fname, 'r')
    field_names = get_name_list(msh, "name_nod_var")
    tau_idx = field_names.index(tauwall_field)
    pressure_idx = field_names.index(pressure_field)

    # Get the coordinates and time
    x = msh.variables['coordx'][:]
    y = msh.variables['coordy'][:]
    time = msh.variables['time_whole'][1:]

    # Get tau_wall and pressure everywhere
    tau_wall_all = msh.variables['vals_nod_var%d' %
                                 (tau_idx + 1)][:]
    pressure_all = msh.variables['vals_nod_var%d' %
                                 (pressure_idx + 1)][:]

    # The wall is wherever tauwall is non-zero,
    wall_idx = tau_wall_all[-1, :] > 1e-16

    # Get the variables on the wall
    tau_wall = tau_wall_all[1:, wall_idx]
    pressure = pressure_all[1:, wall_idx]
    x = x[wall_idx]
    y = y[tau_wall_all[-1, :] > 1e-16]

    # Keep only the last time step in a dataframe
    df = pd.DataFrame()
    df['tau_wall'] = tau_wall[-1, :]
    df['pressure'] = pressure[-1, :]
    df['time'] = time[-1]
    df['x'] = x
    df['y'] = y

    # Calculate coefficients
    u0, rho0, mu = parse_ic(yname)
    dynPres = rho0 * 0.5 * u0 * u0
    df['cf'] = df['tau_wall'] / dynPres
    df['cp'] = df['pressure'] / dynPres

    df.sort_values(by=['x'], inplace=True, )
    return df.reset_index(drop=True)


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
        description='A post-processing tool')
    parser.add_argument(
        '-s', '--show', help='Show the plots', action='store_true')
    args = parser.parse_args()

    # ========================================================================
    # Setup

    fdirs = ['103x28']
    result_root = 'results_p'
    yname_root = 'wallHump_p'
    chord = 420
    separation_exp = 0.665
    reattachement_exp = 1.1
    pname = 'profiles.dat'
    wname = 'wall_vars.dat'
    sname = 'points.dat'

    lst_profiles = []
    lst_wall = []
    points_df = pd.DataFrame()
    for k, fdir in enumerate(fdirs):
        rdirs = sorted(glob.glob(os.path.join(fdir, result_root + '*')))
        orders = [int(re.findall(r'_p(\d+)', d)[0]) for d in rdirs]
        resolution = [int(i) for i in re.findall(r'\d+', fdir)]
        print(rdirs, orders, resolution)

        for order, rdir in zip(orders, rdirs):
            yname = os.path.join(fdir, yname_root +
                                 '{0:d}.i'.format(order))
            u0, rho0, mu = parse_ic(yname)

            # Velocities
            profiles = parse_velocity_probe(rdir, u0)
            profiles['order'] = order
            profiles['nx'] = resolution[0]
            profiles['ny'] = resolution[1]

            # Wall quantities
            fname = os.path.join(rdir, 'wallHump.e')
            walldf = tauwall_hack(fname, yname)
            walldf['order'] = order
            walldf['nx'] = resolution[0]
            walldf['ny'] = resolution[1]

            # Separation and reattachement points
            walldf['abs_cp_grad'] = np.fabs(np.gradient(walldf['cf'],
                                                        walldf['x']))

            # Search for separation point near the experimental one
            interval = 0.1
            ubnd = (1 + interval) * separation_exp
            lbnd = (1 - interval) * separation_exp
            subdf = walldf[(walldf.x >= lbnd) & (walldf.x <= ubnd)]
            idx = subdf['abs_cp_grad'].idxmax()
            separation = subdf.loc[idx].x

            # Search for reattachement point near the experimental one
            interval = 0.2
            ubnd = (1 + interval) * reattachement_exp
            lbnd = (1 - interval) * reattachement_exp
            subdf = walldf[(walldf.x >= lbnd) & (walldf.x <= ubnd)]
            idx = subdf['abs_cp_grad'].idxmax()
            reattachement = subdf.loc[idx].x

            # Append to previous dataframes
            lst_profiles.append(profiles)
            lst_wall.append(walldf)
            points_df = points_df.append({'res': fdir,
                                          'nx': resolution[0],
                                          'ny': resolution[1],
                                          'separation': separation,
                                          'reattachement': reattachement},
                                         ignore_index=True)

    # Concatenate all together and save
    profiles = pd.concat(lst_profiles, ignore_index=True)
    walldf = pd.concat(lst_wall, ignore_index=True)
    profiles.to_csv(pname, index=False)
    walldf.to_csv(wname, index=False)
    points_df.to_csv(sname, index=False)
