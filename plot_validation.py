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
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
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
                           'coordinates[2]': 'z',
                           'velocity_probe[0]': 'u',
                           'velocity_probe[1]': 'v',
                           'velocity_probe[2]': 'w',
                           'reynolds_stress_probe[0]': 'uu',
                           'reynolds_stress_probe[1]': 'uv',
                           'reynolds_stress_probe[2]': 'uw',
                           'reynolds_stress_probe[3]': 'vv',
                           'reynolds_stress_probe[4]': 'vw',
                           'reynolds_stress_probe[5]': 'ww'},
                  inplace=True)
        lst.append(df)
    df = pd.concat(lst)

    # Keep only the last time step
    time = np.unique(df['time'])
    df = df.loc[df['time'] == time[-1]]

    # Non-dimensionalize some quantities
    df['u'] = df['u'] / uref
    df['v'] = df['v'] / uref
    df['w'] = df['w'] / uref
    df['uu'] = df['uu'] / uref**2
    df['uv'] = df['uv'] / uref**2
    df['uw'] = df['uw'] / uref**2
    df['vv'] = df['vv'] / uref**2
    df['vw'] = df['vw'] / uref**2
    df['ww'] = df['ww'] / uref**2
    return df.reset_index(drop=True)


# ========================================================================
def parse_wall_probe(fname, yname):
    """Read a Nalu wall probe file."""

    df = pd.read_csv(fname, delim_whitespace=True)
    df.rename(columns={'Time': 'time',
                       'coordinates[0]': 'x',
                       'coordinates[1]': 'y',
                       'coordinates[2]': 'z',
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
    fdir = os.path.abspath('103x28')
    rdir = os.path.abspath(os.path.join(fdir, 'results'))
    yname = os.path.join(fdir, 'wallHump.i')
    u0, rho0, mu = parse_ic(yname)
    chord = 420

    # ======================================================================
    # Profile data

    # NASA experimental data
    fname = os.path.join(os.path.abspath('nasa_data'), 'profiles_exp.dat')
    nasadf = pd.read_csv(fname, comment='#')
    fname = os.path.join(os.path.abspath('nasa_data'), 'u_inflow_exp.dat')
    inflow = pd.read_csv(fname, comment='#')
    inflow['u'] /= u0
    inflow['y'] /= chord
    inflow['x'] = -2.14
    nasadf = pd.concat([nasadf, inflow])
    nasadf.rename(columns={'y': 'z', 'z': 'y'}, inplace=True)
    nasadf['label'] = 'Exp'

    # NASA CFL3D data
    fname = os.path.join(os.path.abspath('nasa_data'), 'profiles_cfl3d.dat')
    cfl3ddf = pd.read_csv(fname, comment='#')
    cfl3ddf['label'] = 'CFL3D'
    cfl3ddf.rename(columns={'y': 'z', 'z': 'y'}, inplace=True)
    xslices = [-2.14, 0.65, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
    for xslice in xslices:
        cfl3ddf.loc[np.fabs(cfl3ddf['x'] - xslice) < 1e-3, 'x'] = xslice

    # Nalu data
    naludf = parse_velocity_probe(rdir, u0)
    naludf['label'] = 'Nalu'

    # Merge them all
    profiles = pd.concat([nasadf, cfl3ddf, naludf])

    # Explore
    with sns.axes_style("white"):
        sns.set_context("talk")

        p = sns.FacetGrid(col='x',
                          col_wrap=4,
                          hue='label',
                          hue_kws=dict(marker=["s", "", ""],
                                       lw=[0, 2, 2]),
                          sharex=False,
                          sharey=False,
                          data=profiles)
        p = p.map(plt.scatter, 'u', 'z')
        p = p.map(plt.plot, 'u', 'z')

    # ======================================================================
    # Surface quantities

    # NASA experimental data
    fname = os.path.join(os.path.abspath('nasa_data'), 'cp_exp.dat')
    nasacp = pd.read_csv(fname, comment='#')
    fname = os.path.join(os.path.abspath('nasa_data'), 'cf_exp.dat')
    nasacf = pd.read_csv(fname, comment='#')
    nasasurf = pd.concat([nasacp, nasacf])
    nasasurf['label'] = 'Exp'

    # NASA CFL3D data
    fname = os.path.join(os.path.abspath('nasa_data'), 'cp_cfl3d.dat')
    cfl3dcp = pd.read_csv(fname, comment='#')
    fname = os.path.join(os.path.abspath('nasa_data'), 'cf_cfl3d.dat')
    cfl3dcf = pd.read_csv(fname, comment='#')
    cfl3dsurf = pd.concat([cfl3dcp, cfl3dcf])
    cfl3dsurf['label'] = 'CFL3D'

    fname = os.path.join(rdir, 'probe_bottomwall_0.dat')
    nalusurf = parse_wall_probe(fname, yname)
    nalusurf['label'] = 'Nalu'

    # Merge them all
    surface = pd.concat([nasasurf, cfl3dsurf, nalusurf])
    surface['cp'] *= -1

    # Explore
    with sns.axes_style("white"):
        sns.set_context("talk")

        p = sns.FacetGrid(hue='label',
                          hue_kws=dict(marker=["s", "", ""],
                                       lw=[0, 2, 2]),
                          sharex=False,
                          sharey=False,
                          data=surface)
        p = p.map(plt.scatter, 'x', 'cp')
        p = p.map(plt.plot, 'x', 'cp')

        p = sns.FacetGrid(hue='label',
                          hue_kws=dict(marker=["s", "", ""],
                                       lw=[0, 2, 2]),
                          sharex=False,
                          sharey=False,
                          data=surface)
        p = p.map(plt.scatter, 'x', 'cf')
        p = p.map(plt.plot, 'x', 'cf')

    # ======================================================================
    # Save the plots
    plt.figure(1)
    plt.savefig('profiles.png', format='png')

    plt.figure(2)
    plt.savefig('cp.png', format='png')

    plt.figure(3)
    plt.savefig('cf.png', format='png')

    if args.show:
        plt.show()

    #
    # KEEPING IN CASE THIS IS NECESSARY LATER
    #

    # # ======================================================================
    # # NASA experiment
    # fname = os.path.join(os.path.abspath('nasa_data'), 'profiles_exp.dat')
    # df = pd.read_csv(fname, comment='#')

    # x = 0.65
    # subdf = df[df['x'] == x]
    # plt.figure(0)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # plt.figure(1)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # x = 0.8
    # subdf = df[df['x'] == x]
    # plt.figure(2)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # plt.figure(3)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # x = 0.9
    # subdf = df[df['x'] == x]
    # plt.figure(4)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # plt.figure(5)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # x = 1.0
    # subdf = df[df['x'] == x]
    # plt.figure(6)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # plt.figure(7)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # x = 1.1
    # subdf = df[df['x'] == x]
    # plt.figure(8)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # plt.figure(9)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # x = 1.2
    # subdf = df[df['x'] == x]
    # plt.figure(10)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # plt.figure(11)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # x = 1.3
    # subdf = df[df['x'] == x]
    # plt.figure(12)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # plt.figure(13)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # x = -2.14
    # fname = os.path.join(os.path.abspath('nasa_data'), 'u_inflow_exp.dat')
    # df = pd.read_csv(fname, comment='#')
    # plt.figure(14)
    # p = plt.plot(df['u'] / u0,
    #              df['y'] / chord,
    #              ls='None',
    #              lw=1,
    #              color=cmap[-1],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp. x={0:.2f}'.format(x))

    # # Cp
    # fname = os.path.join(os.path.abspath('nasa_data'), 'cp_exp.dat')
    # df = pd.read_csv(fname, comment='#')

    # plt.figure(15)
    # p = plt.plot(df['x'],
    #              -df['cp'],
    #              ls='None',
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='Exp.')

    # # Cf
    # fname = os.path.join(os.path.abspath('nasa_data'), 'cf_exp.dat')
    # df = pd.read_csv(fname, comment='#')

    # plt.figure(16)
    # p = plt.errorbar(df['x'],
    #                  df['cf'],
    #                  yerr=df['error'],
    #                  fmt='o',
    #                  capsize=3,
    #                  color=cmap[-1],
    #                  marker=markertype[0],
    #                  mec=cmap[-1],
    #                  mfc=cmap[-1],
    #                  ms=6,
    #                  label='Exp.')

    # # ======================================================================
    # # NASA CFL3D
    # fname = os.path.join(os.path.abspath('nasa_data'), 'profiles_cfl3d.dat')
    # df = pd.read_csv(fname, comment='#')

    # x = 0.65
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(0)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              marker=markertype[0],
    #              mec=cmap[-1],
    #              mfc=cmap[-1],
    #              ms=6,
    #              label='CFL3D x={0:.2f}'.format(x))

    # plt.figure(1)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # x = 0.8
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(2)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # plt.figure(3)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # x = 0.9
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(4)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # plt.figure(5)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # x = 1.0
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(6)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # plt.figure(7)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # x = 1.1
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(8)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # plt.figure(9)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # x = 1.2
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(10)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # plt.figure(11)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # x = 1.3
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(12)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # plt.figure(13)
    # p = plt.plot(subdf['uv'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # x = -2.14
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(14)
    # p = plt.plot(subdf['u'],
    #              subdf['y'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[0],
    #              label='CFL3D x={0:.2f}'.format(x))

    # # Cp
    # fname = os.path.join(os.path.abspath('nasa_data'), 'cp_cfl3d.dat')
    # df = pd.read_csv(fname, comment='#')

    # plt.figure(15)
    # p = plt.plot(df['x'],
    #              -df['cp'],
    #              color=cmap[0],
    #              label='CFL3D')

    # # Cf
    # fname = os.path.join(os.path.abspath('nasa_data'), 'cf_cfl3d.dat')
    # df = pd.read_csv(fname, comment='#')

    # plt.figure(16)
    # p = plt.plot(df['x'],
    #              df['cf'],
    #              color=cmap[0],
    #              label='CFL3D')

    # # ======================================================================
    # # Nalu output
    # df = parse_velocity_probe(rdir, u0)

    # x = 0.65
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(0)
    # p = plt.plot(subdf['u'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # plt.figure(1)
    # p = plt.plot(subdf['uv'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # x = 0.8
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(2)
    # p = plt.plot(subdf['u'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # plt.figure(3)
    # p = plt.plot(subdf['uv'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # x = 0.9
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(4)
    # p = plt.plot(subdf['u'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # plt.figure(5)
    # p = plt.plot(subdf['uv'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # x = 1.0
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(6)
    # p = plt.plot(subdf['u'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # plt.figure(7)
    # p = plt.plot(subdf['uv'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # x = 1.1
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(7)
    # p = plt.plot(subdf['u'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # plt.figure(8)
    # p = plt.plot(subdf['uv'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # x = 1.2
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(9)
    # p = plt.plot(subdf['u'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # plt.figure(10)
    # p = plt.plot(subdf['uv'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # x = 1.3
    # subdf = df[np.fabs(df['x'] - x) < 1e-3]
    # plt.figure(11)
    # p = plt.plot(subdf['u'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # plt.figure(12)
    # p = plt.plot(subdf['uv'],
    #              subdf['z'],
    #              ls='-',
    #              lw=2,
    #              color=cmap[1],
    #              label='Nalu x={0:.2f}'.format(x))
    # p[0].set_dashes(dashseq[1])

    # # Cp and Cf
    # fname = os.path.join(rdir, 'probe_bottomwall_0.dat')
    # df = parse_wall_probe(fname, yname)

    # plt.figure(15)
    # p = plt.plot(df['x'],
    #              -df['cp'],
    #              color=cmap[1],
    #              label='Nalu')

    # plt.figure(16)
    # p = plt.plot(df['x'],
    #              df['cf'],
    #              color=cmap[1],
    #              label='Nalu')

    # # ======================================================================
    # # Format the plots
    # plt.figure(0)
    # ax = plt.gca()
    # plt.xlabel(r"$u / u_r$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"$y$", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # ax.set_xlim([0, 1.2])
    # ax.set_ylim([1, 5])
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # #plt.savefig('u_m4.pdf', format='pdf')
    # plt.savefig('u_m4.png', format='png')

    # plt.figure(1)
    # ax = plt.gca()
    # plt.xlabel(r"$1000 \times u'v' / u_r$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"$y$", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # ax.set_xlim([-1.5, 0.5])
    # ax.set_ylim([1, 5])
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # #plt.savefig('uv_m4.pdf', format='pdf')
    # plt.savefig('uv_m4.png', format='png')

    # plt.figure(2)
    # ax = plt.gca()
    # plt.xlabel(r"$u / u_r$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"$y$", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # ax.set_xlim([-0.4, 1.2])
    # ax.set_ylim([0, 3])
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # #plt.savefig('u_1.pdf', format='pdf')
    # plt.savefig('u_1.png', format='png')

    # plt.figure(3)
    # ax = plt.gca()
    # plt.xlabel(r"$1000 \times u'v' / u_r$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"$y$", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # ax.set_xlim([-15, 5])
    # ax.set_ylim([0, 3])
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # #plt.savefig('uv_1.pdf', format='pdf')
    # plt.savefig('uv_1.png', format='png')

    # plt.figure(4)
    # ax = plt.gca()
    # plt.xlabel(r"$u / u_r$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"$y$", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # ax.set_xlim([-0.4, 1.2])
    # ax.set_ylim([0, 3])
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # #plt.savefig('u_4.pdf', format='pdf')
    # plt.savefig('u_4.png', format='png')

    # plt.figure(5)
    # ax = plt.gca()
    # plt.xlabel(r"$1000 \times u'v' / u_r$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"$y$", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # ax.set_xlim([-15, 5])
    # ax.set_ylim([0, 3])
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # #plt.savefig('uv_4.pdf', format='pdf')
    # plt.savefig('uv_4.png', format='png')

    # plt.figure(6)
    # ax = plt.gca()
    # plt.xlabel(r"$u / u_r$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"$y$", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # ax.set_xlim([-0.4, 1.2])
    # ax.set_ylim([0, 3])
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # #plt.savefig('u_6.pdf', format='pdf')
    # plt.savefig('u_6.png', format='png')

    # plt.figure(7)
    # ax = plt.gca()
    # plt.xlabel(r"$1000 \times u'v' / u_r$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"$y$", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # ax.set_xlim([-15, 5])
    # ax.set_ylim([0, 3])
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # #plt.savefig('uv_6.pdf', format='pdf')
    # plt.savefig('uv_6.png', format='png')

    # plt.figure(8)
    # ax = plt.gca()
    # plt.xlabel(r"$u / u_r$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"$y$", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # ax.set_xlim([-0.4, 1.2])
    # ax.set_ylim([0, 3])
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # #plt.savefig('u_10.pdf', format='pdf')
    # plt.savefig('u_10.png', format='png')

    # plt.figure(9)
    # ax = plt.gca()
    # plt.xlabel(r"$1000 \times u'v' / u_r$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"$y$", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # ax.set_xlim([-15, 5])
    # ax.set_ylim([0, 3])
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # #plt.savefig('uv_10.pdf', format='pdf')
    # plt.savefig('uv_10.png', format='png')

    # plt.figure(10)
    # ax = plt.gca()
    # plt.xlabel(r"$x$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"$C_p$", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # ax.set_xlim([-5, 30])
    # ax.set_ylim([-0.3, 0.1])
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # #plt.savefig('cp.pdf', format='pdf')
    # plt.savefig('cp.png', format='png')

    # plt.figure(11)
    # ax = plt.gca()
    # plt.xlabel(r"$x$", fontsize=22, fontweight='bold')
    # plt.ylabel(r"$|C_f|$", fontsize=22, fontweight='bold')
    # plt.setp(ax.get_xmajorticklabels(), fontsize=18, fontweight='bold')
    # plt.setp(ax.get_ymajorticklabels(), fontsize=18, fontweight='bold')
    # ax.set_xlim([-5, 30])
    # ax.set_ylim([-0.002, 0.004])
    # legend = ax.legend(loc='best')
    # plt.tight_layout()
    # #plt.savefig('cf.pdf', format='pdf')
    # plt.savefig('cf.png', format='png')
