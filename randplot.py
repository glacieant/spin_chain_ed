#!/usr/bin/env python
import sys
import itertools

import numpy as np

import matplotlib.pyplot as plt

## plot routines

# palette
plt.style.use('tableau-colorblind10')

# marker strings
markers = itertools.cycle(['s','o','^','D','X','p','H','P'])

## mathtext style
plt.rcParams['mathtext.fontset']='cm'
plt.rc('font',family='serif')

DISORDER = sys.argv[1]
YLABEL = sys.argv[2]
NAME = r'$\Delta_h$' if DISORDER == "field" else r'$\Delta$'

sitenum, logavg, relerr, delta = np.loadtxt('randoutput.txt', usecols=(0,1,2,4), unpack=True)

afig, ax = plt.subplots()

for dlt in np.unique(delta):

    sel = delta == dlt

    line = ax.errorbar(sitenum[sel],
                       logavg[sel],
                       yerr=relerr[sel],
                       ms=10,
                       marker=next(markers),
                       mfc='None',
                       lw=4,
                       linestyle='--',
                       capsize=5,
                       elinewidth=2,
                       label=NAME+r'$=$'+str("%g" % dlt),
                       )

ax.tick_params(which='major',width=2,
               labelsize=22,direction='in',
               length=10,
               bottom=True,top=False,
               left=True,right=False)

legend = ax.legend(
        loc='lower left',
        fontsize=18,
        facecolor='w',
        edgecolor='k',
        framealpha=1,
        ncol=1,
        borderpad=0.5,
        handletextpad=0.5,
        handlelength=1.0,
        )

ax.set_xlabel(r'$L$',fontsize=22)
ax.set_ylabel(YLABEL,fontsize=18)
ax.set_xlim(left=np.amin(sitenum)-1,right=np.amax(sitenum)+1)
ax.set_box_aspect(1)

afig.savefig("heisenberg_plot_power_distro.pdf",
             bbox_inches='tight',transparent=True,
             )

plt.close('all')
