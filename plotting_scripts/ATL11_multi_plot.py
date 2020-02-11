#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 10:08:35 2019

@author: ben
"""

import glob
import ATL11
from PointDatabase import ATL06_data
import matplotlib.pyplot as plt
import numpy as np
import pointCollection as pc
import re
pair=2
import sys
import os


def ATL11_multi_plot(ATL11_file, ATL06_wc=None, pair=2, cycles=[3, 4], hemisphere=1, xy0=None, W=5.e4):
    cyc_ind=np.array(cycles).astype(int)-1
    
    if ATL06_wc is None:
        if hemisphere==1:
            ATL06_wc="/Volumes/ice2/ben/scf/GL_06/002/cycle_*/"
            EPSG=3413
        else:
            ATL06_wc="/Volumes/ice2/ben/scf/AA_06/002/cycle_*/"
            EPSG=3031
    ATL11_re=re.compile('ATL11_(\d\d\d\d)(\d\d)')
    rgt, subprod = ATL11_re.search(ATL11_file).groups()
    
    D6_files=[];
    for cycle in range(cycles[0], cycles[1]+1):
        cyc_string='%02d' % cycle
        D6_files += glob.glob(ATL06_wc+'/ATL06_*_'+rgt+cyc_string+subprod+'_*.h5')
    print(D6_files)
    D6=[ ATL06_data(beam_pair=pair, field_dict=ATL11.misc.default_ATL06_fields()).from_file(file) for file in D6_files ]

    D11 = ATL11.data().from_file(ATL11_file, pair=pair)
    if xy0 is not None:
           bounds=[xy0[0]+np.array([-W, W]), xy0[1]+np.array([-W, W])]
           els=(D11.x >= bounds[0][0]) & (D11.x <= bounds[0][1]) & \
                (D11.y >= bounds[1][0]) & (D11.y <= bounds[1][1])
           D11=D11.index(els)
           for ii, D6i in enumerate(D6):
                els=(D6i.x >= bounds[0][0]) & (D6i.x <= bounds[0][1]) & \
                      (D6i.y >= bounds[1][0]) & (D6i.y <= bounds[1][1])
                D6[ii]=D6i.index(els)

    plt.clf()
    cyc_ind=[D11.cycles.index(cycles[0]), D11.cycles.index(cycles[1])]
    dh = D11.corrected_h.h_corr[:,cyc_ind[1]]-D11.corrected_h.h_corr[:,cyc_ind[0]]
    dh_sigma = np.sqrt(D11.corrected_h.h_corr_sigma[:,cyc_ind[1]]**2+D11.corrected_h.h_corr_sigma[:,cyc_ind[0]]**2)
    
    h0=plt.subplot(241)
    colors=['r','b']
    t0=np.zeros(len(D6))
    for ii in range(len(D6)):
        #plt.plot(D6[ii].segment_id, D6[ii].h_li,'k.')
        #good=(D6[ii].snr_significance < 0.005) & (D6[ii].min_along_track_dh < 10)
        #D6[ii].h_li[good==0]=np.NaN
        #plt.plot(D6[ii].segment_id, D6[ii].h_li, colors[ii], marker='.')
        good=D6[ii].atl06_quality_summary==0
        D6[ii].h_li[good==0]=np.NaN
        t0[ii]=np.nanmean(D6[ii].delta_time.ravel())
    for ii in np.argsort(t0):
        plt.plot(D6[ii].segment_id, D6[ii].h_li, colors[ii], marker='.', label=os.path.basename(D6_files[ii]))
    plt.ylabel('ATL06 h_li')
    #plt.legend()
    
    plt.subplot(242, sharex=h0)
    good=(D11.ref_surf.misfit_chi2r < 8000)  & (np.sum(D11.cycle_stats.seg_count[:, cyc_ind[0]:cyc_ind[1]+1]>=2, axis=1)>=1)
    plt.errorbar(D11.corrected_h.ref_pt[good], dh[good], yerr=dh_sigma[good], fmt='ko')
    
    #plt.scatter(xo.)
    #plt.plot(D11.corrected_h.ref_pt[good], dhm[good],'.')
    
    plt.ylabel('ATL11 dh')
    
    plt.subplot(243, sharex=h0)
    for ii in range(len(D6)):
        plt.plot(D6[ii].segment_id, D6[ii].y_atc,colors[ii])
    plt.ylabel('y_atc')
    
    
    plt.subplot(245, sharex=h0)
    plt.plot(D11.corrected_h.ref_pt, D11.ref_surf.deg_x,'r.')
    plt.plot(D11.corrected_h.ref_pt, D11.ref_surf.deg_y,'b.')
    plt.plot(D11.corrected_h.ref_pt, D11.ref_surf.complex_surface_flag,'k')
    plt.ylabel('degree')
    
    plt.subplot(246, sharex=h0)
    plt.plot(D11.corrected_h.ref_pt[good], D11.ref_surf.misfit_chi2r[good],'b.')
    plt.plot(D11.corrected_h.ref_pt[good], D11.ref_surf.quality_summary[good],'r.')
    plt.ylabel('misfit, quality summary')
    
    plt.subplot(247, sharex=h0)
    plt.plot(D11.corrected_h.ref_pt[good], D11.cycle_stats.seg_count[good, cyc_ind[0]],'r.')
    plt.plot(D11.corrected_h.ref_pt[good], D11.cycle_stats.seg_count[good, cyc_ind[1]],'b.')
    plt.ylabel('n_segs')


    if hemisphere==1:
        MOS=pc.grid.data().from_geotif('/Volumes/ice1/ben/MOG/2005/mog_2005_1km.tif')
        plt.subplot(144)
        MOS.show(cmap='gray', vmin=14000, vmax=17000)
        D11.get_xy(EPSG=EPSG)
        plt.plot(D11.x, D11.y,'.')


    return D11, D6
if __name__=='__main__':
    ATL11_multi_plot(sys.argv[1])
