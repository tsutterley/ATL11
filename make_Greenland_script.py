# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 16:24:11 2017

@author: ben
"""
import os
import re
import numpy as np

# script that generates a list of command line calls that can be passed to gnu parallel

# np.seterr(all=raise)
np.seterr(invalid='ignore')

ATL06_base = '/Volumes/ice2/nick/34_IceSat2/SyntheticTracks_Alex/ATL06/'
# Make a list of the track files in the first cycle:
track_file_list=[X for X in os.listdir(os.path.join(ATL06_base,'TrackData_01') if re.match(r'(.*?).h5')]
# loop over the tracks:
for track_file in track_file_list:
    # loop over pairs
    for pair in [1,2,3]:
        # establish output file name
        m = re.search(r"Track_(.*?).h5",track_file)
        track_num=int(m.group(1))
        fileout = 'ATL11_Track{0}_Pair{1:d}.h5'.format(track_num, pair)
        if os.path.isfile(fileout):
            continue
        glob_str=os.path.join(ATL06_base,'TrackData_v4*','Track_{0}.h5'.format(m.group(1)))
        print("python ATL06_to_ATL11.py --ATL06_glob '{0}' -o {1} -v -p {2:d} -t {3:d}".format(glob_str, fileout, pair, track_num))
