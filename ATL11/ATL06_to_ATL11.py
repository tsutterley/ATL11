
import os
os.environ['MKL_NUM_THREADS']="1"
os.environ['NUMEXPR_NUM_THREADS']="1"
os.environ['OMP_NUM_THREADS']="1"
os.environ['OPENBLAS_NUM_THREADS']="1"

import numpy as np
import ATL11
import glob
import sys
import matplotlib.pyplot as plt
from PointDatabase import geo_index

#591 10 -F /Volumes/ice2/ben/scf/AA_06/001/cycle_02/ATL06_20190205041106_05910210_001_01.h5 -b -101. -76. -90. -74.5 -o test.h5 -G "/Volumes/ice2/ben/scf/AA_06/001/cycle*/index/GeoIndex.h5"
#591 10 -F /Volumes/ice2/ben/scf/AA_06/001/cycle_02/ATL06_20190205041106_05910210_001_01.h5 -o test.h5 -G "/Volumes/ice2/ben/scf/AA_06/001/cycle*/index/GeoIndex.h5"

def get_proj4(hemisphere):
    if hemisphere==-1:
        return'+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
    if hemisphere==1:
        return '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +k=1 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs '

def main(argv):
    # account for a bug in argparse that misinterprets negative agruents
    for i, arg in enumerate(argv):
        if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

    # command-line interface: run ATL06_to_ATL11 on a list of ATL06 files
    import argparse
    parser=argparse.ArgumentParser(description='generate an ATL11 file from a collection of ATL06 files.')
    parser.add_argument('rgt', type=int)
    parser.add_argument('subproduct', type=int)
    parser.add_argument('--directory','-d', default=os.getcwd())
    parser.add_argument('--verbose','-v', action='store_true')
    parser.add_argument('--pair','-p', type=int, default=None)
    parser.add_argument('--File', '-F', type=str, default=None)
    parser.add_argument('--GI_file_glob','-G', type=str, default=None)
    parser.add_argument('--out_file','-o', default=None, required=True)
    parser.add_argument('--first_point','-f', type=int, default=None)
    parser.add_argument('--last_point','-l', type=int, default=None)
    parser.add_argument('--num_points','-N', type=int, default=None)
    parser.add_argument('--cycles', '-c', type=int, default=2)
    parser.add_argument('--min_cycle','-m', type=int, default=0)
    parser.add_argument('--Hemisphere','-H', type=int, default=-1)
    parser.add_argument('--bounds', '-b', type=float, nargs=4, default=None, help="latlon bounds: west, south, east, north")

    args=parser.parse_args()

    if args.verbose:
        print("working on :")
        print(args.file)
    if args.File is None:
        glob_str=os.path.join(args.directory,'*ATL06*_*_{0:04d}??{1:02d}_*.h5'.format(args.rgt, args.subproduct))
        files=glob.glob(glob_str)
    else:
        files=[args.File]
    print("found ATL06 files:" + str(files))

    if args.pair is None:
        pairs=[1, 2, 3]
    else:
        pairs=[args.pair]

    if args.GI_file_glob is not None:
        GI_files=glob.glob(args.GI_file_glob)
    else:
        GI_files=None
    print("found GI files:"+str(GI_files))
    for pair in pairs:
        D6 = ATL11.read_ATL06_data(files, beam_pair=pair, min_cycle=args.min_cycle)
        D6, ref_pt_numbers, ref_pt_x = ATL11.select_ATL06_data(D6, first_ref_pt=args.first_point, last_ref_pt=args.last_point, lonlat_bounds=args.bounds, num_ref_pts=args.num_points)
        if D6 is None or len(ref_pt_numbers)==0:
            continue
        D11=ATL11.data().from_ATL06(D6, ref_pt_numbers=ref_pt_numbers, ref_pt_x=ref_pt_x,\
                      N_cycles=args.cycles, beam_pair=pair, verbose=args.verbose, \
                      GI_files=GI_files, hemisphere=args.Hemisphere) # defined in ATL06_to_ATL11
        if D11 is not None:
            D11.write_to_file(args.out_file)

    # create a geo index for the current file.  This gets saved in the '/index' group
    GI=geo_index(SRS_proj4=get_proj4(args.Hemisphere), delta=[1.e4, 1.e4]).for_file(args.out_file, 'ATL11', dir_root=os.path.dirname(args.out_file))
    GI.attrs['bin_root']=None
    # the 'file' attributes of the geo_index are of the form :pair1, :pair2, :pair3, which means that the
    # data for each bin are to be read from the current file
    for file in ['file_1','file_2','file_3']:
        if file in GI.attrs:
            temp = ':'+GI.attrs[file].split(':')[1]
            GI.attrs[file] = temp
    GI.to_file(args.out_file)

    print("ATL06_to_ATL11: done with "+args.out_file)

if __name__=="__main__":
    main(sys.argv)
