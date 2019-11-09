#!/usr/bin/env python
u"""
export_ATL11_pairs.py
Written by Tyler Sutterley (11/2019)
Exports pairs of ATL11 data to an output file format

COMMAND LINE OPTIONS:
    -h, --help: Lists the command line options
    -D X, --directory=X: input directory of ATL11 files
    -O X, --output=X: output directory of combined ATL11 files
    --release=X: ICESat-2 data release to run
    --track=X: ICESat-2 tracks to run separated by commas
    --granule=X: ICESat-2 granule regions to run separated by commas
    --cycle=X: ICESat-2 cycles to run separated by commas
    -P X, --projection=X: Spatial projection (EPSG code)
    -T X, --time=X: Temporal format for output (julian or decimal)
    -F X, --format=X: output data format (ascii, netCDF4, HDF5)
    --bbox=X: Bounding box for spatially subsetting (Xmin,Xmax,Ymin,Ymax)
    --polygon=X: Polygon file (shp, geojson, kml) for subsetting data
    -M X, --mode=X: output data file permissions format

PYTHON DEPENDENCIES:
	numpy: Scientific Computing Tools For Python
		http://www.numpy.org
		http://www.scipy.org/NumPy_for_Matlab_Users
	h5py: Python interface for Hierarchal Data Format 5 (HDF5)
		http://h5py.org
		http://docs.h5py.org/en/stable/mpi.html
	netCDF4: Python interface to the netCDF C library
	 	https://unidata.github.io/netcdf4-python/netCDF4/index.html
	fiona: Python wrapper for vector data access functions from the OGR library
		https://fiona.readthedocs.io/en/latest/manual.html
	geopandas: Python tools for geographic data
		http://geopandas.readthedocs.io/
	shapely: PostGIS-ish operations outside a database context for Python
		http://toblerity.org/shapely/index.html
	pyproj: Python interface to PROJ library
		https://pypi.org/project/pyproj/

UPDATE HISTORY:
	Written 11/2019
"""

import sys
import os
import re
import getopt
import ATL11
import pyproj
import datetime
import numpy as np
import h5py, netCDF4
from shapely.geometry import MultiPoint, Polygon
from PointDatabase import mapData, point_data
from ATL11.convert_julian import convert_julian
from ATL11.count_leap_seconds import count_leap_seconds
from ATL11.convert_calendar_decimal import convert_calendar_decimal
from subsetting.read_geojson_file import read_geojson_file
from subsetting.read_kml_file import read_kml_file
from subsetting.read_shapefile import read_shapefile

#-- PURPOSE: help module to describe the optional input command-line parameters
def usage():
    print('\nHelp: {0}'.format(os.path.basename(sys.argv[0])))
    print(' -D X, --directory=X\tInput directory of ATL11 files')
    print(' -O X, --output=X\tOutput directory of combined ATL11 files')
    print(' --release=X\t\tICESat-2 data release to run')
    print(' --track=X\t\tICESat-2 tracks to run separated by commas')
    print(' --granule=X\t\tICESat-2 granule regions to run separated by commas')
    print(' --cycle=X\t\tICESat-2 cycles to run separated by commas')
    print(' -P X, --projection=X\tSpatial projection (EPSG code)')
    print(' -T X, --time=X\t\tTemporal format for output (julian or decimal)')
    print(' -F X, --format=X\tOutput data format (ascii, netCDF4, HDF5)')
    print(' --bbox=X\t\tBounding box for spatially subsetting')
    print(' --polygon=X\t\tPolygon file for spatially subsetting')
    print('\tshapefile (shp), geojson (json), or Keyhole Markup (kml)')
    print(' -V, --verbose\t\tOutput information for each created file')
    print(' -M X, --mode=X\t\tPermission mode of files created\n')

#-- Main program that calls export_ATL11_pairs()
def main():
    #-- Read the system arguments listed after the program
    long_options = ['help','directory=','output=','release=','track=',
        'granule=','cycle=','projection=','time=','format=','bbox=','polygon=',
        'verbose','mode=']
    optlist,arglist=getopt.getopt(sys.argv[1:],'hD:O:P:T:F:VM:',long_options)

    #-- command line parameters
    input_dir = os.path.join(os.sep,'Volumes','ice2','ben','scf','GL_11')
    output_dir = os.getcwd()
    #-- ICESat-2 data parameters for finding files
    regex_release = '002'
    regex_track = '\d{4}'
    regex_granule = '\d{2}'
    #-- ICESat-2 cycles to merge into combined file
    CYCLES = [3,4]
    #-- output spatial and temporal parameters
    EPSG = 3413
    TIME = 'decimal'
    #-- output file format
    FORMAT = 'ascii'
    #-- spatial subsetting parameters
    BBOX = None
    POLYGON = None
    #-- output information about each input and output file
    VERBOSE = False
    #-- permissions mode of the output files (number in octal)
    MODE = 0o775
    for opt, arg in optlist:
        if opt in ('-h','--help'):
            usage()
            sys.exit()
        elif opt in ("-D","--directory"):
            input_dir = os.path.expanduser(arg)
        elif opt in ("-O","--output"):
            output_dir = os.path.expanduser(arg)
        elif opt in ("--release"):
            RELEASE = np.int(arg)
            regex_release = '{0:03d}'.format(RELEASE)
        elif opt in ("--track"):
            TRACK = np.array(arg.split(','), dtype=np.int)
            regex_track = '|'.join(['{0:04d}'.format(T) for T in TRACK])
        elif opt in ("--granule"):
            GRANULE = np.array(arg.split(','), dtype=np.int)
            regex_granule = '|'.join(['{0:02d}'.format(G) for G in GRANULE])
        elif opt in ("--cycle"):
            CYCLES = np.array(arg.split(','), dtype=np.int)
        elif opt in ("-P","--projection"):
            EPSG = np.int(arg)
        elif opt in ("-T","--time"):
            TIME = arg.lower()
        elif opt in ("-F","--format"):
            FORMAT = arg
        elif opt in ("--bbox"):
            BBOX = np.array(arg.split(','), dtype=np.float)
        elif opt in ("--polygon"):
            POLYGON = os.path.expanduser(arg)
        elif opt in ("-V","--verbose"):
            VERBOSE = True
        elif opt in ("-M","--mode"):
            MODE = int(arg, 8)

    #-- compile regular expression operator for extracting data from files
    args = ('ATL11',regex_track,regex_granule,regex_release)
    regex = re.compile('({0})_({1})({2})_({3})_(\d{{3}}).h5$'.format(*args))
    #-- run program with parameters
    export_ATL11_pairs(input_dir, output_dir, regex, CYCLES=CYCLES, EPSG=EPSG,
        TIME=TIME, BBOX=BBOX, POLYGON=POLYGON, FORMAT=FORMAT, VERBOSE=VERBOSE,
        MODE=MODE)

#-- PURPOSE: read ATL11 pairs and output as a combined file
#-- can spatially subset using a bounding box or a polygon file
def export_ATL11_pairs(input_dir, output_dir, regex, CYCLES=None, EPSG=None,
    TIME=None, BBOX=None, POLYGON=None, FORMAT=None, VERBOSE=False, MODE=0o775):
    #-- recursively create output directory and set permissions mode
    if not os.access(output_dir, os.F_OK):
        os.makedirs(output_dir, MODE)
    #-- output file format and output suffix for each data file format
    output_format = 'processed_{0}_{1}{2}_{3}_{4}.{5}'
    SUFFIX = dict(ascii='txt',netCDF4='nc',HDF5='h5')
    #-- spatially subset data using bounding box or polygon file
    if np.any(BBOX):
        #-- use a bounding box to spatially subset data (Xmin,Xmax,Ymin,Ymax)
        xbox = np.ravel([BBOX[0],BBOX[1],BBOX[1],BBOX[0]])
        ybox = np.ravel([BBOX[2],BBOX[2],BBOX[3],BBOX[3]])
        poly_obj = Polygon(np.concatenate((xbox[:,None],ybox[:,None]),axis=1))
        SUBSETTING = True
    elif POLYGON:
        #-- read shapefile, geojson or kml/kmz file
        fileBasename,fileExtension = os.path.splitext(POLYGON)
        #-- get MultiPolygon object from input spatial file
        if (fileExtension == '.shp'):
            m = read_shapefile(os.path.expanduser(POLYGON), EPSG=EPSG)
        elif (fileExtension == '.zip'):
            m = read_shapefile(os.path.expanduser(POLYGON), EPSG=EPSG, ZIP=True)
        elif (fileExtension == '.kml'):
            m = read_kml_file(os.path.expanduser(POLYGON), EPSG=EPSG)
        elif (fileExtension == '.kmz'):
            m = read_kmz_file(os.path.expanduser(POLYGON), EPSG=EPSG, KMZ=True)
        elif fileExtension in ('.json','.geojson'):
            m = read_geojson_file(os.path.expanduser(POLYGON), EPSG=EPSG)
        #-- calculate the convex hull of the MultiPolygon object for subsetting
        #-- to speed up computational time (versus iterating through polygons)
        poly_obj = m.convex_hull
        SUBSETTING = True
    else:
        SUBSETTING = False

    #-- number of cycles of interest
    ncycles = len(CYCLES)
    #-- find ICESat-2 ATL11 HDF5 files for parameters with regular expression
    input_files = [f for f in os.listdir(input_dir) if bool(regex.match(f))]
    #-- for each input file
    for input_file in sorted(input_files):
        #-- extract parameters from file name
        PRD,TRK,GRAN,RL,VERS = regex.findall(input_file).pop()
        #-- dictionaries with output data
        FLAG = np.zeros((3),dtype=np.bool)
        #-- x and y coordinates (in EPSG)
        D11x = {}
        D11y = {}
        #-- height, error and time
        D11h = {}
        D11e = {}
        D11t = {}
        #-- combined valid mask for cycles
        D11m = {}
        #-- for each pair track
        for i,p in enumerate([1, 2, 3]):
            f = os.path.join(input_dir,input_file)
            D11 = ATL11.data().from_file(f,pair=p).get_xy(None,EPSG=EPSG)
            #-- test that file has the cycles of interest
            npts,nc = np.shape(D11.corrected_h.h_corr)
            cycle_test = np.all([(c in np.arange(1,nc+1)) for c in CYCLES])
            #-- if subsetting with bounding box or polygon
            if SUBSETTING and cycle_test:
                #-- convert coordinates to shapely multipoint object
                D11xy = np.concatenate((D11.x[:,None],D11.y[:,None]),axis=1)
                xy_point = MultiPoint(D11xy)
                #-- finds if points are encapsulated by polygon
                int_test = poly_obj.intersects(xy_point)
                #-- if there are points encapsulated
                if int_test:
                    #-- extract intersected points
                    int_map = list(map(poly_obj.intersects,xy_point))
                    int_indices, = np.nonzero(int_map)
                    nint = np.count_nonzero(int_map)
                    #-- reduce data
                    D11x[p] = D11.x[int_indices]
                    D11y[p] = D11.y[int_indices]
                    #-- extract height and delta time at each cycle of interest
                    D11h[p] = np.zeros((nint,ncycles))
                    D11e[p] = np.zeros((nint,ncycles))
                    D11t[p] = np.zeros((nint,ncycles))
                    D11m[p] = np.ones((nint),dtype=np.bool)
                    for j,c in enumerate(CYCLES):
                        #-- 0-based indexing for cycles
                        D11h[p][:,j] = D11.corrected_h.h_corr[int_indices,c-1]
                        D11e[p][:,j] = D11.corrected_h.h_corr_sigma[int_indices,c-1]
                        D11t[p][:,j] = D11.corrected_h.delta_time[int_indices,c-1]
                        #-- set mask for cycle
                        D11m[p] &= np.isfinite(D11h[p][:,j])
                        #-- output combined data flag
                        FLAG[i] = cycle_test & int_test & np.any(D11m[p])
            elif cycle_test:
                #-- extract data
                D11x[p] = D11.x.copy()
                D11y[p] = D11.y.copy()
                #-- extract height and delta time at each cycle of interest
                D11h[p] = np.zeros((npts,ncycles))
                D11e[p] = np.zeros((npts,ncycles))
                D11t[p] = np.zeros((npts,ncycles))
                #-- for each cycle of interest
                for j,c in enumerate(CYCLES):
                    #-- 0-based indexing for cycles
                    D11h[p][:,j] = D11.corrected_h.h_corr[:,c-1]
                    D11e[p][:,j] = D11.corrected_h.h_corr_sigma[:,c-1]
                    #-- extract delta time and convert to format
                    D11t[p][:,j] = D11.corrected_h.delta_time[:,c-1]
                    #-- set mask for cycle
                    D11m[p] &= np.isfinite(D11h[p][:,j])
                    #-- output combined data flag
                    FLAG[i] = cycle_test & np.any(D11m[p])

        #-- output file if there are any valid data
        if FLAG.any():
            #-- reduce output pair tracks if data
            PAIRS = [p for i,p in enumerate([1, 2, 3]) if FLAG[i]]
            #-- output file name
            args = (PRD,TRK,GRAN,RL,VERS,SUFFIX[FORMAT])
            FILE = os.path.join(output_dir,output_format.format(*args))
            print(FILE) if VERBOSE else None
            #-- output to specified file format
            output_ATL11_file(D11x, D11y, D11h, D11e, D11t, D11m, FILENAME=FILE,
                FORMAT=FORMAT, EPSG=EPSG, TIME=TIME, PAIRS=PAIRS, MODE=MODE)

#-- PURPOSE: convert time from ATLAS SDP seconds into Julian and year-decimal
def convert_delta_time(delta_time, atlas_sdp_gps_epoch=1198800018.0):
    #-- calculate gps time from delta_time
    gps_seconds = atlas_sdp_gps_epoch + delta_time
    time_leaps = count_leap_seconds(gps_seconds)
    #-- calculate julian time
    time_julian = 2444244.5 + (gps_seconds - time_leaps)/86400.0
    #-- convert to calendar date with convert_julian.py
    Y,M,D,h,m,s = convert_julian(time_julian,FORMAT='tuple')
    #-- calculate year-decimal time
    time_decimal = convert_calendar_decimal(Y,M,DAY=D,HOUR=h,MINUTE=m,SECOND=s)
    #-- return both the Julian and year-decimal formatted dates
    return dict(julian=time_julian, decimal=time_decimal)

#-- PURPOSE: output (reduced) pair data for ATL11 files
def output_ATL11_file(D11x, D11y, D11h, D11e, D11t, D11m, FILENAME=None,
    FORMAT=None, EPSG=None, TIME=None, PAIRS=None, MODE=0o775):
    #-- output as format
    if (FORMAT == 'ascii'):
        #-- open output ascii file
        fid = open(FILENAME,'w')
        #-- for each valid pair track
        for p in PAIRS:
            #-- output dimensions (including invalids)
            npts,ncycles = np.shape(D11h[p])
            #-- convert from delta_time into formatted time
            rpt_times = convert_delta_time(D11t[p])[TIME]
            #-- find valid data
            valid, = np.nonzero(D11m[p])
            for i in valid:
                #-- merge cycles
                times = ','.join('{0:0.6f}'.format(t) for t in rpt_times[i,:])
                heights = ','.join('{0:0.6f}'.format(h) for h in D11h[p][i,:])
                errors = ','.join('{0:0.6f}'.format(h) for h in D11e[p][i,:])
                #-- print to file
                args = (D11x[p][i],D11y[p][i],times,heights,errors)
                print('{0:0.6f},{1:0.6f},{2},{3},{4}'.format(*args),file=fid)
        #--- close the file
        fid.close()

#-- run main program
if __name__ == '__main__':
	main()
