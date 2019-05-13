#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 13:16:45 2018a

@author: suzanne
"""
import os.path, csv, sys
import glob
import geopy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import Gsmooth_functions as sm
import time
import argparse
import datetime as dt
from netCDF4 import Dataset

#from mpl_toolkits.basemap import Basemap
# NEED TO get Cartopy working... see pip install Cartopy online.

class smooth_params:
    def __init__(self,args):
        # set domain for smoothing
        self.minlat = args.bbox[1]
        self.maxlat = args.bbox[3]
        self.minlon = args.bbox[0]
        self.maxlon = args.bbox[2]
        self.fwhm = args.fwhm  # full-width, half-max for smoothing, in km
        self.lon_r = np.ceil(self.fwhm/110.) * 3;
        self.lat_r = np.ceil(self.fwhm/110.) * 3;
        self.sigma = 1000. * self.fwhm / 2.355
        self.ATL = args.ATL      # 1 means Atlantic ocean, lons -180 to 180.
        
        days_in_month = {1:31,2:28,3:31,4:30,5:31,6:30,7:31,8:31,9:30,10:31,11:30,12:31}
        yf,mf = [int(item) for item in args.firstdate.split('-')]
        self.firstdate = dt.date(yf,mf,1)  
        yf,mf = [int(item) for item in args.lastdate.split('-')]
        self.lastdate = dt.date(yf,mf,days_in_month[mf])  


if __name__ == "__main__":
    #time_start = time.time()
    # read in path names for input data, putput smoothed fields and output figures.
    try:   
        with open('dataPaths.csv') as csvfile:
            rows = csv.reader(csvfile)
            data_path = {row[0]:row[1] for row in rows}  
    except:
        print('Need dataPaths.csv file containing path info.')
    
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Retrieves domain, parameters for spatial smoothing.')    
    # add parameters to parse
    parser.add_argument('-b', '--boundingbox', dest='bbox', type=float, nargs=4, required=False, help="should be of the form W S E N");
    parser.add_argument('-f', '--fwhm', type=float, help="full-width, half-max in km, default is 300",default=300);
    parser.add_argument('-fd','--firstdate', dest='firstdate', help='first date as yyyy-mm',default='1000-01-01')
    parser.add_argument('-ld','--lastdate', dest='lastdate', help='last date as yyyy-mm',default='2100-01-01')
    parser.add_argument('-a', '--atl', dest='ATL', type=int, help='1 for Atlantic Ocean, default is Pacific Ocean',default=0)
    parser.add_argument('-o', '--outdir', help='output directory for smoothed fields',default='./')
    args = parser.parse_args()
    if len(sys.argv) == 1: # if no command line arguments
        parser.print_help()
        sys.exit(1) 
        
    if args.outdir:
        sm.establish_dir(data_path['smoo_dir']+ args.outdir.strip('/') + '/')

    # user defined parameters, depending on where to smooth and how much smoothing (fwhm) to do.
    params = smooth_params(args)
#    attrs = vars(params)
#    print('\n'.join("%s: %s" % item for item in attrs.items()))
        
    # filename containing the smoothing weights, based on params
    file_wts = data_path['data_dir'] + 'smooth_wts_avisogrid_fwhm' + str(np.int(params.fwhm)) + '.nc'
    print(file_wts)
    print()
    # bring in variable to smooth
    filelist=sorted(glob.glob(data_path['data_dir'] + 'dt*.nc'))
    filelist = sm.files_within_dates(filelist,params)
    [print(item) for item in filelist]

    tt=np.array([])
    sla=np.array([])

    for i,file in enumerate(filelist):
        nc = Dataset(file,mode='r')
#        print(nc.file_format)
#        print(nc.dimensions)  # also nc.dimensions.keys()
#        print(nc.variables)
#        print(nc.variables['time'].units)  
        # append time to tt variable
        tt = np.append(tt,nc.variables['time'][:])
        if i==0:
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            sla=np.squeeze(nc.variables['sla'][:])  # squeeze because was [1,720,1440]
        else:
            # append geophysical variable to sla variable
            sla = np.dstack((sla,np.squeeze(nc.variables['sla'][:]))) 

    # replace FillValue with nans            
    np.place(sla,sla==nc.variables['sla']._FillValue,np.nan)

    # change longitudes to -180 to 180 if working in Atlantic
    if params.ATL == 1:
        np.place(lon,lon>180,lon[lon>180]-360)
        lon=np.roll(lon,len(lon)//2)                        
        sla=np.roll(sla,len(lon)//2,axis=1)  # sla dims: lat, lon, all-months
        sla = np.moveaxis(sla,-1,0)  # move last dim to front: now all-months, lat, lon
    # make land mask from the means of sla over time. Ocean = 1, Land = nan
    mask = np.mean(sla,axis=0)
    mask[~np.isnan(mask)]=1

    levels = np.arange(-0.3,0.31,0.01)
#    plt.figure(1);plt.clf()
#    plt.contourf(lon[np.logical_and(lon>=params.minlon,lon<=params.maxlon)],lat[np.logical_and(lat>=params.minlat,lat<=params.maxlat)],
#                     sla[:,np.logical_and(lat>=params.minlat,lat<=params.maxlat),:][:,:,np.logical_and(lon>=params.minlon,lon<=params.maxlon)][0,:,:],levels);plt.colorbar()
#    plt.title('Sea level anomaly');plt.show()

    # check if smoothing weights file already exists for particular grid and smoothing parameter 'fwhm'
    if os.path.isfile(file_wts):
        print('smoothing weights file exists')
        nc = Dataset(file_wts,mode='r')
        print('nc vars',nc.variables.keys())
        class load_wts:
            def __init__(self,ncDataSet):
                for k in ncDataSet.variables.keys(): 
                    setattr(self,k,ncDataSet.variables[k][:])
        
        Gauss = load_wts(nc)
                
    else:
        print('make smoothing weights file!')
        Gauss = sm.smoothing_weights(lat,lon,params).make_weights()        
        filename = data_path['data_dir'] + 'smooth_wts_avisogrid_fwhm {0%d} .nc'.format(params.fwhm)
        print(filename)
        nc = Dataset(filename,'w',format='NETCDF4_CLASSIC')
        
        # create dimensions
        sm_width = nc.createDimension('smoothing_width',Gauss.lat_wts.shape[0]) # None means unlimitied 
        lat_center = nc.createDimension('lat_center',Gauss.wts.shape[2])
        lon_center = nc.createDimension('lon_center',1)
        print(nc.dimensions.keys())
        
        # create variables
        var=nc.createVariable('lat_ctr',float,('lat_center',))
        var[:]=Gauss.lat_ctr
        var=nc.createVariable('lon_ctr',float,('lon_center',))
        var[:]=Gauss.lon_ctr
        var=nc.createVariable('lat_wts',float,('smoothing_width','lat_center'),fill_value='NaN')
        var[:]=Gauss.lat_wts
        var=nc.createVariable('lon_wts',float,('smoothing_width',))
        var[:]=Gauss.lon_wts
        var=nc.createVariable('wts',float,('smoothing_width','smoothing_width','lat_center'),fill_value='NaN')
        var[:]=Gauss.wts
        nc.close()
    
    
    # time the smoothing
    #t_start = time.perf_counter()
    print('shapes', len(lat),len(lon),len(sla)) #,lat.shape(),lon.shape(),Gauss.shape())
    # smooth the monthly clim anomaly, geo, where geo = var minus monthly mean of var
    sla_sm, lat_sm, lon_sm = sm.Gsmooth(sla,lat,lon,tt,mask,params,Gauss).smooth_it()  # Gsmooth doesn't use tt 
    print('sla_sm size',sla_sm.shape)
    fig = plt.figure(33);plt.clf()
    levels = np.arange(-0.15,0.16,0.01)
    plt.contourf(lon_sm,lat_sm,sla_sm[0,:,:],levels);plt.colorbar()
    plt.title('smoothed AVISO, SSH anomaly minus clim means, Jan 2014')
    fig.savefig('figs/smoothed_SSH.jpg')
    t_stop = time.perf_counter()
    #print('time perf',(t_stop-t_start)/60)
    
    # parse the tt vector ('iterable' in python parlance)
    year  = np.array([ (dt.datetime(1950,1,1) + dt.timedelta(days=np.asscalar(tt[i]))).year for i in range(len(tt))])
    month = np.array([ (dt.datetime(1950,1,1) + dt.timedelta(days=np.asscalar(tt[i]))).month for i in range(len(tt))])
    day   = np.array([ (dt.datetime(1950,1,1) + dt.timedelta(days=np.asscalar(tt[i]))).day for i in range(len(tt))])
    
    # make netcdf file
    filename = data_path['smoo_dir'] + 'test_sla_sm.nc'
    print(filename)
    nc = Dataset(filename,'w',format='NETCDF4_CLASSIC')
    #print(nc.file_format) # CLASSIC data: dims, variables, attributes.
    attrs = vars(params)
    for x in attrs.items():
        if 'date' in x[0]:
            nc.setncattr(x[0],str(x[1]))
        else:
            nc.setncattr(x[0],x[1])
        
    # create dimensions
    latitude = nc.createDimension('latitude',lat_sm.shape[0]) # 'None' indicates unlimitied 
    longitude = nc.createDimension('longitude',lon_sm.shape[0]) 
    time = nc.createDimension('time',tt.shape[0])
    print(nc.dimensions.keys())
    
    # include params in output file
    
    #var=nc.createVariable()
    # create variables
    var=nc.createVariable('lat',np.float32,('latitude',))  # variable id, type, dimension(s)
    var[:]=lat_sm
    var.units = 'degrees_North'
    var=nc.createVariable('lon',np.float32,('longitude',))
    var[:]=lon_sm
    var.units = 'degrees_East'
    var=nc.createVariable('year',np.int32,('time',))
    var[:]=year
    var=nc.createVariable('month',np.int32,('time',))
    var[:]=month
    var=nc.createVariable('day',np.int32,('time',))
    var[:]=day
    
    var=nc.createVariable('ssh_sm',np.float32,('time','latitude','longitude'))
    var[:]=sla_sm
    var.units = 'meters'
    var.land_value = 'NaN'
    var.grid = 'Aviso grid: 1/4 degree in lat and lon, centered on the eighth degree'
    
    # global attribute
    nc.fwhm = 'Full-width, half-max is ' + str(params.fwhm) + ' km'
    nc.close()
    
    #print('this program took', time.time() - time_start(), 'seconds to run')
