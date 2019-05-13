#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 13:16:45 2018a

@author: suzanne
"""
import os.path, csv, sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import Gsmooth_functions as sm
import time
import re
#from mpl_toolkits.basemap import Basemap
# NEED TO get Cartopy working... see pip install Cartopy online.
import datetime as dt
from calendar import monthrange
from netCDF4 import Dataset

data_dir = '/Users/suzanne/luanne/cmip5/FiguresforPaper/py/data/'
figs_dir = '/Users/suzanne/luanne/cmip5/FiguresforPaper/py/figs/'
smoo_dir = '/Users/suzanne/luanne/cmip5/FiguresforPaper/py/smoothed/'

# user defined parameters, depending on where to smooth and how much smoothing (fwhm) to do.
class smooth_params:
    def __init__(self,inPar):
        # set domain for smoothing
        for k,v in inPar.items():
            if 'lat' in k or 'lon' in k or 'fw' in k:
                setattr(self, k, float(v))
            elif 'date' in k:
                yf,mf = [int(item) for item in inPar['firstdate'].split('-')]
                setattr(self, 'firstdate', dt.date(yf,mf,1))
                yf,mf = [int(item) for item in inPar['lastdate'].split('-')]
                setattr(self, 'lastdate', dt.date(yf,mf,monthrange(yf,mf)[1]))                     
#            else:
#                setattr(self, k, v)
        self.lon_r = np.ceil(self.fwhm/110.) * 3;
        self.lat_r = np.ceil(self.fwhm/110.) * 3;
        self.sigma = 1000. * self.fwhm / 2.355
            
if __name__ == "__main__":
    # read in path names for input data, putput smoothed fields and output figures.
    #inputParams = {}
    try:   
        with open('inputParams.csv') as csvfile:
            rows = csv.reader(csvfile)
            inputParams= {row[0]:row[1].strip('\t') for row in rows} 
            #inputParams= {row[0]:row[1].strip('\t') for row in rows[6:]}  
    except:
        print('Need inputParams.csv file to proceed.')

    params = smooth_params(inputParams)
    attrs = vars(params)
    print('\n'.join("%s: %s" % item for item in attrs.items()))

    # filename containing (or to be made of) the smoothing weights
    #print(inputParams['data_dir']  'smooth_wts_' inputParams['datasrc'] + 'grid_fwhm' + str(np.int(params.fwhm)) + '.nc')
    file_wts = inputParams['data_dir'] + 'smooth_wts_' + inputParams['datasrc'] + 'grid_fwhm' + str(np.int(params.fwhm)) + '.nc'
    print(file_wts)
    
    # bring in latent heat to smooth
    filelist=sorted(glob.glob(inputParams['data_dir'] + inputParams['filebase']  + '*.nc'))
    print('data files',filelist)

    filelist = sm.files_within_dates(filelist,params)
    [print(item) for item in filelist]
    month=np.array([])
    year =np.array([])
    day  =np.array([])
    
    for i,file in enumerate(filelist):
        print('file',file)
        nc = Dataset(file,mode='r')
        # append time to mm variable (time is actually month of year)
        if 'oaflux' in inputParams['datasrc']:
            month = np.append(month, nc.variables['time'][:])
            year  = np.append(year, np.ones([12,1],dtype=np.int)*np.int(re.findall(r'\d+',file)[1]))
            day = np.append(day, np.ones([12,1])*15)
            if i==0:
                lat = nc.variables['lat'][:]
                lat_units = nc.variables['lat'].units
                lon = nc.variables['lon'][:]
                lon_units = nc.variables['lon'].units
                origVar = np.empty((len(filelist)*12,len(lat),len(lon)))
                origVar[range(i*12,i*12+12),:,:] = -1 * nc.variables[inputParams['varname']][:]
                units_str = nc.variables[inputParams['varname']].units
            else:
                origVar[range(i*12,i*12+12),:,:] = -1 * nc.variables[inputParams['varname']][:]
            
        elif 'aviso' in inputParams['datasrc']:
            filedate = dt.timedelta(days=nc.variables['time'][:]) + dt.date(1950,1,1)
            year = np.append(year,filedate.year)
            month = np.append(month,filedate.month)
            day = np.append(day,filedate.day)
            if i==0:
                lat = nc.variables['lat'][:]
                lat_units = nc.variables['lat'].units
                lon = nc.variables['lon'][:]
                lon_units = nc.variables['lon'].units
                origVar = np.empty((len(filelist),len(lat),len(lon)))
                origVar[i,:,:]=nc.variables[inputParams['varname']][:] 
                units_str = nc.variables[inputParams['varname']].units
            else:
                origVar[i,:,:]=nc.variables[inputParams['varname']][:] 
    
    # make time iterable
    tt = [dt.date(np.int(yy),np.int(mm),np.int(dd)) for yy,mm,dd in zip(year,month,day)]
    print('tt',origVar.shape)
    # replace fillValues with nans            
    np.place(origVar,origVar==nc.variables[inputParams['varname']].missing_value,np.nan)
    
    # change longitudes to -180 to 180 if working in Atlantic
    if inputParams['ATL']:
        np.place(lon,lon>180,lon[lon>180]-360)
        lon=np.roll(lon,len(lon)//2)                        
        origVar=np.roll(origVar,len(lon)//2,axis=2)  # dims: all-months, lat, lon
        
    
    # make land mask from the means of geovar over time. Ocean = 1, Land = nan
    mask = np.mean(origVar,axis=0)
    mask[mask==-1]=np.nan
    mask[~np.isnan(mask)]=1
    #plt.figure(2);plt.clf()
    #plt.imshow(mask);plt.colorbar()
    #plt.figure(11);plt.clf()
    #plt.imshow(lh[:,np.logical_and(lat>=params.minlat,lat<=params.maxlat),:][:,:,np.logical_and(lon>=params.minlon,lon<=params.maxlon)][0,:,:]*
    #           mask[np.logical_and(lat>=params.minlat,lat<=params.maxlat),:][:,np.logical_and(lon>=params.minlon,lon<=params.maxlon)],vmin=-50,vmax=30);plt.colorbar()
    levels = np.arange(-50,35,5)
    print(params.minlon,params.maxlon,params.minlat,params.maxlat)
#    plt.contourf(lon[np.logical_and(lon>=params.minlon,lon<=params.maxlon)],lat[np.logical_and(lat>=params.minlat,lat<=params.maxlat)],
#                 lh[:,np.logical_and(lat>=params.minlat,lat<=params.maxlat),:][:,:,np.logical_and(lon>=params.minlon,lon<=params.maxlon)][0,:,:]);plt.colorbar()
    #plt.contourf(lon[np.logical_and(lon>=params.minlon,lon<=params.maxlon)],lat[np.logical_and(lat>=params.minlat,lat<=params.maxlat)],
    #             lh[:,np.logical_and(lat>=params.minlat,lat<=params.maxlat),:][:,:,np.logical_and(lon>=params.minlon,lon<=params.maxlon)][0,:,:]*
    #             mask[np.logical_and(lat>=params.minlat,lat<=params.maxlat),:][:,np.logical_and(lon>=params.minlon,lon<=params.maxlon)],levels,colors='k');plt.colorbar()
    
    
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
        #print('shapes',Gauss.lat_wts.shape,Gauss.lon_wts.shape,Gauss.wts.shape)
        
        filename = '{0}smooth_wts_{1}grid_fwhm{2}.nc'.format(inputParams['data_dir'],inputParams['datasrc'],int(params.fwhm))
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
    
    # time the smoothing process
    t_start = time.perf_counter()
    print(origVar.shape)
    # G.smooth smoothes the monthly clim anomaly, geo, where geo = var minus monthly mean of var
    smooVar, lat_sm, lon_sm = sm.Gsmooth(origVar,lat,lon,tt,mask,params,Gauss).smooth_it()  # Gsmooth doesn't use tt 
    
    print('shapes after smoothing',smooVar.shape,lat_sm.shape,lon_sm.shape)
#    fig = plt.figure(33);plt.clf()
#    plt.contourf(lon[np.logical_and(lon>=params.minlon,lon<=params.maxlon)],lat[np.logical_and(lat>=params.minlat,lat<=params.maxlat)],lh_sm[0,:,:],levels);plt.colorbar()
#    plt.title('smoothed OAFLUX, LH anomaly, Jan 2014')
#    fig.savefig(figs_dir + 'smoo_LH.jpg')
    
    t_stop = time.perf_counter()
    print('time perf',(t_stop-t_start)/60)  # prints minutes
    #year  = np.array([ (dt.datetime(1950,1,1) + dt.timedelta(days=np.asscalar(tt[i]))).year for i in range(len(tt))])
    #month = np.array([ (dt.datetime(1950,1,1) + dt.timedelta(days=np.asscalar(tt[i]))).month for i in range(len(tt))])
    #day   = np.array([ (dt.datetime(1950,1,1) + dt.timedelta(days=np.asscalar(tt[i]))).day for i in range(len(tt))])
    #
    # make netcdf file
    filename=inputParams['smoo_dir'] + 'test_' + inputParams['geovar'] + '_sm.nc'
    print(filename)
    nc = Dataset(filename,'w',format='NETCDF4_CLASSIC')
    print(nc.file_format) # CLASSIC data: dims, variables, attributes.
    # include attributes from smoothing parameters as global attributes of the netcdf file
    attrs = vars(params)
    for x in attrs.items():
        print('attrs',x[0],x[1])
        if 'date' in x[0]:
            nc.setncattr(x[0],str(x[1]))
        elif 'fwhm' in x[0]:
            nc.setncattr(x[0],'Full-width, half-max is ' + str(np.int(params.fwhm)) + ' km')   
        else:
            nc.setncattr(x[0],x[1])
    
    # create dimensions
    latitude = nc.createDimension('latitude',lat_sm.shape[0]) # None means unlimitied 
    longitude = nc.createDimension('longitude',lon_sm.shape[0]) 
    time = nc.createDimension('time',len(tt))
    print(nc.dimensions.keys())
    
    # create variables
    var=nc.createVariable('lat',np.float32,('latitude',))  # variable id, type, dimension(s)
    var[:]=lat_sm
    var.units = lat_units
    var=nc.createVariable('lon',np.float32,('longitude',))
    var[:]=lon_sm
    var.units = lon_units
    var=nc.createVariable('year',np.int32,('time',))
    var[:]=year
    var=nc.createVariable('month',np.int32,('time',))
    var[:]=month
    var=nc.createVariable('day',np.int32,('time',))
    var[:]=day
    
    var=nc.createVariable(inputParams['geovar']+'_sm',np.float32,('time','latitude','longitude'))
    var[:]=smooVar
    var.units = units_str
    var.land_value = 'NaN'
    if 'oaflux' in inputParams['datasrc']:
        var.grid = 'OAFlux grid: 1 degree in lat and lon, centered on the half degree'
    elif 'aviso' in inputParams['datasrc']:
        var.grid = 'Aviso grid: 1/4 degree in lat and lon, centered on the eighth degree'
    
    nc.close()
    
    

## bring in sea surface temperature to smooth
#filelist=sorted(glob.glob(data_dir + 'ts_oaflux_*.nc'))
#month=np.array([])
#year =np.array([])
#day  =np.array([])
#sst  =np.array([])
#
#
#for i,file in enumerate(filelist):
#    print('file',file)
#    nc = Dataset(file,mode='r')
#    # append time to mm variable (time is actually month of year)
#    month = np.append(month, nc.variables['time'][:])
#    year  = np.append(year, np.ones([12,1],dtype=np.int)*np.int(re.findall(r'\d+',file)[1]))
#    day = np.append(day, np.ones([12,1])*15)
#    if i==0:
#        lat = nc.variables['lat'][:]
#        lon = nc.variables['lon'][:]
#        sst = nc.variables['tmpsf'][:] 
#    else:
#        # append geophysical variable 
#        sst = np.concatenate((sst,nc.variables['tmpsf'][:]),axis=0) 
#        
## make time iterable
#tt = [dt.date(np.int(yy),np.int(mm),np.int(dd)) for yy,mm,dd in zip(year,month,day)]
#
## replace fillValues with nans            
#np.place(sst,sst==nc.variables['tmpsf'].missing_value,np.nan)
#
## change longitudes to -180 to 180 if working in Atlantic
#if params.ATL == 1:
#    np.place(lon,lon>180,lon[lon>180]-360)
#    lon=np.roll(lon,len(lon)//2)                        
#    sst=np.roll(sst,len(lon)//2,axis=2)  # dims: all-months, lat, lon
#    
## make land mask from the means of geovar over time. Ocean = 1, Land = nan
#mask = np.mean(sst,axis=0)
#mask[mask==-1]=np.nan
#mask[~np.isnan(mask)]=1
#
#plt.figure(11);plt.clf()
#levels = np.arange(24,38,2)
#plt.contourf(lon[np.logical_and(lon>=params.minlon,lon<=params.maxlon)],lat[np.logical_and(lat>=params.minlat,lat<=params.maxlat)],
#             sst[:,np.logical_and(lat>=params.minlat,lat<=params.maxlat),:][:,:,np.logical_and(lon>=params.minlon,lon<=params.maxlon)][0,:,:],levels,colors='k');plt.colorbar()
##plt.contourf(lon[np.logical_and(lon>=params.minlon,lon<=params.maxlon)],lat[np.logical_and(lat>=params.minlat,lat<=params.maxlat)],
##             sst[:,np.logical_and(lat>=params.minlat,lat<=params.maxlat),:][:,:,np.logical_and(lon>=params.minlon,lon<=params.maxlon)][0,:,:]*
##             mask[np.logical_and(lat>=params.minlat,lat<=params.maxlat),:][:,np.logical_and(lon>=params.minlon,lon<=params.maxlon)],levels,colors='k');plt.colorbar()
#
#
## check if smoothing weights file already exists for particular grid and smoothing parameter 'fwhm'
#if os.path.isfile(file_wts):
#    print('smoothing weights file exists')
#    nc = Dataset(file_wts,mode='r')
#    # I don't know how to simply put an attribute on a class that's not define, or do both in one stepp.
#    class Object(object):
#        pass
#    Gauss = Object()        
#    Gauss.wts = nc.variables['wts'][:]
#    Gauss.lat_wts = nc.variables['lat_wts'][:]
#    Gauss.lon_wts = nc.variables['lon_wts'][:]
#    Gauss.lat_ctr = nc.variables['lat_ctr'][:]
#    Gauss.lon_ctr = nc.variables['lon_ctr'][:]
#    
#else:
#    print('make smoothing weights file!')
#    Gauss = sm.smoothing_weights(lat,lon,params).make_weights()
#    
#    filename='/Users/suzanne/luanne/cmip5/FiguresforPaper/py/smooth_wts_oafluxgrid_fwhm' + str(np.int(params.fwhm)) + '.nc'
#    print(filename)
#    nc = Dataset(filename,'w',format='NETCDF4_CLASSIC')
#    #print(nc.file_format) # CLASSIC data: dims, variables, attributes.
#    
#    # create dimensions
#    sm_width = nc.createDimension('smoothing_width',Gauss.lat_wts.shape[0]) # None means unlimitied 
#    lat_center = nc.createDimension('lat_center',Gauss.wts.shape[2])
#    lon_center = nc.createDimension('lon_center',1)
#    print(nc.dimensions.keys())
#    
#    # create variables
#    var=nc.createVariable('lat_ctr',float,('lat_center',))
#    var[:]=Gauss.lat_ctr
#    var=nc.createVariable('lon_ctr',float,('lon_center',))
#    var[:]=Gauss.lon_ctr
#    var=nc.createVariable('lat_wts',float,('smoothing_width','lat_center'),fill_value='NaN')
#    var[:]=Gauss.lat_wts
#    var=nc.createVariable('lon_wts',float,('smoothing_width',))
#    var[:]=Gauss.lon_wts
#    var=nc.createVariable('wts',float,('smoothing_width','smoothing_width','lat_center'),fill_value='NaN')
#    var[:]=Gauss.wts
#    nc.close()
#
## time the smoothing process
#t_start = time.perf_counter()
#
## G.smooth smoothes the monthly clim anomaly, geo, where geo = var minus monthly mean of var
#sst_sm, lat_sm, lon_sm = sm.Gsmooth(sst,lat,lon,tt,mask,params,Gauss).smooth_it()  # Gsmooth doesn't use tt 
#
#print('shapes after smoothing',sst_sm.shape,lat_sm.shape,lon_sm.shape)
#levels = np.arange(-1.5,1,0.25)
#fig = plt.figure(33);plt.clf()
#plt.contourf(lon[np.logical_and(lon>=params.minlon,lon<=params.maxlon)],lat[np.logical_and(lat>=params.minlat,lat<=params.maxlat)],sst_sm[0,:,:],levels);plt.colorbar()
#plt.title('smoothed OAFLUX, SST anomaly, Jan 2014')
#fig.savefig(figs_dir + 'smoo_SST.jpg')
#
#t_stop = time.perf_counter()
#print('time perf',(t_stop-t_start)/60)  # prints minutes
##year  = np.array([ (dt.datetime(1950,1,1) + dt.timedelta(days=np.asscalar(tt[i]))).year for i in range(len(tt))])
##month = np.array([ (dt.datetime(1950,1,1) + dt.timedelta(days=np.asscalar(tt[i]))).month for i in range(len(tt))])
##day   = np.array([ (dt.datetime(1950,1,1) + dt.timedelta(days=np.asscalar(tt[i]))).day for i in range(len(tt))])
## make netcdf file
#filename=smoo_dir + 'test_sst_sm.nc'
#print(filename)
#nc = Dataset(filename,'w',format='NETCDF4_CLASSIC')
#print(nc.file_format) # CLASSIC data: dims, variables, attributes.
#
## create dimensions
#latitude = nc.createDimension('latitude',lat_sm.shape[0]) # None means unlimitied 
#longitude = nc.createDimension('longitude',lon_sm.shape[0]) 
#time = nc.createDimension('time',len(tt))
#print(nc.dimensions.keys())
#
## create variables
#var=nc.createVariable('lat',np.float32,('latitude',))  # variable id, type, dimension(s)
#var[:]=lat_sm
#var.units = 'degrees_North'
#var=nc.createVariable('lon',np.float32,('longitude',))
#var[:]=lon_sm
#var.units = 'degrees_East'
#var=nc.createVariable('year',np.int32,('time',))
#var[:]=year
#var=nc.createVariable('month',np.int32,('time',))
#var[:]=month
#var=nc.createVariable('day',np.int32,('time',))
#var[:]=day
#
#var=nc.createVariable('sst_sm',np.float32,('time','latitude','longitude'))
#var[:]=sst_sm
#var.units = 'deg C'
#var.land_value = 'NaN'
#var.grid = 'OAFlux grid: 1 degree in lat and lon, centered on the half degree'
#
## global attribute
#nc.fwhm = 'Full-width, half-max is ' + str(np.int(params.fwhm)) + ' km'
#nc.close()
