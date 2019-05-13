#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 26 13:16:45 2018a

@author: suzanne
"""
import numpy as np
import scipy as sp
from scipy import stats
from scipy import interpolate
from geopy import distance
import matplotlib.pyplot as plt
import os, re
import datetime as dt
#from mpl_toolkits.basemap import Basemap
# NEED TO get Cartopy working... see pip install Cartopy online.

# object definition containing the smoothing weights if there is no file for it already.
class smoothing_weights:
    def __init__(self,lat,lon,params):
        self.lat = lat
        self.lon = lon
        self.latd= lat #lat[np.logical_and((lat>=params.minlat),(lat<=params.maxlat))]
        #self.lond= lon lon[np.logical_and((lon>=params.minlon),(lon<=params.maxlon))]
        self.params = params
        
    def make_weights(self):
        self.lat_ctr = self.lat
        self.lon_ctr = self.lon[self.lon.shape[0]//2]
        print('lat ctr, lon ctr',self.lat_ctr,self.lon_ctr)
        lon1 = self.lon[np.abs(self.lon_ctr-self.lon)<self.params.lon_r]
        if self.lat.shape[0] % 2 == 0:
            indla=self.lat.shape[0]//2
        else:
            indla=self.lat.shape[0]//2 + 1
        Gauss_wts = np.full([lon1.shape[0],lon1.shape[0],indla],np.nan)
        lat_wts   = np.full([lon1.shape[0],              indla],np.nan)
        lon_wts   = lon1

        print(' in make_weights',self.latd,indla,lat_wts.shape)
        # compute gaussian weights for the first half of the latitude grid 
        for jj,la in enumerate(self.latd[0:indla]):
            print('lat center = {:3.3f}, lon center = {:1.3f}'.format(la,self.lon_ctr),jj,lat_wts.shape)
            lat1 = self.lat[np.abs(la-self.lat)<self.params.lat_r]
            #print(self.params.lat_r,lat1)
            np.put(lat_wts[:,jj],range(lat1.shape[0]),lat1)
            dd = np.array([distance.distance((lat1[j],lon1[i]),(la,self.lon_ctr)).m for j in range(len(lat1)) for i in range(len(lon1))]).reshape(len(lat1),len(lon1))
            Gauss_wts[range(lat1.shape[0]),:,jj] = 1/(2*np.pi*self.params.sigma**2) * np.exp(-dd**2/(2*self.params.sigma**2))
            #dd_wts[range(lat1.shape[0]),:,jj] = dd
#        plt.figure(56);plt.clf()
#        plt.subplot(211)
#        plt.imshow(lat_wts,vmin=-80,vmax=80);plt.colorbar()
        # concatenate with a mirror image for the second hemiphsere
        if self.lat.shape[0] % 2 == 0:
            lat_wts = np.concatenate([lat_wts,-1*np.flipud(np.fliplr(lat_wts))],axis=1)
            Gauss_wts = np.concatenate([Gauss_wts,np.moveaxis(np.array([np.flipud(Gauss_wts[:,:,jj]) for jj in range(indla-1,-1,-1)]),0,-1)],axis=2)
            
        else:  # if lat vector has a length that is odd, HAS NOT BEEN CHECKED!!!!!
            print('you are in uncharted terr')
            lat_wts   = np.concatenate([lat_wts,-1*np.flipud(np.fliplr(lat_wts[0:-1]))],axis=1)
            Gauss_wts = np.concatenate([Gauss_wts,np.moveaxis(np.array([np.flipud(Gauss_wts[:,:,jj]) for jj in range(indla-2,-1,-1)]),0,1)],axis=2)
            #dd_wts    = np.concatenate([dd_wts,np.moveaxis(np.array([np.flipud(dd_wts[:,:,jj]) for jj in range(indla-2,-1,-1)]),0,1)],axis=2)
#        print(Gauss_wts.shape,dd.shape)
#        plt.figure(56)
#        plt.subplot(212)
#        plt.imshow(lat_wts,vmin=-80,vmax=80);plt.colorbar()
#        plt.figure(55);plt.clf()
#        plt.imshow(Gauss_wts[:,:,90]);plt.colorbar()
#        plt.figure(54);plt.clf()
#        plt.imshow(dd);plt.colorbar()
        
#        class Object(object):
#            pass
#        Gauss = Object()        
#        Gauss.wts = Gauss_wts
#        Gauss.lat_wts = lat_wts
#        Gauss.lon_wts = lon_wts
#        Gauss.lat_ctr = self.lat_ctr
#        Gauss.lon_ctr = self.lon_ctr
#        return  Gauss
        self.wts = Gauss_wts
        self.lat_wts = lat_wts
        self.lon_wts = lon_wts
        return self
                
# find dates in directory names and trim to firstdate/lastdate
def files_within_dates(flist,params):
    t0=[]
    for f in flist:
        fb = os.path.basename(f)
        if len(re.findall(r'\d+',fb)) == 1:  # oaflux 
            yr = re.findall(r'\d+',fb)
            t0.append(dt.date(int(yr[0]),1,1))
        elif len(re.findall(r'\d+',fb)) == 2:  # aviso
            yr, mo = re.findall(r'\d+',fb)
            t0.append(dt.date(int(yr),int(mo),1))

    mask = [t0[ii]>=params.firstdate and t0[ii]<=params.lastdate for ii in range(len(t0))]
    return [flist[ii] for ii in range(len(flist)) if mask[ii]]    
        

def establish_dir(download_dir):
   if not os.path.exists(download_dir):
        os.makedirs(download_dir)
        if os.access(download_dir, os.W_OK) is False:  # I haven't checked this one
            print ("WARNING: Cannot write to this path! Check permissions for {0}".format(download_dir))
            exit(-1)
    
# object definiton for smoothing a variable over a specified domain, grid and smoothing distance.
class Gsmooth:
    def __init__(self,geo,lat,lon,time,mask,params,Gauss):
        self.lat = lat
        self.lon = lon
        self.time = time
        #### could check here that time is a factor of 12.
        self.mask = mask
        self.params = params
        self.Gauss = Gauss
        # reshape the total months, to years and months
        self.geo = geo.reshape((geo.shape[0]//12,12,geo.shape[1],geo.shape[2])) # months, years, lats, lons 
        self.nyr, self.nmo = self.geo.shape[0:2]
        # compute climatology, mean over years 
        geo_clim=np.ndarray.mean(self.geo,axis=0)
        # remove clim from geo to get anomalies
        # need to 'moveaxis' after list comprehension, as the iterable axis ends up as the 0th dimension
        #self.geo = np.moveaxis(np.array([self.geo[ii,:,:,:]-geo_clim for ii in range(self.nyr)]),0,1)
        self.geo = np.array([self.geo[ii,:,:,:]-geo_clim for ii in range(self.nyr)])
        # find the domain over which to smooth, get sizes
        self.latd= lat[np.logical_and((lat>=self.params.minlat),(lat<=self.params.maxlat))]
        self.lond= lon[np.logical_and((lon>=self.params.minlon),(lon<=self.params.maxlon))]
        self.nlat = self.latd.shape[0]
        self.nlon = self.lond.shape[0]

#        fig = plt.figure(22);plt.clf()
#        levels = np.arange(-1.5,1,0.25)
#        plt.contourf(self.lond,self.latd,self.geo[:,:,np.in1d(self.lat,self.latd),:][:,:,:,np.in1d(self.lon,self.lond)][0,0,:,:]*
#                     self.mask[np.in1d(self.lat,self.latd),:][:,np.in1d(self.lon,self.lond)],levels);plt.colorbar()
#        plt.title('original OAFLUX, SST anomaly, Jan 2014')
#        fig.savefig('figs/orig_SST.jpg')
#        plt.title('original OAFLUX, LH anomaly, Jan 2014')
#        fig.savefig('figs/orig_LH.jpg')
#        plt.title('original AVISO, SSH anomaly, Jan 2014')
#        fig.savefig('figs/orig_SSH.jpg')

    def smooth_it(self):
        geo_sm = np.full([self.nyr,self.nmo,self.nlat,self.nlon],np.nan)
#        geo_smmask = np.full([self.nmo,self.nyr,self.nlat,self.nlon],np.nan)
        for j in range(self.nlat):
            la = self.Gauss.lat_wts[:,self.latd[j]==self.Gauss.lat_ctr].ravel()
            #print('smoothing along ',self.latd[j])
            for i in range(self.nlon): 
                lo = self.Gauss.lon_wts + self.lond[i] - self.Gauss.lon_ctr
                mask_trunc = self.mask[np.in1d(self.lat,la),:][:,np.in1d(self.lon,lo)]
#                plt.figure(40);plt.clf()
#                plt.imshow(mask_trunc);#plt.colorbar()
   
                trunc = self.geo[:,:,np.in1d(self.lat,la),:][:,:,:,np.in1d(self.lon,lo)]  # likes to slice 1-d at a time. np.in1d booleam mask
#                plt.figure(42);plt.clf()
#                plt.imshow(np.squeeze(trunc[0,0,:,:]),vmin=-0.1,vmax=0.1);plt.colorbar()
                
                wts_trunc = np.squeeze(self.Gauss.wts[:,:,np.in1d(self.Gauss.lat_ctr,self.latd[j])])
#                plt.figure(43);plt.clf()
#                plt.imshow(wts_trunc);plt.colorbar()
                # mask out the land for summing below
                wts_trunc = wts_trunc * mask_trunc
#                plt.figure(44);plt.clf()
#                plt.imshow(wts_trunc);plt.colorbar()
                if np.nansum(wts_trunc.ravel())!=0:
                    geo_sm[:,:,j,i] = np.array([np.nansum(trunc[yr,mo,:,:].ravel() * wts_trunc.ravel()) / np.nansum(wts_trunc.ravel()) for yr in range(self.nyr) for mo in range(self.nmo)]).reshape(self.nyr,self.nmo) 
            
#        plt.figure(3);plt.clf()
#        plt.imshow(geo_sm[0,0,:,:],vmin=-0.1,vmax=0.1);plt.colorbar()
        # apply mask
        geo_sm = np.array([geo_sm[yr,mo,:,:] * self.mask[np.in1d(self.lat,self.latd),:][:,np.in1d(self.lon,self.lond)] for yr in range(self.nyr) for mo in range(self.nmo)]).reshape(self.nyr,self.nmo,self.nlat,self.nlon)

#        plt.figure(13);plt.clf()
#        plt.imshow(geo_sm[0,0,:,:],vmin=-0.1,vmax=0.1);plt.colorbar()
#        plt.figure(23);plt.clf()
#        plt.imshow(self.mask[np.in1d(self.lat,self.latd),:][:,np.in1d(self.lon,self.lond)],vmin=-0.1,vmax=0.1);plt.colorbar()
#        print('wtf',np.isnan(geo_sm).any())
        
        
        return geo_sm.reshape(self.nyr*self.nmo,self.nlat,self.nlon),self.latd,self.lond 
        #return geo_sm,self.latd,self.lond 
        
def correlatem(xx,yy):
    xx = xx.ravel()
    yy = yy.ravel()
    mx = len(xx)
    my = len(yy);print('lengths',mx,my)
    # compute correlation
    xp = xx-np.mean(xx)
    yp = yy-np.mean(yy);
    rho = np.corrcoef(xp,yp)[1,0] # rho[0,0] xp with itself, rho[1,0] xp with yp, rho[0,1] yp with xp, rho[1,1] yp with itself
    # compute covariance
    lags = np.arange(-(mx-1),mx)    
    covx = np.correlate(xp,xp,mode='full')  # get all lags (2*mx -1)
    covx = covx / (mx - np.abs(lags));      # like 'unbiased' scaling in matlab
    covx = covx / np.max(covx)
    covy = np.correlate(yp,yp,mode='full')
    covy = covy / (my - np.abs(lags)); 
    covy = covy / np.max(covy)

    # integrate until zero crossings
    # USE TRAPEZOID RULE INSTEAD OF RECTANGLE
    taux=0.;
    i=mx-1;
    while covx[i]>=0:
    	taux=taux+(covx[i]+covx[i+1])/2;
    	i=i+1;
    tauy=0.;
    i=my-1;
    while covy[i]>=0:
    	tauy=tauy+(covy[i]+covy[i+1])/2;
    	i=i+1;

    # find degrees of freedom
    df=np.min([np.floor(mx/taux), np.floor(my/tauy)]);

    # find 95% confidence correlation
    b=2*1.96/np.sqrt(df-3);
    rho_sig95=(np.exp(b)-1)/(1+np.exp(b));

    return rho, df, rho_sig95

def lagcor(xx,yy,maxlag,fg):
    # remove means
    xp = xx-np.mean(xx);
    mx = len(xp)
    yp = yy-np.mean(yy);
    
    rho0 = np.corrcoef(xp,yp)[1,0] # rho[0,0] xp with itself, rho[1,0] xp with yp, rho[0,1] yp with xp, rho[1,1] yp with itself
    lags = np.arange(-(mx-1),mx)    
    rho = np.correlate(xp,yp,mode='full')  # get all lags (2*mx -1)
    rho = rho/rho[lags==0] * rho0

    # trim to lags
    rho = rho[mx-maxlag-1:mx+maxlag]
    lags = lags[mx-maxlag-1:mx+maxlag]
    if fg==1:
        plt.figure(11);plt.clf()
        plt.plot(lags,rho,'*-');plt.grid(True)
        plt.xlim(-maxlag-1, maxlag+1);plt.ylim(min(rho)-0.1, max(rho)+0.1)
        plt.text(-maxlag/2,min(rho),'xx leads yy')
        plt.text(maxlag/2-2,min(rho),'yy leads xx')
        
    return lags, rho

# detrend temporally (scipy.signal.detrend doesn't like nans or missing data)
def detrendNaN(t,y):
    y_detrend = np.nan*np.ones(y.shape)
    if not np.all(np.isnan(y)):
        m, b, r_val, p_val, std_err = stats.linregress(t[~np.isnan(y)],y[~np.isnan(y)])
        y_detrend = y - (m*t + b)
    
    return y_detrend

# interpolate from one 2d grid to another 2d grid where there are NaNs in the data (geo_in)
def interp2dNaN(lon_in,lat_in,geo_in,lon_out,lat_out):
    nan_map = np.zeros_like(geo_in)
    nan_map[np.isnan(geo_in)] = 1

    filled_z = geo_in.copy()
    filled_z[ np.isnan(geo_in) ] = 0

    f = interpolate.interp2d(lon_in, lat_in, filled_z, kind='linear')
    f_nan = interpolate.interp2d(lon_in, lat_in, nan_map, kind='linear')     

    geo_out = f(lon_out,lat_out)
    nan_new = f_nan(lon_out,lat_out)
    geo_out[nan_new>0.5] = np.nan
    
    return geo_out


