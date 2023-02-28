#import AtmPy.cape_plev as capecalc_plev
import AtmPy.mod_rip_cape_plev_4D as capecalc_plev_4D
#import AtmPy.dcape as dcapecalc_plev
#import AtmPy.dcape_plev_4D as dcapecalc_plev_4D
import numpy as np
import os

__all__ = ['mlcape_plev','mlcape3km_plev','mlcape3km_plev_4D','sfcape_plev','sfcape3km_plev','sfcape3km_plev_4D','mucape_plev','mucape_plev_4D','mlcape_plev','mlcape_plev_4D','sfcape_plev_4D','dcape_plev','dcape_plev_4D']

lookup_file = '/home/data/gefs_reforecasts/psadilookup.dat'


def mucape_plev(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    sfc_t: near-surface (e.g. 2m ) temperature (K)
    sfc_mixr: near-surface (e.g. 2m ) mixing ratio (kg/kg) 
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    mucape: 2D array of 0-3km most unstable parcel cape (J/kg)
    mucin: 2D array of 0-3km most unstable parcel cin (J/kg)
    mulcl: 2D array of 0-3km most unstable parcel lcl height (m)
    mulfc: 2D array of 0-3km most unstable parcel lfc height (m)
    '''
    ter_follow = 0 #pressure level data
    cape_type = 0
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[::-1,:,:]
        tmp = tmp[::-1,:,:]
        mixr = mixr[::-1,:,:]
        hgt = hgt[::-1,:,:]

    mucape,mucin,mulcl,mulfc,muel,mupght,mulclt,muelt = capecalc_plev.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)
    return mucape,mucin,mulcl,mulfc,muel,mupght,mulclt,muelt


def mucape_plev_4D(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (time,plev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (time,plev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (time,plev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (time,plev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 3D array (time,lat,lon) surface pressure (mb)
    sfc_t: 3D array (time,lat,lon) near-surface (e.g. 2m ) temperature (K)
    sfc_mixr: 3D array (time,lat,lon) near-surface (e.g. 2m ) mixing ratio (kg/kg) 
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    mucape: 3D array of 0-3km most unstable parcel cape (J/kg)
    mucin: 3D array of 0-3km most unstable parcel cin (J/kg)
    mulcl: 3D array of 0-3km most unstable parcel lcl height (m)
    mulfc: 3D array of 0-3km most unstable parcel lfc height (m)
    '''
    ter_follow = 0 #pressure level data
    cape_type = 0
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[0,:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[:,::-1,:,:]
        tmp = tmp[:,::-1,:,:]
        mixr = mixr[:,::-1,:,:]
        hgt = hgt[:,::-1,:,:]
    prs_mb = np.swapaxes(prs_mb,0,1)
    tmp = np.swapaxes(tmp,0,1)
    mixr = np.swapaxes(mixr,0,1)
    hgt = np.swapaxes(hgt,0,1)

    mucape,mucin,mulcl,mulfc,muel,mupght,mulclt,muelt = capecalc_plev_4D.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)
    return mucape,mucin,mulcl,mulfc,muel,mupght,mulclt,muelt





def mucape_0to20_plev(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    sfc_t: near-surface (e.g. 2m ) temperature (K)
    sfc_mixr: near-surface (e.g. 2m ) mixing ratio (kg/kg) 
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    mucape: 2D array of 0-3km most unstable parcel cape (J/kg)
    mucin: 2D array of 0-3km most unstable parcel cin (J/kg)
    mulcl: 2D array of 0-3km most unstable parcel lcl height (m)
    mulfc: 2D array of 0-3km most unstable parcel lfc height (m)
    '''
    ter_follow = 0 #pressure level data
    cape_type = 6
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[::-1,:,:]
        tmp = tmp[::-1,:,:]
        mixr = mixr[::-1,:,:]
        hgt = hgt[::-1,:,:]

    mucape,mucin,mulcl,mulfc,muel,mupght,mulclt,muelt = capecalc_plev.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)

    return mucape,mucin,mulcl,mulfc,muel,mupght

def mucape_0to20_plev_4D(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    sfc_t: near-surface (e.g. 2m ) temperature (K)
    sfc_mixr: near-surface (e.g. 2m ) mixing ratio (kg/kg) 
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    mucape: 2D array of 0-3km most unstable parcel cape (J/kg)
    mucin: 2D array of 0-3km most unstable parcel cin (J/kg)
    mulcl: 2D array of 0-3km most unstable parcel lcl height (m)
    mulfc: 2D array of 0-3km most unstable parcel lfc height (m)
    '''
    ter_follow = 0 #pressure level data
    cape_type = 6
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[0,:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[:,::-1,:,:]
        tmp = tmp[:,::-1,:,:]
        mixr = mixr[:,::-1,:,:]
        hgt = hgt[:,::-1,:,:]
    prs_mb = np.swapaxes(prs_mb,0,1)
    tmp = np.swapaxes(tmp,0,1)
    mixr = np.swapaxes(mixr,0,1)
    hgt = np.swapaxes(hgt,0,1)


    mucape,mucin,mulcl,mulfc,muel,mupght,mulclt,muelt = capecalc_plev_4D.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)

    return mucape,mucin,mulcl,mulfc,muel,mupght



def mlcape_plev(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type=3,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    *UPP version of MLCAPE (lifts from mixed-layer level)

    mlcape: 2D array of 100mb mixed layer cape (J/kg)
    mlcin: 2D array of 100mb mixed layer cin (J/kg)
    mllcl: 2D array of 100mb mixed layer lcl height (m)
    mllfc: 2D array of 100mb mixed layer lfc height (m)
    '''
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[::-1,:,:]
        tmp = tmp[::-1,:,:]
        mixr = mixr[::-1,:,:]
        hgt = hgt[::-1,:,:]
    ter_follow = 0 #pressure level data
    #cape_type = 3 #raise from sfc=3 ,raise from midlayer=2
    mlcape,mlcin,mllcl,mllfc,mlel,mlpght,mllclt,mlelt = capecalc_plev.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)

    return mlcape,mlcin,mllcl,mllfc,mlel,mlpght,mllclt,mlelt

def mlcape_plev_4D(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type=3,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    *UPP version of MLCAPE (lifts from mixed-layer level)

    mlcape: 2D array of 100mb mixed layer cape (J/kg)
    mlcin: 2D array of 100mb mixed layer cin (J/kg)
    mllcl: 2D array of 100mb mixed layer lcl height (m)
    mllfc: 2D array of 100mb mixed layer lfc height (m)
    '''
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[0,:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[:,::-1,:,:]
        tmp = tmp[:,::-1,:,:]
        mixr = mixr[:,::-1,:,:]
        hgt = hgt[:,::-1,:,:]

    prs_mb = np.swapaxes(prs_mb,0,1)
    tmp = np.swapaxes(tmp,0,1)
    mixr = np.swapaxes(mixr,0,1)
    hgt = np.swapaxes(hgt,0,1)


    ter_follow = 0 #pressure level data
    #cape_type = 3 #raise from sfc,2 raise from midlyr
    mlcape,mlcin,mllcl,mllfc,mlel,mlpght,mllclt,mlelt = capecalc_plev_4D.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)

    return mlcape,mlcin,mllcl,mllfc,mlel,mlpght,mllclt,mlelt


def mlcape3km_plev(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    *UPP version of MLCAPE (lifts from mixed-layer level)

    mlcape: 2D array of 100mb mixed layer cape (J/kg)
    mlcin: 2D array of 100mb mixed layer cin (J/kg)
    mllcl: 2D array of 100mb mixed layer lcl height (m)
    mllfc: 2D array of 100mb mixed layer lfc height (m)
    '''
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[::-1,:,:]
        tmp = tmp[::-1,:,:]
        mixr = mixr[::-1,:,:]
        hgt = hgt[::-1,:,:]
    ter_follow = 0 #pressure level data
    cape_type = 5
    mlcape,mlcin,mllcl,mllfc,mlel,mlpght,mllclt,mlelt = capecalc_plev.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)
    return mlcape,mlcin,mllcl,mllfc,mlel,mlpght

def mlcape3km_plev_4D(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    *UPP version of MLCAPE (lifts from mixed-layer level)

    mlcape: 2D array of 100mb mixed layer cape (J/kg)
    mlcin: 2D array of 100mb mixed layer cin (J/kg)
    mllcl: 2D array of 100mb mixed layer lcl height (m)
    mllfc: 2D array of 100mb mixed layer lfc height (m)
    '''
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[0,:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[:,::-1,:,:]
        tmp = tmp[:,::-1,:,:]
        mixr = mixr[:,::-1,:,:]
        hgt = hgt[:,::-1,:,:]

    prs_mb = np.swapaxes(prs_mb,0,1)
    tmp = np.swapaxes(tmp,0,1)
    mixr = np.swapaxes(mixr,0,1)
    hgt = np.swapaxes(hgt,0,1)

    ter_follow = 0 #pressure level data
    cape_type = 5
    mlcape,mlcin,mllcl,mllfc,mlel,mlpght,mllclt,mlelt = capecalc_plev_4D.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)
    return mlcape,mlcin,mllcl,mllfc,mlel,mlpght



def sfcape_plev(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    sfc_t: near-surface (e.g. 2m ) temperature (K)
    sfc_mixr: near-surface (e.g. 2m ) mixing ratio (kg/kg) 
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    sfccape: 2D array of surface based cape (J/kg)
    sfccin: 2D array of surface based cin (J/kg)
    sfclcl: 2D array of surface based lcl height (m)
    sfclfc: 2D array of surface based lfc height (m)
    '''
    ter_follow = 0 #pressure level data
    cape_type = 1
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[::-1,:,:]
        tmp = tmp[::-1,:,:]
        mixr = mixr[::-1,:,:]
        hgt = hgt[::-1,:,:]

    sfccape,sfccin,sfclcl,sfclfc,sfcel,sfcpght,sfclclt,sfcelt = capecalc_plev.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)
    return sfccape,sfccin,sfclcl,sfclfc,sfcel,sfcpght,sfclclt,sfcelt

def sfcape_plev_4D(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    sfc_t: near-surface (e.g. 2m ) temperature (K)
    sfc_mixr: near-surface (e.g. 2m ) mixing ratio (kg/kg) 
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    sfccape: 2D array of surface based cape (J/kg)
    sfccin: 2D array of surface based cin (J/kg)
    sfclcl: 2D array of surface based lcl height (m)
    sfclfc: 2D array of surface based lfc height (m)
    '''
    ter_follow = 0 #pressure level data
    cape_type = 1
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[0,:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[:,::-1,:,:]
        tmp = tmp[:,::-1,:,:]
        mixr = mixr[:,::-1,:,:]
        hgt = hgt[:,::-1,:,:]

    prs_mb = np.swapaxes(prs_mb,0,1)
    tmp = np.swapaxes(tmp,0,1)
    mixr = np.swapaxes(mixr,0,1)
    hgt = np.swapaxes(hgt,0,1)

    sfccape,sfccin,sfclcl,sfclfc,sfcel,sfcpght,sfclclt,sfcelt = capecalc_plev_4D.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)
    return sfccape,sfccin,sfclcl,sfclfc,sfcel,sfcpght,sfclclt,sfcelt


def sfcape3km_plev(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    sfc_t: near-surface (e.g. 2m ) temperature (K)
    sfc_mixr: near-surface (e.g. 2m ) mixing ratio (kg/kg) 
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    sfccape: 2D array of 0-3km sfc cape (J/kg)
    sfccin: 2D array of 0-3km most sfc cin (J/kg)
    sfclcl: 2D array of 0-3km most sfc lcl height (m)
    sfclfc: 2D array of 0-3km most sfc lfc height (m)
    '''
    ter_follow = 0 #pressure level data
    cape_type = 4
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[::-1,:,:]
        tmp = tmp[::-1,:,:]
        mixr = mixr[::-1,:,:]
        hgt = hgt[::-1,:,:]

    sfccape,sfccin,sfclcl,sfclfc,sfcel,sfcpght,sfclclt,sfcelt = capecalc_plev.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)

    return sfccape,sfccin,sfclcl,sfclfc,sfcel,sfcpght,sfclclt,sfcelt


def sfcape3km_plev_4D(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    sfc_t: near-surface (e.g. 2m ) temperature (K)
    sfc_mixr: near-surface (e.g. 2m ) mixing ratio (kg/kg) 
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    sfccape: 2D array of 0-3km sfc cape (J/kg)
    sfccin: 2D array of 0-3km most sfc cin (J/kg)
    sfclcl: 2D array of 0-3km most sfc lcl height (m)
    sfclfc: 2D array of 0-3km most sfc lfc height (m)
    '''
    ter_follow = 0 #pressure level data
    cape_type = 4
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[0,:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[:,::-1,:,:]
        tmp = tmp[:,::-1,:,:]
        mixr = mixr[:,::-1,:,:]
        hgt = hgt[:,::-1,:,:]

    prs_mb = np.swapaxes(prs_mb,0,1)
    tmp = np.swapaxes(tmp,0,1)
    mixr = np.swapaxes(mixr,0,1)
    hgt = np.swapaxes(hgt,0,1)


    sfccape,sfccin,sfclcl,sfclfc,sfcel,sfcpght,sfclclt,sfcelt = capecalc_plev_4D.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)

    return sfccape,sfccin,sfclcl,sfclfc,sfcel,sfcpght,sfclclt,sfcelt


def dcape_plev(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    sfc_t: near-surface (e.g. 2m ) temperature (K)
    sfc_mixr: near-surface (e.g. 2m ) mixing ratio (kg/kg) 
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    dcape: 2D array of 0-3km most unstable parcel cape (J/kg)
    dcin: 2D array of 0-3km most unstable parcel cin (J/kg)
    dlcl: 2D array of 0-3km most unstable parcel lcl height (m)
    dlfc: 2D array of 0-3km most unstable parcel lfc height (m)
    '''
    ter_follow = 0 #pressure level data
    cape_type = 0
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[::-1,:,:]
        tmp = tmp[::-1,:,:]
        mixr = mixr[::-1,:,:]
        hgt = hgt[::-1,:,:]

    dcape,dcin,dlcl,dlfc,delh,dpght,dlclt,delt = dcapecalc_plev.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)
    return dcape,dcin

def dcape_plev_4D(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,lookup_file=lookup_file):
    '''
    input
    -----
    prs_mb:3D array (vlev,lat,lon) of pressure on vertical levels (mb)
    tmp: 3D array (vlev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (vlev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 3D array (vlev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 2D array (lat,lon) surface pressure (mb)
    sfc_t: near-surface (e.g. 2m ) temperature (K)
    sfc_mixr: near-surface (e.g. 2m ) mixing ratio (kg/kg) 
    ter_follow: scalar; 0: pressure level data, 1: terrain following data

    ** order needs to be top down for vertical dimension **

    output
    ------
    dcape: 2D array of 0-3km most unstable parcel cape (J/kg)
    dcin: 2D array of 0-3km most unstable parcel cin (J/kg)
    dlcl: 2D array of 0-3km most unstable parcel lcl height (m)
    dlfc: 2D array of 0-3km most unstable parcel lfc height (m)
    '''
    ter_follow = 0 #pressure level data
    cape_type = 0
    #test to make sure vertical pressure is top down
    plevtest = prs_mb[0,:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        prs_mb = prs_mb[:,::-1,:,:]
        tmp = tmp[:,::-1,:,:]
        mixr = mixr[:,::-1,:,:]
        hgt = hgt[:,::-1,:,:]


    dcape,dcin,dlcl,dlfc,delh,dpght,dlclt,delt = dcapecalc_plev_4D.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)
    return dcape,dcin

