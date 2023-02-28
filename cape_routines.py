import AtmPy
from AtmPy import cape_3D_plevs
from AtmPy import cape_4D_plevs
import numpy as np
import os

__all__ = ['cape_plev_3D','cape_plev_4D']


lookup_file = AtmPy.__path__.__dict__["_path"][0] + '/psadilookup.dat'
ter_follow = 0 # Pressure level data for all functions

def cape_plev_3D(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,parcel):
    '''
    Calculate CAPE on Pressure Levels

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
    parcel: string; MU, ML, SB, UL, ML3km, SB3km, : Most unstable, 100mb mixed-layer, surfaced based, most unstable in 0 to -20C layer

    ** order needs to be top down for vertical dimension **
    ** we will check for this below by sorting pressure  **

    output
    ------
    mucape: 2D array of most unstable parcel cape (J/kg)
    mucin: 2D array of most unstable parcel cin (J/kg)
    mulcl: 2D array of most unstable parcel lcl height (m)
    mulfc: 2D array of most unstable parcel lfc height (m)
    muel: 2D array of most unstable parcel equllibrium level height (m)
    mupght: 2D array of most unstable parcel geopotential height (m)
    mulclt: 2D array of most unstable parcel lcl temperature (K)
    muelt: 2D array of most unstable parcel equillibrium temperature (K)
    '''
    ter_follow = 0 #pressure level data
    if parcel == 'SB': # Surface based parcel
        cape_type = 1
    elif parcel == 'MU':# Most Unstable Parcel
        cape_type = 0 
    elif parcel == 'ML': # 100mb MLCAPE
        cape_type = 3
    elif parcel == 'UL': # MUCAPE in the 0 to -20C layer
        cape_type = 6
    elif parcel == 'SB3km': # SBCAPE in lowest 3km 
        cape_type = 4
    elif parcel == 'ML3km': # 100mb MLCAPE in lowest 3km 
        cape_type = 5

    #test to make sure vertical pressure is top down
    plevtest = prs_mb[:,0,0]
    sortedplev = np.sort(plevtest)
    if np.array_equal(plevtest,sortedplev) == False:
        print('Warning, your data are not sorted by pressure increasing (top down). Inverting...')
        prs_mb = prs_mb[::-1,:,:]
        tmp = tmp[::-1,:,:]
        mixr = mixr[::-1,:,:]
        hgt = hgt[::-1,:,:]

    cape,cin,lcl,lfc,el,pght,lclt,elt = cape_3D_plevs.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)
    return cape,cin,lcl,lfc,el,pght,lclt,elt

def cape_plev_4D(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,parcel):
    '''
    input
    -----
    prs_mb:4D array (time,plev,lat,lon) of pressure on vertical levels (mb)
    tmp: 4D array (time,plev,lat,lon) of temperature on vertical levels (K)
    mixr: 3D array (time,plev,lat,lon) of water vapor mixing ratio on vertical levels (kg/kg)
    hgt: 4D array (time,plev,lat,lon) of geopotential height on vertical levels (m)
    ter: 2D array (lat,lon) terrain height/orog (m)
    psfc_mb: 3D array (time,lat,lon) surface pressure (mb)
    sfc_t: 3D array (time,lat,lon) near-surface (e.g. 2m ) temperature (K)
    sfc_mixr: 3D array (time,lat,lon) near-surface (e.g. 2m ) mixing ratio (kg/kg) 

    ** order needs to be top down for vertical dimension **

    output
    ------
    mucape: 3D array of 0-3km most unstable parcel cape (J/kg)
    mucin: 3D array of 0-3km most unstable parcel cin (J/kg)
    mulcl: 3D array of 0-3km most unstable parcel lcl height (m)
    mulfc: 3D array of 0-3km most unstable parcel lfc height (m)
    '''
    ter_follow = 0 #pressure level data
    if parcel == 'SB': # Surface based parcel
        cape_type = 1
    elif parcel == 'MU':# Most Unstable Parcel
        cape_type = 0 
    elif parcel == 'ML': # 100mb MLCAPE
        cape_type = 3
    elif parcel == 'UL': # MUCAPE in the 0 to -20C layer
        cape_type = 6
    elif parcel == 'SB3km': # SBCAPE in lowest 3km 
        cape_type = 4
    elif parcel == 'ML3km': # 100mb MLCAPE in lowest 3km 
        cape_type = 5
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

    cape,cin,lcl,lfc,el,pght,lclt,elt = cape_4D_plevs.dcapecalc3d(prs_mb,tmp,mixr,hgt,ter,psfc_mb,sfc_t,sfc_mixr,cape_type,ter_follow,lookup_file)
    return cape,cin,lcl,lfc,el,pght,lclt,elt

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

