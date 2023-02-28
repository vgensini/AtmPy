#!/usr/bin/env python
import numpy as np
import AtmPy.srh as calc_srh
import AtmPy.SRH_plev as calc_srh_plev
__all__ = ['calc_helicity']

def calc_helicity(u,v,zagl):
    '''
    u: 3d array U wind(m/s) with dimension order level,lat,lon
    v: 3d array of V wind (m/s) with dimension order level,lat,lon
    zagl: 3d array of AGL heights (m) with dimension order level,lat,lon

    returns:
    srh500:
    srh01:
    srh03:
    '''
    depth = [500.,1000.,3000.]
    u2 = np.moveaxis(u,0,-1)
    v2 = np.moveaxis(v,0,-1)
    zagl2 = np.moveaxis(zagl,0,-1)

    helicity = calc_srh.calhel(u2,v2,zagl2,depth)
    srh500 = helicity[:,:,0]
    srh01 = helicity[:,:,1]
    srh03 = helicity[:,:,2]

    return srh500,srh01,srh03

def calc_helicity_plev(u,v,u10,v10,zagl):
    '''
    u: 3d array U wind(m/s) with dimension order level,lat,lon
    v: 3d array of V wind (m/s) with dimension order level,lat,lon
    zagl: 3d array of AGL heights (m) with dimension order level,lat,lon

    returns:
    srh500:
    srh01:
    srh03:
    '''
    depth = [500.,1000.,3000.]
    u2 = np.moveaxis(u,0,-1)
    v2 = np.moveaxis(v,0,-1)

    zagl2 = np.moveaxis(zagl,0,-1)

    helicity = calc_srh_plev.calhelplev(u2,v2,u10,v10,zagl2,depth)
    srh500 = helicity[:,:,0]
    srh01 = helicity[:,:,1]
    srh03 = helicity[:,:,2]

    return srh500,srh01,srh03

