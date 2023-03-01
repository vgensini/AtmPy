import numpy

__all__ = ['srh','bunkers','bulk_wind_shear','mean_wind']

#storm relative helicity
def srh(z, u, v, dz, lower, upper, sm_u, sm_v, topo):

   #Calculates storm relative helicity provided input arrays of height, u, and v; as well as upper and
   #lower limits to the layer and storm motion vectors

   #Adapted from similar sharppy routine (https://github.com/sharppy/)
   #Copyright (c) 2015, Kelton Halbert, Greg Blumberg, Tim Supinie, and Patrick Marsh

   #Halbert, K. T., W. G. Blumberg, and P. T. Marsh, 2015: "SHARPpy: Fueling the Python Cult".
   #Preprints, 5th Symposium on Advances in Modeling and Analysis Using Python, Phoenix AZ.

   #Input:

   #z - 3d numpy array of height AGL (m)
   #u - 3d numpy array of the u component of the wind (m/s)
   #v - 3d numpy array of the v component of the wind (m/s)
   #dz - 3d numpy array [z,y,x] of vertical grid spacing for each grid point (m)
   #lower - float value of the lower boundary of the layer (m AGL)
   #upper - float value of the upper boundary of the layer (m AGL)
   #sm_u - float value of the u component of the storm motion (m/s)
   #sm_v - float value of the v component of the storm motion (m/s)
   #topo - terrain height (m)

   #Returns:

   #srh - 2d numpy array of storm-relative helicity for the layer

   #NOTE: Does NOT currently interpolate to lower and upper values (i.e. SRH is calculated between the nearest vertical 
   #layers to the lower and upper values, not the actual values)

#######################
    z = z - topo
    z[z<0.0]=0.0
    u = numpy.swapaxes(u, 0, 1)
    v = numpy.swapaxes(v, 0, 1)
    dz = numpy.swapaxes(dz, 0, 1)
    z = numpy.swapaxes(z,0,1)
    sr_u = u - sm_u
    sr_v = v - sm_v
    #print(sr_u.shape)
    du = sr_u[1:,:,:,:]-sr_u[:-1,:,:,:]          #du/dz
    dv = sr_v[1:,:,:,:]-sr_v[:-1,:,:,:]          #dv/dz
    layers = (sr_v[:-1,:,:,:]*(du/dz[:-1,:,:,:]) - sr_u[:-1,:,:,:]*(dv/dz[:-1,:,:,:])) * dz[:-1,:,:,:] 
    masked_layers = numpy.ma.masked_where((z[:-1,:,:,:]) > upper, (layers))
    masked_layers = numpy.ma.masked_where((z[:-1,:,:,:]) < lower, (masked_layers))
    srh = numpy.sum(masked_layers, axis=0)
    return srh

#bunkers storm motion
def bunkers(p, z, dz, u, v, topo):
    z = z - topo
    z[z<0.0]=0.0
    d = 7.5 #Empirically-derived deviation value from Bunkers et al. (2014; JOM) 
    upper = 6000. #Upper-limit to storm motion layer (m AGL)
    lower = 0. #Lower-limit to storm motion layer (m AGL)
    mean_u, mean_v = mean_wind(p, z, dz, u, v, lower, upper)
    shear_u, shear_v, shear = bulk_wind_shear(z, u, v, lower, upper)
    modifier = d / numpy.sqrt(shear_u**2 + shear_v**2) 
    bunk_r_u = mean_u + (modifier * shear_v)
    bunk_r_v = mean_v - (modifier * shear_u)
    bunk_l_u = mean_u - (modifier * shear_v)
    bunk_l_v = mean_v + (modifier * shear_u)
    return bunk_r_u, bunk_r_v, bunk_l_u, bunk_l_v, shear

#bulk wind shear
def bulk_wind_shear(z, u, v, lower, upper):
    u = numpy.ma.masked_where(z > upper, (u))
    u = numpy.ma.masked_where(z < lower, (u))
    v = numpy.ma.masked_where(z > upper, (v))
    v = numpy.ma.masked_where(z < lower, (v))
    
    lower_indices, upper_indices = numpy.ma.notmasked_edges(u,axis=1)
    u_upper = u[upper_indices].reshape(u.shape[0],u.shape[2],u.shape[3])
    u_lower = u[lower_indices].reshape(u.shape[0],u.shape[2],u.shape[3])
    v_upper = v[upper_indices].reshape(v.shape[0],v.shape[2],v.shape[3])
    v_lower = v[lower_indices].reshape(v.shape[0],v.shape[2],v.shape[3])

    shear_u = u_upper - u_lower
    shear_v = v_upper - v_lower
    shear = numpy.sqrt(shear_u**2 + shear_v**2)
    return shear_u, shear_v, shear

#pressure weighted mean wind
def mean_wind(p, z, dz, u, v, lower, upper):
    p = numpy.ma.masked_where(z > upper, (p))
    p = numpy.ma.masked_where(z < lower, (p))
    u = numpy.ma.masked_where(z > upper, (u))
    u = numpy.ma.masked_where(z < lower, (u))
    v = numpy.ma.masked_where(z > upper, (v))
    v = numpy.ma.masked_where(z < lower, (v))
    dz = numpy.ma.masked_where(z > upper, (dz))
    dz = numpy.ma.masked_where(z < lower, (dz))
    mean_u = numpy.ma.average(u, axis=1, weights=(dz*p))
    mean_v = numpy.ma.average(v, axis=1, weights=(dz*p)) 
    return mean_u, mean_v 