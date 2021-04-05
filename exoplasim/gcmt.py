import numpy as np
import netCDF4 as nc

class _Dataset:
    def __init__(self,filename):
        if filename[-3:]==".nc":
            self.body = nc.Dataset(filename,"r")
            self.variables=self.body.variables
        elif filename[-4:]==".npy":
            self.body = np.load(filename)
            self.variables = self.body.item()
        else:
            raise DatafileError("Unknown dataset format")
        
    def close(self):
        try:
            self.body.close()
        except:
            pass

class DimensionError(Exception):
    pass

class UnitError(Exception):
    pass

class DatafileError(Exception):
    pass

def parse(file,variable,lat=None,lon=None):
    """Retrieve a variable from a NetCDF file

    Parameters
    ----------
    file : str
        Path to a NetCDF file
    variable : str
        Name of the variable to extract
    lat,lon : str, optional
        If the latitude and longitude arrays have non-standard names, specify them here.

    Returns
    -------
    numpy.ndarray
        Requested output field

    """
    ncd=_Dataset(file)
    variable = ncd.variables[variable][:]
    
    if lat:
        lt = ncd.variables[lat][:]
    elif "lat" in ncd.variables:
        lt = ncd.variables['lat'][:]
    elif "lt" in ncd.variables:
        lt = ncd.variables["lt"][:]
    elif "lats" in ncd.variables:
        lt = ncd.variables['lats'][:]
    elif "latitude" in ncd.variables:
        lt = ncd.variables['latitude'][:]
    elif "latitudes" in ncd.variables:
        lt = ncd.variables['latitudes'][:]
    else:
        raise DatafileError("Unknown datafile format; unsure how to extract latitude")
    
    if lon:
        ln = ncd.variables[lon][:]
    elif "lat" in ncd.variables:
        ln = ncd.variables['lon'][:]
    elif "lt" in ncd.variables:
        ln = ncd.variables["ln"][:]
    elif "lats" in ncd.variables:
        ln = ncd.variables['lons'][:]
    elif "latitude" in ncd.variables:
        ln = ncd.variables['longitude'][:]
    elif "latitudes" in ncd.variables:
        ln = ncd.variables['longitudes'][:]
    else:
        raise DatafileError("Unknown datafile format; unsure how to extract longitude")
    
    ncd.close()
    
    return ln,lt,variable




def make2d(variable,lat=None,lon=None,time=None,lev=None,ignoreNaNs=True,radius=6.371e6,
           latitudes=None,longitudes=None):
    """Compress a variable in two dimensions by slicing or averaging.

    Parameters
    ----------
    variable : numpy.ndarray
        The variable to operate on
    lat,lon,lev : int, str, optional
        Either an index on which to slice, or either of "sum" or "mean", indicating what
        should be done along that axis.
    time : int, optional
        The time index on which to slice. If unspecified, a time average will be returned.
    ignoreNaNs : bool, optional
        If set, will use NaN-safe numpy operators.
    radius : float, optional
        Planet radius in meters (only used for summation)
    latitudes: numpy.ndarray, optional
        Latitude array--required if lat is "mean", or if either lat or lon is "sum"
    longitudes: numpy.ndarray, optional
        Longitude array--required if lon is "mean" or if either lat or lon is "sum"
        
    Returns
    -------
    numpy.ndarray
        A 2-D array
    """
    if ignoreNaNs:
        sumop = np.nansum
        meanop = np.nanmean
    else:
        sumop = np.sum
        meanop = np.mean
    if len(variable.shape)==2:
        return variable
    if time:
        try:
            variable=variable[time,:]
        except:
            raise UnitError("You have probably passed a float time to a variable with no "+
                            "information about what that means. You should pass an integer "+
                            "time index instead")
    elif time==None and len(variable.shape)>2:
        variable=meanop(variable,axis=0)
    elif time==0:
        variable=variable[time,:]
    if len(variable.shape)>2:
        if lev!=None:
            if type(lev)==int:
                variable=variable[lev,:]
            elif lev=="sum":
                variable=sumop(variable,axis=0)
            elif lev=="mean":
                variable=meanop(variable,axis=0)
            else:
                raise UnitError("Unknown level specification")
        elif lat!=None and lon==None:
            if type(lat)==int:
                variable=variable[:,lat,:]
            elif lat=="sum":
                variable=latsum(variable,latitudes,dlon=longitudes[1]-longitudes[0],
                                radius=radius)
            elif lat=="mean":
                variable=latmean(variable,latitudes)
            else:
                raise UnitError("Unknown latitude specification")
        elif lon!=None and lat==None:
            if type(lon)==int:
                variable=variable[:,:,lon]
            elif lon=="sum":
                newvar = np.zeros(variable.shape[:-1])
                gradlat = np.gradient(np.sin(latitudes*np.pi/180.))
                for lt in range(len(latitudes)):
                    newvar[:,lt]=lonsum(variable,longitudes,dsinlat=gradlat[lt],radius=radius)
                variable = newvar
            elif lon=="mean":
                variable=lonmean(variable,longitudes)
            else:
                raise UnitError("Unknown longitude specification")
        else:
            raise DimensionError("Inappropriate or insufficient dimensional constraints")
    
    return variable
    

def spatialmath(variable,lat=None,lon=None,file=None,mean=True,time=None,
               ignoreNaNs=True,lev=None,radius=6.371e6):
    """Compute spatial means or sums of data

    Parameters
    ----------
    variable : str, numpy.ndarray
        The variable to operate on. Can either be a data array, or the name of a variable. If the latter, file must be specified.
    lat,lon : numpy.ndarray, optional
        Latitude and longitude arrays. If file is provided and lat and lon are not, they will be
        extracted from the file.
    file : str, optional
        Path to a NetCDF output file to open and extract data from.
    mean : bool, optional
        If True, compute a global mean. If False, compute a global sum.
    time : int, optional
        The time index on which to slice. If unspecified, a time average will be returned.
    ignoreNaNs : bool, optional
        If True, use NaN-safe numpy operators.
    lev : int, optional
        If set, slice a 3D spatial array at the specified level.
    radius : float, optional
        Radius of the planet in meters. Only used if mean=False.
        
    Returns
    -------
    float

    """
    
    if ignoreNaNs:
        sumop = np.nansum
        meanop = np.nanmean
    else:
        sumop = np.sum
        meanop = np.mean
        
    if file:
        ln,lt,variable = parse(file,variable,lat=lat,lon=lon)
        
    else:
        if type(lat)==type(None) or type(lon)==type(None):
            raise DimensionError("Need to provide latitude and longitude data")
        ln=lon
        lt=lat
    variable = make2d(variable,time=time,lev=lev,ignoreNaNs=ignoreNaNs)
    
    lt1 = np.zeros(len(lt)+1)
    lt1[0] = 90
    for n in range(0,len(lt)-1):
        lt1[n+1] = 0.5*(lt[n]+lt[n+1])
    lt1[-1] = -90
    dln = np.diff(ln)[0]
    ln1 = np.zeros(len(ln)+1)
    ln1[0] = -dln
    for n in range(0,len(ln)-1):
        ln1[n+1] = 0.5*(ln[n]+ln[n+1])
    ln1[-1] = 360.0-dln
    
    lt1*=np.pi/180.0
    ln1*=np.pi/180.0
    
    darea = np.zeros((len(lt),len(ln)))
    for jlat in range(0,len(lt)):
        for jlon in range(0,len(ln)):
            dln = ln1[jlon+1]-ln1[jlon]
            darea[jlat,jlon] = abs(np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*abs(dln)
    
    svar = variable*darea
    if mean:
        outvar = sumop(svar)/sumop(darea)
    else:
        outvar = sumop(svar) * radius**2
    
    return outvar

def latmean(variable,latitudes):
    """Compute meriodional mean (i.e. the variable that changes is latitude).
    
    Compute the area-weighted mean of a latitude array :math:`x`\ , such that:

    .. math::

        \\bar{x} = \\frac{\sum_{i=1}^N |\\sin(\\phi_{i-1/2})-\\sin(\\phi_{i+1/2})|x_i}{\sum_{i=1}^N |\\sin(\\phi_{i-1/2})-\\sin(\\phi_{i+1/2})|}
    
    Parameters
    ----------
    variable : numpy.ndarray
        Array to be averaged. Assumption is that if 2D, lat is the first dimension, if 3D, the second dimension, and if 4D. the 3rd dimension.
    latitudes : array-like
        Array or list of latitudes
        
    Returns
    -------
    scalar or numpy.ndarray
        Depending on the dimensionality of the input array, output may have 0, 1, or 2 dimensions.
    """
        
    lt1 = np.zeros(len(latitudes)+1)
    lt1[0] = 90
    for n in range(0,len(latitudes)-1):
        lt1[n+1] = 0.5*(latitudes[n]+latitudes[n+1])
    lt1[-1] = -90
    
    lt1*=np.pi/180.0
    
    darea = np.zeros(latitudes.shape)
    for jlat in range(0,len(latitudes)):
        darea[jlat] = np.sin(lt1[jlat])-np.sin(lt1[jlat+1])
    
    if len(variable.shape)==1:
        return np.nansum(variable*darea)/np.nansum(darea)
    elif len(variable.shape)==2:
        return np.nansum(variable*darea[:,np.newaxis],axis=0)/np.nansum(darea[:,np.newaxis]*np.ones(variable.shape),axis=0)
    elif len(variable.shape)==3:
        return np.nansum(variable*darea[np.newaxis,:,np.newaxis],axis=1)/np.nansum(darea[np.newaxis,:,np.newaxis]*np.ones(variable.shape),axis=1)
    elif len(variable.shape)==4:
        return np.nansum(variable*darea[np.newaxis,np.newaxis,:,np.newaxis],axis=2)/np.nansum(darea[np.newaxis,np.newaxis,:,np.newaxis]*np.ones(variable.shape),axis=2)
    else:
        raise DimensionError("Variable must have 4 or fewer dimensions. Latitude should be the second-from the right-most dimension if there are 2 or more dimensions.")
        
def latsum(variable,latitudes,dlon=360.0,radius=6.371e6):
    """Compute meriodional sum (i.e. the variable that changes is latitude).
    
    Compute the area-weighted sum of a latitude array :math:`x` given a longitude span :math:`\\Delta\\theta` and planet radius :math:`R`\ , such that:

    .. math::

        X = \sum_{i=1}^N |\\sin(\\phi_{i-1/2})-\\sin(\\phi_{i+1/2})|\\Delta\\theta R^2x_i
    
    Parameters
    ----------
    variable : numpy.ndarray
        Array to be summed. Assumption is that if 2D, lat is the first dimension, if 3D, the second dimension, and if 4D. the 3rd dimension.
    latitudes : array-like
        Array or list of latitudes
    dlon : float, optional
        Longitude span in degrees.
    radius : float, optional
        Planet radius in meters.
        
    Returns
    -------
    scalar or numpy.ndarray
        Depending on the dimensionality of the input array, output may have 0, 1, or 2 dimensions.
    """
        
    lt1 = np.zeros(len(latitudes)+1)
    lt1[0] = 90
    for n in range(0,len(latitudes)-1):
        lt1[n+1] = 0.5*(latitudes[n]+latitudes[n+1])
    lt1[-1] = -90
    
    lt1*=np.pi/180.0
    dlon *= np.pi/180.0
    
    darea = np.zeros(latitudes.shape)
    for jlat in range(0,len(latitudes)):
        darea[jlat] = abs(np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*abs(dlon)*radius**2
    
    if len(variable.shape)==1:
        return np.nansum(variable*darea)
    elif len(variable.shape)==2:
        return np.nansum(variable*darea[:,np.newaxis],axis=0)
    elif len(variable.shape)==3:
        return np.nansum(variable*darea[np.newaxis,:,np.newaxis],axis=1)
    elif len(variable.shape)==4:
        return np.nansum(variable*darea[np.newaxis,np.newaxis,:,np.newaxis],axis=2)
    else:
        raise DimensionError("Variable must have 4 or fewer dimensions. Latitude should be the second-from the right-most dimension if there are 2 or more dimensions.")


def lonmean(variable,longitudes):
    """Compute zonal mean (i.e. the variable that changes is longitude).
    
    Compute the area-weighted mean of a longitude array :math:`x`\ , such that:

    .. math::

        \\bar{x} = \\frac{\sum_{i=1}^N |\\theta_{i-1/2}-\\theta_{i+1/2}|x_i}{\sum_{i=1}^N |\\theta_{i-1/2}-\\theta_{i+1/2}|}
    
    Parameters
    ----------
    variable : numpy.ndarray
        Array to be summed. Assumption is that longitude is always the last dimension.
        
    Returns
    -------
    scalar or numpy.ndarray
        Depending on the dimensionality of the input array, output may be a scalar or have N-1 dimensions.
    """
    
    dlon = np.gradient(longitudes)
    sumlon = np.nansum(dlon)
    dlon = np.broadcast_to(dlon,variable.shape)
    
    return np.nansum(variable*dlon,axis=-1)/sumlon
        
def lonsum(variable,longitudes,dsinlat=2.0,radius=6.371e6):
    """Compute zonal sum (i.e. the variable that changes is longitude).
    
    Compute the area-weighted sum of a longitude array :math:`x` given a latitude span :math:`\\Delta\\sin\\phi` and planet radius :math:`R`\ , such that:

    .. math::

        X = \sum_{i=1}^N |\\theta_{i-1/2}-\\theta_{i+1/2}|\\Delta\\sin\\phi R^2x_i
    
    Parameters
    ----------
    variable : numpy.ndarray
        Array to be summed. Assumption is that longitude is always the last dimension.
    longitudes : array-like
        Array or list of longitudes
    dsinlat : float, optional
        The sine-latitude span for the longitude span considered. The default is 2, corresponding to -90 degrees to 90 degrees.
    radius : float, optional
        Planet radius in meters.
        
    Returns
    -------
    scalar or numpy.ndarray
        Depending on the dimensionality of the input array, output may have 0, 1, or 2 dimensions.
    """
        
    dlon = np.gradient(longitudes)*np.pi/180.0
    
    darea = np.zeros(longitudes.shape)
    for jlon in range(0,len(longitudes)):
        darea[jlon] = abs(dsinlat)*abs(dlon[jlon])*radius**2
    
    darea = np.broadcast_to(darea,variable.shape)
    
    return np.nansum(variable*darea,axis=-1)
    
def cspatialmath(variable,lat=None,lon=None,file=None,mean=True,time=None,
               ignoreNaNs=True,lev=None,radius=6.371e6,poles=False):
    """Compute spatial means or sums of data, but optionally don't go all the way to the poles.

    Sometimes, saying that the latitudes covered go all the way to :math:`\pm90^\circ` results in
    errors, and accurate accounting requires excluding the poles themselves. This function
    is identical to spatialmath, except that it provides that option.

    Parameters
    ----------
    variable : str, numpy.ndarray
        The variable to operate on. Can either be a data array, or the name of a variable. If the latter, file must be specified.
    lat,lon : numpy.ndarray, optional
        Latitude and longitude arrays. If file is provided and lat and lon are not, they will be
        extracted from the file.
    file : str, optional
        Path to a NetCDF output file to open and extract data from.
    mean : bool, optional
        If True, compute a global mean. If False, compute a global sum.
    time : int, optional
        The time index on which to slice. If unspecified, a time average will be returned.
    ignoreNaNs : bool, optional
        If True, use NaN-safe numpy operators.
    lev : int, optional
        If set, slice a 3D spatial array at the specified level.
    radius : float, optional
        Radius of the planet in meters. Only used if mean=False.
    poles : bool, optional
        If False (default), exclude the poles.
        
    Returns
    -------
    float

    """
    
    if ignoreNaNs:
        sumop = np.nansum
        meanop = np.nanmean
    else:
        sumop = np.sum
        meanop = np.mean
        
    if file:
        ln,lt,variable = parse(file,variable,lat=lat,lon=lon)
        
    else:
        if type(lat)==type(None) or type(lon)==type(None):
            raise DimensionError("Need to provide latitude and longitude data")
        ln=lon
        lt=lat
    variable = make2d(variable,time=time,lev=lev,ignoreNaNs=ignoreNaNs)
    
    lt1 = np.zeros(len(lt)+1)
    dlt1 = abs(np.diff(lt)[0])
    dlt2 = abs(np.diff(lt)[-1])
    lt1[0] = lt[0]+0.5*dlt1
    for n in range(0,len(lt)-1):
        lt1[n+1] = 0.5*(lt[n]+lt[n+1])
    lt1[-1] = lt[-1]-0.5*dlt2
    if poles:
        lt1[0] = 90.0
        lt1[-1] = -90.0
    
    dln = np.diff(ln)[0]
    ln1 = np.zeros(len(ln)+1)
    ln1[0] = ln[0]-dln*0.5
    for n in range(0,len(ln)-1):
        ln1[n+1] = 0.5*(ln[n]+ln[n+1])
    ln1[-1] = ln[-1]+dln*0.5
    
    lt1*=np.pi/180.0
    ln1*=np.pi/180.0
    
    darea = np.zeros((len(lt),len(ln)))
    for jlat in range(0,len(lt)):
        for jlon in range(0,len(ln)):
            dln = ln1[jlon+1]-ln1[jlon]
            darea[jlat,jlon] = abs(np.sin(lt1[jlat])-np.sin(lt1[jlat+1]))*abs(dln)
    
    svar = variable*darea
    if mean:
        outvar = sumop(svar)/sumop(darea)
    else:
        outvar = sumop(svar) * radius**2
    
    return outvar

def wrap2d(var):
    '''Add one element to the longitude axis to allow for wrapping'''
    newvar = np.zeros(np.array(var.shape)+np.array((0,1)))
    newvar[:,:-1] = var[:,:]
    newvar[:,-1] = var[:,0]
    return newvar

    
def streamfxn(file,time=None):
    '''Deprecated. Passes args to eqstream().'''
    return eqstream(file,time=time)
    
def eqstream(file,time=None):
    '''Return the streamfunction

    Parameters
    ----------
    file : str
        Path to an ExoPlaSim NetCDF output file.
        
    Returns
    -------
    numpy.ndarray
        The streamfunction for the given file.
    '''
    ln,lt,levs=parse(file,"lev")
    plevs = levs*spatialmath("ps",file=file,time=time)
    ln,lt,va = parse(file,"va")
    va = make2d(va,lon="mean",time=time)
    strf = np.zeros(va.shape)
    pref = 2*np.pi*6.371e6*np.cos(lt*np.pi/180.0)/9.81
    ps = np.array([0.0,]+list(plevs))
    vas = np.zeros(np.array(va.shape)+np.array((1,0)))
    vas[1:,:] = va[:,:]
    for k in range(0,len(plevs)):
        strf[k,:] = pref[:]*np.trapz(vas[0:k+1,:],x=ps[0:k+1],axis=0)
    return strf


def adist(lon1,lat1,lon2,lat2):
    '''Return angular distance(s) in degrees between two points (or sets of points) on a sphere.
    
    Parameters
    ----------
    lon1 : float or numpy.ndarray
        Longitudes of first point(s)
    lat1 : float or numpy.ndarray
        Latitudes of first point(s)
    lon2 : float or numpy.ndarray
        Longitudes of second point(s)
    lat2 : float or numpy.ndarray
        Latitudes of second point(s)
        
    Returns
    -------
    float or numpy.ndarray
        Angular distance(s) between given points
    '''
    rn1 = lon1*np.pi/180.0
    rn2 = lon2*np.pi/180.0
    rt1 = lat1*np.pi/180.0
    rt2 = lat2*np.pi/180.0
    dist = np.arccos(np.sin(rt1)*np.sin(rt2)+np.cos(rt1)*np.cos(rt2)*np.cos(rn1-rn2))
    return dist*180.0/np.pi

def eq2tl_coords(lon,lat,substellar=0.0):
    '''Compute tidally-locked coordinates of a set of equatorial lat-lon coordinates. 
    
    Transforms equatorial coordinates into a tidally-locked coordinate system where 0 degrees longitude is the substellar-south pole-antistellar meridian, and 90 degrees latitude is the substellar point, such that the evening hemisphere is 0-180 degrees longitude, the morning hemisphere is 180-360 degrees longitude, the north equatorial pole is at (0, 180), and easterly flow is counter-clockwise. Note that this differs from the coordinate system introduced in Koll & Abbot (2015) in that theirs is a left-handed coordinate system, with the south pole at (0, 180) and counter-clockwise easterly flow, which represents a south-facing observer inside the sphere, while ours is a right-handed coordinate system, representing a south-facing observer outside the sphere, which is the usual convention for spherical coordinate systems.
    
    Parameters
    ----------
    lon : numpy.ndarray
        Longitudes in equatorial coordinates [degrees]
    lat : numpy.ndarray
        Latitudes in equatorial coordinates [degrees]
    substellar : float, optional
        Longitude of the substellar point. [degrees]
        
    Returns
    -------
    numpy.ndarray, numpy.ndarray
        Transformed longitudes and latitudes [degrees]
        
    '''
    lons,lats = np.meshgrid(lon,lat)
    tlons = np.zeros(lons.shape)
    tlats = np.zeros(lats.shape)
    elons = lons*np.pi/180.0
    elats = lats*np.pi/180.0
    rss = substellar*np.pi/180.0
    tlats = np.arcsin(np.cos(elats)*np.cos(elons-rss))*180.0/np.pi
    tlons = np.arctan2(np.sin(elons-rss),-np.tan(elats))*180.0/np.pi
    tlons[tlons<0] += 360.0
    return tlons,tlats

def tl2eq_coords(lon,lat,substellar=0.0):
    '''Compute equatorial coordinates of a set of tidally-locked lat-lon coordinates. 
    
    Transforms tidally-locked coordinates into the standard equatorial coordinate system. Note that in our tidally-locked coordinate system, 0 degrees longitude is the substellar-south pole-antistellar meridian, and 90 degrees latitude is the substellar point, such that the evening hemisphere is 0-180 degrees longitude, the morning hemisphere is 180-360 degrees longitude, the north equatorial pole is at (0, 180), and easterly flow is counter-clockwise. Note that this differs from the coordinate system introduced in Koll & Abbot (2015) in that theirs is a left-handed coordinate system, with the south pole at (0, 180) and counter-clockwise easterly flow, which represents a south-facing observer inside the sphere, while ours is a right-handed coordinate system, representing a south-facing observer outside the sphere, which is the usual convention for spherical coordinate systems.
    
    Parameters
    ----------
    lon : numpy.ndarray
        Longitudes in tidally-locked coordinates [degrees]
    lat : numpy.ndarray
        Latitudes in tidally-locked coordinates [degrees]
    substellar : float, optional
        Longitude of the substellar point. [degrees]
        
    Returns
    -------
    numpy.ndarray, numpy.ndarray
        Transformed longitudes and latitudes [degrees]
        
    '''
    lons,lats = np.meshgrid(lon,lat)
    qlons = np.zeros(lons.shape)
    qlats = np.zeros(lats.shape)
    tlons = lons*np.pi/180.0
    tlats = lats*np.pi/180.0
    rss = substellar*np.pi/180.0
    qlats = np.arcsin(-np.cos(tlats)*np.cos(tlons))
    qlons = np.arctan2(np.cos(tlats)*np.sin(tlons),np.sin(tlats)) + rss
    qlats *= 180.0/np.pi
    qlons *= 180.0/np.pi
    qlons[qlons<0] += 360.0
    return qlons,qlats


def eq2tl(variable,lon,lat,substellar=0.0, polemethod="interp"):
    '''Transform a variable to tidally-locked coordinates

    Note that in our tidally-locked coordinate system, 0 degrees longitude is the substellar-south pole-antistellar meridian, and 90 degrees latitude is the substellar point, such that the evening hemisphere is 0-180 degrees longitude, the morning hemisphere is 180-360 degrees longitude, the north equatorial pole is at (0, 180), and easterly flow is counter-clockwise. Note that this differs from the coordinate system introduced in Koll & Abbot (2015) in that theirs is a left-handed coordinate system, with the south pole at (0, 180) and counter-clockwise easterly flow, which represents a south-facing observer inside the sphere, while ours is a right-handed coordinate system, representing a south-facing observer outside the sphere, which is the usual convention for spherical coordinate systems.

    Parameters
    ----------
    variable : numpy.ndarray (2D, 3D, or 4D)
        N-D data array to be transformed. Final two dimensions must be (lat,lon)
    lon : numpy.ndarray
        1D array of longitudes [deg]
    lat : numpy.ndarray
        1D array of latitudes [deg]
    substellar : float, optional
        Longitude of the substellar point (defaults to 0 degrees)
    polemethod : str, optional
        Interpolation method for polar latitudes. If "nearest", then instead of inverse-distance linear interpolation, will use nearest-neighbor. This is recommended for vector variables. For scalars, leave as "interp".
        
    Returns
    -------
    numpy.ndarray, numpy.ndarray, numpy.ndarray
        Transformed longitudes, latitudes, and data array.
    '''
    tlon = np.copy(lon)
    tlat = np.copy(lat)
    tlvariable = np.zeros(variable.shape)
    nlon = len(lon)
    nlat = len(lat)
    elons,elats = tl2eq_coords(tlon,tlat,substellar=substellar)
    elons *= np.pi/180.0
    elats *= np.pi/180.0
    for i in range(nlon):
        for j in range(nlat):
            rlon = elons[j,i]
            rlat = elats[j,i]
            dlon = rlon*180.0/np.pi
            dlat = rlat*180.0/np.pi
            if abs(dlat)>abs(lat).max():
                if polemethod!="nearest":
                    jj = abs(lat-dlat).argmin()
                    distances = np.zeros(nlon)
                    for ii in range(nlon):
                        distances[ii] = 1.0/adist(lon[ii],dlon,lat[jj],dlat)
                    tlvariable[...,j,i] = np.average(variable[...,jj,:],weights=distances,axis=len(variable.shape)-2)
                else:
                    ilat=np.argmin(abs(dlat-lat))
                    ilon=np.argmin(abs(dlon-lon))
                    tlvariable[...,j,i] = variable[...,ilat,ilon]
            else:
                latcomparison = abs(lat-dlat)
                loncomparison = abs(lon-dlon)
                colat = (latcomparison.min()==0.0) #We are colatitude
                colon = (loncomparison.min()==0.0) #We are colongitude
                if not colat and not colon:
                    jj = latcomparison.argmin()
                    ii = loncomparison.argmin()
                    if ii==0:
                        ln1 = abs(loncomparison[-1]-360.0-dlon) #if chosen, our indices will be -1 and 0
                        ln2 = loncomparison[ 1] #if chosen, our indices will be 0 and 1
                    elif ii==nlon-1:
                        ln1 = loncomparison[ii-1] #if chosen, our indices will be nlon-2 and nlon-1
                        ln2 = abs(360.0-dlon) #if chosen, our indices will be nlon-1 and 0
                    else:
                        ln1 = loncomparison[ii-1]
                        ln2 = loncomparison[ii+1]
                    if jj==0:
                        lt1 = abs(100.0-dlat) #shouldn't be chosen
                        lt2 = latcomparison[1]
                    elif jj==nlat-1:
                        lt1 = latcomparison[jj-1]
                        lt2 = abs(-100.0-dlat) #shouldn't be chosen
                    else:
                        lt1 = latcomparison[jj-1]
                        lt2 = latcomparison[jj+1]
                    ni = (ii + 2*np.argmin([ln1,ln2])-1)%nlon
                    nj = jj + 2*np.argmin([lt1,lt2])-1
                    if nj>=nlat or nj<0:
                        print(ii,jj,ni,nj,dlon,dlat)
                        raise
                    d11 = 1.0/adist(lon[ii],dlon,lat[jj],dlat)
                    d22 = 1.0/adist(lon[ni],dlon,lat[nj],dlat)
                    d21 = 1.0/adist(lon[ii],dlon,lat[nj],dlat)
                    d12 = 1.0/adist(lon[ni],dlon,lat[jj],dlat)
                    tlvariable[...,j,i] = np.average(np.array([variable[...,jj,ii],variable[...,jj,ni],
                                                            variable[...,nj,ii],variable[...,nj,ni]]),
                                                  weights=[d11,d12,d21,d22],axis=0)
                elif colat: #only changing longitude
                    ii = loncomparison.argmin()
                    jj = latcomparison.argmin()
                    if ii==0:
                        ln1 = abs(loncomparison[-1]-360.0-dlon) #if chosen, our indices will be -1 and 0
                        ln2 = loncomparison[ 1] #if chosen, our indices will be 0 and 1
                    elif ii==nlon-1:
                        ln1 = loncomparison[ii-1] #if chosen, our indices will be nlon-2 and nlon-1
                        ln2 = abs(360.0-dlon) #if chosen, our indices will be nlon-1 and 0
                    else:
                        ln1 = loncomparison[ii-1]
                        ln2 = loncomparison[ii+1]
                    ni = (ii + 2*np.argmin([ln1,ln2])-1)%nlon
                    d1 = 1.0/adist(lon[ii],dlon,lat[jj],dlat)
                    d2 = 1.0/adist(lon[ni],dlon,lat[jj],dlat)
                    tlvariable[...,j,i] = np.average(np.array([variable[...,jj,ii],variable[...,jj,ni]]),
                                                  weights = [d1,d2],axis=0)
                elif colon: #only changing latitude
                    ii = loncomparison.argmin()
                    jj = latcomparison.argmin()
                    if jj==0:
                        lt1 = abs(100.0-dlat) #shouldn't be chosen
                        lt2 = latcomparison[1]
                    elif jj==nlat-1:
                        lt1 = latcomparison[jj-1]
                        lt2 = abs(-100.0-dlat) #shouldn't be chosen
                    else:
                        lt1 = latcomparison[jj-1]
                        lt2 = latcomparison[jj+1]
                    nj = jj + 2*np.argmin([lt1,lt2])-1
                    d1 = 1.0/adist(lon[ii],dlon,lat[jj],dlat)
                    d2 = 1.0/adist(lon[ii],dlon,lat[nj],dlat)
                    tlvariable[...,j,i] = np.average(np.array([variable[...,jj,ii],variable[...,nj,ii]]),
                                                  weights = [d1,d2],axis=0)
                else: #We coincide with a real point
                    ii = loncomparison.argmin()
                    jj = latcomparison.argmin()
                    tlvariable[...,j,i] = variable[...,jj,ii]
    return tlon,tlat,tlvariable

        
def tl2eq(variable,lon,lat,substellar=0.0):
    '''Transform a tidally-locked variable to standard equatorial coordinates

    Note that in our tidally-locked coordinate system, 0 degrees longitude is the substellar-south pole-antistellar meridian, and 90 degrees latitude is the substellar point, such that the evening hemisphere is 0-180 degrees longitude, the morning hemisphere is 180-360 degrees longitude, the north equatorial pole is at (0, 180), and easterly flow is counter-clockwise. Note that this differs from the coordinate system introduced in Koll & Abbot (2015) in that theirs is a left-handed coordinate system, with the south pole at (0, 180) and counter-clockwise easterly flow, which represents a south-facing observer inside the sphere, while ours is a right-handed coordinate system, representing a south-facing observer outside the sphere, which is the usual convention for spherical coordinate systems.

    Parameters
    ----------
    variable : numpy.ndarray (2D, 3D, or 4D)
        N-D data array to be transformed. Final two dimensions must be (lat,lon)
    lon : numpy.ndarray
        1D array of longitudes [deg]
    lat : numpy.ndarray
        1D array of latitudes [deg]
    substellar : float, optional
        Longitude of the substellar point (defaults to 0 degrees)
        
    Returns
    -------
    numpy.ndarray, numpy.ndarray, numpy.ndarray
        Transformed longitudes, latitudes, and data array.
    '''
    qlon = np.copy(lon)
    qlat = np.copy(lat)
    eqvariable = np.zeros(variable.shape)
    nlon = len(lon)
    nlat = len(lat)
    qlons,qlats = eq2tl_coords(qlon,qlat,substellar=substellar)
    qlons *= np.pi/180.0
    qlats *= np.pi/180.0
    for i in range(nlon):
        for j in range(nlat):
            rlon = qlons[j,i]
            rlat = qlats[j,i]
            dlon = rlon*180.0/np.pi
            dlat = rlat*180.0/np.pi
            if abs(dlat)>abs(lat).max():
                jj = abs(lat-dlat).argmin()
                distances = np.zeros(nlon)
                for ii in range(nlon):
                    distances[ii] = 1.0/adist(lon[ii],dlon,lat[jj],dlat)
                eqvariable[...,j,i] = np.average(variable[...,jj,:],weights=distances,axis=len(variable.shape)-2)
            else:
                latcomparison = abs(lat-dlat)
                loncomparison = abs(lon-dlon)
                colat = (latcomparison.min()==0.0) #We are colatitude
                colon = (loncomparison.min()==0.0) #We are colongitude
                if not colat and not colon:
                    jj = latcomparison.argmin()
                    ii = loncomparison.argmin()
                    if ii==0:
                        ln1 = abs(loncomparison[-1]-360.0-dlon) #if chosen, our indices will be -1 and 0
                        ln2 = loncomparison[ 1] #if chosen, our indices will be 0 and 1
                    elif ii==nlon-1:
                        ln1 = loncomparison[ii-1] #if chosen, our indices will be nlon-2 and nlon-1
                        ln2 = abs(360.0-dlon) #if chosen, our indices will be nlon-1 and 0
                    else:
                        ln1 = loncomparison[ii-1]
                        ln2 = loncomparison[ii+1]
                    if jj==0:
                        lt1 = abs(100.0-dlat) #shouldn't be chosen
                        lt2 = latcomparison[1]
                    elif jj==nlat-1:
                        lt1 = latcomparison[jj-1]
                        lt2 = abs(-100.0-dlat) #shouldn't be chosen
                    else:
                        lt1 = latcomparison[jj-1]
                        lt2 = latcomparison[jj+1]
                    ni = (ii + 2*np.argmin([ln1,ln2])-1)%nlon
                    nj = jj + 2*np.argmin([lt1,lt2])-1
                    if nj>=nlat or nj<0:
                        print(ii,jj,ni,nj,dlon,dlat)
                        raise
                    d11 = 1.0/adist(lon[ii],dlon,lat[jj],dlat)
                    d22 = 1.0/adist(lon[ni],dlon,lat[nj],dlat)
                    d21 = 1.0/adist(lon[ii],dlon,lat[nj],dlat)
                    d12 = 1.0/adist(lon[ni],dlon,lat[jj],dlat)
                    eqvariable[...,j,i] = np.average(np.array([variable[...,jj,ii],variable[...,jj,ni],
                                                            variable[...,nj,ii],variable[...,nj,ni]]),
                                                  weights=[d11,d12,d21,d22],axis=0)
                elif colat: #only changing longitude
                    ii = loncomparison.argmin()
                    jj = latcomparison.argmin()
                    if ii==0:
                        ln1 = abs(loncomparison[-1]-360.0-dlon) #if chosen, our indices will be -1 and 0
                        ln2 = loncomparison[ 1] #if chosen, our indices will be 0 and 1
                    elif ii==nlon-1:
                        ln1 = loncomparison[ii-1] #if chosen, our indices will be nlon-2 and nlon-1
                        ln2 = abs(360.0-dlon) #if chosen, our indices will be nlon-1 and 0
                    else:
                        ln1 = loncomparison[ii-1]
                        ln2 = loncomparison[ii+1]
                    ni = (ii + 2*np.argmin([ln1,ln2])-1)%nlon
                    d1 = 1.0/adist(lon[ii],dlon,lat[jj],dlat)
                    d2 = 1.0/adist(lon[ni],dlon,lat[jj],dlat)
                    eqvariable[...,j,i] = np.average(np.array([variable[...,jj,ii],variable[...,jj,ni]]),
                                                  weights = [d1,d2],axis=0)
                elif colon: #only changing latitude
                    ii = loncomparison.argmin()
                    jj = latcomparison.argmin()
                    if jj==0:
                        lt1 = abs(100.0-dlat) #shouldn't be chosen
                        lt2 = latcomparison[1]
                    elif jj==nlat-1:
                        lt1 = latcomparison[jj-1]
                        lt2 = abs(-100.0-dlat) #shouldn't be chosen
                    else:
                        lt1 = latcomparison[jj-1]
                        lt2 = latcomparison[jj+1]
                    nj = jj + 2*np.argmin([lt1,lt2])-1
                    d1 = 1.0/adist(lon[ii],dlon,lat[jj],dlat)
                    d2 = 1.0/adist(lon[ii],dlon,lat[nj],dlat)
                    eqvariable[...,j,i] = np.average(np.array([variable[...,jj,ii],variable[...,nj,ii]]),
                                                  weights = [d1,d2],axis=0)
                else: #We coincide with a real point
                    ii = loncomparison.argmin()
                    jj = latcomparison.argmin()
                    eqvariable[...,j,i] = variable[...,jj,ii]
    return qlon,qlat,eqvariable


def eq2tl_uv(u, v, lon, lat, substellar=0.0):
    '''Transform velocity variables to tidally-locked coordinates

    Parameters
    ----------
    u : numpy.ndarray (2D, 3D, or 4D)
        N-D data array of zonal velocities to be transformed. Final two dimensions must be (lat,lon)
    v : numpy.ndarray (2D, 3D, or 4D)
        N-D data array of meridional velocities to be transformed. Final two dimensions must be (lat,lon)
    lon : numpy.ndarray
        1D array of longitudes [deg]
    lat : numpy.ndarray
        1D array of latitudes [deg]
    substellar : float, optional
        Longitude of the substellar point (defaults to 0 degrees)
        
    Returns
    -------
    numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray
        Transformed longitudes, latitudes, and velocity data arrays.
    '''
    lon_tl,lat_tl,uq_tl = eq2tl(u,lon,lat,substellar=substellar,polemethod='nearest')
    lon_tl,lat_tl,vq_tl = eq2tl(v,lon,lat,substellar=substellar,polemethod='nearest')
    #lon_tl,lat_tl = eq2tl_coords(lon,lat,substellar=substellar)
    lons,lats = np.meshgrid(lon,lat)
    rlons = substellar*np.pi/180.0 - lons*np.pi/180.0
    rlats = lats*np.pi/180.0
    qlons,qlats = tl2eq_coords(lon,lat,substellar=substellar)
    rqlons = substellar*np.pi/180.0 - qlons*np.pi/180.0
    rqlats = qlats*np.pi/180.0
    rtlons,rtlats = np.meshgrid(lon_tl,lat_tl)
    rtlons *= np.pi/180.0
    rtlats *= np.pi/180.0
    if len(u.shape)==3:
        rlons = rlons[np.newaxis,:,:]
        rlats = rlats[np.newaxis,:,:]
        rqlons = rqlons[np.newaxis,:,:]
        rqlats = rqlats[np.newaxis,:,:]
        rtlons = rtlons[np.newaxis,:,:]
        rtlats = rtlats[np.newaxis,:,:]
    elif len(u.shape)==4:
        rqlons = rqlons[np.newaxis,np.newaxis,:,:]
        rqlats = rqlats[np.newaxis,np.newaxis,:,:]
        rlons = rlons[np.newaxis,np.newaxis,:,:]
        rlats = rlats[np.newaxis,np.newaxis,:,:]
        rtlons = rtlons[np.newaxis,np.newaxis,:,:]
        rtlats = rtlats[np.newaxis,np.newaxis,:,:]
        
    ufactor = -np.cos(rtlats)/(np.sin(rqlats)*(1+np.sin(rqlons)**2/np.tan(rqlats)**2))
    vfactor = np.sin(rqlats)/np.sqrt(1-np.cos(rqlats)**2*np.cos(rqlons)**2)
    
    u_tl =  (vq_tl*np.sin(rqlons)/np.sin(rqlats) + uq_tl*np.cos(rqlons))*ufactor
    v_tl =  (uq_tl*np.sin(rqlons)/np.sin(rqlats) - vq_tl*np.cos(rqlons))*vfactor
            
    return lon_tl,lat_tl,u_tl,v_tl    


def tlstream(dataset,plarad=6371.0e3,grav=9.80665,substellar=0.0):
    '''Compute the tidally-locked streamfunction
    
    Parameters
    ----------
    dataset : str or netCDF4.Dataset
        Either path to netCDF4 Dataset of model output or an instance of the dataset.
    plarad : float, optional
        Planetary radius [m]
    grav : float, optional
        Surface gravity [m/s^2]
    substellar : float, optional
        Longitude of the substellar point in degrees.
            
    Returns
    -------
    numpy.ndarray(1D), numpy.ndarray(1D), numpy.ndarray(2D)
        tidally-locked latitude, layer interface pressures, and TL streamfunction
    '''
    from scipy.integrate import cumtrapz
    
    if type(dataset)==str:
        dataset = nc.Dataset(dataset,"r")
    lon = dataset.variables['lon'][:]
    lat = dataset.variables['lat'][:]
    lon_TL,lat_TL,ua_TL,va_TL = eq2tl_uv(dataset.variables['ua'][:],
                                         dataset.variables['va'][:],
                                         lon,lat,substellar=substellar)
    mva = np.nanmean(va_TL,axis=3) #tidally-locked meridional wind
    ps = spatialmath(dataset.variables['ps'][:],lon=lon,lat=lat)
    lev = dataset.variables['lev'][:]
    mvadp = cumtrapz(mva,x=lev[:]*ps*100.0,axis=1) #integrate in Pa not hPa
    prefactor = 2*np.pi*plarad/grav #2piR/g
    stf = prefactor*np.cos(lat_TL[np.newaxis,np.newaxis,:]*np.pi/180.0)*mvadp
    pmid = 0.5*(lev[1:]+lev[:-1])*ps
    return lat_TL,pmid,-stf

#def rhines(U,lat,lon,plarad=6371.0,daylen=15.0,beta=None):
    #'''Return the nondimensional Rhines length scale L_R/a
    
    #Parameters
    #----------
    #U : numpy.ndarray or float
        #Characteristic velocity [m/s] of the atmospheric jets
    #lat : numpy.ndarray
        #1D array of latitudes [deg]
    #lon : numpy.ndarray
        #1D array of longitudes [deg]
    #plarad : float, optional
        #Planetary radius in km 
    #daylen : float, optional
        #Planetary rotation periods in days
       
    #Returns
    #-------
       #float
    #'''
    
    ##if beta is None:
    #Omega = 2*np.pi/(daylen*86400.0)
        ##beta = np.gradient(2*Omega*np.sin(lat*np.pi/180.0),
                           ##plarad*1e3*(lat*np.pi/180.0))
        ##lons,lats = np.meshgrid(lon,lat)
        ##beta = np.ones(U.shape)*beta[np.newaxis,:,np.newaxis]
        ##beta = gt.spatialmath(beta,lat=lat,lon=lon)
    #beta = 2*Omega/(plarad*1e3) #This is the equatorial approximation
    #Lr = np.sqrt(U/beta)/(plarad*1e3)
    #return Lr

#def getheight(T,P,P0,gascon,grav,lat=None,lon=None,lapse=None):
    #ndims = len(T.shape)
    #if lapse is None:
        #altz = -gascon*T/grav*np.log(P/P0)
    #else:
        #if ndims==1:
            #T0 = T[-1]
        #elif ndims==2:
            #T0 = np.nanmean(T[-1,:])
        #elif ndims==3:
            #if lat is None:
                #lat = np.zeros(90,-90,num=T.shape[1])
            #if lon is None:
                #lon = np.zeros(0,360,num=T.shape[2])
            #T0 = gt.spatialmath(T[-1,:,:],lat=lat,lon=lon)
        #elif ndims==4:
            #if lat is None:
                #lat = np.zeros(90,-90,num=T.shape[1])
            #if lon is None:
                #lon = np.zeros(0,360,num=T.shape[2])
            #T0 = gt.spatialmath(T[:,-1,:,:],lat=lat,lon=lon)
        #altz = T0/lapse * ((P/P0)**(-lapse*gascon/grav)-1)
    #return altz

#def Ngradient(var,x,axis=0):
    #'''Assumes x has the same shape as var'''
    #ndims = len(var.shape)
    #if ndims==1: #why are you even using this function
        #return np.gradient(var,x)
    #elif ndims==2:
        #grad = np.zeros(var.shape)
        #if axis==0:
            #for k in range(var.shape[1]):
                #grad[:,k] = np.gradient(var[:,k],x[:,k])
        #else:
            #for k in range(var.shape[0]):
                #grad[k,:] = np.gradient(var[k,:],x[k,:])
    #elif ndims==3:
        #grad = np.zeros(var.shape)
        #if axis==0:
            #for j in range(var.shape[1]):
                #for k in range(var.shape[2]):
                    #grad[:,j,k] = np.gradient(var[:,j,k],x[:,j,k])
        #elif axis==1:
            #for j in range(var.shape[0]):
                #for k in range(var.shape[2]):
                    #grad[j,:,k] = np.gradient(var[j,:,k],x[j,:,k])
        #else:
            #for j in range(var.shape[0]):
                #for k in range(var.shape[1]):
                    #grad[j,k,:] = np.gradient(var[j,k,:],x[j,k,:])
    #else:
        #grad = np.zeros(var.shape)
        #if axis==0:
            #for j in range(var.shape[1]):
                #for k in range(var.shape[2]):
                    #for l in range(var.shape[3]):
                        #grad[:,j,k,l] = np.gradient(var[:,j,k,l],
                                                    #x[:,j,k,l])
        #elif axis==1:
            #for j in range(var.shape[0]):
                #for k in range(var.shape[2]):
                    #for l in range(var.shape[3]):
                        #grad[j,:,k,l] = np.gradient(var[j,:,k,l],
                                                    #x[j,:,k,l])
        #elif axis==2:
            #for j in range(var.shape[0]):
                #for k in range(var.shape[1]):
                    #for l in range(var.shape[3]):
                        #grad[j,k,:,l] = np.gradient(var[j,k,:,l],
                                                    #x[j,k,:,l])
        #else:
            #for j in range(var.shape[0]):
                #for k in range(var.shape[1]):
                    #for l in range(var.shape[2]):
                        #grad[j,k,l,:] = np.gradient(var[j,k,l,:],
                                                    #x[j,k,l,:])
    #return grad
        

#def bruntvasaila(T,P,lat=None,lon=None,lapse=False,gascon=287.0,grav=9.80665,cp=1006.0):
    #'''T and P must have the same shape'''
    #P0 = np.nanmax(P)
    #theta = T*(P0/P)**(gascon/cp)
    #ndims = len(T.shape)
    #if not lapse:
        #altz = -gascon*T/grav*np.log(P/P0)
    #else: #Iteratively find the lapse rates and altitudes
        #altz = getheight(T,P,P0,gascon,grav,lapse=False,lat=lat,lon=lon)
        #if ndims==4:
            #lapse = Ngradient(T,altz,axis=1)
        #else:
            #lapse = Ngradient(T,altz,axis=0)
        #for n in range(40):
            #altz = getheight(T,P,P0,gascon,grav,lapse=lapse,lat=lat,lon=lon)
            #if ndims==4:
                #lapse = Ngradient(T,altz,axis=1)
            #else:
                #lapse = Ngradient(T,altz,axis=0)
    #if ndims==4:
        #N = np.sqrt(grav/theta * Ngradient(theta,altz,axis=1))
    #else:
        #N = np.sqrt(grav/theta * Ngradient(theta,altz,axis=0))
    #return N
    
#def rossby(T,P,plarad=6371.0,daylen=15.0,gascon=287.0,grav=9.80665,
           #lat=None,lon=None,lapse=False,cp=1006.0):
    #'''Return the Rossby deformation number L_R/a
    
    #Options
    #-------
        #T : Air temperature [K] -- numpy.ndarray
        #P : Air pressure -- numpy.ndarray
        #plarad : Planetary radius in km (optional)
        #daylen : Planetary rotation period in days (optional)
        #gascon : Specific gas constant (optional)
        #grav : Surface gravity [m/s^2] (optional)
        #cp : specific heat capacity of air at constant presure [J/kg]
        #lat : Latitudes in degrees (optional)
        #lon : Longitudes in degrees (optional)
        #lapse : Whether to compute the lapse rate {True/False} (optional)
        
    #Returns
    #-------
        #float
    #'''
    #N = bruntvasaila(T,P,lat=lat,lon=lon,lapse=lapse,gascon=gascon,
                    #grav=grav,cp=cp)
    #Omega = 2*np.pi/(daylen*86400.0)
    ##beta = np.gradient(2*Omega*np.sin(lat*np.pi/180.0),
    ##                   plarad*1e3*(lat*np.pi/180.0))
    ##lons,lats = np.meshgrid(lon,lat)
    ##beta = np.ones(U.shape)*beta[np.newaxis,:,np.newaxis]
    ##beta = gt.spatialmath(beta,lat=lat,lon=lon)
    #beta = 2*Omega/(plarad*1e3)
    #scaleH = gascon*T/grav
    #Ld = np.sqrt(N*scaleH/(2*beta))/(plarad*1e3)
    #ndims = len(Ld.shape)
    #return Ld
    #if ndims<3:
        #return np.nanmean(Ld)
    #elif ndims==3:
        #return gt.spatialmath(np.nanmean(Ld,axis=0),lon=lon,lat=lat)
    #else:
        #return gt.spatialmath(np.nanmean(Ld,axis=1),lon=lon,lat=lat)

#def deform(T,lat,lon,plarad=6371.0,daylen=15.0,gascon=287.0,grav=9.80665):
    #'''Return the equatorial (Rossby) deformation length scale L_R/a
    
    #Options
    #-------
       #T : characteristic temperature [K]
       #plarad : Planetary radius in km (optional)
       #daylen : Planetary rotation periods in days (optional)
       #gascon : Specific gas constant (optional)
       #grav : Surface gravity [m/s^2] (optional)
       
    #Returns
    #-------
       #float
    #'''
    
    #Omega = 2*np.pi/(daylen*86400.0)
    ##beta = np.gradient(2*Omega*np.sin(lat*np.pi/180.0),
    ##                   plarad*1e3*(lat*np.pi/180.0))
    ##lons,lats = np.meshgrid(lon,lat)
    ##beta = np.ones(U.shape)*beta[np.newaxis,:,np.newaxis]
    ##beta = gt.spatialmath(beta,lat=lat,lon=lon)
    #beta = 2*Omega/(plarad*1e3)
    #scaleH = gascon*T/grav
    #c = np.sqrt(grav*scaleH)
    #Le = np.sqrt(c/beta)/(plarad*1e3)
    #return Le
    