"""
Read raw exoplasim output files and postprocess them into netCDF output files.
"""
import numpy as np
import struct
import exoplasim.gcmt
import exoplasim.gcmt as gcmt
import exoplasim.filesupport
from exoplasim.filesupport import SUPPORTED
import scipy, scipy.integrate, scipy.interpolate
import os, sys

'''
This module is intended to be a near-replacement for the C++ burn7 utility, which in its present
form relies on deprecated and unsupported versions of the netcdf C and C++ libraries. This utility
relies instead on the netCDF4-python package, whose usage pattern is unlikely to become deprecated
at any point in the near future. Some features of burn7 are not included; notably at present only
netCDF output is supported, rather than the additional SRV and GRaDs formats supported by burn7.
Additionally, at least at first, there are the following restrictions:

    **Vertical levels will only be output in native sigma levels,** and will not be interpolated onto
    a pressure grid. This is unlikely to be supported in the future either, as the conversion of sigma
    to pressure is trivial, and interpolating onto a different pressure grid at the postprocessing step
    risks misleading end users and is liable to run into problems with non-Earthlike surface pressures.
    
    **Horizontal grid will only be Gaussian-spaced latitude-longitude,** rather than being optionally
    spherical harmonics, fourier coefficients, or zonal averages. The additional horizontal options present
    in burn7 will be supported in pyburn in future versions.
    
    **The Mars shortcut for derived variables is unsupported.** This will be supported in future versions.
'''

def _log(destination,string):
    if destination is None:
        if string=="\n":
            print()
        else:
            print(string)
    else:
        with open(destination,"a") as f:
            f.write(string+"\n")

#dictionary that can be searched by string (integer) codes
ilibrary = { "50":["nu"   ,"true_anomaly"                    ,"deg"        ],
             "51":["lambda","ecliptic_longitude"             ,"deg"        ],
             "52":["zdec" ,"solar_declination"               ,"deg"        ],
             "53":["rdist","orbital_distance"                ,"1"          ],
             "54":["rasc" ,"solar right ascension"           ,"deg"        ],
            "110":["mld"  ,"mixed_layer_depth"               ,"m"          ], 
            "129":["sg"   ,"surface_geopotential"            ,"m2 s-2"     ],
            "130":["ta"   ,"air_temperature"                 ,"K"          ],
            "131":["ua"   ,"eastward_wind"                   ,"m s-1"      ],
            "132":["va"   ,"northward_wind"                  ,"m s-1"      ],
            "133":["hus"  ,"specific_humidity"               ,"1"          ],
            "134":["ps"   ,"surface_air_pressure"            ,"hPa"        ],
            "135":["wap"  ,"vertical_air_velocity"           ,"Pa s-1"     ],
            "137":["wa"   ,"upward_wind"                     ,"m s-1"      ],
            "138":["zeta" ,"atm_relative_vorticity"          ,"s-1"        ],
            "139":["ts"   ,"surface_temperature"             ,"K"          ],
            "140":["mrso" ,"lwe_of_soil_moisture_content"    ,"m"          ],
            "141":["snd"  ,"surface_snow_thickness"          ,"m"          ],
            "142":["prl"  ,"lwe_of_large_scale_precipitation","m s-1"      ],
            "143":["prc"  ,"convective_precipitation_rate"   ,"m s-1"      ],
            "144":["prsn" ,"lwe_of_snowfall_amount"          ,"m s-1"      ],
            "145":["bld"  ,"dissipation_in_atmosphere_bl"    ,"W m-2"      ],
            "146":["hfss" ,"surface_sensible_heat_flux"      ,"W m-2"      ],
            "147":["hfls" ,"surface_latent_heat_flux"        ,"W m-2"      ],
            "148":["stf"  ,"streamfunction"                  ,"m2 s-2"     ],
            "149":["psi"  ,"velocity_potential"              ,"m2 s-2"     ],
            "151":["psl"  ,"air_pressure_at_sea_level"       ,"hPa"        ],
            "152":["pl"   ,"log_surface_pressure"            ,"1"          ],
            "155":["d"    ,"divergence_of_wind"              ,"s-1"        ],
            "156":["zg"   ,"geopotential_height"             ,"m"          ],
            "157":["hur"  ,"relative_humidity"               ,"1"          ],
            "158":["tps"  ,"tendency_of_surface_air_pressure","Pa s-1"     ],
            "159":["u3"   ,"ustar"                           ,"m3 s-3"     ],
            "160":["mrro" ,"surface_runoff"                  ,"m s-1"      ],
            "161":["clw"  ,"liquid_water_content"            ,"1"          ],
            "162":["cl"   ,"cloud_area_fraction_in_layer"    ,"1"          ],
            "163":["tcc"  ,"total_cloud_cover"               ,"1"          ],
            "164":["clt"  ,"cloud_area_fraction"             ,"1"          ],
            "165":["uas"  ,"eastward_wind_10m"               ,"m s-1"      ],
            "166":["vas"  ,"northward_wind_10m"              ,"m s-1"      ],
            "167":["tas"  ,"air_temperature_2m"              ,"K"          ],
            "168":["td2m" ,"dew_point_temperature_2m"        ,"K"          ],
            "169":["tsa"  ,"surface_temperature_accumulated" ,"K"          ],
            "170":["tsod" ,"deep_soil_temperature"           ,"K"          ],
            "171":["dsw"  ,"deep_soil_wetness"               ,"1"          ],
            "172":["lsm"  ,"land_binary_mask"                ,"1"          ],
            "173":["z0"   ,"surface_roughness_length"        ,"m"          ],
            "174":["alb1" ,"surface_albedo_SW1"              ,"1"          ],
            "175":["alb"  ,"surface_albedo"                  ,"1"          ],
            "176":["rss"  ,"surface_net_shortwave_flux"      ,"W m-2"      ],
            "177":["rls"  ,"surface_net_longwave_flux"       ,"W m-2"      ],
            "178":["rst"  ,"toa_net_shortwave_flux"          ,"W m-2"      ],
            "179":["rlut" ,"toa_net_longwave_flux"           ,"W m-2"      ],
            "180":["tauu" ,"surface_eastward_stress"         ,"Pa"         ],
            "181":["tauv" ,"surface_northward_stress"        ,"Pa"         ],
            "182":["evap" ,"lwe_of_water_evaporation"        ,"m s-1"      ],
            "183":["tso"  ,"climate_deep_soil_temperature"   ,"K"          ],
            "184":["alb2" ,"surface_albedo_SW2"              ,"1"          ],
            "203":["rsut" ,"toa_outgoing_shortwave_flux"     ,"W m-2"      ],
            "204":["ssru" ,"surface_solar_radiation_upward"  ,"W m-2"      ],
            "205":["stru" ,"surface_thermal_radiation_upward","W m-2"      ],
            "207":["tso2" ,"soil_temperature_level_2"        ,"K"          ],
            "208":["tso3" ,"soil_temperature_level_3"        ,"K"          ],
            "209":["tso4" ,"soil_temperature_level_4"        ,"K"          ],
            "210":["sic"  ,"sea_ice_cover"                   ,"1"          ],
            "211":["sit"  ,"sea_ice_thickness"               ,"m"          ],
            "218":["snm"  ,"snow_melt"                       ,"m s-1"      ],
            "221":["sndc" ,"snow_depth_change"               ,"m s-1"      ],
            "230":["prw"  ,"atmosphere_water_vapor_content"  ,"kg m-2"     ],
            "232":["glac" ,"glacier_cover"                   ,"1"          ],
            "238":["tsn"  ,"snow_temperature"                ,"K"          ],
            "259":["spd"  ,"wind_speed"                      ,"m s-1"      ],
            "260":["pr"   ,"total_precipitation"             ,"m s-1"      ],
            "261":["ntr"  ,"net_top_radiation"               ,"W m-2"      ],
            "262":["nbr"  ,"net_bottom_radiation"            ,"W m-2"      ],
            "263":["hfns" ,"surface_downward_heat_flux"      ,"W m-2"      ],
            "264":["wfn"  ,"net_water_flux"                  ,"m s-1"      ],
            "265":["dqo3" ,"ozone_concentration"             ,"unknown"    ],
            "266":["lwth" ,"local_weathering"                ,"W_earth"    ],
            "267":["grnz" ,"ground_geopotential"             ,"m2 s-2"     ],
            "268":["icez" ,"glacier_geopotential"            ,"m2 s-2"     ],
            "269":["netz" ,"net_geopotential"                ,"m2 s-2"     ],
            "273":["dpdx" ,"d(ps)/dx"                        ,"Pa m-1"     ],
            "274":["dpdy" ,"d(ps)/dy"                        ,"Pa m-1"     ],
            "277":["hlpr" ,"half_level_pressure"             ,"Pa"         ],
            "278":["flpr" ,"full_level_pressure"             ,"Pa"         ],
            "279":["thetah","half_level_potential_temperature","K"         ],
            "280":["theta","full_level_potential_temperature","K"          ],
            "298":["vegf"     ,"forest_cover"                        ,"1"           ],
            "299":["veglai"   ,"vegetation_leaf_area_index"          ,"nondimen"    ],
            "300":["veggpp"   ,"vegetation_gross_primary_production" ,"kg C m-2 s-1"],
            "301":["vegnpp"   ,"vegetation_net_primary_production"   ,"kg C m-2 s-1"],
            "302":["vegppl"   ,"vegetation_light-lim_GPP"            ,"kg C m-2 s-1"],
            "303":["vegppw"   ,"vegetation_h2o-lim_GPP"              ,"kg C m-2 s-1"],
            "304":["vegplantc","vegetation_plant_carbon"             ,"kg C m-2"    ],
            "305":["vegsoilc" ,"vegetation_soil_carbon"              ,"kg C m-2"    ],
            "306":["vegnogrow","vegetation_no_growth_allocation"     ,"kg C m-2 s-1"],
            "307":["veghresp" ,"vegetation_heterotrophic_respiration","kg C m-2 s-1"],
            "308":["veglitter","vegetation_litter_production"        ,"kg C m-2 s-1"],
            "318":["czen" ,"cosine_solar_zenith_angle"       ,"nondimen"   ],
            "319":["wthpr","weatherable_precipitation"       ,"mm day-1"   ],
            "320":["mint" ,"minimum_temperature"             ,"K"          ],
            "321":["maxt" ,"maximum_temperature"             ,"K"          ],
            "322":["cape" ,"convective_available_potential_energy","J kg-1"],
            "323":["lnb"  ,"level_neutral_buoyancy"          ,"hPa"        ],
            "324":["sdef" ,"troposphere_entropy_deficit"     ,"1"          ],
            "325":["absz" ,"sigma=0.85_abs_vorticity"        ,"s-1"        ],
            "326":["umax" ,"maximum_potential_intensity"     ,"m s-1"      ],
            "327":["vent" ,"ventilation_index"               ,"1"          ],
            "328":["vrumax","ventilation_reduced_maximum_wind","m s-1"     ],
            "329":["gpi"  ,"genesis_potential_index"            ,"1"       ],
            "404":["dfu"  ,"shortwave_up"                    ,"W m-2"      ],
            "405":["dfd"  ,"shortwave_down"                  ,"W m-2"      ],
            "406":["dftu"  ,"longwave_up"                    ,"W m-2"      ],
            "407":["dftd"  ,"longwave_down"                  ,"W m-2"      ],
            "408":["dtdt"  ,"rad_heating_rate"               ,"K s-1"      ],
            "409":["dfdz"  ,"flux_convergence"               ,"W m-3"      ]}

#dictionary that can be searched by string codes
slibrary = {}
for key in ilibrary:
    kcode = int(key)
    name = ilibrary[key][0]
    longn = ilibrary[key][1]
    units = ilibrary[key][2]
    slibrary[name] = [kcode,longn,units]
                    
geopotcode  = 129 #done
tempcode    = 130 #done
ucode       = 131 #done
vcode       = 132 #done
humcode     = 133 #done
pscode      = 134 #done
wcode       = 135 #done
wzcode      = 137 #done
vortcode    = 138 #done
tscode      = 139 #done
stfcode     = 148 #done
vpotcode    = 149 #done
slpcode     = 151 #done
lnpscode    = 152 #done
divcode     = 155 #done
geopotzcode = 156 #done
rhumcode    = 157 #done
spdcode     = 259 #done
preccode    = 260 #done
ntopcode    = 261 #done
nbotcode    = 262 #done
nheatcode   = 263 #done
nh2ocode    = 264 #done
swatmcode   = 268 #done
lwatmcode   = 269 #done
natmcode    = 270 #done
sruncode    = 271 #done
dpsdxcode   = 273 #done
dpsdycode   = 274 #done
freshcode   = 275 #done
            
hpresscode  = 277 #done
fpresscode  = 278 #done
thetahcode  = 279 #done
thetafcode  = 280 #done


#Constants

MARS_GRAV   = 3.728
MARS_RADIUS = 3400000.0
MARS_RD     = 189.0
L_TIMES_RHOH2O = -333700000.0
RLAPSE      = 0.0065    #International Standard Atmosphere temperature lapse rate in K/m
RH2O        = 1000 * 1.380658e-23*6.0221367e+23 / 18.0153


def _getEndian(fbuffer):
    '''Determine Endian-ness of the buffer'''
    import sys
    byteorder = sys.byteorder
    if byteorder=="little":
        endian="<"
    else:
        endian=">"
    tag = struct.unpack('i',fbuffer[:4])[0]
    if tag>1024:   #If we get a huge number for what should be an 8-word header, flip the byte order
        if endian==">": #We're in a Big-Endian system, and the buffer was written in little-Endian.
            endian="<"
        else: #We're in a little-Endian system, and the buffer was written in Big-Endian.
            endian=">"
    elif tag==0: #it was an 8-byte record marker
        tag = struct.unpack('l',fbuffer[:8])[0]
        if tag>1024:   #If we get a huge number for what should be an 8-word header, flip the byte order
            if endian==">": #We're in a Big-Endian system, and the buffer was written in little-Endian.
                endian="<"
            else: #We're in a little-Endian system, and the buffer was written in Big-Endian.
                endian=">"
    return endian

def _getmarkerlength(fbuffer,en):
    '''Determine how many bytes were used for record markers'''
    tag = struct.unpack(en+'i',fbuffer[:4])[0]
    markerlength=4
    if tag==0:
        markerlength=8
    return markerlength

def _getwordlength(fbuffer,n,en,fmt='i'):
    '''Determine if we're dealing with 32-bit output or 64-bit.
    
    At the moment, all exoplasim outputs are 32-bit, even if run in 64-bit. However, we should build
    compatibility for this eventuality in the future.
    
    Parameters
    ----------
    fbuffer : bytes
        Binary bytes read from a file opened with ``mode='rb'`` and read with ``file.read()``. 
    n : int
        The index of the byte at which to check. This should be the start of the first word of the 
        variable in question.
    en : str
        Endianness, denoted by ">" or "<"
    fmt : str
        Expected type of the word--'i' for integers (longs), and 'd' for floats (doubles). 
    
    Returns
    -------
    int, str
       Word length in bytes, and format string for a word--4 for 32 bit, and 8 for 64 bit. 'i' for a 
       4-byte int, 'l' for an 8-byte int, 'f' for a 4-byte float, and 'd' for an 8-byte float.
    '''

    tag = struct.unpack(en+fmt,fbuffer[n:n+4])[0]
    wordlength=4
    if tag==0:
        wordlength=8
    if fmt=='i' and wordlength==8:
        fmt='l'
    elif fmt=='f' and wordlength==8:
        fmt='d'
    return wordlength,fmt

def _getknownwordlength(fbuffer,n,en,ml,mf):
    '''Determine word length of a data variable in cases where we know that the header is 8 words and is 32-bit.
    
    Parameters
    ----------
    fbuffer : bytes
        Binary bytes read from a file opened with ``mode='rb'`` and read with ``file.read()``. 
    n : int
        The index of the word at which to start, in bytes. A 32-bit word has length 4, so the current 
        position in words would be 4*n assuming 4-byte words, or 8*n if 64 bits and 8-byte words.
    en : str
        Endianness, denoted by ">" or "<"
    ml : int
        Length of a record marker
    mf : str
        Format of the record marker ('i' or 'l')
    
    Returns
    -------
    int, str
       Word length in bytes, and format string for a word--4 for 32 bit, and 8 for 64 bit. 
       'f' for a 4-byte float, and 'd' for an 8-byte float.
    '''
    
    htag = struct.unpack(en+mf,fbuffer[n:n+ml])
    n+=ml
    header = struct.unpack(en+8*'i',fbuffer[n:n+32])
    n+=32+ml #Add one word for restatement of header length
    dtag = struct.unpack(en+mf,fbuffer[n:n+ml])[0]
    
    dim1 = header[4]
    dim2 = header[5]
    
    length = dim1*dim2
    
    wordlength = dtag//length #dtag tells us the length of the coming record in bytes, while length is the
                              #length of the coming record in words. The ratio of the two is thus wordlength.
    if wordlength==4:
        fmt='f'
    else:
        fmt='d'
    return wordlength,fmt


def readrecord(fbuffer,n,en,ml,mf):
    '''Read a Fortran record from the buffer, starting at index n, and return the header, data, and updated n.
    
    Parameters
    ----------
    fbuffer : bytes
        Binary bytes read from a file opened with ``mode='rb'`` and read with ``file.read()``. 
    n : int
        The index of the word at which to start, in bytes. A 32-bit word has length 4, so the current 
        position in words would be 4*n assuming 4-byte words, or 8*n if 64 bits and 8-byte words.
    en : str
        Endianness, denoted by ">" or "<"
    ml : int
        Length of a record marker
    mf : str
        Format of the record marker ('i' or 'l')
        
    Returns
    -------
    array-like, array-like, int
        A tuple containing first the header, then the data contained in the record, and finally the new
        position in the buffer in bytes.
    '''
    if n<len(fbuffer):
        wl,fmt = _getknownwordlength(fbuffer,n,en,ml,mf)
                
        headerlength = int(struct.unpack(en+mf,fbuffer[n:n+ml])[0]//4)
        n+=ml
        header = struct.unpack(en+headerlength*'i',fbuffer[n:n+headerlength*4])
        n+=headerlength*4+ml #Add one word for restatement of header length (for backwards seeking)
        datalength = int(struct.unpack(en+mf,fbuffer[n:n+ml])[0]//wl)
        n+=ml
        data = struct.unpack(en+datalength*fmt,fbuffer[n:n+datalength*wl])
        n+=datalength*wl+ml #additional 4 for restatement of datalength
        return header,data,n
    else:
        raise Exception("Reached end of buffer!!!")
    
def readvariablecode(fbuffer,kcode,en,ml,mf):
    '''Seek through a binary output buffer and extract all records associated with a variable code.
    
    Note, assembling a variable list piece by piece in this way may be slower than reading **all** variables
    at once, because it requires seeking all the way through the buffer multiple times for each variable.
    This will likely only be faster if you only need a small number of variables.
    
    Parameters
    ----------
    fbuffer : bytes
        Binary bytes read from a file opened with ``mode='rb'`` and read with ``file.read()``.
    kcode : int
        The integer code associated with the variable. For possible codes, refer to the 
        ``Postprocessor Variable Codes. <postprocessor.html#postprocessor-variable-codes>`_
    en : str
        Endianness, denoted by ">" or "<"
    ml : int
        Length of a record marker
    mf : str
        Format of the record marker ('i' or 'l')
    
    Returns
    -------
    array-like, array-like
        A tuple containing first the header, then the variable data, as one concatenated 1D variable.
    '''
    n = 0
    mainheader,zsig,n = readrecord(fbuffer,n,en,ml)
    
    variable = None
    
    while n<len(fbuffer):
        
        recordn0 = n
        
        headerlength = int(struct.unpack(en+mf,fbuffer[n:n+ml])[0]//4)
        n+=ml
        header = struct.unpack(en+headerlength*'i',fbuffer[n:n+headerlength*4])
        n+=headerlength*4+ml
        datalength = struct.unpack(en+mf,fbuffer[n:n+ml])[0]
        n+=ml
        if header[0]==kcode:
            dataheader = header
            wl, fmt = _getknownwordlength(fbuffer,recordn0,en,ml,mf)
            datalength = int(datalength//wl)
            if not variable:
                variable = np.array(struct.unpack(en+datalength*fmt,fbuffer[n:n+datalength*wl]))
            else:
                variable = np.append(variable,struct.unpack(en+datalength*fmt,fbuffer[n:n+datalength*wl]))
            n+=datalength*wl+ml
        else: #Fast-forward past this variable without reading it.
            n+=datalength+ml
    
    return dataheader, variable

def _gettimevar(fbuffer):
    '''Extract the time array, as an array of timesteps'''
    
    en = _getEndian(fbuffer)
    ml,mf = _getwordlength(fbuffer,0,en)
    
    kcode = 139 #Use surface temperature to do this
    time = []
    n = 0
    mainheader,zsig,n = readrecord(fbuffer,n,en,ml)
    
    variable = None
    
    while n<len(fbuffer):
        
        recordn0 = n
        
        headerlength = int(struct.unpack(en+mf,fbuffer[n:n+ml])[0]//4)
        n+=ml
        header = struct.unpack(en+headerlength*'i',fbuffer[n:n+headerlength*4])
        n+=headerlength*4+ml
        datalength = struct.unpack(en+mf,fbuffer[n:n+ml])[0]
        n+=ml
        if header[0]==kcode:
            time.append(header[6]) #nstep-nstep1 (timesteps since start of run)
        n+=datalength+ml
    
    return time
    
def readallvariables(fbuffer):
    '''Extract all variables and their headers from a file byte buffer.
    
    Doing this and then only keeping the codes you want may be faster than extracting variables one by one,
    because it only needs to seek through the file one time.
    
    Parameters
    ----------
    fbuffer : bytes
        Binary bytes read from a file opened with ``mode='rb'`` and read with ``file.read()``.
    
    Returns
    -------
    dict, dict
        A dictionary containing all variable headers (by variable code), and a dictionary containing all
        variables, again by variable code.
    '''
    
    en = _getEndian(fbuffer)
    ml,mf = _getwordlength(fbuffer,0,en)
    
    n=0
    mainheader,zsig,n = readrecord(fbuffer,n,en,ml,mf)
    
    headers= {'main':mainheader}
    variables = {'main':zsig}
    nlev=mainheader[6]
    variables["sigmah"] = zsig[:nlev]
    variables["time"] = []
    
    while n<len(fbuffer):
        header,field,n = readrecord(fbuffer,n,en,ml,mf)
        kcode = str(header[0])
        if int(kcode)==139:
            variables["time"].append(header[6]) #nstep-nstep1 (timesteps since start of run)
        if kcode not in variables:
            variables[kcode] = np.array(field)
            headers[kcode] = header
        else:
            variables[kcode] = np.append(variables[kcode],field)
    
    return headers, variables
    

def refactorvariable(variable,header,nlev=10):
    '''Given a 1D data array extracted from a file with :py:func:`readrecord <exoplasim.pyburn.readrecord>`, reshape it into its appropriate dimensions.
    
    Parameters
    ----------
    variable : array-like
        Data array extracted from an output file using :py:func:`readrecord <exoplasim.pyburn.readrecord>`.
        Can also be the product of a concatenated file assembled with
        :py:func:`readvariable <exoplasim.pyburn.readvariable>`.
    header : array-like
        The header array extracted from the record associated with ``variable``. This header contains
        dimensional information.
    nlev : int, optional
        The number of vertical levels in the variable. If 1, vertical levels will not be a dimension in
        the output variable.
        
    Returns
    -------
    numpy.ndarray
        A numpy array with dimensions inferred from the header.
    '''
    dim1 = max(header[4],header[5])
    dim2 = min(header[4],header[5])
    if header[1]==1:
        nlevs=nlev
        if len(variable)%(float(len(variable))/(dim1*dim2*nlevs))!=0:
            nlevs+=1
    else:
        nlevs=1
    ntimes = int(len(variable)//(dim1*dim2*nlevs))
    if nlevs==1:
        if dim2==1:
            if dim1==1:
                newvar = np.reshape(variable,(ntimes,))
            else:
                newvar = np.reshape(variable,(ntimes,dim1))
        else:
            newvar = np.reshape(variable,(ntimes,dim2,dim1))
    else:
        if dim2==1:
            if dim1==1:
                newvar = np.reshape(variable,(ntimes,nlevs))
            else:
                newvar = np.reshape(variable,(ntimes,nlevs,dim1))
        else:
            newvar = np.reshape(variable,(ntimes,nlevs,dim2,dim1))
            
    return newvar

def readfile(filename):
    '''Extract all variables from a raw plasim output file and refactor them into the right shapes
    
    This routine will only produce what it is in the file; it will not compute derived variables.
    
    Parameters
    ----------
    filename : str
        Path to the output file to read
        
    Returns
    -------
    dict
        Dictionary of model variables, indexed by numerical code
    '''
    
    
    with open(filename,"rb") as fb:
        fbuffer = fb.read()
    
    headers, variables = readallvariables(fbuffer)
    
    nlevs = len(variables['sigmah'])
    sigmah = variables['sigmah']
    
    kcodes = list(variables.keys())
    kcodes.remove('sigmah')
    
    #time = _gettimevar(fbuffer) #This means we do one additional sweep through the file
    time = variables["time"]
    kcodes.remove("time")
    
    data = {}
    
    for key in kcodes:
        data[key] = refactorvariable(variables[key],headers[key],nlev=nlevs)
    
    nlat = min(headers['main'][4],headers['main'][5])
    nlon = max(headers['main'][4],headers['main'][5])
    ntru = headers['main'][7]
    ntimes = len(time)
    
    
    if sys.version[0]=="2":
        if nlat in [192,320]:
            import exoplasim.pyfft991v2 as pyfft
        else:
            import exoplasim.pyfft2 as pyfft
    else:
        if nlat in [192,320]:
            import exoplasim.pyfft991 as pyfft
        else:
            import exoplasim.pyfft as pyfft
    
        
    sid,gwd = pyfft.inigau(nlat)
    rlat = np.arcsin(sid)
    lat = rlat*180.0/np.pi
    
    lon = np.arange(nlon)/float(nlon)*360.0
    rlon = lon*np.pi/180.0
    
    sigmab = np.append([0,],sigmah)
    sigma = 0.5*(sigmab[0:-1]+sigmab[1:]) #Mid-layer sigma
    
    data["lev"] = sigma[:]
    data["lat"] = lat[:]
    data["lon"] = lon[:]
    data["time"]= time[:]

    return data

def _transformvar(lon,lat,variable,meta,nlat,nlon,nlev,ntru,ntime,mode='grid',
                  substellarlon=180.0,physfilter=False,zonal=False,presync=False):
    '''Ensure a variable is in a given horizontal mode.
    
    Parameters
    ----------
    lon : array-like
        Longitude array, in degrees
    lat : array-like
        Latitude array, in degrees
    variable : array-like
        Data array from dataset
    meta : list
        Meta array from dataset
    nlat : int 
        Number of latitudes
    nlon : int 
        Number of longitudes
    nlev : int 
        Number of vertical levels
    ntru : int 
        Truncation wavenumber
    ntime : int 
        Number of output frames
    mode : str, optional
        Horizontal output mode. Can be 'grid', meaning the Gaussian latitude-longitude grid used
        in ExoPlaSim, 'spectral', meaning spherical harmonics, 
        'fourier', meaning Fourier coefficients and latitudes, 'synchronous', meaning a
        Gaussian latitude-longitude grid in the synchronous coordinate system defined in
        Paradise, et al (2021), with the north pole centered on the substellar point, or
        'syncfourier', meaning Fourier coefficients computed along the dipolar meridians in the
        synchronous coordinate system (e.g. the substellar-antistellar-polar meridian, which is 0 degrees,
        or the substellar-evening-antistellar-morning equatorial meridian, which is 90 degrees). Because this
        will get assigned to the original latitude array, that will become -90 degrees for the polar
        meridian, and 0 degrees for the equatorial meridian, identical to the typical equatorial coordinate
        system.
    substellarlon : float, optional
        If mode='synchronous', the longitude of the substellar point in equatorial coordinates,
        in degrees
    physfilter : bool, optional
        Whether or not a physics filter should be used when transforming spectral variables to
        Fourier or grid domains
    zonal : bool, optional
        For grid modes ("grid" and "synchronous"), compute and output zonal means
    presync : bool, optional
        If True, the data have already been rotated into a synchronous coordinate system
        
    Returns
    -------
    numpy.ndarray
        Transformed array
    '''
    
    if sys.version[0]=="2":
        if nlat in [192,320]:
            import exoplasim.pyfft991v2 as pyfft
        else:
            import exoplasim.pyfft2 as pyfft
    else:
        if nlat in [192,320]:
            import exoplasim.pyfft991 as pyfft
        else:
            import exoplasim.pyfft as pyfft
    
    if nlev in variable.shape:
        levd = "lev"
    elif nlev+1 in variable.shape:
        levd = "levp"
        
    if mode=="grid":
        if (ntru+1)*(ntru+2) in variable.shape: #spectral variable
            if len(variable.shape)==3: #Include lev
                nlevs = variable.shape[1]
                ntimes = variable.shape[0]
                spvar = np.asfortranarray(
                           np.transpose(np.reshape(variable,
                                                   (ntimes*nlevs,variable.shape[2])))
                           )
                gridshape = (ntimes,nlevs,nlat,nlon)
                dims = ["time",levd,"lat","lon"]
            else:
                ntimes = variable.shape[0]
                spvar = np.asfortranarray(
                           np.transpose(np.reshape(variable,
                                                   (ntimes,variable.shape[1])))
                           )
                gridshape = (ntimes,nlat,nlon)
                dims = ["time","lat","lon"]
            gridvar = pyfft.sp2gp(spvar,nlat,nlon,ntru,int(physfilter))
            gridvar = np.reshape(np.transpose(gridvar),gridshape)
            gridvar = np.roll(gridvar,nlon//2,axis=-1)
        else: #grid variable
            gridvar = variable
            if len(variable.shape)==3:
                dims = ["time","lat","lon"]
            elif len(variable.shape)==2: #column
                dims = ["time",levd]
            elif len(variable.shape)==1: #scalar
                dims = ["time",]
            else:
                dims = ["time",levd,"lat","lon"]
        if meta[0]=="hus":
            gridvar[gridvar<0] = 0.0
        if zonal and "lon" in dims:
            gridvar = np.nanmean(gridvar,axis=-1)
            dims.remove("lon")
        meta.append(tuple(dims))
        outvar = gridvar
        
    elif mode=="synchronous":
        if (ntru+1)*(ntru+2) in variable.shape: #spectral variable
            if len(variable.shape)==3: #Include lev
                nlevs = variable.shape[1]
                ntimes = variable.shape[0]
                spvar = np.asfortranarray(
                           np.transpose(np.reshape(variable,
                                                   (ntimes*nlevs,variable.shape[2])))
                           )
                gridshape = (ntimes,nlevs,nlat,nlon)
                dims = ["time",levd,"lat","lon"]
            else:
                ntimes = variable.shape[0]
                spvar = np.asfortranarray(
                           np.transpose(np.reshape(variable,
                                                   (ntimes,variable.shape[1])))
                           )
                gridshape = (ntimes,nlat,nlon)
                dims = ["time","lat","lon"]
            gridvar = pyfft.sp2gp(spvar,nlat,nlon,ntru,int(physfilter))
            gridvar = np.reshape(np.transpose(gridvar),gridshape)
            gridvar = np.roll(gridvar,nlon//2,axis=-1)
        else: #grid variable
            gridvar = variable
            if len(variable.shape)==3:
                dims = ["time","lat","lon"]
            elif len(variable.shape)==2: #column
                dims = ["time",levd]
            elif len(variable.shape)==1: #scalar
                dims = ["time",]
            else:
                dims = ["time",levd,"lat","lon"]
        if meta[0]=="hus":
            gridvar[gridvar<0] = 0.0
        if not presync and "lat" in dims:
            lon,lat,tlgridvar = gcmt.eq2tl(gridvar,lon,lat,substellar=substellarlon,
                                           polemethod='interp') #fine bc all vectors are derived
        else:
            tlgridvar = gridvar
        
        if zonal and "lon" in dims:
            tlgridvar = np.nanmean(tlgridvar,axis=-1)
            dims.remove("lon")
        meta.append(tuple(dims))
        outvar = tlgridvar
        
    elif mode=="spectral":
        if (ntru+1)*(ntru+2) in variable.shape: #spectral variable
            specvar = variable
            if len(variable.shape)==3:
                dims = ("time",levd,"modes","complex")
            else:
                dims = ("time","modes","complex")
        else:
            if len(variable.shape)==4: #Include lev
                nlevs = variable.shape[1]
                ntimes = variable.shape[0]
                gpvar = np.asfortranarray(
                           np.transpose(np.reshape(variable,
                                                   (ntimes*nlevs,nlat,nlon)))
                           )
                dims = ("time",levd,"modes","complex")
            elif len(variable.shape)==2: #column
                outvar = variable
                dims = ("time",levd)
                meta.append(dims)
            elif len(vairable.shape)==1: #scalar
                outvar = variable
                dims = ("time",)
                meta.append(dims)
            else:
                ntimes = variable.shape[0]
                gpvar = np.asfortranarray(
                           np.transpose(np.reshape(variable,
                                                   (ntimes,nlat,nlon)))
                           )
                dims = ("time","modes","complex")
            if "modes" in dims:
                spvar = pyfft.gp2sp(gpvar,nlat,nlon,ntru,int(physfilter))
                specvar = np.transpose(spvar)
        if "modes" in dims:
            shape = list(specvar.shape)
            shape[-1] //= 2
            shape.append(2)
            shape = tuple(shape)
            specvar = np.reshape(specvar,shape)
            meta.append(dims)
            outvar = specvar
        
    elif mode=="fourier":
        if (ntru+1)*(ntru+2) in variable.shape: #spectral variable
            if len(variable.shape)==3: #Include lev
                nlevs = variable.shape[1]
                ntimes = variable.shape[0]
                spvar = np.asfortranarray(
                           np.transpose(np.reshape(variable,
                                                   (ntimes*nlevs,variable.shape[2])))
                           )
                fcshape = (ntimes,nlevs,nlat,nlon//2,2)
                dims = ["time",levd,"lat","fourier","complex"]
            else:
                ntimes = variable.shape[0]
                spvar = np.asfortranarray(
                           np.transpose(np.reshape(variable,
                                                   (ntimes,variable.shape[1])))
                           )
                fcshape = (ntimes,nlat,nlon//2,2)
                dims = ["time","lat","fourier","complex"]
            fcvar = pyfft.sp3fc(spvar,nlat,nlon,ntru,int(physfilter))
            fouriervar = np.reshape(np.transpose(fcvar),fcshape)
        else: #grid variable
            if len(variable.shape)==4: #include lev
                nlevs = variable.shape[1]
                ntimes = variable.shape[0]
                gpvar = np.asfortranarray(
                           np.transpose(np.reshape(variable/1.4142135623730951,
                                                   (ntimes*nlevs,nlat,nlon)))
                           )
                fcshape = (ntimes,nlevs,nlat,nlon//2,2)
                dims = ["time",levd,"lat","fourier","complex"]
            elif len(variable.shape)==2: #column
                fouriervar = variable
                dims = ["time",levd]
            elif len(variable.shape)==1: #scalar
                fouriervar = variable
                dims = ["time",]
            else:
                ntimes = variable.shape[0]
                gpvar = np.asfortranarray(
                           np.transpose(np.reshape(variable/1.4142135623730951,
                                                   (ntimes,nlat,nlon)))
                           )
                fcshape = (ntimes,nlat,nlon//2,2)
                dims = ["time","lat","fourier","complex"]
            if "fourier" in dims:
                fcvar = pyfft.gp3fc(gpvar)
                fouriervar = np.reshape(np.transpose(fcvar),fcshape)
        meta.append(tuple(dims))
        outvar = fouriervar
        
    elif mode=="syncfourier":
        #First transform into synchronous coordinate space
        if (ntru+1)*(ntru+2) in variable.shape: #spectral variable
            if len(variable.shape)==3: #Include lev
                nlevs = variable.shape[1]
                ntimes = variable.shape[0]
                spvar = np.asfortranarray(
                           np.transpose(np.reshape(variable,
                                                   (ntimes*nlevs,variable.shape[2])))
                           )
                gridshape = (ntimes,nlevs,nlat,nlon)
                dims = ["time",levd,"lat","lon"]
            else:
                ntimes = variable.shape[0]
                spvar = np.asfortranarray(
                           np.transpose(np.reshape(variable,
                                                   (ntimes,variable.shape[1])))
                           )
                gridshape = (ntimes,nlat,nlon)
                dims = ["time","lat","lon"]
            gridvar = pyfft.sp2gp(spvar,nlat,nlon,ntru,int(physfilter))
            gridvar = np.reshape(np.transpose(gridvar),gridshape)
            gridvar = np.roll(gridvar,nlon//2,axis=-1)
        else: #grid variable
            gridvar = variable
            if len(variable.shape)==3:
                dims = ["time","lat","lon"]
            elif len(variable.shape)==2: #column
                dims = ["time",levd]
            elif len(variable.shape)==1: #scalar
                dims = ["time",]
            else:
                dims = ["time",levd,"lat","lon"]
        if meta[0]=="hus":
            gridvar[gridvar<0] = 0.0
        if zonal and "lon" in dims:
            gridvar = np.nanmean(gridvar,axis=-1)
            dims.remove("lon")
        
        if not presync and "lat" in dims:
            lon,lat,tlgridvar = gcmt.eq2tl(gridvar,lon,lat,substellar=substellarlon,
                                           polemethod='interp') #fine bc all vectors are derived
        else:
            tlgridvar = gridvar
        #Next we want to reshape things so that parallel meridians link up to form circles.
        #This will leave us with half the longitudes, and twice the latitudes.
        
        #0 degrees "latitude" fourier coefficients will be the fourier transform along the 
        #substellar polar meridian, while 90 degrees "latitude" will be the transform along
        #the equatorial meridian.

        if "lat" in dims:
            rottlgridvar = np.zeros(tlgridvar.shape)
            for jlon in range(nlats):
                rottlgridvar[...,jlon,:nlats] = tlgridvar[...,jlon]
                rottlgridvar[...,jlon,nlats:] = tlgridvar[...,jlon+nlats]    
            
            #Compute fourier coefficients along our new "longitudes"
            if len(rottlgridvar.shape)==4: #include lev
                nlevs = rottlgridvar.shape[1]
                ntimes = rottlgridvar.shape[0]
                gpvar = np.asfortranarray(
                        np.transpose(np.reshape(rottlgridvar/1.4142135623730951,
                                                (ntimes*nlevs,nlat,nlon)))
                        )
                fcshape = (ntimes,nlevs,nlat,nlon//2,2)
                dims = ["time",levd,"lat","fourier","complex"]
            else:
                ntimes = rottlgridvar.shape[0]
                gpvar = np.asfortranarray(
                        np.transpose(np.reshape(rottlgridvar/1.4142135623730951,
                                                (ntimes,nlat,nlon)))
                        )
                fcshape = (ntimes,nlat,nlon//2,2)
                dims = ["time","lat","fourier","complex"]
            fcvar = pyfft.gp3fc(gpvar)
            fouriervar = np.reshape(np.transpose(fcvar),fcshape)
        else:
            fouriervar = tlgridvar
        meta.append(tuple(dims))
        outvar = fouriervar
        
    else:
        raise Exception("Invalid output mode selected")
                
    return (outvar,meta)
    

def _transformvectorvar(lon,uvar,vvar,umeta,vmeta,lats,nlon,nlev,ntru,ntime,mode='grid',
                        substellarlon=180.0,physfilter=False,zonal=False,radius=6371220.0):
    '''Ensure a variable is in a given horizontal mode.
    
    Parameters
    ----------
    lon : array-like
        Longitude array, in degrees
    uvar : array-like
        U-axis data array from dataset (e.g. divergence, or u-wind)
    vvar : array-like
        V-axis data array from dataset (e.g. vorticity. or v-wind)
    meta : list
        Meta array from dataset
    lats : array-like 
        Latitude array
    nlon : int 
        Number of longitudes
    nlev : int 
        Number of vertical levels
    ntru : int 
        Truncation wavenumber
    ntime : int 
        Number of output frames
    mode : str, optional
        Horizontal output mode. Can be 'grid', meaning the Gaussian latitude-longitude grid used
        in ExoPlaSim, 'spectral', meaning spherical harmonics, 
        'fourier', meaning Fourier coefficients and latitudes, 'synchronous', meaning a
        Gaussian latitude-longitude grid in the synchronous coordinate system defined in
        Paradise, et al (2021), with the north pole centered on the substellar point, or
        'syncfourier', meaning Fourier coefficients computed along the dipolar meridians in the
        synchronous coordinate system (e.g. the substellar-antistellar-polar meridian, which is 0 degrees,
        or the substellar-evening-antistellar-morning equatorial meridian, which is 90 degrees). Because this
        will get assigned to the original latitude array, that will become -90 degrees for the polar
        meridian, and 0 degrees for the equatorial meridian, identical to the typical equatorial coordinate
        system.
    zonal : bool, optional
        For grid modes ("grid" and "synchronous"), compute and output zonal means
    substellarlon : float, optional
        If mode='synchronous', the longitude of the substellar point in equatorial coordinates,
        in degrees
    physfilter : bool, optional
        Whether or not a physics filter should be used when transforming spectral variables to
        Fourier or grid domains
        
    Returns
    -------
    numpy.ndarray
        Transformed array
    '''
    
    if np.nanmax(lats)>10: #Dealing with degrees, not radians
        rlats = lats*np.pi/180.0
    else:
        rlats = lats[:]
        
    nlat = len(rlats)
    
    if sys.version[0]=="2":
        if nlat in [192,320]:
            import exoplasim.pyfft991v2 as pyfft
        else:
            import exoplasim.pyfft2 as pyfft
    else:
        if nlat in [192,320]:
            import exoplasim.pyfft991 as pyfft
        else:
            import exoplasim.pyfft as pyfft
    
        
    rdcostheta = radius/np.cos(rlats)
    costhetadr = np.cos(rlats)/radius
    
    if nlev in uvar.shape:
        levd = "lev"
    elif nlev+1 in uvar.shape:
        levd = "levp"
    
    if mode=="grid":
        if (ntru+1)*(ntru+2) in uvar.shape: #spectral variable
            if len(uvar.shape)==3: #Include lev
                nlevs = uvar.shape[1]
                ntimes = uvar.shape[0]
                spuvar = np.asfortranarray(
                            np.transpose(np.reshape(uvar,
                                                    (ntimes*nlevs,uvar.shape[2])))
                           )
                spvvar = np.asfortranarray(
                            np.transpose(np.reshape(vvar,
                                                    (ntimes*nlevs,vvar.shape[2])))
                           )
                gridshape = (ntimes,nlevs,nlat,nlon)
                dims = ["time",levd,"lat","lon"]
            else:
                ntimes = uvar.shape[0]
                spuvar = np.asfortranarray(
                            np.transpose(np.reshape(uvar,
                                                    (ntimes,variable.shape[1])))
                           )
                spvvar = np.asfortranarray(
                            np.transpose(np.reshape(vvar,
                                                    (ntimes,variable.shape[1])))
                           )
                gridshape = (ntimes,nlat,nlon)
                dims = ["time","lat","lon"]
            griduvar,gridvvar = pyfft.spvgp(spuvar,spvvar,rdcostheta,nlon,ntru,int(physfilter))
            griduvar = np.reshape(np.transpose(griduvar),gridshape)
            gridvvar = np.reshape(np.transpose(gridvvar),gridshape)
            griduvar = np.roll(griduvar,nlon//2,axis=-1)
            gridvvar = np.roll(gridvvar,nlon//2,axis=-1)
        else: #grid variable
            griduvar = uvar
            gridvvar = vvar
            if len(uvar.shape)==3:
                dims = ["time","lat","lon"]
            else:
                dims = ["time",levd,"lat","lon"]
        if zonal:
            griduvar = np.nanmean(griduvar,axis=-1)
            gridvvar = np.nanmean(gridvvar,axis=-1)
            dims.remove("lon")
        umeta.append(tuple(dims))
        vmeta.append(tuple(dims))
        outuvar = griduvar
        outvvar = gridvvar
        
    elif mode=="synchronous":
        if (ntru+1)*(ntru+2) in uvar.shape: #spectral variable
            if len(uvar.shape)==3: #Include lev
                nlevs = uvar.shape[1]
                ntimes = uvar.shape[0]
                spuvar = np.asfortranarray(
                            np.transpose(np.reshape(uvar,
                                                    (ntimes*nlevs,uvar.shape[2])))
                           )
                spvvar = np.asfortranarray(
                            np.transpose(np.reshape(vvar,
                                                    (ntimes*nlevs,vvar.shape[2])))
                           )
                gridshape = (ntimes,nlevs,nlat,nlon)
                dims = ["time",levd,"lat","lon"]
            else:
                ntimes = uvar.shape[0]
                spuvar = np.asfortranarray(
                            np.transpose(np.reshape(uvar,
                                                    (ntimes,uvar.shape[1])))
                           )
                spvvar = np.asfortranarray(
                            np.transpose(np.reshape(vvar,
                                                    (ntimes,vvar.shape[1])))
                           )
                gridshape = (ntimes,nlat,nlon)
                dims = ["time","lat","lon"]
            griduvar, gridvvar = pyfft.spvgp(spuvar,spvvar,rdcostheta,nlon,ntru,int(physfilter))
            griduvar = np.reshape(np.transpose(griduvar),gridshape)
            gridvvar = np.reshape(np.transpose(gridvvar),gridshape)
            griduvar = np.roll(griduvar,nlon//2,axis=-1)
            gridvvar = np.roll(gridvvar,nlon//2,axis=-1)
        else: #grid uvar
            griduvar = uvar
            gridvvar = vvar
            if len(uvar.shape)==3:
                dims = ["time","lat","lon"]
            else:
                dims = ["time",levd,"lat","lon"]
        lon,lat,tlgriduvar,tlgridvvar = gcmt.eq2tl_uv(griduvar,gridvvar,lon,lats,
                                                      substellar=substellarlon)
        
        if zonal:
            tlgriduvar = np.nanmean(tlgriduvar,axis=-1)
            tlgridvvar = np.nanmean(tlgridvvar,axis=-1)
            dims.remove("lon")
        umeta.append(tuple(dims))
        vmeta.append(tuple(dims))
        outuvar = tlgriduvar
        outvvar = tlgridvvar
        
    elif mode=="spectral":
        if (ntru+1)*(ntru+2) in uvar.shape: #spectral variable
            specuvar = uvar
            specvvar = vvar
            if len(uvar.shape)==3:
                dims = ("time",levd,"modes","complex")
            else:
                dims = ("time","modes","complex")
        else:
            if len(uvar.shape)==4: #Include lev
                nlevs = uvar.shape[1]
                ntimes = uvar.shape[0]
                gpuvar = np.asfortranarray(
                            np.transpose(np.reshape(uvar,
                                                    (ntimes*nlevs,nlat,nlon)))
                           )
                gpvvar = np.asfortranarray(
                            np.transpose(np.reshape(vvar,
                                                    (ntimes*nlevs,nlat,nlon)))
                           )
                dims = ("time",levd,"modes","complex")
            else:
                
                ntimes = uvar.shape[0]
                gpuvar = np.asfortranarray(
                            np.transpose(np.reshape(uvar,
                                                    (ntimes,nlat,nlon)))
                           )
                gpvvar = np.asfortranarray(
                            np.transpose(np.reshape(vvar,
                                                    (ntimes,nlat,nlon)))
                           )
                dims = ("time","modes","complex")
            spuvar,spvvar = pyfft.gpvsp(gpvuar,gpvvar,costhetadr,ntru,int(physfilter))
            specuvar = np.transpose(spuvar)
            specvvar = np.transpose(spvvar)
        shape = list(specuvar.shape)
        shape[-1] //= 2
        shape.append(2)
        shape = tuple(shape)
        specuvar = np.reshape(specuvar,shape)
        specvvar = np.reshape(specvvar,shape)
        
        umeta.append(dims)
        vmeta.append(dims)
        outuvar = specuvar
        outvvar = specvvar
        
    elif mode=="fourier":
        if (ntru+1)*(ntru+2) in uvar.shape: #spectral variable
            if len(uvar.shape)==3: #Include lev
                nlevs = uvar.shape[1]
                ntimes = uvar.shape[0]
                spuvar = np.asfortranarray(
                            np.transpose(np.reshape(uvar,
                                                    (ntimes*nlevs,uvar.shape[2])))
                           )
                spvvar = np.asfortranarray(
                            np.transpose(np.reshape(vvar,
                                                    (ntimes*nlevs,vvar.shape[2])))
                           )
                fcshape = (ntimes,nlevs,nlat,nlon//2,2)
                dims = ["time",levd,"lat","fourier","complex"]
            else:
                ntimes = uvar.shape[0]
                spuvar = np.asfortranarray(
                            np.transpose(np.reshape(uvar,
                                                    (ntimes,uvar.shape[1])))
                           )
                spvvar = np.asfortranarray(
                            np.transpose(np.reshape(vvar,
                                                    (ntimes,vvar.shape[1])))
                           )
                fcshape = (ntimes,nlat,nlon//2,2)
                dims = ["time","lat","fourier","complex"]
            gpuvar, gpvvar = pyfft.spvgp(spuvar,spvvar,rdcostheta,nlon,ntru,int(physfilter))
            fcuvar = pyfft.gp3fc(gpuvar)
            fcvvar = pyfft.gp3fc(gpvvar)
            fourieruvar = np.reshape(np.transpose(fcuvar),fcshape)
            fouriervvar = np.reshape(np.transpose(fcvvar),fcshape)
        else: #grid variable
            if len(uvar.shape)==4: #include lev
                nlevs = uvar.shape[1]
                ntimes = uvar.shape[0]
                gpuvar = np.asfortranarray(
                            np.transpose(np.reshape(uvar/1.4142135623730951,
                                                    (ntimes*nlevs,nlat,nlon)))
                           )
                gpvvar = np.asfortranarray(
                            np.transpose(np.reshape(vvar/1.4142135623730951,
                                                    (ntimes*nlevs,nlat,nlon)))
                           )
                fcshape = (ntimes,nlevs,nlat,nlon//2,2)
                dims = ["time",levd,"lat","fourier","complex"]
            else:
                ntimes = uvar.shape[0]
                gpuvar = np.asfortranarray(
                            np.transpose(np.reshape(uvar/1.4142135623730951,
                                                    (ntimes,nlat,nlon)))
                           )
                gpvvar = np.asfortranarray(
                            np.transpose(np.reshape(vvar/1.4142135623730951,
                                                    (ntimes,nlat,nlon)))
                           )
                fcshape = (ntimes,nlat,nlon//2,2)
                dims = ["time","lat","fourier","complex"]
            fcuvar = pyfft.gp3fc(gpuvar)
            fcvvar = pyfft.gp3fc(gpvvar)
            fourieruvar = np.reshape(np.transpose(fcuvar),fcshape)
            fouriervvar = np.reshape(np.transpose(fcvvar),fcshape)
        umeta.append(tuple(dims))
        vmeta.append(tuple(dims))
        outuvar = fourieruvar
        outvvar = fouriervvar
        
    elif mode=="syncfourier":
        #First transform into synchronous coordinate space
        if (ntru+1)*(ntru+2) in uvar.shape: #spectral variable
            if len(uvar.shape)==3: #Include lev
                nlevs = uvar.shape[1]
                ntimes = uvar.shape[0]
                spuvar = np.asfortranarray(
                            np.transpose(np.reshape(uvar,
                                                    (ntimes*nlevs,uvar.shape[2])))
                           )
                spvvar = np.asfortranarray(
                            np.transpose(np.reshape(vvar,
                                                    (ntimes*nlevs,vvar.shape[2])))
                           )
                gridshape = (ntimes,nlevs,nlat,nlon)
                dims = ["time",levd,"lat","lon"]
            else:
                ntimes = uvar.shape[0]
                spuvar = np.asfortranarray(
                            np.transpose(np.reshape(uvar,
                                                    (ntimes,uvar.shape[1])))
                           )
                spvvar = np.asfortranarray(
                            np.transpose(np.reshape(vvar,
                                                    (ntimes,vvar.shape[1])))
                           )
                gridshape = (ntimes,nlat,nlon)
                dims = ["time","lat","lon"]
            griduvar, gridvvar = pyfft.spvgp(spuvar,spvvar,rdcostheta,nlon,ntru,int(physfilter))
            griduvar = np.reshape(np.transpose(griduvar),gridshape)
            gridvvar = np.reshape(np.transpose(gridvvar),gridshape)
            griduvar = np.roll(griduvar,nlon//2,axis=-1)
            gridvvar = np.roll(gridvvar,nlon//2,axis=-1)
        else: #grid uvar
            griduvar = uvar
            gridvvar = vvar
            if len(uvar.shape)==3:
                dims = ["time","lat","lon"]
            else:
                dims = ["time",levd,"lat","lon"]
        lon,lat,tlgriduvar,tlgridvvar = gcmt.eq2tl_uv(griduvar,gridvvar,lon,lats,
                                                      substellar=substellarlon)
        #Next we want to reshape things so that parallel meridians link up to form circles.
        #This will leave us with half the longitudes, and twice the latitudes.
        
        #0 degrees "latitude" fourier coefficients will be the fourier transform along the 
        #substellar polar meridian, while 90 degrees "latitude" will be the transform along
        #the equatorial meridian.

        rottlgriduvar = np.zeros(tlgriduvar.shape)
        rottlgridvvar = np.zeros(tlgridvvar.shape)
        for jlon in range(nlats):
            rottlgriduvar[...,jlon,:nlats] = tlgriduvar[...,jlon]
            rottlgridvvar[...,jlon,nlats:] = tlgridvvar[...,jlon+nlats]    
        
        #Compute fourier coefficients along our new "longitudes"
        if len(rottlgriduvar.shape)==4: #include lev
            nlevs = rottlgriduvar.shape[1]
            ntimes = rottlgriduvar.shape[0]
            gpuvar = np.asfortranarray(
                        np.transpose(np.reshape(rottlgriduvar/1.4142135623730951,
                                                (ntimes*nlevs,nlat,nlon)))
                       )
            gpvvar = np.asfortranarray(
                        np.transpose(np.reshape(rottlgridvvar/1.4142135623730951,
                                                (ntimes*nlevs,nlat,nlon)))
                       )
            fcshape = (ntimes,nlevs,nlat,nlon//2,2)
            dims = ["time",levd,"lat","fourier","complex"]
        else:
            ntimes = rottlgridvar.shape[0]
            gpuvar = np.asfortranarray(
                        np.transpose(np.reshape(rottlgriduvar/1.4142135623730951,
                                                (ntimes,nlat,nlon)))
                       )
            gpvvar = np.asfortranarray(
                        np.transpose(np.reshape(rottlgridvvar/1.4142135623730951,
                                                (ntimes,nlat,nlon)))
                       )
            fcshape = (ntimes,nlat,nlon//2,2)
            dims = ["time","lat","fourier","complex"]
        fcuvar = pyfft.gp3fc(gpuvar)
        fcvvar = pyfft.gp3fc(gpvvar)
        fourieruvar = np.reshape(np.transpose(fcuvar),fcshape)
        fouriervvar = np.reshape(np.transpose(fcvvar),fcshape)
        umeta.append(tuple(dims))
        vmeta.append(tuple(dims))
        outuvar = fourieruvar
        outvvar = fouriervvar
        
    else:
        raise Exception("Invalid output mode selected")
                
    return (outuvar,outvvar,umeta,vmeta)
    

def dataset(filename, variablecodes, mode='grid', zonal=False, substellarlon=180.0, physfilter=False,
            radius=1.0,gravity=9.80665,gascon=287.0,logfile=None):
    '''Read a raw output file, and construct a dataset.
    
    Parameters
    ----------
    filename : str
        Path to the raw output file
    variablecodes : array-like
        list of variables to include. Can be the integer variable codes from the burn7 postprocessor
        conventions (as either strings or integers), or the short variable name strings 
        (e.g. 'rlut'), or a combination of the two.
    mode : str, optional
        Horizontal output mode. Can be 'grid', meaning the Gaussian latitude-longitude grid used
        in ExoPlaSim, 'spectral', meaning spherical harmonics, 
        'fourier', meaning Fourier coefficients and latitudes, 'synchronous', meaning a
        Gaussian latitude-longitude grid in the synchronous coordinate system defined in
        Paradise, et al (2021), with the north pole centered on the substellar point, or
        'syncfourier', meaning Fourier coefficients computed along the dipolar meridians in the
        synchronous coordinate system (e.g. the substellar-antistellar-polar meridian, which is 0 degrees,
        or the substellar-evening-antistellar-morning equatorial meridian, which is 90 degrees). Because this
        will get assigned to the original latitude array, that will become -90 degrees for the polar
        meridian, and 0 degrees for the equatorial meridian, identical to the typical equatorial coordinate
        system.
    zonal : bool, optional
        For grid modes ("grid" and "synchronous"), compute and output zonal means
    substellarlon : float, optional
        If mode='synchronous', the longitude of the substellar point in equatorial coordinates,
        in degrees
    physfilter : bool, optional
        Whether or not a physics filter should be used when transforming spectral variables to
        Fourier or grid domains
    radius : float, optional
        Planet radius in Earth radii
    gravity : float, optional
        Surface gravity in m/s^2.
    gascon : float, optional
        Specific gas constant for dry gas (R$_d$) in J/kg/K.  
    logfile : str or None, optional
        If None, log diagnostics will get printed to standard output. Otherwise, the log file
        to which diagnostic output should be written.
        
    Returns
    -------
    dict
        Dictionary of extracted variables
    '''
    
    plarad = radius*6371220.0 #convert Earth radii to metres
    
    rawdata = readfile(filename)
    
    
    lat = rawdata["lat"]
    lon = rawdata["lon"]
    lev = rawdata["lev"]
    time = rawdata["time"]
    
    nlat = len(lat)
    nlon = len(lon)
    nlev = len(lev)
    ntime = len(time)
    
    ntru = (nlon-1) // 3
    
    rdataset = {}
    
    windless=True #Once we get ua and va, this should be false
    
    rlat = lat*np.pi/180.0
    rlon = lon*np.pi/180.0
    colat = np.cos(rlat)
    
    gridlnps,lnpsmeta = _transformvar(lon[:],lat[:],rawdata[str(lnpscode)][:],ilibrary[str(lnpscode)][:],nlat,nlon,
                                      nlev,ntru,ntime,mode='grid',substellarlon=substellarlon,
                                      physfilter=physfilter,zonal=False)
    dpsdx = np.zeros(gridlnps.shape)
    for jlat in range(nlat):
        dpsdx[...,jlat,:] = np.gradient(gridlnps[...,jlat,:],rlon*plarad*colat[jlat],axis=-1)
    dpsdy = np.gradient(gridlnps,rlat*plarad,axis=-2)
    gridps = np.exp(gridlnps)
    
    levp = np.zeros(nlev+1)
    levp[-1] = 1.0
    levp[1:-1] = 0.5*(lev[1:]+lev[0:-1])
    levp[0] = 0.5*lev[0]#-(levp[1]-lev[0])
    pa = gridps[:,np.newaxis,:,:] * lev[np.newaxis,:,np.newaxis,np.newaxis]
    hpa = gridps[:,np.newaxis,:,:] * levp[np.newaxis,:,np.newaxis,np.newaxis]
    
    meanpa = np.nanmean(pa,axis=(0,2,3))*1.0e-2
    meanhpa = np.nanmean(hpa,axis=(0,2,3))*1.0e-2
    
    _log(logfile,"Interface Pressure     Mid-Level Pressure")
    _log(logfile,"*****************************************") #%18s
    _log(logfile,"%14f hpa     ------------------"%(meanhpa[0]))
    for jlev in range(nlev):
        _log(logfile,"------------------     %14f hpa"%(meanpa[jlev]))
        _log(logfile,"%14f hpa     ------------------"%(meanhpa[jlev+1]))
    _log(logfile,"*****************************************") #%18s

    specmodes = np.zeros((ntru+1)*(ntru+2))
    
    w=0
    for m in range(ntru+1):
        for n in range(m,ntru+1):
            specmodes[w  ] = n
            specmodes[w+1] = n
            w+=2
    
    for key in variablecodes:
        '''Collect metadata from our built-in list, and extract 
        the variable data if it already exists; if not set a flag
        that we need to derive it.'''
        if type(key)==int:
            try:
                meta = ilibrary[str(key)][:]
            except:
                raise Exception("Unknown variable code requested: %s"%str(key))
            if str(key) in rawdata:
                variable = rawdata[str(key)][:]
                derived=False
            else:
                #_log(logfile,str(key)+" not in rawdata;",rawdata.keys())
                derived=True
        else:
            if key in ilibrary:
                meta = ilibrary[key][:]
            elif key in slibrary:
                meta = slibrary[key][:]
                kcode = meta[0]
                meta[0] = key
                #_log(logfile,"Reassigning key; key was %s and is now %s"%(key,kcode))
                key = str(kcode) #Now key is always the integer code, and meta[0] is always the name
            else:
                raise Exception("Unknown variable code requested: %s"%key)
            if key in rawdata:
                variable = rawdata[key][:]
                derived=False
            else:
                #_log(logfile,key+" not in rawdata;",rawdata.keys())
                derived=True
        #_log(logfile,meta)
        meta.append(key)
        #_log(logfile,meta)
        #_log(logfile,meta,derived,rawdata.keys())
        if not derived:
            #_log(logfile,"Found variable; no need to derive: %s"%meta[0])
            variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
            rdataset[meta[0]]= [variable,meta]
            _log(logfile,"Collected variable: %8s\t.... %3d timestamps"%(meta[0],variable.shape[0]))
        else: #derived=True       
        
            # Add in derived variables
            #_log(logfile,nlat,nlon,ntru,nlevs*ntimes)
            #ta = pyfft.sp2gp(np.asfortranarray(np.transpose(np.reshape(data["130"],(ntimes*nlevs,data["130"].shape[-1])))),
                             #nlat,nlon,ntru,int(physfilter))
            #ta = np.reshape(np.transpose(ta),(ntimes,nlevs,nlat,nlon))
            
            #data["ta"] = ta
            
            if key==str(ucode): #ua
                if windless:
                    div  = rawdata[str(divcode)][:]
                    vort = rawdata[str(vortcode)][:]
                    umeta = ilibrary[key][:]
                    umeta.append(key)
                    vmeta = ilibrary[str(vcode)][:]
                    ua,va,meta,vmeta = _transformvectorvar(lon[:],div,vort,umeta,vmeta,lat,nlon,nlev,ntru,
                                                            ntime,mode=mode,substellarlon=substellarlon,
                                                            physfilter=physfilter,zonal=zonal,
                                                            radius=plarad)
                    windless = False
                else:
                    meta=umeta[:-1]
                    meta.append(key)
                    meta.append(umeta[-1])
                rdataset[meta[0]] = [ua,meta]
            elif key==str(vcode): #va
                if windless:
                    div  = rawdata[str(divcode)][:]
                    vort = rawdata[str(vortcode)][:]
                    umeta = ilibrary[str(ucode)][:]
                    vmeta = ilibrary[key][:]
                    vmeta.append(key)
                    ua,va,umeta,meta = _transformvectorvar(lon[:],div,vort,umeta,vmeta,lat,nlon,nlev,ntru,
                                                            ntime,mode=mode,substellarlon=substellarlon,
                                                            physfilter=physfilter,zonal=zonal,
                                                            radius=plarad)
                    windless = False
                else:
                    meta=vmeta[:-1]
                    meta.append(key)
                    meta.append(vmeta[-1])
                rdataset[meta[0]] = [va,meta]
            elif key==str(spdcode): #spd
                if windless:
                    div  = rawdata[str(divcode)][:]
                    vort = rawdata[str(vortcode)][:]
                    umeta = ilibrary[str(ucode)][:]
                    vmeta = ilibrary[str(vcode)][:]
                    ua,va,umeta,vmeta = _transformvectorvar(lon[:],div,vort,umeta,vmeta,lat,nlon,nlev,ntru,
                                                            ntime,mode=mode,substellarlon=substellarlon,
                                                            physfilter=physfilter,zonal=zonal,
                                                            radius=plarad)
                    windless = False
                meta = ilibrary[key][:]
                meta.append(key)
                meta.append(umeta[-1])
                spd = np.sqrt(ua**2+va**2)
                rdataset[meta[0]] = [spd,meta]
                
            elif key==str(dpsdxcode): #dpsdx
                meta = ilibrary[key][:]
                meta.append(key)
                variable = gridps*dpsdx
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(dpsdycode): #dpsdy
                meta = ilibrary[key][:]
                meta.append(key)
                variable = gridps*dpsdy
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(preccode): #precipiation
                # prc + prl
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["142"][:]+rawdata["143"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(ntopcode): #Net top radiation
                # rst + rlut
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["178"][:]+rawdata["179"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(nbotcode): #Net bottom radiation
                # rss + rls
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["176"][:]+rawdata["177"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(nheatcode): #Net heat flux
                # Melt*L*rho + rss + rls + hfss + hfls
                meta = ilibrary[key][:]
                meta.append(key)
                variable = (rawdata["218"][:]*L_TIMES_RHOH2O +rawdata["176"][:] + rawdata["177"][:]
                           +rawdata["146"][:] + rawdata["147"][:])
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(nh2ocode): #Net water flux
                # evap - mrro + precip
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["182"][:] - rawdata["160"][:] + rawdata["142"][:] + rawdata["143"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(swatmcode): #Shortwave net 
                # rst = rss
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["178"][:] - rawdata["176"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(lwatmcode): #longwave net 
                # rlut - rst
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["179"][:] - rawdata["177"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(natmcode): #Net atmospheric radiation
                # rst + rlut - rss - rst
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["178"][:] - rawdata["176"][:] + rawdata["179"][:] - rawdata["177"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(sruncode): #Precip + Evap - Increase in snow  = water added to bucket
                #Actual runoff should be precip + evap + melt + soilh2o - bucketmax
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["182"][:] - rawdata["221"][:] + rawdata["142"][:] + rawdata["143"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(freshcode): #Precip + Evap
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["142"][:] + rawdata["143"][:] + rawdata["182"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(wcode): #Omega? vertical air velocity in Pa/s
                # w = p(j)*(u(i,j)*dpsdx(i,j)+v(i,j)*dpsdy(i,j)) 
                  #   - deltap(j)*(div(i,j)+u(i,j)*dpsdx(i,j)+v(i,j)*dpsdy(i,j))
                  
                
                #get pa
                #if "ps" in rdataset:
                    #ps,pmeta = _transformvar(lon[:],lat[:],rdataset["ps"][0][:],rdataset["ps"][1][:],nlat,nlon,nlev,
                                             #ntru,ntime,mode='grid',substellarlon=substellarlon,
                                             #physfilter=physfilter,zonal=False)
                #else:
                    #ps,pmeta = _transformvar(lon[:],lat[:],gridps,ilibrary["134"],nlat,nlon,nlev,
                                             #ntru,ntime,mode='grid',substellarlon=substellarlon,
                                             #physfilter=physfilter,zonal=False)
                #pa = ps[:,np.newaxis,:,:]*lev[np.newaxis,:,np.newaxis,np.newaxis]
                #We already got pa
                  
                if not windless:
                    uu,vv,umeta,vmeta = _transformvectorvar(lon[:],div,vort,umeta[:],vmeta[:],lat,
                                                            nlon,nlev,ntru,
                                                            ntime,mode='grid',radius=plarad,
                                                            substellarlon=substellarlon,
                                                            physfilter=physfilter,zonal=False)
                    dv,dmeta = _transformvar(lon[:],lat[:],div,ilibrary[str(divcode)][:],nlat,nlon,nlev,ntru,ntime,
                                             mode='grid',substellarlon=substellarlon,
                                             physfilter=physfilter,zonal=False)
                else:
                    div  = rawdata[str(divcode)][:]
                    vort = rawdata[str(vortcode)][:]
                    umeta = ilibrary[str(ucode)][:]
                    vmeta = ilibrary[str(vcode)][:]
                    uu,vv,umeta,vmeta = _transformvectorvar(lon[:],div,vort,umeta,vmeta,lat,nlon,nlev,ntru,
                                                            ntime,mode='grid',radius=plarad,
                                                            substellarlon=substellarlon,
                                                            physfilter=physfilter,zonal=False)
                    dv,dmeta = _transformvar(lon[:],lat[:],div,ilibrary[str(divcode)][:],nlat,nlon,nlev,ntru,ntime,
                                             mode='grid',substellarlon=substellarlon,
                                             physfilter=physfilter,zonal=False)
                    
                wap = np.zeros(dv.shape)
                for t in range(ntime):
                    for j in range(nlat):
                        for i in range(nlon):
                            wap[t,:,j,i] = (pa[t,:,j,i]*(uu[t,:,j,i]*dpsdx[t,j,i] 
                                                        +vv[t,:,j,i]*dpsdy[t,j,i]) 
                                            - scipy.integrate.cumtrapz(np.append([0,],
                                                                       dv[t,:,j,i]
                                                                       +uu[t,:,j,i]*dpsdx[t,j,i]
                                                                       +vv[t,:,j,i]*dpsdy[t,j,i]),
                                                                       x=np.append([0,],pa[t,:,j,i])))
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],wap,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]] = [variable*1.0e-2,meta] #transform to hPa/s
                
            elif key==str(wzcode): #Vertical wind wa
                # wa = -omega * gascon * ta / (grav * pa)
                
                #get pa
                #if "ps" in rdataset:
                    #ps,pmeta = _transformvar(lon[:],lat[:],rdataset["ps"][0][:],rdataset["ps"][1][:],nlat,nlon,nlev,
                                             #ntru,ntime,mode='grid',substellarlon=substellarlon,
                                             #physfilter=physfilter,zonal=False)
                #else:
                    #ps,pmeta = _transformvar(lon[:],lat[:],gridps,ilibrary["134"],nlat,nlon,nlev,
                                             #ntru,ntime,mode='grid',substellarlon=substellarlon,
                                             #physfilter=physfilter,zonal=False)
                #pa = ps[:,np.newaxis,:,:]*lev[np.newaxis,:,np.newaxis,np.newaxis]
                #We already got pa
                
                if "wap" in rdataset:
                    omega = rdataset["wap"][0]
                else:
                    if not windless:
                        uu,vv,umeta,vmeta = _transformvectorvar(lon[:],div,vort,umeta[:],vmeta[:],
                                                                lat,nlon,nlev,ntru,
                                                                ntime,mode='grid',radius=plarad,
                                                                substellarlon=substellarlon,
                                                                physfilter=physfilter,zonal=False)
                        dv,dmeta = _transformvar(lon[:],lat[:],div,ilibrary[str(divcode)][:],nlat,nlon,nlev,ntru,ntime,
                                                  mode='grid',substellarlon=substellarlon,
                                                  physfilter=physfilter,zonal=False)
                    else:
                        div  = rawdata[str(divcode)][:]
                        vort = rawdata[str(vortcode)][:]
                        umeta = ilibrary[str(ucode)][:]
                        vmeta = ilibrary[str(vcode)][:]
                        uu,vv,umeta,vmeta = _transformvectorvar(lon[:],div,vort,umeta,vmeta,lat,nlon,nlev,ntru,
                                                                ntime,mode='grid',radius=plarad,
                                                                substellarlon=substellarlon,
                                                                physfilter=physfilter,zonal=False)
                        dv,dmeta = _transformvar(lon[:],lat[:],div,ilibrary[str(divcode)][:],nlat,nlon,nlev,ntru,ntime,
                                                  mode='grid',substellarlon=substellarlon,
                                                  physfilter=physfilter,zonal=False)
                    omega = np.zeros(dv.shape)
                    for t in range(ntime):
                        for j in range(nlat):
                            for i in range(nlon):
                                omega[t,:,j,i] = (pa[t,:,j,i]*(uu[t,:,j,i]*dpsdx[t,j,i] 
                                                            +vv[t,:,j,i]*dpsdy[t,j,i]) 
                                                - scipy.integrate.cumtrapz(np.append([0,],
                                                                        dv[t,:,j,i]
                                                                        +uu[t,:,j,i]*dpsdx[t,j,i]
                                                                        +vv[t,:,j,i]*dpsdy[t,j,i]),
                                                                    x=np.append([0,],(pa[t,:,j,i]))))
                omega,wmeta = _transformvar(lon[:],lat[:],omega,vmeta,nlat,nlon,nlev,ntru,ntime,mode='grid',
                                            substellarlon=substellarlon,physfilter=physfilter)
                if "ta" in rdataset:
                    ta,tameta = _transformvar(lon[:],lat[:],rdataset["ta"][0][:],ilibrary[str(tempcode)][:],nlat,nlon,
                                              nlev,ntru,ntime,mode='grid',substellarlon=substellarlon,
                                              physfilter=physfilter,zonal=False)
                else:
                    ta,tameta = _transformvar(lon[:],lat[:],rawdata[str(tempcode)][:],ilibrary[str(tempcode)][:],
                                              nlat,nlon,nlev,ntru,ntime,mode='grid',
                                              substellarlon=substellarlon,
                                              physfilter=physfilter,zonal=False)
                
                
                wa = -omega*gascon*ta / (gravity*pa)
                meta = ilibrary[key][:]
                meta.append(key)
                wa,meta = _transformvar(lon[:],lat[:],wa,meta,nlat,nlon,
                                        nlev,ntru,ntime,mode=mode,
                                        substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]] = [wa,meta]
                
            elif key==str(pscode): #surface pressure (hPa)
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],gridps*1.0e-2,meta,nlat,nlon,
                                              nlev,ntru,ntime,
                                              mode=mode,substellarlon=substellarlon,
                                              physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(vpotcode): #Velocity potential (psi)
                meta = ilibrary[key][:]
                meta.append(key)
                vdiv,vmeta = _transformvar(lon[:],lat[:],rawdata[str(divcode)][:],
                                           meta,nlat,nlon,
                                           nlev,ntru,ntime,mode='spectral',substellarlon=substellarlon,
                                           physfilter=physfilter,zonal=False) #Need it to be spectral
                vdivshape = list(vdiv.shape[:-1])
                vdivshape[-1]*=2
                vdivshape = tuple(vdivshape)
                vdiv = np.reshape(vdiv,vdivshape)
                vpot = np.zeros(vdiv.shape)
                modes = np.resize(specmodes,vdiv.shape)
                vpot[...,2:] = vdiv[...,2:] * plarad**2/(modes**2+modes)[...,2:]

                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],vpot,meta,
                                              nlat,nlon,nlev,ntru,ntime,mode=mode,
                                    substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(stfcode): #Streamfunction (stf)
                if mode in ("fourier","spectral","syncfourier"):
                    if mode in ("fourier","spectral"):
                        tempmode = "grid"
                    else:
                        tempmode = "synchronous"
                else:
                    tempmode = mode
                if windless:
                    div  = rawdata[str(divcode)][:]
                    vort = rawdata[str(vortcode)][:]
                    umeta = ilibrary[str(ucode)][:]
                    vmeta = ilibrary[str(vcode)][:]
                    ua,va,umeta,vmeta = _transformvectorvar(lon[:],div,vort,umeta,vmeta,
                                                            lat,nlon,nlev,ntru,
                                                            ntime,mode=tempmode,
                                                            substellarlon=substellarlon,
                                                            physfilter=physfilter,zonal=False,
                                                            radius=plarad)
                    windless = False
                else:
                    ua,va,umeta2,vmeta2 = _transformvectorvar(lon[:],div,vort,umeta,vmeta,
                                                              lat,nlon,nlev,ntru,
                                                              ntime,mode=tempmode,
                                                              substellarlon=substellarlon,
                                                              physfilter=physfilter,zonal=False,
                                                              radius=plarad)
                                                              
                #svort,smeta = _transformvar(lon[:],lat[:],rawdata[str(vortcode)][:],
                                            #ilibrary[str(vortcode)][:],nlat,
                                            #nlon,nlev,ntru,ntime,mode='spectral',
                                            #substellarlon=substellarlon,physfilter=physfilter,
                                            #zonal=False) #Need it to be spectral
                #svortshape = list(svort.shape[:-1])
                #svortshape[-1]*=2
                #svortshape = tuple(svortshape)
                #svort = np.reshape(svort,svortshape)
                #stf = np.zeros(svort.shape)
                #modes = np.resize(specmodes,svort.shape)
                #stf[...,2:] = svort[...,2:] * plarad**2/(modes**2+modes)[...,2:]
                
                vadp = np.zeros(va.shape)
                for nt in range(ntime):
                    for jlat in range(nlat):
                        for jlon in range(nlon):
                            vadp[nt,:,jlat,jlon] = scipy.integrate.cumtrapz(va[nt,:,jlat,jlon],
                                                                           x=pa[nt,:,jlat,jlon],
                                                                           initial=0.0)
                    
                prefactor = 2*np.pi*plarad*colat/gravity
                sign = 1 - 2*(tempmode=="synchronous") #-1 for synchronous, 1 for equatorial
                stf = sign*prefactor[np.newaxis,np.newaxis,:,np.newaxis]*vadp
            
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],stf,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal,presync=True)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(slpcode): #Sea-level pressure (slp)
                
                if "sg" in rdataset:
                    geopot,gmeta = _transformvar(lon[:],lat[:],rdataset["sg"][0][:],
                                                 ilibrary[str(geopotcode)][:],nlat,
                                                 nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                else:
                    geopot,gmeta = _transformvar(lon[:],lat[:],rawdata[str(geopotcode)][:],
                                                 ilibrary[str(geopotcode)][:],
                                                 nlat,nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                
                #temp should be bottom layer of atmospheric temperature
                
                if "ta" in rdataset:
                    tta,tmeta = _transformvar(lon[:],lat[:],rdataset["ta"][0][:],ilibrary[str(tempcode)][:],nlat,nlon,
                                              nlev,ntru,ntime,mode="grid",substellarlon=substellarlon,
                                              physfilter=physfilter,zonal=False)
                else:
                    tta,tmeta = _transformvar(lon[:],lat[:],rawdata[str(tempcode)][:],ilibrary[str(tempcode)][:],nlat,
                                              nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                temp = tta[:,-1,...]
                
                #aph is half-level pressure
                #apf is full-level pressure
                aph = hpa[:,-1,...] #surface pressure
                apf =  pa[:,-1,...] #mid-layer pressure of bottom layer
                
                slp = np.zeros(geopot.shape)
                slp[abs(geopot)<1.0e-4] = aph[abs(geopot)<1.0e-4]
                
                mask = abs(geopot)>=1.0e-4
                alpha = gascon*RLAPSE/gravity
                tstar = (1 + alpha*(aph[mask]/apf[mask]-1))*temp[mask]
                tstar[tstar<255.0] = 0.5*(255+tstar[tstar<255.0])
                tmsl = tstar + geopot[mask]*RLAPSE/gravity
                ZPRT = geopot[mask] / (gascon*tstar)
                ZPRTAL = np.zeros(ZPRT.shape)
                mask2 = abs(tmsl-tstar)<1.0e-6
                mask3 = abs(tmsl-tstar)>=1.0e-6
                ZPRTAL[mask2] = 0.0
                alpha = gascon * (tmsl[mask3]-tstar[mask3])/geopot[mask][mask3]
                ZPRTAL[mask3] = ZPRT[mask3] * alpha
                slp[mask] = aph[mask] * np.exp(ZPRT*(1.0-ZPRTAL*(0.5-ZPRTAL/3.0)))
                
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],slp,meta,nlat,nlon,
                                              nlev,ntru,ntime,mode=mode,
                                    substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                                              
            elif key==str(geopotzcode): #Geopotential height 
                #we need temperature, humidity, half-level pressure
                if "hus" in rdataset:
                    qq,qmeta = _transformvar(lon[:],lat[:],rdataset["hus"][0][:],
                                             ilibrary[str(humcode)][:],nlat,nlon,
                                             nlev,ntru,ntime,mode="grid",substellarlon=substellarlon,
                                             physfilter=physfilter,zonal=False)
                else:
                    qq,qmeta = _transformvar(lon[:],lat[:],rawdata[str(humcode)][:],ilibrary[str(humcode)][:],nlat,
                                             nlon,nlev,ntru,ntime,mode="grid",
                                             substellarlon=substellarlon,
                                             physfilter=physfilter,zonal=False)
                qq[qq<0] = 0.0
                
                
                if "ta" in rdataset:
                    temp,tmeta = _transformvar(lon[:],lat[:],rdataset["ta"][0][:],ilibrary[str(tempcode)][:],nlat,nlon,
                                              nlev,ntru,ntime,mode="grid",substellarlon=substellarlon,
                                              physfilter=physfilter,zonal=False)
                else:
                    temp,tmeta = _transformvar(lon[:],lat[:],rawdata[str(tempcode)][:],ilibrary[str(tempcode)][:],nlat,
                                               nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                
                if "sg" in rdataset:
                    oro,gmeta = _transformvar(lon[:],lat[:],rdataset["sg"][0][:],ilibrary[str(geopotcode)][:],nlat,
                                                 nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                else:
                    oro,gmeta = _transformvar(lon[:],lat[:],rawdata[str(geopotcode)][:],ilibrary[str(geopotcode)][:],
                                                 nlat,nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                
                gzshape = list(qq.shape)
                gzshape[1] = len(levp)
                gz = np.zeros(gzshape)
                
                gz[:,nlev,...] = oro[:] #bottom layer of geopotential Z is the orographic geopotential
                
                VTMP = RH2O/gascon - 1.0
                twolog2 = 2.0*np.log(2.0)
                
                if np.nanmax(qq)>=1.0e-14: #Non-dry atmosphere
                    for jlev in range(nlev-1,0,-1):
                        gz[:,jlev,...] = (gz[:,jlev+1,...]
                                        + gascon*temp[:,jlev,...]*(1.0+VTMP+qq[:,jlev,...])
                                                *np.log(hpa[:,jlev+1,...])/hpa[:,jlev,...])
                    gz[:,0,...] = gz[:,1,...] + gascon*temp[:,0,...]*(1.0+VTMP+qq[:,0,...])*twolog2
                    
                else: #Dry atmosphere
                    for jlev in range(nlev-1,0,-1):
                        gz[:,jlev,...] = (gz[:,jlev+1,...] + gascon*temp[:,jlev,...]
                                                             *np.log(hpa[:,jlev+1,...])/hpa[:,jlev,...])
                    gz[:,0,...] = gz[:,1,...] + gascon*temp[:,0,...]*twolog2
                
                gz *= 1.0/gravity
                
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],gz,meta,
                                              nlat,nlon,nlev,ntru,ntime,mode=mode,
                                    substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]

                
            elif key==str(rhumcode): #relative humidity (hur)
                
                rv     = 461.51
                TMELT  = 273.16
                ra1    = 610.78
                ra2    =  17.2693882
                ra4    =  35.86
                rdbrv  = gascon / rv
                
                if "ta" in rdataset:
                    temp,tmeta = _transformvar(lon[:],lat[:],rdataset["ta"][0][:],ilibrary[str(tempcode)][:],nlat,nlon,
                                              nlev,ntru,ntime,mode="grid",substellarlon=substellarlon,
                                              physfilter=physfilter,zonal=False)
                else:
                    temp,tmeta = _transformvar(lon[:],lat[:],rawdata[str(tempcode)][:],ilibrary[str(tempcode)][:],nlat,
                                               nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                
                if "hus" in rdataset:
                    qq,qmeta = _transformvar(lon[:],lat[:],rdataset["hus"][0][:],ilibrary[str(humcode)][:],nlat,nlon,
                                             nlev,ntru,ntime,mode="grid",substellarlon=substellarlon,
                                             physfilter=physfilter,zonal=False)
                else:
                    qq,qmeta = _transformvar(lon[:],lat[:],rawdata[str(humcode)][:],ilibrary[str(humcode)][:],nlat,
                                             nlon,nlev,ntru,ntime,mode="grid",
                                             substellarlon=substellarlon,
                                             physfilter=physfilter,zonal=False)
                
                #This is the saturation vapor pressure divided by the local pressure to give saturation
                #specific humidity, but it seems like it must account for the pressure contribution of
                #water.
                zqsat  = rdbrv * ra1 * np.exp(ra2 * (temp-TMELT)/(temp-ra4)) / pa #saturation spec hum
                zqsat *= 1.0 / (1.0 - (1.0/rdbrv-1.0)*zqsat)
                
                rh     = qq/zqsat * 100.0
                
                rh[rh<0.0  ] =   0.0
                rh[rh>100.0] = 100.0
                
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],rh,meta,nlat,nlon,nlev,
                                              ntru,ntime,mode=mode,
                                    substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(hpresscode): #Half-level pressure
                
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],hpa,meta,nlat,nlon,
                                              nlev,ntru,ntime,mode=mode,
                                           substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(fpresscode): #Full-level pressure
                
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],pa,meta,
                                              nlat,nlon,nlev,ntru,ntime,mode=mode,
                                          substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(thetahcode) or key==str(thetafcode): #Potential temperature
                
                if "theta" in rdataset or "thetah" in rdataset:
                    if key==str(thetahcode):
                        meta = ilibrary[key][:]
                        meta.append(key)
                        variable,meta = _transformvar(lon[:],lat[:],thetah,meta,nlat,nlon,nlev,ntru,
                                                      ntime,mode=mode,substellarlon=substellarlon,
                                                      physfilter=physfilter,zonal=zonal)
                        rdataset[meta[0]]= [variable,meta]
                    elif key==str(thetafcode):
                        variable,meta = _transformvar(lon[:],lat[:],theta,meta,nlat,nlon,nlev,ntru,
                                                      ntime,mode=mode,substellarlon=substellarlon,
                                                      physfilter=physfilter,zonal=zonal)
                        rdataset[meta[0]]= [variable,meta]
                else:
                    
                    if "ta" in rdataset:
                        ta,tmeta = _transformvar(lon[:],lat[:],rdataset["ta"][0][:],ilibrary[str(tempcode)][:],nlat,
                                                 nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,
                                                 physfilter=physfilter,zonal=False)
                    else:
                        ta,tmeta = _transformvar(lon[:],lat[:],rawdata[str(tempcode)][:],ilibrary[str(tempcode)][:],
                                                 nlat,nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                    
                    if "ts" in rdataset:
                        tsurf,tsmeta = _transformvar(lon[:],lat[:],rdataset["ts"][0][:],ilibrary[str(tscode)][:],nlat,
                                                     nlon,nlev,ntru,ntime,mode="grid",
                                                     substellarlon=substellarlon,
                                                     physfilter=physfilter,zonal=False)
                    else:
                        tsurf,tsmeta = _transformvar(lon[:],lat[:],rawdata[str(tscode)][:],ilibrary[str(tscode)][:],
                                                     nlat,nlon,nlev,ntru,ntime,mode="grid",
                                                     substellarlon=substellarlon,physfilter=physfilter,
                                                     zonal=False)
                    
                    thetah = np.zeros(hpa.shape)
                    theta  = np.zeros( pa.shape)
                    
                    kappa = 1.0/3.5
                    
                    for jlev in range(nlev-1):
                        thetah[:,jlev+1,...] = (0.5*(ta[:,jlev,...]+ta[:,jlev+1,...])
                                               *(gridps/hpa[:,jlev,...])**kappa)
                    thetah[:,nlev,...] = tsurf[:]
                    theta = 0.5*(thetah[:,:-1,...] + thetah[:,1:,...])
                    
                    meta = ilibrary[key][:]
                    meta.append(key)
                    if key==str(thetahcode):
                        variable,meta = _transformvar(lon[:],lat[:],thetah,meta,
                                                      nlat,nlon,nlev,ntru,
                                                      ntime,mode=mode,substellarlon=substellarlon,
                                                      physfilter=physfilter,zonal=zonal)
                        rdataset[meta[0]]= [variable,meta]
                    elif key==str(thetafcode):
                        variable,meta = _transformvar(lon[:],lat[:],theta,meta,
                                                      nlat,nlon,nlev,ntru,
                                                      ntime,mode=mode,substellarlon=substellarlon,
                                                      physfilter=physfilter,zonal=zonal)
                        rdataset[meta[0]]= [variable,meta]
                        
            _log(logfile,"Collected variable: %8s\t.... %3d timestamps"%(meta[0],variable.shape[0]))
          
    rdataset["lat"] = [np.array(lat),["lat","latitude","deg"] ]
    rdataset["lon"] = [np.array(lon),["lon","longitude","deg"]]
    rdataset["lev"] = [np.array(lev),["lev","sigma_coordinate","nondimensional"]       ]
    rdataset["levp"] = [np.array(levp),["levp","half_sigma_coordinate","nondimensional"]]
    rdataset["time"] = [np.array(time),["time","timestep_of_year","timesteps"]         ]      
    
    return rdataset


def advancedDataset(filename, variablecodes, substellarlon=180.0,
                    radius=1.0,gravity=9.80665,gascon=287.0,logfile=None):
    '''Read a raw output file, and construct a dataset.
    
    Parameters
    ----------
    filename : str
        Path to the raw output file
    variablecodes : dict
        Variables to include. Each member must use the variable name as the key, and contain a sub-dict
        with the horizontal mode, zonal averaging, and physics filtering  options optionall set as 
        members. For example::
            {"ts":{"mode":"grid","zonal":False},
             "stf":{"mode":"grid","zonal":True,"physfilter":True}}
        Options that are not set take on their
        default values from :py:func:`dataset() <exoplasim.pyburn.dataset>`.
    mode : str, optional
        Horizontal output mode. Can be 'grid', meaning the Gaussian latitude-longitude grid used
        in ExoPlaSim, 'spectral', meaning spherical harmonics, 
        'fourier', meaning Fourier coefficients and latitudes, 'synchronous', meaning a
        Gaussian latitude-longitude grid in the synchronous coordinate system defined in
        Paradise, et al (2021), with the north pole centered on the substellar point, or
        'syncfourier', meaning Fourier coefficients computed along the dipolar meridians in the
        synchronous coordinate system (e.g. the substellar-antistellar-polar meridian, which is 0 degrees,
        or the substellar-evening-antistellar-morning equatorial meridian, which is 90 degrees). Because this
        will get assigned to the original latitude array, that will become -90 degrees for the polar
        meridian, and 0 degrees for the equatorial meridian, identical to the typical equatorial coordinate
        system.
    zonal : bool, optional
        For grid modes ("grid" and "synchronous"), compute and output zonal means
    physfilter : bool, optional
        Whether or not a physics filter should be used when transforming spectral variables to
        Fourier or grid domains
    substellarlon : float, optional
        If mode='synchronous', the longitude of the substellar point in equatorial coordinates,
        in degrees
    radius : float, optional
        Planet radius in Earth radii
    gravity : float, optional
        Surface gravity in m/s^2.
    gascon : float, optional
        Specific gas constant for dry gas (R$_d$) in J/kg/K. 
    logfile : str or None, optional
        If None, log diagnostics will get printed to standard output. Otherwise, the log file
        to which diagnostic output should be written.
        
    Returns
    -------
    dict
        Dictionary of extracted variables
    '''
    
    plarad = radius*6371220.0 #convert Earth radii to metres
    
    rawdata = readfile(filename)
    
    
    lat = rawdata["lat"]
    lon = rawdata["lon"]
    lev = rawdata["lev"]
    time = rawdata["time"]
    
    nlat = len(lat)
    nlon = len(lon)
    nlev = len(lev)
    ntime = len(time)
    
    ntru = (nlon-1) // 3
    
    rdataset = {}
    
    windless=True #Once we get ua and va, this should be false
    
    rlat = lat*np.pi/180.0
    rlon = lon*np.pi/180.0
    colat = np.cos(rlat)
    
    gridlnps,lnpsmeta = _transformvar(lon[:],lat[:],rawdata[str(lnpscode)][:],ilibrary[str(lnpscode)][:],nlat,nlon,
                                      nlev,ntru,ntime,mode='grid',substellarlon=substellarlon,
                                      physfilter=physfilter,zonal=False)
    dpsdx = np.zeros(gridlnps.shape)
    for jlat in range(nlat):
        dpsdx[...,jlat,:] = np.gradient(gridlnps[...,jlat,:],rlon*plarad*colat[jlat],axis=-1)
    dpsdy = np.gradient(gridlnps,rlat*plarad,axis=-2)
    gridps = np.exp(gridlnps)
    
    levp = np.zeros(nlev+1)
    levp[-1] = 1.0
    levp[1:-1] = 0.5*(lev[1:]+lev[0:-1])
    levp[0] = 0.5*lev[0]#-(levp[1]-lev[0])
    pa = gridps[:,np.newaxis,:,:] * lev[np.newaxis,:,np.newaxis,np.newaxis]
    hpa = gridps[:,np.newaxis,:,:] * levp[np.newaxis,:,np.newaxis,np.newaxis]
    
    meanpa = np.nanmean(pa,axis=(0,2,3))*1.0e-2
    meanhpa = np.nanmean(hpa,axis=(0,2,3))*1.0e-2
    
    _log(logfile,"Interface Pressure     Mid-Level Pressure")
    _log(logfile,"*****************************************") #%18s
    _log(logfile,"%14f hpa     ------------------"%(meanhpa[0]))
    for jlev in range(nlev):
        _log(logfile,"------------------     %14f hpa"%(meanpa[jlev]))
        _log(logfile,"%14f hpa     ------------------"%(meanhpa[jlev+1]))
    _log(logfile,"*****************************************") #%18s

    specmodes = np.zeros((ntru+1)*(ntru+2))
    
    w=0
    for m in range(ntru+1):
        for n in range(m,ntru+1):
            specmodes[w  ] = n
            specmodes[w+1] = n
            w+=2
    
    for key in variablecodes:
        '''Collect metadata from our built-in list, and extract 
        the variable data if it already exists; if not set a flag
        that we need to derive it.'''
        if type(key)==int:
            try:
                meta = ilibrary[str(key)][:]
            except:
                raise Exception("Unknown variable code requested: %s"%str(key))
            if str(key) in rawdata:
                variable = rawdata[str(key)][:]
                derived=False
            else:
                #_log(logfile,str(key)+" not in rawdata;",rawdata.keys())
                derived=True
        else:
            if key in ilibrary:
                meta = ilibrary[key][:]
            elif key in slibrary:
                meta = slibrary[key][:]
                kcode = meta[0]
                meta[0] = key
                #_log(logfile,"Reassigning key; key was %s and is now %s"%(key,kcode))
                key = str(kcode) #Now key is always the integer code, and meta[0] is always the name
            else:
                raise Exception("Unknown variable code requested: %s"%key)
            if key in rawdata:
                variable = rawdata[key][:]
                derived=False
            else:
                #_log(logfile,key+" not in rawdata;",rawdata.keys())
                derived=True
        #_log(logfile,meta)
        meta.append(key)
        #_log(logfile,meta)
        #_log(logfile,meta,derived,rawdata.keys())
        if not derived:
            #_log(logfile,"Found variable; no need to derive: %s"%meta[0])
            mode = "grid"; zonal=False; physfilter=False
            if "mode" in variablecodes[key]:
                mode=variablecodes[key]["mode"]
            if "zonal" in variablecodes[key]:
                zonal=variablecodes[key]["zonal"]
            if "physfilter" in variablecodes[key]:
                physfilter=variablecodes[key]["physfilter"]
            variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
            rdataset[meta[0]]= [variable,meta]
            _log(logfile,"Collected variable: %8s\t.... %3d timestamps"%(meta[0],variable.shape[0]))
        else: #derived=True       
        
            # Add in derived variables
            #_log(logfile,nlat,nlon,ntru,nlevs*ntimes)
            #ta = pyfft.sp2gp(np.asfortranarray(np.transpose(np.reshape(data["130"],(ntimes*nlevs,data["130"].shape[-1])))),
                             #nlat,nlon,ntru,int(physfilter))
            #ta = np.reshape(np.transpose(ta),(ntimes,nlevs,nlat,nlon))
            
            #data["ta"] = ta
            
            mode = "grid"; zonal=False; physfilter=False
            if "mode" in variablecodes[key]:
                mode=variablecodes[key]["mode"]
            if "zonal" in variablecodes[key]:
                zonal=variablecodes[key]["zonal"]
            if "physfilter" in variablecodes[key]:
                physfilter=variablecodes[key]["physfilter"]
            
            if key==str(ucode): #ua
                if windless:
                    div  = rawdata[str(divcode)][:]
                    vort = rawdata[str(vortcode)][:]
                    umeta = ilibrary[key][:]
                    umeta.append(key)
                    vmeta = ilibrary[str(vcode)][:]
                    ua,va,meta,vmeta = _transformvectorvar(lon[:],div,vort,umeta,vmeta,lat,nlon,nlev,ntru,
                                                            ntime,mode=mode,substellarlon=substellarlon,
                                                            physfilter=physfilter,zonal=zonal,
                                                            radius=plarad)
                    windless = False
                else:
                    meta=umeta[:-1]
                    meta.append(key)
                    meta.append(umeta[-1])
                rdataset[meta[0]] = [ua,meta]
            elif key==str(vcode): #va
                if windless:
                    div  = rawdata[str(divcode)][:]
                    vort = rawdata[str(vortcode)][:]
                    umeta = ilibrary[str(ucode)][:]
                    vmeta = ilibrary[key][:]
                    vmeta.append(key)
                    ua,va,umeta,meta = _transformvectorvar(lon[:],div,vort,umeta,vmeta,lat,nlon,nlev,ntru,
                                                            ntime,mode=mode,substellarlon=substellarlon,
                                                            physfilter=physfilter,zonal=zonal,
                                                            radius=plarad)
                    windless = False
                else:
                    meta=vmeta[:-1]
                    meta.append(key)
                    meta.append(vmeta[-1])
                rdataset[meta[0]] = [va,meta]
            elif key==str(spdcode): #spd
                if windless:
                    div  = rawdata[str(divcode)][:]
                    vort = rawdata[str(vortcode)][:]
                    umeta = ilibrary[str(ucode)][:]
                    vmeta = ilibrary[str(vcode)][:]
                    ua,va,umeta,vmeta = _transformvectorvar(lon[:],div,vort,umeta,vmeta,lat,nlon,nlev,ntru,
                                                            ntime,mode=mode,substellarlon=substellarlon,
                                                            physfilter=physfilter,zonal=zonal,
                                                            radius=plarad)
                    windless = False
                meta = ilibrary[key][:]
                meta.append(key)
                meta.append(umeta[-1])
                spd = np.sqrt(ua**2+va**2)
                rdataset[meta[0]] = [spd,meta]
                
            elif key==str(dpsdxcode): #dpsdx
                meta = ilibrary[key][:]
                meta.append(key)
                variable = gridps*dpsdx
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(dpsdycode): #dpsdy
                meta = ilibrary[key][:]
                meta.append(key)
                variable = gridps*dpsdy
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(preccode): #precipiation
                # prc + prl
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["142"][:]+rawdata["143"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(ntopcode): #Net top radiation
                # rst + rlut
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["178"][:]+rawdata["179"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(nbotcode): #Net bottom radiation
                # rss + rls
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["176"][:]+rawdata["177"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(nheatcode): #Net heat flux
                # Melt*L*rho + rss + rls + hfss + hfls
                meta = ilibrary[key][:]
                meta.append(key)
                variable = (rawdata["218"][:]*L_TIMES_RHOH2O +rawdata["176"][:] + rawdata["177"][:]
                           +rawdata["146"][:] + rawdata["147"][:])
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(nh2ocode): #Net water flux
                # evap - mrro + precip
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["182"][:] - rawdata["160"][:] + rawdata["142"][:] + rawdata["143"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(swatmcode): #Shortwave net 
                # rst = rss
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["178"][:] - rawdata["176"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(lwatmcode): #longwave net 
                # rlut - rst
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["179"][:] - rawdata["177"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(natmcode): #Net atmospheric radiation
                # rst + rlut - rss - rst
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["178"][:] - rawdata["176"][:] + rawdata["179"][:] - rawdata["177"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(sruncode): #Precip + Evap - Increase in snow  = water added to bucket
                #Actual runoff should be precip + evap + melt + soilh2o - bucketmax
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["182"][:] - rawdata["221"][:] + rawdata["142"][:] + rawdata["143"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(freshcode): #Precip + Evap
                meta = ilibrary[key][:]
                meta.append(key)
                variable = rawdata["142"][:] + rawdata["143"][:] + rawdata["182"][:]
                variable,meta = _transformvar(lon[:],lat[:],variable,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(wcode): #Omega? vertical air velocity in Pa/s
                # w = p(j)*(u(i,j)*dpsdx(i,j)+v(i,j)*dpsdy(i,j)) 
                  #   - deltap(j)*(div(i,j)+u(i,j)*dpsdx(i,j)+v(i,j)*dpsdy(i,j))
                  
                
                #get pa
                #if "ps" in rdataset:
                    #ps,pmeta = _transformvar(lon[:],lat[:],rdataset["ps"][0][:],rdataset["ps"][1][:],nlat,nlon,nlev,
                                             #ntru,ntime,mode='grid',substellarlon=substellarlon,
                                             #physfilter=physfilter,zonal=False)
                #else:
                    #ps,pmeta = _transformvar(lon[:],lat[:],gridps,ilibrary["134"],nlat,nlon,nlev,
                                             #ntru,ntime,mode='grid',substellarlon=substellarlon,
                                             #physfilter=physfilter,zonal=False)
                #pa = ps[:,np.newaxis,:,:]*lev[np.newaxis,:,np.newaxis,np.newaxis]
                #We already got pa
                  
                if not windless:
                    uu,vv,umeta,vmeta = _transformvectorvar(lon[:],div,vort,umeta[:],vmeta[:],lat,
                                                            nlon,nlev,ntru,
                                                            ntime,mode='grid',radius=plarad,
                                                            substellarlon=substellarlon,
                                                            physfilter=physfilter,zonal=False)
                    dv,dmeta = _transformvar(lon[:],lat[:],div,ilibrary[str(divcode)][:],nlat,nlon,nlev,ntru,ntime,
                                             mode='grid',substellarlon=substellarlon,
                                             physfilter=physfilter,zonal=False)
                else:
                    div  = rawdata[str(divcode)][:]
                    vort = rawdata[str(vortcode)][:]
                    umeta = ilibrary[str(ucode)][:]
                    vmeta = ilibrary[str(vcode)][:]
                    uu,vv,umeta,vmeta = _transformvectorvar(lon[:],div,vort,umeta,vmeta,lat,nlon,nlev,ntru,
                                                            ntime,mode='grid',radius=plarad,
                                                            substellarlon=substellarlon,
                                                            physfilter=physfilter,zonal=False)
                    dv,dmeta = _transformvar(lon[:],lat[:],div,ilibrary[str(divcode)][:],nlat,nlon,nlev,ntru,ntime,
                                             mode='grid',substellarlon=substellarlon,
                                             physfilter=physfilter,zonal=False)
                    
                wap = np.zeros(dv.shape)
                for t in range(ntime):
                    for j in range(nlat):
                        for i in range(nlon):
                            wap[t,:,j,i] = (pa[t,:,j,i]*(uu[t,:,j,i]*dpsdx[t,j,i] 
                                                        +vv[t,:,j,i]*dpsdy[t,j,i]) 
                                            - scipy.integrate.cumtrapz(np.append([0,],
                                                                       dv[t,:,j,i]
                                                                       +uu[t,:,j,i]*dpsdx[t,j,i]
                                                                       +vv[t,:,j,i]*dpsdy[t,j,i]),
                                                                       x=np.append([0,],pa[t,:,j,i])))
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],wap,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]] = [variable*1.0e-2,meta] #transform to hPa/s
                
            elif key==str(wzcode): #Vertical wind wa
                # wa = -omega * gascon * ta / (grav * pa)
                
                #get pa
                #if "ps" in rdataset:
                    #ps,pmeta = _transformvar(lon[:],lat[:],rdataset["ps"][0][:],rdataset["ps"][1][:],nlat,nlon,nlev,
                                             #ntru,ntime,mode='grid',substellarlon=substellarlon,
                                             #physfilter=physfilter,zonal=False)
                #else:
                    #ps,pmeta = _transformvar(lon[:],lat[:],gridps,ilibrary["134"],nlat,nlon,nlev,
                                             #ntru,ntime,mode='grid',substellarlon=substellarlon,
                                             #physfilter=physfilter,zonal=False)
                #pa = ps[:,np.newaxis,:,:]*lev[np.newaxis,:,np.newaxis,np.newaxis]
                #We already got pa
                
                if "wap" in rdataset:
                    omega = rdataset["wap"][0]
                else:
                    if not windless:
                        uu,vv,umeta,vmeta = _transformvectorvar(lon[:],div,vort,umeta[:],vmeta[:],
                                                                lat,nlon,nlev,ntru,
                                                                ntime,mode='grid',radius=plarad,
                                                                substellarlon=substellarlon,
                                                                physfilter=physfilter,zonal=False)
                        dv,dmeta = _transformvar(lon[:],lat[:],div,ilibrary[str(divcode)][:],nlat,nlon,nlev,ntru,ntime,
                                                  mode='grid',substellarlon=substellarlon,
                                                  physfilter=physfilter,zonal=False)
                    else:
                        div  = rawdata[str(divcode)][:]
                        vort = rawdata[str(vortcode)][:]
                        umeta = ilibrary[str(ucode)][:]
                        vmeta = ilibrary[str(vcode)][:]
                        uu,vv,umeta,vmeta = _transformvectorvar(lon[:],div,vort,umeta,vmeta,lat,nlon,nlev,ntru,
                                                                ntime,mode='grid',radius=plarad,
                                                                substellarlon=substellarlon,
                                                                physfilter=physfilter,zonal=False)
                        dv,dmeta = _transformvar(lon[:],lat[:],div,ilibrary[str(divcode)][:],nlat,nlon,nlev,ntru,ntime,
                                                  mode='grid',substellarlon=substellarlon,
                                                  physfilter=physfilter,zonal=False)
                        omega = np.zeros(dv.shape)
                        for t in range(ntime):
                            for j in range(nlat):
                                for i in range(nlon):
                                    omega[t,:,j,i] = (pa[t,:,j,i]*(uu[t,:,j,i]*dpsdx[t,j,i] 
                                                                  +vv[t,:,j,i]*dpsdy[t,j,i]) 
                                                      - scipy.integrate.cumtrapz(np.append([0,],
                                                                             dv[t,:,j,i]
                                                                             +uu[t,:,j,i]*dpsdx[t,j,i]
                                                                             +vv[t,:,j,i]*dpsdy[t,j,i]),
                                                                        x=np.append([0,],(pa[t,:,j,i]))))
                omega,wmeta = _transformvar(lon[:],lat[:],omega,vmeta,nlat,nlon,nlev,ntru,ntime,mode='grid',
                                            substellarlon=substellarlon,physfilter=physfilter)
                if "ta" in rdataset:
                    ta,tameta = _transformvar(lon[:],lat[:],rdataset["ta"][0][:],ilibrary[str(tempcode)][:],nlat,nlon,
                                              nlev,ntru,ntime,mode='grid',substellarlon=substellarlon,
                                              physfilter=physfilter,zonal=False)
                else:
                    ta,tameta = _transformvar(lon[:],lat[:],rawdata[str(tempcode)][:],ilibrary[str(tempcode)][:],
                                              nlat,nlon,nlev,ntru,ntime,mode='grid',
                                              substellarlon=substellarlon,
                                              physfilter=physfilter,zonal=False)
                
                
                wa = -omega*gascon*ta / (gravity*pa)
                meta = ilibrary[key][:]
                meta.append(key)
                wa,meta = _transformvar(lon[:],lat[:],wa,meta,nlat,nlon,
                                        nlev,ntru,ntime,mode=mode,
                                        substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]] = [wa,meta]
                
            elif key==str(pscode): #surface pressure
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],gridps*1.0e-2,meta,nlat,nlon,
                                              nlev,ntru,ntime,
                                              mode=mode,substellarlon=substellarlon,
                                              physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(vpotcode): #Velocity potential (psi)
                meta = ilibrary[key][:]
                meta.append(key)
                vdiv,vmeta = _transformvar(lon[:],lat[:],rawdata[str(divcode)][:],
                                           meta,nlat,nlon,
                                           nlev,ntru,ntime,mode='spectral',substellarlon=substellarlon,
                                           physfilter=physfilter,zonal=False) #Need it to be spectral
                vdivshape = list(vdiv.shape[:-1])
                vdivshape[-1]*=2
                vdivshape = tuple(vdivshape)
                vdiv = np.reshape(vdiv,vdivshape)
                vpot = np.zeros(vdiv.shape)
                modes = np.resize(specmodes,vdiv.shape)
                vpot[...,2:] = vdiv[...,2:] * plarad**2/(modes**2+modes)[...,2:]

                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],vpot,meta,
                                              nlat,nlon,nlev,ntru,ntime,mode=mode,
                                    substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(stfcode): #Streamfunction (stf)
                if mode in ("fourier","spectral","syncfourier"):
                    if mode in ("fourier","spectral"):
                        tempmode = "grid"
                    else:
                        tempmode = "synchronous"
                else:
                    tempmode = mode
                if windless:
                    div  = rawdata[str(divcode)][:]
                    vort = rawdata[str(vortcode)][:]
                    umeta = ilibrary[str(ucode)][:]
                    vmeta = ilibrary[str(vcode)][:]
                    ua,va,umeta,vmeta = _transformvectorvar(lon[:],div,vort,umeta,vmeta,
                                                            lat,nlon,nlev,ntru,
                                                            ntime,mode=tempmode,
                                                            substellarlon=substellarlon,
                                                            physfilter=physfilter,zonal=False,
                                                            radius=plarad)
                    windless = False
                else:
                    ua,va,umeta2,vmeta2 = _transformvectorvar(lon[:],div,vort,umeta,vmeta,
                                                              lat,nlon,nlev,ntru,
                                                              ntime,mode=tempmode,
                                                              substellarlon=substellarlon,
                                                              physfilter=physfilter,zonal=False,
                                                              radius=plarad)
                                                              
                #svort,smeta = _transformvar(lon[:],lat[:],rawdata[str(vortcode)][:],
                                            #ilibrary[str(vortcode)][:],nlat,
                                            #nlon,nlev,ntru,ntime,mode='spectral',
                                            #substellarlon=substellarlon,physfilter=physfilter,
                                            #zonal=False) #Need it to be spectral
                #svortshape = list(svort.shape[:-1])
                #svortshape[-1]*=2
                #svortshape = tuple(svortshape)
                #svort = np.reshape(svort,svortshape)
                #stf = np.zeros(svort.shape)
                #modes = np.resize(specmodes,svort.shape)
                #stf[...,2:] = svort[...,2:] * plarad**2/(modes**2+modes)[...,2:]
                
                vadp = np.zeros(va.shape)
                for nt in range(ntime):
                    for jlat in range(nlat):
                        for jlon in range(nlon):
                            vadp[nt,:,jlat,jlon] = scipy.integrate.cumtrapz(va[nt,:,jlat,jlon],
                                                                           x=pa[nt,:,jlat,jlon],
                                                                           initial=0.0)
                    
                prefactor = 2*np.pi*plarad*colat/gravity
                sign = 1 - 2*(tempmode=="synchronous") #-1 for synchronous, 1 for equatorial
                stf = sign*prefactor[np.newaxis,np.newaxis,:,np.newaxis]*vadp
                
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],stf,meta,nlat,nlon,nlev,ntru,ntime,mode=mode,
                                              substellarlon=substellarlon,physfilter=physfilter,
                                              zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(slpcode): #Sea-level pressure (slp)
                
                if "sg" in rdataset:
                    geopot,gmeta = _transformvar(lon[:],lat[:],rdataset["sg"][0][:],
                                                 ilibrary[str(geopotcode)][:],nlat,
                                                 nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                else:
                    geopot,gmeta = _transformvar(lon[:],lat[:],rawdata[str(geopotcode)][:],
                                                 ilibrary[str(geopotcode)][:],
                                                 nlat,nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                
                #temp should be bottom layer of atmospheric temperature
                
                if "ta" in rdataset:
                    tta,tmeta = _transformvar(lon[:],lat[:],rdataset["ta"][0][:],ilibrary[str(tempcode)][:],nlat,nlon,
                                              nlev,ntru,ntime,mode="grid",substellarlon=substellarlon,
                                              physfilter=physfilter,zonal=False)
                else:
                    tta,tmeta = _transformvar(lon[:],lat[:],rawdata[str(tempcode)][:],ilibrary[str(tempcode)][:],nlat,
                                              nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                temp = tta[:,-1,...]
                
                #aph is half-level pressure
                #apf is full-level pressure
                aph = hpa[:,-1,...] #surface pressure
                apf =  pa[:,-1,...] #mid-layer pressure of bottom layer
                
                slp = np.zeros(geopot.shape)
                slp[abs(geopot)<1.0e-4] = aph[abs(geopot)<1.0e-4]
                
                mask = abs(geopot)>=1.0e-4
                alpha = gascon*RLAPSE/gravity
                tstar = (1 + alpha*(aph[mask]/apf[mask]-1))*temp[mask]
                tstar[tstar<255.0] = 0.5*(255+tstar[tstar<255.0])
                tmsl = tstar + geopot[mask]*RLAPSE/gravity
                ZPRT = geopot[mask] / (gascon*tstar)
                ZPRTAL = np.zeros(ZPRT.shape)
                mask2 = abs(tmsl-tstar)<1.0e-6
                mask3 = abs(tmsl-tstar)>=1.0e-6
                ZPRTAL[mask2] = 0.0
                alpha = gascon * (tmsl[mask3]-tstar[mask3])/geopot[mask][mask3]
                ZPRTAL[mask3] = ZPRT[mask3] * alpha
                slp[mask] = aph[mask] * np.exp(ZPRT*(1.0-ZPRTAL*(0.5-ZPRTAL/3.0)))
                
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],slp,meta,nlat,nlon,
                                              nlev,ntru,ntime,mode=mode,
                                    substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                                              
            elif key==str(geopotzcode): #Geopotential height 
                #we need temperature, humidity, half-level pressure
                if "hus" in rdataset:
                    qq,qmeta = _transformvar(lon[:],lat[:],rdataset["hus"][0][:],
                                             ilibrary[str(humcode)][:],nlat,nlon,
                                             nlev,ntru,ntime,mode="grid",substellarlon=substellarlon,
                                             physfilter=physfilter,zonal=False)
                else:
                    qq,qmeta = _transformvar(lon[:],lat[:],rawdata[str(humcode)][:],ilibrary[str(humcode)][:],nlat,
                                             nlon,nlev,ntru,ntime,mode="grid",
                                             substellarlon=substellarlon,
                                             physfilter=physfilter,zonal=False)
                qq[qq<0] = 0.0
                
                
                if "ta" in rdataset:
                    temp,tmeta = _transformvar(lon[:],lat[:],rdataset["ta"][0][:],ilibrary[str(tempcode)][:],nlat,nlon,
                                              nlev,ntru,ntime,mode="grid",substellarlon=substellarlon,
                                              physfilter=physfilter,zonal=False)
                else:
                    temp,tmeta = _transformvar(lon[:],lat[:],rawdata[str(tempcode)][:],ilibrary[str(tempcode)][:],nlat,
                                               nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                
                if "sg" in rdataset:
                    oro,gmeta = _transformvar(lon[:],lat[:],rdataset["sg"][0][:],ilibrary[str(geopotcode)][:],nlat,
                                                 nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                else:
                    oro,gmeta = _transformvar(lon[:],lat[:],rawdata[str(geopotcode)][:],ilibrary[str(geopotcode)][:],
                                                 nlat,nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                
                gzshape = list(qq.shape)
                gzshape[1] = len(levp)
                gz = np.zeros(gzshape)
                
                gz[:,nlev,...] = oro[:] #bottom layer of geopotential Z is the orographic geopotential
                
                VTMP = RH2O/gascon - 1.0
                twolog2 = 2.0*np.log(2.0)
                
                if np.nanmax(qq)>=1.0e-14: #Non-dry atmosphere
                    for jlev in range(nlev-1,0,-1):
                        gz[:,jlev,...] = (gz[:,jlev+1,...]
                                        + gascon*temp[:,jlev,...]*(1.0+VTMP+qq[:,jlev,...])
                                                *np.log(hpa[:,jlev+1,...])/hpa[:,jlev,...])
                    gz[:,0,...] = gz[:,1,...] + gascon*temp[:,0,...]*(1.0+VTMP+qq[:,0,...])*twolog2
                    
                else: #Dry atmosphere
                    for jlev in range(nlev-1,0,-1):
                        gz[:,jlev,...] = (gz[:,jlev+1,...] + gascon*temp[:,jlev,...]
                                                             *np.log(hpa[:,jlev+1,...])/hpa[:,jlev,...])
                    gz[:,0,...] = gz[:,1,...] + gascon*temp[:,0,...]*twolog2
                
                gz *= 1.0/gravity
                
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],gz,meta,
                                              nlat,nlon,nlev,ntru,ntime,mode=mode,
                                    substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]

                
            elif key==str(rhumcode): #relative humidity (hur)
                
                rv     = 461.51
                TMELT  = 273.16
                ra1    = 610.78
                ra2    =  17.2693882
                ra4    =  35.86
                rdbrv  = gascon / rv
                
                if "ta" in rdataset:
                    temp,tmeta = _transformvar(lon[:],lat[:],rdataset["ta"][0][:],ilibrary[str(tempcode)][:],nlat,nlon,
                                              nlev,ntru,ntime,mode="grid",substellarlon=substellarlon,
                                              physfilter=physfilter,zonal=False)
                else:
                    temp,tmeta = _transformvar(lon[:],lat[:],rawdata[str(tempcode)][:],ilibrary[str(tempcode)][:],nlat,
                                               nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                
                if "hus" in rdataset:
                    qq,qmeta = _transformvar(lon[:],lat[:],rdataset["hus"][0][:],ilibrary[str(humcode)][:],nlat,nlon,
                                             nlev,ntru,ntime,mode="grid",substellarlon=substellarlon,
                                             physfilter=physfilter,zonal=False)
                else:
                    qq,qmeta = _transformvar(lon[:],lat[:],rawdata[str(humcode)][:],ilibrary[str(humcode)][:],nlat,
                                             nlon,nlev,ntru,ntime,mode="grid",
                                             substellarlon=substellarlon,
                                             physfilter=physfilter,zonal=False)
                
                #This is the saturation vapor pressure divided by the local pressure to give saturation
                #specific humidity, but it seems like it must account for the pressure contribution of
                #water.
                zqsat  = rdbrv * ra1 * np.exp(ra2 * (temp-TMELT)/(temp-ra4)) / pa #saturation spec hum
                zqsat *= 1.0 / (1.0 - (1.0/rdbrv-1.0)*zqsat)
                
                rh     = qq/zqsat * 100.0
                
                rh[rh<0.0  ] =   0.0
                rh[rh>100.0] = 100.0
                
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],rh,meta,nlat,nlon,nlev,
                                              ntru,ntime,mode=mode,
                                    substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(hpresscode): #Half-level pressure
                
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],hpa,meta,nlat,nlon,
                                              nlev,ntru,ntime,mode=mode,
                                           substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(fpresscode): #Full-level pressure
                
                meta = ilibrary[key][:]
                meta.append(key)
                variable,meta = _transformvar(lon[:],lat[:],pa,meta,
                                              nlat,nlon,nlev,ntru,ntime,mode=mode,
                                          substellarlon=substellarlon,physfilter=physfilter,zonal=zonal)
                rdataset[meta[0]]= [variable,meta]
                
            elif key==str(thetahcode) or key==str(thetafcode): #Potential temperature
                
                if "theta" in rdataset or "thetah" in rdataset:
                    if key==str(thetahcode):
                        meta = ilibrary[key][:]
                        meta.append(key)
                        variable,meta = _transformvar(lon[:],lat[:],thetah,meta,nlat,nlon,nlev,ntru,
                                                      ntime,mode=mode,substellarlon=substellarlon,
                                                      physfilter=physfilter,zonal=zonal)
                        rdataset[meta[0]]= [variable,meta]
                    elif key==str(thetafcode):
                        variable,meta = _transformvar(lon[:],lat[:],theta,meta,nlat,nlon,nlev,ntru,
                                                      ntime,mode=mode,substellarlon=substellarlon,
                                                      physfilter=physfilter,zonal=zonal)
                        rdataset[meta[0]]= [variable,meta]
                else:
                    
                    if "ta" in rdataset:
                        ta,tmeta = _transformvar(lon[:],lat[:],rdataset["ta"][0][:],ilibrary[str(tempcode)][:],nlat,
                                                 nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,
                                                 physfilter=physfilter,zonal=False)
                    else:
                        ta,tmeta = _transformvar(lon[:],lat[:],rawdata[str(tempcode)][:],ilibrary[str(tempcode)][:],
                                                 nlat,nlon,nlev,ntru,ntime,mode="grid",
                                                 substellarlon=substellarlon,physfilter=physfilter,
                                                 zonal=False)
                    
                    if "ts" in rdataset:
                        tsurf,tsmeta = _transformvar(lon[:],lat[:],rdataset["ts"][0][:],ilibrary[str(tscode)][:],nlat,
                                                     nlon,nlev,ntru,ntime,mode="grid",
                                                     substellarlon=substellarlon,
                                                     physfilter=physfilter,zonal=False)
                    else:
                        tsurf,tsmeta = _transformvar(lon[:],lat[:],rawdata[str(tscode)][:],ilibrary[str(tscode)][:],
                                                     nlat,nlon,nlev,ntru,ntime,mode="grid",
                                                     substellarlon=substellarlon,physfilter=physfilter,
                                                     zonal=False)
                    
                    thetah = np.zeros(hpa.shape)
                    theta  = np.zeros( pa.shape)
                    
                    kappa = 1.0/3.5
                    
                    for jlev in range(nlev-1):
                        thetah[:,jlev+1,...] = (0.5*(ta[:,jlev,...]+ta[:,jlev+1,...])
                                               *(gridps/hpa[:,jlev,...])**kappa)
                    thetah[:,nlev,...] = tsurf[:]
                    theta = 0.5*(thetah[:,:-1,...] + thetah[:,1:,...])
                    
                    meta = ilibrary[key][:]
                    meta.append(key)
                    if key==str(thetahcode):
                        variable,meta = _transformvar(lon[:],lat[:],thetah,meta,
                                                      nlat,nlon,nlev,ntru,
                                                      ntime,mode=mode,substellarlon=substellarlon,
                                                      physfilter=physfilter,zonal=zonal)
                        rdataset[meta[0]]= [variable,meta]
                    elif key==str(thetafcode):
                        variable,meta = _transformvar(lon[:],lat[:],theta,meta,
                                                      nlat,nlon,nlev,ntru,
                                                      ntime,mode=mode,substellarlon=substellarlon,
                                                      physfilter=physfilter,zonal=zonal)
                        rdataset[meta[0]]= [variable,meta]
                        
            _log(logfile,"Collected variable: %8s\t.... %3d timestamps"%(meta[0],variable.shape[0]))
          
    rdataset["lat"] = [np.array(lat),["lat","latitude","deg"] ]
    rdataset["lon"] = [np.array(lon),["lon","longitude","deg"]]
    rdataset["lev"] = [np.array(lev),["lev","sigma_coordinate","nondimensional"]       ]
    rdataset["levp"] = [np.array(levp),["levp","half_sigma_coordinate","nondimensional"]]
    rdataset["time"] = [np.array(time),["time","timestep_of_year","timesteps"]         ]      
    
    return rdataset

                
def netcdf(rdataset,filename="most_output.nc",append=False,logfile=None):
    '''Write a dataset to a netCDF file.
    
    Parameters
    ----------
    rdataset : dict
        A dictionary of outputs as generated from :py:func:`pyburn.dataset()<exoplasim.pyburn.dataset>`
    filename : str, optional
        Path to the output file that should be written.
    append : bool, optional
        Whether the file should be opened in "append" mode, or overwritten (default).
    logfile : str or None, optional
        If None, log diagnostics will get printed to standard output. Otherwise, the log file
        to which diagnostic output should be written.
        
    Returns
    -------
    object
        A netCDF object corresponding to the file that has been written.
    '''
    import netCDF4 as nc
    
    latitude  = rdataset["lat" ]
    longitude = rdataset["lon" ]
    level     = rdataset["lev" ]
    levelp    = rdataset["levp"]
    timestamp = rdataset["time"]
    
    nlats  = len( latitude[0])
    nlons  = len(longitude[0])
    nlevs  = len(    level[0])
    ntimes = len(timestamp[0])
    ntru = (nlons-1)//3
    nmodes = ((ntru+1)*(ntru+2))//2
    
    complex64   = np.dtype([("real",np.float32),("imag",np.float32)]) 
    
    
    if append:
        ncd = nc.Dataset(filename, "a", format="NETCDF4")
        
        complex64_t = ncd.cmptypes["complex64"]
        
        ntimes0 = len(ncd.variables["time"][:])
        t0 = ntimes0
        t1 = t0+ntimes
        times = ncd.variables['time']
        times[t0:t1] = np.array(timestamp[0]).astype("float32")
    else:
        ncd = nc.Dataset(filename, "w", format="NETCDF4")
        
        complex64_t = ncd.createCompoundType(complex64,"complex64")
        
        t0 = 0
        t1 = ntimes
        
        lat = ncd.createDimension("lat",   nlats)
        lon = ncd.createDimension("lon",   nlons)
        lev = ncd.createDimension("lev",   nlevs)
        levp= ncd.createDimension("levp",  nlevs+1)
        ttime = ncd.createDimension("time",None)
        #cmplx = ncd.createDimension("complex",2)
        fourier = ncd.createDimension("fourier",nlons//2)
        sphmods = ncd.createDimension("modes",nmodes)
        
        latitudes   = ncd.createVariable("lat", "f4",("lat", ),zlib=True,least_significant_digit=6)
        longitudes  = ncd.createVariable("lon", "f4",("lon", ),zlib=True,least_significant_digit=6)
        levels      = ncd.createVariable("lev", "f4",("lev", ),zlib=True,least_significant_digit=6)
        hlevels     = ncd.createVariable("levp","f4",("levp",),zlib=True,least_significant_digit=6)
        times       = ncd.createVariable("time","f4",("time",),zlib=True,least_significant_digit=6)
        #complexn    = ncd.createVariable("complex",complex64_t,("complex",),
                                        #zlib=True,least_significant_digit=6)
        fourierc    = ncd.createVariable("fharmonic","f4",("fourier",),zlib=True,least_significant_digit=6)
        spharmonics = ncd.createVariable("modes","f4",("modes",),zlib=True,least_significant_digit=6)
    
        ncd.set_auto_mask(False)
        latitudes.set_auto_mask(False)  
        longitudes.set_auto_mask(False)   
        levels.set_auto_mask(False)  
        hlevels.set_auto_mask(False)
        times.set_auto_mask(False)
        #complexn.set_auto_mask(False)
        fourierc.set_auto_mask(False)
        spharmonics.set_auto_mask(False)
        
        latitudes.units  =  latitude[1][2]
        longitudes.units = longitude[1][2]
        levels.units     =     level[1][2]
        hlevels.units     =   levelp[1][2]
        times.units      = timestamp[1][2]
        #complexn.units   = "n/a"
        fourierc.units   = "n/a"
        spharmonics.units= "n/a"
        
        latitudes[:]   = np.array( latitude[0]).astype("float32")
        longitudes[:]  = np.array(longitude[0]).astype("float32")
        levels[:]      = np.array(    level[0]).astype("float32")
        hlevels[:]     = np.array(   levelp[0]).astype("float32")
        times[:]       = np.array(timestamp[0]).astype("float32")
        #complexn[0]["real"] = 1.0; complexn[0]["imag"] = 0.0
        #complexn[0]["real"] = 0.0; complexn[0]["imag"] = 1.0
        fourierc[:] = np.arange(nlons//2,dtype='float32')
        
        specmodes = np.zeros(nmodes)
        w=0
        for m in range(ntru+1):
            for n in range(m,ntru+1):
                specmodes[w  ] = n
                w+=1
        
        spharmonics[:] = np.array(specmodes).astype('float32')
        
        longitudes.axis = 'X'
        fourierc.axis   = 'X'
        spharmonics.axis    = 'X'
        latitudes.axis  = 'Y'
        levels.axis     = 'Z'
        hlevels.axis     = 'Z'
        
        latitudes.standard_name = latitude[1][1]
        longitudes.standard_name= longitude[1][1]
        levels.standard_name    = level[1][1]
        hlevels.standard_name    = levelp[1][1]
        #complexn.standard_name  = "complex_plane"
        fourierc.standard_name  = "fourier_coefficients"
        spharmonics.standard_name   = "spherical_real_modes"
        times.standard_name     = timestamp[1][1]
        
        latitudes.long_name = latitude[1][1]
        longitudes.long_name= longitude[1][1]
        levels.long_name    = "sigma at layer midpoints"
        hlevels.long_name   = "sigma at layer interfaces"
        #complexn.long_name  = "complex coefficients"
        fourierc.long_name  = "Fourier coefficients"
        spharmonics.long_name   = "Spherical harmonic real global modes"
        times.long_name     = timestamp[1][1]
        
        levels.positive  = "down"
        hlevels.positive = "down"
    
    keyvars = list(rdataset.keys())
    keyvars.remove("time")
    keyvars.remove("lat" )
    keyvars.remove("lon" )
    keyvars.remove("lev" )
    keyvars.remove("levp")
    
    for key in keyvars:
        datavar,meta = rdataset[key]
        try:
            dims = meta[4]
        except:
            _log(logfile,meta)
            raise
        shape = datavar.shape
        if "complex" in dims: #Complex dtype
            dims = dims[:-1]
            if not append:
                dvarmask = datavar[np.isfinite(datavar)]
                dvarmask = dvarmask[dvarmask!=0]
                if len(dvarmask)==0:
                    lsd = None
                else:
                    lsd = max(int(round(abs(np.log10(abs(dvarmask).min()))+0.5))+6,6) #get decimal place of smallest value
                    if abs(dvarmask).min()>=1.0:
                        lsd=6
                _log(logfile,"%s Decimal precision: "%key+str(lsd))
                try:
                    variable = ncd.createVariable(key,complex64_t,dims,zlib=True,
                                                  least_significant_digit=lsd)
                except:
                    _log(logfile,meta)
                    raise
                variable.set_auto_mask(False)
            else:
                variable = ncd.variables[key]
            data = np.empty(shape[:-1],complex64)
            data["real"] = datavar[...,0]; data["imag"] = datavar[...,1]
            variable[t0:t1,...] = data
        else:
            if not append:
                dvarmask = datavar[np.isfinite(datavar)]
                dvarmask = dvarmask[dvarmask!=0]
                if len(dvarmask)==0:
                    lsd = None
                else:
                    lsd = max(int(round(abs(np.log10(abs(dvarmask).min()))+0.5))+6,6) #get decimal place of smallest value
                    if abs(dvarmask).min()>=1.0:
                        lsd=6
                _log(logfile,"%s Decimal precision: "%key+str(lsd))
                try:
                    variable = ncd.createVariable(key,"f4",dims,zlib=True,least_significant_digit=lsd)
                except:
                    _log(logfile,meta)
                    raise
                variable.set_auto_mask(False)
                try:
                    variable[:] = datavar[:]
                except:
                    _log(logfile,dims,variable.shape,datavar.shape)
                    _log(logfile,meta)
                    raise
            else:
                variable = ncd.variables[key]
                variable[t0:t1,...] = datavar[:]
        if not append:
            variable.units = meta[2]
            variable.standard_name = meta[1]
            variable.long_name = meta[1]
            variable.units = meta[2]
            variable.code = meta[3]
        if "fourier" not in dims and "modes" not in dims:
            variable.grid_type = "gaussian"
        _log(logfile,"Packing %8s in %s\t....... %d timestamps"%(key,filename,ntimes))
        
    ncd.sync()
    return ncd
    
def npsavez(rdataset,filename="most_output.npz",logfile=None):
    '''Write a dataset to a NumPy compressed .npz file.
    
    Two output files will be created: filename as specified (e.g. most_output.npz), which contains the
    data variables, and a metadata file (e.g. most_output_metadata.npz), which contains the metadata
    headers associated with each variable.
    
    Parameters
    ----------
    rdataset : dict
        A dictionary of outputs as generated from :py:func:`pyburn.dataset()<exoplasim.pyburn.dataset>`
    filename : str, optional
        Path to the output file that should be written.
    logfile : str or None, optional
        If None, log diagnostics will get printed to standard output. Otherwise, the log file
        to which diagnostic output should be written.
        
    Returns
    -------
    tuple
        A 2-item tuple containing (variables, meta), each of which is a 
        dictionary with variable names as keys.
    '''
    variables = {}
    meta = {}
    variables["lat" ] = rdataset["lat" ][0]
    variables["lon" ] = rdataset["lon" ][0]
    variables["lev" ] = rdataset["lev" ][0]
    variables["levp"] = rdataset["levp"][0]
    variables["time"] = rdataset["time"][0]
    meta["lat" ] = rdataset["lat" ][1]
    meta["lon" ] = rdataset["lon" ][1]
    meta["lev" ] = rdataset["lev" ][1]
    meta["levp"] = rdataset["levp"][1]
    meta["time"] = rdataset["time"][1]
    keyvars = list(rdataset.keys())
    keyvars.remove("lat" )
    keyvars.remove('lon' ) 
    keyvars.remove("lev" )
    keyvars.remove("levp")
    keyvars.remove("time")
    for key in keyvars:
        variables[key] = rdataset[key][0].astype("float32")
        vmeta = rdataset[key][1]
        for n in range(len(vmeta)):
            if type(vmeta[n])==tuple or type(vmeta[n])==list:
                vmeta[n] = '/'.join(vmeta[n])
        meta[key] = rdataset[key][1]
        _log(logfile,"Packing %8s in %s\t....... %d timestamps"%(key,filename,rdataset[key][0].shape[0]))
    metafilename = filename[:-4]+"_metadata.npz"
    
    np.savez_compressed(metafilename,**meta)
    np.savez_compressed(filename,**variables)
    return (variables,meta)
    
def _writecsvs(filename,variables,meta,extension=None,logfile=None):
    '''Write CSV output files
    
    Files are placed in a subdirectory named from the filename naming pattern (stripping off the extension).
    
    Parameters
    ----------
    filename : str
        Filename pattern
    variables : dict
        Dictionary of variable data arrays
    meta : dict
        Dictionary of metadata fields for associated variables
    extension : str, optional
        File extension to use for individual files
    logfile : str or None, optional
        If None, log diagnostics will get printed to standard output. Otherwise, the log file
        to which diagnostic output should be written.
        
    Returns
    -------
    list, str 
        List of paths to output files, and the containing directory.
    '''
    idx = filename[::-1].find(".")+1 #Index of last period separator (negative)
    dirname = filename[:-idx]
    if dirname[-4:]==".tar":
        dirname = dirname[:-4]
        idx+=4
    if extension is None:
        extension = filename[-idx:]
    fname = dirname
    if "/" in dirname:
        fname = dirname.split("/")[-1]
    os.system("mkdir %s"%dirname) #Create a subdirectory that just omits the file extension
    files = []
    dimvars = ["lat","lon","lev","levp","time"]
    keyvars = list(variables.keys())
    for var in dimvars:
        outname =  "%s/%s_%s%s"%(dirname,fname,var,extension)
        np.savetxt(outname,variables[var].astype("float32"),
                   header=(str(len(variables[var]))+",|||,"
                          +','.join(meta[var])),delimiter=',')
        keyvars.remove(var)
        files.append(outname)
    maxlen = 0
    for var in keyvars:
        maxlen = max(maxlen,len("%s/%s_%s%s"%(dirname,fname,var,extension)))
    maxlen+=1
    for var in keyvars:
        #This creates e.g. most_output/most_output_ts.csv if filename was most_output.csv
        shape = variables[var].shape
        dim2 = shape[-1]
        dim1 = int(np.prod(shape[:-1]))
        var2d = np.reshape(variables[var],(dim1,dim2)) #np.savetxt can only handle 2 dimensions
        outname =  "%s/%s_%s%s"%(dirname,fname,var,extension)
        try:
            np.savetxt(outname,var2d.astype("float32"),
                   header=(','.join(np.array(shape).astype(str))+",|||,"
                          +','.join(meta[var][:-1])+",".join(meta[var][-1])),delimiter=',')
        except:
            print(",".join(np.array(shape).astype(str)))
            print(",|||,")
            print(meta[var])
            print(','.join(meta[var]))
            print(','.join(np.array(shape).astype(str))+",|||,"
                          +','.join(meta[var]))
            raise
        #The original shape of the array to which it should be reshaped on unpacking is in the header,
        #with the actual metadata separated from the shape by '|||'
        files.append(outname)
        try:
            writeline = "Writing %8s to %"+str(maxlen)+"s\t....... %d timestamps"
            _log(logfile,writeline%(var,outname,variables[var].shape[0]))
        except:
            _log(logfile,"Writing %8s to %s"%(var,filename))
    return files,dirname+"/"
    
def csv(rdataset,filename="most_output.tar.gz",logfile=None,extracompression=False):
    '''Write a dataset to CSV/TXT-type output, optionally compressed.
    
    If a tarball format (e.g. \*.tar or \*.tar.gz) is used, output files will be packed into a tarball.
    gzip (.gz), bzip2 (.bz2), and lzma (.xz) compression types are supported. If a tarball format is 
    not used, then accepted file extensions are .csv, .txt, or .gz. All three will produce a directory
    named following the filename pattern, with one file per variable in the directory. If the .gz extension
    is used, NumPy will compress each output file using gzip compression. 
    
    Files will only contain 2D
    variable information, so the first N-1 dimensions will be flattened. The original variable shape is
    included in the file header (prepended with a # character) as the first items in a comma-separated
    list, with the first non-dimension item given as the '|||' placeholder. On reading variables from these
    files, they should be reshaped according to these dimensions. This is true even in tarballs (which 
    contain CSV files).
    
    Parameters
    ----------
    rdataset : dict
        A dictionary of outputs as generated from :py:func:`pyburn.dataset()<exoplasim.pyburn.dataset>`
    filename : str, optional
        Path to the output file that should be written. This will be parsed to determine output type.
    logfile : str or None, optional
        If None, log diagnostics will get printed to standard output. Otherwise, the log file
        to which diagnostic output should be written.
    extracompression : bool, optional
        If True, then component files in tarball outputs will be compressed individually with gzip, 
        instead of being plain-text CSV files.
        
    Returns
    -------
    tuple or str
        If non-tarball output was used, a tuple containing a list of paths to output files, and a string
        giving the name of the output directory. If tarball output was used, a relative path to the tarball.
    '''
    
    variables = {}
    meta = {}
    for key in rdataset:
        variables[key] = rdataset[key][0]
        meta[key] = rdataset[key][1]
    
    fileparts = filename.split('.')

    if fileparts[-2]=="tar": #We will be doing compression
        import tarfile
        ext = ".csv"
        if extracompression:
            ext = ".gz"
        files,dirname = _writecsvs(filename,variables,meta,extension=ext,logfile=logfile)
        namelen = len(filename)
        with tarfile.open(filename,"w:%s"%fileparts[-1]) as tarball:
            
            maxlen = 0
            for tfile in files:
                maxlen = max(maxlen,len(tfile))
            for var in files:
                varname = var
                if len(varname.split("/"))>2:
                    varname = "/".join(varname.split("/")[-2:])
                tarball.add(var,arcname=varname)
                writeline = "Packing %"+str(maxlen)+"s in %"+str(namelen)+"s"
                try:
                    _log(logfile,(writeline+" ....... %d timestamps")%(var,filename,
                                                                  rdataset[var][0].shape[0]))
                except:
                    _log(logfile,writeline%(var,filename))
        os.system("rm -rf %s"%dirname)
        return filename
        
    elif fileparts[-1]=="tar": #We're still making a tarball, but it won't be compressed
        import tarfile
        ext = ".csv"
        if extracompression:
            ext = ".gz"
        files,dirname = _writecsvs(filename,variables,meta,extension=ext,logfile=logfile)
        namelen = len(filename)
        
        with tarfile.open(filename,"w") as tarball:
            
            maxlen = 0
            for tfile in files:
                maxlen = max(maxlen,len(tfile))
            for var in files:
                tarball.add(var)
                writeline = "Packing %"+str(maxlen)+"s in %"+str(namelen)+"s"
                try:
                    _log(logfile,(writeline+"\t....... %d timestamps")%(var,filename,
                                                                  rdataset[var][0].shape[0]))
                except:
                    _log(logfile,writeline%(var,filename))
        os.system("rm -rf %s"%dirname)
        return filename
        
    else: #Just a collection of CSV/TXT-type files in a subdirectory, which may be individually-compressed.
        #These files can have .txt, .csv, or .gz file extensions.
        files,dirname = _writecsvs(filename,variables,meta,logfile=logfile)
        return files,dirname
     
def hdf5(rdataset,filename="most_output.hdf5",append=False,logfile=None):
    '''Write a dataset to HDF5 output.
    
    Note: HDF5 files are opened in append mode. This means that this format can be used to create
    a single output dataset for an entire simulation.
    
    HDF5 files here are generated with gzip compression at level 9, with chunk rearrangement and
    Fletcher32 checksum data protection.
    
    Parameters
    ----------
    rdataset : dict
        A dictionary of outputs as generated from :py:func:`pyburn.dataset()<exoplasim.pyburn.dataset>`
    filename : str, optional
        Path to the output file that should be written.
    append : bool, optional
        Whether or not the file should be opened in append mode, or overwritten (default). 
    logfile : str or None, optional
        If None, log diagnostics will get printed to standard output. Otherwise, the log file
        to which diagnostic output should be written.
        
    Returns
    -------
    object
        An HDF5 object corresponding to the file that has been written.
    '''
    
    import h5py
    
    mode = "w"
    if append:
        mode = "a"
    hdfile = h5py.File(filename,mode)
    
    latitude  = rdataset["lat" ]
    longitude = rdataset["lon" ]
    level     = rdataset["lev" ]
    levelp    = rdataset["levp"]
    time      = rdataset["time"]
    
    keyvars = list(rdataset.keys())
    keyvars.remove("lat" )
    keyvars.remove("lon" )
    keyvars.remove("lev" )
    keyvars.remove("levp")
    keyvars.remove("time")
    
    #We only add lat, lon, and lev once to the file, so we do it here if the file appears to be new
    if "lat" not in hdfile:
        hdfile.create_dataset("lat",data=latitude[0].astype('float32'),compression='gzip',
                              compression_opts=9,shuffle=True,fletcher32=True)
        hdfile.attrs["lat"] = np.array(latitude[1]).astype('S') #Store metadata
    if "lon" not in hdfile:
        hdfile.create_dataset("lon",data=longitude[0].astype('float32'),compression='gzip',
                              compression_opts=9,shuffle=True,fletcher32=True)
        hdfile.attrs["lon"] = np.array(longitude[1]).astype('S') #Store metadata
    if "lev" not in hdfile:
        hdfile.create_dataset("lev",data=level[0].astype('float32'),compression='gzip',
                              compression_opts=9,shuffle=True,fletcher32=True)
        hdfile.attrs["lev"] = np.array(level[1]).astype('S') #Store metadata
    if "levp" not in hdfile:
        hdfile.create_dataset("levp",data=levelp[0].astype('float32'),compression='gzip',
                              compression_opts=9,shuffle=True,fletcher32=True)
        hdfile.attrs["levp"] = np.array(levelp[1]).astype('S') #Store metadata
    if "time" not in hdfile:
        hdfile.create_dataset("time",data=time[0].astype('float32'),compression='gzip',
                              maxshape=(len(time[0]),),compression_opts=9,
                              shuffle=True,fletcher32=True)
        hdfile.attrs["time"] = np.array(time[1]).astype('S') #Store metadata
    else:
        hdfile["time"].resize((hdfile["time"].shape[0]+len(time[0])),axis=0)
        hdfile["time"][-len(time[0]):] = time[0].astype("float32'")
        
    for var in keyvars:
        if var not in hdfile:
            maxshape = [None,]
            for dim in rdataset[var][0].shape[1:]:
                maxshape.append(dim)
            maxshape=tuple(maxshape)
            hdfile.create_dataset(var,data=rdataset[var][0].astype("float32"),compression="gzip",
                                  maxshape=maxshape,compression_opts=9,shuffle=True,fletcher32=True)
            #_log(logfile,rdataset[var][1])
            #_log(logfile,np.array(rdataset[var][1]).astype('S'))
            meta = rdataset[var][1]
            meta[-1] = ','.join(meta[-1])
            try:
                hdfile.attrs[var] = np.array(meta).astype('S') #Store metadata
            except:
                _log(logfile,meta)
        else:
            hdfile[var].resize((hdfile[var].shape[0]+rdataset[var][0].shape[0]),axis=0) #Expand file
            hdfile[var][-rdataset[var][0].shape[0]:] = rdataset[var][0].astype("float32") #Append
        _log(logfile,"Packing %8s in %s\t....... %d timestamps"%(var,filename,rdataset[var][0].shape[0]))
    return hdfile
             

def postprocess(rawfile,outfile,logfile=None,namelist=None,variables=None,mode='grid',
                zonal=False, substellarlon=180.0, physfilter=False,timeaverage=True,stdev=False,
                times=12,interpolatetimes=True,radius=1.0,gravity=9.80665,gascon=287.0,mars=False):
    '''Convert a raw output file into a postprocessed formatted file.
    
    Output format is determined by the file extension of outfile. Current supported formats are 
    NetCDF (\*.nc), HDF5 (\*.hdf5, \*.he5, \*.h5), numpy's ``np.savez_compressed`` format (\*.npz), and CSV format. If NumPy's 
    single-array .npy extension is used, .npz will be substituted--this is a compressed ZIP archive 
    containing .npy files. Additionally, the CSV output format can be used in compressed form either
    individually by using the .gz file extension, or collectively via tarballs (compressed or uncompressed).
    
    If a tarball format (e.g. \*.tar or \*.tar.gz) is used, output files will be packed into a tarball.
    gzip (.gz), bzip2 (.bz2), and lzma (.xz) compression types are supported. If a tarball format is 
    not used, then accepted file extensions are .csv, .txt, or .gz. All three will produce a directory
    named following the filename pattern, with one file per variable in the directory. If the .gz extension
    is used, NumPy will compress each output file using gzip compression. 
    
    CSV-type files will only contain 2D
    variable information, so the first N-1 dimensions will be flattened. The original variable shape is
    included in the file header (prepended with a # character) as the first items in a comma-separated
    list, with the first non-dimension item given as the '|||' placeholder. On reading variables from these
    files, they should be reshaped according to these dimensions. This is true even in tarballs (which 
    contain CSV files).
    
    A T21 model output with 10 vertical levels, 12 output times, all supported variables in grid 
    mode,and no standard deviation computation will have the following sizes for each format:
    
        +----------------+-----------+
        |Format          | Size      |
        +================+===========+
        |netCDF          | 12.8 MiB  |
        +----------------+-----------+
        |HDF5            | 17.2 MiB  |
        +----------------+-----------+
        |NumPy (default) | 19.3 MiB  |
        +----------------+-----------+
        |tar.xz          | 33.6 MiB  |
        +----------------+-----------+
        |tar.bz2         | 36.8 MiB  |
        +----------------+-----------+
        |gzipped         | 45.9 MiB  |
        +----------------+-----------+
        |uncompressed    | 160.2 MiB |
        +----------------+-----------+
            
    Using the NetCDF (.nc) format requires the netCDF4 python package.
    
    Using the HDF5 format (.h5, .hdf5, .he5) requires the h5py python package.
    
    Parameters
    ----------
    rawfile : str
        Path to the raw output file
    outfile : str
        Path to the destination output file. The file extension determines the format. Currently,
        netCDF (\*.nc). numpy compressed (\*.npz), HDF5 (\*.hdf5, \*.he5, \*.h5), or CSV-type (\*.csv, \*.txt, \*.gz, \*.tar, \*.tar.gz,
        \*.tar.bz2, \*.tar.xz) are supported. If a format (such as npz) that requires
        that metadata be placed in a separate file is chosen, a second file with a '_metadata' suffix will be
        created.
    append : bool, optional
        If True, and outfile already exists, then append to outfile rather than overwriting it. Currently
        only supported for netCDF and HDF5 formats. Support for other formats coming soon.
    logfile : str or None, optional
        If None, log diagnostics will get printed to standard output. Otherwise, the log file
        to which diagnostic output should be written.
    namelist : str, optional
        Path to a burn7 postprocessor namelist file. If not given, then ``variables`` must be set. 
    variables : list or dict, optional
        If a list is given, a list of either variable keycodes (integers or strings), or the abbreviated
        variable name (e.g. 'ts' for surface temperature). If a dict is given, each item in the dictionary
        should have the keycode or variable name as the key, and the desired horizontal mode and additional
        options for that variable as a sub-dict. Each member of the subdict should be passable as **kwargs 
        to :py:func:`advancedDataset() <exoplasim.pyburn.advancedDataset>`. If None, then ``namelist`` must be set.
    mode : str, optional
        Horizontal output mode, if modes are not specified for individual variables. Options are 
        'grid', meaning the Gaussian latitude-longitude grid used
        in ExoPlaSim, 'spectral', meaning spherical harmonics, 
        'fourier', meaning Fourier coefficients and latitudes, 'synchronous', meaning a
        Gaussian latitude-longitude grid in the synchronous coordinate system defined in
        Paradise, et al (2021), with the north pole centered on the substellar point, or
        'syncfourier', meaning Fourier coefficients computed along the dipolar meridians in the
        synchronous coordinate system (e.g. the substellar-antistellar-polar meridian, which is 0 degrees,
        or the substellar-evening-antistellar-morning equatorial meridian, which is 90 degrees). Because this
        will get assigned to the original latitude array, that will become -90 degrees for the polar
        meridian, and 0 degrees for the equatorial meridian, identical to the typical equatorial coordinate
        system.
    zonal : bool, optional
        Whether zonal means should be computed for applicable variables.
    substellarlon : float, optional
        Longitude of the substellar point. Only relevant if a synchronous coordinate output mode is chosen.
    physfilter : bool, optional
        Whether or not a physics filter should be used in spectral transforms.
    times : int or array-like or None, optional
        Either the number of timestamps by which to divide the output, or a list of times given as a fraction
        of the output file duration (which enables e.g. a higher frequency of outputs during periapse of an
        eccentric orbit, when insolation is changing more rapidly). If None, the timestamps in the raw output will be written directly to file.
    timeaverage : bool, optional
        Whether or not timestamps in the output file should be averaged to produce the requested number of 
        output timestamps. Timestamps for averaged outputs will correspond to the middle of the averaged time period.
    stdev : bool, optional
        Whether or not standard deviations should be computed. If timeaverage is True, this will be the 
        standard deviation over the averaged time period; if False, then it will be the standard deviation
        over the whole duration of the output file
    interpolatetimes : bool, optional
        If true, then if the times requested don't correspond to existing timestamps, outputs will be
        linearly interpolated to those times. If false, then nearest-neighbor interpolation will be used.
    radius : float, optional
        Planet radius in Earth radii
    gravity : float, optional
        Surface gravity in m/s^2.
    gascon : float, optional
        Specific gas constant for dry gas (R$_d$) in J/kg/K.  
    mars : bool, optional
        If True, use Mars constants
    
    '''
    #Check output format legality
    
    fileparts = outfile.split('.')
    if fileparts[-1] == "nc":
        pass #OK
    elif fileparts[-1] == "npz" or fileparts[-1] == "npy":
        pass #OK
    elif (fileparts[-1] in ("csv","txt","gz","tar") or \
          (fileparts[-2]+"."+fileparts[-1]) in ("tar.gz","tar.bz2","tar.xz")):
        pass #OK
    elif fileparts[-1] in ("hdf5","h5","he5"):
        pass #OK
    else:
        raise Exception("Unsupported output format detected. Supported formats are:\n\t\n\t%s"%("\n\t".join(SUPPORTED)))
    
    _log(logfile,"==================================")
    _log(logfile,"| PYBURN EXOPLASIM POSTPROCESSOR |")
    _log(logfile,"|  v1.0, Adiv Paradise (C) 2021  |")
    _log(logfile,"==================================")
    _log(logfile,"\n")
    _log(logfile,("--------" +"-"*len(rawfile) + "----"))
    _log(logfile,("Reading %"+"%d"%len(rawfile)+"s ...")%rawfile)
    _log(logfile,("--------" +"-"*len(rawfile) + "----"))
    _log(logfile,"\n")
    
    if namelist is None:
        if variables is None:
            variables = list(ilibrary.keys())
        if mars:
            radius = MARS_RADIUS/6371220.0
            gascon = MARS_RD
            gravity = MARS_GRAV
        if type(variables)==tuple or type(variables)==list:
            data = dataset(rawfile, variables, mode=mode,radius=radius,gravity=gravity,gascon=gascon, 
                           zonal=zonal, substellarlon=substellarlon, 
                           physfilter=physfilter,logfile=logfile)
        elif type(variables)==dict: #for advancedDataset
            data = advancedDataset(rawfile, variables, mode=mode,radius=radius,gravity=gravity,
                                   gascon=gascon, zonal=zonal, substellarlon=substellarlon, 
                                   physfilter=physfilter, logfile=logfile)
        
    else:
        #Scrape namelist
        variables=None
        with open(namelist,"r") as fname:
            fnamelist = fname.read().split('\n')
        for line in fnamelist:
            parts = line.split('=')
            if len(parts)>1:
                key = parts[0]
                arg = parts[1]
                if key=="code" or key=="CODE":
                    variables = arg.split(',')
                elif key=="htype" or key=="HYTPE":
                    if arg=="g" or arg=="G":
                        mode="grid"
                    elif arg=="z" or arg=="Z":
                        mode="grid"
                        zonal=True
                    elif arg=="s" or arg=="S":
                        mode="spectral"
                    elif arg=="f" or arg=="F":
                        mode="fourier"
                elif key=="mean" or key=="MEAN":
                    if int(arg)==0:
                        timeaverage=False
                        stdev = False
                    elif int(arg)==1:
                        timeaverage=True
                        stdev = False
                    elif int(arg)==2:
                        stdev = True
                    elif int(arg)==3:
                        timeaverage=True
                        stdev = True
                elif key=="radius" or key=="RADIUS":
                    radius = float(arg)/6371220.0      #Radii in the namelist are given in metres
                elif key=="gravity" or key=="GRAVITY":
                    gravity = float(arg)
                elif (key=="mars" or key=="MARS") and int(arg)==1:
                    gravity = MARS_GRAV
                    radius  = MARS_RADIUS/6371220.0    #We want to start off with radii in Earth radii
                    gascon  = MARS_RD      #This is called RD in burn7, not gascon
        data = dataset(rawfile, variables, mode=mode,radius=radius,gravity=gravity,gascon=gascon, 
                       zonal=zonal, substellarlon=substellarlon, physfilter=physfilter,logfile=logfile)
        
    # Compute time averages, binning, stdev, etc
    
    dtimes = data["time"][0]
    ntimes = len(dtimes)
    
    if times is None:
        times = ntimes
    
    if type(times)==int: #A number of time outputs was specified
        times = max(times,1)
        if ntimes==times: #The number of outputs exactly equals the number provided
            timeaverage = False #This way stdev still gets computed, but over the whole file.
        
        else:
            if timeaverage:
               _log(logfile,"\nComputing averages, going from %d timestamps to %d ..."%(ntimes,times))
               if times>ntimes and interpolatetimes:
                   _log(logfile,
                        "Interpolating by a factor of 10 to compute averages at super-resolution....")
                   ttimes = np.linspace(dtimes[0],dtimes[-1],num=10*ntimes)
                   indices = np.digitize(np.linspace(dtimes[0],dtimes[-1],num=times+1),ttimes)-1
                   indices[-1] = len(ttimes)
                   counts = np.diff(indices)
                   newtimes = np.linspace(dtimes[0],dtimes[-1],num=times+1)
                   newtimes = 0.5*(newtimes[:-1]+newtimes[1:])
                   varkeys = list(data.keys())
                   varkeys.remove("time")
                   varkeys.remove("lat")
                   varkeys.remove("lon")
                   varkeys.remove("lev")
                   varkeys.remove("levp")
                   if stdev:
                       _log(logfile,"Computing standard deviations ....")
                   for var in varkeys:
                        odata = data[var][0][:]
                        interpfunc = scipy.interpolate.interp1d(dtimes,odata,axis=0,kind='linear')
                        tempdata = interpfunc(ttimes)
                        newshape = list(odata.shape)
                        newshape[0] = times
                        newshape = tuple(newshape)
                        newcounts = np.transpose(np.resize(counts,newshape[::-1]))
                        data[var][0] = np.add.reduceat(tempdata,indices[:-1],axis=0)/ newcounts
                        if stdev:  #Compute standard deviation
                            stdvar = np.zeros(data[var][0].shape)
                            for nt in range(stdvar.shape[0]):
                                stdvar[nt,...] = np.std(tempdata[indices[nt]:indices[nt+1]],axis=0)
                            stdmeta = list(data[var][1][:])
                            stdmeta[0]+="_std"
                            stdmeta[1]+="_standard_deviation"
                            data[var+"_std"] = (stdvar,stdmeta)
               else:
                   indices = np.linspace(0,ntimes,times+1,True).astype(int)
                   counts = np.diff(indices)
                   varkeys = list(data.keys())
                   varkeys.remove("time")
                   varkeys.remove("lat")
                   varkeys.remove("lon")
                   varkeys.remove("lev")
                   varkeys.remove("levp")
                   newtimes = np.add.reduceat(dtimes,indices[:-1]) / counts
                   if stdev:
                       _log(logfile,"Computing standard deviations ....")
                   for var in varkeys:
                       odata = data[var][0][:]
                       newshape = list(odata.shape)
                       newshape[0] = times
                       newshape = tuple(newshape)
                       newcounts = np.transpose(np.resize(counts,newshape[::-1]))
                       data[var][0] = np.add.reduceat(odata,indices[:-1],axis=0) / newcounts
                       if stdev:  #Compute standard deviation
                           stdvar = np.zeros(data[var][0].shape)
                           for nt in range(stdvar.shape[0]):
                               stdvar[nt,...] = np.std(odata[indices[nt]:indices[nt+1]],axis=0)
                           stdmeta = list(data[var][1][:])
                           stdmeta[0]+="_std"
                           stdmeta[1]+="_standard_deviation"
                           data[var+"_std"] = (stdvar,stdmeta)
               data["time"][0] = newtimes
            else:
               if interpolatetimes:
                   interpolation = "linear"
                   _log(logfile,"\nInterpolating from %d timestamps to %d ..."%(ntimes,times))
               else:
                   interpolation = "nearest"
                   _log(logfile,"\nSelecting %d timestamps from %d ..."%(times,ntimes))
               newtimes = np.linspace(dtimes[0],dtimes[-1],num=times) #We can't extrapolate; only interpolate
               varkeys = list(data.keys())
               varkeys.remove("time")
               varkeys.remove("lat")
               varkeys.remove("lon")
               varkeys.remove("lev")
               varkeys.remove("levp")
               for var in varkeys:
                   odata = data[var][0][:]
                   interpfunc = scipy.interpolate.interp1d(dtimes,odata,axis=0,kind=interpolation)
                   data[var][0] = interpfunc(newtimes)
               if not interpolatetimes:
                   interpfunc = scipy.interpolate.interp1d(dtimes,dtimes,kind="nearest")
                   data["time"][0] = interpfunc(newtimes)
               else:
                   data["time"][0] = newtimes
                   
    else: #A list of times was specified
        
        if timeaverage: #times values are assumed to be bin edges
            _log(logfile,"\nComputing averages, going from %d timestamps to %d ..."%(ntimes,len(times)-1))
            if not interpolatetimes: #We will always round down to the nearest neighbor
                _log(logfile,
                     "Interpolation disabled, so bin edges are being selected via nearest-neighbor.")
                indices = np.digitize(np.array(times)*(dtimes[-1]-dtimes[0])+dtimes[0],dtimes)-1
                counts = np.diff(indices)
                varkeys = list(data.keys())
                varkeys.remove("time")
                varkeys.remove("lat")
                varkeys.remove("lon")
                varkeys.remove("lev")
                varkeys.remove("levp")
                newtimes = np.add.reduceat(dtimes,indices[:-1]) / counts
                if stdev:
                    _log(logfile,"Computing standard deviations ....")
                for var in varkeys:
                    odata = data[var][0][:]
                    newshape = list(odata.shape)
                    newshape[0] = len(times)-1
                    newshape = tuple(newshape)
                    newcounts = np.transpose(np.resize(counts,newshape[::-1]))
                    data[var][0] = np.add.reduceat(odata,indices[:-1],axis=0) / newcounts
                    if stdev:  #Compute standard deviation
                        stdvar = np.zeros(data[var][0].shape)
                        for nt in range(stdvar.shape[0]):
                            stdvar[nt,...] = np.std(odata[indices[nt]:indices[nt+1]],axis=0)
                        stdmeta = list(data[var][1][:])
                        stdmeta[0]+="_std"
                        stdmeta[1]+="_standard_deviation"
                        data[var+"_std"] = (stdvar,stdmeta)
                data["time"][0] = newtimes
            else: #First we'll interpolate to high time resolution, then compute average via binning
                _log(logfile,"Interpolating to find bin edges....")
                ttimes = np.linspace(dtimes[0],dtimes[-1],num=10*ntimes)
                indices = np.digitize(np.array(times)*(dtimes[-1]-dtimes[0])+dtimes[0],ttimes)-1
                counts = np.diff(indices)
                varkeys = list(data.keys())
                varkeys.remove("time")
                varkeys.remove("lat")
                varkeys.remove("lon")
                varkeys.remove("lev")
                varkeys.remove("levp")
                if stdev:
                    _log(logfile,"Computing standard deviations ....")
                for var in varkeys:
                    odata = data[var][0][:]
                    interpfunc = scipy.interpolate.interp1d(dtimes,odata,axis=0,kind='linear')
                    tempdata = interpfunc(ttimes)
                    newshape = list(odata.shape)
                    newshape[0] = len(times)-1
                    newshape = tuple(newshape)
                    newcounts = np.transpose(np.resize(counts,newshape[::-1]))
                    data[var][0] = np.add.reduceat(tempdata,indices[:-1],axis=0)/ newcounts
                    if stdev:  #Compute standard deviation
                        stdvar = np.zeros(data[var][0].shape)
                        for nt in range(stdvar.shape[0]):
                            stdvar[nt,...] = np.std(tempdata[indices[nt]:indices[nt+1]],axis=0)
                        stdmeta = list(data[var][1][:])
                        stdmeta[0]+="_std"
                        stdmeta[1]+="_standard_deviation"
                        data[var+"_std"] = (stdvar,stdmeta)
                newtimes = np.array(times)
                newtimes = 0.5*(newtimes[:-1]+newtimes[1:])*dtimes[-1]
                data["time"][0] = newtimes
                    
        else:
            newtimes = np.array(times)*(dtimes[-1]-dtimes[0])+dtimes[0] #Convert to timestamps
            if interpolatetimes:
                interpolation = "linear"
                _log(logfile,"\nInterpolating from %d timestamps to %d ..."%(ntimes,len(times)))
            else:
                interpolation = "nearest"
                _log(logfile,"\nSelecting %d timestamps from %d ..."%(len(times),ntimes))
            varkeys = list(data.keys())
            varkeys.remove("time")
            varkeys.remove("lat")
            varkeys.remove("lon")
            varkeys.remove("lev")
            varkeys.remove("levp")
            for var in varkeys:
                odata = data[var][0][:]
                interpfunc = scipy.interpolate.interp1d(dtimes,odata,axis=0,kind=interpolation)
                data[var][0] = interpfunc(newtimes)
            if not interpolatetimes:
                interpfunc = scipy.interpolate.interp1d(dtimes,dtimes,kind="nearest")
                data["time"][0] = interpfunc(newtimes)
            else:
                data["time"][0] = newtimes
        
    #Standard deviation if timeaverage is False
    if not timeaverage and stdev:
        _log(logfile,"\nComputing standard deviations ....")
        for var in varkeys:
            stdvar = np.std(data[var],axis=0)
            stdmeta = list(data[var][1][:])
            stdmeta[0]+="_std"
            stdmeta[1]+="_standard_deviation"
            data[var+"_std"] = (stdvar,stdmeta)
    
        
    # Write to output
    
    _log(logfile,"\n")
    _log(logfile,("--------" +"-"*len(outfile) + "----"))
    _log(logfile,("Writing %"+"%d"%len(outfile)+"s ...")%outfile)
    _log(logfile,("--------" +"-"*len(outfile) + "----"))
    _log(logfile,"\n")
    
    fileparts = outfile.split('.')
    if fileparts[-1] == "nc":
        output=netcdf(data,filename=outfile,logfile=logfile)
        output.close()
    elif fileparts[-1] == "npz" or fileparts[-1] == "npy":
        output=npsavez(data,filename=outfile,logfile=logfile)
    elif (fileparts[-1] in ("csv","txt","gz","tar") or \
          (fileparts[-2]+"."+fileparts[-1]) in ("tar.gz","tar.bz2","tar.xz")):
        output=csv(data,filename=outfile,logfile=logfile)
    elif fileparts[-1] in ("hdf5","h5","he5"):
        output=hdf5(data,filename=outfile,logfile=logfile)
        output.close()
    else:
        raise Exception("Unsupported output format detected. Supported formats are:\n\t\n\t%s"%("\n\t".join(SUPPORTED)))
    
    _log(logfile,"\n")
    _log(logfile,"%s closed."%outfile)
    _log(logfile,"\n")
    
    _log(logfile,"================================")
    _log(logfile,"| PYBURN FINISHED SUCCESSFULLY |")
    _log(logfile,"================================")
    
    
    
        
    
    

        
        
        
    
    
