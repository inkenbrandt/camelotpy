import math as math

def FindSpecificStorage(density_fluid, compress_soil, \
                        compress_fluid, porosity, gravity=9.81):
    '''Finds the specific storage (Ss) from soil and fluid
       compressibility, fluid density, and porosity.
       Ss = density_fluid * gravity * (compress_soil +
            compress_fluid * porosity)
       density_fluid     density of fluid [M/L^3]
       compress_soil     compressibility of soil [L^2/F]
       compress_fluid    compressibility of fluid [L^2/F]
       porosity          porosity of aquifer [-]
       gravity           gravitational accleration [L/T^2]
                         (default = 9.81)'''
    #return specific storage
    return density_fluid * gravity * (compress_soil + \
                         compress_fluid * porosity)

def FindStorativity(specific_storage, thickness):
    '''Finds the storativity (S) from specific storage and
       aquifer thickness.
       S = specific_storage * thickness
       specific_storage   specific storage of aquifer [1/L]
       thickness          aquifer thickness [L]'''
    #return storativity
    return specific_storage * thickness

def FindHydraulicConductivity(permeability, density_fluid, \
                              viscosity_fluid, gravity=9.81):
    '''Finds the hydraulic conductivity (K) from the permeability
       (also known as intrinsic conductivity), fluid density,
       fluid viscosity, and gravity
       K = permeability * density_fluid * gravity / viscosity_fluid
       permeability      permeability of aquifer [L^2]
       density_fluid     density of fluid [M/L^3]
       viscosity_fluid   dynamic viscosity of the fluid [FT/L^2]
       gravity           gravitational accleration [L/T^2]
                         (default = 9.81)'''
    #return hydraulic conductivity
    return permeability * density_fluid * gravity / viscosity_fluid

def FindVolumetricInjectionRate(massrate, density):
    '''Finds the volumetric injection rate (Q) from mass injection
       rate (M) and density (rho)
       Q = M / rho
       massrate   mass injection rate [M/T]
       density    density of injected fluid [M/L^3]'''
    #return volumetric injection rate
    return massrate / density

def Interpolate(x1, y1, x2, y2, yknown):
    '''Linear interpolation between two points (x1,y1) and (x2,y2).
       yknown is the point between points 1 and 2 for which the
       interpolated value is sought.'''
    delx = x2 - x1
    dely = y1- y2
    return x2 - (yknown - y2) * delx / dely
