import pynbody
import numpy as np
import pynbody.plot.sph as sph
from pynbody import units

#Physical constants used in the derrivations of number of particles.
HI_molarMass = 1.007940
OVI_molarMass = 15.99940
g_in_msol = 1.989e33
cm_in_kpc = 3.086e21

mol = 6.022e23
#number of particles in one Msol for HI and OVI
HI_N_in_Msol = g_in_msol * mol / HI_molarMass
OVI_N_in_Msol= g_in_msol * mol / OVI_molarMass  

#two factors used for converting surface density from MSol/kpc**2 (units initially given in defining covering fractions)
#to N/cm**2 (units used in literature)
HI_factor = HI_N_in_Msol/(cm_in_kpc**2) 
OVI_factor = OVI_N_in_Msol/(cm_in_kpc**2) 







def calculate_gas_mass(s):
    """
    calculates the oxygen six and hydrogren one gas masses and adds them to the respective sim values. 
    """
    # Calculate for each gas particle the fraction of the hydrogen that is neutral hydrogen
    s.gas['hiif'] = pynbody.analysis.ionfrac.calculate(s.gas,ion='hi')
    # Calculate the fraction of oxygen that is OVI
    s.gas['oviif'] = pynbody.analysis.ionfrac.calculate(s.gas,ion='ovi')
    
    #holder variable for mass, so that if "mass" is changed for particles (such as for HI and OVI rho calculations),
    #no dependencies are affected
    @pynbody.snapshot.simsnap.SimSnap.stable_derived_array
    def mass_holder(sim):
        return sim.g['mass']
    
    ##fractions of HI and OVI mass relative to the entire particle
    
    @pynbody.derived_array
    def HI_frac(sim):
        return sim.g['hiif']*sim.gas['hydrogen']
        
    @pynbody.derived_array
    def OVI_frac(sim):
        
        return sim.g['oviif']*sim.gas['OxMassFrac']    
    #mass of HI and OVI in a given particle
    @pynbody.derived_array
    def HI_mass(sim):
        return sim.gas['mass_holder']* sim.gas['HI_frac']
    
    @pynbody.derived_array
    def OVI_mass(sim):
        return np.multiply(sim.gas['mass_holder'], sim.gas['OVI_frac'])
    
    
    # number of atoms of HI and OVI in a given particle.
    @pynbody.derived_array
    def HI_N(sim):
        return (sim.gas['HI_mass']*HI_N_in_Msol)
    
    @pynbody.derived_array
    def OVI_N(sim):
        return (sim.gas['OVI_mass']*OVI_N_in_Msol)
    
    #set units as relevant for mass which needs it.
    s.g['OVI_mass'].units = units.Unit('Msol')
    s.g['HI_mass'].units = units.Unit('Msol')

def calculate_gas_rhos(s):
    """
    in order to take advantage of pynbodys SPH modeling to generate our covering fractions, the pynbody.image function was used.
    this allows us to get a visual image of density(rho), which takes into account sph modeling to avoid having all the mass from particles located
    on their center.
    
    further detail on how this is done can be found in the "get CD" function, however in order for this to work, the rho for OVI and HI must also 
    be recalculated, as while each particles mass is independent of its surrounding (allowing us to simply multiply it by the respective ion fraction
    to derrive this) density depends on the mass of the surrounding particles.
    
    no built in pynbody function exists that does this directly, however, by setting the mass of the particle to equal the mass of OVI or HI 
    the call to sph.rho will return an array of particle masses based on the mass of the surrounding elements, hence while janky, the above is done.
    
    
    tldr: to use pynbodys sph functionality in calculating covering fractions, the density(rho) of OVI  and HI particles is needed
    """
    
    #set the mass of all particles to their HI mass, (np.copyto is needed as to not mess with units)
    np.copyto(s.g['mass'], s.g['HI_mass'])
    s.gas['HI_rho'] = pynbody.sph.rho(s.gas) #get an array of each particles respective "HI density (sph.rho automatically pulls the ['mass'] in
    
    #repeat for OVI
    np.copyto(s.g['mass'], s.g['OVI_mass'])
    s.gas['OVI_rho'] = pynbody.sph.rho(s.gas)
    
    #as the original mass values are now OVI, pull in from the stable "mass holder" and reset them. Reset rho as well.
    np.copyto(s.g['mass'], s.g['mass_holder'])
    s.gas['rho'] = pynbody.sph.rho(s.gas)