import numpy as np
import pandas as pd
import h5py

tracked_filepath = "/home/christenc/Code/python/Justice_League_Code/Data/tracked_particles.hdf5"
def find_halo_keys(simulations):

    """
    pre-processing to aid in aprticle tracking -- takes in the filepath to all tracked particles, and divides them by relevant halos
    ie identifies which halos belong in which simulations
    """

    #load pandas file with particle tracking for all four sims
    halo_keys= np.array(list(h5py.File(tracked_filepath).keys()))#load in tracked particle data for all the halos and sims
    halo_keys = sorted(halo_keys, key=lambda x: int(x.split('_', 1)[1]))
    #sort halo keys so that it appears in an order of decending halo size for all simulations


    ##make arrays to store the respective halos in each simulation
    h148 = []
    h229 = []
    h329 = []
    h242 = []

    for key in halo_keys:#create seperate lists with only the keys (ie halos) pertaining to each simulation. 
        #As the halos where added in decending size in the original sorted keys arays, each of these will also be sorted in decending halo size
        
        if key.startswith('h148'):
            h148.append(key)
        elif key.startswith('h229'):
            h229.append(key)
        elif key.startswith('h242'):
            h242.append(key)
        elif key.startswith('h329'):
            h329.append(key)
    halos = pd.Series({'Sandra': h148, 'Ruth': h229, 'Sonia': h242, 'Elena':h329})
    simulations['Halo keys'] = halos
    return()

def find_halo_particles(h1, simulations, filename = 'Sandra'):
    """
    for each halo in the list of halos, find all tracked particles, and add them toa  dictionary where they can be pulled from
    """
    
    halo_subsims = {}#create dictionary
    all_halos = np.array([])#create a list for particles in any halo
    
    halolist = simulations['Halo keys'][filename]#get the list of halos for the given simulation
    
    for halo in halolist:
        # particle_list = [] ##array to store all particle ids within it
        halo_of_origin = "h"+(halo.split('_')[1]) #string denoting the halo of origin
    
        tracked_particles = pd.read_hdf(tracked_filepath, key = halo)#find all information on particles tracked in given halo
        halo_particle_IDs = (tracked_particles['pid'].to_numpy())#isolate and create a np list of PID's 
        bools = np.isin(h1.g['iord'], halo_particle_IDs)#Using the particle ID's defined above, for every particle in the central
        #galaxy make a boolean list of wether or not it is in the satelite currently observed
    
        halo_subsims[halo_of_origin] =  h1.g[bools]#create a dictionary entry with the
        #list of all particles ID's(particles), and turn the boolean list to have a callable pynbody data file
        
        
        all_halos = np.append(all_halos, halo_particle_IDs)#add PID's of this specific halo, to the list of the particles accreted
        #from all satelites
        
    all_bool = np.isin(h1.g['iord'], all_halos)
    halo_subsims['all_halos'] =  h1.g[all_bool]#repeat steps above for list of all acreted particles
    
    local_bool = np.invert(all_bool) 
    halo_subsims['local'] = h1.g[local_bool]#flip the boolean array to get a callable list 
    return halo_subsims
    #of just the particles in the central halo