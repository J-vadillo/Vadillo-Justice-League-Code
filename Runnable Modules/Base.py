import pandas as pd
import pynbody
import numpy as np
import h5py
import pynbody.plot.sph as sph
import matplotlib as mpl
from matplotlib import colors
from pynbody import units

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

# set up matplotlib preferences
mpl.rc('font',**{'family':'serif','monospace':['Palatino']})
mpl.rc('text', usetex=True)
mpl.rcParams.update({'figure.dpi': 200,
                     'font.size': 9,
                     'xtick.direction': 'in',
                     'ytick.direction': 'in',
                     'legend.frameon': False,
                     'figure.constrained_layout.use': True,
                     'xtick.top': True,
                     'ytick.right': True})



def Polar_Profile(image_array, edge, averege=True, radial=True, bin_width=5):
    
    """
    Makes a polar profile of the given image array! it can either make a radial profoile (the default) or a angular profile,  and  

    Works by taking the total sum of bins that fall within a set amount of radial or angular profiles, and if averege = true: returns the averege value, wherass if not
    it returns the total sum of values that fall within that bin. 
    
    Mainly used in this project to calculate covering fractions by passing a 2d image array of where the density crosses a detectability threshold and where it does not: and returns
    the averege number of bins that cross the threshold at each radius, or angular marker.

    Parameters:
    image_array (np.ndarray): A 2D array representing the image data from which the polar profile is calculated.
    edge (int): The maximum radius (in kpc) to consider for the calculation. Treaitionally will be "edge" ie np.max
    averege (bool): If True, returns the averege in each bin; otherwise, returns raw counts. Default is True.
    radial (bool): If True, calculates radial metrics; otherwise, caulculates angular metrics. Default is True.
    bin_width (float): The width(in kpc) of each bin in the histogram. Default is 5.

    Returns:
    tuple: A tuple containing:
        - Covering_fraction (np.ndarray): The fraction of bins that cross the threshold at each radius.
        - np.sqrt(Threshold_crossing[1]) (np.ndarray): The square root of the bin edges for the histogram.
    """
    
    profile_values = []  # an array to store the equivalent radial or anguilar distance of each bin within the given array
    values = []  # an array to store the corresponding value of if said bin is covered or not
    if len(image_array) != len(image_array[0]):
        raise Exception("Sorry, image array must be square") 
    ncells = len(image_array) #the number of cells in the image array along a given axis. Image is assumed to be square.
    cell_width = 2*edge / ncells  # defines how many kpc's each cell in the image currently possesses.
   

    edge_squared = edge**2 ##useful for calculating if a given bin is within the virial radius.
    midpoint = len(image_array) / 2 ## the midpoitn of the image array useful for calculating the radial distance of a given bin.

    for column in range(ncells):
        for row in range(ncells):
            x_position = (column - midpoint) * cell_width
            y_position = (row - midpoint) * cell_width
            rad_sq = int((x_position**2 + y_position**2))

            if rad_sq <= edge_squared:  # only add in values within the virial radius.
                if radial:
                    profile_values.append(rad_sq)
                else:
                    angle = np.arctan2(y_position, x_position)
                    if angle < 0:
                        angle += 2 * np.pi
                    profile_values.append(angle)
                    
                values.append(image_array[column, row])
    if radial:
        bins = np.square(np.linspace(0, edge, int(edge / bin_width)))
        bin_edges = np.sqrt(bins)  # convert to radial distances
    else:
        bins = np.linspace(0, 2 * np.pi, int(2 * np.pi / bin_width))
        bin_edges = bins
    values = np.histogram(profile_values, bins=bins, weights=values)
    
    if averege:
        Total = np.histogram(profile_values, bins=bins)
        return_values = np.where(Total[0] != 0, values[0] / Total[0], 0)
        values = return_values
    return (values, bin_edges)
