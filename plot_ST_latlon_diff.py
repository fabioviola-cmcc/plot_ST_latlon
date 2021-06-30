#!/usr/bin/python3

# ignore warnings
import warnings
warnings.filterwarnings("ignore")

# global requirements
import os
import sys
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap


##############################################
#
# Read input param
#
##############################################

# Open the NetCDF file
try:
    filename_A = sys.argv[1]
    b_filename_A = os.path.basename(filename_A)    
except:
    print("ERROR: you must provide a NetCDF file as input!")
    sys.exit(255)

try:
    filename_B = sys.argv[2]
    b_filename_B = os.path.basename(filename_B)    
except:
    print("ERROR: you must provide a NetCDF file as input!")
    sys.exit(255)    


##############################################
#
# Determine area, regions and model for FILE A
#
##############################################

# nea or nwp?
filetype = None
regions = {}
if "NEA_re" in b_filename_A.split("-"):
    filetype = 'NEA'
    lats = [40, 34]
    regions = {
        1: { "lat": [40], "lon": [-45, 27], "min": 33, "max": 40 },
        2: { "lat": [34], "lon": [-45, 36], "min": 33, "max": 40 },
        3: { "lat": [31, 67], "lon": [-35], "min": -35, "max": 40 }
    }    
else:
    filetype = 'NWP'
    regions = {
        1: { "lat": [10], "lon": [106, 154], "min": 32, "max": 36 },
        2: { "lat": [20], "lon": [105, 154], "min": 32, "max": 36 },
        3: { "lat": [6, 40], "lon": [123], "min": 33, "max": 40 }
    }

datavar = None
if "TEMP" in b_filename_A.split("-"):
    datavar = "TEMP"
else:
    datavar = "PSAL"
    
model_A = None
if "CGLORS" in b_filename_A:
    model_A = "C-GLORSv7"
elif "GREPv2" in b_filename_A:
    model_A = "GREPv2e"
else:
    model_A = "GLORYSv12"

    
##############################################
#
# Determine area, regions and model for FILE B
#
##############################################

# nea or nwp?
filetype = None
regions = {}
if "NEA_re" in b_filename_B.split("-"):
    filetype = 'NEA'
    lats = [40, 34]
    regions = {
        1: { "lat": [40], "lon": [-45, 27], "min": 33, "max": 40 },
        2: { "lat": [34], "lon": [-45, 36], "min": 33, "max": 40 },
        3: { "lat": [31, 67], "lon": [-35], "min": -35, "max": 40 }
    }    
else:
    filetype = 'NWP'
    regions = {
        1: { "lat": [10], "lon": [106, 154], "min": 32, "max": 36 },
        2: { "lat": [20], "lon": [105, 154], "min": 32, "max": 36 },
        3: { "lat": [6, 40], "lon": [123], "min": 33, "max": 40 }
    }

datavar = None
if "TEMP" in b_filename_B.split("-"):
    datavar = "TEMP"
else:
    datavar = "PSAL"
    
model_B = None
if "CGLORS" in b_filename_B:
    model_B = "C-GLORSv7"
elif "GREPv2" in b_filename_B:
    model_B = "GREPv2e"
else:
    model_B = "GLORYSv12"
    

##############################################
#
# Open and process the NetCDF file
#
##############################################

ncfile_A = xr.open_dataset(filename_A)
ncfile_B = xr.open_dataset(filename_B)

# iterate over lats
for region in regions:

    ############ FILE A
    
    # Extract the desired data for filename A
    indata_A = None    
    if len(regions[region]["lon"]) == 2:

        if datavar == "TEMP":
            indata_A = ncfile_A.thetao.sel(lon=slice(regions[region]["lon"][0], regions[region]["lon"][1])).sel(lat=regions[region]["lat"][0], method='nearest')[0,:,:]
        else:
            indata_A = ncfile_A.so.sel(lon=slice(regions[region]["lon"][0], regions[region]["lon"][1])).sel(lat=regions[region]["lat"][0], method='nearest')[0,:,:]
                  
        # get indexes
        minlon=regions[region]["lon"][0] 
        maxlon=regions[region]["lon"][1]
        lons = ncfile_A.lon
        lons_indexes = np.where((lons >= minlon) & (lons <= maxlon))[0]
        lons_sel = lons[lons_indexes]
                
    else:
        if datavar == "TEMP":
            indata_A = ncfile_A.thetao.sel(lat=slice(regions[region]["lat"][0], regions[region]["lat"][1])).sel(lon=regions[region]["lon"][0], method='nearest')[0,:,:]
        else:
            indata_A = ncfile_A.so.sel(lat=slice(regions[region]["lat"][0], regions[region]["lat"][1])).sel(lon=regions[region]["lon"][0], method='nearest')[0,:,:]

        # get indexes
        minlat=regions[region]["lat"][0] 
        maxlat=regions[region]["lat"][1]
        lats = ncfile_A.lat
        lats_indexes = np.where((lats >= minlat) & (lats <= maxlat))[0] 
        lats_sel = lats[lats_indexes]            

    # depths
    depths = ncfile_A.depth
    
    minl=regions[region]["min"]
    maxl=regions[region]["max"]
    
    ############## FILE B

    # Extract the desired data for filename B
    indata_B = None    
    if len(regions[region]["lon"]) == 2:

        if datavar == "TEMP":
            indata_B = ncfile_B.thetao.sel(lon=slice(regions[region]["lon"][0], regions[region]["lon"][1])).sel(lat=regions[region]["lat"][0], method='nearest')[0,:,:]
        else:
            indata_B = ncfile_B.so.sel(lon=slice(regions[region]["lon"][0], regions[region]["lon"][1])).sel(lat=regions[region]["lat"][0], method='nearest')[0,:,:]
                
    else:
        if datavar == "TEMP":
            indata_B = ncfile_B.thetao.sel(lat=slice(regions[region]["lat"][0], regions[region]["lat"][1])).sel(lon=regions[region]["lon"][0], method='nearest')[0,:,:]
        else:
            indata_B = ncfile_B.so.sel(lat=slice(regions[region]["lat"][0], regions[region]["lat"][1])).sel(lon=regions[region]["lon"][0], method='nearest')[0,:,:]

            
    #######################################################################
    #
    # PLOT
    #
    #######################################################################

    # calculate data
    indata = indata_A - indata_B
    
    # Create the figure
    fig = plt.figure(figsize=(12,6))
    ax = plt.axes()

    # Set the indexes of coordinates
    depths_indexes = np.where((depths >= depths.min()) & (depths <= depths.max()))[0]     
    time_indexes = 0

    # Extract the coordinates subsets
    depths_sel = -depths[depths_indexes]
    
    # Create the meshgrid for the plot
    if len(regions[region]["lon"]) == 2:
        xx, yy = np.meshgrid(lons_sel, depths_sel)
    else:
        xx, yy = np.meshgrid(lats_sel, depths_sel)
    
    # contour fill
    min_value = regions[region]["min"]
    max_value = regions[region]["max"]
    step_value = 0.02
    # contour_levels = np.arange(min_value, max_value, step_value)
    contour_levels = np.arange(-2, 2, step_value)
    plt.contourf(xx, yy, indata, contour_levels, vmin=-2, vmax=2, cmap="bwr") #, vmin=min_value, vmax=max_value, extend='both')

    # Set the color bar
    plt.colorbar(extend='both')

    # Set the x-axis and  y-axis labels and figure title
    if len(regions[region]["lat"]) == 2:
        ax.set_xlabel("Longitude", fontsize=12)
        degUnit = "N"
    else:
        ax.set_xlabel("Latitude", fontsize=12)
        degUnit = "E"
    ax.set_ylabel("Depth (m)", fontsize=12)
    
    if datavar is "TEMP":
        plt.title("%s [2] Temperature difference (degC) (%s%s, 2004-2019 mean, %s VS %s rean)" % (filetype.split("_")[0], regions[region]["lat"], degUnit, model_A, model_B))
    else:
        plt.title("%s [2] Salinity (psu) (%s%s, 2004-2019 mean, %s VS %s rean)" % (filetype.split("_")[0], regions[region]["lat"], degUnit, model_A, model_B))

    # show the figure
    if not os.path.exists("./images"):
        os.mkdir("./images")
    if len(regions[region]["lat"]) == 2:
        plt.savefig("./images/diff_%s_%s_%s%s_%s_VS_%s_2004_2019.png" % (filetype.split("_")[0],
                                                                   datavar,
                                                                   regions[region]["lon"][0], degUnit, model_A, model_B))
        print("Generated ./images/diff_%s_%s_%s%s_%s_VS_%s_2004_2019.png" % (filetype.split("_")[0],
                                                                       datavar,
                                                                       regions[region]["lon"][0], degUnit, model_A, model_B))
    else:        
        plt.savefig("./images/diff_%s_%s_%s%s_%s_VS_%s_2004_2019.png" % (filetype.split("_")[0],
                                                                         datavar,
                                                                         regions[region]["lat"][0], degUnit, model_A, model_B))
        print("Generated ./images/diff_%s_%s_%s%s_%s_VS_%s2004_2019.png" % (filetype.split("_")[0],
                                                                            datavar,
                                                                            regions[region]["lat"][0], degUnit, model_A, model_B))


