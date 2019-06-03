import os, sys
import time
import mne
import numpy as np
import scipy.io as io
import scipy.signal as signal
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import spectralevents_functions as tse

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Top-level run script for finding spectral events in test data 
# provided by original SpectralEvents GitHub repo.

####################################################################
# Main Code
####################################################################

# Set dataset and analysis parameters
numSubj = 10

eventBand = [15,29]      # Frequency range of spectral events
fVec = np.arange(1,60+1)            # Vector of fequency values over which to calculate TFR
Fs = 600                            # Sampling rate of time-series
findMethod = 1                      # Event-finding method (1 allows for maximal overlap while 2 limits overlap in each respective suprathreshold region)
width = 7

thrFOM = 6; #Factors of Median threshold (see Shin et al. eLife 2017 for details concerning this value)

footprintFreq = 4
footprintTime = 80
threshold = 0.00
neighbourhood_size = (footprintFreq,footprintTime)

vis = True

################################
# Processing starts here
subjectIDs = np.arange(numSubj)+1

# Load data sessions/subjects from the same experimental setup so that 
# spectral event features are differentially characterized only between the 
# desired trial classification labels: in this case, detection vs. 
# non-detection
#
# Note: each .mat file contains: 
#   'prestim_raw_yes_no' - 200 trials x 600 time samples matrix of time series data
#   'YorN' - 200 trials x 1 matrix of 1s or 0s to indicate trial label
x = []
for s in subjectIDs:
    testFile = os.path.join('test_data', "".join(['prestim_humandetection_600hzMEG_subject',
        str(s), '.mat']))
    a = io.loadmat(testFile)
    x.append( a )

numTrials, numSamples = a['prestim_raw_yes_no'].shape

# Validate fVec input
Fn = Fs/2                   # Nyquist frequency
dt = 1/Fs                   # Sampling time interval
Fmin = 1/(numSamples*dt)    # Minimum resolvable frequency

if fVec[0] < Fmin: 
    sys.exit('Frequency vector includes values outside the resolvable/alias-free range.')
elif fVec[-1] > Fn: 
    sys.exit('Frequency vector includes values outside the resolvable/alias-free range.')
elif np.abs(fVec[1]-fVec[0]) < Fmin:
    sys.exit('Frequency vector includes values outside the resolvable/alias-free range.')

# Run spectral event analysis per dataset (Dataset loop)
TFR = []
specEvents = []
ctr = 1
for thisX in x:

    # Convert data to TFR 
    thisData = thisX['prestim_raw_yes_no']
    thisClassLabels = thisX['YorN']
    thisTFR, tVec, fVec = tse.spectralevents_ts2tfr( thisData.T, fVec, Fs, width )
    TFR.append( thisTFR )

    # Normalize the TFR data [tr x f x t] to the median value per frequency band 
    numTrials, numFreqBins, numSamples = thisTFR.shape
    TFR_order = np.transpose(thisTFR, axes=[1,0,2]) # [f x tr x t]
    TFR_reshape = np.reshape(TFR_order, (numFreqBins, numTrials*numSamples))
    TFRmeds = np.median(TFR_reshape, axis=1)        # f vector
    TFRmeds_expanded = np.transpose(np.tile(TFRmeds, (numSamples,numTrials,1)), axes=[1,2,0])
    thisTFR_norm = thisTFR/TFRmeds_expanded

    # Find local maxima in TFR
    thisSpecEvents = tse.spectralevents_find (findMethod, thrFOM, tVec, fVec, thisTFR, thisClassLabels, 
        neighbourhood_size, threshold, Fs)
    thisSpecEvents = pd.DataFrame( thisSpecEvents )
    specEvents.append( thisSpecEvents )

    # Extract event attributes for this test data 
    classes = np.unique( thisSpecEvents['Hit/Miss'].tolist() )

    # Plot results?
    if vis:

        # Plot results for each class of trial
        for clss in classes:

            # Get TFR, time course, and trial IDs for this class of trials only 
            trial_inds = np.where(thisClassLabels == clss)[0]
            classTFR = thisTFR[trial_inds,:,:]
            classTFR_norm = thisTFR_norm[trial_inds,:,:]
            classData = thisData[trial_inds,:]

            # Get events data for this class only, and update trial indices to be consecutive
            #   starting at 0
            df = thisSpecEvents[thisSpecEvents['Hit/Miss']==clss]
            classEvents = df.copy()
            classEvents = classEvents.replace(trial_inds, np.arange(len(trial_inds)))

            # Drop events that are low threshold (below 6 FOM)
            classEvents = classEvents[classEvents['Outlier Event']==1]

            # Make figure
            fig, axs = tse.spectralevents_vis( classEvents, classData, classTFR, classTFR_norm, 
                tVec, fVec, eventBand )
            # Add title
            axs[0,0].set_title( 'DataSet ' + str(ctr) + ', Trial class ' + str(clss) )

            # Save figure
            figName = os.path.join('test_results', 'python', 
                    "".join(['prestim_humandetection_600hzMEG_subject', str(ctr), '_class_', 
                    str(clss), '.png']))
            fig.savefig(figName)
            plt.close()
    
    ctr = ctr + 1
    




