'''Code Tests'''

# Authors: Ryan Thorpe <ryvthorpe@gmail.com>
# Authors: Darcy Diesburg <darcy.diesburg@gmail.com>

import sys
import os
import os.path as op

import numpy as np
from scipy.io import loadmat

# sys.path.append('/Users/darcy/Desktop/SpectralEvents/')
import spectralevents as se

def test_event_comparison():
    '''Test whether the output from MATLAB and Python event detection scripts align

    Parameters
    ----------
    method : float
    The event detection method (1, 2, or 3) used to detect spectral events.

    Notes
    -----
    Currently only supports testing method 1 and 3, with method 2 to be added
    to the Python event detection functions in the future. Python output is
    tested against the output of MATLAB scripts from Shin et al., 2017. Detection
    of the same number of beta events within each trial is required, and tolerable
    margins of deviation for the other parameters.
    '''

    # Define tolerable deviance for spectral event characteristics
    dev_max_time = 5  # +-5ms allowed for timing latency diff
    dev_count_avg = 0.001 # 1 of 1000 trials, on average, may have different event count

    # Load MATLAB data
    data_dir = os.getcwd() 
    matlab_ev_struct = loadmat(
        op.join(data_dir, 'tests', 'beta_events_shin_2017.mat'))

    # Get python events by running find_events on demo data
    subj_ids = [str(id) for id in range(1, 10 + 1)]  # subject IDs 1-10

    # set parameters
    samp_freq = 600
    n_times = 600
    # fequency values (Hz) over which to calculate TFR
    freqs = list(range(1, 60 + 1))
    times = np.arange(n_times) / samp_freq  # seconds
    event_band = [15, 29]  # beta band (Hz)
    thresh_fom = 6.0  # factor-of-the-median threshold

    # get python evs by running find_events()
    py_ev_dict = list()
    for _, id_val in enumerate(subj_ids):
        fname = op.join(
            data_dir, 'data','prestim_humandetection_600hzMEG_subject' + id_val + '.mat')
        raw_data = loadmat(fname)
        data = raw_data['prestim_raw_yes_no']  # MEG time series (trials x samples)
        # calculate TFR
        tfrs = se.tfr(data, freqs, samp_freq)

        # find events
        subj_evs = se.find_events(tfr=tfrs, times=times, freqs=freqs, event_band=event_band,
                                thresholds=None, threshold_FOM=thresh_fom)
        py_ev_dict.append(subj_evs)

    # take ev counts and timing latencies from respective dictionaries and structures 
    # to put into a subj x trial matrix for comparison
    py_ev_count_mat = np.zeros((10,200))
    for subj_idx, id_val in enumerate(subj_ids):
        py_ev_count_mat[subj_idx,:] = np.array([len(trial) for trial in py_ev_dict[subj_idx]])

    # matlab - these are already in an array, so just reshape
    matlab_ev_count_mat = np.array([trial.shape[1] for trial in matlab_ev_struct['event_times'][:, 0]])
    matlab_ev_count_mat = np.reshape(matlab_ev_count_mat,(10,200))

    # Overall check that same # of evs detected (or a tolerable # of differences)
    assert np.mean(np.abs(matlab_ev_count_mat-py_ev_count_mat))<=dev_count_avg, "Should be mostly same number of events across all trials"

    # Once assured number of evs are within tolerable comparison limits, extract event characteristics
    # matlab evs
    same_ev_count_bool = (matlab_ev_count_mat==py_ev_count_mat)
    same_ev_count_bool_reshape = np.reshape(same_ev_count_bool,(2000,))

    # get event timing in trials in which the event count was equal
    matlab_ev_timing_list = list()
    py_ev_timing_list = list()
    # extract latencies from matlab structure, which has all subject trials in array of arrays within dictionary key
    for trial in range(len(same_ev_count_bool_reshape)):
        if same_ev_count_bool_reshape[trial]:
            matlab_ev_timing_list = [np.append(matlab_ev_timing_list, np.array(event)) for event in matlab_ev_struct['event_times'][trial, 0]][0]
    # extract latencies from py dictionaries for each subject
    for subj_idx in range(10):
        subj_evs = py_ev_dict[subj_idx]
        for trial_idx in range(200):
            trial_evs = list()
            if same_ev_count_bool[subj_idx, trial_idx]:
                # append sorted event latencies
                if len(subj_evs[trial_idx])>0:
                    for event in subj_evs[trial_idx]:
                        trial_evs = np.append(trial_evs, event['Peak Time'])
                    sorted = np.argsort(trial_evs)
                    py_ev_timing_list = np.append(py_ev_timing_list, trial_evs[sorted])

    # ensure that, in the events that were the same detected in each method, timing is within a tolerable limit of difference
    assert ((matlab_ev_timing_list - py_ev_timing_list)*1000 <= dev_max_time).all(
        ), "Timing of events should be within tolerable limits of deviation"

    pass
