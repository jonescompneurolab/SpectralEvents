'''Code Tests'''

# Authors: Ryan Thorpe <ryvthorpe@gmail.com>
#          Darcy Diesburg <darcy.diesburg@gmail.com>

import os.path as op

import numpy as np
from scipy.io import loadmat

import SpectralEvents.spectralevents as se


def test_event_comparison():
    '''Test if the output from MATLAB and Python event detection scripts align

    Python output is tested against the output of MATLAB scripts from Shin et
    al., 2017. Detection of the same number of beta events within each trial
    is required, and tolerable margins of deviation for the other parameters.
    '''

    # number of subjects, trials in demo data
    n_subjs = 10
    n_trials = 200

    # Define tolerable deviance for spectral event characteristics
    dev_event_time = 2.0  # +-2ms allowed for timing latency diff
    dev_count_avg = 0.0  # 0/1000 trials, on average

    # Load MATLAB data
    root_dir = op.dirname(se.__file__)
    matlab_ev_struct = loadmat(op.join(root_dir, 'tests',
                                       'beta_events_shin_2017.mat'))

    # Get python events by running find_events on demo data
    subj_ids = [str(id) for id in range(1, n_subjs + 1)]  # subject IDs 1-10

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
            root_dir, 'data', 'prestim_humandetection_600hzMEG_subject' +
            id_val + '.mat')
        raw_data = loadmat(fname)
        data = raw_data['prestim_raw_yes_no']
        # calculate TFR
        tfrs = se.tfr(data, freqs, samp_freq)

        # find events
        subj_evs = se.find_events(tfr=tfrs, times=times, freqs=freqs,
                                  event_band=event_band, thresholds=None,
                                  threshold_FOM=thresh_fom)
        py_ev_dict.append(subj_evs)

    # take ev counts and timing latencies from respective dictionaries
    # and structures to put into a subj x trial matrix for comparison
    py_ev_count_mat = np.zeros((n_subjs, n_trials))
    for subj_idx, id_val in enumerate(subj_ids):
        py_ev_count_mat[subj_idx, :] = [
            len(trial) for trial in py_ev_dict[subj_idx]]

    # matlab - these are already in an array, so just reshape
    matlab_ev_count_mat = [
        trial.shape[1] for trial in matlab_ev_struct['event_times'][:, 0]]
    matlab_ev_count_mat = np.reshape(matlab_ev_count_mat, (n_subjs, n_trials))

    # Overall check that same # of evs detected
    # (or a tolerable # of differences)
    assert np.mean(np.abs(
        matlab_ev_count_mat - py_ev_count_mat
        )) <= dev_count_avg, "Should be mostly same number of evs"

    # Once assured number of evs are within tolerable limits, extract evs
    # matlab evs
    same_ev_count_bool = (matlab_ev_count_mat == py_ev_count_mat)
    same_ev_count_bool_reshape = np.reshape(
        same_ev_count_bool, (n_subjs*n_trials,))

    # get event timing in trials in which the event count was equal
    matlab_ev_timing_list = list()
    py_ev_timing_list = list()
    # extract latencies from matlab structure, which has all subject trials
    # in array of arrays within dictionary key
    # first, extract the events
    trial_indices = np.where(same_ev_count_bool_reshape)[0]
    event_times = [
        matlab_ev_struct['event_times'][trial, 0].flatten()
        for trial in trial_indices
    ]

    # then, concatenate the event times into a single list
    matlab_ev_timing_list = np.concatenate(event_times)

    # extract latencies from py dictionaries for each subject
    # First, create a list of sorted event times for each trial
    trial_event_times = []
    for subj_idx, subj_evs in enumerate(py_ev_dict):
        for trial_idx in range(n_trials):
            if same_ev_count_bool[subj_idx, trial_idx]:
                event_times = [event['Peak Time']
                               for event in subj_evs[trial_idx]]
                sorted_event_times = np.sort(event_times)
                trial_event_times.append(sorted_event_times)

    # Then, concatenate the lists into a single array
    py_ev_timing_list = np.concatenate(trial_event_times)

    # ensure that, in the events that were the same detected in each method,
    # timing is within a tolerable limit of difference
    assert (
        abs(matlab_ev_timing_list - py_ev_timing_list)*1000 <= dev_event_time
    ).all(), "Timing of events should be within tolerable limits"
    pass
