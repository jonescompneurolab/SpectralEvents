'''Code Tests'''

# Authors: Ryan Thorpe <ryvthorpe@gmail.com>
# Authors: Darcy Diesburg <darcy.diesburg@gmail.com>

import sys
import os
import os.path as op

import pytest
import numpy as np
from scipy.io import loadmat

sys.path.append('/gpfs/home/ddiesbur/SpectralEvents')
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
    dev_max_time = 5  # +-5ms
    dev_max_freq = 1  # +-1Hz
    dev_low_freq = 1  # +-1Hz
    dev_high_freq = 1  # +-1Hz
    dev_onset = 10  # +-10ms
    dev_offset = 10  # +-10ms
    dev_fom_pow = 0.5  # +-0.5 FOM

    # Load MATLAB data
    data_dir = os.getcwd()
    if method == 1:
        matlab_file = loadmat(
            op.join(data_dir, 'tests', 'beta_events_shin_2017.mat'))
    #elif method == 3:
    #    matlab_file = loadmat(
    #        op.join(data_dir, 'tests', 'beta_events_shin_2017_method3.mat'))
    else:
        raise ValueError('Unsupported method! Please use 1') #or 3')
    mat_all_evs = matlab_file['event_times'][:, 0]

    # Get python events by running find_events on demo data
    # <----------SIMPLIFY BY GETTING PRE-SAVED PICKLE FILES FOR THIS INSTEAD----------->
    subj_ids = [str(id) for id in range(1, 10 + 1)]  # subject IDs 1-10

    # set parameters
    samp_freq = 600
    n_times = 600
    # fequency values (Hz) over which to calculate TFR
    freqs = list(range(1, 60 + 1))
    times = np.arange(n_times) / samp_freq  # seconds
    event_band = [15, 29]  # beta band (Hz)
    thresh_fom = 6.0  # factor-of-the-median threshold

    # get python evs
    py_all_evs = list()
    for _, id_val in enumerate(subj_ids):
        fname = op.join(
            data_dir, 'data','prestim_humandetection_600hzMEG_subject' + id_val + '.mat')
        raw_data = loadmat(fname)
        data = raw_data['prestim_raw_yes_no']  # MEG time series (trials x samples)
        # calculate TFR
        tfrs = se.tfr(data, freqs, samp_freq)

        # find events
        subj_evs = se.find_events(tfr=tfrs, times=times, freqs=freqs, event_band=event_band,
                                  thresholds=None, threshold_FOM=thresh_fom, find_method=method)
        py_all_evs.append(subj_evs)

    # put from cell array/dict into sets
    py_evs = np.array([])
    for subj_idx, id_val in enumerate(subj_ids):
        py_evs = np.append(py_evs,np.array([len(trial) for trial in py_all_evs[subj_idx]]))
    mat_evs = np.array([trial.shape[1] for trial in mat_all_evs])

    # Overall check that same # of evs detected (or we adapt this to a tolerable # of difference)
    assert np.array_equal(
        mat_evs, py_evs), "Should be same number of events across all trials"

    # Once assured same num of evs, extract event characteristics into np arrays
    # matlab evs
    mat_times = np.array([])
    mat_max_freq = np.array([])
    mat_low_freq = np.array([])
    mat_high_freq = np.array([])
    mat_onsets = np.array([])
    mat_offsets = np.array([])
    mat_power = np.array([])
    trial_ind = -1
    for trial in mat_all_evs:
        trial_ind += 1
        mat_times = [np.append(mat_times, event) for event in trial]
        mat_max_freq = [np.append(mat_max_freq, event)
                      for event in matlab_file['event_maxfreqs'][trial_ind]]
        mat_low_freq = [np.append(mat_low_freq, event)
                      for event in matlab_file['event_lowfreqs'][trial_ind]]
        mat_high_freq = [np.append(mat_high_freq, event)
                       for event in matlab_file['event_highfreqs'][trial_ind]]
        mat_onsets = [np.append(mat_onsets, event)
                     for event in matlab_file['event_onsets'][trial_ind]]
        mat_offsets = [np.append(mat_offsets, event)
                      for event in matlab_file['event_offsets'][trial_ind]]
        mat_power = [np.append(mat_power, event)
                    for event in matlab_file['event_powers'][trial_ind]]

    # python evs
    py_times = np.array([])
    py_onsets = np.array([])
    py_offsets = np.array([])
    py_max_freq = np.array([])
    py_low_freq = np.array([])
    py_high_freq = np.array([])
    py_power = np.array([])
    for trial in py_all_evs:
        if len(trial) > 0:
            trial_times = np.array([])
            trial_onsets = np.array([])
            trial_offsets = np.array([])
            trial_max_freq = np.array([])
            trial_low_freq = np.array([])
            trial_high_freq = np.array([])
            trial_power = np.array([])
            for event in trial:
                trial_times = np.append(trial_times, event['Peak Time'])
                trial_onsets = np.append(trial_onsets, event['Event Onset Time'])
                trial_offsets = np.append(trial_offsets, event['Event Offset Time'])
                trial_max_freq = np.append(trial_max_freq, event['Peak Frequency'])
                trial_low_freq = np.append(
                    trial_low_freq, event['Lower Frequency Bound'])
                trial_high_freq = np.append(
                    trial_high_freq, event['Upper Frequency Bound'])
                trial_power = np.append(trial_power, event['Normalized Peak Power'])
                sorter = np.argsort(trial_times)
            py_times = np.append(py_times, trial_times[sorter])
            py_onsets = np.append(py_onsets, trial_onsets[sorter])
            py_offsets = np.append(py_offsets, trial_offsets[sorter])
            py_max_freq = np.append(py_max_freq, trial_max_freq[sorter])
            py_low_freq = np.append(py_low_freq, trial_low_freq[sorter])
            py_high_freq = np.append(py_high_freq, trial_high_freq[sorter])
            py_power = np.append(py_power, trial_power[sorter])

    # Timing - check max timing
    assert ((mat_times - py_times)*1000 <= dev_max_time).all(
    ), "Timing of events should be within tolerable limits of deviation"

    # Frequency - check max freq
    assert ((mat_max_freq - py_max_freq) <= dev_max_freq).all(
    ), "Peak frequency of events should be within tolerable limits of deviation"

    # check lower freq bound
    assert ((mat_low_freq - py_low_freq) <= dev_low_freq).all(
    ), "Lower frequency of events should be within tolerable limits of deviation"

    # check higher freq bound
    assert ((mat_high_freq - py_high_freq) <= dev_high_freq).all(
    ), "Higher frequency of events should be within tolerable limits of deviation"

    # Duration - check onset
    assert ((mat_onsets - py_onsets)*1000 <= dev_onset).all(
    ), "Onset of events should be within tolerable limits of deviation"

    # check offset
    assert ((mat_offsets - py_offsets)*1000 <= dev_offset).all(
    ), "Offset of events should be within tolerable limits of deviation"

    # Power - check maximum power in FOM
    assert ((mat_power - py_power) <= dev_fom_pow).all(
    ), "Power of events should be within tolerable limits of deviation"

    pass
