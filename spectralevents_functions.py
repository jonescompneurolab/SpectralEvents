'''Spectral Event Analysis functions'''

# Authors: Tim Bardouille <tim.bardouille@dal.ca>
#          Ryan Thorpe <ryvthorpe@gmail.com>

import numpy as np
import scipy.signal as signal
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt


def spectralevents_ts2tfr(S, freqs, Fs, width):
    '''
    Calculates the TFR (in spectral power) of a time-series waveform by 
    convolving in the time-domain with a Morlet wavelet.                            

    Input
    -----
    S    : signals = Trials x Time
    freqs    : frequencies over which to calculate TF energy        
    Fs   : sampling frequency
    width: number of cycles in wavelet (> 5 advisable)  

    Output
    ------
    t    : time
    f    : frequency
    B    : phase-locking factor = frequency x time

    Adapted from Ole Jensen's traces2TFR in the 4Dtools toolbox.

    See also SPECTRALEVENTS, SPECTRALEVENTS_FIND, SPECTRALEVENTS_VIS.
    '''

    S = np.array(S)
    n_trials = S.shape[0]
    n_samps = S.shape[1]
    n_freqs = len(freqs)
    # Validate freqs input

    Fn = Fs / 2  # Nyquist frequency
    dt = 1 / Fs  # Sampling time interval
    min_freq = 1 / (n_samps * dt)  # Minimum resolvable frequency

    if freqs[0] < min_freq:
        raise ValueError('Frequency vector includes values outside the resolvable/alias-free range.')
    elif freqs[-1] > Fn:
        raise ValueError('Frequency vector includes values outside the resolvable/alias-free range.')
    elif np.abs(freqs[1] - freqs[0]) < min_freq:
        raise ValueError('Frequency vector includes values outside the resolvable/alias-free range.')

    TFR = np.zeros((n_trials, n_freqs, n_samps))

    # Trial Loop
    for trial_idx in np.arange(n_trials):
        # Frequency loop
        for freq_idx in np.arange(n_freqs):
            TFR[trial_idx, freq_idx, :] = energyvec(freqs[freq_idx], signal.detrend(S[trial_idx, :]), Fs, width)

    return TFR


def tfr_normalize(tfr):
    '''Normalize the power in each frequency band of a time-frequency response

    Parameters
    ----------
    tfr : array, shape ([n_trials,] n_freqs, n_times)
        The time-frequency response (TFR) to be normalized.

    Returns
    -------
    tfr_norm : array
        The normalized TFR calculated by dividing the power values in each
        frequency bin by the median power across all trials and time samples.
    '''

    if len(tfr.shape) == 3:
        n_trials, _, n_times = tfr.shape
        med_powers = np.median(tfr, axis=(0, 2))
        med_powers_tiled = np.tile(med_powers, (n_trials, n_times, 1))
        med_powers = np.transpose(med_powers_tiled, axes=(0, 2, 1))

    elif len(tfr.shape) == 2:
        _, n_times = tfr.shape
        med_powers = np.median(tfr, axis=1)
        med_powers_tiled = np.tile(med_powers, (n_times, 1))
        med_powers = np.transpose(med_powers_tiled, axes=(1, 0))

    else:
        raise ValueError(f'TFR must be an array of at least 2 dimensions. Got '
                         f'{tfr.shape}')

    return tfr / med_powers


def find_events(event_band, threshold_fom, times, freqs, TFR, Fs):
    '''
    SPECTRALEVENTS_FIND Algorithm for finding and calculating spectral 
      events on a trial-by-trial basis of of a single subject/session. Uses 
      one of three methods before further analyzing and organizing event 
      features:

      1) (Primary event detection method in Shin et al. eLife 2017): Find 
          spectral events by first retrieving all local maxima in 
          un-normalized TFR using imregionalmax, then selecting suprathreshold
          peaks within the frequency band of interest. This method allows for 
          multiple, overlapping events to occur in a given suprathreshold 
          region and does not guarantee the presence of within-band, 
          suprathreshold activity in any given trial will render an event.


    specEv_struct = spectralevents_find(event_band,threshold,times,freqs,TFR)

    Inputs:
      event_band - range of frequencies ([Fmin_event Fmax_event]; Hz) over 
          which above-threshold spectral power events are classified.
      thrFOM - factors of median threshold; positive real number used to
          threshold local maxima and classify events (see Shin et al. eLife 
          2017 for discussion concerning this value).
      times - time vector (s) over which the time-frequency response (TFR) is 
          calcuated.
      freqs - frequency vector (Hz) over which the time-frequency response 
          (TFR) is calcuated.
      TFR - time-frequency response (TFR) (trial-frequency-time) for a
          single subject/session.

    Outputs:
      specEv_struct - event feature structure with three main sub-structures:
          TrialSummary (trial-level features), Events (individual event 
          characteristics), and IEI (inter-event intervals from all trials 
          and those associated with only a given class label).
    '''

    # Initialize general data parameters
    # Number of elements in discrete frequency spectrum
    flength = TFR.shape[1]
    # Number of point in time
    tlength = TFR.shape[2]
    # Number of trials
    numTrials = TFR.shape[0]

    # Median power at each frequency across all trials
    TFRpermute = np.transpose(TFR, [1, 2, 0])  # freq x time x trial
    TFRreshape = np.reshape(TFRpermute, (flength, tlength * numTrials))
    medianPower = np.median(TFRreshape, axis=1)

    # Spectral event threshold for each frequency value
    eventThresholdByFrequency = threshold_fom*medianPower

    # Validate consistency of parameter dimensions
    if flength != len(freqs):
        raise ValueError('Mismatch in frequency dimensions!')
    if tlength != len(times):
        raise ValueError('Mismatch in time dimensions!')

    # Find spectral events using appropriate method
    #    Implementing find_method=1 for now
    spectralEvents = find_localmax_method_1(TFR, freqs, times, event_band,
                                            eventThresholdByFrequency,
                                            medianPower, Fs)

    return spectralEvents


def energyvec(f, s, Fs, width):
    '''
    Return a vector containing the energy as a
    function of time for frequency f. The energy
    is calculated using Morlet's wavelets. 
    s : signal
    Fs: sampling frequency
    width : width of Morlet wavelet (>= 5 suggested).
    '''

    dt = 1 / Fs
    sf = f / width
    st = 1 / (2 * np.pi * sf)

    t = np.arange(-3.5 * st, 3.5 * st, dt)
    m = morlet(f, t, width)

    y = np.convolve(s, m)
    y = 2 * (dt * np.abs(y)) ** 2
    lowerLimit = int(np.ceil(len(m) / 2))
    upperLimit = int(len(y) - np.floor(len(m) / 2) + 1)
    y = y[lowerLimit:upperLimit]

    return y


def morlet(f, t, width):
    '''
    Morlet's wavelet for frequency f and time t. 
    The wavelet will be normalized so the total energy is 1.
    width defines the ``width'' of the wavelet. 
    A value >= 5 is suggested.

    Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
    '''

    sf = f / width
    st = 1 / (2 * np.pi * sf)
    A = 1 / (st * np.sqrt(2 * np.pi))
    y = A * np.exp(-t ** 2 / (2 * st ** 2)) * np.exp(1j * 2 * np.pi * f * t)

    return y


def fwhm_lower_upper_bound1(vec, peakInd, peakValue):
    '''
    Function to find the lower and upper indices within which the vector is less than the FWHM
      with some rather complicated boundary rules (Shin, eLife, 2017)
    '''

    halfMax = peakValue/2

    # Extract data before the peak only (data should be rising at the end of the new array)
    vec1 = vec[0:peakInd]
    # Find indices less than half the max
    vec1_underThreshold = np.where(vec1 < halfMax)[0]
    if len(vec1_underThreshold) == 0:
        # There are no indices less than half the max, so we have to estimate the lower edge
        estimateLowerEdge = True
    else:
        # There are indices less than half the max, take the last one under halfMax as the lower edge
        estimateLowerEdge = False
        lowerEdgeIndex = vec1_underThreshold[-1]

    # Extract data following the peak only (data should be falling at the start of the new array)
    vec2 = vec[peakInd:]
    # Find indices less than half the max
    vec2_underThreshold = np.where(vec2 < halfMax)[0]
    if len(vec2_underThreshold) == 0:
        # There are no indices less than half the max, so we have to estimate the upper edge
        estimateUpperEdge = True
    else:
        # There are indices less than half the max, take the first one under halfMax as the upper edge
        estimateUpperEdge = False
        upperEdgeIndex = vec2_underThreshold[0] + len(vec1)

    if not estimateLowerEdge:
        if not estimateUpperEdge:
            # FWHM fits in the range, so pick off the edges of the FWHM
            lowerInd = lowerEdgeIndex
            upperInd = upperEdgeIndex
            FWHM = upperInd - lowerInd
        if estimateUpperEdge:
            # FWHM fits in on the low end, but hits the edge on the high end
            lowerInd = lowerEdgeIndex
            upperInd = len(vec)-1
            FWHM = 2 * (peakInd - lowerInd + 1)
    else:
        if not estimateUpperEdge:
            # FWHM hits the edge on the low end, but fits on the high end
            lowerInd = 0
            upperInd = upperEdgeIndex
            FWHM = 2 * (upperInd - peakInd + 1)
        if estimateUpperEdge:
            # FWHM hits the edge on the low end and the high end
            lowerInd = 0
            upperInd = len(vec)-1
            FWHM = 2*len(vec)

    return lowerInd, upperInd, FWHM


def find_localmax_method_1(TFR, freqs, times, event_band,
                           eventThresholdByFrequency,
                           medianPower, Fs):
    '''
    1st event-finding method (primary event detection method in Shin et
    al. eLife 2017): Find spectral events by first retrieving all local
    maxima in un-normalized TFR using imregionalmax, then selecting
    suprathreshold peaks within the frequency band of interest. This
    method allows for multiple, overlapping events to occur in a given
    suprathreshold region and does not guarantee the presence of
    within-band, suprathreshold activity in any given trial will render
    an event.

    spectralEvents: 12 column matrix for storing local max event metrics:
            trial index,            hit/miss,         maxima frequency,
            lowerbound frequency,     upperbound frequency,
            frequency span,         maxima timing,     event onset timing,
            event offset timing,     event duration, maxima power,
            maxima/median power
    '''
    # Number of elements in discrete frequency spectrum
    flength = TFR.shape[1]
    # Number of point in time
    tlength = TFR.shape[2]
    # Number of trials
    numTrials = TFR.shape[0]

    spectralEvents = list()

    # Retrieve all local maxima in TFR using python equivalent of imregionalmax
    for ti in range(numTrials):

        # Get TFR data for this trial [frequency x time]
        thisTFR = TFR[ti, :, :]

        # Find local maxima in the TFR data
        data = thisTFR
        # Find maximum amoung adjacent pixels (3x3 footprint) for each pixel
        data_max = filters.maximum_filter(data, size=(3, 3))
        maxima = (data == data_max)
        data_min = filters.minimum_filter(data, size=(3, 3))
        # Rule out pixels with footprints that have flatlined
        maxima[data_max == data_min] = False
        labeled, num_objects = ndimage.label(maxima)
        yx = ndimage.center_of_mass(data, labels=labeled,
                                    index=range(1, num_objects + 1))

        #numPeaks = len(yx)

        peakF = list()
        peakT = list()
        peakPower = list()
        for f_idx, t_idx in yx:
            f_idx = int(round(f_idx))
            t_idx = int(round(t_idx))
            event_freq = freqs[f_idx]
            if (event_freq >= event_band[0] and event_freq <= event_band[1] and
               thisTFR[f_idx, t_idx] > eventThresholdByFrequency[f_idx]):
                peakF.append(f_idx)
                peakT.append(t_idx)
                peakPower.append(thisTFR[f_idx, t_idx])

        numPeaks = len(peakF)

        # Find local maxima lowerbound, upperbound, and full width at half max
        #    for both frequency and time
        for lmi in range(numPeaks):
            thisPeakF = peakF[lmi]
            thisPeakT = peakT[lmi]
            thisPeakPower = peakPower[lmi]

            # Indices of TFR frequencies < half max power at the time of a given local peak
            TFRFrequencies = thisTFR[:, thisPeakT]
            lowerInd, upperInd, FWHM = fwhm_lower_upper_bound1(TFRFrequencies,
                                                               thisPeakF, thisPeakPower)
            lowerEdgeFreq = freqs[lowerInd]
            upperEdgeFreq = freqs[upperInd]
            FWHMFreq = FWHM * (freqs[1] - freqs[0])

            # Indices of TFR times < half max power at the frequency of a given local peak
            TFRTimes = thisTFR[thisPeakF, :]
            lowerInd, upperInd, FWHM = fwhm_lower_upper_bound1(TFRTimes,
                                                               thisPeakT, thisPeakPower)
            lowerEdgeTime = times[lowerInd]
            upperEdgeTime = times[upperInd]
            FWHMTime = FWHM / Fs

            # Put peak characteristics to a dictionary
            #        trial index,            hit/miss,         maxima frequency,
            #        lowerbound frequency,     upperbound frequency,
            #        frequency span,         maxima timing,     event onset timing,
            #        event offset timing,     event duration, maxima power,
            #        maxima/median power
            peakParameters = {
                'Trial': ti,
                'Peak Frequency': freqs[thisPeakF],
                'Lower Frequency Bound': lowerEdgeFreq,
                'Upper Frequency Bound': upperEdgeFreq,
                'Frequency Span': FWHMFreq,
                'Peak Time': times[thisPeakT],
                'Event Onset Time': lowerEdgeTime,
                'Event Offset Time': upperEdgeTime,
                'Event Duration': FWHMTime,
                'Peak Power': thisPeakPower,
                'Normalized Peak Power': thisPeakPower / medianPower[thisPeakF]
            }

            # Build a list of dictionaries
            spectralEvents.append(peakParameters)

    return np.array(spectralEvents)


def plot_events(TFR, times, freqs, event_band, spec_events=None,
                timeseries=None, ax=None, vlim=None, ylim_ts=None, label=None):

    if TFR.shape != (len(freqs), len(times)):
        raise ValueError(f'TFR must be an array of shape (n_freqs, n_times) '
                         f'got TFR: {TFR.shape}, freqs: ({len(freqs)},), '
                         f'times: ({len(times)},)')

    if vlim is None:
        vlim = [None, None]

    # convert to numpy array if not already
    freqs = np.array(freqs)

    # frequencies within the band of interest
    # band_mask = np.logical_and(freqs >= event_band[0], freqs <= event_band[1])

    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        fig = ax.get_figure()

    # plot TFR
    im = ax.pcolormesh(times, freqs, TFR, cmap='jet', vmin=vlim[0],
                       vmax=vlim[1], shading='nearest')
    fig.colorbar(im, ax=ax)
    ax.axhline(y=event_band[0], c='w', linewidth=2., linestyle=':',
               alpha=.7)
    ax.axhline(y=event_band[1], c='w', linewidth=2., linestyle=':',
               alpha=.7)
    ax.set_yticks([freqs[0], event_band[0], event_band[1], freqs[-1]])
    ax.set_xlim(times[0], times[-1])

    # overlay with timecourse
    # twin axis needs to be created and yticks set regardless of if a
    # timeseries is provided to ensure consistency with seaborn tick formatting
    ax_twin = ax.twinx()
    ax_twin.set_yticks([])
    if timeseries is not None:
        ax_twin.plot(times, timeseries, 'w', linewidth=1., alpha=0.8)

        if ylim_ts is not None:
            ax_twin.set_ylim(ylim_ts[0], ylim_ts[1])

    # plot event locations
    if spec_events is not None:
        event_times = [event['Peak Time'] for event in spec_events]
        event_freqs = [event['Peak Frequency'] for event in spec_events]
        # reverse sigmoid: make scatter markers more transparent when there are
        # more of them
        alpha = lambda x : (0.6 + 0.4 * np.exp(-0.2 * (x - 30))  # noqa
                            / (1 + np.exp(-0.2 * (x - 30))))  # noqa
        ax.scatter(event_times, event_freqs, s=20, c='w', marker='x',
                   alpha=alpha(len(spec_events)))

    # add label on top of spectrogram plot
    if label is not None:
        ax.annotate(label, xy=(0.01, 0.95), xycoords='axes fraction',
                    va='top', ha='left', color='w', size=10, fontweight='bold')

    return fig


def plot_avg_spectrogram(TFR, times, freqs, event_band, spec_events=None,
                         timeseries=None, example_epochs=None, vlim=None):
    '''
    Function to plot spectral events on test data (to check against Matlab code)

    spec_events = spectral event characteristics 
    timeseries = trials x time electrophysiological data
    TFR = trials x frequency x time TFR of timeseries
    TFR = trials x frequency x time normalized TFR of timeseries
    times = vector of time samples
    freqs = vector of frequency bins
    event_band = vector with min and max frequency for spectral event mapping
    '''

    if example_epochs is not None:
        trial_idx_set = set(range(TFR.shape[0]))
        trial_idx_subset = set(example_epochs)
        if trial_idx_subset.intersection(trial_idx_set) != trial_idx_subset:
            raise ValueError('One or more of the specified example trial '
                             'indices does not exist in the provided TFR '
                             'array.')
    else:
        # set to empty list
        example_epochs = list()

    fig, axs = plt.subplots(nrows=len(example_epochs) + 1, ncols=1,
                            sharex=True)

    # plot trial-average TFR
    TFR_avg = np.mean(TFR, axis=0).squeeze()
    plot_events(TFR=TFR_avg, times=times, freqs=freqs,
                event_band=event_band, spec_events=spec_events,
                ax=axs[0], vlim=vlim, label='trial avg.')

    # plot TFR + events for example trials
    if timeseries is not None and example_epochs is not None:
        max_ts_amplitude = np.max(timeseries[example_epochs])
        min_ts_amplitude = np.min(timeseries[example_epochs])
        ylim_ts = [max_ts_amplitude, min_ts_amplitude]

        for count_idx, trial_idx in enumerate(example_epochs):
            # get spectral events for the current trial
            trial_events = [event for event in spec_events
                            if event['Trial'] == trial_idx]

            # plot trial TFR
            TFR_trial = TFR[trial_idx, :, :].squeeze()
            timeseries_trial = timeseries[trial_idx, :]

            plot_events(TFR=TFR_trial, times=times, freqs=freqs,
                        event_band=event_band, spec_events=trial_events,
                        timeseries=timeseries_trial, ax=axs[count_idx + 1],
                        vlim=vlim, ylim_ts=ylim_ts, label=f'epoch {trial_idx}')

    axs[-1].set_xlabel('time (s)')
    axs[0].set_ylabel('freq. (Hz)')

    return fig
