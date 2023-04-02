'''Spectral Event Analysis functions'''

# Authors: Tim Bardouille <tim.bardouille@dal.ca>
#          Ryan Thorpe <ryvthorpe@gmail.com>

import numpy as np
import scipy.signal as signal
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt


def tfr(timeseries, freqs, samp_freq, width=7.):
    '''Calculate the time-frequency response by convolving with Morlet wavelets

    Parameters
    ----------
    timeseries : array-like, shape ([n_epochs,] n_times)
        The timeseries signal, one for each epoch of data.
    freqs : array-like, shape (n_freqs,)
        Frequency domain of the TFR (Hz).
    samp_freq : float
        Sampling frequency (Hz).
    width : int
        Number of cycles in each Morlet wavelet (>5 advisable).

    Returns
    -------
    tfr : array, shape (n_epochs, n_freqs, n_times)
        A collection of time-frequency responses (TFRs), one for each epoch.

    Notes
    -----
    Adapted from Ole Jensen's traces2tfr in the 4Dtools toolbox.
    '''

    ts = np.atleast_2d(timeseries)
    n_epochs = ts.shape[0]
    n_samps = ts.shape[1]
    n_freqs = len(freqs)

    # Validate freqs input
    f_nyquist = samp_freq / 2  # Nyquist frequency
    dt = 1 / samp_freq  # Sampling time interval
    min_freq = 1 / (n_samps * dt)  # Minimum resolvable frequency

    if freqs[0] < min_freq:
        raise ValueError('Frequency vector includes values outside the '
                         'resolvable/alias-free range.')
    elif freqs[-1] > f_nyquist:
        raise ValueError('Frequency vector includes values outside the '
                         'resolvable/alias-free range.')
    elif np.abs(freqs[1] - freqs[0]) < min_freq:
        raise ValueError('Frequency vector includes values outside the '
                         'resolvable/alias-free range.')

    tfr = np.zeros((n_epochs, n_freqs, n_samps))

    # Trial Loop
    for trial_idx in np.arange(n_epochs):
        ts_detrended = signal.detrend(ts[trial_idx, :])
        # Frequency loop
        for freq_idx in np.arange(n_freqs):
            tfr[trial_idx, freq_idx, :] = _energyvec(freqs[freq_idx],
                                                     ts_detrended, samp_freq,
                                                     width)

    return tfr


def _get_power_thresholds(tfr, FOM_threshold=6.):
    '''Get the power threshold for each frequency band of a TFR'''

    med_powers = np.median(tfr, axis=1)
    return med_powers * FOM_threshold, med_powers


def tfr_normalize(tfr):
    '''Normalize the power in each frequency band of a TFR

    Parameters
    ----------
    tfr : array, shape ([n_epochs,] n_freqs, n_times)
        The time-frequency response (TFR) to be normalized.

    Returns
    -------
    tfr_norm : array
        The normalized TFR calculated by dividing the power values in each
        frequency bin by the median power across all trials and time samples.
    '''

    if len(tfr.shape) == 3:
        n_epochs, _, n_times = tfr.shape
        med_powers = np.median(tfr, axis=(0, 2))
        med_powers_tiled = np.tile(med_powers, (n_epochs, n_times, 1))
        med_powers = np.transpose(med_powers_tiled, axes=(0, 2, 1))

    elif len(tfr.shape) == 2:
        _, n_times = tfr.shape
        med_powers = np.median(tfr, axis=1)
        med_powers_tiled = np.tile(med_powers, (n_times, 1))
        med_powers = np.transpose(med_powers_tiled, axes=(1, 0))

    else:
        raise ValueError(f'TFR must be an array of at least 2 dimensions. Got '
                         f'{tfr.shape}.')

    return tfr / med_powers


def find_events(tfr, times, freqs, event_band, thresholds=None,
                threshold_FOM=6.):
    '''Locate spectral events in a time-frequency response.

    Parameters
    ----------
    tfr : array, shape ([n_epochs,] n_freqs, n_times)
        The time-frequency response (TFR) in which to search for high power
        spectral events.
    times : array-like, shape (n_times,)
        Time domain of the TFR in seconds.
    freqs : array-like, shape (n_freqs,)
        Frequency domain of the TFR in Hertz.
    event_band : list
        Lower and upper bounds (inclusive, Hz) of the frequency band-of-
        interest in which to search for spectral events.
    thresholds : None | array-like, shape (n_freqs,)
        If not None, these frequency-specific threshold values (in units of
        spectral power) will be used to identify suprathreshold spectral
        events.
    threshold_FOM : float | int
        Factor-of-the-median threshold (a.u.) with which to identify
        suprathreshold spectral events (default: 6). This threshold value is
        applied across frequencies of the TFR when thresholds=None. See Shin et
        al. eLife 2017 for more details concerning this value.

    Returns
    -------
    events : list of list of dict
        A nested list with n_epochs elements in the outer list and n_events in
        the inner list. Each inner element comprises an event that is
        characterized by dictionary items including it's location in time,
        location in frequency, duration, and frequency span, and more.

    Notes
    -----
    This version only supports find-method #1 at the moment. As outlined in
    Shin et al. eLife (2017), it isolates spectral events by first retrieving
    all local maxima in un-normalized TFR using imregionalmax, then selecting
    suprathreshold peaks within the frequency band of interest. This method
    allows for multiple, overlapping events to occur in a given suprathreshold
    region and does not guarantee that the presence of within-band,
    suprathreshold activity in any given trial will render an event.
    '''

    # ensure tfr has 3 dimensions (epochs x freqs x time)
    if len(tfr.shape) < 3:
        tfr = tfr[np.newaxis, ...]
    n_epochs = tfr.shape[0]
    n_freqs = tfr.shape[1]
    n_times = tfr.shape[2]

    # some time steps might be slightly different due to rounding error, etc.
    samp_freq = 1 / np.unique(np.diff(times).round(10))
    if len(samp_freq) > 1:
        raise ValueError('Sampling rate is not consistent across time '
                         'samples.')
    samp_freq = samp_freq[0]

    # concatenate trials together to make one big spectrogram, then find thresh
    tfr_permute = np.transpose(tfr, [1, 2, 0])  # freq x time x trial
    tfr_concat_trials = np.reshape(tfr_permute, (n_freqs, n_times * n_epochs))
    thresholds, med_powers = _get_power_thresholds(tfr_concat_trials,
                                                   FOM_threshold=threshold_FOM)

    # Validate consistency of parameter dimensions
    if n_freqs != len(freqs):
        raise ValueError('Mismatch in frequency dimensions!')
    if n_times != len(times):
        raise ValueError('Mismatch in time dimensions!')

    # Find spectral events using appropriate method
    #    Implementing find_method=1 for now
    events = _find_localmax_method_1(tfr, freqs, times, event_band,
                                     thresholds, med_powers, samp_freq)

    return events


def _energyvec(f, s, Fs, width=7.):
    '''Spectral energy as a function of frequency using Morlet wavelets.'''

    dt = 1 / Fs
    sf = f / width
    st = 1 / (2 * np.pi * sf)

    t_single_side = np.arange(0., 3.5 * st, dt)
    # mirror about 0; always an odd # of elements
    times = np.r_[-t_single_side[-1:0:-1], t_single_side]
    wavelet = _morlet(f, times, width)

    fourier_coeffs = np.convolve(s, wavelet)
    spec_energy = 2 * (dt * np.abs(fourier_coeffs)) ** 2
    lower_idx = int(len(wavelet) // 2)
    upper_idx = int(len(spec_energy) - (len(wavelet) // 2))

    return spec_energy[lower_idx:upper_idx]


def _morlet(f, t, width):
    '''
    Morlet's wavelet for frequency f and time t. The wavelet will be normalized
    so the total energy is 1. width defines the ``width'' of the wavelet. A
    value >= 5 is suggested.

    Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
    '''

    sf = f / width
    st = 1 / (2 * np.pi * sf)
    A = 1 / (st * np.sqrt(2 * np.pi))
    y = A * np.exp(-t ** 2 / (2 * st ** 2)) * np.exp(1j * 2 * np.pi * f * t)

    return y


def _fwhm_lower_upper_bound1(vec, peakInd, peakValue):
    '''
    Function to find the lower and upper indices within which the vector is
    less than the FWHM with some rather complicated boundary rules
    (Shin, eLife, 2017).
    '''

    halfMax = peakValue/2

    # Extract data before the peak only (data should be rising at the end of
    # the new array)
    vec1 = vec[0:peakInd]
    # Find indices less than half the max
    vec1_underThreshold = np.where(vec1 < halfMax)[0]
    if len(vec1_underThreshold) == 0:
        # There are no indices less than half the max, so we have to estimate
        # the lower edge
        estimateLowerEdge = True
    else:
        # There are indices less than half the max, take the last one under
        # halfMax as the lower edge
        estimateLowerEdge = False
        lowerEdgeIndex = vec1_underThreshold[-1]

    # Extract data following the peak only (data should be falling at the start
    # of the new array)
    vec2 = vec[peakInd:]
    # Find indices less than half the max
    vec2_underThreshold = np.where(vec2 < halfMax)[0]
    if len(vec2_underThreshold) == 0:
        # There are no indices less than half the max, so we have to estimate
        # the upper edge
        estimateUpperEdge = True
    else:
        # There are indices less than half the max, take the first one under
        # halfMax as the upper edge
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


def _find_localmax_method_1(tfr, freqs, times, event_band,
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
            hit/miss,         maxima frequency,
            lowerbound frequency,     upperbound frequency,
            frequency span,         maxima timing,     event onset timing,
            event offset timing,     event duration, maxima power,
            maxima/median power
    '''

    n_epochs = tfr.shape[0]

    all_epochs_events = list()

    # Retrieve all local maxima in tfr using python equivalent of imregionalmax
    for trial_idx in range(n_epochs):

        epoch_events = list()

        # Get tfr data for this trial [frequency x time]
        thistfr = tfr[trial_idx, :, :]

        # Find local maxima in the tfr data
        data = thistfr
        # Find maximum amoung adjacent pixels (3x3 footprint) for each pixel
        data_max = filters.maximum_filter(data, size=(3, 3))
        maxima = (data == data_max)
        data_min = filters.minimum_filter(data, size=(3, 3))
        # Rule out pixels with footprints that have flatlined
        maxima[data_max == data_min] = False
        labeled, num_objects = ndimage.label(maxima)
        yx = ndimage.center_of_mass(data, labels=labeled,
                                    index=range(1, num_objects + 1))

        peakF = list()
        peakT = list()
        peakPower = list()
        for f_idx, t_idx in yx:
            f_idx = int(round(f_idx))
            t_idx = int(round(t_idx))
            event_freq = freqs[f_idx]
            if (event_freq >= event_band[0] and event_freq <= event_band[1] and
               thistfr[f_idx, t_idx] > eventThresholdByFrequency[f_idx]):
                peakF.append(f_idx)
                peakT.append(t_idx)
                peakPower.append(thistfr[f_idx, t_idx])

        numPeaks = len(peakF)

        # Find local maxima lowerbound, upperbound, and full width at half max
        #    for both frequency and time
        for lmi in range(numPeaks):
            thisPeakF = peakF[lmi]
            thisPeakT = peakT[lmi]
            thisPeakPower = peakPower[lmi]

            # Indices of tfr frequencies < half max power at the time of a
            # given local peak
            tfrFrequencies = thistfr[:, thisPeakT]
            lowerInd, upperInd, FWHM = _fwhm_lower_upper_bound1(tfrFrequencies,
                                                                thisPeakF,
                                                                thisPeakPower)
            lowerEdgeFreq = freqs[lowerInd]
            upperEdgeFreq = freqs[upperInd]
            FWHMFreq = FWHM * (freqs[1] - freqs[0])

            # Indices of tfr times < half max power at the frequency of a given
            # local peak
            tfrTimes = thistfr[thisPeakF, :]
            lowerInd, upperInd, FWHM = _fwhm_lower_upper_bound1(tfrTimes,
                                                                thisPeakT,
                                                                thisPeakPower)
            lowerEdgeTime = times[lowerInd]
            upperEdgeTime = times[upperInd]
            FWHMTime = FWHM / Fs

            # Put peak characteristics to a dictionary
            peakParameters = {
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
            epoch_events.append(peakParameters)

        all_epochs_events.append(epoch_events)

    return all_epochs_events


def plot_events(tfr, times, freqs, event_band, spec_events=None,
                timeseries=None, ax=None, vlim=None, ylim_ts=None, label=None):
    '''Plot a single TFR spectrogram overlayed with spectral events.

    Parameters
    ----------
    tfr : array, shape (n_freqs, n_times)
        A time-frequency response (TFR).
    times : array-like, shape (n_times,)
        Time domain of the TFR in seconds.
    freqs : array-like, shape (n_freqs,)
        Frequency domain of the TFR in Hertz.
    event_band : list
        Lower and upper bounds (inclusive, Hz) of the frequency band-of-
        interest in which spectral events were detected.
    spec_events : list of dict
        A single-level list where each element comprises an event that is
        characterized by dictionary items including it's location in time,
        location in frequency, duration, and frequency span, and more. If None,
        no events are plotted.
    timeseries : array, shape (n_times,) | None
        The timeseries associated with the TFR.
    ax : None | Matplotlib.axes.Axes instance
        The Axes instance with wich to plot the spectrogram on. If None, a new
        Axes object is created.
    vlim : list, shape (2,) | None
        If not None, sets the colorbar lower and upper bounds for the
        spectrogram across all axes of the returned figure. If None, each
        spectrogram (including example epochs) sets its vlim independently.
    ylim_ts : list | None
        The lower and upper limits of the timeseries (if applicable) overlayed
        on the spectrogram. If None, these values are automatically assigned.
    label : str | None
        A label to overlay in the upper left corner of this plot.

    Return
    ------
    fig : Matplotlib.figure.Figure instance
        The Figure instance containing the spectrogram.
    '''

    if tfr.shape != (len(freqs), len(times)):
        raise ValueError(f'tfr must be an array of shape (n_freqs, n_times), '
                         f'got tfr: {tfr.shape}, freqs: ({len(freqs)},), '
                         f'times: ({len(times)},)')

    if vlim is None:
        vlim = [None, None]

    # convert to numpy array if not already
    freqs = np.array(freqs)

    # frequencies within the band of interest
    # band_mask = np.logical_and(freqs >= event_band[0],
    #                            freqs <= event_band[1])

    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        fig = ax.get_figure()

    # plot tfr
    im = ax.pcolormesh(times, freqs, tfr, cmap='jet', vmin=vlim[0],
                       vmax=vlim[1], shading='nearest')
    fig.colorbar(im, ax=ax)
    ax.axhline(y=event_band[0], c='w', linewidth=2., linestyle=':', alpha=.7)
    ax.axhline(y=event_band[1], c='w', linewidth=2., linestyle=':', alpha=.7)
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
        ax.annotate(label, xy=(0.01, 0.98), xycoords='axes fraction',
                    va='top', ha='left', color='w', size=10, fontweight='bold')

    return fig


def plot_avg_spectrogram(tfr, times, freqs, event_band, spec_events=None,
                         timeseries=None, example_epochs=None, vlim=None,
                         show_events=False):
    '''Plot the average spectrogram over all TFR epochs/trials.

    Parameters
    ----------
    tfr : array, shape (n_epochs, n_freqs, n_times)
        A stack of time-frequency response (TFR) epochs/trials.
    times : array-like, shape (n_times,)
        Time domain of each TFR in seconds.
    freqs : array-like, shape (n_freqs,)
        Frequency domain of each TFR in Hertz.
    event_band : list
        Lower and upper bounds (inclusive, Hz) of the frequency band-of-
        interest in which spectral events were detected.
    spec_events : list of list of dict
        A nested list with n_epochs elements in the outer list and n_events in
        the inner list. Each inner element comprises an event that is
        characterized by dictionary items including it's location in time,
        location in frequency, duration, and frequency span, and more. If None,
        no events are plotted.
    timeseries : array, shape (n_epochs, n_times) | None
        The stack of timeseries associated with the TFR epochs/trials.
    example_epochs : list | None
        The epoch/trial indices that will be used to plot example spectrogram
        trials alonside the average spectrogram.
    vlim : list, shape (2,) | None
        If not None, sets the colorbar lower and upper bounds for the
        spectrogram across all axes of the returned figure. If None, each
        spectrogram (including example epochs) sets its vlim independently.
    show_events : boolean
        Whether or not to overlay spectral event locations on the average
        spectrogram.

    Return
    ------
    fig : Matplotlib.figure.Figure instance
        The Figure instance containing the average spectrogram.
    '''

    if len(tfr.shape) != 3:
        raise ValueError(f'tfr should be a 3D array of shape (n_epochs, '
                         f'n_freqs, n_times), got {tfr.shape}.')

    if example_epochs is not None:
        trial_idx_set = set(range(tfr.shape[0]))
        trial_idx_subset = set(example_epochs)
        if trial_idx_subset.intersection(trial_idx_set) != trial_idx_subset:
            raise ValueError('One or more of the specified example trial '
                             'indices does not exist in the provided tfr '
                             'array.')
    else:
        # set to empty list
        example_epochs = list()

    if show_events:
        spec_events_agg = sum(spec_events, [])
    else:
        spec_events_agg = None

    fig, axs = plt.subplots(nrows=len(example_epochs) + 1, ncols=1,
                            sharex=True)

    # if example_epochs is None, ensure that axs is still subscriptable
    axs = np.atleast_1d(axs)

    # plot trial-average tfr
    tfr_avg = np.mean(tfr, axis=0).squeeze()
    plot_events(tfr=tfr_avg, times=times, freqs=freqs,
                event_band=event_band, spec_events=spec_events_agg,
                ax=axs[0], vlim=vlim, label='epoch avg.')

    # plot tfr + events for example trials
    if timeseries is not None and example_epochs is not None:
        max_ts_amplitude = np.max(timeseries[example_epochs])
        min_ts_amplitude = np.min(timeseries[example_epochs])
        ylim_ts = [max_ts_amplitude, min_ts_amplitude]

        for plot_idx, trial_idx in enumerate(example_epochs):
            # get spectral events for the current trial
            if spec_events is not None:
                trial_events = spec_events[trial_idx]
            else:
                trial_events = None

            # plot trial tfr
            tfr_trial = tfr[trial_idx, :, :].squeeze()
            timeseries_trial = timeseries[trial_idx, :]

            plot_events(tfr=tfr_trial, times=times, freqs=freqs,
                        event_band=event_band, spec_events=trial_events,
                        timeseries=timeseries_trial, ax=axs[plot_idx + 1],
                        vlim=vlim, ylim_ts=ylim_ts, label=f'epoch {trial_idx}')

    axs[-1].set_xlabel('time (s)')
    axs[0].set_ylabel('freq. (Hz)')
    fig.tight_layout()

    return fig
