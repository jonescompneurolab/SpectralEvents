[![tests](https://github.com/jonescompneurolab/SpectralEvents/actions/workflows/tests.yml/badge.svg?branch=master)](https://github.com/jonescompneurolab/SpectralEvents/actions/?query=branch:master+event:push)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jonescompneurolab/SpectralEvents/HEAD?labpath=tutorial.ipynb)

# SpectralEvents Toolbox
This toolbox is composed of a series of functions that find and analyze transient spectral activity (spectral events) in a time-series dataset, allowing for spectral event feature comparison between trial/epoch outcomes/conditions. Spectral events are defined as local maxima above a power threshold of a specified band in the non-averaged time-frequency response (TFR).

This repository contains two versions written in Python and Matlab, respectively. While both versions render similar
results, we encourage new users to use the Python version since it contains the most recent features and is more actively maintained.

## Usage (Python)
After cloning this repository, you can import the `spectralevents` module by
first setting the path to the `SpectralEvents` directory.

```python
import sys
sys.path.append('/PATH/TO/FOLDER/SpectralEvents')
import spectralevents
```

You'll need to load your timeseries datadata as an array of epochs-by-time
samples. With your data loaded, the general workflow for event detection
follows this progression: **TIMESERIES** -> **TFR** -> **SPECTRAL EVENTS**.
This can be accomplished using the `spectralevents.tfr` and 
`spectralevents.find_events` functions, which return an array of epoch TFRs and nested
list of spectral events, respectively. Keep in mind that there are function
arguments that you must set that influence how Spectral Event detection will be
conducted. In particular, these include
- `freqs`, the frequency values over which you will calculate your **t**ime-**f**requency **r**esponse (TFR)
- `times`, the time values at which your signal was sampled relative to each epoch or trial
- `event_band`, the bounds of the frequency band in which you will look for Spectral Events
- `thresh_FOM`, the factor-of-the-median threshold that will be used to find suprathreshold spectral power in each frequency bin across the spectrogram(s)

Note that the Python version currently only supports event Find Method 1,
which is the same method introduced in Shin et al. eLife 2017 that sets no 
explicit restrictions on event overlap in the time-frequency domain.
Spectral Event Analysis results can be visualized by passing the results of
your TFR and Spectral Events tranformations into
`spectralevents.plot_avg_spectrogram` and `spectralevents.plot_events`. For a
more detailed introduction to Spectral Event Analysis, see the Jupyter Notebook
tutorial, `tutorial.ipynb`.

### Dependencies
- Python 3
- numpy
- scipy
- matplotlib

## Usage (Matlab)
The Matlab version follows a similar workflow as the Python version, albeit
with an additional wrapper function that can be used to streamline your entire
analysis workflow, `spectralevents.m`. Depending on your data and research
goals, `spectralevents.m` might be helpful as it allows you to analyze data
from multiple subjects/sessions, each containing multiple epochs or trials, and
return results all with a single line of code (see `example.m`). Alternatively,
some use-cases might be better served by calling the constituent sub-functions
`spectralevents_ts2tfr.m`, `spectralevents_find.m`, and `spectralevents_vis.m`
individually as they follow a similar logical structure as the Python
version described above. These functions are described in more detail as
follows.

### spectralevents
```
[specEv_struct, TFRs, X] = spectralevents(eventBand, fVec, Fs, findMethod, vis, X, classLabels)
```
or
```
[specEv_struct, TFRs, X] = spectralevents(eventBand, fVec, Fs, findMethod, vis, X{1}, classLabels{1}, X{2}, classLabels{2}, ...)
```
Imports time-series dataset and finds the TFR for each trial using `spectralevents_ts2tfr`, calls the event-finding function `spectralevents_find` for each subject/session within the dataset, and runs `spectralevents_vis` in order to capture and view spectral event features.

Returns a structure array of spectral event features (`specEv_struct`), cell array containing the time-frequency responses (`TFRs`), and cell array of all time-series trials (`X`) for each subject/session within the dataset comparing various experimental conditions or outcome states corresponding to each trial. By default, this function sets the factors of median threshold at 6.

*IMPORTANT*: the `findMethod` specified in `spectralevents_find` will bias the results. It’s important to understand the different methods and to state clearly in any presentation or publication which method was used. In the foundational paper for this toolbox, Shin et al. eLife 2017, `findMethod=1` was applied.

Inputs:
* `eventBand` - range of frequencies ([Fmin_event, Fmax_event]; Hz) over which above-threshold spectral power events are classified.
* `fVec` - frequency vector (Hz) over which the time-frequency response (TFR) is calculated. Note that this set must fall within the range of resolvable/alias-free frequency values (i.e. Fmin>=1/(trial duration), Fmax<=(Nyquist freq)).
* `Fs` - sampling frequency (Hz).
* `findMethod` - integer value specifying which event-finding method to use (`1`, `2`, or `3`). Note that the method specifies how much overlap exists between events. Use `1` to replicate the method used in Shin et al. eLife 2017.
* `vis` - logical value that determines whether to run basic feature analysis and output standard figures.
* `X{a}` - m-by-n matrix (of the a-th subject/session cell in cell array X) representing the time-series trials of the given subject. m is the number of timepoints and n is the number of trials. Note that m timepoints must be uniform across all trials and subjects.
* `classLabels{a}` - numeric or logical trial classification labels (of the a-th subject/session cell in cell array classLabels); associates each trial of the given subject/session to an experimental condition/outcome/state (e.g., hit or miss, detect or non-detect, attend-to or attend away). If classLabels{a} is entered as a single value (e.g., 0 or 1), all trials in the a-th subject/session are associated with that label. Alternatively, classLabels{a} can be entered as a vector of n elements, each corresponding to a trial within the a-th subject/session.

Outputs:
* `specEv_struct` - array of event feature structures, each corresponding with a subject/session, respectively.
* `TFRs` - cell array with each cell containing the time-frequency response (freq-by-time-by-trial) for a given subject/session.
* `X` - cell array with each cell containing the time-series trials for a given subject/session.

### spectralevents_find
```
specEv_struct = spectralevents_find(findMethod,eventBand,thrFOM,tVec,fVec,TFR,classLabels)
```
Algorithm for finding and calculating spectral events on a trial-by-trial basis of of a single subject/session. Uses one of three methods before further analyzing and organizing event features:

1. (Primary event detection method in Shin et al. eLife 2017): Find spectral events by first retrieving all local maxima in un-normalized TFR using imregionalmax, then selecting suprathreshold peaks within the frequency band of interest. This method allows for multiple, overlapping events to occur in a given suprathreshold region and does not guarantee the presence of within-band, suprathreshold activity in any given trial will render an event.
2. Find spectral events by first thresholding entire normalize TFR (over all frequencies), then finding local maxima. Discard those of lesser magnitude in each suprathreshold region, respectively, s.t. only the greatest local maximum in each region survives (when more than one local maxima in a region have the same greatest value, their respective event timing, frequency location, and boundaries at full-width half-max are calculated separately and averaged). This method does not allow for overlapping events to occur in a given suprathreshold region and does not guarantee the presence of within-band, suprathreshold activity in any given trial will render an event.
3. Find spectral events by first thresholding normalized TFR in frequency band of interest, then finding local maxima. Discard those of lesser magnitude in each suprathreshold region, respectively, s.t. only the greatest local maximum in each region survives (when more than one local maxima in a region have the same greatest value, their respective event timing, freq. location, and boundaries at full-width half-max are calculated separately and averaged). This method does not allow for overlapping events to occur in a given suprathreshold region and ensures the presence of within-band, suprathreshold activity in any given trial will render an event.

Inputs:
* `findMethod` - integer value specifying which event-finding method to use (`1`, `2`, or `3`). Note that the method specifies how much overlap exists between events. Use `1` to replicate the method used in Shin et al. eLife 2017.
* `eventBand` - range of frequencies ([Fmin_event Fmax_event]; Hz) over which above-threshold spectral power events are classified.
* `thrFOM` - factors of median threshold; positive real number used to threshold local maxima and classify events (see Shin et al. eLife 2017 for discussion concerning this value).
* `tVec` - time vector (s) over which the time-frequency response (TFR) is calculated.
* `fVec` - frequency vector (Hz) over which the time-frequency response (TFR) is calculated.
* `TFR` - time-frequency response (TFR) (frequency-by-time-trial) for a single subject/session.
* `classLabels` - numeric or logical 1-row array of trial classification labels; associates each trial of the given subject/session to an experimental condition/outcome/state (e.g., hit or miss, detect or non-detect, attend-to or attend away).

Outputs:
* `specEv_struct` - event feature structure with three main sub-structures: TrialSummary (trial-level features), Events (individual event characteristics), and IEI (inter-event intervals from all trials and those associated with only a given class label).

### spectralevents_ts2tfr
```
[TFR, tVec, fVec] = spectralevents_ts2tfr(S, fVec, Fs, width)
```
Calculates the TFR (in spectral power) of a time-series waveform by convolving in the time-domain with a Morlet wavelet. **Note that this version does not currently control for convolution edge effects where the Morlet wavelet is not completely overlapping the time-series. The purpose of this TFR calculation is to give an approximation of transient activity, though a more thorough analysis should crop edge effects off of the TFR.**

Inputs:
* `S` - column vector of the time-series signal.
* `fVec` - frequency vector (Hz) over which the time-frequency response (TFR) is calculated.
* `Fs` - sampling frequency (Hz).
* `width` - number of cycles in wavelet.

Outputs:
* `TFR` - frequency-by-time time-frequency response (TFR).
* `tVec` - time vector (s) over which the time-frequency response (TFR) is calculated.
* `fVec` - frequency vector (Hz) over which the time-frequency response (TFR) is calculated.

### spectralevents_vis
```
spectralevents_vis(specEv_struct, timeseries, TFRs, tVec, fVec)
```
Conducts basic analysis for the purpose of visualizing dataset spectral event features and generates spectrogram and probability histogram plots.

Inputs:
* `specEv_struct` - spectralevents structure array.
* `timeseries` - cell array containing time-series trials by subject/session.
* `TFRs` - cell array containing time-frequency responses by subject/session.
* `tVec` - time vector (s) over which the time-frequency responses are shown.
* `fVec` - frequency vector (Hz) over which the time-frequency responses are shown.

### Dependencies
- MATLAB R2019a
- MATLAB Image Processing Toolbox

## Contributors:
* Hyeyoung Shin, Department of Neuroscience, Brown University, shinehyeyoung@gmail.com
* Ryan Thorpe, Department of Neuroscience, Brown University, ryvthorpe@gmail.com
* Tim Bardouille, Department of Physics and Atmospheric Science, Dalhousie University, tim.bardouille@dal.ca
* Darcy Diesburg, Department of Neuroscience, Brown University, darcy.diesburg@gmail.com
