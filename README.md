# SpectralEvents Toolbox

This toolbox is composed of a series of functions that find and analyze spectral events (local maxima above a power threshold within a specified band of the non-averaged time-frequency response) on a trial-by-trial basis allowing for spectral event feature comparison between trial outcomes/conditions.

## Functions
### spectralevents
```
[specEv_struct,TFRs,X] = spectralevents(eventBand,fVec,Fs,findMethod,vis,X,classLabels)
```
or
```
[specEv_struct,TFRs,X] = spectralevents(eventBand,fVec,Fs,findMethod,vis,X{1},classLabels{1},X{2},classLabels{2},...)
```
Imports time-series dataset and finds the TFR for each trial using `traces2TFR` from 4DToolbox, calls event-finding method, `spectralevents_find_1` or `spectralevent_find_2`, for each subject/session within the dataset, and runs `spectralevents_vis` in order to capture and view spectral event features.

Returns a structure array of spectral event features (specEv_struct), cell array containing the time-frequency responses (TFRs), and cell array of all time-series trials (X) for each subject/session within the dataset comparing various experimental conditions or outcome states corresponding to each trial. By default, this function sets the factors of median threshold at 6.

Inputs:
* `eventBand` - range of frequencies ([Fmin_event,Fmax_event]; Hz) over which above-threshold spectral power events are classified.
* `fVec` - frequency vector (Hz) over which the time-frequency response (TFR) is calculated. Note that this set must fall within the range of resolvable/alias-free frequency values (i.e. Fmin>=1/(trial duration), Fmax<=(Nyquist freq)).
* `Fs` - sampling frequency (Hz).
* `findMethod` - integer value specifying which event-finding method
function (i.e., `spectralevents_find_1` or `spectralevents_find_2`) to run. Note that the method specifies how much overlap exists between events.
* `vis` - logical value that determines whether to run basic feature analysis and output standard figures.
* `X{a}` - m-by-n matrix (of the a-th subject/session cell in cell array X) representing the time-series trials of the given subject. m is the number of timepoints and n is the number of trials. Note that m timepoints must be uniform across all trials and subjects.
* `classLabels{a}` - numeric or logical trial classification labels (of the a-th subject/session cell in cell array classLabels); associates each trial of the given subject/session to an experimental condition/outcome/state (e.g., hit or miss, detect or non-detect, attend-to or attend away). If classLabels{a} is entered as a single value (e.g., 0 or 1), all trials in the a-th subject/session are associated with that label. Alternatively, classLabels{a} can be entered as a vector of n elements, each corresponding to a trial within the a-th subject/session.

Outputs:
* `specEv_struct` - array of event feature structures, each corresponding with a subject/session, respectively.
* `TFRs` - cell array with each cell containing the time-frequency response (freq-by-time-by-trial) for a given subject/session.
* `X` - cell array with each cell containing the time-series trials for a given subject/session.

### spectralevents_find_1
```
specEv_struct = spectralevents_find_1(eventBand,thrFOM,tVec,fVec,TFR,classLabels)
```
Algorithm for finding and calculating spectral events on a trial-by-trial basis of of a single subject/session. As the primary event detection method in Shin et al. eLife 2017, events are found as follows:
1. Retrieve all local maxima in TFR using imregionalmax
2. Pick out maxima above threshold and within the frequency band of interest
3. Identify and organize event features

Inputs:
* `eventBand` - range of frequencies ([Fmin_event Fmax_event]; Hz) over which above-threshold spectral power events are classified.
* `thrFOM` - factors of median threshold; positive real number used to threshold local maxima and classify events (see Shin et al. eLife 2017 for discussion concerning this value).
* `tVec` - time vector (s) over which the time-frequency response (TFR) is calcuated.
* `fVec` - frequency vector (Hz) over which the time-frequency response (TFR) is calcuated.
* `TFR` - time-frequency response (TFR) (frequency-by-time-trial) for a single subject/session.
* `classLabels` - numeric or logical 1-row array of trial classification labels; associates each trial of the given subject/session to an experimental condition/outcome/state (e.g., hit or miss, detect or non-detect, attend-to or attend away).

Outputs:
* `specEv_struct` - event feature structure with three main sub-structures: TrialSummary (trial-level features), Events (individual event characteristics), and IEI (inter-event intervals from all trials and those associated with only a given class label).

### spectralevents_find_2
```
specEv_struct = spectralevents_find_2(eventBand,thrFOM,tVec,fVec,TFR,classLabels)
```
Algorithm for finding and calculating spectral events on a trial-by-trial basis of of a single subject/session. Uses the following method:
1. Retrieve all suprathreshold local maxima in normalized TFR using imregionalmax
2. Discard those of lesser magnitude in each suprathreshold region, respectively, s.t. only the greatest local maximum in each region survives (when more than one local maxima in a region have the same greatest value, their respective event timing, freq. location, and boundaries at full-width half-max are calculated separately and averaged)
3. Of the remaining maxima, pick out maxima within the frequency band of interest
4. Identify and organize event features

Inputs:
* `eventBand` - range of frequencies ([Fmin_event Fmax_event]; Hz) over which above-threshold spectral power events are classified.
* `thrFOM` - factors of median threshold; positive real number used to threshold local maxima and classify events (see Shin et al. eLife 2017 for discussion concerning this value).
* `tVec` - time vector (s) over which the time-frequency response (TFR) is calcuated.
* `fVec` - frequency vector (Hz) over which the time-frequency response (TFR) is calcuated.
* `TFR` - time-frequency response (TFR) (frequency-by-time-trial) for a single subject/session.
* `classLabels` - numeric or logical 1-row array of trial classification labels; associates each trial of the given subject/session to an experimental condition/outcome/state (e.g., hit or miss, detect or non-detect, attend-to or attend away).

Outputs:
* `specEv_struct` - event feature structure with three main sub-structures: TrialSummary (trial-level features), Events (individual event characteristics), and IEI (inter-event intervals from all trials and those associated with only a given class label).

### spectralevents_vis (script)
```
spectralevents_vis(specEv_struct,timeseries,TFRs, tVec,fVec)
```
Conducts basic analysis for the purpose of visualizing dataset spectral event features and generates spectrogram and probability histogram plots.

Inputs:
* `specEv_struct` - spectralevents structure array.
* `timeseries` - cell array containing time-series trials by subject/session.
* `TFRs` - cell array containing time-frequency responses by subject/session.
* `tVec` - time vector (s) over which the time-frequency responses are shown.
* `fVec` - frequency vector (Hz) over which the time-frequency responses are shown.

## Example
See test script `test.m`.

## Dependencies
* 4DToolbox by Ole Jensen

## Contributors:
* Hyeyoung Shin, Department of Neuroscience, Brown University, shinehyeyoung@gmail.com
* Ryan Thorpe, School of Engineering, Brown University, ryvthorpe@gmail.com

Written for MATLAB R2017b.
