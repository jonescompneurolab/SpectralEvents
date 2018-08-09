# SpectralEvents Toolbox

Find and analyze spectral events (local maxima above a power threshold) of a specified band in the non-averaged time-frequency response (TFR) in a dataset.

```
[specEv_struct,TFRs,X] = SPECTRALEVENTS(eventBand,analyze,fVec,Fs,X,classLabels)
```
or
```
[specEv_struct,TFRs,X] = SPECTRALEVENTS(eventBand,analyze,fVec,Fs,X{1},classLabels{1},X{2},classLabels{2},...)
```
returns a structure array of spectral event features (specEv_struct), cell array containing the time-frequency responses (TFRs), and cell array of all time-series trials (X) for each subject/session within the dataset comparing various experimental conditions or outcome states corresponding to each trial. By default, this function sets the factors of median threshold at 6.

### Dependencies

  * 4DToolbox by Ole Jensen, Helsinki University of Technology

### Contributors:
  * Hyeyoung Shin, Department of Neuroscience, Brown University, shinehyeyoung@gmail.com
  * Ryan Thorpe, School of Engineering, Brown University, ryvthorpe@gmail.com

Written for MATLAB R2017b.
