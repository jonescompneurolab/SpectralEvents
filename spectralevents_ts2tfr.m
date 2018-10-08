function [TFR,tVec,fVec] = spectralevents_ts2tfr(S,fVec,Fs,width)
% function [TFR,tVec,fVec] = spectralevents_ts2tfr(S,fVec,Fs,width);
%
% Calculates the TFR (in spectral power) of a time-series waveform by 
% convolving in the time-domain with a Morlet wavelet.                            
%
% Input
% -----
% S    : signals = time x Trials      
% fVec    : frequencies over which to calculate TF energy        
% Fs   : sampling frequency
% width: number of cycles in wavelet (> 5 advisable)  
%
% Output
% ------
% t    : time
% f    : frequency
% B    : phase-locking factor = frequency x time
%     
% Adapted from Ole Jensen's traces2TFR in the 4Dtools toolbox.
%
% See also SPECTRALEVENTS, SPECTRALEVENTS_FIND, SPECTRALEVENTS_VIS.

%   -----------------------------------------------------------------------
%   SpectralEvents::spectralevents_ts2tfr
%   Copyright (C) 2018  Ryan Thorpe
%
%   This file is part of the SpectralEvents toolbox.
% 
%   SpectralEvents is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   SpectralEvents is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.
%   -----------------------------------------------------------------------

S = S';

tVec = (1:size(S,2))/Fs;  

B = zeros(length(fVec),size(S,2)); 

for i=1:size(S,1)          
    for j=1:length(fVec)
        B(j,:) = energyvec(fVec(j),detrend(S(i,:)),Fs,width) + B(j,:);
    end
end
TFR = B/size(S,1);     


function y = energyvec(f,s,Fs,width)
% Return a vector containing the energy as a
% function of time for frequency f. The energy
% is calculated using Morlet's wavelets. 
% s : signal
% Fs: sampling frequency
% width : width of Morlet wavelet (>= 5 suggested).

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);

y = conv(s,m);
y = 2*(dt*abs(y)).^2;
y = y(ceil(length(m)/2):length(y)-floor(length(m)/2));
end

function y = morlet(f,t,width)
% Morlet's wavelet for frequency f and time t. 
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet. 
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)

sf = f/width;
st = 1/(2*pi*sf);
A = 1/(st*sqrt(2*pi));
y = A*exp(-t.^2/(2*st^2)).*exp(1i*2*pi*f.*t);
end

end
