function [abs_trans, pseudofreq, scales] = wavtrans(record, t, fsample, rowsPerOct, freqSpan, padmode, wavelet, plottype)
% Calculate and show the wavelet transform for a single eeg channel
%
% record: 1x[time] Vector of recording
% fsample: sampling frequency
% voicesPerOct: defines the transform's resolution in scales
% freqSpan: frequency span to compute transform
% padmode: see dwtmode for other choices
% mwave: 'morl' (morlet) or 'mexh' (mexican hat)

% Get adjusted scales
scales= helperCWTTimeFreqVector(freqSpan(1),freqSpan(2),centfrq(wavelet),1/fsample,rowsPerOct);
pseudofreq= mean([freqSpan(1)*scales(end),freqSpan(2)*scales(1)])./scales;

% Transform
trans = cwtft({record,1/fsample},'wavelet',wavelet,'scales',scales,'padmode',padmode);
abs_trans = abs(trans.cfs);

% Plot
if plottype == false
else
    figure
    wscalogram(plottype,flip(abs_trans,1),'scales',flip(round(pseudofreq),2),'xdata',t,'ydata',record);
    xlabel('Time (s)'); ylabel('Freq (Hz)');
    title('Pseudofrequency-time transform');

end

end