% general filter and time-frequency considerations

% segmentsFilt@100: lowpasspar = 0.01; highpasspar = 100; firorder = 65;

display = 1;   % 1 if want to display
show = [];
samplefreq = 1450;
waveFrq = [4,350];       % Transform frequency range
rowsPerOct = 32;
padmode = 'zpd';
wavelet = 'morl';          % Mother wavelet; must be either 'morl' | 'mexh'
t= (0:size(segment,2)-1)/samplefreq*1000 - 210;

filtersignal = 0; % if to filter set 1
lowpasspar = 0.01; % was = 0.3
highpasspar = 100; % was 40
firorder = 65;  % was = 60

figure
plot(t,segment')
if filtersignal
    segment = eegfilt(segment,samplefreq,lowpasspar,highpasspar,0,firorder);
end
figure
plot(t,segment')

if display
    show = 'image';
end

 % 'wavePlot',false...
chan = 61;
[wcf, pfreq, scales]= wavtrans(segment,chan,t,samplefreq,rowsPerOct,waveFrq,padmode,wavelet,show);