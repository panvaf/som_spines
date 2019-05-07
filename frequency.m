tic
% load data

load('cellg110319.mat')
volt = Ch3.values; % in mV
t = Ch3.times;     % in s
curr = Ch6.values; % in nA
samplefreq = 1/Ch3.interval;  % in Hz

% remove segments where stimulation occurs

seg = [1 100; 550 750;];  % segments to be removed, in seconds
seg = seg*samplefreq;
for i = 1:size(seg,2)
    volt(seg(i,1):seg(i,2)) = [];
    t(seg(i,1):seg(i,2)) = [];
    curr(seg(i,1):seg(i,2)) = [];
end

% general filter and time-frequency considerations

display = 1;   % 1 if want to display
show = [];
waveFrq = [40,350];       % Transform frequency range (from Logothetis paper)
rowsPerOct = 32;
padmode = 'zpd';
wavelet = 'mexh';          % Mother wavelet; must be either 'morl' | 'mexh'

filtersignal = 0; % if to filter set 1
lowpasspar = 0.01; % was = 0.3
highpasspar = 100; % was 40
firorder = 65;  % was = 60

figure
plot(t,volt)
if filtersignal
    %y = eegfilt(y,samplefreq,lowpasspar,highpasspar,0,firorder);
end
figure
plot(t,volt)

if display
    show = 'image';
end

% 'wavePlot',false...
% Compute wavelet transform

signal = volt(50000:100000); time = Ch3.times(50000:100000);
[wcf, pfreq, scales]= wavtrans(signal,time,samplefreq,rowsPerOct,waveFrq,padmode,wavelet,show);
toc