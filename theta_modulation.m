% quantify theta modulation around stimulation events using wavelet
% transform

% load data

load('cellg110319.mat')
load('cellg110319stim.mat')
volt = Ch3.values; % in mV
t = Ch3.times;     % in s
curr = Ch6.values; % in nA
samplefreq = 1/Ch3.interval;  % in Hz

% parameters

stimulus = 'call';
times = NonSibCall;
display = 1;   % 1 if want to display
waveFrq = [6,10];       % Transform frequency range
rowsPerOct = 32;
padmode = 'zpd';
wavelet = 'mexh';          % Mother wavelet; must be either 'morl' | 'mexh'

if strcmp(stimulus,'odor') 
    bef = 2; % in s, baseline before stimulus presentation
    stim = 5; % in s, stimulation duration
    aft = 5; % in s, rebound activity after stimulation
elseif strcmp(stimulus,'call') 
    bef = 2;
    stim = 1;
    aft = 1;
end

bef_size = floor(bef*samplefreq);
stim_size = floor(stim*samplefreq);
aft_size = floor(aft*samplefreq);

if display
    show = 'image';
else
    show = false;
end

amps = zeros(size(times,1),3);

for i=1:size(times,1)
    index = times(i)*samplefreq;
    signal = volt((index-bef_size):(index+stim_size+aft_size));
    t = (times(i) - bef):Ch3.interval:(times(i) + stim + aft);
    [wcf, pfreq, scales] = wavtrans(signal,t,samplefreq,rowsPerOct,waveFrq,padmode,wavelet,show);
    amps(i,1) = sum(wcf(:,1:bef_size),'all');   % amplitude before stimulation
    amps(i,2) = sum(wcf(:,(bef_size+1):(bef_size+stim_size)),'all');   % amplitude during stimulation
    amps(i,3) = sum(wcf(:,(bef_size+stim_size+1):end),'all');   % amplitude after stimulation
    w = waitforbuttonpress;
end

figure
plot(amps)
xlabel('Event #')
ylabel('Intensity of theta rhythm')
legend('Baseline','Stimulation','Rebound')
title('Modulation of theta rhythm by stimulation (Non-Sibling Call)')