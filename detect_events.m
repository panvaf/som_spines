% receive and manipulate wavelet transform of physiological recordings in
% order to detect various types of events

tic
% load data

load('cella270319.mat')
volt = Ch3.values; % in mV
curr = Ch6.values; % in nA
t = Ch3.times;     % in s
sampleint = Ch3.interval;
samplefreq = 1/sampleint;  % in Hz

% remove segments that are not wanted (stimulation etc)
% not needed since the wavelet algorithm is ultra fast
% better remove events in the end

% segments to be removed, in seconds

seg = [1040 1090; 1589 1610;];          % for cella270319
%seg = [110 150; 1630 1660;];          % for cella260319
%seg = [368 490;];                     % for cellb250319
%seg = [];                             % for cella250319
%seg = [0 60; 105 140; 180 210;];      % for cellb220319
%seg = [10 90; 160 190;];              % for cella220319
%seg = [0 100;];                       % for cellb080319
%seg = [];                             % for cellg040219
%seg = [30 45;];                       % for cellb110319
%seg = [0 40; 390 450; 1030 1180;];    % for cellb070319
%seg = [0 50; 1600 1750;];             % for celle070319
%seg = [290 360; 450 520; 1280 1320;]; % for cellf110319
%seg = [0 30; 1620 1730;];             % for cellf070319
%seg = [1280 1380;];                   % for cellc020419
%seg = [650 750;];                     % for celll230119
%seg = [0 100; 550 750;];              % for cellg110319
%{
seg = seg*samplefreq;
for i = 1:size(seg,2)
    volt(seg(i,1):seg(i,2)) = [];
    t(seg(i,1):seg(i,2)) = [];
    curr(seg(i,1):seg(i,2)) = [];
end
%}

% parameters

batch = 1e6;  % batch size to break computation in parts
win = .1; % in s, maximum expected size of event
win_size = floor(win*samplefreq);
threshold = 300; % for detection of events in general
threshold_sp = 400; % events above this threshold are classified as spikes
isclose = win_size/10; % collate events that are close enough
len = win_size/2;  % lenght of window of integration for detection

display = 0;   % 1 if want to display
show = [];
waveFrq = [50,100];       % Transform frequency range
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
else
    show = false;
end

% Compute wavelet transform in batches

wavs = [];
times = [];
amps = [];

for i=1:ceil(size(volt,1)/batch)
    disp(i)
    fin = min(i*batch,size(volt,1));
    signal = volt(((i-1)*batch+1):fin); t = Ch3.times(((i-1)*batch+1):fin);
    [wcf, pfreq, scales] = wavtrans(signal,t,samplefreq,rowsPerOct,waveFrq,padmode,wavelet,show);
    
    % extract events by the energy in a moving window
    amp = movmean(sum(wcf),len);
    trans = diff(amp>threshold);
    
    % merge events that are too close
    start = find(trans==1); fin = find(trans==-1);
    
    if size(start,2) > size(fin,2)
        fin(size(start,2)) = size(trans,2);
    elseif size(start,2) < size(fin,2)
        start = [1 start];
    end
    
    for j=2:size(fin,2)
        if start(j)-fin(j-1)<isclose
            start(j) = 0;
            fin(j-1) = 0;
        end
    end

    start(start==0) = [];
    fin(fin==0) = [];
    
    % center around peaks of events
    event = false(1,size(trans,2));
    for j=1:size(start,2)
        [argvalue, argmax] = max(amp(start(j):fin(j)));
        event(argmax+start(j)-1) = true;
    end
    
    index = find(event); amp = amp(index);
    
    % get absolute times of events
    time = t(index); % + Ch3.times(((i-1)*batch+1));
    
    % extract wavelet transform around events
    wav = zeros(size(wcf,1),2*win_size,length(time));
    for j=1:length(time)
        start = index(j)-.5*win_size; fin = index(j)+1.5*win_size-1;
        if start > 0 && fin < size(wcf,2)
            wav(:,:,j) = wcf(:,start:fin);
        end
    end
    wavs = cat(3,wavs,wav);
    times = [times; time];
    amps = [amps amp];
end

% remove events that are in times we are not interested in

keep = true(size(times));
for i=1:size(seg,1)
    logical_temp = times < seg(i,1) | times > seg(i,2);
    keep = keep & logical_temp;
end

keep = keep & (amps < threshold_sp)';
delete = not(keep);
wavs(:,:,delete) = [];
times = times(keep);
amps = amps(keep);

toc

% ~4 min for the recording to run