% receive and manipulate wavelet transform of physiological recordings

tic
% load data

load('cellg110319.mat')
volt = Ch3.values; % in mV
t = Ch3.times;     % in s
curr = Ch6.values; % in nA
samplefreq = 1/Ch3.interval;  % in Hz

% remove segments that are not wanted (stimulation etc)
% not needed since the wavelet algorithm is ultra fast
% better remove events in the end

seg = [1 100; 550 750;];  % segments to be removed, in seconds
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
win = .1; % in s, maximum expected size of postsynaptic event
win_size = floor(win*samplefreq);
threshold = 1000; % for detection of events in general

display = 0;   % 1 if want to display
show = [];
waveFrq = [10,300];       % Transform frequency range (from Logothetis paper)
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
    amp = movmean(sum(wcf),win_size);
    event = diff(amp>threshold);
    event(event<0) = 0; event = logical(event);
    
    % required space between events
    for j=1:length(event)
        if event(j)
            event((j+1):(j+1+win_size/5)) = false;
            j = j + win_size;
        end
    end
    
    index = find(event); amp = amp(event);
    
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
for i=1:size(seg,2)
    logical_temp = times < seg(i,1) | times > seg(i,2);
    keep = keep & logical_temp;
end

delete = not(keep);
wavs(:,:,delete) = [];
times = times(keep);
amps = amps(keep);

toc

% ~4 min for the recording to run