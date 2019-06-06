% quantify theta modulation around stimulation events using wavelet
% transform

% load data

cellname = 'cella270319';
recording = strcat(cellname,'.mat');
load(recording)
stimulation = strcat(cellname,'stim.mat');
load(stimulation)
volt = Ch3.values; % in mV
t = Ch3.times;     % in s
sampleint = Ch3.interval;
samplefreq = 1/sampleint;  % in Hz

% parameters

stimulus = 'call';
times = NonSibCall;
display = 0;   % 1 if want to display
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
P_stim = 0; P_base = 0;

for i=1:size(times,1)
    index = times(i)*samplefreq;
    signal = volt((index-bef_size):(index+stim_size+aft_size));
    t = (times(i) - bef):sampleint:(times(i) + stim + aft);
    [wcf, pfreq, scales] = wavtrans(signal,t,samplefreq,rowsPerOct,waveFrq,padmode,wavelet,show);
    amps(i,1) = sum(wcf(:,1:bef_size),'all')/bef;   % amplitude before stimulation
    amps(i,2) = sum(wcf(:,(bef_size+1):(bef_size+stim_size)),'all')/stim;   % amplitude during stimulation
    amps(i,3) = sum(wcf(:,(bef_size+stim_size+1):end),'all')/aft;   % amplitude after stimulation
    
    % plot modulation in fourier transform
    
    base = signal(1:bef_size);
    stm = signal((bef_size+1):(bef_size+stim_size));
    
    ft_temp_stim = fft(stm);
    L_stim = size(stm,1);
    P2 = abs(ft_temp_stim/L_stim);
    P1_temp = P2(1:L_stim/2+1);
    P1_temp(2:end-1) = 2*P1_temp(2:end-1);
    P_stim = P_stim + P1_temp;
    
    ft_temp_base = fft(base);
    L_base = size(base,1);
    P2 = abs(ft_temp_base/L_base);
    P1_temp = P2(1:L_base/2+1);
    P1_temp(2:end-1) = 2*P1_temp(2:end-1);
    P_base = P_base + P1_temp;
    
    w = waitforbuttonpress;
end

f = samplefreq*(0:(L_stim/2))/L_stim;
f_base = samplefreq*(0:(L_base/2))/L_base;
P_stim = P_stim/size(times,1); P_base = P_base/size(times,1);

figure
plot(amps)
xlabel('Event #')
ylabel('Intensity of theta rhythm')
legend('Baseline','Stimulation','Rebound')
title('Modulation of theta rhythm by stimulation (Non-Sibling call)')

% use one-sample Kolmogorov-Smirnov test to test if distributions are normal

norm_bef = ~kstest(amps(:,1));
norm_stim = ~kstest(amps(:,2));
norm_aft = ~kstest(amps(:,3));

% use paired-sample t-test if distributions are normal

if norm_bef && norm_stim
    disp('Normally distributed:')
    if ttest(amps(:,1),amps(:,2))
        disp('Different mean for before and after stimulation')
    else
        disp('Same mean for before and after stimulation')
    end
end

if norm_bef && norm_aft
    disp('Normally distributed:')
    if ttest(amps(:,1),amps(:,3))
        disp('Different mean for before and after stimulation')
    else
        disp('Same mean for before and after stimulation')
    end
end

% use Wilcoxon signed rank test if distributions are not normal

[p1,h1] = signrank(amps(:,1),amps(:,2));

if h1
    disp('Different mean for before and during stimulation')
else
    disp('Same mean for before and during stimulation')
end
disp(p1)

[p2,h2] = signrank(amps(:,1),amps(:,3));

if h2
    disp('Different mean for before and after stimulation')
else
    disp('Same mean for before and after stimulation')
end
disp(p2)


% see difference in fourier transform

tsin = timeseries(P_base,f_base);
tsout = resample(tsin,f);
P_base = tsout.Data;

figure
semilogy(f,P_base)
hold on
semilogy(f,P_stim)
xlim([1 50])
xlabel('Frequency (Hz)')
ylabel('Power')
legend('Baseline','Stimulation')
title('Fourier transform (Sibling call)')