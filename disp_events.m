% Display wavelet transform around events from a given recording

% load data

load('cellf070319.mat')
volt = Ch3.values; % in mV
t = Ch3.times;     % in s
curr = Ch6.values; % in nA
sampleint = Ch3.interval;
samplefreq = 1/sampleint;  % in Hz
load('cellf070319_th_w5th500.mat')
win = 5; % in s, maximum expected size of postsynaptic event
win_size = floor(win*samplefreq);

% sort with descending order
[amps, order] = sort(amps,'descend');
wavs = wavs(:,:,order); times = times(order);

% plots wavelet transform: click on current figure to get the next

for i = 1:size(wavs,3)
    t = (times(i) - 0.5*win):sampleint:(times(i) + 1.5*win - sampleint);
    record = volt(((times(i) - 0.5*win)*samplefreq):((times(i) - 0.5*win)*samplefreq + 2*win_size - 1));
    figure
    wscalogram('image',flip(wavs(:,:,i),1),'scales',flip(round(pfreq,2),2),'xdata',t,'ydata',record);
    xlabel('Time (s)'); ylabel('Freq (Hz)');
    title('Pseudofrequency-time transform');
    w = waitforbuttonpress;
end