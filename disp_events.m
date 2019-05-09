% Display wavelet transform around events from a given recording

% load data

load('cellg110319.mat')
volt = Ch3.values; % in mV
t = Ch3.times;     % in s
curr = Ch6.values; % in nA
samplefreq = 1/Ch3.interval;  % in Hz
load('wav_trans.mat')
win = .1; % in s, maximum expected size of postsynaptic event
win_size = floor(win*samplefreq);

% plots wavelet transform: click on current figure to get the next

for i = 1:size(wavs,3)
    t = (times(i) - 0.05):Ch3.interval:(times(i) + 0.15 - Ch3.interval);
    record = volt(((times(i) - 0.05)*samplefreq):((times(i) - 0.05)*samplefreq + 2*win_size - 1));
    figure
    wscalogram('image',flip(wavs(:,:,i),1),'scales',flip(round(pfreq,2),2),'xdata',t,'ydata',record);
    xlabel('Time (s)'); ylabel('Freq (Hz)');
    title('Pseudofrequency-time transform');
    w = waitforbuttonpress;
end