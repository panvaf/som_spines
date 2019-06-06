% Compute coherence between respiratory rhythm and rhythms in single cell
% recordings
tic

% load data

cellname = 'cella270319';
recording = strcat(cellname,'.mat');
load(recording)
stimulation = strcat(cellname,'stim.mat');
load(stimulation)
volt = Ch3.values(1:end-1); % in mV
t = Ch3.times(1:end-1);     % in s
sampleint = Ch3.interval;
samplefreq = 1/sampleint;  % in Hz
resp = Ch7.values; % in mV

% parameters
odor_dur = 5; % in s
call_dur = 1; % in s
f = 1:10;    % frequencies for which to compute coherence
win_size = 5;  % window in which to compute coherence
overlap = 4;  % 80 % overlap = moves 1 s per computation

window = floor(win_size*samplefreq);
noverlap = floor(overlap*samplefreq);

% Look for coherence during odor presentations, calls and baseline

Odor = cat(1,LemOd,MomOd,NonMomOd,NonSibOd,SibOd);
Call = cat(1,MomCall,NonMomCall,NonSibCall,SibCall);

t_od = []; t_call = [];

for i=1:size(Odor,1)
    temp = t(t>Odor(i) & t<(Odor(i)+odor_dur));
    t_od = [t_od; temp];
end

for i=1:size(Call,1)
    temp = t(t>Call(i) & t<(Call(i)+call_dur));
    t_call = [t_call; temp];
end

t_rem = cat(1,t_od,t_call);
t_base = setdiff(t,t_rem);

% convert times to indices
odor = round(t_od*samplefreq);
call = round(t_call*samplefreq);
base = round(t_base*samplefreq); base = base(2:end);

% receive averaged coherences

[cxy_odor,f] = mscohere(resp(odor),volt(odor),window,noverlap,f,samplefreq);
[cxy_call,f] = mscohere(resp(call),volt(call),window,noverlap,f,samplefreq);
[cxy_base,f] = mscohere(resp(base),volt(base),window,noverlap,f,samplefreq);

% plot

figure
plot(f,cxy_odor)
hold on
plot(f,cxy_call)
plot(f,cxy_base)
legend('Odors','Calls','Baseline')
title('Respiration - single cell recording coherence')
ylabel('Coherence')
xlabel('Frequency (Hz)')

toc