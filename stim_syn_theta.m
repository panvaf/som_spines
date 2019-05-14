% check simultaneous occurence of theta oscillations and synaptic events 
% during stimulation

% load data

cellname = 'cellg110319';
stimulation = strcat(cellname,'stim.mat');
load(stimulation)
load('theta_w5th300.mat')
t_th = times;
load('wav_trans_w100th150r50-100.mat')
t_syn = times;

% parameters
odor_dur = 5; % in s
call_dur = 1; % in s

Odor = cat(1,LemOd,MomOd,NonMomOd,NonSibOd,SibOd);
Call = cat(1,MomCall,NonMomCall,NonSibCall,SibCall);

t_od_th = []; t_call_th = []; t_od_syn = []; t_call_syn = [];

for i=1:size(Odor,1)
    temp_th = t_th(t_th>Odor(i) & t_th<(Odor(i)+odor_dur));
    temp_syn = t_syn(t_syn>Odor(i) & t_syn<(Odor(i)+odor_dur));
    t_od_th = [t_od_th; temp_th]; t_od_syn = [t_od_syn; temp_syn];
end

for i=1:size(Call,1)
    temp_th = t_th(t_th>Call(i) & t_th<(Call(i)+call_dur));
    temp_syn = t_syn(t_syn>Call(i) & t_syn<(Call(i)+call_dur));
    t_call_th = [t_call_th; temp_th]; t_call_syn = [t_call_syn; temp_syn];
end

savefile = strcat(cellname,'_th_syn.mat');
save(savefile,'t_od_th','t_call_th','t_od_syn','t_call_syn')