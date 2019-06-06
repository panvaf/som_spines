% check for preferential occurence of events after odor and call
% presentations

% load files

cellname = 'cella260319';
recording = strcat(cellname,'.mat');
load(recording)
stimulation = strcat(cellname,'stim.mat');
load(stimulation)
events = strcat(cellname,'_th_syn.mat');
load(events)

% parameters
dur = Ch3.times(end);
odor = 5; % in s
call = 1; % in s
Odor = cat(1,LemOd,MomOd,NonMomOd,NonSibOd,SibOd);
Call = cat(1,MomCall,NonMomCall,NonSibCall,SibCall);

% subtract time of segments that were not analyzed
t_seg = sum(diff(seg,1,2));
dur = dur - t_seg;

% total duration of presentations
odor_dur = odor*size(Odor,1);
call_dur = call*size(Call,1);

% proportion of recording during a certain presentation was active
odor_t_prop = odor_dur/dur;
call_t_prop = call_dur/dur;

% proportion of events happening during a presentation
odor_th_prop = size(t_od_th,1)/size(Odor,1);
odor_syn_prop = size(t_od_syn,1)/size(Odor,1);
call_th_prop = size(t_call_th,1)/size(Call,1);
call_syn_prop = size(t_call_syn,1)/size(Call,1);

% relative tendency of events to happen during presentations
odor_th_pref = odor_th_prop/odor_t_prop
odor_syn_pref = odor_syn_prop/odor_t_prop
call_th_pref = call_th_prop/call_t_prop
call_syn_pref = call_syn_prop/call_t_prop