% load data from excel file of events during single cell recordings

filename = 'SiblingTimes.xls';
cellname = 'celle070319';
directory = '';
filepath = strcat(directory,'/',filename);
savefile = strcat(cellname,'stim.mat');

data = readtable(filepath);
SibOd = table2array(data(:,2)); SibOd(find(isnan(SibOd))) = [];
NonSibOd = table2array(data(:,4)); NonSibOd(find(isnan(NonSibOd))) = [];
LemOd = table2array(data(:,6)); LemOd(find(isnan(LemOd))) = [];
MomOd = table2array(data(:,8)); MomOd(find(isnan(MomOd))) = [];
NonMomOd = table2array(data(:,10)); NonMomOd(find(isnan(NonMomOd))) = [];
SibCall = table2array(data(:,12)); SibCall(find(isnan(SibCall))) = [];
NonSibCall = table2array(data(:,14)); NonSibCall(find(isnan(NonSibCall))) = [];
MomCall = table2array(data(:,16)); MomCall(find(isnan(MomCall))) = [];
NonMomCall = table2array(data(:,18)); NonMomCall(find(isnan(NonMomCall))) = [];

save(savefile,'SibOd','NonSibOd','LemOd','MomOd','NonMomOd','SibCall','NonSibCall','MomCall','NonMomCall')