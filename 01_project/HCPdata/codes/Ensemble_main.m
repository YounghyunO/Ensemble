%% setup 
% go to the run directory where the necessary files are stored
cd /local_raid1/03_user/younghyun/01_project/HCPdata/run

% delete previously stored files for variables
delete BOLD* FASK* FAS*

% Add paths 
addpath(genpath('/local_raid1/01_software/toolboxes/toolboxes/cifti-matlab/')); %Add path to the gifti read
addpath(genpath('/local_raid1/01_software/toolboxes/toolboxes/matlab_util/')); % Add path to the utils

 
uplist = load('caselist_replication.txt'); % load subject list based on exclusion criteria
load('EnsembleMdl(210112).mat') % load Ensemble matrix
load('sorted.mat')

listsize = size(uplist,1);
predicted = cell(1,listsize);
rDCMs = cell(1,listsize); VARs = cell(1,listsize);
FASKs = cell(1,listsize); GCs = cell(1,listsize);
FCs = cell(1,listsize); BOLDs = cell(1,listsize);

% load atlas
parcel = 'MMP'; 
switch parcel
    case 'Schaefer'
        rois = 454;
        % ctx = 400;
    case 'MMP'
        rois = 414;
end
atlas = atlasinvest(uplist,parcel);

%% run Ensemble

for i = 102:listsize
    [pred,rDCM,FASK,VAR,GC,FC,BOLD] = Ensemble(i,uplist,EnsembleMdl,atlas,rois,sorted);
     predicted{i} = pred; rDCMs{i} = rDCM; FASKs{i} = FASK;
     VARs{i} = VAR; GCs{i} = GC; FCs{i} = FC; BOLDs{i} = BOLD;
     
end

cd /local_raid1/03_user/younghyun/01_project/HCPdata/results/MMP414_rep

save('Ensembles.mat','predicted')
save('FCs.mat','FCs')
save('rDCMs.mat','rDCMs')
save('FASKs.mat','FASKs')
save('VARs.mat','VARs')
save('GCs.mat','GCs')
save('BOLDs.mat','BOLDs')