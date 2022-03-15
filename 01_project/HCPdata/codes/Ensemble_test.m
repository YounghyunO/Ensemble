%% Add paths 
addpath(genpath('/local_raid1/01_software/toolboxes/cifti-matlab/')); %Add path to the gifti read
addpath(genpath('/local_raid1/01_software/toolboxes/matlab_util/')); % Add path to the utils

% Add caselist
listDir = fullfile('/local_raid1/03_user/younghyun/01_project/HCPdata');
list = load(fullfile(listDir,'/caselist_HCP.txt'));
list = censorsub(list); list(18) = 0; list(27) = 0; list(37) = 0; list(52) = 0; list(75) = 0; list = nonzeros(list);

load('EnsembleMdl2.mat')
load('F2390.mat')
% Fs = 1/0.72; % Sampling frequency (1/TR)
% T = 1/Fs; % Sampling period (TR)
% N = 1200; % Number of timepoints
% t = (0:N-1)*T; % Total sampling time, i.e., time-vector
% % Lower bound is DC (0 Hz) and the Upper bound is Nyquist Frequency,
% % that is, heoretically minimum number of points that one needs to sample in
% % order to see the real frequency of the sample. --> 
% hz = linspace(0,Fs/2,N/2+1); %The formula to convert indices to Hz 
% F = filtermumford(200,0.72,2390);


%%
tic
predicteds = [];
FCs = [];
rDCMs = []; 
VARs = [];
FASKs = [];
Patels = [];

for sub =1:size(list,1)
    Data = [];

    sub
    %%
    
%     % Glasser et al., 2016 atlas
%     dir_atlas = '/local_raid1/02_data/99_parcellation/Glasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k';
%     aname = fullfile(dir_atlas, '/Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
%     atlas = ft_read_cifti(aname).indexmax;
%     
%     dir_atals = '/local_raid1/03_user/younghyun/01_project/HCPdata/data';
%     aname = fullfile(dir_atals, 'Gordon333_FreesurferSubcortical.10k_fs_LR.dlabel.nii');
%     atlas = ft_read_cifti(aname)

    % Schaefer 2017 atlas
    dir_atlas = '/local_raid1/02_data/00_PARCELLATION/HCP/fslr32k/Schaefer2018_CerebCort_LocalGlobal';
    aname = fullfile(dir_atlas, '/Schaefer2018_100Parcels_7Networks_order.dlabel.nii');
    atlas = ft_read_cifti(aname).parcels;

    % Set routes for the files
    dir_struct = sprintf('/local_raid1/03_user/younghyun/01_project/HCPdata/Data/Structural/%d/MNINonLinear/fsaverage_LR32k',list(sub));
    dir_func1 = sprintf('/local_raid1/03_user/younghyun/01_project/HCPdata/Data/Functional/%d/MNINonLinear/Results/rfMRI_REST1_LR',list(sub));
    dir_func2 = sprintf('/local_raid1/03_user/younghyun/01_project/HCPdata/Data/Functional/%d/MNINonLinear/Results/rfMRI_REST2_RL',list(sub));
    
    % Load timeseries
    fname1 = fullfile(dir_func1, '/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii');
    fname2 = fullfile(dir_func2, '/rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii');
    
    ts_LR = read_cifti(fname1);
    ts_RL = read_cifti(fname2);
    
    % concat LR and RL
    ts_all = [zscore(ts_LR.cdata(:,6:end), 0, 2), zscore(ts_RL.cdata(:,6:end), 0, 2)];
    nvols = size(ts_all,2);
    
    % map to parcellation
    list_struct = ts_LR.diminfo{1, 1}.models;
    
    % index of  both cortices
    ctx_idx = list_struct{1}.start:list_struct{2}.start+list_struct{2}.count-1;
    
   % timeseries of cortices
    ts_ctx = ts_all(ctx_idx,:);


    % Load Atlas
    atlasroi = [gifti(fullfile(dir_struct, sprintf('%d.L.atlasroi.32k_fs_LR.shape.gii',list(sub)))).cdata; ...
                gifti(fullfile(dir_struct, sprintf('%d.R.atlasroi.32k_fs_LR.shape.gii',list(sub)))).cdata];

    % mean of the cortex timeseries within a parcel
    atlas(atlasroi == 0) = 0; atlas = nonzeros(atlas);
    
    ctx_roi =100;
    tsmean_ctx = zeros(ctx_roi, nvols);
    for roi = 1:ctx_roi
        tsmean_ctx(roi, :) = mean(ts_ctx(atlas==roi, :), 1);
    end
    
    % subcortex
    sub_roi = 19;
    tsmean_sub = zeros(sub_roi,nvols);
    for roi = 3:21
        start_vox = list_struct{roi}.start;
        end_vox = start_vox + list_struct{roi}.count -1;
        tsmean_sub(roi-2,:) = mean(ts_all(start_vox:end_vox,:),1);
    end
    tsmean_sub(5,:) = []; tsmean_sub(7:10,:) = []; 
    tsmean = [tsmean_ctx; tsmean_sub];

    % temporal filtering
    BOLD = F*tsmean'; 
    
    % Ensemble method
    
    writematrix(BOLD,'BOLD5');   
    FC = corr(BOLD); FC = FC - diag(diag(FC));
    FCs = cat(3,FCs,FC);
    FC = FC(:);
    
    [VAR] = VAR_run(BOLD');
    VAR = VAR - diag(diag(VAR));
    VARs = cat(3,VARs,VAR);
       
    Y.y = BOLD;
    Y.dt = 0.72;
    rdcm = tapas_rdcm_model_specification(Y,[],[]);
    [output,~] = tapas_rdcm_estimate(rdcm,'r',[],1);
    rDCM = output.Ep.A; 
    rDCMs = cat(3,rDCMs,rDCM);
    rDCM = rDCM - diag(diag(rDCM));
    
    oAlpha = 0.5;
    command = sprintf('java  -jar causal-cmd-1.3.0-jar-with-dependencies.jar --algorithm fask --test fisher-z-test --orientationAlpha %d --faskAdjacencyMethod 1 --alpha 0.01 --data-type continuous --dataset BOLD5.txt --delimiter comma --no-header --prefix FASK5', oAlpha);
    system(command)
    FASK = Tetrad2Matrix('FASK5.txt','directed')';
    FASKs = cat(3,FASKs,FASK);
    
    alpha = 1*10^-3; 
    command = sprintf('java -jar causal-cmd-1.3.0-jar-with-dependencies.jar --algorithm fas --stableFAS --test fisher-z-test --alpha %d --data-type continuous --dataset BOLD5.txt --delimiter comma --no-header --prefix FAS5',alpha);
    system(command)
    FAS = Tetrad2Matrix('FAS5.txt','undirected');
    [patel] = patel_tau(BOLD,FAS);
    Patels = cat(3,Patels,patel);
    
    % vectorize output matrices
    rDCM = rDCM(:);
    VAR = VAR(:);
    FASK = FASK(:);
    patel = patel(:);        

    % output matrix for Ensemble learning
    Data = [Data; rDCM FASK VAR FC patel];
    temp = predict(EnsembleMdl,Data);
    predicted = reshape(temp,114,114);
    
    predicteds = cat(3,predicteds,predicted);
    
end
toc

% save('Ensemble_MMP360_test.mat')
% EC = sum(predicteds,3)/85;
% save('EC_MMP360_session1.mat','EC');