% % Add paths 
% addpath(genpath('/local_raid1/01_software/toolboxes/cifti-matlab/')); %Add path to the gifti read
% addpath(genpath('/local_raid1/01_software/toolboxes/matlab_util/')); % Add path to the utils
% % addpath('/local_raid1/03_user/younghyun/03_software/fieldtrip-master');
% addpath('/local_raid1/03_user/younghyun/03_software/surfstat');

% Add caselist
listDir = fullfile('/local_raid1/03_user/younghyun/01_project/HCPdata');
list = load(fullfile(listDir,'/caselist_HCP.txt'));
list = censorsub(list);

% PSD_total = [];
% FC_total = [];
BOLD_total = cell(1,85);
% Glasser et al., 2016 atlas
dir_atlas = '/local_raid1/02_data/99_parcellation/Glasser2016/Q1-Q6_RelatedParcellation210/MNINonLinear/fsaverage_LR32k';
aname = fullfile(dir_atlas, '/Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
atlas = ft_read_cifti(aname).indexmax;

% %Schaefer400 atlas
% dir_atlas = '/local_raid1/02_data/00_PARCELLATION/HCP/fslr32k/Schaefer2018_CerebCort_LocalGlobal';
% aname = fullfile(dir_atlas, '/Schaefer2018_100Parcels_7Networks_order.dlabel.nii');
% atlas = ft_read_cifti(aname).parcels;

Fs = 1/0.72; % Sampling frequency (1/TR)
T = 1/Fs; % Sampling period (TR)
N = 1200; % Number of timepoints
t = (0:N-1)*T; % Total sampling time, i.e., time-vector
% Lower bound is DC (0 Hz) and the Upper bound is Nyquist Frequency,
% that is, heoretically minimum number of points that one needs to sample in
% order to see the real frequency of the sample. --> 
hz = linspace(0,Fs/2,N/2+1); %The formula to convert indices to Hz 
F = filtermumford(200,0.72,1200);

for sub = 1:size(list,1)
    Data = [];

    sub
    %%
    % Set routes for the files
    dir_struct = sprintf('/local_raid1/03_user/younghyun/01_project/HCPdata/Data/Structural/%d/MNINonLinear/fsaverage_LR32k',list(sub));
    dir_func1 = sprintf('/local_raid1/03_user/younghyun/01_project/HCPdata/Data/Functional/%d/MNINonLinear/Results/rfMRI_REST1_LR',list(sub));
    
    % Load timeseries
    fname1 = fullfile(dir_func1, '/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii');
%     fname2 = fullfile(dir_func2, '/rfMRI_REST1_RL_Atlas_hp2000_clean.dtseries.nii');

    ts1 = ft_read_cifti(fname1);
    nvols = length(ts1.time);
    ts1 = ts1.dtseries((ts1.brainstructure == 1) | (ts1.brainstructure == 2), :);

    % Load Atlas

    atlasroi = [gifti(fullfile(dir_struct, sprintf('%d.L.atlasroi.32k_fs_LR.shape.gii',list(sub)))).cdata; ...
                gifti(fullfile(dir_struct, sprintf('%d.R.atlasroi.32k_fs_LR.shape.gii',list(sub)))).cdata];

    % Averaging
    atlas(atlasroi == 0) = 0;

    tsmean = zeros(180, nvols);
    for roi = 1:180
        tsmean(roi, :) = mean(ts1(atlas==roi, :), 1);
    end
    tsmean = tsmean';
    BOLD = F*tsmean;

    fhat = fft(tsmean);
    P2 = abs(fhat);
    P1 = P2(1:1200/2+1,:);
    P1(2:end-1,:) = P1(2:end-1,:).^2;   
%     figure; plot(hz,P1)
    PSD_total = cat(3,PSD_total,P1);


end

%%
total_dissiml = [];
total_dissimr = [];
for i = 1:85
    P1l = squeeze(PSD_total(:,:,i));
%     P1r = squeeze(PSD_total(:,201:400,i));
    corrdistl = squareform(pdist(P1l','correlation'));
%     corrdistr = squareform(pdist(P1r','correlation'));
    total_dissiml = cat(3,total_dissiml,corrdistl);
%     total_dissimr = cat(3,total_dissimr,corrdistr);
end

avg_dissiml = mean(total_dissiml,3);
avg_dissimr = mean(total_dissimr,3);
Zl = linkage(avg_dissiml,'complete');
Zr = linkage(avg_dissimr,'complete');
Tl = cluster(Zl,'maxclust',7);
Tr = cluster(Zr,'maxclust',6);
cutoff = median(Zl(end-8,3));

figure; dendrogram(Zl,200,'ColorThreshold',cutoff); xtickangle(45); set(gca,'FontSize',5);
figure; dendrogram(Zr,0); xtickangle(45)
figure; imagesc(avg_dissiml);
figure; plot(hz,P1(:,316)); ylim([0 5*10^4]); title('parcel316 Cluster 10')

% Mapping the values to the brain

tempL = gifti('/local_raid1/01_software/HCPpipelines/global/templates/standard_mesh_atlases/cortical_surface/S900.L.midthickness_MSMAll.10k_fs_LR.surf.gii');
tempR = gifti('/local_raid1/01_software/HCPpipelines/global/templates/standard_mesh_atlases/cortical_surface/S900.R.midthickness_MSMAll.10k_fs_LR.surf.gii');

surfL.coord = tempL.vertices';
surfL.tri = tempL.faces;
surfR.coord = tempR.vertices';
surfR.tri = tempR.faces;
surf.coord = [ tempL.vertices' tempR.vertices' ];
surf.tri = [ tempL.faces; tempR.faces+32492; ];

% figure; project_detection_community(lmetric, 0, atlas(1:32492)-180, surfL, 1);
% figure; project_detection_community(rmetric, 0, atlas(32493:end), surfR, 1);
project_detection_community(Tl', 0, atlas', surf, 1);

