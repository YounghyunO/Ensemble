dir_atlas = '/local_raid1/03_user/younghyun/02_data/parcellations';
aname = fullfile(dir_atlas,'Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
atlas = ft_read_cifti(aname).indexmax;


% dir_struct = sprintf('/local_raid1/03_user/younghyun/02_data/structure/%d/MNINonLinear/fsaverage_LR32k',uplist(1));

% Load Atlas
% atlasroi = [gifti(fullfile(dir_struct, sprintf('%d.L.atlasroi.32k_fs_LR.shape.gii',uplist(1)))).cdata; ...
%             gifti(fullfile(dir_struct, sprintf('%d.R.atlasroi.32k_fs_LR.shape.gii',uplist(1)))).cdata];
atlasroi2 = [gifti(fullfile('sub-NDARAD615WLJ.L.atlasroi.32k_fs_LR.shape.gii')).cdata; ...
                gifti(fullfile('sub-NDARAD615WLJ.R.atlasroi.32k_fs_LR.shape.gii')).cdata];

% mean of the cortex timeseries within a parcel
atlas(atlasroi == 0) = 0; 
rois = 360;

fid = fopen('log.txt');
txt = textscan(fid,'%s','delimiter','\n');
temp = txt{1};
FCs = cell(246,1); rDCMs = cell(246,1);

parfor i = 1:246
    [FC,rDCM] = movie_fmri(temp,i,atlas);
    FCs{i} = FC; rDCMs{i} = rDCM;
end

save('group_rDCM.mat','rDCMs')
save('group_FC.mat','FCs')