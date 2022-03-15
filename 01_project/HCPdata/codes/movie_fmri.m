function [FC,rDCM] =movie_fmri(temp,i,atlas)

sub = temp{i};
dir_func = sprintf('/local_raid1/03_user/younghyun/02_data/mri/%s/ciftify/%s/MNINonLinear/Results/task-movieDM',sub,sub);
fname = fullfile(dir_func,'task-movieDM_Atlas_s3.dtseries.nii');

tstemp = ft_read_cifti(fname);
ts = zscore(tstemp.dtseries(:,6:end),0,2);

nvols = size(ts,2);
tsmean = zeros(rois,nvols);
for roi = 1:rois
    tsmean(roi, :) = mean(ts(atlas==roi, :), 1);
end

BOLD = zeros(nvols,rois);
for i = 1:360
    id = sorted(i);
    BOLD(:,i) = tsmean(id,:);
end

BOLD = detrend(BOLD);

FC = corr(BOLD); FC = FC - diag(diag(FC));
figure; imagesc(FC);

Y.y = BOLD;
Y.dt = 0.8;
rdcm = tapas_rdcm_model_specification(Y,[],[]);
[output,~] = tapas_rdcm_estimate(rdcm,'r',[],1);
rDCM = output.Ep.A; 
rDCM = rDCM - diag(diag(rDCM)); 