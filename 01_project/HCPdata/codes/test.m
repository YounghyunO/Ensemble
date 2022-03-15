%% Add paths 
addpath(genpath('/local_raid1/01_software/toolboxes/toolboxes/cifti-matlab/')); %Add path to the gifti read


% %  Load surface gifti file
L_surf = gifti('/local_raid1/03_user/younghyun/02_data/parcellations/Q1-Q6_RelatedParcellation210.L.midthickness_MSMAll_2_d41_WRN_DeDrift.32k_fs_LR.surf.gii');
R_surf = gifti('/local_raid1/03_user/younghyun/02_data/parcellations/Q1-Q6_RelatedParcellation210.R.midthickness_MSMAll_2_d41_WRN_DeDrift.32k_fs_LR.surf.gii');

% L_surf = gifti('/local_raid1/03_user/younghyun/02_data/parcellations/S900.L.midthickness_MSMAll.32k_fs_LR.surf.gii');
% R_surf = gifti('/local_raid1/03_user/younghyun/02_data/parcellations/S900.R.midthickness_MSMAll.32k_fs_LR.surf.gii');

surf.coord = [L_surf.vertices', R_surf.vertices'];
surf.tri = [L_surf.faces; R_surf.faces+32492];
% surf.tri = [L_surf.faces; R_surf.faces+10242];

% % Load shaefer parcellation
% parcel = ft_read_cifti('/local_raid1/03_user/younghyun/02_data/parcellations/Schaefer2018_100Parcels_7Networks_order.dlabel.nii').parcels;


% % Glasser et al., 2016 atlas
% dir_atlas = '/local_raid1/03_user/younghyun/02_data/parcellations/';
% aname = fullfile(dir_atlas, '/Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
% parcel = ft_read_cifti(aname).indexmax;

%%
% Load effective connectiviy
% EC = load('EC_MMP374.mat').EC;

EC = EC - diag(diag(EC)); 
ECvec = EC(:); 
% ECvec(ECvec<-0.02) = -1; ECvec(ECvec>0.02) = 1;
% ECvec(abs(ECvec)<0.01) = 0;
% ECmod = reshape(ECvec,114,114);
EC_ctx = ECmod(1:100,1:100);
EC_sub = ECmod(101:end,101:end);

% for cortex
ctx_pos_of = sum(max(EC1p, 0), 1);  
ctx_neg_of = -sum(min(EC1p, 0), 1);    
ctx_pos_if = sum(max(EC1p, 0), 2);
ctx_neg_if = -sum(min(EC1p, 0), 2);

ctx_of = sum(abs(EC), 1);  
ctx_if = sum(abs(EC), 2);


% for subcortex
sub_pos_of = sum(max(EC_sub,0), 1);
sub_neg_of = -sum(min(EC_sub,0), 1);
sub_pos_if = sum(max(EC_sub,0), 2);
sub_neg_if = -sum(min(EC_sub,0), 2);

pos_info = zeros(size(EC,1),1);
neg_info = zeros(size(EC,1),1);
% for information flow
for i = 1:size(EC,1)
    pos_info(i,1) = sum(max(ECmod(i,:),0))/sum(max(ECmod(:,i),0));
    neg_info(i,1) = sum(min(ECmod(i,:),0))/sum(min(ECmod(:,i),0));
end

pos_info_ctx = pos_info(1:100,1);
% pos_info_sub = pos_info(361:end,1);
neg_info_ctx = neg_info(1:100,1);
% neg_info_sub = neg_info(361:end,1);

% Visualization
close all
% project parcellation to vertex
% target = 112;
% target_of = ECmod(:,target); 
% target_if = ECmod(target,:);
% target_if(1,101) = 1;
% target_of(target,1) = 1;
% target_if(1,target) = 1;

% feature = target_if(1,1:360);
% feature2 = target_if(1,361:374);
% 
% feature = target_of(1:100,1);
% feature2 = target_of(101:114,1);

% % feature for visualization
feature = en_neg_of;
% feature = neg_info_ctx;
% feature2 = sub_pos_of;

ref = zeros(1, 64984);
for roi = 1:100
    ref(Schaefer == roi) = feature(roi);
end

% ref = zeros(1, 64984);
% for roi = 1:100
%     ref(parcel == roi) = feature(roi);
% end

figure; plot_cortical(ref, surf); 
colormap(parula); caxis([0 2.5])
saveas(gca,'neg_of.png')
close
% saveas(gca,'mmp_pos_if_ctx(sum).png')
% figure; plot_subcortical(feature2,'cmap','jet','ventricles','False'); 
% saveas(gca,'mmp_pos_if_sub(sum).png')


Schaefer(32493:end) = Schaefer(1:32492) +50;

figure;
SurfStatView(Schaefer, surf);
colormap(lines(400))
% color = 'jet';   % 위의 colormap 중에 선택
% SurfStatColormap( color );
% saveas(gca','pos_if.png')