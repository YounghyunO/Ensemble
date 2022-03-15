function atlas = atlasinvest(uplist,parcellation)

switch parcellation
    case 'MMP'
    % Glasser et al., 2016 atlas
    dir_atlas = '/local_raid1/03_user/younghyun/02_data/parcellations';
%     aname = fullfile(dir_atlas,'Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
    aname = fullfile(dir_atlas, '/Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S4.dlabel.nii');
    atlas = ft_read_cifti(aname).indexmax;
    case 'Schaefer'
    %     Schaefer 2017 atlas
    dir_atlas = '/local_raid1/03_user/younghyun/02_data/parcellations';
    aname = fullfile(dir_atlas, '/Schaefer2018_400Parcels_7Networks_order_Tian_Subcortex_S4.dlabel.nii');
    atlas = ft_read_cifti(aname).parcels;
%     atlas(32493:end) = atlas(1:32492) +50;
    case 'Yeo'
    %Yeo 7 atlas
    dir_atlas = '/local_raid1/03_user/younghyun/02_data/parcellations';
    aname = fullfile(dir_atlas, 'Yeo2011_7Networks_N1000.dlabel.nii');
    atlas = ft_read_cifti(aname).parcels;
    temp = atlas(32493:end);
    indx = temp > 0;
    temp(indx,1) = temp(indx,1) +7;
    atlas = [atlas(1:32492,1);temp];
end

dir_struct = sprintf('/local_raid1/03_user/younghyun/02_data/structure/%d/MNINonLinear/fsaverage_LR32k',uplist(1));

% Load Atlas
atlasroi = [gifti(fullfile(dir_struct, sprintf('%d.L.atlasroi.32k_fs_LR.shape.gii',uplist(1)))).cdata; ...
            gifti(fullfile(dir_struct, sprintf('%d.R.atlasroi.32k_fs_LR.shape.gii',uplist(1)))).cdata];

% mean of the cortex timeseries within a parcel
atlas(atlasroi == 0) = 0; 
% atlas = nonzeros(atlas);

end