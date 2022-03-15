%% Read dtseries & concat

% phase: L->R file
% ts_temp = cifti_read('data/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii');
list_struct = ts_RL.diminfo{1, 1}.models;

EC = load('data/EC119.mat').temp5;
pos_os = sum(max(EC, 0), 1);    % EC의 excitatory 부분만 row sum
neg_os = sum(min(EC, 0), 1);    % EC의 inhibitory 부분만 row sum


%% Visualization subcortical areas

% [Map to SubCortex]
% template = load_nii('data/MNI152_T1_2mm_brain.nii');

subcortex = zeros(91, 109, 91);
for roi=3:21
     disp(list_struct{roi}.struct)
     vox_idx = list_struct{roi}.voxlist+1;
     for vox = vox_idx
        subcortex(vox(1), vox(2), vox(3)) = pos_os(100 + roi - 2);
     end
end

volumeViewer


    
