function [predicted,rDCM,FASK,VAR,GC,FC,BOLD] = Ensemble(i,uplist,EnsembleMdl,atlas,rois,sorted)

    sub = i;
    
    % Set routes for the files
    dir_func1 = sprintf('/local_raid1/03_user/younghyun/02_data/resting/%d/MNINonLinear/Results/rfMRI_REST1_LR',uplist(sub));
    dir_func2 = sprintf('/local_raid1/03_user/younghyun/02_data/resting/%d/MNINonLinear/Results/rfMRI_REST1_RL',uplist(sub));
%     dir_func3 = sprintf('/local_raid1/03_user/younghyun/02_data/resting/%d/MNINonLinear/Results/rfMRI_REST2_LR',uplist(sub));
%     dir_func4 = sprintf('/local_raid1/03_user/younghyun/02_data/resting/%d/MNINonLinear/Results/rfMRI_REST2_RL',uplist(sub));
    
    % Load timeseries
    fname1 = fullfile(dir_func1, '/rfMRI_REST1_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii');
    fname2 = fullfile(dir_func2, '/rfMRI_REST1_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii');
%     fname3 = fullfile(dir_func3, '/rfMRI_REST2_LR_Atlas_MSMAll_hp2000_clean.dtseries.nii');
%     fname4 = fullfile(dir_func4, '/rfMRI_REST2_RL_Atlas_MSMAll_hp2000_clean.dtseries.nii');

    ts1temp = ft_read_cifti(fname1);
    ts2temp = ft_read_cifti(fname2);


    ts1 = zscore(ts1temp.dtseries(:,6:end),0,2);
    ts2 = zscore(ts2temp.dtseries(:,6:end),0,2);

    
    ts = [ts1, ts2];
    nvols = size(ts,2);
    tsmean = zeros(rois,nvols);
    for roi = 1:rois
        tsmean(roi, :) = mean(ts(atlas==roi, :), 1);
    end
    
    if rois == 414
        BOLD = zeros(nvols,rois);
        tsmean_ctx = tsmean(55:414,:);
        for i = 1:360
            id = sorted(i);
            BOLD(:,i) = tsmean_ctx(id,:);
        end
        tsmean_sub = tsmean(1:54,:);
        for i = 361:414
            BOLD(:,i) = tsmean_sub(i-360,:);
        end
    else
        BOLD = zeros(nvols,rois);
        BOLD(:,1:400) = tsmean(55:end,:)';
        BOLD(:,401:end) = tsmean(1:54,:)';
    end
    
    FC = corr(BOLD); FC = FC - diag(diag(FC));
%     figure; imagesc(FC);
    % Ensemble method
%     BOLD = CBIG_bandpass_matrix(BOLD,0.008,0.08,0.72);
    
    writematrix(BOLD,sprintf('BOLD%d',sub));   
    [VAR,GC] =StateSpaceGC(BOLD',0);
    
    Y.y = BOLD;
    Y.dt = 0.72;
    rdcm = tapas_rdcm_model_specification(Y,[],[]);
    [output,~] = tapas_rdcm_estimate(rdcm,'r',[],1);
    rDCM = output.Ep.A; 
    rDCM = rDCM - diag(diag(rDCM)); 

    FASK = subsample_wrapper(BOLD,i,0.05,0.00005);
        
    % vectorize output matrices
    rDCMvec = rDCM(:);
    VARvec = VAR(:);
    FASKvec = FASK(:);
    GCvec = GC(:);
    
    % output matrix for Ensemble learning
    Data = [rDCMvec VARvec GCvec FASKvec];
    temp = predict(EnsembleMdl,Data);
    predicted = reshape(temp,rois,rois);
        