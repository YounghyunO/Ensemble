%% script for information flow mapping 
infall = []; pos_infall = []; pos_outfall = [];
outfall = []; neg_infall = []; neg_outfall = [];
    vis = [1:9,51:58];
    sm = [10:15,59:66];
    da = [16:23,67:73];
    va = [24:30,74:78];
    lim = [31:33,79:80];
    fp = [34:37,81:89];
    dm = [38:50,90:100];

    seed = [119, 299];
for i = 1:202
    EC = cell2mat(Ensembles(i));
    EC = EC(1:360,1:360);
    ECpos = max(EC,0);
    ECpos1 = binarize_conn(ECpos);
    ECneg = abs(min(EC,0));
    ECneg1 = binarize_conn(ECneg);
    ECmod = ECpos1+ -1*ECneg1;
%     EC(vis,vis) = 0;
%     EC(sm,sm) = 0;
%     EC(da,da) = 0;
%     EC(va,va) = 0;
%     EC(lim,lim) = 0;
%     EC(fp,fp) = 0;
%     EC(dm,dm) = 0;
    
    inflow = sum(abs(ECmod),2);
    outflow = sum(abs(ECmod),1);
    
    pos_inflow = sum(max(ECmod,0),2);
    pos_outflow = sum(max(ECmod,0),1);
    neg_inflow = sum(-min(ECmod,0),2);
    neg_outflow = sum(-min(ECmod,0),1);
    
    infall = cat(2,infall,inflow);
    outfall = cat(1,outfall,outflow);
    
    pos_infall = cat(2,pos_infall,pos_inflow);
    neg_infall = cat(2,neg_infall,neg_inflow);
    pos_outfall = cat(1,pos_outfall,pos_outflow);
    neg_outfall = cat(1,neg_outfall,neg_outflow);
end

% investigate the information flows
figure; hist(outfall(:,50),50);

% set path
addpath(genpath('/local_raid1/01_software/toolboxes/toolboxes/cifti-matlab/')); %Add path to the gifti read
dir_atlas = '/local_raid1/03_user/younghyun/02_data/parcellations';

% Schaefer 100
aname = fullfile(dir_atlas, '/Schaefer2018_100Parcels_7Networks_order.dlabel.nii');
Schaefer = ft_read_cifti(aname).parcels;
[surf_lh,surf_rh] = load_conte69;


% calculate average
ave_inf = mean(infall,2);
ave_outf = mean(outfall,1);

ave_inf(:,2) = ave_outf';
contrastmap = zeros(size(ave_inf,1),1);
for i = 1:size(ave_inf,1)
    [temp1,ind1] = max(ave_inf(i,:)); temp2 = min(ave_inf(i,:));
%     if temp1-temp2 < 1
%         temp1 = 0;
%     end
    if ind1 == 1
        temp1 = -temp1;
    end
    contrastmap(i,1) = temp1;
end


% plot net inflow & outflow
ifandof = [ave_outf' ave_inf];
plot_hemispheres(ifandof,{surf_lh,surf_rh},...
                         'parcellation', atlas, ...
                         'labeltext', {'outflow','inflow'});
                     
% calculate average
ave_pos_inf = mean(pos_infall,2);
ave_pos_outf = mean(pos_outfall,1);
ave_neg_inf = mean(neg_infall,2);
ave_neg_outf = mean(neg_outfall,1);

ifandof = [ave_pos_outf' ave_pos_inf ave_neg_outf' ave_neg_inf];
plot_hemispheres(ifandof,{surf_lh,surf_rh},...
                         'parcellation', atlas, ...
                         'labeltext', {'pos outflow','pos inflow','neg outflow','neg inflow'});
                     
                     
                     
                     