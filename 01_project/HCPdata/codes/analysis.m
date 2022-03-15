%% rDCM
temp1 = rDCMs(~cellfun(@isempty,rDCMs));

rDCMtemp = [];
for i = 1:size(temp1,2)
    rDCMtemp = cat(3,rDCMtemp,cell2mat(temp1(i)));
end

EC = mean(rDCMtemp,3);
EC_rDCM = EC_rDCM(1:50,1:50);
std_rDCM = std(rDCMtemp,1,3);


figure; imagesc(EC_rDCM); title('Group rDCM Result'); colorbar
saveas(gca,'group_rdcm_result.png')
figure; imagesc(std_rDCM); title('standard deviation'); colorbar
saveas(gca,'standard_deviation_rdcm.png')
save('EC_rDCM.mat','EC_rDCM')

%% FASK
temp2 = FASKs(~cellfun(@isempty,FASKs));

FASKtemp = [];
for i = 1:size(temp2,2)
    FASKtemp = cat(3,FASKtemp,cell2mat(temp2(i)));
end

EC_FASK = mean(FASKtemp,3);
EC_FASK = EC_FASK(1:50,1:50);
std_FASK = std(FASKtemp,1,3);

figure; imagesc(EC_FASK); title('Group FASK Result (Right hemi)'); colorbar
saveas(gca,'group_fask_result.png')
figure; imagesc(std_FASK); title('standard deviation FASK'); colorbar
saveas(gca,'standard_deviation_fask.png')
save('EC_FASK.mat','EC_FASK')

%% VAR
temp3 = VARs(~cellfun(@isempty,VARs));


VARtemp = [];
for i = 1:size(temp3,2)
    VARtemp = cat(3,VARtemp,cell2mat(temp3(i)));
end

EC_VAR = mean(VARtemp,3);
EC_VAR = EC_VAR(1:50,1:50);
% std_VAR = std(VARtemp,1,3);

figure; imagesc(EC_VAR); title('Group VAR Result (Left hemi)'); colorbar
saveas(gca,'group_var_result.png')
figure; imagesc(std_VAR); title('VAR standard deviation'); colorbar
saveas(gca,'standard_deviation_var.png')

save('EC_VAR.mat','EC_VAR')

%% GC
temp4 = GCs(~cellfun(@isempty,GCs));

GCtemp = [];
for i = 1:size(temp4,2)
    GCtemp = cat(3,GCtemp,cell2mat(temp4(i)));
end

EC_GC = mean(GCtemp,3);
EC_GC = EC_GC(1:50,1:50);
std_GC = std(GCtemp,1,3);

figure; imagesc(EC_GC); title('Group GC Result (Right hemi)'); colorbar
saveas(gca,'group_gc_result.png')
figure; imagesc(std_GC); title('GC standard deviation'); colorbar
saveas(gca,'standard_deviation_gc.png')

save('EC_GC.mat','EC_GC')

%% Ensemble
temp5 = Ensembles(~cellfun(@isempty,Ensembles));

Ensembletemp = [];
for i = 1:size(temp5,2)
    Ensembletemp = cat(3,Ensembletemp,cell2mat(temp5(i)));
end

EC = mean(Ensembletemp,3);


EC = zeros(400,400);
EC_std = zeros(400,400);
EC_var_coeff = zeros(400,400);
for i = 1:400
    for j = 1:400
        if i ~= j
            temp = squeeze(Ensembletemp(i,j,:));
%             std_temp = std(nonzeros(temp));
            temp(abs(temp)<std(temp))=0;
            EC(i,j) = mean(nonzeros(temp));
            EC_std(i,j) = std(nonzeros(temp));
            temp2 = EC_std(i,j)^2/EC(i,j);
            if abs(temp2) > 10
                temp2 =0;
            end
            EC_var_coeff(i,j) = temp2;
        end
    end
end
EC(isnan(EC)) = 0;
 [~,BOLD,~] = Network_sim_210110(EC/4);
FC = corr(BOLD); FC = FC-diag(diag(FC));
figure; imagesc(FC)

figure; imagesc(EC); title('Group Ensemble Result'); colorbar
saveas(gca,'group_average_ensemble.png')
figure; imagesc(std_En); title('Ensemble standarnd deviation'); colorbar
saveas(gca,'standard_deviation_ensemble.png')

save('EC_En.mat','EC_En')

%% FC
temp6 = FCs(~cellfun(@isempty,FCs));

FC_temp = [];
for i = 1:size(temp6,2)
    FC_temp = cat(3,FC_temp,cell2mat(temp6(i)));
end

FC = mean(FC_temp,3);
FC = FC(51:100,51:100);
std_FC = std(FC_temp,1,3);


figure; imagesc(FC); title('Group FC Result'); colorbar
saveas(gca,'group_average_fc.png')
figure; imagesc(std_FC); title('FC standarnd deviation'); colorbar
saveas(gca,'fc_standard_deviation.png')

save('FC.mat','FC')




%% Extract only cotrext
EC_En = EC_En(1:100,1:100);
EC_En = EC_En - diag(diag(EC_En));

EC_rDCM = EC_rDCM(1:100,1:100);
EC_rDCM = EC_rDCM - diag(diag(EC_rDCM));

EC_GC = EC_GC(1:100,1:100);
for i = 1:100
    EC_GC(i,i) = 0;
end

EC_VAR = EC_VAR(1:100,1:100);
EC_VAR = EC_VAR - diag(diag(EC_VAR));

FC = FC(1:360,1:360);

%% First principal gradient of all algirhtms
addpath(genpath('/local_raid1/01_software/toolboxes/toolboxes/cifti-matlab/')); %Add path to the gifti read

% gradient analysis using principal component analysis
dir_atlas = '/local_raid1/03_user/younghyun/02_data/parcellations';

%Yeo 7 atlas
aname = fullfile(dir_atlas, 'Yeo2011_7Networks_N1000.dlabel.nii');
Yeo7 = ft_read_cifti(aname).parcels;
aname = fullfile(dir_atlas, 'Yeo2011_17Networks_N1000.dlabel.nii');
Yeo17 = ft_read_cifti(aname).parcels;
% Yeo_parcel = ft_read_cifti(aname2).split_components;
% Glasser et al., 2016 atlas
aname = fullfile(dir_atlas, '/Schaefer2018_100Parcels_7Networks_order.dlabel.nii');
Schaefer = ft_read_cifti(aname).parcels;
[surf_lh,surf_rh] = load_conte69;
gm = GradientMaps('approach','pca');


% gradient calculation
gm_en = gm.fit(EC_En);
gm_rdcm = gm.fit(EC_rDCM);
gm_gc = gm.fit(EC_GC);
gm_var = gm.fit(EC_VAR);
gm_fc = gm.fit(FC);

% first gradients
gradients = [gm_en.gradients{1}(:,1) gm_rdcm.gradients{1}(:,1) ...
                     gm_gc.gradients{1}(:,1) gm_var.gradients{1}(:,1) gm_fc.gradients{1}(:,1)];

%plot
 
plot_hemispheres(Yeo_parcel,{surf_lh,surf_rh},...
                             'parcellation', labeling, ...
                             'labeltext', {'Yeo 7'});
saveas(gca,'first principal gradients.png')

% three principal gradients of Ensemble 
gradients = gm_en.gradients{1}(:,1:3);

%plot
plot_hemispheres(gradients,{surf_lh,surf_rh},...
                             'parcellation', labeling, ...
                             'labeltext', {'gradient1','gradient2','gradient3'});
saveas(gca,'Ensemble_three_gradients.png')


%% three principal gradients of ensemble

% first show gradient plot of Ensemble
scree_plot(gm_en.lambda{1}); title('Ensemble')
saveas(gca,'Ensemble eigenvalues.png')

gradient_in_euclidean(gm_en.gradients{1}(:,1:3),{surf_lh,surf_rh},labeling); title('Ensemble')
saveas(gca,'Ensemble gradient distance plot.png')

% now show gradient plot of FC to elucidate the difference
scree_plot(gm_fc.lambda{1}); title('FC')
saveas(gca,'FC eigenvalues.png')

gradient_in_euclidean(gm_fc.gradients{1}(:,1:3),{surf_lh, surf_rh},labeling); title('FC')
saveas(gca,'FC gradient distance plot.png')

%% inflow and outflow 

% optional
% extract column-wise top 10 percent
EC_pos = max(EC_En, 0);
EC_pos_10p = EC_pos;
EC_pos_10p(EC_pos < prctile(EC_pos,95)) = 0; 

EC_neg = min(EC_En, 0);
EC_neg_10p = EC_neg;
EC_neg_10p(EC_neg>prctile(EC_neg,5)) = 0;

EC_En_10p = EC_pos_10p + EC_neg_10p;
 
EC_out = EC_En_10p; % outflow -- column wise 10percent

EC_pos = max(EC_En, 0);
EC_pos_10p = EC_pos;
EC_pos_10p(EC_pos < prctile(EC_pos,95,2)) = 0; 

EC_neg = min(EC_En, 0);
EC_neg_10p = EC_neg;
EC_neg_10p(EC_neg>prctile(EC_neg,5,2)) = 0;

EC_En_10p = EC_pos_10p + EC_neg_10p;

EC_in = EC_En_10p; % inflow -- row wise 10 percent

% otherwise

EC = EC_En;

target = 2;
target_of = EC(:,target);  
target_of(target,1) = 1;
target_if = EC(target,:);

ifandof = [target_of target_if'];
plot_hemispheres(ifandof,{surf_lh,surf_rh},...
                         'parcellation', labeling2, ...
                         'labeltext', {'outflow','inflow'});

% inflow and outflow of rDCM, GC & Ensemble
en_pos_of = rescale(sum(max(EC, 0), 1));  
en_neg_of = rescale(-sum(min(EC, 0), 1));    
en_pos_if = rescale(sum(max(EC, 0), 2));
en_neg_if = rescale(-sum(min(EC, 0), 2));

% inflow and outflow of rDCM, GC & Ensemble
en_pos_of = sum(max(EC, 0), 1);  
en_neg_of = -sum(min(EC, 0), 1);    
en_pos_if = sum(max(EC, 0), 2);
en_neg_if = -sum(min(EC, 0), 2);



ifandof = [en_pos_of' en_pos_if en_neg_of' en_neg_if];

plot_hemispheres(ifandof,{surf_lh,surf_rh},...
                             'parcellation', Schaefer, ...
                             'labeltext', {'pos of','pos if','neg of','neg if'});
% colormap(mymap)
saveas(gca,'positive&negative flows.png')

% net inflow and outflow
en_of = rescale(sum(abs(EC),1));
en_if = rescale(sum(abs(EC),2));
 
ifandof  = [en_of' en_if];

plot_hemispheres(ifandof,{surf_lh,surf_rh},...
                             'parcellation', Schaefer, ...
                             'labeltext', {'net of','net if'});
saveas(gca,'NET positive&negative flows.png')
 
%% target inflow and outflow
close all
EC = Ensembletemp(1:360,1:360,1);
EC_pos = max(EC, 0);
EC_pos_10p = EC_pos;
EC_pos_10p(EC_pos < prctile(EC_pos,90)) = 0; 

EC_neg = min(EC, 0);
EC_neg_10p = EC_neg;
EC_neg_10p(EC_neg>prctile(EC_neg,10)) = 0;

% inflow and outflow of rDCM, GC & Ensemble
en_pos_of = sum(EC_pos_10p, 1);  
en_neg_of = -sum(EC_neg_10p, 1);    
en_pos_if = sum(EC_pos_10p, 2);
en_neg_if = -sum(EC_neg_10p, 2);



%[1:9,51:58] = visual
%[10:15,59:66] = somato-motor
%[16:23,67:73] = Dorsal attention
%[24:30,74:78] = Ventral attention
%[31:33,79:80] = Limbic
%[34:37,81:89] = Froto parietal
%[38:50,90:100] = Default
target_of = sum(EC(:,target),2);
target_if = sum(EC(target,:),1)';
% target_of(target,1) = 1;
% target_if(1,target) = 1;

% plot_hemispheres([EC(:,24),EC(:,25),EC(:,26),EC(:,27)],{surf_lh,surf_rh},...
%                              'parcellation', Schaefer);
plot_hemispheres([target_of target_if],{surf_lh,surf_rh},...
                             'parcellation', Schaefer);
 
target = 24:30;    
figure('Position',[100 100 600 800]); title('seed based flow plot')
subplot(1,4,1)
imagesc(sum(max(EC(:,target),0),2)); title('positive outflow')
subplot(1,4,2)
imagesc(sum(-min(EC(:,target),0),2)); title('negative outflow')
subplot(1,4,3)
imagesc(sum(max(EC(target,:),0),1)'); title('positive inflow')
subplot(1,4,4)
 imagesc(sum(-min(EC(target,:),0),1)'); title('negative inflow')
 

%%




scree_plot(gm_gc.lambda{1}); title('GC')
saveas(gca,'GC eigenvalues.png')
scree_plot(gm_fask.lambda{1}); title('FASK')
saveas(gca,'FASK eigenvalues.png')
scree_plot(gm_rdcm.lambda{1}); title('rDCM')
saveas(gca,'rDCM eigenvalues.png')
scree_plot(gm_var.lambda{1}); title('VAR')
saveas(gca,'VAR eigenvalues.png')


gradient_in_euclidean(gm_rdcm.gradients{1}(:,1:2),{surf_lh,surf_rh},labeling); title('rDCM')
saveas(gca,'rDCM gradient distance plot.png')
gradient_in_euclidean(gm_fask.gradients{1}(:,1:2),{surf_lh,surf_rh},labeling); title('FASK')
saveas(gca,'FASK gradient distance plot.png')
gradient_in_euclidean(gm_var.gradients{1}(:,1:2),{surf_lh,surf_rh},labeling); title('VAR')
saveas(gca,'VAR gradient distance plot.png')
gradient_in_euclidean(gm_gc.gradients{1}(:,1:2),{surf_lh,surf_rh},labeling); title('GC')
saveas(gca,'GC gradient distance plot.png')

%% network analysis

%modularity

module_rdcm = modularity_dir(EC_rDCM);
module_fask = modularity_dir(EC_FASK);
module_var = modularity_dir(EC_VAR);
module_gc = modularity_dir(EC_GC);
module_ec = modularity_dir(EC_En);

modules = [module_gc module_fask ...
                   module_rdcm module_var module_ec];

plot_hemispheres(modules,{surf_lh,surf_rh},...
                             'parcellation', Schaefer, ...
                             'labeltext', {'GC','FASK','rDCM','VAR','EC'});

%Louvain community

community_rdcm = community_louvain(EC_rDCM,[],[],'negative_asym');
community_fask = community_louvain(EC_FASK);
community_var = community_louvain(EC_VAR,[],[],'negative_asym');
community_gc = community_louvain(EC_GC);
community_ec = community_louvain(EC_En,[],[],'negative_asym');

communities = [community_gc community_fask ...
                   community_rdcm community_var community_ec];

plot_hemispheres(communities,{surf_lh,surf_rh},...
                             'parcellation', labeling, ...
                             'labeltext', {'GC','FASK','rDCM','VAR','EC'});
 saveas(gca,'louvain community detection.png')
 
 %% plot in Yeo 7network
 
 Yeo_atlas = cell(7,1);

for i = 1:100
    idx = Schaefer == i;
    temp = Yeo7(idx,1);
    temp2 = mode(temp);
    Yeo_atlas{temp2}(end+1) = i;
end

data = [];
for i = 1:28
    if i<8
        temp1 = en_pos_of(1,cell2mat(Yeo_atlas(i)))';
        data = cat(1,data,temp1);
    elseif i>7 && i<15
        temp1 = en_neg_of(1,cell2mat(Yeo_atlas(i-7)))';
        data = cat(1,data,temp1);
    elseif i>14 && i<22
        temp1 = en_pos_if(cell2mat(Yeo_atlas(i-14)),1);
        data = cat(1,data,temp1);
    elseif i>21
        temp1 = en_neg_if(cell2mat(Yeo_atlas(i-21)),1);
        data = cat(1,data,temp1);
    end
end

idx = [];
for i = 1:7
    idx(end+1) =  length(cell2mat(Yeo_atlas(i)));
end
Yeo_idx = cumsum(idx);

origin = cell(100,1);
for i = 1:100
    if i < Yeo_idx(1)+1
        origin(i) = {'Visual'};
    elseif i < Yeo_idx(2)+1
        origin(i) = {'Somatomotor'};
    elseif i < Yeo_idx(3)+1
            origin(i) = {'Dorsal Attention'};
    elseif i < Yeo_idx(4)+1
            origin(i) = {'Ventral Attention'};
    elseif i < Yeo_idx(5)+1
            origin(i) = {'Limbic'};
    elseif i < Yeo_idx(6)+1
            origin(i) = {'Frontoparietal'};
    else
            origin(i) = {'Default'};
    end
end
origincopy = [origin; origin; origin; origin];

color = cell(400,1);
for i = 1:400
    if i<101
        color(i) = {'positive outflow'};
    elseif i<201
        color(i) = {'negative outflow'};
    elseif i<301
        color(i) = {'positive inflow'};
    else
        color(i) = {'negative inflow'};
    end
end

g = gramm('x',color,'y',data,'color',origincopy);
g.stat_boxplot('dodge',0.8)
figure('Position',[100 100 800 550])
g.draw()


% PCA analysis
eigs = zeros(100,1);
for i = 1: 100 % Schaefer 100 parcellation

    x = ts(atlas == i,:)'; % ts is z-scored time series (64984*1195) observation*variable
    [coeff,score,latent,~,lambda,mu] = pca(x);
    eigs(i) = lambda(1);
    
end

