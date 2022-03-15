% set path
addpath(genpath('/local_raid1/01_software/toolboxes/toolboxes/cifti-matlab/')); %Add path to the gifti read
dir_atlas = '/local_raid1/03_user/younghyun/02_data/parcellations';

% MMP 370
aname = fullfile(dir_atlas, '/Q1-Q6_RelatedParcellation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
MMP = ft_read_cifti(aname).indexmax;

%Schaefer 100
aname = fullfile(dir_atlas,'/Schaefer2018_400Parcels_7Networks_order.dlabel.nii');
Schaefer = ft_read_cifti(aname).parcels;
dir_struct = '/local_raid1/03_user/younghyun/02_data/structure/100206/MNINonLinear/fsaverage_LR32k';

% Load Atlas
atlasroi = [gifti(fullfile(dir_struct, '100206.L.atlasroi.32k_fs_LR.shape.gii')).cdata; ...
            gifti(fullfile(dir_struct, '100206.R.atlasroi.32k_fs_LR.shape.gii')).cdata];


%%
close all

% load cortical surfaces 
subject_name = 'S900';
surf_type = 'inflated_MSMAll';
% surf_type = 'flat';
try
    [vertices, faces] = connRSMreadGII(dir_atlas, subject_name, surf_type); % get surfaces
catch
    error('Error reading surface. Perhaps gifti() function is not working. Try using gifti toolbox provided by HCP, https://wiki.humanconnectome.org/display/PublicData/HCP+Users+FAQ, Question 2 for help')
end

Cmap = zeros(64984,1);
for i =1:360
    ind = MMP == sorted(i);
%     seed = 27;
%     tp = max(sfc(:,seed,1),0);
%     inflow(seed) = 100;
    areas = [att; att];
%     temp = areas(i);
    temp = ave_inf(i);
    Cmap(ind,1) = temp;
end
Cmap(atlasroi == 0) = 0;
if strcmp(surf_type,'inflated_MSMAll')
    nsub = 7;
else
    nsub = 3;
end
f = figure;
% options for subtightplot
gapleft = 0.035;
gapright = 0.035;
margin_height = [0.03 0.03]; % lower, upper
margin_height_flat = [0.15 0.15]; % lower, upper
margin_width = [0.04 0.06]; % left, right


clrs = connMapVibModes2CmapRainbow_v2(min(Cmap),max(Cmap),360);
% clrs = lines(22);
clr_IDs = round((Cmap-min(Cmap))/max(Cmap-min(Cmap))*359)+1;
plot_colors = clrs(clr_IDs,:);


plot_data_mf_w_borders = plot_colors;
plot_data_mf_w_borders = plot_HCP_boundaries(plot_colors,'all','black');

 if strcmp(surf_type,'flat')
            subtightplot(1,nsub,1,[gapleft gapright], margin_height_flat, margin_width);
            options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
            plot_mesh(vertices.left, faces.left, options);
            view([5 90])
            material('dull')
            
            subtightplot(1,nsub,2,[gapleft gapright], margin_height_flat, margin_width);
            options.face_vertex_color = plot_data_mf_w_borders(length(plot_data_mf_w_borders)/2+1:end,:);
            plot_mesh(vertices.right, faces.right, options);
            view([-10 90])
            material('dull')
elseif strcmp(surf_type,'inflated_MSMAll')

    % left lateral
    subtightplot(1,nsub,5,[gapleft gapright], margin_height, [0.04 0.07]);
    options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
    plot_mesh(vertices.left, faces.left, options);
    view([270 0])
    material('dull')
    camlight(190,-180)

    % right lateral
    subtightplot(1,nsub,2,[gapleft gapright], [0.2 0.2], [0.04 0.02]);
    options.face_vertex_color = plot_data_mf_w_borders(length(plot_data_mf_w_borders)/2+1:end,:);
    plot_mesh(vertices.right, faces.right, options);
    view([90 0])
    material('dull')
    camlight(145,215)

    % both from the top
    subtightplot(1,nsub,3,[gapleft gapright], [0.1 0.1], [0.02 0]);
    options.face_vertex_color = plot_data_mf_w_borders;
    plot_mesh(vertices.all, faces.all, options);
    view([0 90])
    material('dull')
    camlight(0,270)

    % both from below
    subtightplot(1,nsub,4,[gapleft gapright], [0.1 0.1], [0 0.04]);
    options.face_vertex_color = plot_data_mf_w_borders;
    plot_mesh(vertices.all, faces.all, options);
    view([180 -90])
    material('dull')
    camlight(180,180)


    % left medial
    subtightplot(1,nsub,1,[gapleft gapright], margin_height, margin_width);
    options.face_vertex_color = plot_data_mf_w_borders(1:length(plot_data_mf_w_borders)/2,:);
    plot_mesh(vertices.left, faces.left, options);
    view([90 0])
    material('dull')
    camlight(-190,240)
    
    % right medial
    subtightplot(1,nsub,6,[gapleft gapright], margin_height, margin_width);
    options.face_vertex_color = plot_data_mf_w_borders(length(plot_data_mf_w_borders)/2+1:end,:);
    plot_mesh(vertices.right, faces.right, options);
    view([270 0])
    material('dull')
    camlight(-190,-180)


 end
 
 if min(Cmap)<0
    subtightplot(1,nsub,7,[0.04 0.04], [0.025 0], [0.05 0.01]);
    colormap(clrs)
    axis off
    h=colorbar();
    set(h,'XTick',[0,abs(min(Cmap))/(abs(min(Cmap))+abs(max(Cmap))),1])
    % set labels
    xlabs = [min(Cmap),0,max(Cmap)];
    xlabs = round(xlabs);
    set(h,'XTickLabel',xlabs)
    set(h, 'Position', [0.8542    0.2373    0.0179    0.5763]);
 else
    subtightplot(1,nsub,7,[0.04 0.04], [0.025 0], [0.05 0.01]);
    colormap(clrs)
    axis off
    h=colorbar();
    set(h,'XTick',[0,1])
    xlabs = [0,max(Cmap)];
    set(h,'XTickLabel',xlabs)
    set(h, 'Position', [0.8542    0.2373    0.0179    0.5763]);
 end
 
 if strcmp(surf_type,'flat')
        set(h,'AxisLocation','out')
    else
        set(h,'AxisLocation','in')
end   
 set(gcf,'Color','white')
 f.Units='centimeters';    

    if strcmp(surf_type,'inflated_MSMAll')
        f.Position = [26.4319   19.0235   31.8823    6.8263];
    else
        f.Position = [26.4319   20.2406   20.3465    6.8263];
    end