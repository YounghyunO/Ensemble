%% Load Data
% Load toolbox
addpath('toolbox/surfstat');
addpath('toolbox/BCT/2019_03_03_BCT')

% Load surface gifti file
L_surf = gifti('data/S900.L.midthickness_MSMAll.10k_fs_LR.surf.gii');
R_surf = gifti('data/S900.R.midthickness_MSMAll.10k_fs_LR.surf.gii');
avg_surf.coord = [L_surf.vertices', R_surf.vertices'];
avg_surf.tri = [L_surf.faces; R_surf.faces+10242];

% Load shaefer parcellation
parcel = cifti_read('data/Schaefer2018_100Parcels_7Networks_order.dlabel.nii').cdata;

% Load effective connectiviy
EC = load('data/EC.mat').temp;

%% Analysis
% See: https://sites.google.com/site/bctnet/Home/functions


% 말씀하신 positive / negative outstrength (row sum)
pos_os = sum(max(EC, 0), 2);    % EC의 excitatory 부분만 row sum
neg_os = sum(min(EC, 0), 2);    % EC의 inhibitory 부분만 row sum

% BCT functions
% [id, od, deg] = degrees_dir(EC); % node degree
% [is, os, str] = strengths_dir(EC); % node strengths


%% Color List
% cyan to yellow
bs = 64; % block size
r = [zeros(1, bs * 2), linspace(0, 0.5, bs * 2), linspace(0.5, 0, bs * 2), zeros(1, bs * 2), linspace(0, 1, bs*5), ones(1, bs*3)];
g = [ones(1, bs), linspace(1, 0, bs*3), zeros(1, bs * 9), linspace(0, 1, bs * 3)];
b = [linspace(1, 0, bs), zeros(1, bs * 2), linspace(0, 2/3, bs), linspace(2/3, 0.5, bs), linspace(0.5, 2/3, bs), linspace(0.5, 0, bs*2), zeros(1, bs*8)];
mycol.shella = [r; g; b]';

mycol.blackblue =   flipud( [zeros(1,3)*0.8; zeros(127,1) ...
    (0:126)'/127 ones(127,1)]);
mycol.blue      =   flipud( [ones(1,3)*0.8; ...
    zeros(127,1) (0:126)'/127 ones(127,1)]);
mycol.red       =   flipud( [ones(1,3)*0.8; ...
    ones(500,1) linspace(0,253,500)'/254 zeros(500,1);...
    ones(64,1) ones(64,1) linspace(0,200,64)'/254]);

mycol.bluered = flipud( [[1, 0, 0]; ...
    [(126:-1:0)'/127 (0:126)'/127 zeros(127, 1)]; ...
    [zeros(127, 1) (126:-1:0)'/127 (0:126)'/127]; ...
    zeros(1,3)]);
yeo_colormap = [ 0 0 0;
    120 18 134;
    70 130 180;
    0 118 14;
    196 58 250;
    220 248 164;
    230 148 34;
    205 62 78 ]/255;

%% Visualization

% project parcellation to vertex

% feature for visualization
feature = pos_os;

ref = zeros(1, 20484);
for roi = 1:100
    ref(parcel == roi) = feature(roi);
end

SurfStatView(ref, avg_surf);
color = mycol.shella;   % 위의 colormap 중에 선택
SurfStatColormap( mycol.red  );