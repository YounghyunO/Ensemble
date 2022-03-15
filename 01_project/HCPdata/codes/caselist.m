

for j = 1:6
    switch j 
        case 1; task = 'emotion'; case 2; task = 'gambling'; case 3; task = 'relational'; ...
        case 4; task = 'wm'; case 5; task = 'social'; case 6; task = 'language'; 
    end
    list = load(sprintf('%s_caselist.txt',task));
    temp = [];
    for i = 1:size(list,1)

        temp2 = list(i);
        path = sprintf('/local_raid1/03_user/younghyun/02_data/task/%s/%d/MNINonLinear/Results',task,temp2);
        all_files = dir(path);
        num_dir = numel(all_files);
        
        if num_dir == 4
            temp = cat(1,temp,temp2);
        end
    end
    writematrix(temp,sprintf('%s_caselist.txt',task));
end

load('emotion_caselist.txt')
load('gambling_caselist.txt')
load('relational_caselist.txt')
load('wm_caselist.txt')
load('social_caselist.txt')
load('language_caselist.txt')


list = [];
for i = 1:1072
    temp = gambling_caselist(i,1);
    if sum(temp == emotion_caselist(:,1)) && sum(temp == relational_caselist(:,1)) ...
            && sum(temp == wm_caselist(:,1)) && sum(temp == social_caselist(:,1)) ...
            && sum(temp == language_caselist(:,1))
        list = [list; temp];
    end
end

% list = load('caselist.txt');
% idx = list ==668361;
idx = list == 146533; % movement doesn't exist
idx = list == 196952; % MSMA doesn't exist
list(idx) = [];

writematrix(temp,'caselist.txt')


% currently not using as the caselist is already computed
% Add caselist
listDir = fullfile('/local_raid1/03_user/younghyun/01_project/HCPdata/run');
%%
list = load(fullfile(listDir,'/caselist_discovery.txt'));
idx = list == 227432; 
list(idx) = [];
idx = list == 833148;
list(idx) = [];
idx = list == 987074;
list(idx) = [];
idx = list == 749361;
list(idx) = [];
idx = list == 837560;
list(idx) = [];
idx = list == 992673;
list(idx) = [];
idx = list == 989987;
list(idx) = [];
writematrx(list,'caselist_discovery.txt')

parfor i = 1:size(list,1)
    uplist(i,1) = censorsub(list(i))
end
uplist = nonzeros(uplist);
writematrix(uplist,'replication_caselist_wo_outlier.txt')


% for discovery caselist
% for replication caselist
% idx = list == 135225;
% list(idx) = [];
% idx = list == 620434;
% list(idx) =[];
% idx = list == 668361;
% list(idx) = [];
% idx = list == 751348;
% list(idx) = [];
% idx = list == 737960;
% list(idx) = [];
% idx = list == 987983;
% list(idx) = [];
% idx = list == 978578;
% list(idx) = [];
% idx = list == 993675;
% list(idx) = [];
% idx = list == 990366;
% list(idx) = [];
