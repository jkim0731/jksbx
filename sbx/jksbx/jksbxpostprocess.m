%% list of post processing

%% basic information
base_dir = 'D:\2p\JK';
mice = {'025','027','030'};
% sessions = 1:19;
% mice = {'650'};

%% split into trials
for mi = 1 : length(mice)
    cd([base_dir filesep mice{mi}])
    jksbxsplittrial_dir
end

%% alignment
for mi = 1 : length(mice)
    cd([base_dir filesep mice{mi}])
    jksbxaligndir_optotune_top    
end

%% align the whole imaging sessions
% for mi = 1 : length(mice)
%     cd([base_dir filesep mice{mi}])
%     flist = [];
%     for sn = sessions
%         flist = [flist;dir([sprintf('%s_%03d_',mice{mi},sn),'*.align'])];
%     end    
%     mg = cell(length(flist),1);
%     for i = 1 : length(flist)
%         load(flist(i).name,'-mat');
%         mg{i} = m;
%     end
%     [mm, TT] = jkaligncell(mg);
%     save([mice{mi}, '_align.mat'], 'mm', 'TT')
%     for i = 1 : length(flist)
%         load(flist(i).name,'-mat');
%         T = T + TT(i,:);
%         m = circshift(m,TT(i,:));
%         save(flist(i).name,'-mat','m','T')
%     end
% end


%% segment
% jksbxsegmenttool


