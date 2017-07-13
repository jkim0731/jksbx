%% Alignment
mice = {'648','649'};
base_dir = 'D:\2p\JK';
%

for i = 1 : size(mice,2)    
    cd([base_dir, filesep, mice{i}])
    sbxaligndir
end

%% Split the image into each trials
mice = {'650', '653'};
base_dir = 'D:\2p\JK';
for i = 1 : length(mice)
    cd([base_dir, filesep, mice{i}])   
    jksbxsplittrial_dir
end
