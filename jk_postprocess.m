%% Split the image into each trials
mice = {'036', '037', '039','041'};
base_dir = 'D:\2p\JK';
for i = 1 : length(mice)
    cd([base_dir, filesep, mice{i}])   
    jksbxsplittrial_dir
end


%% Alignment
mice = {'036','037','039','041'};
base_dir = 'D:\2p\JK';

for i = 1 : size(mice,2)    
    cd([base_dir, filesep, mice{i}])
    jksbxaligndir({},'no')
end

