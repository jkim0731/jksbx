%% Alignment
mice = {'045','046','047','048'};
base_dir = 'D:\2p\JK';

for i = 1 : size(mice,2)    
    cd([base_dir, filesep, mice{i}])
    jksbxaligndir({},'no')
end

%% Split the image into each trials
mice = {'025', '027', '030'};
base_dir = 'D:\TPM\JK';
for i = 1 : length(mice)
    cd([base_dir, filesep, mice{i}])   
    jksbxsplittrial_dir
end
