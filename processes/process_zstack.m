dBase = 'D:\TPM\JK\';
% fnList = {'025_1000', '027_1000', '030_1000', '030_2000', '036_997', '037_998', '038_998', '039_998', '041_998', '052_994', '052_995'};
fnList = {'053_990', '054_995', '056_994'};

for i = 1 : length(fnList)
    fnStr = strsplit(fnList{i}, '_');
    cd([dBase, fnStr{1}])
    make_zstack(fnList{i}, 1);
end


%%
dBase = 'D:\TPM\JK\';
fnList = {'025_1000', '027_1000', '030_1000', '030_2000', '036_997', '037_998', '038_998', '039_998', '041_998', '052_994', '052_995','053_990', '054_995', '056_994'};

for i = 1 : length(fnList)
    fnStr = strsplit(fnList{i}, '_');
    cd([dBase, fnStr{1}])
    draw_dura(['zstack_', fnList{i}, '.mat'])
end

%%
dBase = 'D:\TPM\JK\';
% s2pList = {'025_004', '025_019', '027_003', '027_016', '030_003', '030_021', '036_001', '036_017', '037_007', '038_002', '039_001', '039_022', '039_022', ...
%     '041_003', '052_003', '052_004', '052_021', '052_022', '053_003', '054_003', '054_004', '056_003', '056_004', '056_005'};
s2pList = {'052_003', '052_004', '052_021', '052_022', '053_003', '054_003', '054_004', '056_003', '056_004', '056_005'};
% zstackList = {'025_1000', '027_1000', '030_1000', '036_997', '037_998', '038_998', '039_998', '041_998', '052_994','053_990', '054_995', '056_994'};
zstackList = {'052_994','053_990', '054_995', '056_994'};
zstackMice = zeros(length(zstackList),1);
for i = 1 : length(zstackMice)
    temp = strsplit(zstackList{i},'_');
    zstackMice(i) = str2double(temp{1});
end

for i = 1 : length(s2pList)
    temp = strsplit(s2pList{i},'_');
    cd([dBase, temp{1}])
    s2pMouse = str2double(temp{1});
    zlistInd = find(zstackMice == s2pMouse, 1);
    zstackFn = dir(['zstack_', zstackList{zlistInd}, '*.mat']);    
    if length(zstackFn) > 1
        error('Too many zstack files')
    end
    zstackDuraFn = dir(['zstack_dura_', zstackList{zlistInd}, '*.mat']);
    if length(zstackDuraFn) > 1
        error('Too many zstackDura files')
    end
%     draw_dura(zstackFn(1).name)
%     draw_C2(zstackFn(1).name)
    tempFlist = dir(['F_', s2pList{i}, '_plane*_proc.mat']); % assuming everything is processed already
    for j = 1 : length(tempFlist)
        depth_and_C2(tempFlist(j).name, zstackFn(1).name, zstackDuraFn(1).name);
    end
end