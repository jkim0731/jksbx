dBase = 'D:\TPM\JK\';
% fnList = {'025_1000', '027_1000', '030_1000', '030_2000', '036_997', '037_998', '038_998', '039_998', '041_998', '052_994', '052_995', '053_990', '054_995', '056_994'};
fnList = {'038_998', '041_998'};

for i = 1 : length(fnList)
    fnStr = strsplit(fnList{i}, '_');
    cd([dBase, fnStr{1}])
    make_zstack(fnList{i}, 1);
end


%%
dBase = 'D:\TPM\JK\';
% fnList = {'025_1000', '027_1000', '030_1000', '030_2000', '036_997', '037_998', '038_998', '039_998', '041_998', '052_994', '052_995','053_990', '054_995', '056_994'};
fnList = {'041_998'};

for i = 1 : length(fnList)
    fnStr = strsplit(fnList{i}, '_');
    cd([dBase, fnStr{1}])
    draw_dura(['zstack_', fnList{i}, '.mat'])
end

%%
% dBase = 'Y:\Whiskernas\JK\suite2p\';
% s2pList = {'025_004', '025_019', '025_022', '027_003', '027_016', '027_017', '030_003', '030_021', '030_022', '036_001', '036_017', '036_018','037_007', '038_002', '039_001', '039_022', '039_023','039_024','039_025', ...
%     '041_003', '052_003', '052_004', '052_021', '052_022', '053_003', '054_003', '054_004', '056_003', '056_004', '056_005'};
% zstackList = {'025_1000', '027_1000', '030_1000', '036_997', '037_998', '038_998', '039_998', '041_998', '052_994','053_990', '054_995', '056_994'};

dBase = 'D:\TPM\JK\suite2p\';
s2pList = {'027_010'};
zstackList = {'027_1000'};

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
        assign_depth(tempFlist(j).name, zstackFn(1).name, zstackDuraFn(1).name);
    end
end


%%
%% depth assigment (manual)
%%
% dBase = 'Y:\Whiskernas\JK\suite2p\';

% mice = [25,27,30,36,37,38,39,41,52,53,54,56];
% sessions = {[4,19,22],[3,16,17],[3,21,22],[1,17,18],[7],[2],[1,23,24],[3],[3,21,26],[3],[3],[3]};
% depths = {[110 : 165/3 : 110+165,  244 : 165/3 : 244+165], ... % 025
%     [110 : 165/3 : 110+165,  269 : 165/3 : 269+165], ... % 027
%     [98  : 165/3 : 98+165,   244 : 165/3 : 244+165], ... % 030
%     [98  : 165/3 : 98+165,   280 : 165/3 : 280+165], ... % 036    
%     [98  : 165/3 : 98+165,   299 : 165/3 : 299+165], ... % 037
%     [98  : 165/3 : 98+165,   269 : 165/3 : 269+165], ... % 038
%     [98  : 165/3 : 98+165,   299 : 165/3 : 299+165], ... % 039
%     [98  : 165/3 : 98+165,   299 : 165/3 : 299+165], ... % 041
%     [110 : 131/3 : 110+131,  269 : 131/3 : 269+131] ... % 052
%     [110 : 131/3 : 110+131,  269 : 131/3 : 269+131], ... % 053
%     [110 : 131/3 : 110+131,  269 : 131/3 : 269+131], ... % 054
%     [110 : 131/3 : 110+131,  269 : 131/3 : 269+131]};    % 056

dBase = 'D:\TPM\JK\suite2p\';
mice = 27;
sessions = {[10]};
depths = {[110 : 165/3 : 110+165,  269 : 165/3 : 269+165]};
for i = 1 : length(mice)
    mouse = mice(i);
    cd(sprintf('%s%03d',dBase,mice(i)))
    for si = 1 : length(sessions{i})
        session = sessions{i}(si);
        for j = 5 : 8
            fn = sprintf('F_%03d_%03d_plane%d_proc.mat',mouse, session,j);
    %         for k = 1 : length(flist)
                load(fn) % loading dat and spk

                dat.depth.imagingDepth = depths{i}(j);
                dat.depth.FOV = dat.depth.imagingDepth + round(dat.depth.depthComp);
                for kk = 1 : length(dat.stat)
                    % for readability, extend out actual x, y position of the
                    % cell, from the whole FOV of tpm image
                    yind = dat.ops.useY(dat.ops.yrange(round(dat.stat(kk).med(1))));
                    xind = dat.ops.useX(dat.ops.xrange(round(dat.stat(kk).med(2))));
                    dat.stat(kk).depth = dat.depth.FOV(yind, xind) * cosd(dat.depth.windowAngle);
                end
                save(fn, 'dat')
    %         end
        end
    end
end
