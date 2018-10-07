clear
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  
rollingWindowForBaseF = 100; % in s
baseFprctile = 5;
baseDir = 'C:\Data\suite2p\';
for mi = 1 : length(mice)
    cd([baseDir, sprintf('%03d',mice(mi))])
    fList = dir('F_*_proc_final.mat');
    for i = 1 : length(fList)
        load(fList(i).name)
        dat.rollingWindowForBaseF = rollingWindowForBaseF;
        dat.baseFprctile = baseFprctile;
        dat.dF = zeros(size(dat.Fcell{1}), 'like', dat.Fcell{1});
        for j = 1 : size(dat.Fcell{1},1)
            tempF = dat.Fcell{1}(j,:) - dat.FcellNeu{1}(j,:) * dat.stat(j).neuropilCoefficient;
            window = round(dat.rollingWindowForBaseF*(dat.ops.imageRate/dat.ops.num_plane));
            tempFF = [tempF, tempF(end-window+1:end)];
            baseF = zeros(size(tempF));
            for k = 1 : length(baseF)
                baseF(k) = prctile(tempFF(k:k+window),dat.baseFprctile);
            end
            dat.dF(j,:) = (tempF-baseF)./baseF;
        end
        save(fList(i).name, 'dat')
    end
end