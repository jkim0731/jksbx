clear
mice = [36,37,39,52,53,54,56];
sessions = {[1,17],[7],[1,22],[3,21],[3],[3],[3]};  
rollingWindowForBaseF = 100; % in s
baseFprctile = 5;
baseDir = 'D:\TPM\JK\suite2p\';
for mi = 1 : length(mice)
    cd([baseDir, sprintf('%03d',mice(mi))])
    fList = dir('F_*_proc_final.mat');
    for i = 1 : length(fList)
        fprintf('mouse %03d session %d\n', mice(mi), i)
        load(fList(i).name)
        dat.rollingWindowForBaseF = rollingWindowForBaseF;
        dat.baseFprctile = baseFprctile;
        window = round(dat.rollingWindowForBaseF*(dat.ops.imageRate/dat.ops.num_plane));
        dF = zeros(size(dat.Fcell{1}), 'like', dat.Fcell{1});
        parfor j = 1 : size(dat.Fcell{1},1)
            if (dat.stat(j).iscell)
                fprintf('cell #%d/%d\n',j,size(dat.Fcell{1},1))
                tempF = dat.Fcell{1}(j,:) - dat.FcellNeu{1}(j,:) * dat.stat(j).neuropilCoefficient;
                window = round(dat.rollingWindowForBaseF*(dat.ops.imageRate/dat.ops.num_plane));
                tempFF = [tempF, tempF(end-window+1:end)];
                baseF = zeros(size(tempF));
                for k = 1 : length(baseF)
                    baseF(k) = prctile(tempFF(k:k+window),dat.baseFprctile);
                end
                dF(j,:) = (tempF-baseF)./baseF;
            end
        end
        dat.dF = dF;
        save(fList(i).name, 'dat')
    end
end