clear
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  
rollingWindowForBaseF = 100; % in s
baseFprctile = 5;
baseDir = 'D:\TPM\JK\suite2p\';
for mi = 1 : length(mice)
% for mi = 1
    cd([baseDir, sprintf('%03d',mice(mi))])
    fList = dir('F_*_proc_final.mat');
    for i = 1 : length(fList)
%     for i = 1
        fprintf('mouse %03d session %d\n', mice(mi), i)
        load(fList(i).name)
        dat.rollingWindowForBaseF = rollingWindowForBaseF;
        dat.baseFprctile = baseFprctile;
        
        dF = zeros(size(dat.Fcell{1}), 'like', dat.Fcell{1});
        inds = find([dat.stat.iscell]);
        npcoeffs = min(min(dat.Fcell{1}(inds,:) ./ dat.FcellNeu{1}(inds,:), [], 2), ones(length(inds),1)*0.7);
                            
        len = size(dat.Fcell{1},2);
        window = round(rollingWindowForBaseF*(dat.ops.imageRate/dat.ops.num_plane));

        tempF = dat.Fcell{1}(inds,:) - dat.FcellNeu{1}(inds,:) .* npcoeffs;
        tempFF = [tempF, tempF(:,end-window+1:end)];
        ttbase = zeros(length(inds),size(tempF,2));
        parfor k = 1 : len
            ttbase(:,k) = prctile(tempFF(:,k:k+window),baseFprctile,2);
        end        
        dat.dF = (tempF - ttbase)./ttbase;
        dat.npcoeffs = npcoeffs;
        save(fList(i).name, 'dat')
    end
end

%% Split just for now. Integrate into the upper for loops later (2018/10/15)
%% Error estimation. Using 5th percentile of std of smoothed n (5?) frames window.
% Gather all, then merge, and take unique frames. Lower bound = 5 % of
% total frames + (n-1) frames. Upper bound = 5 % of total frames X n.
% -> No. This cannot average out slow drift. When decided to subtract the mean
% rather than using it's values. Do not care redundancy. 
clear
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  
stdwindow = 5;
lowerprct = 5; % 5th percentile
iteration = 10000;
baseDir = 'D:\TPM\JK\suite2p\';
for mi = 1 : length(mice)
    cd([baseDir, sprintf('%03d',mice(mi))])
    fList = dir('F_*_proc_final.mat');
    for i = 1 : length(fList)
        fprintf('mouse %03d session %d\n', mice(mi), i)
        load(fList(i).name) % loading dat
        frameNum = round(dat.ops.imageRate/dat.ops.num_plane);
        errorDiffs = cell(size(dat.dF,1),1);
        errorResponse = zeros(size(dat.dF,1),2);
        for j = 1 : size(dat.dF,1)
            tempstd = zeros(size(dat.dF,2) - (stdwindow-1), 1);
            temp = dat.dF(j,:);
            tempsmooth = smooth(temp,5);
            parfor k = 1 : length(tempstd)
                tempstd(k) = std(tempsmooth(k+2:k+(stdwindow-1)/2+2)); % chopping off flanking half-windows
            end
            tempstdInds = find(tempstd<=prctile(tempstd,lowerprct));
            baseInds = zeros(length(tempstdInds),stdwindow);
            diffs = zeros(length(tempstdInds),stdwindow);
            parfor k = 1 : length(tempstdInds)
                baseInds(k,:) = tempstdInds(k) : tempstdInds(k)+(stdwindow-1);
                diffs(k,:) = temp(baseInds(k,:)) - mean(temp(baseInds(k,:)));
            end
%             errorInds = unique(baseInds(:)); % collect all frames that were within threshold from smoothed std.
%             % making chunks before diff
%             chunkPoints = [1, find(diff(errorInds')>1) + 1, length(errorInds')+1]; % +1 because of diff. first 1 for the first chunk. So this is actually start points of each chunk.
%             chunks = cell(1,length(chunkPoints)-1); %
%             diffs = [];
%             for ci = 1 : length(chunks)
%                 chunks{ci} = errorInds(chunkPoints(ci) : chunkPoints(ci+1)-1);
%                 diffs = [diffs, temp(chunks{ci}) - mean(chunks{ci})]; % maybe this is better representation of errors.
%                 % diffs = [diffs, diff(temp(chunks{ci}))]; % diverges a lot. 5th really low and 95th really high. Again, maybe because the constraint is gone. 
%             end
            
            diffs = diffs(:);
            errorDiffs{j} = diffs;
            errorResponseDist = zeros(iteration,1); 
            parfor ii = 1 : iteration
                inds = randperm(length(diffs),frameNum);
                errorResponseDist(ii) = mean(cumsum(diffs(inds)));                
            end
            errorResponse(j,1) = prctile(errorResponseDist,5);
            errorResponse(j,2) = prctile(errorResponseDist,95);
        end
        dat.noise = cellfun(@(x) sqrt(sum(x.^2)/length(x)), errorDiffs);
        dat.errorResponse = errorResponse;
        dat.errorDiffs = errorDiffs;
        dat.stdwindow = stdwindow;
        dat.lowerprct = lowerprct;
        save(fList(i).name, 'dat')
    end
end
