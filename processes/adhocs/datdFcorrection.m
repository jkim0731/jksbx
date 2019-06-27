clear
% baseDir = 'Y:\Whiskernas\JK\suite2p\';
% mice = [25,27,30,36,37,38,39,41,52,53,54,56, 70,74,75,76];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3], [6],[4],[4],[4]};  
baseDir = 'D:\TPM\JK\suite2p\';
mice = 27;
sessions = {[9,10]};



rollingWindowForBaseF = 100; % in s
baseFprctile = 5;


for mi = 1 : length(mice)
% for mi = 1
    mouse = mice(mi);
    cd([baseDir, sprintf('%03d',mouse)])
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);
    
        fList = dir(sprintf('F_%03d_%03d_plane*_proc_final.mat', mouse, session));
        for i = 1 : length(fList)
%         for i = 1
            fprintf('mouse %03d session %d plane %d\n', mouse, session, i)
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
end

%% Split just for now. Integrate into the upper for loops later (2018/10/15)
%% Error estimation. Using 5th percentile of std of smoothed n (5?) frames window.
% Gather all, then merge, and take unique frames. Lower bound = 5 % of
% total frames + (n-1) frames. Upper bound = 5 % of total frames X n.
% -> No. This cannot average out slow drift. When decided to subtract the mean
% rather than using it's values. Do not care redundancy. 
clear
% baseDir = 'Y:\Whiskernas\JK\suite2p\';
% mice = [25,27,30,36,37,38,39,41,52,53,54,56, 70,74,75,76];
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3], [6],[4],[4],[4]};  
baseDir = 'D:\TPM\JK\suite2p\';
mice = 27;
sessions = {[9,10]};

stdwindow = 5;
lowerprct = 5; % 5th percentile
iteration = 10000;


for mi = 1 : length(mice)
    mouse = mice(mi);
    cd([baseDir, sprintf('%03d',mouse)])
    
    for si = 1 : length(sessions{mi})
        session = sessions{mi}(si);

        fList = dir(sprintf('F_%03d_%03d_plane*_proc_final.mat', mouse, session));
        for i = 1 : length(fList)
%         for i = 1
            fprintf('mouse %03d session %d plane %d\n', mouse, session, i)
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
            dat.noise = cellfun(@(x) sqrt(sum(diffs.^2)/length(diffs)), errorDiffs);
            dat.errorResponse = errorResponse;
            dat.errorDiffs = errorDiffs;
            dat.stdwindow = stdwindow;
            dat.lowerprct = lowerprct;
            save(fList(i).name, 'dat')
        end
    end
end
