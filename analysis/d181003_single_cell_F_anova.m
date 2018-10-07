% # of trials, # of touches in each angle in 2 different volumes
% Considering just unimodal (or monotonic) tuning for now 2018/06/27 JK

% load u and, if it exists, ANOVA tuning file
mice = [25,27,30,36,37,39,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[1,22],[3,21],[3],[3],[3]};  

showImages = 0;

% for mi = 1 : length(mice)
for mi = 1
%     for si = 1 : length(sessions{mi})
    for si = 1

        u = Uber.buildUberArray(mice(mi), sessions{mi}(si));

        savefn = [u.mouseName,u.sessionName,'singleCell_F_ANOVA.mat'];
        if exist(savefn, 'file')
            load(savefn)
        else

        % settings
        angles = 45:15:135;
        baseFrameNum = 3;
        afterFrameNum = ceil(u.frameRate);
        thresholdAnova = 0.05;
        thresholdTtestSharpness = 0.05;
        thresholdTtestResponse = 0.05;

        excludeDrinkingTime = 1;
        onlyBeforeDecision = 1;
        if onlyBeforeDecision 
            onlyAfterDecision = 0;
        else
            onlyAfterDecision = 1;
        end
        allowOverlap = 0;

        cellsTuned = [];
        tuneAngle = []; % Max response angle, in case of multimodal response
        tuneDirection = []; % 1 up, 2 down, 3 bipolar
        tuneBipolarity = []; % Ratio between |max| and |min|. Calculated only when tuneDirection is 3(bipolar)
        tuneAmplitude = []; % response to the tuned angle.
        tuneModulationMaxmin = []; % dF/F0 (%p). Always positive. (max - min)
        tuneSharpness = []; % based on t-test to neighbors
        tuneSharpnessFWHM = []; % # of neighbors with response larger than half of the modulation. For now it will be useful for only within-angle comparison
        tuneSharpnessTtest = []; % # of neighbors with same tune direction. For now it will be useful for only within-angle comparison
        tuneResponseProb = []; % proportion of response to the tuned angle
        tuneReliabilty = []; % To the max response angle. Average of Pearson's correlation to the average time series
        
        tuneSingle = [];
        tuneBroad = []; 
        tuneLOO = []; % leave-one-out
        tuneMM = []; % multimodal. Including bipolar.
        tuneCateg = [] ; % categorical (>= 90 or <= 90)
        tuneRamp = []; % ramping up or down

        cellsNTResponse = [];
        NTRdirection = []; % 1 up, 2 down
        NTRresponseProb = []; % proportion of response to the touch
        NTRreliability = []; % average of Pearson's correlation to the average time series (from all the angles)
        cellsNTNR = [];
        
        % for confirmation
        dFtotal = cell(length(u.cellNums),1); % for all cells. ordered by cell number
        cellsTotal = u.cellNums;
        
        

        % exit = 0;
        for cellid = 1:length(u.cellNums)
        % figure, set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.2, 0.3, 0.7, 0.5])
        % ci = 1;
        % while true
        %     cellid = ci;
            disp(['Processing cell id ', num2str(cellid)])
            cellNum = u.cellNums(cellid);
            touchNumTrial = cell(length(angles),1);
            touchNumChunk = cell(length(angles),1);
            countedTouchs = cell(length(angles),1);
            dF = cell(length(angles),1);
            F = cell(length(angles),1);
            for ai = 1 : length(angles)
                angle = angles(ai);
                plane = floor(cellNum/1000);
                trialPlaneInd = find(cellfun(@(x) ~isempty(find(x.planes == plane, 1)), u.trials));
                trialAngleInd = find(cellfun(@(x) x.angle == angle, u.trials));
                trialInd = intersect(trialPlaneInd, trialAngleInd);

                touchNum = 0;
                touchNumTrial{ai} = [];
                touchNumChunk{ai} = [];
                for i = 1 : length(trialInd)
                    if ~isempty(u.trials{trialInd(i)}.touchChunks)                
                        if excludeDrinkingTime && ~isempty(u.trials{trialInd(i)}.drinkingTime)
                            tempChunksDrinking = u.trials{trialInd(i)}.touchChunks;
                            chunkIndDrinking = 1:length(tempChunksDrinking);
                            drinkInd = [];
                            for j = 1 : length(tempChunksDrinking)
                                if tempChunksDrinking{j}(end) > u.trials{trialInd(i)}.drinkingTime(1) && tempChunksDrinking{j}(1) < u.trials{trialInd(i)}.drinkingTime(2)
                                    drinkInd = [drinkInd, j];
                                end
                            end
                            chunkIndDrinking = setdiff(chunkIndDrinking, drinkInd);
                            touchChunks = cell(length(chunkIndDrinking),1);
                            for j = 1 : length(touchChunks)
                                touchChunks{j} = tempChunksDrinking{chunkIndDrinking(j)};                        
                            end
                            chunkIndSave = chunkIndDrinking;
                        else
                            touchChunks = u.trials{trialInd(i)}.touchChunks;
                            chunkIndSave = 1 : length(u.trials{trialInd(i)}.touchChunks);
                        end

                        if onlyBeforeDecision && ~isempty(u.trials{trialInd(i)}.answerLickTime)
                            decisionChunkInd = length(touchChunks)+1;
                            for j = 1 : length(touchChunks)
                                if touchChunks{j}(1) > u.trials{trialInd(i)}.answerLickTime
                                    decisionChunkInd = j;
                                    break
                                end
                            end
                            if decisionChunkInd == 1
                                touchChunks = {};
                                chunkIndSave = [];
                            else
                                decisionInd = 1 : decisionChunkInd - 1;
                                tempTouchChunks = touchChunks;
                                touchChunks = cell(length(decisionInd),1);
                                for j = 1 : length(decisionInd)
                                    touchChunks{j} = tempTouchChunks{decisionInd(j)};
                                end
                                chunkIndSave = chunkIndSave(decisionInd);
                            end
                        end

                        if onlyAfterDecision && ~isempty(u.trials{trialInd(i)}.answerLickTime)
                            decisionChunkInd = length(touchChunks)+1;
                            for j = 1 : length(touchChunks)
                                if touchChunks{j}(1) > u.trials{trialInd(i)}.answerLickTime
                                    decisionChunkInd = j;
                                    break
                                end
                            end
                            if decisionChunkInd == length(touchChunks)
                                touchChunks = {};
                                chunkIndSave = [];
                            else
                                decisionInd = decisionChunkInd : length(touchChunks);
                                tempTouchChunks = touchChunks;
                                touchChunks = cell(length(decisionInd),1);
                                for j = 1 : length(decisionInd)
                                    touchChunks{j} = tempTouchChunks{decisionInd(j)};
                                end
                                chunkIndSave = chunkIndSave(decisionInd);
                            end
                        end

                        if ~isempty(touchChunks)
                            touchNum = touchNum + length(touchChunks);
                            touchNumTrial{ai} = [touchNumTrial{ai}; ones(length(touchChunks),1) * trialInd(i)];
                            touchNumChunk{ai} = [touchNumChunk{ai}; chunkIndSave'];
                        end
                    end
                end
                %%
                countedTouchs{ai} = 0;
                dF{ai} = nan(touchNum, baseFrameNum + afterFrameNum);
                F{ai} = nan(touchNum, baseFrameNum + afterFrameNum);
                for i = 1 : touchNum    
                    trial = touchNumTrial{ai}(i);
                    chunk = touchNumChunk{ai}(i);
                    time = u.trials{trial}.touchChunks{chunk}(1); % for now, just consider the first touch point
                    frameInd = find(u.trials{trial}.tpmTime{mod(plane-1,4)+1} >= time, 1, 'first');
                    if i > 1 && ~allowOverlap
                        if trial == trialOld
                            if frameInd - frameOld < baseFrameNum
                                continue
                            end
                        end
                    end
                    countedTouchs{ai} = countedTouchs{ai} + 1;
                    Find = find(u.trials{trial}.neuindSession == cellNum);
                    if frameInd < baseFrameNum
                        tempBaseNum = frameInd;
                        F{ai}(i,baseFrameNum-frameInd+1:baseFrameNum) = u.trials{trial}.F(Find, 1:frameInd);
                    else
                        F{ai}(i,1:baseFrameNum) = u.trials{trial}.F(Find, frameInd-baseFrameNum + 1 : frameInd);
                    end

                    maxafternum = max(size(u.trials{trial}.F,2)-frameInd, afterFrameNum);

                    if size(u.trials{trial}.F,2)-frameInd < afterFrameNum
                        tempAfter = size(u.trials{trial}.F,2)-frameInd;
                        F{ai}( i, baseFrameNum + 1 : baseFrameNum + tempAfter ) = u.trials{trial}.F(Find, frameInd + 1 : frameInd + tempAfter);
                    else
                        F{ai}( i, baseFrameNum + 1 : end) = u.trials{trial}.F(Find, frameInd + 1 : frameInd + afterFrameNum);
                    end

                    tempBase = nanmean(F{ai}(i, 1 : baseFrameNum));
                    dF{ai}(i,:) = (F{ai}(i,:) - tempBase)/tempBase * 100;
                    frameOld = frameInd;
                    trialOld = trial;
                end
                dF{ai}(isnan(mean(dF{ai},2)),:) = []; % removing all the NaN rows
            end
            dFtotal{cellid} = dF;
            
            % ANOVA
            touchNumGroups = cellfun(@(x) size(x,1), dF);
            timeAveragedF = zeros(sum(touchNumGroups),1);
            anovaGroupsF = zeros(sum(cellfun(@(x) size(x,1), dF)),1);
            for ai = 1 : length(angles)
                timeAveragedF( sum(touchNumGroups(1:(ai-1))) + 1 : sum(touchNumGroups(1:ai)) ) = mean(dF{ai}(:,baseFrameNum+1:end),2);
                anovaGroupsF( sum(touchNumGroups(1:(ai-1))) + 1 : sum(touchNumGroups(1:ai)) ) = deal(ai);
            end

            [anovaP, anovaTable, anovaStat] = anova1(timeAveragedF, anovaGroupsF, 'off');
            [pairComp, meanNsem] = multcompare(anovaStat, 'Display', 'off');
            [~, ind] = min(abs(anovaStat.means));
            meanSign = sign(anovaStat.means);
            if anovaP <= thresholdAnova        
                tTestH = zeros(length(angles),1);
                for ai = 1 : length(angles)
                    tTestH(ai) = ttest(timeAveragedF( sum(touchNumGroups(1:(ai-1))) + 1 : sum(touchNumGroups(1:ai)) ));
                end
                tTestInd = find(tTestH);
                if isempty(tTestInd) % happens sometimes, when the response level is low
                    cellsNTNR = [cellsNTNR; cellNum];
                else
                    cellsTuned = [cellsTuned; cellNum];
                    [maxAmp, maxInd] = max(abs(anovaStat.means(tTestInd)));
                    tunedAngleInd = tTestInd(maxInd); % for readability
                    tuneAngle = [tuneAngle; angles(tunedAngleInd)];
                    tuneAmplitude = [tuneAmplitude; maxAmp];
                    tuneModulationMaxmin = [tuneModulationMaxmin; max(anovaStat.means) - min(anovaStat.means)];
                    if min(anovaStat.means(tTestInd)) * max(anovaStat.means(tTestInd)) < 0
                        tuneDirection = [tuneDirection; 3]; % bipolar                
                        if anovaStat.means(tunedAngleInd) > 0 % Response to major tuned angle is positive
                            tuneBipolarity = [tuneBipolarity; abs(min(anovaStat.means(tTestInd))) / abs(max(anovaStat.means(tTestInd)))];
                        else
                            tuneBipolarity = [tuneBipolarity; abs(max(anovaStat.means(tTestInd))) / abs(min(anovaStat.means(tTestInd)))];
                        end
                    elseif min(anovaStat.means(tTestInd)) > 0
                        tuneDirection = [tuneDirection; 1]; % increase
                    elseif max(anovaStat.means(tTestInd)) < 0
                        tuneDirection = [tuneDirection; 2]; % decrease
                    else
                        tuneDirection = [tuneDirection; 0]; % error
                        tuneModulationMaxmin = [tuneModulationMaxmin(1:end-1); 0]; % error
                    end

                    % Sharpness
                    % Using pairwise comparison p-value. 
                    % (Failed because even though anova p < 0.05, all the pair-wise mean differences can be insignificant.??)
                    lowerThanThis = find(pairComp(:,2) == tunedAngleInd);
                    tempNumNeighbors = 1; % itself
                    if ~isempty(lowerThanThis)
                        for li = 1 : length(lowerThanThis)                            
                            if pairComp(lowerThanThis(end-li+1),6) <= thresholdTtestSharpness % going backwards in "lowerThanThis"
                                break
                            else
                                tempNumNeighbors = tempNumNeighbors + 1;
                            end
                        end
                    end
                    higherThanThis = find(pairComp(:,1) == tunedAngleInd);
                    if ~isempty(higherThanThis)
                        for hi = 1 : length(higherThanThis)
                            if pairComp(higherThanThis(hi),6) <= thresholdTtestSharpness
                                break
                            else
                                tempNumNeighbors = tempNumNeighbors + 1;
                            end
                        end
                    end        
                    tuneSharpness = [tuneSharpness; tempNumNeighbors];

                    % Using FWHM
                    tempNumNeighbors = 1; % include itself
                    if anovaStat.means(tunedAngleInd) > 0 % Positive tuning
                        FWHM = anovaStat.means(tunedAngleInd) - tuneModulationMaxmin(end)/2;
                        for li = tunedAngleInd-1 : -1 : 1
                            if anovaStat.means(li) > FWHM
                                tempNumNeighbors = tempNumNeighbors + 1;
                            else
                                break
                            end
                        end
                        for hi = tunedAngleInd+1 : length(angles)
                            if anovaStat.means(hi) > FWHM
                                tempNumNeighbors = tempNumNeighbors + 1;
                            else
                                break
                            end
                        end
                    else % Negative tuning
                        FWHM = anovaStat.means(tunedAngleInd) + tuneModulationMaxmin(end)/2;
                        for li = tunedAngleInd-1 : -1 : 1
                            if anovaStat.means(li) < FWHM
                                tempNumNeighbors = tempNumNeighbors + 1;
                            else
                                break
                            end
                        end
                        for hi = tunedAngleInd+1 : length(angles)
                            if anovaStat.means(hi) < FWHM
                                tempNumNeighbors = tempNumNeighbors + 1;
                            else
                                break
                            end
                        end
                    end
                    tuneSharpnessFWHM = [tuneSharpnessFWHM; tempNumNeighbors];

                    % tune sharpness using neighbors of significant response in the same direction
                    % How many bins are tuned?? (Not testing whether tuning strength is the same 
                    tempNumNeighbors = 1; % include itself                    
                    if anovaStat.means(tunedAngleInd) > 0 % Positive tuning
                        posTuning = find(anovaStat.means > 0);                        
                    else % Negative tuning
                        posTuning = find(anovaStat.means < 0);
                    end
                    simTuning = posTuning .* tTestH;
                    for li = tunedAngleInd-1:-1:1
                        if simTuning(li)
                            tempNumNeighbors = tempNumNeighbors + 1;
                        else
                            break
                        end
                    end
                    for hi = tunedAngleInd+1:length(tTestH)
                        if simTuning(hi)
                            tempNumNeighbors = tempNumNeighbors + 1;
                        else
                            break
                        end
                    end
                    tuneSharpnessTtest = [tuneSharpnessTtest; tempNumNeighbors];                    
                    
                    % Response probability
                    dFtuned = dF{tunedAngleInd}(:,baseFrameNum+1 : end); % for readability
                    responseNum = 0;
                    for ri = 1:size(dFtuned,1)
                        tempTS = dFtuned(ri,:);
                        [~, p] = ttest(tempTS');
                        if p <= thresholdTtestResponse
                            if anovaStat.means(tunedAngleInd) * mean(tempTS) > 0
                                responseNum = responseNum + 1;
                            end
                        end                        
                    end
                    tuneResponseProb = [tuneResponseProb; responseNum/ri*100];

                    % Response reliability
                    template = mean(dFtuned);
                    rho = zeros(size(dFtuned,1),1);
                    for ri = 1 : length(rho)
                        tempTS = dFtuned(ri,:);
                        rho(ri) = corr(template',tempTS');
                    end
                    tuneReliabilty = [tuneReliabilty; mean(rho)];
                    
                    % Categorization
                    if sum(tTestH) == 1
                        tuneSingle = [tuneSingle; cellNum];
                    elseif sum(tTestH) > 1
                        tuneBroad = [tuneBroad; cellNum];
                        if sum(tTestH) == length(tTestH) - 1
                            tuneLOO = [tuneLOO; cellNum]; % leave-one-out
                        elseif sum(tTestH) == length(tTestH) && abs(length(find(meanSign > 0)) - length(find(meanSign < 0))) == length(tTestH) - 2
                            tuneLOO = [tuneLOO; cellNum]; % leave-one-out
                        elseif find(tTestH>0,1) < length(tTestH) 
                            if sum(abs(diff(tTestH(find(tTestH>0,1) : end)))) > 1
                                tuneMM = [tuneMM; cellNum]; % multimodal. Including bipolar.
                            end
                        elseif sum(abs( diff( tTestH(1 : (end-1)/2) ) )) == 0 && sum(abs( diff( tTestH( (end+1)/2 : end ) ) )) == 0 % all values < 90 are the same in t-test, AND all values > 90 are the same                             
                            if tTestH(1) == 1  && tTestH(end) == 0 % only first half is responding
                                if sum(abs( diff( meanSign(1 : (end-1)/2) ) )) == 0 % all the first half direction is the same
                                    tuneCateg = [tuneCateg; cellNum] ; % categorical (>= 90 or <= 90)
                                end
                            elseif tTestH(1) == 0  && tTestH(end) == 1 % only second half is responding
                                if sum(abs( diff( meanSign( (end+1)/2 : end ) ) )) == 0 % all the second half direction is the same
                                    tuneCateg = [tuneCateg; cellNum] ; % categorical (>= 90 or <= 90)
                                end
                            elseif tTestH(1) == 1 && tTestH(end) == 1 % both first half and second half is responding
                                if sum(abs( diff( meanSign(1 : (end-1)/2) ) )) == 0 && ... % all the first half direction is the same
                                    sum(abs( diff( meanSign( (end+1)/2 : end ) ) )) == 0 && ... % all the second half direction is the same
                                    meanSign(1) ~= meanSign(end) % direction between first half and second half is different
                                    tuneCateg = [tuneCateg; cellNum] ; % categorical (>= 90 or <= 90)
                                end
                            end 
                        elseif isempty(find(diff(sign(diff(anovaStat.means))),1)) % everything is going up or down 
                            tuneRamp = [tuneRamp; cellNum]; % ramping up or down
                        end
                    end
                end
            else
                responseH = ttest(timeAveragedF);
                if responseH
                    cellsNTResponse = [cellsNTResponse; cellNum];
                    if mean(timeAveragedF) > 0
                        NTRdirection = [NTRdirection; 1]; % increase
                    elseif mean(timeAveragedF) < 0
                        NTRdirection = [NTRdirection; 2]; % decrease
                    else
                        NTRdirection = [NTRdirection; 0]; % error
                    end

                    % Response probability
                    dFall = zeros(sum(cellfun(@(x) size(x,1), dF)), afterFrameNum); % for readability
                    for ai = 1 : length(angles)
                        dFall( sum(touchNumGroups(1:(ai-1))) + 1 : sum(touchNumGroups(1:ai)), : ) = dF{ai}(:,baseFrameNum+1 : end);
                    end
                    responseNum = 0;
                    for ri = 1 : size(dFall,1)
                        tempTS = dFall(ri,:);
                        [~, p] = ttest(tempTS');
                        if p <= thresholdTtestResponse
                            if mean(timeAveragedF) * mean(tempTS) > 0
                                responseNum = responseNum + 1;
                            end
                        end
                    end
                    NTRresponseProb = [NTRresponseProb; responseNum/ri*100];            

                    % Response reliability
                    template = mean(dFall);
                    rho = zeros(size(dFall,1),1);
                    for ri = 1 : size(dFall,1)
                        tempTS = dFall(ri,:);                
                        rho(ri) = corr(template', tempTS');
                    end
                    NTRreliability = [NTRreliability; mean(rho)];
                else
                    cellsNTNR = [cellsNTNR; cellNum];
                end
            end
        end

        save(savefn, 'cells*','tune*','NTR*', 'dFtotal')
        end
    end
end