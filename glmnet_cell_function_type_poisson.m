% Extracting input matrices for GLM analysis in each neuron from an Uber_2padArray u
% 
% Dependency: 
%     - Uber class
%     - jkWhiskerOnsetNAmplitude
% 
% 
% inputs: 
%     - mouse (as in number)
%     - session (as in number) 
%     - cellnum (1~length(total number of cells)),
%     - nShift (number of frames to shift, either forward or backward. Default: 3)
% 
%
% outputs: 
%     - cid: cell id (1000~8999)
%     - frameRate
%     - spk: spikes (vector. Padded with NaN's of length nShift before and after each trial)
%     
%     % sensory variables: shift backward only
%     % Same length as spk. 
%     - tTouchCount: # of total touches within the frame (parameter)
%     - pTouchCount: # of protraction touches within the frame (parameter)
%     - rTouchCount: # of retraction touches within the frame (parameter)
%     
%     - tTouchFrames: total touch frames (binary)
%     - pTouchFrames: protraction touch frames (binary)
%     - rTouchFrames: retraction touch frames (binary)
%
%     - tTouchDuration: total touch duration within each tpm frame (parameter, in ms)
%     - pTouchDuration: protraction touch duration within each tpm frame (parameter, in ms)
%     - rTouchDuration: retraction touch duration within each tpm frame (parameter, in ms)     

%       Up to here, have one as all angles and add each angles (total number of predictors: 1 + length(angles)
%
%     - scPiezo: piezo sound cue onset (binary)
%     - scPoleup: pole up sound cue onset (binary)
%     - scPoledown: pole down sound cue onset (binary)
%     
%
%     - drinkOnset: drinking onset (binary)
%     
%     
%     % motor variables: shift both backward and forward
%     - whiskingOnset: whisking onset (parameter; # of onset in each frame)
%     - whiskingAmp: whisking amplitude (parameter; from whisker decomposition; maximum of the frame)
%     - whiskingOA: whisking onset & amplitude. maximum amplitude where there was whisking onset (>= 1)
%     - whiskingMidpoint: whisking midpoint(parameter; from whisker decomposition)
%
%
%     - bLick: total licks within the frame (parameter; # of licks in each frame)
%     - lLick: left licks within the frame (parameter)
%     - rLick: right licks within the frame (parameter)
%     
%     - bLickOnset: the frame where lick onset happened (binary; each bout is calculated as 1 s interval)
%     - bLickOffset: the frame where lick offset happened (binary; each bout is calculated as 1 s interval)
%     - lLickOnset
%     - lLickOffset
%     - rLickOnset
%     - rLickoffset
% 
%     - firstLick: the frame where the first lick of the trial happened (binary)
%     - lastLick: the frame where the last lick of the trial happened (binary)
%     - firstLeftLick
%     - lastLeftLick
%     - firstRightLick
%     - lastRightLick

baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 

% myCluster = parcluster('local');
% delete(myCluster.Jobs)
% poolobj = gcp('nocreate'); % If no pool, do not create new one.
% if isempty(poolobj)    
%     parpool(34, 'SpmdEnabled', false);
% elseif poolobj.NumWorkers < poolobj.Cluster.NumWorkers - 2
%     nw = poolobj.Cluster.NumWorkers-2;
%     delete(poolobj)
%     parpool(nw, 'SpmdEnabled', false);
% end
%%
% mice = [36,37,38,39,41,52,53,54,56];
% sessions = {[17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 


for mi = 1:length(mice)
    % for mi = 1
    for si = 1:length(sessions{mi})
%     for si = 2

        mouse = mice(mi);
        session = sessions{mi}(si);

        posShift = 4;
        negShift = 2;
        testPortion = 0.3; % 30 % test set
        pThresholdNull = 0.05;
        pThresholdPartial = 0.05;
        lickBoutInterval = 1; % licks separated by 1 s regarded as different licking bouts
        numShuffle = 1000; % number of shuffling for testing r for ridge regression

        glmnetOpt = glmnetSet;
        glmnetOpt.standardize = 0;
        glmnetOpt.alpha = 0;
        
        partialGlmOpt = glmnetOpt;
        partialGlmOpt.alpha = 0;
        % lambdaCV = 3; % cross-validation fold number

        dn = sprintf('%s%03d',baseDir,mouse);
        ufn = sprintf('UberJK%03dS%02d.mat', mouse, session);
        cd(dn)
        if exist(ufn,'file')
            load(ufn)
        else
            u = Uber.buildUberArray(mouse, session);
        end
        frameRate = u.frameRate;

        savefnResult = sprintf('glmResponseType_JK%03dS%02d_glmnet_m11',mouse, session); % m(n) meaining method(n)
%         savefnPredictors = sprintf('glmResponseType_JK%03dS%02d_predictors_m7',mouse, session); % just for the input. Run once and use them for test

%         if exist(savefnPredictors, 'file')
%             load(savefnPredictors)
%         else
            %% pre-processing for lick onset and offset
            % regardless of licking alternating, each l and r has it's own lick onset and offset. both licking, just take the union

            v = cell(length(u.trials),1);

            for ui = 1 : length(u.trials)
                bothLickTime = union(u.trials{ui}.leftLickTime, u.trials{ui}.rightLickTime);
                if length(bothLickTime) == 1
                    v{ui}.bothLickOnset = bothLickTime;
                    v{ui}.bothLickOffset = bothLickTime;
                    v{ui}.firstLick = bothLickTime;
                    if u.trials{ui}.response < 1
                        v{ui}.lastLick = bothLickTime;
                    else
                        v{ui}.lastLick = [];
                    end
                elseif length(bothLickTime) > 1
                    onsets = find(diff(bothLickTime) > lickBoutInterval);
                    if isempty(onsets)
                        v{ui}.bothLickOnset = bothLickTime(1);
                        v{ui}.bothLickOffset = bothLickTime(end);
                    else
                        v{ui}.bothLickOnset = bothLickTime([1; onsets+1]);
                        v{ui}.bothLickOffset = bothLickTime([onsets; end]);
                    end
                    v{ui}.firstLick = bothLickTime(1);
                    if u.trials{ui}.response < 1
                        v{ui}.lastLick = bothLickTime(end);
                    else
                        v{ui}.lastLick = [];
                    end
                else
                    v{ui}.bothLickOnset = [];
                    v{ui}.bothLickOffset = [];
                    v{ui}.firstLick = [];
                    v{ui}.lastLick = [];
                end


                if length(u.trials{ui}.leftLickTime) == 1
                    v{ui}.leftLickOnset = u.trials{ui}.leftLickTime;
                    v{ui}.leftLickOffset = u.trials{ui}.leftLickTime;
                    v{ui}.firstLeftLick = u.trials{ui}.leftLickTime;
                    if u.trials{ui}.response < 1 
                        if isempty(u.trials{ui}.rightLickTime)
                            v{ui}.lastLeftLick = u.trials{ui}.leftLickTime;
                        elseif u.trials{ui}.rightLickTime(end) < u.trials{ui}.leftLickTime(end)
                            v{ui}.lastLeftLick = u.trials{ui}.leftLickTime;
                        else
                            v{ui}.lastLeftLick = [];
                        end
                    else
                        v{ui}.lastLeftLick = [];
                    end
                elseif length(u.trials{ui}.leftLickTime) > 1        
                    onsets = find(diff(u.trials{ui}.leftLickTime) > lickBoutInterval);
                    if isempty(onsets)
                        v{ui}.leftLickOnset = u.trials{ui}.leftLickTime(1);
                        v{ui}.leftLickOffset = u.trials{ui}.leftLickTime(end);                    
                    else
                        v{ui}.leftLickOnset = u.trials{ui}.leftLickTime([1; onsets+1]);
                        v{ui}.leftLickOffset = u.trials{ui}.leftLickTime([onsets; end]);
                    end
                    v{ui}.firstLeftLick = u.trials{ui}.leftLickTime(1);
                    if u.trials{ui}.response < 1 
                        if isempty(u.trials{ui}.rightLickTime)
                            v{ui}.lastLeftLick = u.trials{ui}.leftLickTime;
                        elseif u.trials{ui}.rightLickTime(end) < u.trials{ui}.leftLickTime(end)
                            v{ui}.lastLeftLick = u.trials{ui}.leftLickTime;
                        else
                            v{ui}.lastLeftLick = [];
                        end
                    else
                        v{ui}.lastLeftLick = [];
                    end
                else
                    v{ui}.leftLickOnset = [];
                    v{ui}.leftLickOffset = [];
                    v{ui}.firstLeftLick = [];
                    v{ui}.lastLeftLick = [];
                end


                if length(u.trials{ui}.rightLickTime) == 1
                    v{ui}.rightLickOnset = u.trials{ui}.rightLickTime;
                    v{ui}.rightLickOffset = u.trials{ui}.rightLickTime;
                    v{ui}.firstRightLick = u.trials{ui}.rightLickTime;
                    if u.trials{ui}.response < 1 
                        if isempty(u.trials{ui}.leftLickTime)
                            v{ui}.lastRightLick = u.trials{ui}.rightLickTime;
                        elseif u.trials{ui}.leftLickTime(end) < u.trials{ui}.rightLickTime(end)
                            v{ui}.lastRightLick = u.trials{ui}.rightLickTime;
                        else
                            v{ui}.lastRightLick = [];
                        end
                    else
                        v{ui}.lastRightLick = [];
                    end
                elseif length(u.trials{ui}.rightLickTime) > 1        
                    onsets = find(diff(u.trials{ui}.rightLickTime) > lickBoutInterval);
                    if isempty(onsets)
                        v{ui}.rightLickOnset = u.trials{ui}.rightLickTime(1);
                        v{ui}.rightLickOffset = u.trials{ui}.rightLickTime(end);
                    else
                        v{ui}.rightLickOnset = u.trials{ui}.rightLickTime([1; onsets+1]);
                        v{ui}.rightLickOffset = u.trials{ui}.rightLickTime([onsets; end]);
                    end
                    v{ui}.firstRightLick = u.trials{ui}.rightLickTime(1);
                    if u.trials{ui}.response < 1 
                        if isempty(u.trials{ui}.leftLickTime)
                            v{ui}.lastRightLick = u.trials{ui}.rightLickTime;
                        elseif u.trials{ui}.leftLickTime(end) < u.trials{ui}.rightLickTime(end)
                            v{ui}.lastRightLick = u.trials{ui}.rightLickTime;
                        else
                            v{ui}.lastRightLick = [];
                        end
                    else
                        v{ui}.lastRightLick = [];
                    end
                else
                    v{ui}.rightLickOnset = [];
                    v{ui}.rightLickOffset = [];
                    v{ui}.firstRightLick = [];
                    v{ui}.lastRightLick = [];
                end
                v{ui}.tpmTime = u.trials{ui}.tpmTime;
            end


%             %% repetition test
%     %         division = 20;
%             repetition = 1;
%             rtest = struct;
%             for ri = 1 : repetition % repetition index
                %% divide into training set and test set (70%, 30%)
                % based on the animal touched or not, the choice (same as the result since I'm going to mix the pole angles, so right, wrong, and miss), pole angles (2 or 7), and the distance (if there were multiple distances)
                % in this order, make trees, and take 30% of the leaves (or equivalently, take all the possible intersections and take 30%)

                angles = unique(cellfun(@(x) x.angle, u.trials));
                distances = unique(cellfun(@(x) x.distance, u.trials));

                touchGroup = cell(2,1);
                choiceGroup = cell(3,1);
                angleGroup = cell(length(angles),1);
                distanceGroup = cell(length(distances),1);

                ptouchGroup{1} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) length(x.protractionTouchChunks), u.trials))));
                ptouchGroup{2} = setdiff(u.trialNums, ptouchGroup{1});

                rtouchGroup{1} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) length(x.retractionTouchChunks), u.trials))));
                rtouchGroup{2} = setdiff(u.trialNums, rtouchGroup{1});

                choiceGroup{1} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.response == 1, u.trials))));
                choiceGroup{2} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.response == 0, u.trials))));
                choiceGroup{3} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.response == -1, u.trials))));

                for i = 1 : length(angles)
                    angleGroup{i} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.angle == angles(i), u.trials))));
                end

                for i = 1 : length(distances)
                    distanceGroup{i} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.distance == distances(i), u.trials))));
                end
                %%
                testTn = [];
                for pti = 1 : length(ptouchGroup)
                    for rti = 1 : length(rtouchGroup)
                        for ci = 1 : length(choiceGroup)
                            for ai = 1 : length(angleGroup)
                                for di = 1 : length(distanceGroup)
                                    tempTn = intersect(ptouchGroup{pti}, intersect(rtouchGroup{rti}, intersect(choiceGroup{ci}, intersect(angleGroup{ai}, distanceGroup{di}))));
                                    if ~isempty(tempTn)
                                        tempTn = tempTn(randperm(length(tempTn)));
                                        testTn = [testTn; tempTn(1:round(length(tempTn)*0.3))];
                                    end
                                end
                            end
                        end
                    end
                end
                %
                [~,testInd] = ismember(testTn, u.trialNums);

                trainingTn = setdiff(u.trialNums, testTn);
                [~,trainingInd] = ismember(trainingTn, u.trialNums);
                %% design matrix for training lambda
                trainingInputMat = cell(8,1);
                for cgi = 1:2 % cell group index
                % for cgi = 1
                    tindcell = find(cellfun(@(x) ismember(1001+(cgi-1)*4000, x.neuindSession), u.trials));

                    tind = intersect(tindcell, trainingInd);
                    for plane = 1 : 4    
                %     for plane = 1
                        pTouchCount = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(cellfun(@(y) y(1), x.protractionTouchChunks), [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
%                         rTouchCount = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(cellfun(@(y) y(1), x.retractionTouchChunks), [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
%                         tTouchCount = pTouchCount + rTouchCount;
%                         pTouchCount(pTouchCount>0) = deal(1);
% 
                        whiskerVideoFrameDuration = u.trials{tind(1)}.frameDuration * 1000; % in ms
% 
%                         pTouchDuration = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.protractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime{plane}]), nan(1,posShift)], ...
%                             u.trials(tind)','uniformoutput',false));
%                         rTouchDuration = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.retractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime{plane}]), nan(1,posShift)], ...
%                             u.trials(tind)','uniformoutput',false));
%                         tTouchDuration = pTouchDuration + rTouchDuration;
% 
%                         pTouchFrames = pTouchDuration;
%                         pTouchFrames(pTouchDuration > 0) = 1;
%                         rTouchFrames = rTouchDuration;
%                         rTouchFrames(rTouchDuration > 0) = 1;
%                         tTouchFrames = pTouchFrames + rTouchFrames;

                        pTouchCountAngles = cell(length(angles)+1,1);
%                         rTouchCountAngles = cell(length(angles)+1,1);
%                         tTouchCountAngles = cell(length(angles)+1,1);
%                         pTouchDurationAngles = cell(length(angles)+1,1);
%                         rTouchDurationAngles = cell(length(angles)+1,1);
%                         tTouchDurationAngles = cell(length(angles)+1,1);
%                         pTouchFramesAngles = cell(length(angles)+1,1);
%                         rTouchFramesAngles = cell(length(angles)+1,1);
%                         tTouchFramesAngles = cell(length(angles)+1,1);
                        for ai = 1 : length(angles)
                            tempAngleBinary = cell2mat(cellfun(@(x) ones(length(x.tpmTime{plane}) + 2 * posShift, 1) * (x.angle == angles(ai)), u.trials(tind), 'uniformoutput', false));
                            pTouchCountAngles{ai} = pTouchCount .* tempAngleBinary';
%                             rTouchCountAngles{ai} = rTouchCount .* tempAngleBinary';
%                             tTouchCountAngles{ai} = tTouchCount .* tempAngleBinary';
%                             pTouchDurationAngles{ai} = pTouchDuration .* tempAngleBinary';
%                             rTouchDurationAngles{ai} = rTouchDuration .* tempAngleBinary';
%                             tTouchDurationAngles{ai} = tTouchDuration .* tempAngleBinary';
%                             pTouchFramesAngles{ai} = pTouchFrames .* tempAngleBinary';
%                             rTouchFramesAngles{ai} = rTouchFrames .* tempAngleBinary';
%                             tTouchFramesAngles{ai} = tTouchFrames .* tempAngleBinary';
                        end
                        pTouchCountAngles{end} = pTouchCount;
%                         rTouchCountAngles{end} = rTouchCount;
%                         tTouchCountAngles{end} = tTouchCount;
%                         pTouchDurationAngles{end} = pTouchDuration;
%                         rTouchDurationAngles{end} = rTouchDuration;
%                         tTouchDurationAngles{end} = tTouchDuration;
%                         pTouchFramesAngles{end} = pTouchFrames;
%                         rTouchFramesAngles{end} = rTouchFrames;
%                         tTouchFramesAngles{end} = tTouchFrames;

%                         scPiezo = cell2mat(cellfun(@(x) [nan(1,posShift), 1, zeros(1,length(x.tpmTime{plane})-1), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        scPoleup = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.poleUpOnsetTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        scPoledown = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.poleDownOnsetTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        drinkOnset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.drinkingOnsetTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));

                        bLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(union(x.leftLickTime, x.rightLickTime), [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        lLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.leftLickTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        rLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.rightLickTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));

                        bLickOnset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.bothLickOnset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        lLickOnset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.leftLickOnset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        rLickOnset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.rightLickOnset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));

                        bLickOffset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.bothLickOffset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        lLickOffset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.leftLickOffset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        rLickOffset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.rightLickOffset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));

                        firstLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.firstLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        lastLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.lastLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        firstLeftLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.firstLeftLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        lastLeftLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.lastLeftLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        firstRightLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.firstRightLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        lastRightLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.lastRightLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        %%
                        whiskingOnsetCell = cell(1,length(tind));
                        whiskingAmpCell = cell(1,length(tind));
                        whiskingMaxAmpCell = cell(1,length(tind));
                        whiskingMidpointCell = cell(1,length(tind));

                        for ti = 1 : length(tind)
                            currTrial = u.trials{tind(ti)};
                            time = [0, currTrial.tpmTime{plane}];
                            wtimes = currTrial.whiskerTime;
                            [onsetFrame, amplitude, midpoint] = jkWhiskerOnsetNAmplitude(currTrial.theta, 5);
                            onsetTimes = onsetFrame*whiskerVideoFrameDuration/1000; % back to s
                            whiskingOnsetCell{ti} = [nan(1,posShift), histcounts(onsetTimes, time), nan(1,posShift)];

                            tempAmp = zeros(1,length(time)-1);        
                            tempMid = zeros(1,length(time)-1);
                            for i = 1 : length(tempAmp)
                                startInd = find(wtimes >= time(i), 1, 'first');
                                endInd = find(wtimes < time(i+1), 1, 'last');
                                tempAmp(i) = max(amplitude(startInd:endInd));
                                tempMid(i) = mean(midpoint(startInd:endInd));
                            end
                            tempMid(isnan(tempMid)) = deal(mode(tempMid(isfinite(tempMid))));
                            whiskingAmpCell{ti} = [nan(1,posShift), tempAmp, nan(1,posShift)];        
                            whiskingMidpointCell{ti} = [nan(1,posShift), tempMid, nan(1,posShift)];    
                        end
                        whiskingOnset = cell2mat(whiskingOnsetCell);
                        whiskingAmp = cell2mat(whiskingAmpCell);
                        whiskingOA = (whiskingOnset > 0).*whiskingAmp;
                        whiskingMidpoint = cell2mat(whiskingMidpointCell);    

                        %%
%                         tTouchCountMat = zeros(length(tTouchCount), (posShift + 1) * (length(angles)+1));
                        pTouchCountMat = zeros(length(pTouchCount), (posShift + 1) * (length(angles)+1));
%                         rTouchCountMat = zeros(length(rTouchCount), (posShift + 1) * (length(angles)+1));
%                         tTouchFramesMat = zeros(length(tTouchFrames),(posShift + 1) * (length(angles)+1));
%                         pTouchFramesMat = zeros(length(pTouchFrames),(posShift + 1) * (length(angles)+1));
%                         rTouchFramesMat = zeros(length(rTouchFrames), (posShift + 1) * (length(angles)+1));
%                         tTouchDurationMat = zeros(length(tTouchDuration), (posShift + 1) * (length(angles)+1));
%                         pTouchDurationMat = zeros(length(pTouchDuration), (posShift + 1) * (length(angles)+1));
%                         rTouchDurationMat = zeros(length(rTouchDuration), (posShift + 1) * (length(angles)+1));
% 
%                         scPiezoMat = zeros(length(scPiezo), posShift + 1);
                        scPoleUpMat = zeros(length(scPoleup), posShift + 1);
                        scPoleDownMat = zeros(length(scPoledown), posShift + 1);
                        drinkOnsetMat = zeros(length(drinkOnset), posShift + 1);
                        for i = 1 : posShift + 1
                            for ai = 1 : length(angles) + 1
%                                 tTouchCountMat(:,(i-1)*(length(angles)+1) + ai) = circshift(tTouchCountAngles{ai}, [0 i-1])';
                                pTouchCountMat(:,(i-1)*(length(angles)+1) + ai) = circshift(pTouchCountAngles{ai}, [0 i-1])';
%                                 rTouchCountMat(:,(i-1)*(length(angles)+1) + ai) = circshift(rTouchCountAngles{ai}, [0 i-1])';
%                                 tTouchFramesMat(:,(i-1)*(length(angles)+1) + ai) = circshift(tTouchFramesAngles{ai}, [0 i-1])';
%                                 pTouchFramesMat(:,(i-1)*(length(angles)+1) + ai) = circshift(pTouchFramesAngles{ai}, [0 i-1])';
%                                 rTouchFramesMat(:,(i-1)*(length(angles)+1) + ai) = circshift(rTouchFramesAngles{ai}, [0 i-1])';
%                                 tTouchDurationMat(:,(i-1)*(length(angles)+1) + ai) = circshift(tTouchDurationAngles{ai}, [0 i-1])';
%                                 pTouchDurationMat(:,(i-1)*(length(angles)+1) + ai) = circshift(pTouchDurationAngles{ai}, [0 i-1])';
%                                 rTouchDurationMat(:,(i-1)*(length(angles)+1) + ai) = circshift(rTouchDurationAngles{ai}, [0 i-1])';
                            end
%                             scPiezoMat(:,i) = circshift(scPiezo, [0 i-1])';
                            scPoleUpMat(:,i) = circshift(scPoleup, [0 i-1])';
                            scPoleDownMat(:,i) = circshift(scPoledown, [0 i-1])';
                            drinkOnsetMat(:,i) = circshift(drinkOnset, [0 i-1])';
                        end

                        whiskingOnsetMat = zeros(length(whiskingOnset), negShift + posShift + 1);
                        whiskingAmpMat = zeros(length(whiskingAmp), negShift + posShift + 1);
                        whiskingOAMat = zeros(length(whiskingOA), negShift + posShift + 1);
                        whiskingMidpointMat = zeros(length(whiskingMidpoint), negShift + posShift + 1);

                        bLickMat = zeros(length(bLick), negShift + posShift + 1);
                        lLickMat = zeros(length(lLick), negShift + posShift + 1);
                        rLickMat = zeros(length(rLick), negShift + posShift + 1);
                        bLickOnsetMat = zeros(length(bLickOnset), negShift + posShift + 1);
                        lLickOnsetMat = zeros(length(lLickOnset), negShift + posShift + 1);
                        rLickOnsetMat = zeros(length(rLickOnset), negShift + posShift + 1);
                        bLickOffsetMat = zeros(length(bLickOffset), negShift + posShift + 1);
                        lLickOffsetMat = zeros(length(lLickOffset), negShift + posShift + 1);
                        rLickOffsetMat = zeros(length(rLickOffset), negShift + posShift + 1);
                        firstLickMat = zeros(length(firstLick), negShift + posShift + 1);
                        lastLickMat = zeros(length(lastLick), negShift + posShift + 1);
                        firstLeftLickMat = zeros(length(firstLeftLick), negShift + posShift + 1);
                        lastLeftLickMat = zeros(length(lastLeftLick), negShift + posShift + 1);
                        firstRightLickMat = zeros(length(firstRightLick), negShift + posShift + 1);
                        lastRightLickMat = zeros(length(lastRightLick), negShift + posShift + 1);
                        for i = 1 : negShift + posShift + 1
                            whiskingOnsetMat(:,i) = circshift(whiskingOnset, [0 -negShift + i - 1])';
                            whiskingAmpMat(:,i) = circshift(whiskingAmp, [0 -negShift + i - 1])';
                            whiskingOAMat(:,i) = circshift(whiskingOA, [0 -negShift + i - 1])';
                            whiskingMidpointMat(:,i) = circshift(whiskingMidpoint, [0 -negShift + i - 1])';

                            bLickMat(:,i) = circshift(bLick, [0 -negShift + i - 1])';
                            lLickMat(:,i) = circshift(lLick, [0 -negShift + i - 1])';
                            rLickMat(:,i) = circshift(rLick, [0 -negShift + i - 1])';
                            bLickOnsetMat(:,i) = circshift(bLickOnset, [0 -negShift + i - 1])';
                            lLickOnsetMat(:,i) = circshift(lLickOnset, [0 -negShift + i - 1])';
                            rLickOnsetMat(:,i) = circshift(rLickOnset, [0 -negShift + i - 1])';
                            bLickOffsetMat(:,i) = circshift(bLickOffset, [0 -negShift + i - 1])';
                            lLickOffsetMat(:,i) = circshift(lLickOffset, [0 -negShift + i - 1])';
                            rLickOffsetMat(:,i) = circshift(rLickOffset, [0 -negShift + i - 1])';
                            firstLickMat(:,i) = circshift(firstLick, [0 -negShift + i - 1])';
                            lastLickMat(:,i) = circshift(lastLick, [0 -negShift + i - 1])';
                            firstLeftLickMat(:,i) = circshift(firstLeftLick, [0 -negShift + i - 1])';
                            lastLeftLickMat(:,i) = circshift(lastLeftLick, [0 -negShift + i - 1])';
                            firstRightLickMat(:,i) = circshift(firstRightLick, [0 -negShift + i - 1])';
                            lastRightLickMat(:,i) = circshift(lastRightLick, [0 -negShift + i - 1])';
                        end
%                         touchMat = [tTouchCountMat, pTouchCountMat, rTouchCountMat, tTouchFramesMat, pTouchFramesMat, rTouchFramesMat, tTouchDurationMat, pTouchDurationMat, rTouchDurationMat];
                        touchMat = [pTouchCountMat];                        
                        soundMat = [scPoleUpMat, scPoleDownMat];
                        drinkMat = drinkOnsetMat;
                        whiskingMat = [whiskingOnsetMat, whiskingAmpMat, whiskingOAMat, whiskingMidpointMat];
                        lickingMat = [bLickMat, lLickMat, rLickMat, bLickOnsetMat, lLickOnsetMat, rLickOnsetMat, bLickOffsetMat, lLickOffsetMat, rLickOffsetMat, ...
                                firstLickMat, lastLickMat, firstLeftLickMat, lastLeftLickMat, firstRightLickMat, lastRightLickMat];
%                         tempIM = [touchMat, soundMat, drinkMat, whiskingMat, lickingMat];
                        tempIM = touchMat;
                        trainingInputMat{(cgi-1)*4 + plane} = (tempIM - ones(size(tempIM,1),1)*min(tempIM)) ./ (ones(size(tempIM,1),1)*(max(tempIM) - min(tempIM)));
                        
                    end
                end

                %% design matrix for test 
                testInputMat = cell(8,1);
                for cgi = 1:2
                % for cgi = 1
                    tindcell = find(cellfun(@(x) ismember(1001+(cgi-1)*4000, x.neuindSession), u.trials));

                    tind = intersect(tindcell, testInd);
                    for plane = 1 : 4    
                %     for plane = 1
                        pTouchCount = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(cellfun(@(y) y(1), x.protractionTouchChunks), [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
%                         rTouchCount = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(cellfun(@(y) y(1), x.retractionTouchChunks), [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
%                         tTouchCount = pTouchCount + rTouchCount;
%                         pTouchCount(pTouchCount>0) = deal(1);
% 
%                         whiskerVideoFrameDuration = u.trials{tind(1)}.frameDuration * 1000; % in ms
% 
%                         pTouchDuration = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.protractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime{plane}]), nan(1,posShift)], ...
%                             u.trials(tind)','uniformoutput',false));
%                         rTouchDuration = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.retractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime{plane}]), nan(1,posShift)], ...
%                             u.trials(tind)','uniformoutput',false));
%                         tTouchDuration = pTouchDuration + rTouchDuration;
% 
%                         pTouchFrames = pTouchDuration;
%                         pTouchFrames(pTouchDuration > 0) = 1;
%                         rTouchFrames = rTouchDuration;
%                         rTouchFrames(rTouchDuration > 0) = 1;
%                         tTouchFrames = pTouchFrames + rTouchFrames;
% 
%                         pTouchCountAngles = cell(length(angles)+1,1);
%                         rTouchCountAngles = cell(length(angles)+1,1);
%                         tTouchCountAngles = cell(length(angles)+1,1);
%                         pTouchDurationAngles = cell(length(angles)+1,1);
%                         rTouchDurationAngles = cell(length(angles)+1,1);
%                         tTouchDurationAngles = cell(length(angles)+1,1);
%                         pTouchFramesAngles = cell(length(angles)+1,1);
%                         rTouchFramesAngles = cell(length(angles)+1,1);
%                         tTouchFramesAngles = cell(length(angles)+1,1);
                        for ai = 1 : length(angles)
                            tempAngleBinary = cell2mat(cellfun(@(x) ones(length(x.tpmTime{plane}) + 2 * posShift, 1) * (x.angle == angles(ai)), u.trials(tind), 'uniformoutput', false));
                            pTouchCountAngles{ai} = pTouchCount .* tempAngleBinary';
%                             rTouchCountAngles{ai} = rTouchCount .* tempAngleBinary';
%                             tTouchCountAngles{ai} = tTouchCount .* tempAngleBinary';
%                             pTouchDurationAngles{ai} = pTouchDuration .* tempAngleBinary';
%                             rTouchDurationAngles{ai} = rTouchDuration .* tempAngleBinary';
%                             tTouchDurationAngles{ai} = tTouchDuration .* tempAngleBinary';
%                             pTouchFramesAngles{ai} = pTouchFrames .* tempAngleBinary';
%                             rTouchFramesAngles{ai} = rTouchFrames .* tempAngleBinary';
%                             tTouchFramesAngles{ai} = tTouchFrames .* tempAngleBinary';
                        end
                        pTouchCountAngles{end} = pTouchCount;
%                         rTouchCountAngles{end} = rTouchCount;
%                         tTouchCountAngles{end} = tTouchCount;
%                         pTouchDurationAngles{end} = pTouchDuration;
%                         rTouchDurationAngles{end} = rTouchDuration;
%                         tTouchDurationAngles{end} = tTouchDuration;
%                         pTouchFramesAngles{end} = pTouchFrames;
%                         rTouchFramesAngles{end} = rTouchFrames;
%                         tTouchFramesAngles{end} = tTouchFrames;

%                         scPiezo = cell2mat(cellfun(@(x) [nan(1,posShift), 1, zeros(1,length(x.tpmTime{plane})-1), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        scPoleup = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.poleUpOnsetTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        scPoledown = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.poleDownOnsetTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        drinkOnset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.drinkingOnsetTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));

                        bLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(union(x.leftLickTime, x.rightLickTime), [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        lLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.leftLickTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        rLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.rightLickTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));

                        bLickOnset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.bothLickOnset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        lLickOnset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.leftLickOnset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        rLickOnset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.rightLickOnset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));

                        bLickOffset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.bothLickOffset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        lLickOffset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.leftLickOffset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        rLickOffset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.rightLickOffset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));

                        firstLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.firstLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        lastLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.lastLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        firstLeftLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.firstLeftLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        lastLeftLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.lastLeftLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        firstRightLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.firstRightLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        lastRightLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.lastRightLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        %%
                        whiskingOnsetCell = cell(1,length(tind));
                        whiskingAmpCell = cell(1,length(tind));
                        whiskingMaxAmpCell = cell(1,length(tind));
                        whiskingMidpointCell = cell(1,length(tind));

                        for ti = 1 : length(tind)
                            currTrial = u.trials{tind(ti)};
                            time = [0, currTrial.tpmTime{plane}];
                            wtimes = currTrial.whiskerTime;
                            [onsetFrame, amplitude, midpoint] = jkWhiskerOnsetNAmplitude(currTrial.theta, 5);
                            onsetTimes = onsetFrame*whiskerVideoFrameDuration/1000; % back to s
                            whiskingOnsetCell{ti} = [nan(1,posShift), histcounts(onsetTimes, time), nan(1,posShift)];

                            tempAmp = zeros(1,length(time)-1);        
                            tempMid = zeros(1,length(time)-1);
                            for i = 1 : length(tempAmp)
                                startInd = find(wtimes >= time(i), 1, 'first');
                                endInd = find(wtimes < time(i+1), 1, 'last');
                                tempAmp(i) = max(amplitude(startInd:endInd));
                                tempMid(i) = mean(midpoint(startInd:endInd));
                            end
                            tempMid(isnan(tempMid)) = deal(mode(tempMid(isfinite(tempMid))));
                            whiskingAmpCell{ti} = [nan(1,posShift), tempAmp, nan(1,posShift)];        
                            whiskingMidpointCell{ti} = [nan(1,posShift), tempMid, nan(1,posShift)];    
                        end
                        whiskingOnset = cell2mat(whiskingOnsetCell);
                        whiskingAmp = cell2mat(whiskingAmpCell);
                        whiskingOA = (whiskingOnset > 0).*whiskingAmp;
                        whiskingMidpoint = cell2mat(whiskingMidpointCell);

                        %%
%                         tTouchCountMat = zeros(length(tTouchCount), (posShift + 1) * (length(angles)+1));
                        pTouchCountMat = zeros(length(pTouchCount), (posShift + 1) * (length(angles)+1));
%                         rTouchCountMat = zeros(length(rTouchCount), (posShift + 1) * (length(angles)+1));
%                         tTouchFramesMat = zeros(length(tTouchFrames),(posShift + 1) * (length(angles)+1));
%                         pTouchFramesMat = zeros(length(pTouchFrames),(posShift + 1) * (length(angles)+1));
%                         rTouchFramesMat = zeros(length(rTouchFrames), (posShift + 1) * (length(angles)+1));
%                         tTouchDurationMat = zeros(length(tTouchDuration), (posShift + 1) * (length(angles)+1));
%                         pTouchDurationMat = zeros(length(pTouchDuration), (posShift + 1) * (length(angles)+1));
%                         rTouchDurationMat = zeros(length(rTouchDuration), (posShift + 1) * (length(angles)+1));

%                         scPiezoMat = zeros(length(scPiezo), posShift + 1);
                        scPoleUpMat = zeros(length(scPoleup), posShift + 1);
                        scPoleDownMat = zeros(length(scPoledown), posShift + 1);
                        drinkOnsetMat = zeros(length(drinkOnset), posShift + 1);
                        for i = 1 : posShift + 1
                            for ai = 1 : length(angles) + 1
%                                 tTouchCountMat(:,(i-1)*(length(angles)+1) + ai) = circshift(tTouchCountAngles{ai}, [0 i-1])';
                                pTouchCountMat(:,(i-1)*(length(angles)+1) + ai) = circshift(pTouchCountAngles{ai}, [0 i-1])';
%                                 rTouchCountMat(:,(i-1)*(length(angles)+1) + ai) = circshift(rTouchCountAngles{ai}, [0 i-1])';
%                                 tTouchFramesMat(:,(i-1)*(length(angles)+1) + ai) = circshift(tTouchFramesAngles{ai}, [0 i-1])';
%                                 pTouchFramesMat(:,(i-1)*(length(angles)+1) + ai) = circshift(pTouchFramesAngles{ai}, [0 i-1])';
%                                 rTouchFramesMat(:,(i-1)*(length(angles)+1) + ai) = circshift(rTouchFramesAngles{ai}, [0 i-1])';
%                                 tTouchDurationMat(:,(i-1)*(length(angles)+1) + ai) = circshift(tTouchDurationAngles{ai}, [0 i-1])';
%                                 pTouchDurationMat(:,(i-1)*(length(angles)+1) + ai) = circshift(pTouchDurationAngles{ai}, [0 i-1])';
%                                 rTouchDurationMat(:,(i-1)*(length(angles)+1) + ai) = circshift(rTouchDurationAngles{ai}, [0 i-1])';
                            end
%                             scPiezoMat(:,i) = circshift(scPiezo, [0 i-1])';
                            scPoleUpMat(:,i) = circshift(scPoleup, [0 i-1])';
                            scPoleDownMat(:,i) = circshift(scPoledown, [0 i-1])';
                            drinkOnsetMat(:,i) = circshift(drinkOnset, [0 i-1])';
                        end

                        whiskingOnsetMat = zeros(length(whiskingOnset), negShift + posShift + 1);
                        whiskingAmpMat = zeros(length(whiskingAmp), negShift + posShift + 1);
                        whiskingOAMat = zeros(length(whiskingOA), negShift + posShift + 1);
                        whiskingMidpointMat = zeros(length(whiskingMidpoint), negShift + posShift + 1);

                        bLickMat = zeros(length(bLick), negShift + posShift + 1);
                        lLickMat = zeros(length(lLick), negShift + posShift + 1);
                        rLickMat = zeros(length(rLick), negShift + posShift + 1);
                        bLickOnsetMat = zeros(length(bLickOnset), negShift + posShift + 1);
                        lLickOnsetMat = zeros(length(lLickOnset), negShift + posShift + 1);
                        rLickOnsetMat = zeros(length(rLickOnset), negShift + posShift + 1);
                        bLickOffsetMat = zeros(length(bLickOffset), negShift + posShift + 1);
                        lLickOffsetMat = zeros(length(lLickOffset), negShift + posShift + 1);
                        rLickOffsetMat = zeros(length(rLickOffset), negShift + posShift + 1);
                        firstLickMat = zeros(length(firstLick), negShift + posShift + 1);
                        lastLickMat = zeros(length(lastLick), negShift + posShift + 1);
                        firstLeftLickMat = zeros(length(firstLeftLick), negShift + posShift + 1);
                        lastLeftLickMat = zeros(length(lastLeftLick), negShift + posShift + 1);
                        firstRightLickMat = zeros(length(firstRightLick), negShift + posShift + 1);
                        lastRightLickMat = zeros(length(lastRightLick), negShift + posShift + 1);
                        for i = 1 : negShift + posShift + 1
                            whiskingOnsetMat(:,i) = circshift(whiskingOnset, [0 -negShift + i - 1])';
                            whiskingAmpMat(:,i) = circshift(whiskingAmp, [0 -negShift + i - 1])';
                            whiskingOAMat(:,i) = circshift(whiskingOA, [0 -negShift + i - 1])';
                            whiskingMidpointMat(:,i) = circshift(whiskingMidpoint, [0 -negShift + i - 1])';

                            bLickMat(:,i) = circshift(bLick, [0 -negShift + i - 1])';
                            lLickMat(:,i) = circshift(lLick, [0 -negShift + i - 1])';
                            rLickMat(:,i) = circshift(rLick, [0 -negShift + i - 1])';
                            bLickOnsetMat(:,i) = circshift(bLickOnset, [0 -negShift + i - 1])';
                            lLickOnsetMat(:,i) = circshift(lLickOnset, [0 -negShift + i - 1])';
                            rLickOnsetMat(:,i) = circshift(rLickOnset, [0 -negShift + i - 1])';
                            bLickOffsetMat(:,i) = circshift(bLickOffset, [0 -negShift + i - 1])';
                            lLickOffsetMat(:,i) = circshift(lLickOffset, [0 -negShift + i - 1])';
                            rLickOffsetMat(:,i) = circshift(rLickOffset, [0 -negShift + i - 1])';
                            firstLickMat(:,i) = circshift(firstLick, [0 -negShift + i - 1])';
                            lastLickMat(:,i) = circshift(lastLick, [0 -negShift + i - 1])';
                            firstLeftLickMat(:,i) = circshift(firstLeftLick, [0 -negShift + i - 1])';
                            lastLeftLickMat(:,i) = circshift(lastLeftLick, [0 -negShift + i - 1])';
                            firstRightLickMat(:,i) = circshift(firstRightLick, [0 -negShift + i - 1])';
                            lastRightLickMat(:,i) = circshift(lastRightLick, [0 -negShift + i - 1])';
                        end
%                         touchMat = [tTouchCountMat, pTouchCountMat, rTouchCountMat, tTouchFramesMat, pTouchFramesMat, rTouchFramesMat, tTouchDurationMat, pTouchDurationMat, rTouchDurationMat];
                        touchMat = [pTouchCountMat];
                        soundMat = [scPoleUpMat, scPoleDownMat];
                        drinkMat = drinkOnsetMat;
                        whiskingMat = [whiskingOnsetMat, whiskingAmpMat, whiskingOAMat, whiskingMidpointMat];
                        lickingMat = [bLickMat, lLickMat, rLickMat, bLickOnsetMat, lLickOnsetMat, rLickOnsetMat, bLickOffsetMat, lLickOffsetMat, rLickOffsetMat, ...
                                firstLickMat, lastLickMat, firstLeftLickMat, lastLeftLickMat, firstRightLickMat, lastRightLickMat];
%                         tempIM = [touchMat, soundMat, drinkMat, whiskingMat, lickingMat];
                        tempIM = touchMat;
                        testInputMat{(cgi-1)*4 + plane} = (tempIM - ones(size(tempIM,1),1)*min(tempIM)) ./ (ones(size(tempIM,1),1)*(max(tempIM) - min(tempIM)));

                    end
                end

                %%
                touchInd = 1 : size(touchMat,2);
                soundInd = max(touchInd) + 1 : max(touchInd) + size(soundMat,2);
                rewardInd = max(soundInd) + 1 : max(soundInd) + size(drinkMat,2);
                whiskingInd = max(rewardInd) + 1 : max(rewardInd) + size(whiskingMat,2);
                lickInd = max(whiskingInd) + 1 : max(whiskingInd) + size(lickingMat,2);

                indPartial{1} = touchInd;
                indPartial{2} = soundInd;
                indPartial{3} = rewardInd;
                indPartial{4} = whiskingInd;
                indPartial{5} = lickInd;
%         end
%         if ~exist(savefnPredictors, 'file')                
%             save(savefnPredictors, '*InputMat', 'indPartial', '*Group', '*Tn');
%         end
        %%
%             rtest(ri).fitInd = cell(length(u.cellNums),1); % parameters surviving lasso in training set
%             rtest(ri).fitCoeffs = cell(length(u.cellNums),1); % intercept + coefficients of the parameters in training set
%             rtest(ri).fitCoeffInds = nan(length(u.cellNums),6); % first column is dummy
%             
%             rtest(ri).fitResults = zeros(length(u.cellNums), 6);
%                 % fitResult(:,1) if full fitting is significant (compared to null model), 0 if not
%                 % fitResult(:,2) for touchInd, compared to full fitting. if excluding touch
%                 % is significantly less fit, then 1, 0 otherwise
%                 % fitResult(:,3) for sound, (:,4) for reward, (:,5) for whisking, and (:,6) for licking    
%             rtest(ri).devExplained = zeros(length(u.cellNums),1);
            
    
            fitInd = cell(length(u.cellNums),1); % parameters surviving lasso in training set
            fitCoeffs = cell(length(u.cellNums),1); % intercept + coefficients of the parameters in training set
            fitCoeffInds = nan(length(u.cellNums),6); % first column is a dummy
            
            fitResults = zeros(length(u.cellNums), 6); % fitting result from test set
            devExplained = zeros(length(u.cellNums),1); % deviance explained from test set
            cvDev = zeros(length(u.cellNums),1); % deviance explained from training set
                    
            parfor cellnum = 1 : length(u.cellNums)
    %         for cellnum = 1
    %         ci = 0;
    %         for cellnum = 1:division:length(u.cellNums)
    %             ci = ci + 1;
    %         for cellnum = 1
            %     cellnum = 1;
                fitCoeffInd = zeros(1,6);
                
%                 fprintf('Mouse JK%03d session S%02d Loop %d: Running cell %d/%d \n', mouse, session, ri,cellnum, length(u.cellNums));

                fprintf('Mouse JK%03d session S%02d: Running cell %d/%d \n', mouse, session,cellnum, length(u.cellNums));
                
                cID = u.cellNums(cellnum);
    
                % find out trial indices for this specific cell
                tindcell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
    
                tind = intersect(tindcell, trainingInd);
                % find out row number of this cell
                cind = find(u.trials{tind(1)}.neuindSession == cID);
                planeInd = floor(cID/1000);
    
                spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cind,:), nan(1,posShift)], u.trials(tind)','uniformoutput',false));                
                %%
                finiteInd = intersect(find(isfinite(spkTrain)), find(isfinite(sum(trainingInputMat{planeInd},2))));
                input = trainingInputMat{planeInd}(finiteInd,:);
                spk = spkTrain(finiteInd)';
%                 spk = (spk - min(spk)) / (max(spk) - min(spk)); % min-max normalization
    
                cv = cvglmnet(input, spk, 'poisson', glmnetOpt);
            %     [BFull, FitInfoFull] = lassoglm(trainingInputMat{planeInd}, spkTrain, 'poisson', 'Lambda', logspace(-4,2), 'CV', lambdaCV, 'Options', opt);
                %% survived coefficients
                iLambda = find(cv.lambda == cv.lambda_1se);
                fitCoeffs{cellnum} = [cv.glmnet_fit.a0(iLambda);cv.glmnet_fit.beta(:,iLambda)];
                coeffInds = find(cv.glmnet_fit.beta(:,iLambda));                
    %             rtest(ri).fitInd{cellnum} = coeffInds;
                fitInd{cellnum} = coeffInds;
                for i = 1 : length(indPartial)
                    if sum(ismember(indPartial{i},coeffInds)>0)
                        fitCoeffInd(i + 1) = 1;
                    else
                        fitCoeffInd(i + 1) = 0;
                    end
                end
                    
            % %     find(ismember(touchInd, mincoefs))
            % %     find(ismember(soundInd, mincoefs))
            % %     find(ismember(rewardInd, mincoefs))
            % %     find(ismember(whiskingInd, mincoefs))
            % %     find(ismember(lickInd, mincoefs))
    
                %% test
                cID = u.cellNums(cellnum);
    
                % find out trial indices for this specific cell
                tindcell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
    
                tind = intersect(tindcell, testInd);
                % find out row number of this cell
                cind = find(u.trials{tind(1)}.neuindSession == cID);
                planeInd = floor(cID/1000);
    
                spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cind,:), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                spkTest = spkTest';
                finiteInd = intersect(find(isfinite(spkTest)), find(isfinite(sum(testInputMat{planeInd},2))));
                spkTest = spkTest(finiteInd)';
%                 spkTest = (spkTest - min(spkTest)) / (max(spkTest) - min(spkTest));
                %% (1) if the full model is significant
                fitResult = zeros(1,6);
    
                model = exp([ones(length(finiteInd),1),testInputMat{planeInd}(finiteInd,:)]*[cv.glmnet_fit.a0(iLambda); cv.glmnet_fit.beta(:,iLambda)]);
                mu = mean(spkTest); % null poisson parameter
                nullLogLikelihood = sum(log(poisspdf(spkTest,mu)));
                fullLogLikelihood = sum(log(poisspdf(spkTest',exp(model))));
                saturatedLogLikelihood = sum(log(poisspdf(spkTest,spkTest)));
                devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
                dfFullNull = length(coeffInds);
                devExplained(cellnum) = 1 - (saturatedLogLikelihood - fullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
                cvDev(cellnum) = cv.glmnet_fit.dev(iLambda);
                if devianceFullNull > chi2inv(1-pThresholdNull, dfFullNull)
                    fitResult(1) = 1;
%                     %% (2) test without each parameter (as a group)                
%                     for pi = 1 : 5
%                         if find(ismember(coeffInds, indPartial{pi}))
%                             if all(ismember(coeffInds, indPartial{pi}))
%                                 fitResult(pi+1) = 1;
%                             else
%                                 tempTrainInput = trainingInputMat{planeInd}(:,setdiff(coeffInds,indPartial{pi}));
%                                 tempTestInput = testInputMat{planeInd}(finiteInd,setdiff(coeffInds,indPartial{pi}));
%                                 cvPartial = cvglmnet(tempTrainInput(finiteInd,:), spk, 'poisson', partialGlmOpt);
%                                 iLambda = find(cvPartial.lambda == cvPartial.lambda_1se);
%                                 partialLogLikelihood = sum(log(poisspdf(spkTest', exp([ones(length(finiteInd),1), tempTestInput] * [cvPartial.glmnet_fit.a0(iLambda); cvPartial.glmnet_fit.beta(:,iLambda)]))));
%                                 devianceFullPartial = 2*(fullLogLikelihood - partialLogLikelihood);
%                                 dfFullPartial = dfFullNull - cvPartial.glmnet_fit.df(iLambda);
%                                 if devianceFullPartial > chi2inv(1-pThresholdPartial, dfFullPartial)
%                                     fitResult(pi+1) = 1;
%                                 end
%                             end
%                         end
%                     end
                end
                
                
                corVal = corr(model, spkTest');
                
                shuffleCorVals = zeros(numShuffle,1);
                for iShuffle = 1 : numShuffle
                    spkTestShuffled = spkTest(randperm(length(spkTest)));
                    shuffleCorVals(iShuffle) = corr(model,spkTestShuffled');
                end
                
                if corVal > prctile(shuffleCorVals,95)
                    fitResult(2) = 1;
                end
                
                fitResults(cellnum,:) = fitResult;
                fitCoeffInds(cellnum,:) = fitCoeffInd;
            end % end of parfor cellnum
            
%             rtest(ri).fitInd = fitInd; % parameters surviving lasso in training set
%             rtest(ri).fitCoeffs = fitCoeffs;
%             rtest(ri).fitCoeffInds = fitCoeffInds; % first column is dummy
%             
%             rtest(ri).fitResults = fitResults;        
%             % fitResult(:,1) if full fitting is significant (compared to null model), 0 if not
%             % fitResult(:,2) for touchInd, compared to full fitting. if excluding touch
%             % is significantly less fit, then 1, 0 otherwise
%             % fitResult(:,3) for sound, (:,4) for reward, (:,5) for whisking, and (:,6) for licking
%     
%             rtest(ri).devExplained = devExplained;
%             rtest(ri).cvDev = cvDev;
            
            save(savefnResult, 'fitInd', 'fitCoeffs', 'fitCoeffInds', 'fitResults', 'devExplained', 'cvDev', '*InputMat', 'indPartial', '*Group', '*Tn');

%         end % of ri. random group selection index
        push_myphone(sprintf('GLM done for JK%03d S%02d', mouse, session))

    end
end



