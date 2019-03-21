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
%     - pTouchCount: # of protraction touches within the frame (parameter).
%     Retraction touches removed due to limited number of trials (most times none).
%     
%     - pTouchFrames: protraction touch frames (binary)
%
%     - pTouchDuration: protraction touch duration within each tpm frame (parameter, in ms)

%       Up to here, have one as all angles and add each angles (total number of predictors: 1 + length(angles)
%
%     - scPoleup: pole up sound cue onset (binary)
%     - scPoledown: pole down sound cue onset (binary)
%     Piezo sound cue removed because it is always at the first frame, and cannot be dealt with NaN paddings
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
%     - lLick: left licks within the frame (parameter)
%     - rLick: right licks within the frame (parameter)
%     
%     - lLickOnset: the frame where left lick onset happened (binary; each bout is calculated as 1 s interval)
%     - lLickOffset: the frame where left lick offset happened (binary; each bout is calculated as 1 s interval)
%     - rLickOnset
%     - rLickoffset
% 
%     - firstLick: 
%     - firstLeftLick: the frame where the first lick of the trial happened (binary)
%     - lastLeftLick: the frame where the last lick of the trial happened (binary)
%     - firstRightLick
%     - lastRightLick

baseDir = 'D:\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 

            repetition = 3;
            startRepetition = 1;
% sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]};
% errorCell = {{[],[224]},{[],[]},{[],[]},{[],[]},{[]},{[]},{[1211,1972],[1286]},{[]},{[],[605, 676, 740, 755, 811]},{[]},{[]},{[]}};
% errorCell = {{[],[]},{[],[]},{[],[]},{[],[]},{[]},{[]},{[2042,2059],[]},{[]},{[],[]},{[]},{[]},{[]}};
errorCell = {{[],[]},{[],[]},{[],[]},{[],[]},{[]},{[]},{[],[]},{[]},{[],[]},{[]},{[]},{[]}};
%%
% mice = [25,27,30];
% sessions = {[17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 
% sessions = {[19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 
% for mi = 1 : length(mice)
for mi = [1:4,7,9]
    for si = 1:length(sessions{mi})
%     for si = 2
        errorCellSession = errorCell{mi}{si};
    
        poolobj = gcp('nocreate');
        if poolobj.SpmdEnabled == 0
            error('SpmdEnabled turned to false at #1');
        end
        
        mouse = mice(mi);
        session = sessions{mi}(si);
        
        posShiftTouch = 2;
        posShiftSound = 3;
        posShiftReward = 4;
        posShiftWhisking = 4;
        posShiftLicking = 1;
        posShift = 4; % maximum posShift
        negShift = 2;
        testPortion = 0.3; % 30 % test set
        pThresholdNull = 0.05;
        pThresholdPartial = 0.05;
%         lickBoutInterval = 1; % licks separated by 1 s regarded as different licking bouts

        glmnetOpt = glmnetSet;
        glmnetOpt.standardize = 0; % do the standardization at the level of predictors, including both training and test
        glmnetOpt.alpha = 0.95;
        
        partialGlmOpt = glmnetOpt;
        partialGlmOpt.alpha = 0;
        lambdaCV = 5; % cross-validation fold number

        dn = sprintf('%s%03d',baseDir,mouse);
        ufn = sprintf('UberJK%03dS%02d.mat', mouse, session);
        cd(dn)
        if exist(ufn,'file')
            load(ufn)
        else
            u = Uber.buildUberArray(mouse, session);
        end
        frameRate = u.frameRate;

        savefnResult = sprintf('glmResponseType_JK%03dS%02d_glmnet_m40',mouse, session); % m(n) meaining method(n)

%             %% pre-processing for lick onset and offset
%             % regardless of licking alternating, each l and r has it's own lick onset and offset. both licking, just take the union
% 
%             v = cell(length(u.trials),1);
% 
%             for ui = 1 : length(u.trials)
%                 bothLickTime = union(u.trials{ui}.leftLickTime, u.trials{ui}.rightLickTime);
%                 if length(bothLickTime) == 1
%                     v{ui}.bothLickOnset = bothLickTime;
%                     v{ui}.bothLickOffset = bothLickTime;
%                     v{ui}.firstLick = bothLickTime;
%                     if u.trials{ui}.response < 1
%                         v{ui}.lastLick = bothLickTime;
%                     else
%                         v{ui}.lastLick = [];
%                     end
%                 elseif length(bothLickTime) > 1
%                     onsets = find(diff(bothLickTime) > lickBoutInterval);
%                     if isempty(onsets)
%                         v{ui}.bothLickOnset = bothLickTime(1);
%                         v{ui}.bothLickOffset = bothLickTime(end);
%                     else
%                         v{ui}.bothLickOnset = bothLickTime([1; onsets+1]);
%                         v{ui}.bothLickOffset = bothLickTime([onsets; end]);
%                     end
%                     v{ui}.firstLick = bothLickTime(1);
%                     if u.trials{ui}.response < 1
%                         v{ui}.lastLick = bothLickTime(end);
%                     else
%                         v{ui}.lastLick = [];
%                     end
%                 else
%                     v{ui}.bothLickOnset = [];
%                     v{ui}.bothLickOffset = [];
%                     v{ui}.firstLick = [];
%                     v{ui}.lastLick = [];
%                 end
% 
% 
%                 if length(u.trials{ui}.leftLickTime) == 1
%                     v{ui}.leftLickOnset = u.trials{ui}.leftLickTime;
%                     v{ui}.leftLickOffset = u.trials{ui}.leftLickTime;
%                     v{ui}.firstLeftLick = u.trials{ui}.leftLickTime;
%                     if u.trials{ui}.response < 1 
%                         if isempty(u.trials{ui}.rightLickTime)
%                             v{ui}.lastLeftLick = u.trials{ui}.leftLickTime;
%                         elseif u.trials{ui}.rightLickTime(end) < u.trials{ui}.leftLickTime(end)
%                             v{ui}.lastLeftLick = u.trials{ui}.leftLickTime;
%                         else
%                             v{ui}.lastLeftLick = [];
%                         end
%                     else
%                         v{ui}.lastLeftLick = [];
%                     end
%                 elseif length(u.trials{ui}.leftLickTime) > 1        
%                     onsets = find(diff(u.trials{ui}.leftLickTime) > lickBoutInterval);
%                     if isempty(onsets)
%                         v{ui}.leftLickOnset = u.trials{ui}.leftLickTime(1);
%                         v{ui}.leftLickOffset = u.trials{ui}.leftLickTime(end);                    
%                     else
%                         v{ui}.leftLickOnset = u.trials{ui}.leftLickTime([1; onsets+1]);
%                         v{ui}.leftLickOffset = u.trials{ui}.leftLickTime([onsets; end]);
%                     end
%                     v{ui}.firstLeftLick = u.trials{ui}.leftLickTime(1);
%                     if u.trials{ui}.response < 1 
%                         if isempty(u.trials{ui}.rightLickTime)
%                             v{ui}.lastLeftLick = u.trials{ui}.leftLickTime;
%                         elseif u.trials{ui}.rightLickTime(end) < u.trials{ui}.leftLickTime(end)
%                             v{ui}.lastLeftLick = u.trials{ui}.leftLickTime;
%                         else
%                             v{ui}.lastLeftLick = [];
%                         end
%                     else
%                         v{ui}.lastLeftLick = [];
%                     end
%                 else
%                     v{ui}.leftLickOnset = [];
%                     v{ui}.leftLickOffset = [];
%                     v{ui}.firstLeftLick = [];
%                     v{ui}.lastLeftLick = [];
%                 end
% 
% 
%                 if length(u.trials{ui}.rightLickTime) == 1
%                     v{ui}.rightLickOnset = u.trials{ui}.rightLickTime;
%                     v{ui}.rightLickOffset = u.trials{ui}.rightLickTime;
%                     v{ui}.firstRightLick = u.trials{ui}.rightLickTime;
%                     if u.trials{ui}.response < 1 
%                         if isempty(u.trials{ui}.leftLickTime)
%                             v{ui}.lastRightLick = u.trials{ui}.rightLickTime;
%                         elseif u.trials{ui}.leftLickTime(end) < u.trials{ui}.rightLickTime(end)
%                             v{ui}.lastRightLick = u.trials{ui}.rightLickTime;
%                         else
%                             v{ui}.lastRightLick = [];
%                         end
%                     else
%                         v{ui}.lastRightLick = [];
%                     end
%                 elseif length(u.trials{ui}.rightLickTime) > 1        
%                     onsets = find(diff(u.trials{ui}.rightLickTime) > lickBoutInterval);
%                     if isempty(onsets)
%                         v{ui}.rightLickOnset = u.trials{ui}.rightLickTime(1);
%                         v{ui}.rightLickOffset = u.trials{ui}.rightLickTime(end);
%                     else
%                         v{ui}.rightLickOnset = u.trials{ui}.rightLickTime([1; onsets+1]);
%                         v{ui}.rightLickOffset = u.trials{ui}.rightLickTime([onsets; end]);
%                     end
%                     v{ui}.firstRightLick = u.trials{ui}.rightLickTime(1);
%                     if u.trials{ui}.response < 1 
%                         if isempty(u.trials{ui}.leftLickTime)
%                             v{ui}.lastRightLick = u.trials{ui}.rightLickTime;
%                         elseif u.trials{ui}.leftLickTime(end) < u.trials{ui}.rightLickTime(end)
%                             v{ui}.lastRightLick = u.trials{ui}.rightLickTime;
%                         else
%                             v{ui}.lastRightLick = [];
%                         end
%                     else
%                         v{ui}.lastRightLick = [];
%                     end
%                 else
%                     v{ui}.rightLickOnset = [];
%                     v{ui}.rightLickOffset = [];
%                     v{ui}.firstRightLick = [];
%                     v{ui}.lastRightLick = [];
%                 end
%                 v{ui}.tpmTime = u.trials{ui}.tpmTime;
%             end


%             %% repetition test
%     %         division = 20;


        for ri = startRepetition : repetition % repetition index
                %% divide into training set and test set (70%, 30%)
                % based on the animal touched or not, the choice (same as the result since I'm going to mix the pole angles, so right, wrong, and miss), pole angles (2 or 7), and the distance (if there were multiple distances)
                % in this order, make trees, and take 30% of the leaves (or equivalently, take all the possible intersections and take 30%)
                
                angles = unique(cellfun(@(x) x.angle, u.trials));
                distances = unique(cellfun(@(x) x.distance, u.trials));

                touchGroup = cell(2,1);
                choiceGroup = cell(3,1);
                angleGroup = cell(length(angles),1);
                distanceGroup = cell(length(distances),1);
                timeGroup = cell(3,1); % dividing whole session into 5 different time points

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
%                     for rti = 1 : length(rtouchGroup)
                        for ci = 1 : length(choiceGroup)
                            for ai = 1 : length(angleGroup)
                                for di = 1 : length(distanceGroup)
                                        tempTn = intersect(ptouchGroup{pti}, intersect(choiceGroup{ci}, intersect(angleGroup{ai}, distanceGroup{di})));
                                        if ~isempty(tempTn)
                                            tempTn = tempTn(randperm(length(tempTn)));
                                            testTn = [testTn; tempTn(1:round(length(tempTn)*0.3))];
                                        end
                                end
                            end
                        end
%                     end
                end
                %
                totalTn = u.trialNums;
                [~,testInd] = ismember(testTn, totalTn);

                trainingTn = setdiff(totalTn, testTn);
                [~,trainingInd] = ismember(trainingTn, totalTn);
                %% Design matrices
                % standardized using all the trials
                allPredictors = cell(8,1);
                allPredictorsMean = cell(8,1);
                allPredictorsStd = cell(8,1);
                nani = cell(8,1);
                trainingPredictorInd = cell(8,1);
                testPredictorInd = cell(8,1);
                trainingInputMat = cell(8,1);
                testInputMat = cell(8,1);
                
                for cgi = 1:2 % cell group index
                % for cgi = 1
                    tindcell = find(cellfun(@(x) ismember(1001+(cgi-1)*4000, x.neuindSession), u.trials));

                    tind = tindcell;
                    for plane = 1 : 4    
                %     for plane = 1
                        trainingPredictorInd{(cgi-1)*4 + plane} = cell2mat(cellfun(@(x) (ones(1,length(x.tpmTime{plane})+posShift*2)) * ismember(x.trialNum, trainingTn), u.trials(tind)','uniformoutput',false));
                        testPredictorInd{(cgi-1)*4 + plane} = cell2mat(cellfun(@(x) (ones(1,length(x.tpmTime{plane})+posShift*2)) * ismember(x.trialNum, testTn), u.trials(tind)','uniformoutput',false));
                        
                        pTouchCount = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(cellfun(@(y) y(1), x.protractionTouchChunks), [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        pTouchFrame = pTouchCount;
                        pTouchFrame(pTouchFrame > 0) = 1;
%                         pTouchDuration = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.protractionTouchChunks, 'uniformoutput', false)), [0, x.tpmTime{plane}]), nan(1,posShift)], ...
%                             u.trials(tind)','uniformoutput',false));
                        
%                         pTouchCountAngles = cell(length(angles)+1,1);
%                         pTouchDurationAngles = cell(length(angles)+1,1);
                        pTouchFrameAngles = cell(length(angles)+1,1);
                        for ai = 1 : length(angles)
                            tempAngleBinary = cell2mat(cellfun(@(x) ones(length(x.tpmTime{plane}) + 2 * posShift, 1) * (x.angle == angles(ai)), u.trials(tind), 'uniformoutput', false));
                            pTouchFrameAngles{ai} = pTouchFrame .* tempAngleBinary';
%                             pTouchCountAngles{ai} = pTouchCount .* tempAngleBinary';
%                             pTouchDurationAngles{ai} = pTouchDuration .* tempAngleBinary';
                        end
                        pTouchFrameAngles{end} = pTouchFrame;
%                         pTouchCountAngles{end} = pTouchCount;
%                         pTouchDurationAngles{end} = pTouchDuration;

%                         scPiezo = cell2mat(cellfun(@(x) [nan(1,posShift), 1, zeros(1,length(x.tpmTime{plane})-1), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        scPoleup = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.poleUpOnsetTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
%                         scPoledown = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.poleDownOnsetTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
%                         drinkOnset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.drinkingOnsetTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        drink = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.drinkingOnsetTime:0.03:x.drinkingOnsetTime+1, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        drink(drink>0) = 1;
                        
%                         drinkL = cell2mat(cellfun(@(x) [nan(1,posShift), (x.trialType(1) == 'l') * histcounts(x.drinkingOnsetTime:0.03:x.drinkingOnsetTime+1, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
%                         drinkL(drinkL>0) = 1;
%                         
%                         drinkR = cell2mat(cellfun(@(x) [nan(1,posShift), (x.trialType(1) == 'r') * histcounts(x.drinkingOnsetTime:0.03:x.drinkingOnsetTime+1, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
%                         drinkR(drinkR>0) = 1;
                        
                        drinkAngles = cell(length(angles)+1,1);
                        for ai = 1 : length(angles)
                            tempAngleBinary = cell2mat(cellfun(@(x) ones(length(x.tpmTime{plane}) + 2 * posShift, 1) * (x.angle == angles(ai)), u.trials(tind), 'uniformoutput', false));
                            drinkAngles{ai} = drink .* tempAngleBinary';
                        end
                        drinkAngles{end} = drink;
                        
                        
                        
                        lLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.leftLickTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        rLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.rightLickTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                        
                        lLickFrame = lLick;
                        lLickFrame(lLickFrame>0) = 1;
                        
                        rLickFrame = rLick;
                        rLickFrame(rLickFrame>0) = 1;
                        
                        
                        
                        
%                         drinkL = drink + lLickFrame;
%                         drinkL(drinkL<2) = 0;
%                         drinkL(drinkL>0) = 1;
%                         
%                         drinkR = drink + rLickFrame;
%                         drinkR(drinkR<2) = 0;
%                         drinkR(drinkR>0) = 1;
%                         
%                         lLick(drink>0) = 0;
%                         rLick(drink>0) = 0;

%                         lLickOnset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.leftLickOnset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
%                         rLickOnset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.rightLickOnset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        
%                         lLickOffset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.leftLickOffset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
%                         rLickOffset = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.rightLickOffset, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        
%                         firstLeftLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.firstLeftLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
%                         firstRightLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.firstRightLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        
%                         lastLeftLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.lastLeftLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
%                         lastRightLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.lastRightLick, [0, x.tpmTime{plane}]), nan(1,posShift)], v(tind)','uniformoutput',false));
                        %%
                        whiskingOnsetCell = cell(1,length(tind));
%                         whiskingFrameCell = cell(1,length(tind));
                        whiskingAmplitudeCell = cell(1,length(tind));
                        whiskingMidpointCell = cell(1,length(tind));

                        for ti = 1 : length(tind)
                            currTrial = u.trials{tind(ti)};
                            time = [0, currTrial.tpmTime{plane}];
                            wtimes = currTrial.whiskerTime;
                            [onsetFrame, amplitude, midpoint] = jkWhiskerOnsetNAmplitude(currTrial.theta);
                            whiskerVideoFrameDuration = u.trials{tind(1)}.frameDuration; % in s
                            onsetTimes = onsetFrame*whiskerVideoFrameDuration; % back to s
                            tempOnset = histcounts(onsetTimes, time);
%                             tempFrame = tempOnset;
%                             tempFrame(tempFrame > 0) = 1;
                            whiskingOnsetCell{ti} = [nan(1,posShift), tempOnset, nan(1,posShift)];
%                             whiskingFrameCell{ti} = [nan(1,posShift), tempFrame, nan(1,posShift)];

                            tempMid = zeros(1,length(time)-1);
                            tempAmp = zeros(1,length(time)-1);
                            for i = 1 : length(tempMid)
                                startInd = find(wtimes >= time(i), 1, 'first');
                                endInd = find(wtimes < time(i+1), 1, 'last');
                                tempMid(i) = mean(midpoint(startInd:endInd));
                                tempAmp(i) = max(amplitude(startInd:endInd));
                            end
                            tempMid(isnan(tempMid)) = deal(mode(tempMid(isfinite(tempMid))));
                            tempAmp(isnan(tempAmp)) = deal(mode(tempAmp(isfinite(tempAmp))));
                            whiskingMidpointCell{ti} = [nan(1,posShift), tempMid, nan(1,posShift)];
                            whiskingAmplitudeCell{ti} = [nan(1,posShift), tempAmp, nan(1,posShift)];
                        end
                        whiskingOnset = cell2mat(whiskingOnsetCell);
%                         whiskingFrame = cell2mat(whiskingFrameCell);
                        whiskingMidpoint = cell2mat(whiskingMidpointCell);
                        whiskingAmplitude = cell2mat(whiskingAmplitudeCell);

                        %%
%                         pTouchCountMat = zeros(length(pTouchCount), (posShiftTouch + 1) * (length(angles)+1));
                        pTouchFrameMat = zeros(length(pTouchFrame), (posShiftTouch + 1) * (length(angles)+1));
%                         pTouchDurationMat = zeros(length(pTouchDuration), (posShiftTouch + 1) * (length(angles)+1));

%                         scPiezoMat = zeros(length(scPiezo), posShiftSound + 1);
                        scPoleUpMat = zeros(length(scPoleup), posShiftSound + 1);
%                         scPoleDownMat = zeros(length(scPoledown), posShiftSound + 1);
%                         drinkOnsetMat = zeros(length(drinkOnset), posShiftReward + 1);
%                         drinkLMat = zeros(length(drinkL), posShiftReward + 1);
%                         drinkRMat = zeros(length(drinkR), posShiftReward + 1);
%                         drinkMat = zeros(length(drink), posShiftReward + 1);
                        drinkAnglesMat = zeros(length(drink), (posShiftReward + 1) * (length(angles)+1));
                        for i = 1 : posShiftTouch + 1
                            for ai = 1 : length(angles) + 1
%                                 pTouchCountMat(:,(i-1)*(length(angles)+1) + ai) = circshift(pTouchCountAngles{ai}, [0 i-1])';
                                pTouchFrameMat(:,(i-1)*(length(angles)+1) + ai) = circshift(pTouchFrameAngles{ai}, [0 i-1])';
%                                 pTouchDurationMat(:,(i-1)*(length(angles)+1) + ai) = circshift(pTouchDurationAngles{ai}, [0 i-1])';
                            end
                        end
                        for i = 1 : posShiftSound + 1
%                             scPiezoMat(:,i) = circshift(scPiezo, [0 i-1])';
                            scPoleUpMat(:,i) = circshift(scPoleup, [0 i-1])';
%                             scPoleDownMat(:,i) = circshift(scPoledown, [0 i-1])';
                        end
                        for i = 1 : posShiftReward + 1
%                             drinkOnsetMat(:,i) = circshift(drinkOnset, [0 i-1])';
%                             drinkLMat(:,i) = circshift(drinkL, [0 i-1])';
%                             drinkRMat(:,i) = circshift(drinkR, [0 i-1])';
%                             drinkMat(:,i) = circshift(drink, [0 i-1])';
                            for ai = 1 : length(angles) + 1
                                drinkAnglesMat(:,(i-1)*(length(angles)+1) + ai) = circshift(drinkAngles{ai}, [0 i-1])';
                            end                                
                        end

                        whiskingOnsetMat = zeros(length(whiskingOnset), negShift + posShiftWhisking + 1);
%                         whiskingFrameMat = zeros(length(whiskingOnset), negShift + posShiftMotor + 1);
                        
                        whiskingAmplitudeMat = zeros(length(whiskingAmplitude), negShift + posShiftWhisking + 1);
                        whiskingMidpointMat = zeros(length(whiskingMidpoint), negShift + posShiftWhisking + 1);

                        lLickMat = zeros(length(lLick), negShift + posShiftLicking + 1);
                        rLickMat = zeros(length(rLick), negShift + posShiftLicking + 1);
%                         lLickOnsetMat = zeros(length(lLickOnset), negShift + posShiftMotor + 1);
%                         rLickOnsetMat = zeros(length(rLickOnset), negShift + posShiftMotor + 1);
%                         lLickOffsetMat = zeros(length(lLickOffset), negShift + posShiftMotor + 1);
%                         rLickOffsetMat = zeros(length(rLickOnset), negShift + posShiftMotor + 1);
%                         firstLeftLickMat = zeros(length(firstLeftLick), negShift + posShiftMotor + 1);
%                         firstRightLickMat = zeros(length(firstRightLick), negShift + posShiftMotor + 1);
%                         lastLeftLickMat = zeros(length(lastLeftLick), negShift + posShiftMotor + 1);
%                         lastRightLickMat = zeros(length(lastRightLick), negShift + posShiftMotor + 1);
                        for i = 1 : negShift + posShiftWhisking + 1
                            whiskingOnsetMat(:,i) = circshift(whiskingOnset, [0 -negShift + i - 1])';
%                             whiskingFrameMat(:,i) = circshift(whiskingFrame, [0 -negShift + i - 1])';
                            whiskingMidpointMat(:,i) = circshift(whiskingMidpoint, [0 -negShift + i - 1])';
                            whiskingAmplitudeMat(:,i) = circshift(whiskingAmplitude, [0 -negShift + i - 1])';                            
                        end
                        
                        for i = 1 : negShift + posShiftLicking + 1
                            lLickMat(:,i) = circshift(lLick, [0 -negShift + i - 1])';
                            rLickMat(:,i) = circshift(rLick, [0 -negShift + i - 1])';
                            
%                             lLickOnsetMat(:,i) = circshift(lLickOnset, [0 -negShift + i - 1])';
%                             rLickOnsetMat(:,i) = circshift(rLickOnset, [0 -negShift + i - 1])';
%                             
%                             lLickOffsetMat(:,i) = circshift(lLickOffset, [0 -negShift + i - 1])';
%                             rLickOffsetMat(:,i) = circshift(rLickOffset, [0 -negShift + i - 1])';
%                             
%                             firstLeftLickMat(:,i) = circshift(firstLeftLick, [0 -negShift + i - 1])';
%                             firstRightLickMat(:,i) = circshift(firstRightLick, [0 -negShift + i - 1])';
%                             
%                             lastLeftLickMat(:,i) = circshift(lastLeftLick, [0 -negShift + i - 1])';
%                             lastRightLickMat(:,i) = circshift(lastRightLick, [0 -negShift + i - 1])';
                        end
                        
%                         touchMat = [tTouchCountMat, pTouchCountMat, rTouchCountMat, tTouchFramesMat, pTouchFramesMat, rTouchFramesMat, tTouchDurationMat, pTouchDurationMat, rTouchDurationMat];
                        touchMat = [pTouchFrameMat];
%                         soundMat = [scPiezoMat, scPoleUpMat, scPoleDownMat];
                        soundMat = [scPoleUpMat];
%                         drinkMat = drinkOnsetMat;
%                         drinkMat = [drinkLMat, drinkRMat];
                        drinkMat = drinkAnglesMat;
                        whiskingMat = [whiskingOnsetMat, whiskingAmplitudeMat, whiskingMidpointMat];
%                         lickingMat = [lLickMat, rLickMat, lLickOnsetMat, rLickOnsetMat, lLickOffsetMat, rLickOffsetMat, firstLeftLickMat, firstRightLickMat, lastLeftLickMat, lastRightLickMat];
                        lickingMat = [lLickMat, rLickMat];
%                         allPredictors{(cgi-1)*4 + plane} = [touchMat, soundMat, drinkMat, whiskingMat, lickingMat];
                        allPredictors{(cgi-1)*4 + plane} = [touchMat, soundMat, drinkMat, whiskingMat, lickingMat];
                        nani{(cgi-1)*4 + plane} = find(nanstd(allPredictors{(cgi-1)*4 + plane})==0);
                        allPredictorsMean{(cgi-1)*4 + plane} = nanmean(allPredictors{(cgi-1)*4 + plane});
                        allPredictorsStd{(cgi-1)*4 + plane} = nanstd(allPredictors{(cgi-1)*4 + plane});
                        % normalization of all predictors
                        allPredictors{(cgi-1)*4 + plane} = (allPredictors{(cgi-1)*4 + plane} - nanmean(allPredictors{(cgi-1)*4 + plane})) ./ nanstd(allPredictors{(cgi-1)*4 + plane});
                        allPredictors{(cgi-1)*4 + plane}(:,nani{(cgi-1)*4 + plane}) = deal(0);
                        trainingInputMat{(cgi-1)*4 + plane} = allPredictors{(cgi-1)*4 + plane}(find(trainingPredictorInd{(cgi-1)*4 + plane}),:);
                        testInputMat{(cgi-1)*4 + plane} = allPredictors{(cgi-1)*4 + plane}(find(testPredictorInd{(cgi-1)*4 + plane}),:);
                    end
                end

                %%
                touchInd = 1 : size(touchMat,2);
                soundInd = max(touchInd) + 1 : max(touchInd) + size(soundMat,2);
                rewardInd = max(soundInd) + 1 : max(soundInd) + size(drinkMat,2);
%                 rewardInd = max(touchInd) + 1 : max(touchInd) + size(drinkMat,2);
                whiskingInd = max(rewardInd) + 1 : max(rewardInd) + size(whiskingMat,2);
                lickInd = max(whiskingInd) + 1 : max(whiskingInd) + size(lickingMat,2);

                indPartial{1} = touchInd;
                indPartial{2} = soundInd;
                indPartial{3} = rewardInd;
                indPartial{4} = whiskingInd;
                indPartial{5} = lickInd;
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
    
            cIDAll = u.cellNums;
            numCell = length(cIDAll); 
            fitInd = cell(numCell,1); % parameters surviving lasso in training set
            fitCoeffs = cell(numCell,1); % intercept + coefficients of the parameters in training set
            fitCoeffInds = nan(numCell,6); % first column is a dummy
            
            fitResults = zeros(numCell, 6); % fitting result from test set
            fitDeviance = zeros(numCell,1);
            fitCorrelation = zeros(numCell,1);
            fitCorrPval = zeros(numCell,1);
                
            fitDevExplained = zeros(numCell,1); % deviance explained from test set
            fitCvDev = zeros(numCell,1); % deviance explained from training set
            fitLambda = zeros(numCell,1);
            fitDF = zeros(numCell,1);
            started = zeros(numCell,1);
            done = zeros(numCell,1);
            cellTime = zeros(numCell,1);
          
            tindcellAll = cell(numCell,1);
            cindAll = zeros(numCell,1);
            planeIndAll = zeros(numCell,1);
            iTrainAll = cell(numCell,1);
            iTestAll = cell(numCell,1);
            for i = 1 : numCell
                tindcellAll{i} = find(cellfun(@(x) ismember(cIDAll(i), x.neuindSession), u.trials));
                cindAll(i) = find(u.trials{tindcellAll{i}(1)}.neuindSession == cIDAll(i));
                planeIndAll(i) = floor(cIDAll(i)/1000);
                iTrainAll{i} = intersect(tindcellAll{i}, trainingInd);
                iTestAll{i} = intersect(tindcellAll{i}, testInd);
            end
            spikeAll = cellfun(@(x) x.spk, u.trials, 'uniformoutput', false);
            
            
            
            poolobj = gcp('nocreate');
            if poolobj.SpmdEnabled == 0
                error('SpmdEnabled turned to false at #2');
            end
            
            parfor cellnum = 1 : numCell
%             for cellnum = 102, 127, (212 convergence error), 221, 658
    %         ci = 0;
    %         for cellnum = 1:division:length(u.cellNums)
    %             ci = ci + 1;
    %         for cellnum = 1
            %     cellnum = 1;
            if ~ismember(numCell, errorCellSession)
                cellTimeStart = tic;
                fitCoeffInd = zeros(1,6);
                
                fprintf('Mouse JK%03d session S%02d Loop %d: Running cell %d/%d \n', mouse, session, ri, cellnum, numCell);
%                 fprintf('Mouse JK%03d session S%02d: Running cell %d/%d \n', mouse, session,cellnum, numCell);
                started(cellnum) = cellnum;
                
                iTrain = iTrainAll{cellnum};
                cind = cindAll(cellnum);
                planeInd = planeIndAll(cellnum);
    
                spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTrain)','uniformoutput',false));

                finiteIndTrain = intersect(find(isfinite(spkTrain)), find(isfinite(sum(trainingInputMat{planeInd},2))));
                input = trainingInputMat{planeInd}(finiteIndTrain,:);
                spkTrain = spkTrain(finiteIndTrain)';
    
                cv = cvglmnet(input, spkTrain, 'poisson', glmnetOpt, [], lambdaCV);
                %% survived coefficients
                fitLambda(cellnum) = cv.lambda_1se;
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

                %% test
                iTest = iTestAll{cellnum};                
                spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x(cind,:), nan(1,posShift)], spikeAll(iTest)','uniformoutput',false));
                spkTest = spkTest';
                finiteIndTest = intersect(find(isfinite(spkTest)), find(isfinite(sum(testInputMat{planeInd},2))));
                spkTest = spkTest(finiteIndTest)';
                %% (1) if the full model is significant
                fitResult = zeros(1,6);
    
                model = exp([ones(length(finiteIndTest),1),testInputMat{planeInd}(finiteIndTest,:)]*[cv.glmnet_fit.a0(iLambda); cv.glmnet_fit.beta(:,iLambda)]);
                mu = mean(spkTest); % null poisson parameter
                nullLogLikelihood = sum(log(poisspdf(spkTest,mu)));
                fullLogLikelihood = sum(log(poisspdf(spkTest',model)));
                saturatedLogLikelihood = sum(log(poisspdf(spkTest,spkTest)));
                devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
                fitDeviance(cellnum) = devianceFullNull;
                [fitCorrelation(cellnum), fitCorrPval(cellnum)] = corr(spkTest', model);            
                dfFullNull = length(coeffInds);                
                fitDF(cellnum) = dfFullNull;
                fitDevExplained(cellnum) = 1 - (saturatedLogLikelihood - fullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
                fitCvDev(cellnum) = cv.glmnet_fit.dev(iLambda);

                if devianceFullNull > chi2inv(1-pThresholdNull, dfFullNull)
                    fitResult(1) = 1;
                    
                    %% (2) test without each parameter (as a group)                
%                     for pi = 1 : 5
%                         if find(ismember(coeffInds, indPartial{pi}))
%                             if all(ismember(coeffInds, indPartial{pi}))
%                                 fitResult(pi+1) = 1;
%                                 break
%                             else
%                                 tempTrainInput = trainingInputMat{planeInd}(:,setdiff(coeffInds,indPartial{pi}));
%                                 tempTestInput = testInputMat{planeInd}(finiteIndTest,setdiff(coeffInds,indPartial{pi}));
%                                 cvPartial = cvglmnet(tempTrainInput(finiteIndTrain,:), spkTrain, 'poisson', partialGlmOpt, [], lambdaCV);
%                                 iLambda = find(cvPartial.lambda == cvPartial.lambda_1se);
%                                 partialLogLikelihood = sum(log(poisspdf(spkTest', exp([ones(length(finiteIndTest),1), tempTestInput] * [cvPartial.glmnet_fit.a0(iLambda); cvPartial.glmnet_fit.beta(:,iLambda)]))));
%                                 devianceFullPartial = 2*(fullLogLikelihood - partialLogLikelihood);
%                                 dfFullPartial = dfFullNull - length(setdiff(coeffInds, indPartial{pi}));
%                                 if devianceFullPartial > chi2inv(1-pThresholdPartial, dfFullPartial)
%                                     fitResult(pi+1) = 1;
%                                 end
%                             end
%                         end
%                     end
                end
                
                fitResults(cellnum,:) = fitResult;
                fitCoeffInds(cellnum,:) = fitCoeffInd;
                done(cellnum) = cellnum;
                cellTime(cellnum) = toc(cellTimeStart);
            end
            end % end of parfor cellnum
%%
            save(sprintf('%s_R%02d',savefnResult, ri), 'fit*', 'allPredictors', '*InputMat', 'indPartial', '*Group', '*Tn', 'lambdaCV', '*Opt', 'done', 'pThreshold*', '*Shift', 'cellTime', 'testInd', 'trainingInd', 'cIDAll');
%             push_myphone(sprintf('Lasso GLM done for JK%03d S%02d Loop #%03d', mouse, session, ri))
        end % of ri. random group selection index
%         push_myphone(sprintf('Lasso GLM done for JK%03d S%02d', mouse, session))

    end
end



