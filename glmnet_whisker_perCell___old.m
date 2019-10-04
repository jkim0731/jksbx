% whisker glm 
% because of frequence crashes in glmnet, it's better to have each cell saved within parfor,
% and collect only the ones that are calculated correctly.
% Run it cell by cell, with parfor.
% This is because it's a few cell stuck that makes the whole process really long.
% 10 iterations of 10 reps, if the saved file (success) is less than 10, then the cell is an error cell.
% Depending on the iteration and the number of cores, one can run this code in multiple MATLAB windows.
% To reduce I/O burden in NAS, use local drive for temporary save.

% Modified from previous version 'glmnet_whisker_lasso'


% 2019/10/03 JK


%%
%% 
% For now, follow file saving format as before. 
% Need to make it more efficient later.
%%
%%

% Extracting input matrices for GLM analysis in each neuron from an Uber_2padArray u
% 
% Dependency: 
%     - Uber class
%     - jkWhiskerOnsetNAmplitude
%     - parfun_glmnet_whisker_perCell
% 
% inputs: 
%     - mouse (as in number)
%     - session (as in number) 
%     - cellnum (1~length(total number of cells)),
%     - nShift (number of frames to shift, either forward or backward. Default: 3)
%
% outputs: 
%     - cid: cell id (1000~8999)
%     - frameRate
%     - spk: spikes (vector. Padded with NaN's of length nShift before and after each trial)
%     
%     % sensory variables: shift backward only
%     % Same length as spk. 
% whisker touch variables upon touch. 2019/04/09 JK
%     - maxDtheta
%     - maxDphi
%     - maxDkappaV
%     - maxDkappaH
%     - maxSlideDistance
%     - maxDuration
% 
%     - thetaAtTouch
%     - phiAtTouch
%     - kappaHAtTouch
%     - kappaVAtTouch
%     - arcLengthAtTouch
%     - touchCounts
%
%     - scPoleup: pole up sound cue onset (binary)
%     - scPoledown: pole down sound cue onset (binary)
%     Piezo sound cue removed because it is always at the first frame, and cannot be dealt with NaN paddings
%
%     - drinkOnset: drinking onset (binary)
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
% 
% For all cells
% load 'spk' from "JKoooSooangle_tuning.mat"
% 
% Use whisker touch variables instead of touch variables
% To compare their explaining strenght 
% 2019/04/10 JK

%% basic settings
clear

baseDir = 'Y:\Whiskernas\JK\suite2p\';
% baseDir = 'D:\TPM\JK\suite2p\';

localDir = 'C:\JK\tempDataForGLM\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,10],[3,21],[1,17],[7],[2],[1,23],[3],[3,21],[3],[3],[3]}; 
repetition = 10;



poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool(repetition, 'SpmdEnabled', true);
elseif poolobj.SpmdEnabled == 0
    myCluster = parcluster('local');
    delete(myCluster.Jobs)
    clear myCluster
    pause(1) % just in case...
    parpool(repetition, 'SpmdEnabled', true);
end



% errorCell: optional
% errorCell = {{[92,103,219,220],[],[]},...
%     {[],[],[]},...
%     {[],[391],[]},...
%     {[],[],[]},...
%     {[]},...
%     {[316,631]},...
%     {[],[12;26;35;73;78;79;83;87;89;93;98;108;109;110;112;113;114;115;116;117;119;123;130;132;133;139;141;142;143;145;148;150;151;152;153;164;170;171;172;176;177;179;181;185;190;201;205;206;211;212;219;223;231;233;246;249;264;292;293;298;299;300;305;310;313;316;319;322;339;346;350;360;361;364;366;368;376;377;381;385;391;393;394;398;399;405;408;409;411;413;417;420;422;430;434;435;436;438;441;442;446;542;593;596;608;609;615;616;619;620;623;627;827;834;845;894;905;957;981;1097;1102;1122;1131;1140;1155;1166;1209;1238;1372;1822;1846;1912],...
%         [117,139,163,437]},...
%     {[2017,2103,2121]},...
%     {[],[176,763],[160,966]},...
%     {[]},...
%     {[]},...
%     {[]}};

glmnetOpt = glmnetSet;
glmnetOpt.standardize = 0; % do the standardization at the level of predictors, including both training and test
glmnetOpt.alpha = 0.95;
lambdaCV = 5; % cross-validation fold number


        
posShiftTouch = 2;
posShiftSound = 3;
posShiftReward = 4;
posShiftWhisking = 4;
posShiftLicking = 1;
posShift = 4; % maximum posShift
negShift = 2;
testPortion = 0.3; % 30 % test set


for mi = 1
    for si = 1
        mouse = mice(mi);
        session = sessions{mi}(si);
        
        dn = sprintf('%s%03d',baseDir,mouse);
        ufn = sprintf('UberJK%03dS%02d_NC.mat', mouse, session);
        cd(dn)
        load(ufn, 'u')
        frameRate = u.frameRate;

        savefnResult = sprintf('glmWhisker_lasso_allCell_NC_JK%03dS%02d',mouse, session); % m(n) meaining method(n)

        %% building universal data set (including both training set and test set)
        %% divide into training set and test set (70%, 30%)
        % based on the animal touched or not, the choice (same as the result since I'm going to mix the pole angles, so right, wrong, and miss), pole angles (2 or 7), and the distance (if there were multiple distances)
        % in this order, make trees, and take 30% of the leaves (or equivalently, take all the possible intersections and take 30%)

        angles = unique(cellfun(@(x) x.angle, u.trials));
        distances = unique(cellfun(@(x) x.distance, u.trials));

        touchGroup = cell(2,1);
        choiceGroup = cell(3,1);
        angleGroup = cell(length(angles),1);
        distanceGroup = cell(length(distances),1);

        ptouchGroup{1} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) length(x.protractionTouchChunksByWhisking), u.trials))));
        ptouchGroup{2} = setdiff(u.trialNums, ptouchGroup{1});

        choiceGroup{1} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.response == 1, u.trials))));
        choiceGroup{2} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.response == 0, u.trials))));
        choiceGroup{3} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.response == -1, u.trials))));

        for i = 1 : length(angles)
            angleGroup{i} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.angle == angles(i), u.trials))));
        end

        for i = 1 : length(distances)
            distanceGroup{i} = cellfun(@(x) x.trialNum, u.trials(find(cellfun(@(x) x.distance == distances(i), u.trials))));
        end

        stratificationGroups = {ptouchGroup, choiceGroup, angleGroup, distanceGroup};
        %% Design matrices
        % standardized using all the trials
        allPredictors = cell(8,1);
        allPredictorsRaw = cell(8,1);
        allPredictorsMean = cell(8,1);
        allPredictorsStd = cell(8,1);
        nani = cell(8,1);
        for cgi = 1:2 % cell group index
            tindcell = find(cellfun(@(x) ismember(1001+(cgi-1)*4000, x.neuindSession), u.trials));

            tind = tindcell;
            for plane = 1 : 4    
                pTouchCount = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(cellfun(@(y) x.whiskerTime(y(1)), x.protractionTouchChunksByWhisking), [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                pTouchFrame = pTouchCount;
                pTouchFrame(pTouchFrame > 0) = 1;
                pTouchFrameAngles = cell(length(angles)+1,1);
                for ai = 1 : length(angles)
                    tempAngleBinary = cell2mat(cellfun(@(x) ones(length(x.tpmTime{plane}) + 2 * posShift, 1) * (x.angle == angles(ai)), u.trials(tind), 'uniformoutput', false));
                    pTouchFrameAngles{ai} = pTouchFrame .* tempAngleBinary';
                end
                pTouchFrameAngles{end} = pTouchFrame;
                scPoleup = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.poleUpOnsetTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                drinkL = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.drinkingOnsetTime, [0, x.tpmTime{plane}]) * strcmp(x.choice, 'l'), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                drinkR = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.drinkingOnsetTime, [0, x.tpmTime{plane}]) * strcmp(x.choice, 'r'), nan(1,posShift)], u.trials(tind)','uniformoutput',false));                        
                lLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.leftLickTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
                rLick = cell2mat(cellfun(@(x) [nan(1,posShift), histcounts(x.rightLickTime, [0, x.tpmTime{plane}]), nan(1,posShift)], u.trials(tind)','uniformoutput',false));

                %%
                % whisking
                whiskingOnsetCell = cell(1,length(tind));
                whiskingAmplitudeCell = cell(1,length(tind));
                whiskingMidpointCell = cell(1,length(tind));

                % whisker touch variables during touch
                maxDkappaHCell = cell(1,length(tind));
                maxDkappaVCell = cell(1,length(tind));
                maxDthetaCell = cell(1,length(tind));
                maxDphiCell = cell(1,length(tind));
                maxSlideDistanceCell = cell(1,length(tind));
                maxDurationCell = cell(1,length(tind));
                % whisker touch onset variables
                thetaAtTouchCell = cell(1,length(tind));
                phiAtTouchCell = cell(1,length(tind));
                kappaHAtTouchCell = cell(1,length(tind));
                kappaVAtTouchCell = cell(1,length(tind));
                arcLengthAtTouchCell = cell(1,length(tind));
                % + touchCount later

                % use this to confirm it matches with previous protractionTouchFrames calculation
                touchFrameConfirmCell = cell(1,length(tind));
                for ti = 1 : length(tind)
                    currTrial = u.trials{tind(ti)};
                    time = [0, currTrial.tpmTime{plane}];                            
                    wtimes = currTrial.whiskerTime;
                    if iscell(wtimes)
                        wtimes = wtimes{1};
                    end

                    % whisking allocation
                    [onsetFrame, amplitude, midpoint] = jkWhiskerOnsetNAmplitude(currTrial.theta);
                    whiskerVideoFrameDuration = u.trials{tind(1)}.frameDuration; % in s
                    onsetTimes = onsetFrame*whiskerVideoFrameDuration; % back to s
                    tempOnset = histcounts(onsetTimes, time);
                    whiskingOnsetCell{ti} = [nan(1,posShift), tempOnset, nan(1,posShift)];

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

                    % whisker variables allocation  
                    tempDtheta = zeros(1,length(time)-1);
                    tempDphi = zeros(1,length(time)-1);
                    tempDkappaH = zeros(1,length(time)-1);
                    tempDkappaV = zeros(1,length(time)-1);
                    tempSlideDistance = zeros(1,length(time)-1);
                    tempDuration = zeros(1,length(time)-1);
                    tempThetaAtTouch = zeros(1,length(time)-1);
                    tempPhiAtTouch = zeros(1,length(time)-1);
                    tempKappaHAtTouch = zeros(1,length(time)-1);
                    tempKappaVAtTouch = zeros(1,length(time)-1);
                    tempArcLengthAtTouch = zeros(1,length(time)-1);

                    tempTouchFramesForConfirm = zeros(1,length(time)-1);
                    if ~isempty(currTrial.protractionTouchChunksByWhisking)
                        % assign all whisker touch variables to the
                        % touch onset timepoints in tpm time
                        touchOnsetTimes = cellfun(@(x) wtimes(x(1)), currTrial.protractionTouchChunksByWhisking);
                        touchHistCounts = histcounts(touchOnsetTimes, time);
                        touchFrames = find(touchHistCounts); % use this to confirm it matches with previous protractionTouchFrames calculation
                        tempTouchFramesForConfirm(touchFrames) = deal(1);
                        cumsumTouchFrames = [0, cumsum(touchHistCounts(touchFrames))];
                        for tfi = 1 : length(touchFrames)
                            tempInds = cumsumTouchFrames(tfi)+1:cumsumTouchFrames(tfi+1);
                            tempDtheta(touchFrames(tfi)) = nansum(currTrial.protractionTouchDThetaByWhisking(tempInds));
                            tempDphi(touchFrames(tfi)) = nansum(currTrial.protractionTouchDPhiByWhisking(tempInds));
                            tempDkappaH(touchFrames(tfi)) = nansum(currTrial.protractionTouchDKappaHByWhisking(tempInds));
                            tempDkappaV(touchFrames(tfi)) = nansum(currTrial.protractionTouchDKappaVByWhisking(tempInds));
                            tempSlideDistance(touchFrames(tfi)) = nansum(currTrial.protractionTouchSlideDistanceByWhisking(tempInds));
                            tempDuration(touchFrames(tfi)) = nansum(currTrial.protractionTouchDurationByWhisking(tempInds));

                            tempThetaAtTouch(touchFrames(tfi)) = nansum(cellfun(@(x) currTrial.theta(x(1)), currTrial.protractionTouchChunksByWhisking(tempInds)));
                            tempPhiAtTouch(touchFrames(tfi)) = nansum(cellfun(@(x) currTrial.phi(x(1)), currTrial.protractionTouchChunksByWhisking(tempInds)));
                            tempKappaHAtTouch(touchFrames(tfi)) = nansum(cellfun(@(x) currTrial.kappaH(x(1)), currTrial.protractionTouchChunksByWhisking(tempInds)));
                            tempKappaVAtTouch(touchFrames(tfi)) = nansum(cellfun(@(x) currTrial.kappaV(x(1)), currTrial.protractionTouchChunksByWhisking(tempInds)));
                            tempArcLengthAtTouch(touchFrames(tfi)) = nansum(cellfun(@(x) currTrial.arcLength(x(1)), currTrial.protractionTouchChunksByWhisking(tempInds)));
                        end
                    end
                    maxDthetaCell{ti} = [nan(1,posShift), tempDtheta, nan(1,posShift)];
                    maxDphiCell{ti} = [nan(1,posShift), tempDphi, nan(1,posShift)];
                    maxDkappaHCell{ti} = [nan(1,posShift), tempDkappaH, nan(1,posShift)];
                    maxDkappaVCell{ti} = [nan(1,posShift), tempDkappaV, nan(1,posShift)];
                    maxSlideDistanceCell{ti} = [nan(1,posShift), tempSlideDistance, nan(1,posShift)];
                    maxDurationCell{ti} = [nan(1,posShift), tempDuration, nan(1,posShift)];
                    thetaAtTouchCell{ti} = [nan(1,posShift), tempThetaAtTouch, nan(1,posShift)];
                    phiAtTouchCell{ti} = [nan(1,posShift), tempPhiAtTouch, nan(1,posShift)];
                    kappaHAtTouchCell{ti} = [nan(1,posShift), tempKappaHAtTouch, nan(1,posShift)];
                    kappaVAtTouchCell{ti} = [nan(1,posShift), tempKappaVAtTouch, nan(1,posShift)];
                    arcLengthAtTouchCell{ti} = [nan(1,posShift), tempArcLengthAtTouch, nan(1,posShift)];

                    touchFrameConfirmCell{ti} = [nan(1,posShift), tempTouchFramesForConfirm, nan(1,posShift)];
                end
                whiskingOnset = cell2mat(whiskingOnsetCell);
                whiskingMidpoint = cell2mat(whiskingMidpointCell);
                whiskingAmplitude = cell2mat(whiskingAmplitudeCell);

                maxDtheta = cell2mat(maxDthetaCell);
                maxDphi = cell2mat(maxDphiCell);
                maxDkappaH = cell2mat(maxDkappaHCell);
                maxDkappaV = cell2mat(maxDkappaVCell);
                maxSlideDistance = cell2mat(maxSlideDistanceCell);
                maxDuration = cell2mat(maxDurationCell);
                thetaAtTouch = cell2mat(thetaAtTouchCell);
                phiAtTouch = cell2mat(phiAtTouchCell);
                kappaHAtTouch = cell2mat(kappaHAtTouchCell);
                kappaVAtTouch = cell2mat(kappaVAtTouchCell);
                arcLengthAtTouch = cell2mat(arcLengthAtTouchCell);

                touchFrameConfirm = cell2mat(touchFrameConfirmCell);

                %% check validity of whisker touch variable allocation
                if touchFrameConfirm ~= pTouchFrame
                    error('Whisker touch variable allocation is wrong')
                end
                %%
                pTouchFrameMat = zeros(length(pTouchFrame), (posShiftTouch + 1) * (length(angles)+1)); % leave this for now, just in case

                scPoleUpMat = zeros(length(scPoleup), posShiftSound + 1);
                drinkLMat = zeros(length(drinkL), posShiftReward + 1);
                drinkRMat = zeros(length(drinkR), posShiftReward + 1);
                for i = 1 : posShiftTouch + 1
                    for ai = 1 : length(angles) + 1
                        pTouchFrameMat(:,(i-1)*(length(angles)+1) + ai) = circshift(pTouchFrameAngles{ai}, [0 i-1])';
                    end
                end
                for i = 1 : posShiftSound + 1
                    scPoleUpMat(:,i) = circshift(scPoleup, [0 i-1])';
                end
                for i = 1 : posShiftReward + 1
                    drinkLMat(:,i) = circshift(drinkL, [0 i-1])';
                    drinkRMat(:,i) = circshift(drinkR, [0 i-1])';
                end

                whiskingOnsetMat = zeros(length(whiskingOnset), negShift + posShiftWhisking + 1);
                whiskingAmplitudeMat = zeros(length(whiskingAmplitude), negShift + posShiftWhisking + 1);
                whiskingMidpointMat = zeros(length(whiskingMidpoint), negShift + posShiftWhisking + 1);

                lLickMat = zeros(length(lLick), negShift + posShiftLicking + 1);
                rLickMat = zeros(length(rLick), negShift + posShiftLicking + 1);

                for i = 1 : negShift + posShiftWhisking + 1
                    whiskingOnsetMat(:,i) = circshift(whiskingOnset, [0 -negShift + i - 1])';
                    whiskingMidpointMat(:,i) = circshift(whiskingMidpoint, [0 -negShift + i - 1])';
                    whiskingAmplitudeMat(:,i) = circshift(whiskingAmplitude, [0 -negShift + i - 1])';                            
                end

                for i = 1 : negShift + posShiftLicking + 1
                    lLickMat(:,i) = circshift(lLick, [0 -negShift + i - 1])';
                    rLickMat(:,i) = circshift(rLick, [0 -negShift + i - 1])';
                end

                maxDthetaMat = zeros(length(pTouchFrame), (posShiftTouch + 1));
                maxDphiMat = zeros(length(pTouchFrame), (posShiftTouch + 1));
                maxDkappaHMat = zeros(length(pTouchFrame), (posShiftTouch + 1));
                maxDkappaVMat = zeros(length(pTouchFrame), (posShiftTouch + 1));
                maxSlideDistanceMat = zeros(length(pTouchFrame), (posShiftTouch + 1));
                maxDurationMat = zeros(length(pTouchFrame), (posShiftTouch + 1));
                thetaAtTouchMat = zeros(length(pTouchFrame), (posShiftTouch + 1));
                phiAtTouchMat = zeros(length(pTouchFrame), (posShiftTouch + 1));
                kappaHAtTouchMat = zeros(length(pTouchFrame), (posShiftTouch + 1));
                kappaVAtTouchMat = zeros(length(pTouchFrame), (posShiftTouch + 1));
                arcLengthAtTouchMat = zeros(length(pTouchFrame), (posShiftTouch + 1));
                touchCountMat = zeros(length(pTouchFrame), (posShiftTouch + 1));
                for i = 1 : posShiftTouch + 1
                    maxDthetaMat(:,i) = circshift(maxDtheta, [0 i-1]);
                    maxDphiMat(:,i) = circshift(maxDphi, [0 i-1]);
                    maxDkappaHMat(:,i) = circshift(maxDkappaH, [0 i-1]);
                    maxDkappaVMat(:,i) = circshift(maxDkappaV, [0 i-1]);
                    maxSlideDistanceMat(:,i) = circshift(maxSlideDistance, [0 i-1]);
                    maxDurationMat(:,i) = circshift(maxDuration, [0 i-1]);
                    thetaAtTouchMat(:,i) = circshift(thetaAtTouch, [0 i-1]);
                    phiAtTouchMat(:,i) = circshift(phiAtTouch, [0 i-1]);
                    kappaHAtTouchMat(:,i) = circshift(kappaHAtTouch, [0 i-1]);
                    kappaVAtTouchMat(:,i) = circshift(kappaVAtTouch, [0 i-1]);
                    arcLengthAtTouchMat(:,i) = circshift(arcLengthAtTouch, [0 i-1]);
                    touchCountMat(:,i) = circshift(pTouchCount, [0 i-1]);
                end

        %                         touchMat = [pTouchFrameMat];
                soundMat = [scPoleUpMat];
                drinkMat = [drinkLMat, drinkRMat];
                whiskingMat = [whiskingOnsetMat, whiskingAmplitudeMat, whiskingMidpointMat];
                lickingMat = [lLickMat, rLickMat];
                whiskerTouchMat = [maxDthetaMat, maxDphiMat, maxDkappaHMat, maxDkappaVMat, maxSlideDistanceMat, maxDurationMat, ...    
                    thetaAtTouchMat, phiAtTouchMat, kappaHAtTouchMat, kappaVAtTouchMat, arcLengthAtTouchMat, touchCountMat];
                %%
                %%
                allPredictorsRaw{(cgi-1)*4 + plane} = [whiskerTouchMat, soundMat, drinkMat, whiskingMat, lickingMat];
                nani{(cgi-1)*4 + plane} = find(nanstd(allPredictorsRaw{(cgi-1)*4 + plane})==0);
                allPredictorsMean{(cgi-1)*4 + plane} = nanmean(allPredictorsRaw{(cgi-1)*4 + plane});
                allPredictorsStd{(cgi-1)*4 + plane} = nanstd(allPredictorsRaw{(cgi-1)*4 + plane});
                % normalization of all predictors
                allPredictors{(cgi-1)*4 + plane} = (allPredictorsRaw{(cgi-1)*4 + plane} - nanmean(allPredictorsRaw{(cgi-1)*4 + plane})) ./ nanstd(allPredictorsRaw{(cgi-1)*4 + plane});
                allPredictors{(cgi-1)*4 + plane}(:,nani{(cgi-1)*4 + plane}) = deal(0);
            end
        end

        %%
        whiskerTouchInd = 1 : size(whiskerTouchMat,2);
        soundInd = max(whiskerTouchInd) + 1 : max(whiskerTouchInd) + size(soundMat,2);
        rewardInd = max(soundInd) + 1 : max(soundInd) + size(drinkMat,2);
        whiskingInd = max(rewardInd) + 1 : max(rewardInd) + size(whiskingMat,2);
        lickInd = max(whiskingInd) + 1 : max(whiskingInd) + size(lickingMat,2);

        indPartial{1} = whiskerTouchInd;
        indPartial{2} = soundInd;
        indPartial{3} = rewardInd;
        indPartial{4} = whiskingInd;
        indPartial{5} = lickInd;
        %% settings for saving

        cIDAll = u.cellNums;
        numCell = length(cIDAll); 
        fitCoeffs = cell(numCell,1); % intercept + coefficients of the parameters in training set

        fitDeviance = zeros(numCell,1);
        fitCorrelation = zeros(numCell,1);
        fitCorrPval = zeros(numCell,1);

        fitDevExplained = zeros(numCell,1); % deviance explained from test set
        fitCvDev = zeros(numCell,1); % deviance explained from training set
        fitLambda = zeros(numCell,1);

        cellTime = zeros(numCell,1);
        errorCell = zeros(numCell,2); % 1 for ci, 2 for number of saved files

        tindcellAll = cell(numCell,1);
        cindAll = zeros(numCell,1);
        planeIndAll = zeros(numCell,1);
        for i = 1 : numCell
            tindcellAll{i} = find(cellfun(@(x) ismember(cIDAll(i), x.neuindSession), u.trials));
            cindAll(i) = find(u.trials{tindcellAll{i}(1)}.neuindSession == cIDAll(i));
            planeIndAll(i) = floor(cIDAll(i)/1000);
        end
        spikeAll = cellfun(@(x) x.spk, u.trials, 'uniformoutput', false);

        testTn = cell(numCell,1);
        trainingTn = cell(numCell,1);
        ratioi = zeros(numCell,1);
        ratioInd = zeros(numCell,1);

        %% run it cell by cell
        for ci = 1 : numCell
            tic
            tempPredictor = allPredictors{planeIndAll(ci)};
            tempPlane = mod(planeIndAll(ci),4);
            if tempPlane ==0
                tempPlane = 4;
            end
            tindCell = tindcellAll{ci};
            cind = cindAll(ci);
            
            info.localDir = localDir;
            info.mouse = mouse;
            info.session = session;
            info.ci = ci;
            info.repetition = repetition;
            info.trialNums = u.trialNums;
            info.posShift = posShift;
            info.numFrames = cellfun(@(x) length(x.tpmTime{tempPlane}), u.trials(tindCell));
            info.numCell = numCell;
            info.glmnetOpt = glmnetOpt;
            info.lambdaCV = lambdaCV;
            
            tempSpike = cellfun(@(x) x(cind,:), spikeAll(tindCell), 'un', 0);

            totalNumSpikeFrames = sum(cellfun(@(x) length(find(x)), tempSpike));
            if totalNumSpikeFrames < 10 % from the results of d190415_error_and_nonfit_cells_threshold_finding
                fprintf('Abort because of too little spike frames in cell index %d\n', cellnum)
            else
                flagRun = 0;
                while flagRun < 10
                    flagRun = flagRun + 1;
%                     try
%                         parfun_glmnet_whisker_perCell(info, tempSpike, tempPredictor, stratificationGroups, tindCell); % no particular output since it will save each data into files
                        parfun_glmnet_whisker_perRep(info, allSpikes, allPredictors, stratificationGroups)
%                     catch
%                         myCluster = parcluster('local');
%                         delete(myCluster.Jobs)
%                         clear myCluster
%                         pause(1) % just in case...
%                         parpool(repetition, 'SpmdEnabled', true);
%                     end
                    % IMPORTANT %
                    % This saved file name format should match with that in
                    % parfun_glmnet_whisker_perCell.m
                    savedFnList = dir([sprintf('%sJK%03dS%02dci%04d_save_R',localDir,mouse,session,ci), '*']);
                    if length(savedFnList) >= 10
                        flagRun = 100;
                    end
                end
                savedFnList = dir([sprintf('%sJK%03dS%02dci%04d_save_R',localDir,mouse,session,ci), '*']);
                if length(savedFnList) < 10
                    errorCell(ci,1) = ci;
                    if ~isempty(savedFnList)
                        errorCell(ci,2) = length(savedFnList);                        
                    end
                end
            end
            cellTime(ci) = toc;
        end % end of for ci = 1 : length(cIDAll)

        %% summarize and save the results
        done = find(errorCell(:,1) == 0);
        %%
        for ri = 1 : repetition
            for dci = 1 : length(done)
                ci = done(dci);
                tempFn = sprintf('%sJK%03dS%02dci%04d_save_R%02d',localDir,mouse,session,ci,ri);
                dat = load(tempFn);
                
                fitCoeffs{ci} = dat.fitCoeffs;
                fitDeviance(ci) = dat.fitDeviance;
                fitCorrelation(ci) = dat.fitCorrelation;
                fitCorrPval(ci) = dat.fitCorrPval;
                fitDevExplained(ci) = dat.fitDevExplained;
                fitCvDev(ci) = dat.fitCvDev;
                fitLambda(ci) = dat.fitLambda;
                
                testTn{ci} = dat.testTn;
                trainingTn{ci} = dat.trainingTn;
                ratioi(ci) = dat.ratioi;
                ratioInd(ci) = dat.ratioInd;
            end
            save(sprintf('%s%d\\%s_R%02d',baseDir, mouse, savefnResult, ri), 'fit*', 'allPredictors*', 'indPartial', '*Group', 'testTn', 'trainingTn', 'lambdaCV', '*Opt', 'done', '*Shift', 'cellTime', 'cIDAll', 'ratio*', 'errorCell');
        end
        
    end
end
