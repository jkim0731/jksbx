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
%     - tTouchOnset: total touch onset (parameter)
%     - pTouchOnset: protraction touch onset (parameter)
%     - rTouchOnset: retraction touch onset (parameter)
%     
%     - tTouchFrames: total touch frames (binary)
%     - pTouchFrames: protraction touch frames (binary)
%     - rTouchFrames: retraction touch frames (binary)
%
%     - tTouchDuration: total touch duration within each tpm frame (parameter, in ms)
%     - pTouchDuration: protraction touch duration within each tpm frame (parameter, in ms)
%     - rTouchDuration: retraction touch duration within each tpm frame (parameter, in ms)
%     
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
%     - tLick: total licks (parameter; # of licks in each frame)
%     - lLick: left licks (parameter)
%     - rLick: right licks (parameter)
%     
%     - lickOnset
%     - lickOffset
    
baseDir = 'Y:\Whiskernas\JK\suite2p\';

mice = [25,27,30,36,37,38,39,41,52,53,54,56];
sessions = {[4,19],[3,16],[3,21],[1,17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 

poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isempty(poolobj)    
    parpool('SpmdEnabled', false);
end

% mice = [36,37,38,39,41,52,53,54,56];
% sessions = {[17],[7],[2],[1,22],[3],[3,21],[3],[3],[3]}; 


for mi = 1:length(mice)
% for mi = 1
    for si = 1:length(sessions{mi})
%     for si = 1

        mouse = mice(mi);
        session = sessions{mi}(si);

        posShift = 1;
        negShift = 1;
        testPortion = 0.3; % 30 % test set
        pThresholdNull = 0.05;
        pThresholdPartial = 0.05;
        lickBoutInterval = 1; % licks separated by 1 s regarded as different licking bouts


        glmnetOpt = glmnetSet;
        glmnetOpt.standardize = 0;
        glmnetOpt.alpha = 0.95;
        
        partialGlmOpt = glmnetOpt;
        partialGlmOpt.alpha = 0;
        % lambdaCV = 3; % cross-validation fold number

        
        dn = sprintf('%s%03d',baseDir,mouse);
        ufn = sprintf('UberJK%03dS%02d.mat', mouse, session);
        cd(dn)
%         if exist(ufn,'file')
%             load(ufn)
%         else
            u = Uber.buildUberArray(mouse, session);
%         end
        frameRate = u.frameRate;

        savefn = sprintf('glmResponseType_JK%03dS%02d_glmnet_m6',mouse, session); % m(n) meaining method(n)

        
        
        %% pre-processing for lick onset and offset
        % regardless of licking alternating, each l and r has it's own lick onset and offset. both licking, just take the union

        v = cell(length(u.trials),1);
        % v.bothLickOnset = cell(length(u.trials),1);
        % v.bothLickOffset = cell(length(u.trials),1);
        % v.leftLickOnset = cell(length(u.trials),1);
        % v.leftLickOffset = cell(length(u.trials),1);
        % v.rightLickOnset = cell(length(u.trials),1);
        % v.rightLickOffset = cell(length(u.trials),1);

        for ui = 1 : length(u.trials)
            bothLickTime = union(u.trials{ui}.leftLickTime, u.trials{ui}.rightLickTime);
            if length(bothLickTime) == 1
                v{ui}.bothLickOnset = bothLickTime;
                v{ui}.bothLickOffset = bothLickTime;
            elseif length(bothLickTime) > 1        
                onsets = find(diff(bothLickTime) > lickBoutInterval);
                if isempty(onsets)
                    v{ui}.bothLickOnset = bothLickTime(1);
                    v{ui}.bothLickOffset = bothLickTime(end);
                else
                    v{ui}.bothLickOnset = bothLickTime([1; onsets+1]);
                    v{ui}.bothLickOffset = bothLickTime([onsets; end]);
                end
            else
                v{ui}.bothLickOnset = [];
                v{ui}.bothLickOffset = [];
            end

            if length(u.trials{ui}.leftLickTime) == 1
                v{ui}.leftLickOnset = u.trials{ui}.leftLickTime;
                v{ui}.leftLickOffset = u.trials{ui}.leftLickTime;
            elseif length(u.trials{ui}.leftLickTime) > 1        
                onsets = find(diff(u.trials{ui}.leftLickTime) > lickBoutInterval);
                if isempty(onsets)
                    v{ui}.leftLickOnset = u.trials{ui}.leftLickTime(1);
                    v{ui}.leftLickOffset = u.trials{ui}.leftLickTime(end);
                else
                    v{ui}.leftLickOnset = u.trials{ui}.leftLickTime([1; onsets+1]);
                    v{ui}.leftLickOffset = u.trials{ui}.leftLickTime([onsets; end]);
                end
            else
                v{ui}.leftLickOnset = [];
                v{ui}.leftLickOffset = [];
            end

            if length(u.trials{ui}.rightLickTime) == 1
                v{ui}.rightLickOnset = u.trials{ui}.rightLickTime;
                v{ui}.rightLickOffset = u.trials{ui}.rightLickTime;
            elseif length(u.trials{ui}.rightLickTime) > 1        
                onsets = find(diff(u.trials{ui}.rightLickTime) > lickBoutInterval);
                if isempty(onsets)
                    v{ui}.rightLickOnset = u.trials{ui}.rightLickTime(1);
                    v{ui}.rightLickOffset = u.trials{ui}.rightLickTime(end);
                else
                    v{ui}.rightLickOnset = u.trials{ui}.rightLickTime([1; onsets+1]);
                    v{ui}.rightLickOffset = u.trials{ui}.rightLickTime([onsets; end]);
                end
            else
                v{ui}.rightLickOnset = [];
                v{ui}.rightLickOffset = [];
            end
            v{ui}.tpmTime = u.trials{ui}.tpmTime;
        end
        
        
        %% repetition test

%         division = 20;
        repetition = 1;
        rtest = struct;
        for ri = 1 : repetition
        
            
            
            
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
        for ci = 1:2
        % for ci = 1
            tindcell = find(cellfun(@(x) ismember(1001+(ci-1)*4000, x.neuindSession), u.trials));

            tind = intersect(tindcell, trainingInd);
            for plane = 1 : 4    
        %     for plane = 1
                pTouchOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(cellfun(@(y) y(1), x.protractionTouchChunks), [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
                rTouchOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(cellfun(@(y) y(1), x.retractionTouchChunks), [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
                
                whiskerVideoFrameDuration = u.trials{tind(1)}.frameDuration * 1000; % in ms
                
                pTouchDuration = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.protractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime{plane}]), zeros(1,posShift)], ...
                    u.trials(tind)','uniformoutput',false));
                rTouchDuration = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.retractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime{plane}]), zeros(1,posShift)], ...
                    u.trials(tind)','uniformoutput',false));

                pTouchFrames = pTouchDuration;
                pTouchFrames(pTouchDuration > 0) = 1;
                rTouchFrames = rTouchDuration;
                rTouchFrames(rTouchDuration > 0) = 1;
                
                pTouchOnsetAngles = cell(length(angles),1);
                rTouchOnsetAngles = cell(length(angles),1);
                pTouchDurationAngles = cell(length(angles),1);
                rTouchDurationAngles = cell(length(angles),1);
                pTouchFramesAngles = cell(length(angles),1);
                rTouchFramesAngles = cell(length(angles),1);                
                for ai = 1 : length(angles)
                    tempAngleBinary = cell2mat(cellfun(@(x) ones(length(x.tpmTime{plane}) + 2 * posShift, 1) * (x.angle == angles(ai)), u.trials(tind), 'uniformoutput', false));
                    pTouchOnsetAngles{ai} = pTouchOnset .* tempAngleBinary';
                    rTouchOnsetAngles{ai} = rTouchOnset .* tempAngleBinary';
                    pTouchDurationAngles{ai} = pTouchDuration .* tempAngleBinary';
                    rTouchDurationAngles{ai} = rTouchDuration .* tempAngleBinary';
                    pTouchFramesAngles{ai} = pTouchFrames .* tempAngleBinary';
                    rTouchFramesAngles{ai} = rTouchFrames .* tempAngleBinary';                    
                end
                                
                scPiezo = cell2mat(cellfun(@(x) [zeros(1,posShift), 1, zeros(1,length(x.tpmTime{plane})-1), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
                scPoleup = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.poleUpOnsetTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
                scPoledown = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.poleDownOnsetTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
                drinkOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.drinkingOnsetTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));

                lLick = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.leftLickTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
                rLick = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.rightLickTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));

                lLickOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.leftLickOnset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
                rLickOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.rightLickOnset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));

                lLickOffset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.leftLickOffset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
                rLickOffset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.rightLickOffset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));

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
                    whiskingOnsetCell{ti} = [zeros(1,posShift), histcounts(onsetTimes, time), zeros(1,posShift)];

                    tempAmp = zeros(1,length(time)-1);        
                    tempMid = zeros(1,length(time)-1);
                    for i = 1 : length(tempAmp)
                        startInd = find(wtimes >= time(i), 1, 'first');
                        endInd = find(wtimes < time(i+1), 1, 'last');
                        tempAmp(i) = max(amplitude(startInd:endInd));
                        tempMid(i) = mean(midpoint(startInd:endInd));
                    end
                    whiskingAmpCell{ti} = [zeros(1,posShift), tempAmp, zeros(1,posShift)];        
                    whiskingMidpointCell{ti} = [nan(1,posShift), tempMid, nan(1,posShift)];    
                end
                whiskingOnset = cell2mat(whiskingOnsetCell);
                whiskingAmp = cell2mat(whiskingAmpCell);
                whiskingOA = (whiskingOnset > 0).*whiskingAmp;
                whiskingMidpoint = cell2mat(whiskingMidpointCell);    
                whiskingMidpoint(isnan(whiskingMidpoint)) = deal(mode(whiskingMidpoint(isfinite(whiskingMidpoint))));
                                
                %%
                pTouchOnsetMat = zeros(length(pTouchOnset), (posShift + 1) * length(angles));
                rTouchOnsetMat = zeros(length(rTouchOnset), (posShift + 1) * length(angles));
                pTouchFramesMat = zeros(length(pTouchFrames),(posShift + 1) * length(angles));
                rTouchFramesMat = zeros(length(rTouchFrames), (posShift + 1) * length(angles));
                pTouchDurationMat = zeros(length(pTouchDuration), (posShift + 1) * length(angles));
                rTouchDurationMat = zeros(length(rTouchDuration), (posShift + 1) * length(angles));
                
                scPiezoMat = zeros(length(scPiezo), posShift + 1);
                scPoleUpMat = zeros(length(scPoleup), posShift + 1);
                scPoleDownMat = zeros(length(scPoledown), posShift + 1);
                drinkOnsetMat = zeros(length(drinkOnset), posShift + 1);
                for i = 1 : posShift + 1
                    for ai = 1 : length(angles)                        
                        pTouchOnsetMat(:,(i-1)*length(angles) + ai) = circshift(pTouchOnsetAngles{ai}, [0 i-1])';
                        rTouchOnsetMat(:,(i-1)*length(angles) + ai) = circshift(rTouchOnsetAngles{ai}, [0 i-1])';
                        pTouchFramesMat(:,(i-1)*length(angles) + ai) = circshift(pTouchFramesAngles{ai}, [0 i-1])';
                        rTouchFramesMat(:,(i-1)*length(angles) + ai) = circshift(rTouchFramesAngles{ai}, [0 i-1])';                        
                        pTouchDurationMat(:,(i-1)*length(angles) + ai) = circshift(pTouchDurationAngles{ai}, [0 i-1])';
                        rTouchDurationMat(:,(i-1)*length(angles) + ai) = circshift(rTouchDurationAngles{ai}, [0 i-1])';
                    end
                    scPiezoMat(:,i) = circshift(scPiezo, [0 i-1])';
                    scPoleUpMat(:,i) = circshift(scPoleup, [0 i-1])';
                    scPoleDownMat(:,i) = circshift(scPoledown, [0 i-1])';
                    drinkOnsetMat(:,i) = circshift(drinkOnset, [0 i-1])';
                end

                whiskingOnsetMat = zeros(length(whiskingOnset), negShift + posShift + 1);
                whiskingAmpMat = zeros(length(whiskingAmp), negShift + posShift + 1);
                whiskingOAMat = zeros(length(whiskingOA), negShift + posShift + 1);
                whiskingMidpointMat = zeros(length(whiskingMidpoint), negShift + posShift + 1);
                lLickMat = zeros(length(lLick), negShift + posShift + 1);
                rLickMat = zeros(length(rLick), negShift + posShift + 1);
                lLickOnsetMat = zeros(length(lLickOnset), negShift + posShift + 1);
                rLickOnsetMat = zeros(length(rLickOnset), negShift + posShift + 1);
                lLickOffsetMat = zeros(length(lLickOffset), negShift + posShift + 1);
                rLickOffsetMat = zeros(length(rLickOffset), negShift + posShift + 1);
                for i = 1 : negShift + posShift + 1
                    whiskingOnsetMat(:,i) = circshift(whiskingOnset, [0 -negShift + i - 1])';
                    whiskingAmpMat(:,i) = circshift(whiskingAmp, [0 -negShift + i - 1])';
                    whiskingOAMat(:,i) = circshift(whiskingOA, [0 -negShift + i - 1])';
                    whiskingMidpointMat(:,i) = circshift(whiskingMidpoint, [0 -negShift + i - 1])';
                    lLickMat(:,i) = circshift(lLick, [0 -negShift + i - 1])';
                    rLickMat(:,i) = circshift(rLick, [0 -negShift + i - 1])';
                    lLickOnsetMat(:,i) = circshift(lLickOnset, [0 -negShift + i - 1])';
                    rLickOnsetMat(:,i) = circshift(rLickOnset, [0 -negShift + i - 1])';
                    lLickOffsetMat(:,i) = circshift(lLickOffset, [0 -negShift + i - 1])';
                    rLickOffsetMat(:,i) = circshift(rLickOffset, [0 -negShift + i - 1])';
                end
                tempIM = [pTouchOnsetMat, rTouchOnsetMat, pTouchFramesMat, rTouchFramesMat, pTouchDurationMat, rTouchDurationMat, ...
                    scPiezoMat, scPoleUpMat, scPoleDownMat, ...
                    drinkOnsetMat, ...
                    whiskingOnsetMat, whiskingAmpMat, whiskingOAMat, whiskingMidpointMat, ...
                    lLickMat, rLickMat, lLickOnsetMat, rLickOnsetMat, lLickOffsetMat, rLickOffsetMat];
                trainingInputMat{(ci-1)*4 + plane} = (tempIM - ones(size(tempIM,1),1)*min(tempIM)) ./ (ones(size(tempIM,1),1)*(max(tempIM) - min(tempIM)));
                trainingInputMat{(ci-1)*4 + plane}(isnan(trainingInputMat{(ci-1)*4 + plane})) = deal(0); % the ones with NaN, from all of the entries being 0
                                
            end
        end

        %% design matrix for test 
        testInputMat = cell(8,1);
        for ci = 1:2
        % for ci = 1
            tindcell = find(cellfun(@(x) ismember(1001+(ci-1)*4000, x.neuindSession), u.trials));

            tind = intersect(tindcell, testInd);
            for plane = 1 : 4    
        %     for plane = 1
                pTouchOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(cellfun(@(y) y(1), x.protractionTouchChunks), [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
                rTouchOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(cellfun(@(y) y(1), x.retractionTouchChunks), [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
                
                whiskerVideoFrameDuration = u.trials{tind(1)}.frameDuration * 1000; % in ms
                
                pTouchDuration = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.protractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime{plane}]), zeros(1,posShift)], ...
                    u.trials(tind)','uniformoutput',false));
                rTouchDuration = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.retractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime{plane}]), zeros(1,posShift)], ...
                    u.trials(tind)','uniformoutput',false));

                pTouchFrames = pTouchDuration;
                pTouchFrames(pTouchDuration > 0) = 1;
                rTouchFrames = rTouchDuration;
                rTouchFrames(rTouchDuration > 0) = 1;
                
                pTouchOnsetAngles = cell(length(angles),1);
                rTouchOnsetAngles = cell(length(angles),1);
                pTouchDurationAngles = cell(length(angles),1);
                rTouchDurationAngles = cell(length(angles),1);
                pTouchFramesAngles = cell(length(angles),1);
                rTouchFramesAngles = cell(length(angles),1);                
                for ai = 1 : length(angles)
                    tempAngleBinary = cell2mat(cellfun(@(x) ones(length(x.tpmTime{plane}) + 2 * posShift, 1) * (x.angle == angles(ai)), u.trials(tind), 'uniformoutput', false));
                    pTouchOnsetAngles{ai} = pTouchOnset .* tempAngleBinary';
                    rTouchOnsetAngles{ai} = rTouchOnset .* tempAngleBinary';
                    pTouchDurationAngles{ai} = pTouchDuration .* tempAngleBinary';
                    rTouchDurationAngles{ai} = rTouchDuration .* tempAngleBinary';
                    pTouchFramesAngles{ai} = pTouchFrames .* tempAngleBinary';
                    rTouchFramesAngles{ai} = rTouchFrames .* tempAngleBinary';                    
                end
                                
                scPiezo = cell2mat(cellfun(@(x) [zeros(1,posShift), 1, zeros(1,length(x.tpmTime{plane})-1), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
                scPoleup = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.poleUpOnsetTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
                scPoledown = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.poleDownOnsetTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
                drinkOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.drinkingOnsetTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));

                lLick = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.leftLickTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
                rLick = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.rightLickTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));

                lLickOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.leftLickOnset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
                rLickOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.rightLickOnset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));

                lLickOffset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.leftLickOffset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
                rLickOffset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.rightLickOffset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));

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
                    whiskingOnsetCell{ti} = [zeros(1,posShift), histcounts(onsetTimes, time), zeros(1,posShift)];

                    tempAmp = zeros(1,length(time)-1);        
                    tempMid = zeros(1,length(time)-1);
                    for i = 1 : length(tempAmp)
                        startInd = find(wtimes >= time(i), 1, 'first');
                        endInd = find(wtimes < time(i+1), 1, 'last');
                        tempAmp(i) = max(amplitude(startInd:endInd));            
                        tempMid(i) = mean(midpoint(startInd:endInd));
                    end
                    whiskingAmpCell{ti} = [zeros(1,posShift), tempAmp, zeros(1,posShift)];        
                    whiskingMidpointCell{ti} = [nan(1,posShift), tempMid, nan(1,posShift)];    
                end
                whiskingOnset = cell2mat(whiskingOnsetCell);
                whiskingAmp = cell2mat(whiskingAmpCell);
                whiskingOA = (whiskingOnset > 0).*whiskingAmp;
                whiskingMidpoint = cell2mat(whiskingMidpointCell);    
                whiskingMidpoint(isnan(whiskingMidpoint)) = deal(mode(whiskingMidpoint(isfinite(whiskingMidpoint))));
                                
                %%
                pTouchOnsetMat = zeros(length(pTouchOnset), (posShift + 1) * length(angles));
                rTouchOnsetMat = zeros(length(rTouchOnset), (posShift + 1) * length(angles));
                pTouchFramesMat = zeros(length(pTouchFrames),(posShift + 1) * length(angles));
                rTouchFramesMat = zeros(length(rTouchFrames), (posShift + 1) * length(angles));
                pTouchDurationMat = zeros(length(pTouchDuration), (posShift + 1) * length(angles));
                rTouchDurationMat = zeros(length(rTouchDuration), (posShift + 1) * length(angles));
                
                scPiezoMat = zeros(length(scPiezo), posShift + 1);
                scPoleUpMat = zeros(length(scPoleup), posShift + 1);
                scPoleDownMat = zeros(length(scPoledown), posShift + 1);
                drinkOnsetMat = zeros(length(drinkOnset), posShift + 1);
                for i = 1 : posShift + 1
                    for ai = 1 : length(angles)                        
                        pTouchOnsetMat(:,(i-1)*length(angles) + ai) = circshift(pTouchOnsetAngles{ai}, [0 i-1])';
                        rTouchOnsetMat(:,(i-1)*length(angles) + ai) = circshift(rTouchOnsetAngles{ai}, [0 i-1])';
                        pTouchFramesMat(:,(i-1)*length(angles) + ai) = circshift(pTouchFramesAngles{ai}, [0 i-1])';
                        rTouchFramesMat(:,(i-1)*length(angles) + ai) = circshift(rTouchFramesAngles{ai}, [0 i-1])';                        
                        pTouchDurationMat(:,(i-1)*length(angles) + ai) = circshift(pTouchDurationAngles{ai}, [0 i-1])';
                        rTouchDurationMat(:,(i-1)*length(angles) + ai) = circshift(rTouchDurationAngles{ai}, [0 i-1])';
                    end
                    scPiezoMat(:,i) = circshift(scPiezo, [0 i-1])';
                    scPoleUpMat(:,i) = circshift(scPoleup, [0 i-1])';
                    scPoleDownMat(:,i) = circshift(scPoledown, [0 i-1])';
                    drinkOnsetMat(:,i) = circshift(drinkOnset, [0 i-1])';
                end

                whiskingOnsetMat = zeros(length(whiskingOnset), negShift + posShift + 1);
                whiskingAmpMat = zeros(length(whiskingAmp), negShift + posShift + 1);
                whiskingOAMat = zeros(length(whiskingOA), negShift + posShift + 1);
                whiskingMidpointMat = zeros(length(whiskingMidpoint), negShift + posShift + 1);
                lLickMat = zeros(length(lLick), negShift + posShift + 1);
                rLickMat = zeros(length(rLick), negShift + posShift + 1);
                lLickOnsetMat = zeros(length(lLickOnset), negShift + posShift + 1);
                rLickOnsetMat = zeros(length(rLickOnset), negShift + posShift + 1);
                lLickOffsetMat = zeros(length(lLickOffset), negShift + posShift + 1);
                rLickOffsetMat = zeros(length(rLickOffset), negShift + posShift + 1);
                for i = 1 : negShift + posShift + 1
                    whiskingOnsetMat(:,i) = circshift(whiskingOnset, [0 -negShift + i - 1])';
                    whiskingAmpMat(:,i) = circshift(whiskingAmp, [0 -negShift + i - 1])';
                    whiskingOAMat(:,i) = circshift(whiskingOA, [0 -negShift + i - 1])';
                    whiskingMidpointMat(:,i) = circshift(whiskingMidpoint, [0 -negShift + i - 1])';
                    lLickMat(:,i) = circshift(lLick, [0 -negShift + i - 1])';
                    rLickMat(:,i) = circshift(rLick, [0 -negShift + i - 1])';
                    lLickOnsetMat(:,i) = circshift(lLickOnset, [0 -negShift + i - 1])';
                    rLickOnsetMat(:,i) = circshift(rLickOnset, [0 -negShift + i - 1])';
                    lLickOffsetMat(:,i) = circshift(lLickOffset, [0 -negShift + i - 1])';
                    rLickOffsetMat(:,i) = circshift(rLickOffset, [0 -negShift + i - 1])';
                end
                tempIM = [pTouchOnsetMat, rTouchOnsetMat, pTouchFramesMat, rTouchFramesMat, pTouchDurationMat, rTouchDurationMat, ...
                    scPiezoMat, scPoleUpMat, scPoleDownMat, ...
                    drinkOnsetMat, ...
                    whiskingOnsetMat, whiskingAmpMat, whiskingOAMat, whiskingMidpointMat, ...
                    lLickMat, rLickMat, lLickOnsetMat, rLickOnsetMat, lLickOffsetMat, rLickOffsetMat];
                testInputMat{(ci-1)*4 + plane} = (tempIM - ones(size(tempIM,1),1)*min(tempIM)) ./ (ones(size(tempIM,1),1)*(max(tempIM) - min(tempIM)));
                testInputMat{(ci-1)*4 + plane}(isnan(testInputMat{(ci-1)*4 + plane})) = deal(0); % the ones with NaN, from all of the entries being 0
            end
        end

        %%
        touchInd = 1 : 2 * (posShift + 1) * (1 + length(angles));
        soundInd = max(touchInd) + 1 : max(touchInd) + 3 * (posShift + 1);
        rewardInd = max(soundInd) + 1 : max(soundInd) + posShift + 1;
        whiskingInd = max(rewardInd) + 1 : max(rewardInd) + 4 * (negShift + posShift + 1);
        lickInd = max(whiskingInd) + 1 : max(whiskingInd) + 2 * 3 * (negShift + posShift + 1);

        indPartial{1} = touchInd;
        indPartial{2} = soundInd;
        indPartial{3} = rewardInd;
        indPartial{4} = whiskingInd;
        indPartial{5} = lickInd;

        %%
%         rtest(ri).fitInd = cell(length(u.cellNums),1); % parameters surviving lasso in training set
%         rtest(ri).fitCoeffs = nan(length(u.cellNums),6); % first column is dummy
%         
%         rtest(ri).fitResults = zeros(length(u.cellNums), 6);        
%         % fitResult(:,1) if full fitting is significant (compared to null model), 0 if not
%         % fitResult(:,2) for touchInd, compared to full fitting. if excluding touch
%         % is significantly less fit, then 1, 0 otherwise
%         % fitResult(:,3) for sound, (:,4) for reward, (:,5) for whisking, and (:,6) for licking
% 
%         rtest(ri).devExplained = zeros(length(u.cellNums),1);

        
        fitInd = cell(length(u.cellNums),1); % parameters surviving lasso in training set
        fitCoeffs = nan(length(u.cellNums),6); % first column is a dummy
        
        fitResults = zeros(length(u.cellNums), 6);        
        devExplained = zeros(length(u.cellNums),1);
                
        parfor cellnum = 1 : length(u.cellNums)
%         for cellnum = 1
%         ci = 0;
%         for cellnum = 1:division:length(u.cellNums)
%             ci = ci + 1;
%         for cellnum = 1
        %     cellnum = 1;
            fitCoeff = zeros(1,6);
            fprintf('Loop %d: Running cell %d/%d \n', ri,cellnum, length(u.cellNums));
            cID = u.cellNums(cellnum);

            % find out trial indices for this specific cell
            tindcell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));

            tind = intersect(tindcell, trainingInd);
            % find out row number of this cell
            cind = find(u.trials{tind(1)}.neuindSession == cID);
            planeInd = floor(cID/1000);

            spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cind,:), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
            %%
            finiteInd = find(isfinite(spkTrain));
            input = trainingInputMat{planeInd}(finiteInd,:);
            spk = spkTrain(finiteInd)';

            cv = cvglmnet(input, spk, 'poisson', glmnetOpt);
        %     [BFull, FitInfoFull] = lassoglm(trainingInputMat{planeInd}, spkTrain, 'poisson', 'Lambda', logspace(-4,2), 'CV', lambdaCV, 'Options', opt);
            %% survived coefficients
            iLambda = find(cv.lambda == cv.lambda_1se);
            coeffs = find(cv.glmnet_fit.beta(:,iLambda));
%             rtest(ri).fitInd{cellnum} = coeffs;
            fitInd{cellnum} = coeffs;
            for i = 1 : length(indPartial)
                if sum(ismember(indPartial{i},coeffs)>0)
                    fitCoeff(i + 1) = 1;
                else
                    fitCoeff(i + 1) = 0;
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
        %     spkTest = spkTest';
        %     finiteInd = find(isfinite(spkTest));
        %     spkTest = spkTest(finiteInd)';
            %% (1) if the full model is significant
            fitResult = zeros(1,6);

            mu = nanmean(spkTest); % null poisson parameter
            nullLogLikelihood = nansum(log(poisspdf(spkTest,mu)));
            fullLogLikelihood = nansum(log(poisspdf(spkTest',exp([ones(size(testInputMat{planeInd},1),1),testInputMat{planeInd}]*[cv.glmnet_fit.a0(iLambda); cv.glmnet_fit.beta(:,iLambda)]))));
            saturatedLogLikelihood = nansum(log(poisspdf(spkTest,spkTest)));
            devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
            dfFullNull = length(coeffs);
            devExplained(cellnum) = 1 - (saturatedLogLikelihood - fullLogLikelihood)/(saturatedLogLikelihood - nullLogLikelihood);
            if devianceFullNull > chi2inv(1-pThresholdNull, dfFullNull)
                fitResult(1) = 1;
                %% (2) test without each parameter (as a group)                
                for pi = 1 : 5
                    if find(ismember(coeffs, indPartial{pi}))
                        if all(ismember(coeffs, indPartial{pi}))
                            fitResult(pi+1) = 1;
                        else
                            tempTrainInput = trainingInputMat{planeInd}(:,setdiff(coeffs,indPartial{pi}));
                            tempTestInput = testInputMat{planeInd}(:,setdiff(coeffs,indPartial{pi}));
                            cvPartial = cvglmnet(tempTrainInput(finiteInd,:), spk, 'poisson', partialGlmOpt);
                            iLambda = find(cvPartial.lambda == cvPartial.lambda_1se);
                            partialLogLikelihood = nansum(log(poisspdf(spkTest', exp([ones(size(tempTestInput,1),1), tempTestInput] * [cvPartial.glmnet_fit.a0(iLambda); cvPartial.glmnet_fit.beta(:,iLambda)]))));
                            devianceFullPartial = 2*(fullLogLikelihood - partialLogLikelihood);
                            dfFullPartial = dfFullNull - cvPartial.glmnet_fit.df(iLambda);
                            if devianceFullPartial > chi2inv(1-pThresholdPartial, dfFullPartial)
                                fitResult(pi+1) = 1;
                            end
                        end
                    end
                end
            end
            fitResults(cellnum,:) = fitResult;
            fitCoeffs(cellnum,:) = fitCoeff;
        end % end of parfor cellnum
        
        rtest(ri).fitInd = fitInd; % parameters surviving lasso in training set
        rtest(ri).fitCoeffs = fitCoeffs; % first column is dummy
        
        rtest(ri).fitResults = fitResults;        
        % fitResult(:,1) if full fitting is significant (compared to null model), 0 if not
        % fitResult(:,2) for touchInd, compared to full fitting. if excluding touch
        % is significantly less fit, then 1, 0 otherwise
        % fitResult(:,3) for sound, (:,4) for reward, (:,5) for whisking, and (:,6) for licking

        rtest(ri).devExplained = devExplained;
        
        save(savefn, 'rtest');
%         save('poisson_binary_rtest.m', 'rtest')
        end % of ri. random group selection index
    end
end



