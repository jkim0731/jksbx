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
%     - scPiezo: piezo sound cue onset (binary)
%     - scPoleup: pole up sound cue onset (binary)
%     - scPoledown: pole down sound cue onset (binary)
%     
%     - drinkOnset: drinking onset (binary)
%     
%     
%     % motor variables: shift both backward and forward
%     - whiskingOnset: whisking onset (parameter; # of onset in each frame)
%     - whiskingAmp: whisking amplitude (parameter; from whisker decomposition)
%     - whiskingMidpoint: whisking midpoint(parameter; from whisker decomposition)
%
%     - tLick: total licks (parameter; # of licks in each frame)
%     - lLick: left licks (parameter)
%     - rLick: right licks (parameter)
%     
    
baseDir = 'Y:\Whiskernas\JK\suite2p\';

mouse = 25;
session = 4;


lambdaCV = 3; % cross-validation fold number
posShift = 4;
negShift = 2;
testPortion = 0.3; % 30 % test set
pThresholdNull = 0.05;
pThresholdPartial = 0.05;
lickBoutInterval = 1; % licks separated by 1 s regarded as different licking bouts
opt = statset('UseParallel', true);

dn = sprintf('%s%03d',baseDir,mouse);
ufn = sprintf('UberJK%03dS%02d.mat', mouse, session);
cd(dn)
load(ufn)
% u = Uber.buildUberArray(mouse, session);
frameRate = u.frameRate;

savefn = sprintf('glmResponseType_JK%03dS%02d',mouse, session);

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
        tTouchOnset = pTouchOnset + rTouchOnset;

        whiskerVideoFrameDuration = u.trials{tind(1)}.frameDuration * 1000; % in ms
        pTouchDuration = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.protractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime{plane}]), zeros(1,posShift)], ...
            u.trials(tind)','uniformoutput',false));
        rTouchDuration = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.retractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime{plane}]), zeros(1,posShift)], ...
            u.trials(tind)','uniformoutput',false));
        tTouchDuration = pTouchDuration + rTouchDuration;

        pTouchFrames = pTouchDuration;
        pTouchFrames(pTouchDuration > 0) = 1;
        rTouchFrames = rTouchDuration;
        rTouchFrames(rTouchDuration > 0) = 1;
        tTouchFrames = tTouchDuration;
        tTouchFrames(tTouchDuration > 0) = 1;

        scPiezo = cell2mat(cellfun(@(x) [zeros(1,posShift), 1, zeros(1,length(x.tpmTime{plane})-1), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
        scPoleup = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.poleUpOnsetTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
        scPoledown = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.poleDownOnsetTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
        drinkOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.drinkingOnsetTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));

        lLick = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.leftLickTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
        rLick = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.rightLickTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
        tLick = lLick + rLick;
        tLick(tLick>0) = 1;
        
        lLickOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.leftLickOnset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
        rLickOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.rightLickOnset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
        tLickOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.bothLickOnset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
        
        lLickOffset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.leftLickOffset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
        rLickOffset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.rightLickOffset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
        tLickOffset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.bothLickOffset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));

        %%
        whiskingOnsetCell = cell(1,length(tind));
        whiskingAmpCell = cell(1,length(tind));
        whiskingMaxAmpCell = cell(1,length(tind));
        whiskingMidpointCell = cell(1,length(tind));

        for ti = 1 : length(tind)
            currTrial = u.trials{tind(ti)};
            time = [0, currTrial.tpmTime{plane}];
            wtimes = [0:currTrial.nof-1] * currTrial.frameDuration;
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
        % corr(spk(isfinite(spk))', tTouchOnset(isfinite(tTouchOnset))')
        % corr(spk(isfinite(spk))', pTouchOnset(isfinite(pTouchOnset))')
        % corr(spk(isfinite(spk))', rTouchOnset(isfinite(rTouchOnset))')
        % corr(spk(isfinite(spk))', tTouchFrames(isfinite(tTouchFrames))')
        % corr(spk(isfinite(spk))', pTouchFrames(isfinite(pTouchFrames))')
        % corr(spk(isfinite(spk))', rTouchFrames(isfinite(rTouchFrames))')
        % corr(spk(isfinite(spk))', tTouchDuration(isfinite(tTouchDuration))')
        % corr(spk(isfinite(spk))', pTouchDuration(isfinite(pTouchDuration))')
        % corr(spk(isfinite(spk))', rTouchDuration(isfinite(rTouchDuration))')
        % corr(spk(isfinite(spk))', scPiezo(isfinite(scPiezo))')
        % corr(spk(isfinite(spk))', scPoleup(isfinite(scPoleup))')
        % corr(spk(isfinite(spk))', scPoledown(isfinite(scPoledown))')
        % corr(spk(isfinite(spk))', drinkOnset(isfinite(drinkOnset))')
        % corr(spk(isfinite(spk))', whiskingOnset(isfinite(whiskingOnset))')
        % corr(spk(isfinite(spk))', whiskingAmp(isfinite(whiskingAmp))')
        % corr(spk(isfinite(spk))', whiskingMidpoint(isfinite(whiskingMidpoint))')
        % corr(spk(isfinite(spk))', tLick(isfinite(tLick))')
        % corr(spk(isfinite(spk))', lLick(isfinite(lLick))')
        % corr(spk(isfinite(spk))', rLick(isfinite(rLick))')

        %%
        % [b, dev, stat] = glmfit([tTouchOnset', pTouchOnset', rTouchOnset', tTouchFrames', pTouchFrames', rTouchFrames', ...
        %     tTouchDuration', pTouchDuration', rTouchDuration', scPiezo', scPoleup', scPoledown', drinkOnset', ...
        %     whiskingOnset', whiskingAmp', whiskingMidpoint', tLick', lLick', rLick'], spk', 'poisson');
        
        %%
        tTouchOnsetMat = zeros(length(tTouchOnset), posShift + 1);
        pTouchOnsetMat = zeros(length(pTouchOnset), posShift + 1);
        rTouchOnsetMat = zeros(length(rTouchOnset), posShift + 1);
        tTouchFramesMat = zeros(length(tTouchFrames), posShift + 1);
        pTouchFramesMat = zeros(length(pTouchFrames), posShift + 1);
        rTouchFramesMat = zeros(length(rTouchFrames), posShift + 1);
        tTouchDurationMat = zeros(length(tTouchDuration), posShift + 1);
        pTouchDurationMat = zeros(length(pTouchDuration), posShift + 1);
        rTouchDurationMat = zeros(length(rTouchDuration), posShift + 1);
        scPiezoMat = zeros(length(scPiezo), posShift + 1);
        scPoleUpMat = zeros(length(scPoleup), posShift + 1);
        scPoleDownMat = zeros(length(scPoledown), posShift + 1);
        drinkOnsetMat = zeros(length(drinkOnset), posShift + 1);
        for i = 1 : posShift + 1
            tTouchOnsetMat(:,i) = circshift(tTouchOnset, [0 i-1])';
            pTouchOnsetMat(:,i) = circshift(pTouchOnset, [0 i-1])';
            rTouchOnsetMat(:,i) = circshift(rTouchOnset, [0 i-1])';
            tTouchFramesMat(:,i) = circshift(tTouchFrames, [0 i-1])';
            pTouchFramesMat(:,i) = circshift(pTouchFrames, [0 i-1])';
            rTouchFramesMat(:,i) = circshift(rTouchFrames, [0 i-1])';
            tTouchDurationMat(:,i) = circshift(tTouchDuration, [0 i-1])';
            pTouchDurationMat(:,i) = circshift(pTouchDuration, [0 i-1])';
            rTouchDurationMat(:,i) = circshift(rTouchDuration, [0 i-1])';
            scPiezoMat(:,i) = circshift(scPiezo, [0 i-1])';
            scPoleUpMat(:,i) = circshift(scPoleup, [0 i-1])';
            scPoleDownMat(:,i) = circshift(scPoledown, [0 i-1])';
            drinkOnsetMat(:,i) = circshift(drinkOnset, [0 i-1])';
        end
        
        whiskingOnsetMat = zeros(length(whiskingOnset), negShift + posShift + 1);
        whiskingAmpMat = zeros(length(whiskingAmp), negShift + posShift + 1);
        whiskingOAMat = zeros(length(whiskingOA), negShift + posShift + 1);
        whiskingMidpointMat = zeros(length(whiskingMidpoint), negShift + posShift + 1);
        tLickMat = zeros(length(tLick), negShift + posShift + 1);
        lLickMat = zeros(length(lLick), negShift + posShift + 1);
        rLickMat = zeros(length(rLick), negShift + posShift + 1);
        tLickOnsetMat = zeros(length(tLickOnset), negShift + posShift + 1);
        lLickOnsetMat = zeros(length(lLickOnset), negShift + posShift + 1);
        rLickOnsetMat = zeros(length(rLickOnset), negShift + posShift + 1);
        tLickOffsetMat = zeros(length(tLickOffset), negShift + posShift + 1);
        lLickOffsetMat = zeros(length(lLickOffset), negShift + posShift + 1);
        rLickOffsetMat = zeros(length(rLickOffset), negShift + posShift + 1);
        for i = 1 : negShift + posShift + 1
            whiskingOnsetMat(:,i) = circshift(whiskingOnset, [0 -negShift + i - 1])';
            whiskingAmpMat(:,i) = circshift(whiskingAmp, [0 -negShift + i - 1])';
            whiskingOAMat(:,i) = circshift(whiskingOA, [0 -negShift + i - 1])';
            whiskingMidpointMat(:,i) = circshift(whiskingMidpoint, [0 -negShift + i - 1])';
            tLickMat(:,i) = circshift(tLick, [0 -negShift + i - 1])';
            lLickMat(:,i) = circshift(lLick, [0 -negShift + i - 1])';
            rLickMat(:,i) = circshift(rLick, [0 -negShift + i - 1])';
            tLickOnsetMat(:,i) = circshift(tLickOnset, [0 -negShift + i - 1])';
            lLickOnsetMat(:,i) = circshift(lLickOnset, [0 -negShift + i - 1])';
            rLickOnsetMat(:,i) = circshift(rLickOnset, [0 -negShift + i - 1])';
            tLickOffsetMat(:,i) = circshift(tLickOffset, [0 -negShift + i - 1])';
            lLickOffsetMat(:,i) = circshift(lLickOffset, [0 -negShift + i - 1])';
            rLickOffsetMat(:,i) = circshift(rLickOffset, [0 -negShift + i - 1])';
        end
        tempIM = [tTouchOnsetMat, pTouchOnsetMat, rTouchOnsetMat, tTouchFramesMat, pTouchFramesMat, rTouchFramesMat, ...
            tTouchDurationMat, pTouchDurationMat, rTouchDurationMat, scPiezoMat, scPoleUpMat, scPoleDownMat, ...
            drinkOnsetMat, ...
            whiskingOnsetMat, whiskingAmpMat, whiskingOAMat, whiskingMidpointMat, ...
            tLickMat, lLickMat, rLickMat, tLickOnsetMat, lLickOnsetMat, rLickOnsetMat, tLickOffsetMat, lLickOffsetMat, rLickOffsetMat];
        trainingInputMat{(ci-1)*4 + plane} = (tempIM - min(tempIM)) ./ (max(tempIM) - min(tempIM));
        trainingInputMat{(ci-1)*4 + plane}(isnan(trainingInputMat{(ci-1)*4 + plane})) = deal(0); % the ones with NaN, from all of the entries being 0
        
    end
end



%% correlation between predictors, and with spike


cellnum = 400

cID = u.cellNums(cellnum);

% find out trial indices for this specific cell
tindcell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));

tind = intersect(tindcell, trainingInd);
% find out row number of this cell
cind = find(u.trials{tind(1)}.neuindSession == cID);
planeInd = floor(cID/1000);

spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cind,:), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
input = trainingInputMat{planeInd};

corrIn = [input, spkTrain'];
corrIn = corrIn(isfinite(sum(corrIn,2)),:);

corrMat = corrcoef(corrIn);
figure, imagesc(corrMat)


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
        tTouchOnset = pTouchOnset + rTouchOnset;

        whiskerVideoFrameDuration = u.trials{tind(1)}.frameDuration * 1000; % in ms
        
        pTouchDuration = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.protractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime{plane}]), zeros(1,posShift)], ...
            u.trials(tind)','uniformoutput',false));
        if ~isempty(find(isnan(pTouchDuration)))
                error('nan pTouchDruation')
        end
        rTouchDuration = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(cell2mat(cellfun(@(y) y', x.retractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime{plane}]), zeros(1,posShift)], ...
            u.trials(tind)','uniformoutput',false));
        if ~isempty(find(isnan(pTouchDuration)))
            error('nan rTouchDruation')
        end
        tTouchDuration = pTouchDuration + rTouchDuration;

        pTouchFrames = pTouchDuration;
        pTouchFrames(pTouchDuration > 0) = 1;
        rTouchFrames = rTouchDuration;
        rTouchFrames(rTouchDuration > 0) = 1;
        tTouchFrames = tTouchDuration;
        tTouchFrames(tTouchDuration > 0) = 1;

        scPiezo = cell2mat(cellfun(@(x) [zeros(1,posShift), 1, zeros(1,length(x.tpmTime{plane})-1), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
        scPoleup = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.poleUpOnsetTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
        scPoledown = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.poleDownOnsetTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
        drinkOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.drinkingOnsetTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));

        lLick = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.leftLickTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
        rLick = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.rightLickTime, [0, x.tpmTime{plane}]), zeros(1,posShift)], u.trials(tind)','uniformoutput',false));
        tLick = lLick + rLick;
        tLick(tLick>0) = 1;        
        
        lLickOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.leftLickOnset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
        rLickOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.rightLickOnset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
        tLickOnset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.bothLickOnset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
        
        lLickOffset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.leftLickOffset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
        rLickOffset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.rightLickOffset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
        tLickOffset = cell2mat(cellfun(@(x) [zeros(1,posShift), histcounts(x.bothLickOffset, [0, x.tpmTime{plane}]), zeros(1,posShift)], v(tind)','uniformoutput',false));
        %%
        whiskingOnsetCell = cell(1,length(tind));
        whiskingAmpCell = cell(1,length(tind));
        whiskingMaxAmpCell = cell(1,length(tind));
        whiskingMidpointCell = cell(1,length(tind));

        for ti = 1 : length(tind)
            currTrial = u.trials{tind(ti)};
            time = [0, currTrial.tpmTime{plane}];
            wtimes = [0:currTrial.nof-1] * currTrial.frameDuration;
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
        tTouchOnsetMat = zeros(length(tTouchOnset), posShift + 1);
        pTouchOnsetMat = zeros(length(pTouchOnset), posShift + 1);
        rTouchOnsetMat = zeros(length(rTouchOnset), posShift + 1);
        tTouchFramesMat = zeros(length(tTouchFrames), posShift + 1);
        pTouchFramesMat = zeros(length(pTouchFrames), posShift + 1);
        rTouchFramesMat = zeros(length(rTouchFrames), posShift + 1);
        tTouchDurationMat = zeros(length(tTouchDuration), posShift + 1);
        pTouchDurationMat = zeros(length(pTouchDuration), posShift + 1);
        rTouchDurationMat = zeros(length(rTouchDuration), posShift + 1);
        scPiezoMat = zeros(length(scPiezo), posShift + 1);
        scPoleUpMat = zeros(length(scPoleup), posShift + 1);
        scPoleDownMat = zeros(length(scPoledown), posShift + 1);
        drinkOnsetMat = zeros(length(drinkOnset), posShift + 1);
        for i = 1 : posShift + 1
            tTouchOnsetMat(:,i) = circshift(tTouchOnset, [0 i-1])';
            pTouchOnsetMat(:,i) = circshift(pTouchOnset, [0 i-1])';
            rTouchOnsetMat(:,i) = circshift(rTouchOnset, [0 i-1])';
            tTouchFramesMat(:,i) = circshift(tTouchFrames, [0 i-1])';
            pTouchFramesMat(:,i) = circshift(pTouchFrames, [0 i-1])';
            rTouchFramesMat(:,i) = circshift(rTouchFrames, [0 i-1])';
            tTouchDurationMat(:,i) = circshift(tTouchDuration, [0 i-1])';
            pTouchDurationMat(:,i) = circshift(pTouchDuration, [0 i-1])';
            rTouchDurationMat(:,i) = circshift(rTouchDuration, [0 i-1])';
            scPiezoMat(:,i) = circshift(scPiezo, [0 i-1])';
            scPoleUpMat(:,i) = circshift(scPoleup, [0 i-1])';
            scPoleDownMat(:,i) = circshift(scPoledown, [0 i-1])';
            drinkOnsetMat(:,i) = circshift(drinkOnset, [0 i-1])';
        end
        
        whiskingOnsetMat = zeros(length(whiskingOnset), negShift + posShift + 1);
        whiskingAmpMat = zeros(length(whiskingAmp), negShift + posShift + 1);
        whiskingOAMat = zeros(length(whiskingOA), negShift + posShift + 1);
        whiskingMidpointMat = zeros(length(whiskingMidpoint), negShift + posShift + 1);
        tLickMat = zeros(length(tLick), negShift + posShift + 1);
        lLickMat = zeros(length(lLick), negShift + posShift + 1);
        rLickMat = zeros(length(rLick), negShift + posShift + 1);
        tLickOnsetMat = zeros(length(tLickOnset), negShift + posShift + 1);
        lLickOnsetMat = zeros(length(lLickOnset), negShift + posShift + 1);
        rLickOnsetMat = zeros(length(rLickOnset), negShift + posShift + 1);
        tLickOffsetMat = zeros(length(tLickOffset), negShift + posShift + 1);
        lLickOffsetMat = zeros(length(lLickOffset), negShift + posShift + 1);
        rLickOffsetMat = zeros(length(rLickOffset), negShift + posShift + 1);
        for i = 1 : negShift + posShift + 1
            whiskingOnsetMat(:,i) = circshift(whiskingOnset, [0 -negShift + i - 1])';
            whiskingAmpMat(:,i) = circshift(whiskingAmp, [0 -negShift + i - 1])';
            whiskingOAMat(:,i) = circshift(whiskingOA, [0 -negShift + i - 1])';
            whiskingMidpointMat(:,i) = circshift(whiskingMidpoint, [0 -negShift + i - 1])';
            tLickMat(:,i) = circshift(tLick, [0 -negShift + i - 1])';
            lLickMat(:,i) = circshift(lLick, [0 -negShift + i - 1])';
            rLickMat(:,i) = circshift(rLick, [0 -negShift + i - 1])';
            tLickOnsetMat(:,i) = circshift(tLickOnset, [0 -negShift + i - 1])';
            lLickOnsetMat(:,i) = circshift(lLickOnset, [0 -negShift + i - 1])';
            rLickOnsetMat(:,i) = circshift(rLickOnset, [0 -negShift + i - 1])';
            tLickOffsetMat(:,i) = circshift(tLickOffset, [0 -negShift + i - 1])';
            lLickOffsetMat(:,i) = circshift(lLickOffset, [0 -negShift + i - 1])';
            rLickOffsetMat(:,i) = circshift(rLickOffset, [0 -negShift + i - 1])';
        end
        tempIM = [tTouchOnsetMat, pTouchOnsetMat, rTouchOnsetMat, tTouchFramesMat, pTouchFramesMat, rTouchFramesMat, ...
            tTouchDurationMat, pTouchDurationMat, rTouchDurationMat, scPiezoMat, scPoleUpMat, scPoleDownMat, ...
            drinkOnsetMat, ...
            whiskingOnsetMat, whiskingAmpMat, whiskingOAMat, whiskingMidpointMat, ...
            tLickMat, lLickMat, rLickMat, tLickOnsetMat, lLickOnsetMat, rLickOnsetMat, tLickOffsetMat, lLickOffsetMat, rLickOffsetMat];
        testInputMat{(ci-1)*4 + plane} = (tempIM - min(tempIM)) ./ (max(tempIM) - min(tempIM));
        testInputMat{(ci-1)*4 + plane}(isnan(testInputMat{(ci-1)*4 + plane})) = deal(0); % the ones with NaN, from all of the entries being 0
    end
end

%%
touchInd = 1:3 * 3 * (posShift + 1);
soundInd = max(touchInd) + 1 : max(touchInd) + 3 * (posShift + 1);
rewardInd = max(soundInd) + 1 : max(soundInd) + posShift + 1;
whiskingInd = max(rewardInd) + 1 : max(rewardInd) + 4 * (negShift + posShift + 1);
lickInd = max(whiskingInd) + 1 : max(whiskingInd) + 3 * 3 * (negShift + posShift + 1);

%%
fitResults = zeros(length(u.cellNums), 6);
% fitResult(:,1) if full fitting is significant (compared to null model), 0 if not
% fitResult(:,2) for touchInd, compared to full fitting. if excluding touch
% is significantly less fit, then 1, 0 otherwise
% fitResult(:,3) for sound, (:,4) for reward, (:,5) for whisking, and (:,6) for licking

fitInd = cell(length(u.cellNums),1);

% parameters surviving lasso in training set

for cellnum = 1 : length(u.cellNums)
% for cellnum = 1
%     cellnum = 1;    
    fprintf('Running cell %d/%d \n', cellnum, length(u.cellNums));
    cID = u.cellNums(cellnum);
    
    % find out trial indices for this specific cell
    tindcell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
    
    tind = intersect(tindcell, trainingInd);
    % find out row number of this cell
    cind = find(u.trials{tind(1)}.neuindSession == cID);
    planeInd = floor(cID/1000);

    spkTrain = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cind,:), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
    %%
    [BFull, FitInfoFull] = lassoglm(trainingInputMat{planeInd}, spkTrain, 'poisson', 'Lambda', logspace(-4,2), 'CV', lambdaCV, 'Options', opt);
    %% survived coefficients
%     lassoPlot(B,FitInfoFull,'plottype','CV'); 
%     legend('show') % Show legend
%     idxLambdaMinDeviance = FitInfoFull.IndexMinDeviance;
    idxLambda1SE = FitInfoFull.Index1SE;
    mincoefs = find(BFull(:,idxLambda1SE));
    fitInd{cellnum} = mincoefs;

%     find(ismember(touchInd, mincoefs))
%     find(ismember(soundInd, mincoefs))
%     find(ismember(rewardInd, mincoefs))
%     find(ismember(whiskingInd, mincoefs))
%     find(ismember(lickInd, mincoefs))
        
    %% test
    cID = u.cellNums(cellnum);
    
    % find out trial indices for this specific cell
    tindcell = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
    
    tind = intersect(tindcell, testInd);
    % find out row number of this cell
    cind = find(u.trials{tind(1)}.neuindSession == cID);
    planeInd = floor(cID/1000);

    spkTest = cell2mat(cellfun(@(x) [nan(1,posShift), x.spk(cind,:), nan(1,posShift)], u.trials(tind)','uniformoutput',false));
    
    %% (1) if the full model is significant
    fitResult = zeros(1,6);
    
    mu = nanmean(spkTest); % null poisson parameter
    nullLogLikelihood = nansum(log(poisspdf(spkTest,mu)));
    fullLogLikelihood = nansum(log(poisspdf(spkTest',exp([ones(size(testInputMat{planeInd},1),1),testInputMat{planeInd}]*[FitInfoFull.Intercept(idxLambda1SE); BFull(:,idxLambda1SE)]))));
    saturatedLogLikelihood = nansum(log(poisspdf(spkTest,spkTest)));
    devianceFullNull = 2*(fullLogLikelihood - nullLogLikelihood);
    dfFullNull = length(mincoefs);
    if devianceFullNull > chi2inv(1-pThresholdNull, dfFullNull)
        fitResult(1) = 1;
        %% (2) test without each parameter (as a group)
        tempTrainInput = trainingInputMat{planeInd}(:,setdiff(1:size(trainingInputMat{planeInd},2),touchInd));
        tempTestInput = testInputMat{planeInd}(:,setdiff(1:size(testInputMat{planeInd},2),touchInd));
        [BPartialTouch, FitInfoPartialTouch] = lassoglm(tempTrainInput, spkTrain, 'poisson', 'Lambda', logspace(-4,2), 'CV', lambdaCV, 'Options', opt);        
        partialTouchLL = nansum(log(poisspdf(spkTest', exp([ones(size(tempTestInput,1),1), tempTestInput] * [FitInfoPartialTouch.Intercept(FitInfoPartialTouch.Index1SE); BPartialTouch(:,FitInfoPartialTouch.Index1SE)]))));
        devianceFullTouch = 2*(fullLogLikelihood - partialTouchLL);
        dfFullTouch = dfFullNull - FitInfoPartialTouch.DF(FitInfoPartialTouch.Index1SE);
        if devianceFullTouch > chi2inv(1-pThresholdPartial, dfFullTouch)
            fitResult(2) = 1;
        end
        
        tempTrainInput = trainingInputMat{planeInd}(:,setdiff(1:size(trainingInputMat{planeInd},2),soundInd));
        tempTestInput = testInputMat{planeInd}(:,setdiff(1:size(testInputMat{planeInd},2),soundInd));
        [BPartialSound, FitInfoPartialSound] = lassoglm(tempTrainInput, spkTrain, 'poisson', 'Lambda', logspace(-4,2), 'CV', lambdaCV, 'Options', opt);        
        partialSoundLL = nansum(log(poisspdf(spkTest', exp([ones(size(tempTestInput,1),1), tempTestInput] * [FitInfoPartialSound.Intercept(FitInfoPartialSound.Index1SE); BPartialSound(:,FitInfoPartialSound.Index1SE)]))));
        devianceFullSound = 2*(fullLogLikelihood - partialSoundLL);
        dfFullSound = dfFullNull - FitInfoPartialSound.DF(FitInfoPartialSound.Index1SE);
        if devianceFullSound > chi2inv(1-pThresholdPartial, dfFullSound)
            fitResult(3) = 1;
        end
        
        tempTrainInput = trainingInputMat{planeInd}(:,setdiff(1:size(trainingInputMat{planeInd},2),rewardInd));
        tempTestInput = testInputMat{planeInd}(:,setdiff(1:size(testInputMat{planeInd},2),rewardInd));
        [BPartialReward, FitInfoPartialReward] = lassoglm(tempTrainInput, spkTrain, 'poisson', 'Lambda', logspace(-4,2), 'CV', lambdaCV, 'Options', opt);        
        partialRewardLL = nansum(log(poisspdf(spkTest', exp([ones(size(tempTestInput,1),1), tempTestInput] * [FitInfoPartialReward.Intercept(FitInfoPartialReward.Index1SE); BPartialReward(:,FitInfoPartialReward.Index1SE)]))));
        devianceFullReward = 2*(fullLogLikelihood - partialRewardLL);
        dfFullReward = dfFullNull - FitInfoPartialReward.DF(FitInfoPartialReward.Index1SE);
        if devianceFullReward > chi2inv(1-pThresholdPartial, dfFullReward)
            fitResult(4) = 1;
        end
        
        tempTrainInput = trainingInputMat{planeInd}(:,setdiff(1:size(trainingInputMat{planeInd},2),whiskingInd));
        tempTestInput = testInputMat{planeInd}(:,setdiff(1:size(testInputMat{planeInd},2),whiskingInd));
        [BPartialWhisking, FitInfoPartialWhisking] = lassoglm(tempTrainInput, spkTrain, 'poisson', 'Lambda', logspace(-4,2), 'CV', lambdaCV, 'Options', opt);        
        partialWhiskingLL = nansum(log(poisspdf(spkTest', exp([ones(size(tempTestInput,1),1), tempTestInput] * [FitInfoPartialWhisking.Intercept(FitInfoPartialWhisking.Index1SE); BPartialWhisking(:,FitInfoPartialWhisking.Index1SE)]))));
        devianceFullWhisking = 2*(fullLogLikelihood - partialWhiskingLL);
        dfFullWhisking = dfFullNull - FitInfoPartialWhisking.DF(FitInfoPartialWhisking.Index1SE);
        if devianceFullWhisking > chi2inv(1-pThresholdPartial, dfFullWhisking)
            fitResult(5) = 1;
        end
        
        tempTrainInput = trainingInputMat{planeInd}(:,setdiff(1:size(trainingInputMat{planeInd},2),lickInd));
        tempTestInput = testInputMat{planeInd}(:,setdiff(1:size(testInputMat{planeInd},2),lickInd));
        [BPartialLick, FitInfoPartialLick] = lassoglm(tempTrainInput, spkTrain, 'poisson', 'Lambda', logspace(-4,2), 'CV', lambdaCV, 'Options', opt);        
        partialLickLL = nansum(log(poisspdf(spkTest', exp([ones(size(tempTestInput,1),1), tempTestInput] * [FitInfoPartialLick.Intercept(FitInfoPartialLick.Index1SE); BPartialLick(:,FitInfoPartialLick.Index1SE)]))));
        devianceFullLick = 2*(fullLogLikelihood - partialLickLL);
        dfFullLick = dfFullNull - FitInfoPartialLick.DF(FitInfoPartialLick.Index1SE);
        if devianceFullLick > chi2inv(1-pThresholdPartial, dfFullLick)
            fitResult(6) = 1;
        end
    end
    fitResults(cellnum,:) = fitResult;
end


save(savefn, 'fitResults', 'fitInd');
%%
% nanmean(spk)
