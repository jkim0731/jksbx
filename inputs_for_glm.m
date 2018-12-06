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
cellnum = 1;
nShift = 3;

dn = sprintf('%s%03d',baseDir,mouse);
ufn = sprintf('UberJK%03dS%02d.mat', mouse, session);
cd(dn)
% load(ufn)
u = Uber.buildUberArray(mouse, session);

cID = u.cellNums(cellnum);
frameRate = u.frameRate;


% find out trial indices for this specific cell
tind = find(cellfun(@(x) ismember(cID, x.neuindSession), u.trials));
% find out row number of this cell
cind = find(u.trials{tind(1)}.neuindSession == cID);



spk = cell2mat(cellfun(@(x) [nan(1,nShift), x.spk(cind,:), nan(1,nShift)], u.trials(tind)','uniformoutput',false));
pTouchOnset = cell2mat(cellfun(@(x) [nan(1,nShift), histcounts(cellfun(@(y) y(1), x.protractionTouchChunks), [0, x.tpmTime]), nan(1,nShift)], u.trials(tind)','uniformoutput',false));
rTouchOnset = cell2mat(cellfun(@(x) [nan(1,nShift), histcounts(cellfun(@(y) y(1), x.retractionTouchChunks), [0, x.tpmTime]), nan(1,nShift)], u.trials(tind)','uniformoutput',false));
tTouchOnset = pTouchOnset + rTouchOnset;

whiskerVideoFrameDuration = mean(diff(u.trials{tind(1)}.whiskerTime)) * 1000; % in ms
pTouchDuration = cell2mat(cellfun(@(x) [nan(1,nShift), histcounts(cell2mat(cellfun(@(y) y', x.protractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime]), nan(1,nShift)], ...
    u.trials(tind)','uniformoutput',false));
rTouchDuration = cell2mat(cellfun(@(x) [nan(1,nShift), histcounts(cell2mat(cellfun(@(y) y', x.retractionTouchChunks, 'uniformoutput', false)) * whiskerVideoFrameDuration, [0, x.tpmTime]), nan(1,nShift)], ...
    u.trials(tind)','uniformoutput',false));
tTouchDuration = pTouchDuration + rTouchDuration;

pTouchFrames = pTouchDuration;
pTouchFrames(pTouchDuration > 0) = 1;
rTouchFrames = rTouchDuration;
rTouchFrames(rTouchDuration > 0) = 1;
tTouchFrames = tTouchDuration;
tTouchFrames(tTouchDuration > 0) = 1;

scPiezo = cell2mat(cellfun(@(x) [nan(1,nShift), 1, zeros(1,length(x.tpmTime)-1), nan(1,nShift)], u.trials(tind)','uniformoutput',false));
scPoleup = cell2mat(cellfun(@(x) [nan(1,nShift), histcounts(x.poleUpOnsetTime, [0, x.tpmTime]), nan(1,nShift)], u.trials(tind)','uniformoutput',false));
scPoledown = cell2mat(cellfun(@(x) [nan(1,nShift), histcounts(x.poleDownOnsetTime, [0, x.tpmTime]), nan(1,nShift)], u.trials(tind)','uniformoutput',false));
drinkOnset = cell2mat(cellfun(@(x) [nan(1,nShift), histcounts(x.drinkingOnsetTime, [0, x.tpmTime]), nan(1,nShift)], u.trials(tind)','uniformoutput',false));

lLick = cell2mat(cellfun(@(x) [nan(1,nShift), histcounts(x.leftLickTime, [0, x.tpmTime]), nan(1,nShift)], u.trials(tind)','uniformoutput',false));
rLick = cell2mat(cellfun(@(x) [nan(1,nShift), histcounts(x.rightLickTime, [0, x.tpmTime]), nan(1,nShift)], u.trials(tind)','uniformoutput',false));
tLick = lLick + rLick;
tLick(tLick>0) = 1;
%%
whiskingOnsetCell = cell(1,length(tind));
whiskingAmpCell = cell(1,length(tind));
whiskingMidpointCell = cell(1,length(tind));

for ti = 1 : length(tind)
    currTrial = u.trials{tind(ti)};
    time = [0, currTrial.tpmTime];
    theta = nan(currTrial.nof,1);
    wtimes = [0:currTrial.nof-1] * currTrial.frameDuration;
    wInds = round(currTrial.whiskerTime/currTrial.frameDuration)+1;
    theta(wInds) = currTrial.theta;
    [onsetFrame, amplitude, midpoint] = jkWhiskerOnsetNAmplitude(theta, 5);
    onsetTimes = onsetFrame*whiskerVideoFrameDuration/1000; % back to s
    whiskingOnsetCell{ti} = [nan(1,nShift), histcounts(onsetTimes, time), nan(1,nShift)];
    
    tempAmp = zeros(1,length(time)-1);
    tempMid = zeros(1,length(time)-1);
    for i = 1 : length(tempAmp)
        startInd = find(wtimes >= time(i), 1, 'first');
        endInd = find(wtimes < time(i+1), 1, 'last');
        tempAmp(i) = mean(amplitude(startInd:endInd));
        tempMid(i) = mean(midpoint(startInd:endInd));
    end
    whiskingAmpCell{ti} = [nan(1,nShift), tempAmp, nan(1,nShift)];
    whiskingMidpointCell{ti} = [nan(1,nShift), tempMid, nan(1,nShift)];    
end
whiskingOnset = cell2mat(whiskingOnsetCell);
whiskingAmp = cell2mat(whiskingAmpCell);
whiskingMidpoint = cell2mat(whiskingMidpointCell);

%%
corr(spk(isfinite(spk))', tTouchOnset(isfinite(tTouchOnset))')
corr(spk(isfinite(spk))', pTouchOnset(isfinite(pTouchOnset))')
corr(spk(isfinite(spk))', rTouchOnset(isfinite(rTouchOnset))')
corr(spk(isfinite(spk))', tTouchFrames(isfinite(tTouchFrames))')
corr(spk(isfinite(spk))', pTouchFrames(isfinite(pTouchFrames))')
corr(spk(isfinite(spk))', rTouchFrames(isfinite(rTouchFrames))')
corr(spk(isfinite(spk))', tTouchDuration(isfinite(tTouchDuration))')
corr(spk(isfinite(spk))', pTouchDuration(isfinite(pTouchDuration))')
corr(spk(isfinite(spk))', rTouchDuration(isfinite(rTouchDuration))')
corr(spk(isfinite(spk))', scPiezo(isfinite(scPiezo))')
corr(spk(isfinite(spk))', scPoleup(isfinite(scPoleup))')
corr(spk(isfinite(spk))', scPoledown(isfinite(scPoledown))')
corr(spk(isfinite(spk))', drinkOnset(isfinite(drinkOnset))')
corr(spk(isfinite(spk))', whiskingOnset(isfinite(whiskingOnset))')
corr(spk(isfinite(spk))', whiskingAmp(isfinite(whiskingAmp))')
corr(spk(isfinite(spk))', whiskingMidpoint(isfinite(whiskingMidpoint))')
corr(spk(isfinite(spk))', tLick(isfinite(tLick))')
corr(spk(isfinite(spk))', lLick(isfinite(lLick))')
corr(spk(isfinite(spk))', rLick(isfinite(rLick))')

%%
