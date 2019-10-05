function onFrames = laser_on_frames (fn, varargin)
if nargin > 1
    if isinteger(varargin{1}) && varargin{1} > 1
        numThresh = varargin{1};
    else
        disp('Frame threshold number is wrong. Second input should be integer > 1. numThresh = 1000.')
        numThresh = 1000;
    end
else
    numThresh = 1000;
end

maxIdx = round(jkget_maxidx(fn));
if maxIdx > numThresh
    numChunks = ceil(maxIdx/numThresh); % dividing into chunks just because of memory issue.
end
msignal = zeros(maxIdx,1);
%%
global info
for i = 1 : numChunks-1
    a = sbxread(fn,(i-1)*numThresh,numThresh);
    a = squeeze(a(1,:,:,:));
    msignal((i-1)*numThresh+1:i*numThresh) = mean(mean(a));
end

a = sbxread(fn, (numChunks-1)*numThresh, maxIdx - (numChunks-1) * numThresh);
a = squeeze(a(1,:,:,:));
msignal((numChunks-1)*numThresh+1:maxIdx) = mean(mean(a));

laserOffInd = find(msignal < min(msignal) + 50); % Start by considering total blanked frames. 50 is arbitrary (safe threshold for dark noise).
if ~isempty(find(diff(laserOffInd)>1, 1))
    laserOffStartInds = [laserOffInd(1); laserOffInd(find(diff(laserOffInd)>1)+1)];
    laserOffEndInds = [laserOffInd(diff(laserOffInd)>1); laserOffInd(end)];
else
    laserOffStartInds = laserOffInd(1);
    laserOffEndInds = laserOffInd(end);
end

if info.volscan
    numPlanes = length(info.otwave);
else
    numPlanes = 1;
end

for i = 1 : length(laserOffStartInds)
    if laserOffStartInds(i) > 1
        laserOffStartInds(i) = laserOffStartInds(i)-1; % 1 frame reduction, to include half-blanking. 
        % More conservative than directly calculate partly blanked frame by intensity.
        % There can be false positive (prob. ~ 1/length(lines), usually 1/512)
    end
    if numPlanes > 1 % if there are more than one planes (volumetric scanning), 
        % treat each volume as one. Always start at the beginning of each volume.
        if mod(laserOffStartInds(i)-1, numPlanes)            
            laserOffStartInds(i) = laserOffStartInds(i) - mod(laserOffStartInds(i)-1, numPlanes);          
            if laserOffStartInds(i) < 1
                laserOffStartInds(i) = 1;
            end
        end
    end
end
for i = 1 : length(laserOffEndInds)
    if laserOffEndInds(i) < maxIdx
        laserOffEndInds(i) = laserOffEndInds(i)+1; % similar treatment as in laserOffStartInds
    end
    if numPlanes > 1 % if there are more than one planes (volumetric scanning), 
        % treat each volume as one. Always end at the end of each volume.
        if mod(laserOffEndInds(i), numPlanes)
            laserOffEndInds(i) = laserOffEndInds(i) + numPlanes - mod(laserOffEndInds(i), numPlanes);
            if laserOffEndInds(i) > maxIdx
                laserOffEndInds(i) = maxIdx;
            end
        end
    end
end

onFrames = 1:maxIdx;
for i = 1 : length(laserOffStartInds)
    onFrames = setdiff(onFrames,laserOffStartInds(i):laserOffEndInds(i));
end
