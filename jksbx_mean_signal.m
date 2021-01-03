function msignal = jksbx_mean_signal(fn, varargin)

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