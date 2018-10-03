function [tiltAngle, tiltGrad] = dura_tilt_calc(fn)
load(fn) % loading zstackDura. 1's are cortex, 0's are outside of the cortex
height = size(zstackDura,1); 
width = size(zstackDura,2);
boxLine = [...
    1   1       width   width   1; ...
    1   height  height  1       1];
tilt = sum(zstackDura,3); % lower the value, lower the surface (shallower in the cortex)
numTrials = 1000;
angles = zeros(numTrials,1, 'uint8');
gvals = zeros(numTrials,1);
r = sqrt(width^2 + height^2);
for i = 1 : numTrials
    % randomize starting point. Limit this to the center quarter.
    subi = randi(round(height/2)) + round(height/2);
    subj = randi(round(width/2)) + round(width/2);
    startPoint = [subj; subi];
    grad = zeros(180,1);

    for j = 1 : 180
        tempTipUp = [subj + r * cosd(j);     subi + r * sind(j)]; % going upward in plot coord, so going downward in image 
        tempTipDown = [subj + r * cosd(j+180); subi + r * sind(j+180)]; % downward, so upward in image 
        edgeIndUp = InterX([startPoint,tempTipUp], boxLine);
        edgeIndDown = InterX([startPoint,tempTipDown], boxLine);
        lineLength = sqrt((edgeIndUp(1,1) - edgeIndDown(1,1)) ^2 + (edgeIndUp(2,1) - edgeIndDown(2,1)) ^2);
        edgeValDiff = tilt(round(edgeIndUp(2,1)), round(edgeIndUp(1,1))) - tilt(round(edgeIndDown(2,1)), round(edgeIndDown(1,1)));
        grad(j) = edgeValDiff/lineLength;
    end
    [gvals(i), deg] = max(abs(grad));
    if grad(deg) > 0 
        angles(i) = deg+180;
    else
        angles(i) = deg;
    end
end

tiltAngle = round(mean(angles));
inds = find(angles == tiltAngle);
if length(inds) < 10
    inds = find(angles - tiltAngle < 2);
end
tiltGrad = mean(gvals(inds));
save(fn, 'zstackDura', 'tiltAngle', 'tiltGrad')