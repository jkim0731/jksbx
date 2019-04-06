function assign_depth(suite2pFn, zstackFn, zstackDuraFn)
% % input: suite2p file, z-stack file
% % output: create a file with combination of suite2p filename and z-stack
% % file name (since there can be different z-stack file. Save all data from
% % _proc suite2p file and assign cortical depth and whether it is in C2 or
% % not.
% 
% % Cortical surface (determined by autofluorescence from dura mater) and C2
% % will be drawn manually. 
% 
% % Based on the assumption that width of the reference image is longer that
% % its height. (size(refIm,2) > size(refIm,1))
% 
% % %% Dura marking (regions where cortical depth is 0)
% % [temp, rest] = strtok(zstackFile,'_');
% % savefn = [temp, '_dura', rest];
% % tempflist = dir([savefn, '*']);
% % if isempty(tempflist)
% %     draw_dura(zstackFile)
% % end
% 
load(suite2pFn) % loading dat
load(zstackFn) % loading zstack, knobbyInfo, and zstack_maxproj

% 
% 
% 
% 
% 
% 
% 
% %%
% zstackDuraFn = 'zstack_dura_030_1000';
% assign depth allocation
load(zstackDuraFn)
% height = size(zstackDura,1); width = size(zstackDura,2);
if ~exist('tiltAngle', 'var') || ~exist('tiltGrad', 'var')
    [tiltAngle, tiltGrad] = dura_tilt_calc(zstackDuraFn); % it has save function at the end. Takes more than a minute.     
end

zstepum = abs(mean(diff([knobbyInfo(:).z]))) / cosd(knobbyInfo(1).a); % assuming knobby angle did not change during z-stack imaging (must be so)
tiltGradUmperpix = tiltGrad*zstepum; % um/pix
centerPoint = round(dat.ops.info.sz/2); % (1) is for y, (2) is for x.
depthCompensation = zeros(dat.ops.info.sz);
height = dat.ops.info.sz(1); width = dat.ops.info.sz(2);

isoAngle = tiltAngle-90;
if isoAngle ~= 0 && isoAngle ~= 180 % move the line left or right
    initLineY = 1 : height;
    if abs(isoAngle) ~= 90
        initLineX = round( centerPoint(2) + (initLineY-centerPoint(1)) / tand(isoAngle) );
    else
        initLineX = ones(size(initLineY)) * centerPoint(2);
    end
%     zperpix = -tiltGradUmperpix / cosd(tiltAngle); % for moving to right, meaning increasing x
    zperpix = -tiltGradUmperpix * cosd(tiltAngle); % for moving to right, meaning increasing x
    for i = 1 : height
        if initLineX(i) >= 1 % left of the initial line
            depthCompensation(i,1:min(initLineX(i),width)) = ([1:min(initLineX(i),width)] - initLineX(i)) * zperpix;
        end
        if initLineX(i) <= width % right of the initial line
            depthCompensation(i,max(initLineX(i),1):width) = ([max(initLineX(i),1):width] - initLineX(i)) * zperpix;
        end
    end
else % move up or down
    if tiltAngle == 90
        zperpix = -tiltGradUmperpix;
    else % tiltAngle == 270
        zperpix = tiltGradUmperpix;
    end
    for i = 1 : height
        depthCompensation(i,:) = (i - centerPoint(1)) * zperpix;
    end        
end
%%
% figure, subplot(1,2,1), imagesc((sum(zstackDura,3) - sum(zstackDura(centerPoint(1),centerPoint(2),:),3)) * zstepum), axis image,
% subplot(1,2,2), imagesc(depthCompensation), axis image

dat.depth.depthComp = depthCompensation;
dat.depth.tiltAngle = tiltAngle;
xyUmperpix = 1.4 / str2double(dat.ops.info.config.magnification_list(dat.ops.info.config.magnification,:));
dat.depth.xyUmperpix = xyUmperpix;
dat.depth.tiltGradUmperpix = tiltGradUmperpix;
% calculate the window angle (relative orthogonal to the focal plane (35 degrees))
dat.depth.windowAngle = atand(tiltGradUmperpix/xyUmperpix);
% 
% 
% 
% 
% %%
% 
% 
% 
% 
% 
% 
% % assign C2 in the plane
% load(C2Fn)
% 
% save
% if contains(suite2pFn,'.mat')
%     suite2pFn = suite2pFn(1:end-4);
% end
% saveFn = [suite2pFn, '_final.mat'];
saveFn = suite2pFn;
save(saveFn, 'dat', 'spk')