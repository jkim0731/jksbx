function ops1 = jk_reg2P(ops)
%% find the mean frame after aligning a random subset
ops.doRegistration  = getOr(ops, {'doRegistration'}, 1); % register tiffs?

if ops.doRegistration
    disp('running registration');
else
    disp('skipping registration, but assembling binary file');
end

ops = buildRegOps(ops);

% --- number of planes in recording and number of channels --- %
nplanes            = ops.nplanes;
% numPlanes removed because it is redundant with nplanes, unless planesToProcess ~= 1:nplanes, but cannot think of this occasion for now 2018/02/21 JK
% which channel is the functional channel
ichannel           = getOr(ops, {'gchannel'}, 1);
% which channel is the non-functional channel
rchannel           = getOr(ops, {'rchannel'}, 2);

% --- red channel options ---%
red_align          = getOr(ops, {'AlignToRedChannel'}, 0); % register planes to red channel
red_binary         = getOr(ops, {'REDbinary'}, 0); % write red channel to a binary file
% extract mean red channel from blocks of recording with two channels
red_mean           = getOr(ops, {'redMeanImg'}, 0);
if red_mean 
    disp('computing mean RED image');
else
    disp('no red channel')
end

BiDiPhase          = ops.BiDiPhase;
%% find the mean frame after aligning a random subset

ops.Ly = length(ops.useY);
ops.Lx = length(ops.useX);

% split into subsets (for high scanning resolution recordings)
[xFOVs, yFOVs] = get_xyFOVs(ops);

targetSession      = getOr(ops, {'targetSession'}, []);

if ops.doRegistration
    % get frames for initial registration
    IMG = jkGetRandFrames(ops); % ops.useX and useY are already treated inside jkGetRandFrames
    % compute phase shifts from bidirectional scanning    
    if ops.dobidi
        ops.BiDiPhase = BiDiPhaseOffsets(IMG(:,:,1,:)); % BiDi is same across all planes in each session        
    end
    BiDiPhase = ops.BiDiPhase;    
    fprintf('bi-directional scanning offset = %d pixels for session %d\n', BiDiPhase, ops.session);
    if abs(BiDiPhase) > 0
       IMG = ShiftBiDi(BiDiPhase, IMG, ops.Ly, ops.Lx);
    end

    % makes blocks (number = numBlocks) and masks for smoothing registration offsets across blocks
    % this is for nonrigid
    if ops.nonrigid
        ops = MakeBlocks(ops);
    end
    
    % for each plane: align chosen frames to average to generate target image
    ops1 = cell(nplanes,size(xFOVs,2));
    for iplane = 1 : nplanes
        if ops.nonrigid && ~(ops.alignAcrossPlanes || ops.interpolateAcrossPlanes)
            ops1{iplane} = nonrigidAlignIterative(squeeze(IMG(:,:,ops.planesToProcess(iplane),:)), ops);
        else
            for iFOV = 1:size(xFOVs,2)
                ops1{iplane,iFOV} = alignIterative(single(squeeze(IMG(yFOVs(:,iFOV),xFOVs(:,iFOV),ops.planesToProcess(iplane),:))), ops);
            end
        end
        fprintf('target image acquired for plane %d/%d of mouse %s session %d\n', ops.planesToProcess(iplane), nplanes, ops.mouse_name, ops.session);
    end
    
    if ops.alignTargetImages % align target images of all planes to each other
        % (reasonable only if interplane distance during imaging was small)        
        newTargets = getImgOverPlanes(ops1, ops.planesToInterpolate);
        for i = 1:length(ops.planesToInterpolate)
            for l = 1:size(xFOVs,2)
                ops1{ops.planesToInterpolate(i),l}.mimg = newTargets(:,:,i,l);
            end
        end
% JK This is aligning between different depths (e.g., different optotune planes)
    end
    
    % display target image
    if ops.fig && ops.showTargetRegistration
        jkPlotRegMean(ops1,ops);
        drawnow
    end
    clear IMG
    
else   % don't recompute mean image
    ops1 = cell(nplanes,1);
    for iplane = 1 : nplanes
        ops1{iplane} = ops;
        ops1{iplane}.mimg = zeros(ops.Ly, ops.Lx);
    end
end

%% initialize mean imgs and shifts

for iplane = 1 : nplanes
    if ops.nonrigid
        ops1{iplane}.mimgB = cell(prod(ops.numBlocks),1);
        for ib = 1:ops.numBlocks(1)*ops.numBlocks(2)
            ops1{iplane}.mimgB{ib} = ops1{iplane}.mimg(ops1{iplane}.yBL{ib}, ops1{iplane}.xBL{ib});
        end
    else
        for iFOV = 1:size(xFOVs,2)
            if red_mean || red_align
                ops1{iplane,iFOV}.mimgRED       = zeros(ops1{iplane,iFOV}.Ly, ops1{iplane,iFOV}.Lx);
            end
            ops1{iplane,iFOV}.DS                = [];
            ops1{iplane,iFOV}.CorrFrame         = [];
            ops1{iplane,iFOV}.mimg1             = zeros(ops1{iplane,iFOV}.Ly, ops1{iplane,iFOV}.Lx);
        end
    end
end


%% open files for registration

[ops1, fid, fidRED, fidIntpol] = jkOpenBinFiles(ops, ops1, red_binary);

%%
tic
% compute registration offsets and align using offsets

xyValid = true(ops.Ly,ops.Lx);
if red_align
    reg_channel = rchannel;
else
    reg_channel = ichannel;
end

% in case different experiments have different numbers of channels     
for iplane = 1 : nplanes
    % only load frames of registration channel
    data = jkLoadFramesBuff(ops, iplane, reg_channel);
    data = data(ops.useY, ops.useX, :, :);
    if abs(BiDiPhase) > 0; data = ShiftBiDi(BiDiPhase, data, ops.Ly, ops.Lx); end

    % get the registration offsets for each frame
    if ops.doRegistration
        if ops.nonrigid
            [dsall, ops1{iplane}] = jkNonrigidOffsets(data, ops, ops1{iplane});
        else                
            dsall = zeros(size(data,3),2,size(xFOVs,2));
            for iFOV = 1 : size(xFOVs,2)
                [dsall(:,:,iFOV), ops1{iplane,iFOV}] = jkRigidOffsets(data(yFOVs(:,iFOV),xFOVs(:,iFOV),:), ops, ops1{iplane,iFOV});
            end                
        end
    end

    if ops.doRegistration
        % if aligning by the red channel, data needs to be reloaded as the green channel
        if red_align
            data_red = data;
            data = jkLoadFramesBuff(ops, iplane, ichannel);
            data = data(ops.useY, ops.useX, :, :);
            if abs(BiDiPhase) > 0; data = ShiftBiDi(BiDiPhase, data, ops.Ly, ops.Lx);  end
        end
        % load red channel for red binary file if not already loaded for red alignment
        if red_mean && ~red_align
            data_red = jkLoadFramesBuff(ops, iplane, rchannel);
            if abs(BiDiPhase) > 0; data_red = ShiftBiDi(BiDiPhase, data_red, ops.Ly, ops.Lx); end
        end
        % align the frames according to the registration offsets
        if ops.nonrigid
            [dreg, xyValid] = nonrigidMovie(data, ops, dsall, xyValid);
        else
            dreg = rigidMovie(data, ops1, dsall, yFOVs, xFOVs); % Just need to pass ops.useGPU, but to comply with previous code, pass all of ops1.
        end

        if red_mean || red_align
            if ops.nonrigid
                [dreg2, xyValid] = nonrigidMovie(data_red, ops, dsall, xyValid);
            else
                dreg2 = rigidMovie(data_red, ops1, dsall, yFOVs, xFOVs);
            end
        end
    else
        dreg = data;
    end

    %%
    % write dreg to bin file+
    for iFOV = 1:size(xFOVs,2)
        dwrite = dreg(yFOVs(:,iFOV),xFOVs(:,iFOV),:);
        fwrite(fid{iplane,iFOV}, dwrite, class(data));
        ops1{iplane,iFOV}.Nframes = size(dwrite,3);
        ops1{iplane,iFOV}.mimg1 = ops1{iplane,iFOV}.mimg1 + sum(dwrite,3);

        if red_mean || red_align
            dwrite = dreg2(yFOVs(:,iFOV),xFOVs(:,iFOV),:);
            ops1{iplane,iFOV}.mimgRED = ops1{iplane,iFOV}.mimgRED + sum(dwrite,3);
        end
        if red_binary
            fwrite(fidRED{iplane,iFOV}, dwrite, class(data));
        end
    end

    fprintf('Plane #%d/%d of session %d done in time %2.2f \n', iplane, nplanes, ops.session, toc)    
end

if ops.alignAcrossPlanes && ops.doRegistration % align each frame with the best matching target image
    if ops.nonrigid
        [ops1, planesToInterp] = registerBlocksAcrossPlanes(ops1, ops, fid, fidIntpol);
    else
        [ops1, planesToInterp] = registerAcrossPlanes(ops1, ops, fid, fidIntpol);
    end
end

for i = 1:numel(ops1)
    ops1{i}.mimg1 = ops1{i}.mimg1/sum(ops1{i}.Nframes);
    if red_mean || red_align
        red_expts = ismember(ops.session, getOr(ops, 'expred', []));
        ops1{i}.mimgRED = ops1{i}.mimgRED/sum(ops1{i}.Nframes(red_expts));
        if red_binary
            fclose(fidRED{i});
        end
    end
    ops1{i}.badframes = false(1, size(ops1{i}.DS,1));
    if isfield(ops, 'badframes0') && ~isempty(ops.badframes0)
        ops1{i}.badframes(ops.badframes0) = true;
    end
    fclose(fid{i});
    if ops.interpolateAcrossPlanes
        fclose(fidIntpol{i});
    end
end

%%
% get mean of first and last frames in block (to check for drift)
if ~isempty(ops.nimgbegend) && ops.nimgbegend>0
    for i = 1:numel(ops1)
        fid{i}      = fopen(ops1{i}.RegFile, 'r');
        ops1{i}     = getBlockBegEnd(fid{i}, ops1{i}); % get mean of first and last frames in block (to check for drift)
        fclose(fid{i});
    end
end


% write registered tiffs to disk if ~isempty(ops.RegFileTiffLocation)
if ~isempty(ops.RegFileTiffLocation)
    for i = 1:numel(ops1)
        fid{i}      = fopen(ops1{i}.RegFile, 'r');
        ops1{i}     = jk_write_reg_to_tiff(fid{i}, ops1{i}, i, 0);
        fclose(fid{i});

        if red_binary
            fidRED{i}   = fopen(ops1{i}.RegFile2, 'r');
            ops1{i}     = jk_write_reg_to_tiff(fidRED{i}, ops1{i}, i, 1);
            fclose(fidRED{i});
        end   
    end
end

% copy binary file if ~isempty(ops.RegFileBinLocation)
if ~isempty(ops.RegFileBinLocation)
    copyBinFile(ops1);
end

if ops.interpolateAcrossPlanes == 1 && ~isempty(RegFileBinLocation)
    for i = 1 : numel(ops1)
        folder = fullfile(ops1{i}.RegFileBinLocation, ops1{i}.mouse_name, ...
            ops1{i}.session, 'interpolated');
        filename = fullfile(folder, ...
            sprintf('%s_%03d_plane%d.bin', ops.mouse_name, ops.session, i));
        if ismember(i, planesToInterp)
            fidCopy = fopen(ops1{i}.RegFile, 'w');
            fidOrig = fopen(filename, 'r');
            sz = ops1{i}.Lx * ops1{i}.Ly;
            parts = ceil(sum(ops1{i}.Nframes) / 2000);
            for p = 1:parts
                toRead = 2000;
                if p == parts
                    toRead = sum(ops1{i}.Nframes) - 2000 * (parts-1);
                end
                data = fread(fidOrig,  sz*toRead, '*int16');
                fwrite(fidCopy, data, class(data));
            end
            fclose(fidCopy);
            fclose(fidOrig);
        else
            delete(filename)
        end
    end
end

%%
% compute outlier frames and xrange, yrange of registered frames
if ops.doRegistration
    ops1 = getRangeNoOutliers(ops, ops1);
else
    for i = 1:numel(ops1)
        ops1{i}.yrange = 1:ops.Ly;
        ops1{i}.xrange = 1:ops.Lx;
    end
end


savepath = sprintf('%s/', ops.ResultsSavePath);
if ~exist(savepath, 'dir')
    mkdir(savepath)
end
save(sprintf('%s/regops_%s_%03d.mat', ops.ResultsSavePath, ...
    ops.mouse_name, ops.session),  'ops1')
end
