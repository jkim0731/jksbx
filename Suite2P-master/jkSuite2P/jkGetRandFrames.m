% gets random frames (in total NimgFirstRegistration)
% frames are used for initializing registration
function IMG = jkGetRandFrames(ops)

curr_dir = pwd;
cd(ops.RootDir)

nplanes            = getOr(ops, {'nplanes'}, 1);
ichannel           = getOr(ops, {'gchannel'}, 1);
rchannel           = getOr(ops, {'rchannel'}, 2);
red_align          = getOr(ops, {'AlignToRedChannel'}, 0);

% % if not empty, generate target image from specified experiment (useful if there is slow drift in data)
% targetSession      = getOr(ops, {'targetSession'}, []); 
% Treat this after BiDiPhase for every session. Each session might have different bidirection phase setting JK 2018/02/22
min_nframes = min(cellfun(@(x) length(x), ops.frame_to_use));
ops.NimgFirstRegistration = min(min_nframes, ops.NimgFirstRegistration);
IMG = zeros(ops.Ly, ops.Lx, nplanes, ops.NimgFirstRegistration, 'single');
% grab frames from all files
for iplane = 1 : nplanes
    img_plane = permute(jksbxreadrandframes_multifile(ops.sbxfnlist, ops.NimgFirstRegistration, ops.frame_to_use{iplane}),[2 3 1 4]);
    img_plane = img_plane(ops.useY, ops.useX, :, :);
    if red_align
        IMG(:,:,iplane,:) = img_plane(:,:,rchannel,:);
    else
        IMG(:,:,iplane,:) = img_plane(:,:,ichannel,:);
    end
end
cd(curr_dir)