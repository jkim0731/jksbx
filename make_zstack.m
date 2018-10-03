function [zstack, zstack_maxproj] = make_zstack(fn_base, varargin)

% output zstack (image height, image width, z-stack slices, number of
% channels)
% squeezed if number of channels == 1
% save the file

savefn = ['zstack_', fn_base];
fnlist = dir([fn_base, '*.sbx']);
load([fnlist(1).name(1:end-4), '.mat']) % loading an example info to check the number of channels (and which channel it was)

if nargin > 1
    channels = varargin{1}; % either 1, 2, or [1,2]
    if channels == 1 || channels == 2 || channels == [1,2]
    else
        error('Wrong input of channels: should be 1, 2, or [1,2]')
    end
else
    switch info.channels
        case 1
            channels = [1,2];
        case 2
            channels = 1;
        case 3
            channels = 2;
        otherwise
            error('error in info.channels')
    end    
end

if channels == 1 && info.channels == 3
    error('no green channel in this z-stack')
end
if channels == 2 && info.channels == 2
    error('no red channel in this z-stack')
end
if length(channels) == 2 && info.channels ~= 1
    error('no two channels in this z-stack')
end

zstack = zeros([info.sz(1),info.sz(2),length(fnlist),length(channels)],'uint16');
knobbyInfo = struct;

for i = 1 : length(fnlist)
    load([fnlist(i).name(1:end-4), '.mat']) % loading info
    fname = strsplit(fnlist(i).name(1:end-4),'_');
    knobbyInfo(i).fname = str2double(fname{end});
    knobbyInfo(i).x = info.config.knobby.pos.x;
    knobbyInfo(i).y = info.config.knobby.pos.y;
    knobbyInfo(i).z = info.config.knobby.pos.z;
    knobbyInfo(i).a = info.config.knobby.pos.a;
    align_fn = [fnlist(i).name(1:end-4), '.align']; 
    if ~exist(align_fn,'file')
        jksbxaligndir({align_fn(1:end-4)}) % running alignment if not done yet.
    end
    load(align_fn, '-mat') % loading m1 and m2
    if length(channels) == 2
        zstack(:,:,i,1) = m1{1};
        zstack(:,:,i,2) = m2{1};
    elseif channels == 1
        zstack(:,:,i,1) = m1{1};
    else % channels == 2
        zstack(:,:,i,1) = m2{1};
    end
end
zstack_maxproj = zeros([size(zstack,1),size(zstack,2),length(channels)], 'like', zstack);
for i = 1 : size(zstack,4)
    zstack_maxproj(:,:,i) = align_zstack(zstack(:,:,:,i));
end
zstack = squeeze(zstack);
zstack_maxproj = squeeze(zstack_maxproj);
save(savefn, 'zstack', 'zstack_maxproj', 'knobbyInfo')
        
