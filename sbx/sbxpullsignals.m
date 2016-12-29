function sig = sbxpullsignals(fname, varargin)

z = sbxread(fname,1,1);
global info;

load([fname '.segment'],'-mat'); % load segmentation

if nargin > 1
    im_plane = varargin{1};
    plane_num = length(info.aligned.T);
    if exist([fname '.signals'],'file')
        load([fname '.signals'],'-mat')
    else
        sig = cell(1,plane_num);
    end
else
    im_plane = 1;
end


if iscell(mask)
    temp_mask = mask{im_plane};
else
    temp_mask = mask; % make them temp_ because mask is global in sbxsegmenttool
end

ncell = max(temp_mask(:));

for(i=1:ncell)
    idx{i} = find(temp_mask==i);
end

if info.volscan
    sig{im_plane} = zeros(size(info.aligned.T{im_plane},1),ncell); % include the first frame of optotune imaging
else
    sig = zeros(info.max_idx, ncell);
end

h = waitbar(0,sprintf('Pulling %d signals out...',ncell));

if info.volscan    
    for i = 1 : size(sig{im_plane},1)
        waitbar(i/size(sig{im_plane},1),h);
        z = sbxread(fname, (i-1)*plane_num + im_plane, 1); % include the first frame of optotune imaging
        z = squeeze(z(1,:,:));
        z = circshift(z,info.aligned.T{im_plane}(i,:));
        for j = 1 : ncell
            sig{im_plane}(i+1,j) = mean(z(idx{j}));
        end
    end
else
    for(i=0:info.max_idx-1)
        waitbar(i/(info.max_idx-1),h);          % update waitbar...
        z = sbxread(fname,i,1);
        z = squeeze(z(1,:,:));
       z = circshift(z,info.aligned.T(i+1,:)); % align the image
        for(j=1:ncell)                          % for each cell
            sig(i+1,j) = mean(z(idx{j}));       % pull the mean signal out...
        end
    end
end
delete(h);

save([fname '.signals'],'sig');     % append the motion estimate data...


