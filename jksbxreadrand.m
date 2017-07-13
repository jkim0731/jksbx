function z = jksbxreadrand(fname,N,frames,varargin)

% read a set of N images skip frames apart and correct for rigid motion

global info;
z = sbxread(fname,1,1);

if nargin < 3
    error('too low input arguments')
elseif nargin == 3    
elseif nargin == 4 % optotune volumetric imaging
    if info.volscan % optotune volumetric imaging
        im_plane = varargin{1};
        plane_num = length(info.otwave);
        frames = unique(floor(frames/plane_num));        
    end
else
    error('too many input arguments')
end

z = zeros([size(z,2) size(z,3) N]);
idx = frames(randperm(length(frames),N));

if info.volscan
    idx = idx(2:end); % to ignore first frame in optotune imaging
    h = waitbar(0,sprintf('Reading %d frames',length(idx)));
    for(j=1:length(idx))
        waitbar(j/length(idx),h);
        q = sbxread(fname,idx(j)*plane_num+im_plane-1,1);
        z(:,:,j) = circshift(squeeze(q(1,:,:)),info.aligned.T{im_plane}(idx(j)+1,:));        
    end
else
    h = waitbar(0,sprintf('Reading %d frames',length(idx)));
    for(j=1:length(idx))
        waitbar(j/length(idx),h);
        q = sbxread(fname,idx(j),1);
        z(:,:,j) = circshift(squeeze(q(1,:,:)),info.aligned.T(idx(j)+1,:));        
    end
end
delete(h);
