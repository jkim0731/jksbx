function [z, maxidx_list] = jksbxreadrandframes_multifile(cellfn,N,frames,vargin)

% read a set of N random frames from the list of "frames"
% input should be a form of cell containing filenames (excluding
% extensions)

if ~iscell(cellfn) && ischar(cellfn)
    tempfn = cellfn;
    cellfn = cell(1);
    cellfn{1} = tempfn;
end

if nargin < 3
    error('need at least 3 input arguments')
end
    
if nargin == 4
    maxidx_list = vargin{1};
else
    maxidx_list = 0;
    for i = 1 : length(cellfn)        
        maxidx_list = [maxidx_list, maxidx_list(end) + sbx_maxidx(cellfn{i}) + 1]; % each idx starts with 0
    end
end
maxidx_list(1) = [];
z = sbxread(cellfn{1},0,1);
z = zeros([size(z,1) size(z,2) size(z,3) N],'uint16');
idx = frames(randperm(length(frames),N));

for j = 1 : length(idx)     
    ifile = length(find(maxidx_list < idx(j))) + 1;
    if ifile == 1
        z(:,:,:,j) = sbxread(cellfn{1}, idx(j), 1);    
    else
        z(:,:,:,j) = sbxread(cellfn{ifile}, idx(j) - maxidx_list(ifile-1) - 1 , 1);    
    end
end