function z = jksbxreadrandframes_multifile(cellfn,N,frames,maxidx_list)

% read a set of N random frames from the list of "frames"
% input should be a form of cell containing filenames (excluding
% extensions)

if ~iscell(cellfn) && ischar(cellfn)
    tempfn = cellfn;
    cellfn = cell(1);
    cellfn{1} = tempfn;
end

if nargin < 4
    error('need at least 4 input arguments')
end 


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