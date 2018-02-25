function z = jksbxreadframes_multifile(cellfn, frames, vargin)

% read desired frames (in a vector form) from sbx file
% input should be a form of cell containing filenames (excluding
% extensions) or just a single file name

if ~iscell(cellfn) && ischar(cellfn)
    tempfn = cellfn;
    cellfn = cell(1);
    cellfn{1} = tempfn;
end
    
if nargin < 2
    error('need at least 2 input arguments')
end

if nargin == 3
    maxidx_list = vargin{1};
else
    maxidx_list = [];
    for i = 1 : length(cellfn)        
        maxidx_list = [maxidx_list, maxidx_list(end) + sbx_maxidx(cellfn{i}) + 1]; % each idx starts with 0
    end
end

z = sbxread(cellfn{1},0,1);
z = zeros([size(z,1) size(z,2) size(z,3) length(frames)]);

for i = 1 : length(frames)     
    ifile = length(maxidx_list < frames(i)) + 1;
    if ifile == 1
        z(:,:,:,i) = sbxread(cellfn{1}, frames(i), 1);    
    else
        z(:,:,:,i) = sbxread(cellfn{ifile}, frames(i) - maxidx_list(ifile-1) - 1 , 1);    
    end
end
