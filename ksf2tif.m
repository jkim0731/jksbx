function ksf2tif(fname,varargin)

% sbx2tif
% Generates tif file from sbx files
% Argument is the number of frames to convert
% If no argument is passed the whole file is written

z = sbxread(fname,1,1);
global info;

if(nargin>1)
    N = min(varargin{1},info.max_idx);
else
    N = info.max_idx;
end

load([fname, '.ksf'], '-mat', 'ksf');

tifname = [fname '_ksf.tif'];
options.ask = true;
options.message = true;
% options.big = true; % Use BigTIFF format
saveastiff(ksf, tifname, options);

end