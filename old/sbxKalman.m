function sbxKalman(fname, varargin)
% appply Kalman Stack Filtering to .sbx file 
% (How about .align? Think later...)
% from start frame with N consecutive frames
% save as .ksf file '-mat'

% if length(varargin) == 1 then that arg is the total number of frames to be 
% applied, starting from the first frame (0th frame)
% if length(varargin) == 2 then the first argument is the first frame (starting
% from 0) and the second argument is the total number of frames

% 02/11/2016 JK

a = sbxread(fname,0,1);
global info

switch nargin
    case 1
        k = 0;
        N = info.max_idx;
    case 2
        k = 0;
        N = min(varargin{1}, info.max_idx);
    case 3
        k = min(varargin{1}, info.max_idx);
        N = min(varargin{2}, info.max_idx - varargin{1});
    otherwise
        error('wrong input argument')
end
q = squeeze(sbxread(fname,k,N));
ksf = Kalman_Stack_Filter(q);
save([fname, '.ksf'], 'ksf')
end


