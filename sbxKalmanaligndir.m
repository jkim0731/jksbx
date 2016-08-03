function sbxKalmanaligndir(varargin)
% appply Kalman Stack Filtering to all .align file in the directory
% save as .ksf file '-mat'
% 02/21/2016 JK

if(nargin>=1) % cell with filenames to be aligned
    for(i=1:length(varargin{1}))
        d(i).name = varargin{1}{i};
    end
else
    d = dir('*.align');
end

% Kalman all *.align files in the list

for(i=1:length(d))
    fn = strtok(d(i).name,'.');
    a = sbxread(fn,0,1);
    global info
    k = 0;
    N = info.max_idx;

    z = squeeze(sbxread(fn,k,N));
    z = zeros(size(z,1), size(z,2), size(z,3));
    for(j=1:length(idx))
        z_j = squeeze(sbxread(fname,idx(j),1));
        z(:,:,j) = circshift(z_j,info.aligned.T(idx(i)+1,:));
    end
    ksf = Kalman_Stack_Filter(z);
    save([fn, '.ksf'], 'ksf')
end