function z = jksbxreadrandframes(fname,N,frames)

% read a set of N random frames from the list of "frames"

if nargin < 3
    error('too low input arguments')
end
z = sbxread(fname);
z = zeros([size(z,1) size(z,2) size(z,3) N]);
idx = frames(randperm(length(frames),N));

for j = 1 : length(idx) 
    z(:,:,:,j) = sbxread(fname,idx(j),1);    
end