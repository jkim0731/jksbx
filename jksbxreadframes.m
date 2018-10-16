function z = jksbxreadframes(fn, frames)

% read desired frames (in a vector form) from sbx file

z = sbxread(fn,0,1);
z = zeros([size(z,2) size(z,3) length(frames)]);

for i = 1:length(frames)
    z(:,:,i) = sbxread(fn,frames(i),1);    
end