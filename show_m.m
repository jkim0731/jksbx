function show_m(fn)
% load variable m from .align and show the image
load([fn,'.align'],'-mat')
% imtool(m(:,52:end-52))
imtool(m)
end