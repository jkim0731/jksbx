function [] = tiffstack_jk(fname, frames)
% read frames from fname.sbx and make it a stack
% and save to fname.tif
% if fname.tif exists, then quit


tiffname = strcat(fname,'.tif');
if exist(tiffname)
    error('tif file exists');
    return
else
    
% sz = size(squeeze(sbxread(fname,0,1)));
% resultstack = zeros(sz(1),sz(2),length(frames));
for i = 1:length(frames)
%     resultstack(:,:,i) = squeeze(sbxread(fname,frames(i)-1,1)); % sbxread reads from 0 (the first frame is numbered as 0)
    temp = squeeze(sbxread(fname,frames(i)-1,1));
    imwrite(temp,tiffname,'WriteMode','append','Compression','none');
    
end
return
end
