function [] = tiffaligned_jk(fname, varargin)
% read frames from fname.sbx and make it a stack
% and save to fname_aligned.tif
% if fname_aligned.tif exists, then quit
global info;
temp = sbxread(fname,0,1);
if nargin > 1
    frames = varargin{1};
else 
    frames = 1:info.max_idx;
end

tiffname = strcat(fname,'_aligned.tif');
if exist(tiffname)
    error('tif file exists');
    return
elseif ~exist(strcat(fname,'.align'))
    display('not aligned. align using sbxalignx.');
    tic
    [m,T] = sbxalignx(fn,0:info.max_idx-1);   
    save([fn '.align'],'m','T');
    display(sprintf('Done %s: Aligned %d images in %d min',fn,info.max_idx,round(toc/60)));
end

info.aligned = load([fname ,'.align'],'-mat');
for i = 1:length(frames)
    temp = squeeze(sbxread(fname,frames(i)-1,1));
    temp = circshift(temp,info.aligned.T(frames(i),:));
    temp = temp(:,70:end-70);
    imwrite(temp,tiffname,'WriteMode','append','Compression','none');
    
end
return
end
