function [optimizer, metric] = align_jk(fname)
% % align using intensity-based registration
% % For now, let's use tiff
% % output as another tiff
% tiffname = strcat(fname,'.tif');
% info = imfinfo(tiffname);
% num_images = numel(info);
% 
% tiffalign_name = strcat(fname,'_align.tif');
% if exist(tiffalign_name)
%     error('already aligned.');
%     return
% else
%     fixed = imread(tiffname,1);
%     imwrite(fixed,tiffalign_name,'Compression','none');
% 
%     for i = 2 : num_images
%         temp = imread(tiffname,i);
%         [optimizer, metric]  = imregconfig('monomodal');
%         registered = imregister(temp,fixed,'rigid',optimizer, metric);
%         imwrite(registered,tiffalign_name,'WriteMode','append','Compression','none');
%     end
%     return
% end

% changed to using sbx file 2016/05/14
global info
fixed = squeeze(sbxread(fname,0,1)); % Assume now for having just 1 color;
for i = 1 %: info.max_idx -1
    temp = squeeze(sbxread(fname,i,1));
    [optimizer, metric] = imregconfig('monomodal');
    registered = imregister(temp,fixed,'rigid',optimizer,metric);
end