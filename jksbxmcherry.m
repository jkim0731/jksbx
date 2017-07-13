function jksbxmcherry(fn)

a = sbxread(fn,0,1);
if size(a,1) == 1
    error('the file should have 2 channels recorded simultaneously');
    return
end
maxidx = jkget_maxidx(fn);
temp = sbxread(fn,0,maxidx);
output = zeros(size(temp,2),size(temp,3)-99,2);
output(:,:,1) = squeeze(sum(uint32(temp(1,:,100:end,:)),4));
output(:,:,2) = squeeze(sum(uint32(temp(2,:,100:end,:)),4));

merge = imfuse(output(:,:,1),output(:,:,2),'falsecolor','ColorChannels',[2 1 0]);
ronly = imfuse(output(:,:,1),zeros(size(output,1),size(output,2)),'falsecolor','ColorChannels',[2 1 0]);
gonly = imfuse(zeros(size(output,1),size(output,2)),output(:,:,2),'falsecolor','ColorChannels',[2 1 0]);
f = figure; 
imshow(merge)

r = 1;
g = 1;

while(1)
    waitforbuttonpress
    fpos = get(f,'Position');
    ha = axes('Parent',f,'Units','pixels',...
            'Position',fpos);
    click_type=get(f,'SelectionType');
    if strcmp(click_type,'normal') %right click
        if g == 1
            imshow(ronly, 'Parent', ha)
            g = 0;
        else
            imshow(merge, 'Parent', ha)            
            g = 1;
        end
    elseif strcmp(click_type,'alt') %left click
        if r == 1
            imshow(gonly, 'Parent', ha)
            r = 0;
        else
            imshow(merge, 'Parent', ha)
            r = 1;
        end
    end
end

end