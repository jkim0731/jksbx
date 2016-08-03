function manual = roi_manual(mouse_num, reffn)

sn = [num2str(mouse_num), '_roi.mat'];
if (exist([reffn,'.alignref'],'file'))
    load([reffn, '.alignref'],'-mat');
elseif (exist([reffn,'.align'],'file'))
    load([reffn, '.align'],'-mat');
else
    aligndir({reffn});
    load([reffn, '.align'],'-mat');
end
figure, imshow(m), axis off, axis image, hold all, imcontrast

manual = []; i = 1;
if exist([pwd, '\', sn])
    load(sn, 'manual')
    i = length(manual) + 1;
    scatter(manual(:,2), manual(:,1), 'go');
end

while (1)
    [y, x] = ginput(1);
    if y < 0 || x < 0
        break
    else
        scatter(y, x, 'mo');
    end
    manual(i,2)= y; manual(i,1) = x;
    i = i + 1;
end

if exist([pwd, '\', sn])
    save(sn,'manual', 'append')
else
    save(sn,'manual')
end
