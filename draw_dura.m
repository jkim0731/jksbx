function draw_dura(zstackFn)

zstack = [];
load(zstackFn,'zstack','knobbyInfo') % loading zstack
figure('units','normalized','position',[0.1 0.1 0.8 0.8]);
slhan = uicontrol('style','slider','units', 'normalized', 'position', [0.85 0.2 0.01 0.6], 'min',1,'max',size(zstack,3), 'value', 1, ...
    'sliderstep', [1/size(zstack,3), 5/size(zstack,3)],'callback',@slidercb);
hsttext = uicontrol('style','text','units', 'normalized', 'position', [0.4 0.85 0.2 0.05]);
uicontrol('style','pushbutton', 'units', 'normalized', 'position', [0.22 0.82 0.2/3 0.1/3], 'string', 'Clear', 'callback', @clearcb);
uicontrol('style','pushbutton', 'units', 'normalized', 'position', [0.32 0.82 0.2/3 0.1/3], 'string', 'Set', 'callback', @setcb);

uicontrol('style','pushbutton', 'units', 'normalized', 'position', [0.63 0.82 0.2/3 0.1/3], 'string', 'Confirm', 'callback', @confirmcb);
uicontrol('style','pushbutton', 'units', 'normalized', 'position', [0.73 0.82 0.2/3 0.1/3], 'string', 'Save', 'callback', @savecb);
num = round(get(slhan,'value'));
set(hsttext,'visible','on','string',num2str(num))
axes('units','normalized','position',[0.2 0.2 0.6 0.6]);
imshow(mat2gray(zstack(:,:,num))), axis image
points = []; intpoints = cell(size(zstack,3),1);
zstackDura = ones(size(zstack),'logical');
set(gcf, 'WindowButtonDownFcn', @durabdcb);

    function slidercb(source,eventdata)
        num=round(get(slhan,'value'));
        set(hsttext,'visible','on','string',num2str(num))
        imshow(mat2gray(zstack(:,:,num))), axis image
        if ~isempty(points)
            points = [];
        end
        if ~isempty(intpoints{num})
            hold on, plot(intpoints{num}(1,:), intpoints{num}(2,:), 'r.', 'markersize', 5), hold off
        end
    end

    function durabdcb(source,eventdata)
        action_type = get(gcbf, 'SelectionType');
        if strcmp(action_type, 'normal') % left click
            hold on
            p = gca;
            z = round(p.CurrentPoint); 
            z = z(1,1:2);
            plot(z(1),z(2),'r.', 'markersize', 5)
            points = [points, [z(1);z(2)]];
            hold off
        end 
    end

    function clearcb(source, eventdata)
        if ~isempty(points)
            points = [];
            intpoints{num} = [];
        end
        imshow(mat2gray(zstack(:,:,num))), axis image
    end

    function setcb(source, eventdata)
        x = [points(1,:), points(1,1)];
        y = [points(2,:), points(2,1)];
        t = 1:length(x);
        intpoints{num} = zeros(2,length(1:.01:length(x)));
        intpoints{num}(1,:) = interp1(t,x,1:.01:length(x),'pchip');
        intpoints{num}(2,:) = interp1(t,y,1:.01:length(x),'pchip');

        hold on; plot(intpoints{num}(1,:),intpoints{num}(2,:),'g'); hold off        
    end
    
    function confirmcb(source, eventdata)
        if ~isempty(intpoints{num})
        % image compartmentalization
            BW = poly2mask(intpoints{num}(1,:), intpoints{num}(2,:), size(zstack,1), size(zstack,2));
            imshow(BW)
            zstackDura(:,:,num) = BW;
        end
    end


    function savecb(source, eventdata)
        [temp, rest] = strtok(zstackFn,'_');
        savefn = [temp, '_dura', rest];
        tempflist = dir([savefn, '*']);
        if isempty(tempflist)
            save(savefn, 'zstackDura')
        end
    end
end